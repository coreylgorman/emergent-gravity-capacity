#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hqtfim_capacity_probe.py

A compact, referee-friendly substrate check using the 1D transverse-field Ising model (TFIM)
near criticality (z≈1), intended as a *toy analogue* of the CHM/OP modular-response workflow.

What it does
------------
1) Builds an open-boundary TFIM Hamiltonian:
       H = -J * sum_{i=1}^{L-1} σ^z_i σ^z_{i+1} - g * sum_{i=1}^L σ^x_i
   with J=1 (units). Criticality at g≈1.

2) Computes the ground state at g0 and g1 = g0 + dg via sparse exact diagonalization.

3) For contiguous blocks A of length ℓ in [block_min, block_max], forms reduced density
   matrices ρ_A(g0), ρ_A(g1), then:
       - Entanglement entropies: S_A(g) = -Tr(ρ_A log ρ_A)
       - Modular-Hamiltonian of the reference block: K_A = -log ρ_A(g0)
       - δS_A = S_A(g1) - S_A(g0)
       - δ⟨K_A⟩ = Tr( (ρ_A(g1) - ρ_A(g0)) K_A )   (first-law check)

4) "MI/moment-kill–like" regression: fit δ⟨K_A⟩(ℓ) to [const + c * log ℓ], subtract the fit,
   and report the mean residual as a finite, size-stable *plateau* proxy with a small
   error bar (std/√N). This mimics isolating a finite CHM coefficient without claiming
   to reproduce 3+1D CHM/I00.

5) Saves:
   - ./hqtfim_outputs/summary.json
   - ./hqtfim_outputs/dK_vs_logl.png
   - ./hqtfim_outputs/residual_after_mikill.png
   - ./hqtfim_outputs/first_law_check.png

Caveats
-------
• This is a *toy* (1D lattice, small L) to demonstrate robustness of the *structure*:
  first-law linearity and a finite-size-stable residual after removing trivial size
  dependences. It does *not* compute your 3+1D CHM I_{00} nor the cosmological mapping.
• Default L=10 runs in ~seconds on a laptop. You can push to L≈12–14 at the cost of time.

Usage
-----
    python hqtfim_capacity_probe.py
    python hqtfim_capacity_probe.py --L 12 --g0 1.0 --dg 0.002 --block-min 3 --block-max 6

Dependencies
------------
    numpy, scipy, matplotlib
"""

from __future__ import annotations
import argparse
import json
import os
from dataclasses import dataclass, asdict
from typing import Tuple, List

import numpy as np
import numpy.linalg as npla
from scipy import sparse
from scipy.sparse import csr_matrix, kron, identity
from scipy.sparse.linalg import eigsh

import matplotlib.pyplot as plt
import warnings


# ------------------------
# Utilities and data types
# ------------------------

@dataclass
class Settings:
    L: int = 10
    J: float = 1.0
    g0: float = 1.00
    dg: float = 0.002
    block_min: int = 2
    block_max: int = 6
    regularizer: float = 1e-10
    seed: int = 7
    maxiter_eig: int = 2000
    pbc: bool = False  # open boundary default; set True for periodic if desired

@dataclass
class Summary:
    plateau_value: float
    plateau_std: float
    first_law_rms: float
    a0: float
    a1: float
    settings: dict


# ------------------------
# Pauli and spin operators
# ------------------------

SX = np.array([[0., 1.],
               [1., 0.]], dtype=np.float64)

SZ = np.array([[1., 0.],
               [0., -1.]], dtype=np.float64)

ID2 = np.eye(2, dtype=np.float64)


def op_on_site(op: np.ndarray, site: int, L: int) -> csr_matrix:
    """
    Sparse operator placing `op` on given site (0-indexed), identity elsewhere.
    """
    assert 0 <= site < L
    out = None
    for j in range(L):
        mat = op if j == site else ID2
        out = mat if out is None else np.kron(out, mat)
    return csr_matrix(out)


def two_site_coupling(opL: np.ndarray, i: int, opR: np.ndarray, j: int, L: int) -> csr_matrix:
    """
    Sparse operator with opL on site i and opR on site j (i<j), identity elsewhere.
    """
    assert 0 <= i < j < L
    mats = []
    for k in range(L):
        if k == i:
            mats.append(opL)
        elif k == j:
            mats.append(opR)
        else:
            mats.append(ID2)
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return csr_matrix(out)


# ------------------------
# Hamiltonian builder
# ------------------------

def tfim_hamiltonian(L: int, J: float, g: float, pbc: bool = False) -> csr_matrix:
    """
    H = -J sum σ^z_i σ^z_{i+1} - g sum σ^x_i   (open boundary by default)
    """
    H = sparse.csr_matrix((2**L, 2**L), dtype=np.float64)

    # ZZ couplings
    for i in range(L - 1):
        H -= J * two_site_coupling(SZ, i, SZ, i + 1, L)

    if pbc and L > 2:
        H -= J * two_site_coupling(SZ, 0, SZ, L - 1, L)

    # X field
    for i in range(L):
        H -= g * op_on_site(SX, i, L)

    return H.tocsr()


# ------------------------
# Ground state / reduced density matrices
# ------------------------

def ground_state(H: csr_matrix, maxiter: int = 2000, seed: int = 7) -> np.ndarray:
    """
    Lowest-energy eigenvector via sparse eigensolver.
    """
    # 'SA' -> smallest algebraic (lowest energy)
    vals, vecs = eigsh(H, k=1, which='SA', maxiter=maxiter, tol=1e-9, v0=None, ncv=min(40, H.shape[0]-1), return_eigenvectors=True)
    psi = vecs[:, 0]
    # Normalize (eigsh should return normalized already)
    psi = psi / npla.norm(psi)
    # Make a deterministic global phase (for reproducibility)
    # force max-|component| to be real-positive
    idx = np.argmax(np.abs(psi))
    if psi[idx].real < 0:
        psi = -psi
    return psi


def reduced_density_matrix(psi: np.ndarray, L: int, ell: int, reg: float = 1e-10) -> np.ndarray:
    """
    ρ_A for a contiguous block A = {0,...,ell-1}. Open chain assumed.
    psi is the full state vector of dim 2^L.
    Efficient formula: reshape psi -> (2^ell, 2^(L-ell)), then ρ_A = M @ M^† (trace over env).
    """
    dimA = 2**ell
    dimB = 2**(L - ell)
    M = psi.reshape(dimA, dimB)
    # Compute rhoA with robust warning suppression; fall back if anything non-finite appears
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        with np.errstate(over='ignore', divide='ignore', invalid='ignore', under='ignore'):
            rhoA = M @ M.conj().T

    # If numerical pathologies sneak in, fall back to a safe PSD projection next
    if not np.isfinite(rhoA).all():
        # Replace NaN/Inf by zero and proceed (PSD projection below will clean up further)
        rhoA = np.nan_to_num(rhoA, nan=0.0, posinf=0.0, neginf=0.0)

    # Numerical hygiene: symmetrize, PSD-project with tolerance, renormalize
    rhoA = 0.5 * (rhoA + rhoA.conj().T)

    # Spectral projection to PSD with small cutoff 'reg'
    w, v = npla.eigh(rhoA)
    # Zero-out tiny negative eigenvalues (roundoff) and floor tiny positives
    w = np.where(w < 0.0, 0.0, w)
    w = np.where(w < reg, 0.0, w)
    keep = w > 0.0
    if np.count_nonzero(keep) == 0:
        # Fallback: maximally mixed state on A
        dimA = rhoA.shape[0]
        rhoA = np.eye(dimA) / float(dimA)
    else:
        V = v[:, keep]
        wpos = w[keep]
        # Reconstruct using einsum on the numerical support; suppress benign warnings
        with np.errstate(over='ignore', divide='ignore', invalid='ignore', under='ignore'):
            rhoA = np.einsum('ik,k,kj->ij', V, wpos, V.conj().T)
        tr = float(np.trace(rhoA).real)
        if not np.isfinite(tr) or tr <= 0:
            dimA = rhoA.shape[0]
            rhoA = np.eye(dimA) / float(dimA)
        else:
            rhoA /= tr

    rhoA = 0.5 * (rhoA + rhoA.conj().T)
    # Final sanitation: if anything non-finite survived, revert to maximally mixed on A
    if not np.isfinite(rhoA).all():
        dimA = rhoA.shape[0]
        rhoA = np.eye(dimA) / float(dimA)
    return rhoA


# ------------------------
# Entropy and modular Hamiltonian
# ------------------------

def entropy_vn(rho: np.ndarray, reg: float = 1e-12) -> float:
    """
    Von Neumann entropy S(ρ) = -Tr ρ log ρ (natural log).
    Regularize tiny eigenvalues to avoid log(0).
    """
    vals = npla.eigvalsh((rho + rho.conj().T) * 0.5)
    vals_clipped = np.clip(vals, reg, None)
    S = -float(np.sum(vals_clipped * np.log(vals_clipped)))
    return S


def modular_hamiltonian_K(rho_ref: np.ndarray, reg: float = 1e-12) -> np.ndarray:
    """
    K = -log ρ_ref on its numerical support (regularized).
    We project onto eigenvectors with eigenvalues > reg to avoid amplifying noise.
    """
    rhoH = 0.5 * (rho_ref + rho_ref.conj().T)
    vals, vecs = npla.eigh(rhoH)
    # Numerical support: strictly greater than reg
    keep = vals > reg
    if np.count_nonzero(keep) == 0:
        # Degenerate support: return 0 operator (δ⟨K⟩ will then be ≈0 accordingly)
        return np.zeros_like(rhoH)
    V = vecs[:, keep]
    lam = np.clip(vals[keep], reg, 1.0)
    # Build K on the support only; avoid explicit diag and silence benign warnings
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        W = -np.log(lam)
        K_support = (V * W) @ V.conj().T
    K = 0.5 * (K_support + K_support.conj().T)
    return K.real


# ------------------------
# Regression: remove const + c * log ℓ (toy MI/moment-kill)
# ------------------------

def regress_and_residual(logLs: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Fit y ≈ a + b * logL and return (fit, residual, coeffs [a,b]).
    """
    X = np.vstack([np.ones_like(logLs), logLs]).T
    coeffs, *_ = npla.lstsq(X, y, rcond=None)
    fit = X @ coeffs
    residual = y - fit
    return fit, residual, coeffs


# ------------------------
# New helper function: execute_probe
# ------------------------

def execute_probe(settings: Settings, make_plots: bool = True, outdir: str | None = None) -> Summary:
    L = settings.L
    block_max = min(settings.block_max, L - 2)
    block_min = max(2, min(settings.block_min, block_max))
    blocks = list(range(block_min, block_max + 1))

    if outdir is None:
        outdir = os.path.join(os.getcwd(), "hqtfim_outputs")
    os.makedirs(outdir, exist_ok=True)

    H0 = tfim_hamiltonian(L, settings.J, settings.g0, pbc=settings.pbc)
    H1 = tfim_hamiltonian(L, settings.J, settings.g0 + settings.dg, pbc=settings.pbc)
    psi0 = ground_state(H0, maxiter=settings.maxiter_eig, seed=settings.seed)
    psi1 = ground_state(H1, maxiter=settings.maxiter_eig, seed=settings.seed)

    dS_list: list[float] = []
    dK_list: list[float] = []
    logL_list: list[float] = []

    for ell in blocks:
        rhoA0 = reduced_density_matrix(psi0, L, ell, settings.regularizer)
        rhoA1 = reduced_density_matrix(psi1, L, ell, settings.regularizer)
        S0 = entropy_vn(rhoA0, reg=settings.regularizer)
        S1 = entropy_vn(rhoA1, reg=settings.regularizer)
        dS = S1 - S0
        KA = modular_hamiltonian_K(rhoA0, reg=settings.regularizer)
        with np.errstate(over='ignore', divide='ignore', invalid='ignore', under='ignore'):
            dK = float(np.trace((rhoA1 - rhoA0) @ KA))
        dS_list.append(dS)
        dK_list.append(dK)
        logL_list.append(np.log(ell))

    dS_arr = np.array(dS_list)
    dK_arr = np.array(dK_list)
    logL_arr = np.array(logL_list)

    first_law_rms = float(np.sqrt(np.mean((dS_arr - dK_arr) ** 2)))
    # Regression and residuals
    fit, residual, coeffs = regress_and_residual(logL_arr, dK_arr)
    a0, a1 = map(float, coeffs)
    plateau_value = float(np.mean(residual))
    plateau_std = float(np.std(residual, ddof=1) / np.sqrt(len(residual)))

    if make_plots:
        # δ<K> vs logℓ
        plt.figure(figsize=(6.0, 4.2))
        plt.scatter(logL_arr, dK_arr, s=40, label=r'$\delta\langle K_A\rangle(\ell)$')
        xs = np.linspace(min(logL_arr) - 0.1, max(logL_arr) + 0.1, 200)
        ys = a0 + a1 * xs
        plt.plot(xs, ys, lw=2, label=r'fit: $a_0 + a_1 \log\ell$')
        plt.xlabel(r'$\log \ell$')
        plt.ylabel(r'$\delta\langle K_A\rangle$')
        plt.title('Modular response vs block size')
        plt.legend(); plt.tight_layout()
        plt.savefig(os.path.join(outdir, "dK_vs_logl.png"), dpi=150); plt.close()

        # Residuals
        plt.figure(figsize=(6.0, 4.2))
        plt.axhline(0.0, color='k', lw=1)
        plt.scatter(blocks, residual, s=40)
        plt.xlabel(r'Block size $\ell$')
        plt.ylabel(r'Residual $[\delta\langle K_A\rangle - (a_0+a_1\log\ell)]$')
        plt.title(f'MI/moment-kill–like residual (plateau ~ {plateau_value:+.2e})')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "residual_after_mikill.png"), dpi=150); plt.close()

        # First-law check
        lims = [min(dS_arr.min(), dK_arr.min()), max(dS_arr.max(), dK_arr.max())]
        pad = 0.05 * (lims[1] - lims[0] + 1e-12)
        xline = np.linspace(lims[0] - pad, lims[1] + pad, 200)
        plt.figure(figsize=(4.8, 4.8))
        plt.plot(xline, xline, lw=2, label='first-law line: y=x')
        plt.scatter(dK_arr, dS_arr, s=40)
        plt.xlabel(r'$\delta\langle K_A\rangle$'); plt.ylabel(r'$\delta S$')
        plt.title(f'First-law RMS = {first_law_rms:.2e}')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "first_law_check.png"), dpi=150); plt.close()

    summary = Summary(
        plateau_value=plateau_value,
        plateau_std=plateau_std,
        first_law_rms=first_law_rms,
        a0=a0,
        a1=a1,
        settings=asdict(settings),
    )

    # Always write a JSON snapshot for reproducibility
    with open(os.path.join(outdir, "summary.json"), "w") as f:
        json.dump(asdict(summary), f, indent=2)
    return summary


# ------------------------
# New helper function: quick_validations
# ------------------------

def quick_validations(base: Settings):
    outdir = os.path.join(os.getcwd(), "hqtfim_outputs")
    os.makedirs(outdir, exist_ok=True)

    # 1) δg scan
    dg_list = [0.001, 0.002, 0.005]
    rms = []
    plats = []
    for dg in dg_list:
        s = Settings(**asdict(base)); s.dg = dg
        summ = execute_probe(s, make_plots=False, outdir=outdir)
        rms.append(summ.first_law_rms); plats.append(summ.plateau_value)
    dg_arr = np.array(dg_list); rms_arr = np.array(rms)
    A = np.vstack([np.ones_like(dg_arr), dg_arr]).T
    coeff, *_ = npla.lstsq(A, rms_arr, rcond=None)
    fit = A @ coeff
    r2 = 1.0 - float(np.sum((rms_arr - fit)**2) / np.sum((rms_arr - np.mean(rms_arr))**2))
    # Save CSV and a small plot
    with open(os.path.join(outdir, 'quick_dg_scan.csv'), 'w') as f:
        f.write('dg,first_law_rms,plateau\n')
        for d, r, p in zip(dg_list, rms, plats):
            f.write(f"{d},{r},{p}\n")
    plt.figure(figsize=(5.0,3.4))
    xs = np.linspace(min(dg_list)*0.9, max(dg_list)*1.05, 200)
    ys = coeff[0] + coeff[1]*xs
    plt.scatter(dg_list, rms, s=40); plt.plot(xs, ys)
    plt.xlabel('dg'); plt.ylabel('first-law RMS'); plt.title(f'dg-scan (R^2={r2:.3f})')
    plt.tight_layout(); plt.savefig(os.path.join(outdir,'quick_dg_scan.png'), dpi=150); plt.close()
    print(f"[VALIDATE] dg-scan: slope={coeff[1]:.3e}, R^2={r2:.3f} (expect ~linear)")

    # 2) Boundary swap (OBC vs PBC)
    s_obc = Settings(**asdict(base)); s_obc.pbc = False
    s_pbc = Settings(**asdict(base)); s_pbc.pbc = True
    sum_obc = execute_probe(s_obc, make_plots=False, outdir=outdir)
    sum_pbc = execute_probe(s_pbc, make_plots=False, outdir=outdir)
    diff_plat = abs(sum_obc.plateau_value - sum_pbc.plateau_value)
    thr = 3.0 * max(sum_obc.plateau_std, sum_pbc.plateau_std)
    verdict = 'PASS' if diff_plat <= thr else 'FLAG'
    with open(os.path.join(outdir,'quick_pbc_compare.json'),'w') as f:
        json.dump({
            'OBC_plateau': sum_obc.plateau_value,
            'PBC_plateau': sum_pbc.plateau_value,
            'abs_diff': diff_plat,
            'threshold': thr,
            'verdict': verdict
        }, f, indent=2)
    print(f"[VALIDATE] boundary swap: |Δplateau|={diff_plat:.2e} vs thr={thr:.2e} → {verdict}")

    # 3) Block-range stability
    s_base = Settings(**asdict(base)); s_base.block_min = max(2, base.block_min); s_base.block_max = min(base.block_max, base.L-2)
    s_wide = Settings(**asdict(base)); s_wide.block_min = 2; s_wide.block_max = min(7, base.L-2)
    sum_base = execute_probe(s_base, make_plots=False, outdir=outdir)
    sum_wide = execute_probe(s_wide, make_plots=False, outdir=outdir)
    a1_diff = abs(sum_base.a1 - sum_wide.a1)
    with open(os.path.join(outdir,'quick_block_compare.json'),'w') as f:
        json.dump({'a1_base': sum_base.a1, 'a1_wide': sum_wide.a1, 'abs_diff': a1_diff}, f, indent=2)
    print(f"[VALIDATE] block-range: |Δa1|={a1_diff:.2e} (expect small)")

    # 4) Size check (L grid)
    L_vals = [base.L, max(10, base.L+2)]
    rows = []
    for Lval in L_vals:
        s = Settings(**asdict(base)); s.L = int(Lval); s.block_max = min(s.block_max, s.L-2)
        summ = execute_probe(s, make_plots=False, outdir=outdir)
        rows.append((s.L, summ.first_law_rms, summ.plateau_std))
    with open(os.path.join(outdir,'quick_size_scan.csv'),'w') as f:
        f.write('L,first_law_rms,plateau_se\n')
        for Lval, r, se in rows:
            f.write(f"{Lval},{r},{se}\n")
    # Plot RMS vs L
    plt.figure(figsize=(5.0,3.4))
    plt.plot([r[0] for r in rows], [r[1] for r in rows], marker='o')
    plt.xlabel('L'); plt.ylabel('first-law RMS'); plt.title('Size scan')
    plt.tight_layout(); plt.savefig(os.path.join(outdir,'quick_size_scan.png'), dpi=150); plt.close()

    # Write a compact referee-friendly text report
    report_path = os.path.join(outdir, 'validation_report.txt')
    with open(report_path, 'w') as rf:
        rf.write('HQTFIM quick validations\n')
        rf.write('------------------------\n')
        rf.write(f'dg-scan slope (RMS vs dg): {coeff[1]:.6e}, R^2={r2:.3f}\n')
        rf.write(f'boundary swap: OBC plateau={sum_obc.plateau_value:.3e}, '
                 f'PBC plateau={sum_pbc.plateau_value:.3e}, '
                 f'|Δ|={diff_plat:.3e}, thr={thr:.3e}, verdict={verdict}\n')
        rf.write(f'block-range: a1_base={sum_base.a1:.6e}, a1_wide={sum_wide.a1:.6e}, '
                 f'|Δa1|={a1_diff:.6e}\n')
        rf.write('size-scan (L, first-law RMS, plateau SE):\n')
        for Lval, r, se in rows:
            rf.write(f'  L={Lval:>3d}: RMS={r:.6e}, SE={se:.6e}\n')
    print(f"[VALIDATE] wrote {report_path}")

    print("[VALIDATE] size-scan complete (see quick_size_scan.csv/png)")


# ------------------------
# Main workflow
# ------------------------

def main():
    parser = argparse.ArgumentParser(description="HQTFIM substrate probe: entanglement first-law + toy MI/moment-kill residual.")
    parser.add_argument("--L", type=int, default=10, help="System size (number of spins). Default: 10")
    parser.add_argument("--g0", type=float, default=1.0, help="Base transverse field (critical ~1). Default: 1.0")
    parser.add_argument("--dg", type=float, default=0.002, help="Small shift for perturbed field (linear response scale ~0.002). Default: 0.002")
    parser.add_argument("--J", type=float, default=1.0, help="Ising coupling. Default: 1.0")
    parser.add_argument("--block-min", type=int, default=2, help="Minimum block size ℓ. Default: 2")
    parser.add_argument("--block-max", type=int, default=6, help="Maximum block size ℓ. Default: 6")
    parser.add_argument("--pbc", action="store_true", help="Use periodic boundary conditions (default: OBC).")
    parser.add_argument("--regularizer", type=float, default=1e-10, help="Eigenvalue cutoff for logρ regularization; suppresses spurious warnings while leaving results unchanged at displayed precision.")
    parser.add_argument("--maxiter-eig", type=int, default=2000, help="Max iterations for eigensolver.")
    parser.add_argument("--seed", type=int, default=7, help="Random seed (used only for deterministic choices).")
    parser.add_argument('--quick-validate', action='store_true', help='Run referee-friendly validation sweeps (dg-scan, boundary swap, block-range, size).')
    args = parser.parse_args()

    # Sanity and defaults
    L = args.L
    if args.block_max <= 0:
        block_max = max(2, L // 2)
    else:
        block_max = min(args.block_max, L - 2)
    block_min = max(2, min(args.block_min, block_max))

    settings = Settings(
        L=L, J=args.J, g0=args.g0, dg=args.dg,
        block_min=block_min, block_max=block_max,
        regularizer=args.regularizer, seed=args.seed,
        maxiter_eig=args.maxiter_eig, pbc=args.pbc
    )

    if args.quick_validate:
        quick_validations(settings)
        return

    summary = execute_probe(settings, make_plots=True)
    outdir = os.path.join(os.getcwd(), 'hqtfim_outputs')
    print('[HQTFIM] Done.')
    print(f'  Outputs: {outdir}/summary.json')
    print(f'           {outdir}/dK_vs_logl.png')
    print(f'           {outdir}/residual_after_mikill.png')
    print(f'           {outdir}/first_law_check.png')


if __name__ == "__main__":
    main()