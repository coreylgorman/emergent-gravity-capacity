#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gaussian_capacity_probe.py

Referee-ready Gaussian substrate probe for the "capacity-throttling" modular response:
- 1D free fermion chain (tight-binding)
- 1D harmonic (scalar) chain

For each model, we compute δS and δ⟨K⟩ across a range of subregion sizes ℓ,
fit δ⟨K⟩ ≈ a0 + a1 log ℓ, subtract the [1, log ℓ] trend, and report the
residual "plateau" (mean ± SE). We also report the first-law RMS: RMS(δS - δ⟨K⟩).
For fermions, δ⟨K⟩ is computed exactly from the quadratic modular kernel.
For scalars, δ⟨K⟩ is taken to equal δS to first order (first-law linear response),
and this is clearly labeled in outputs.

Quick validations (--quick-validate) run small, fast checks that referees asked for:
- δg (deformation) scan → first-law RMS ~ linear in deformation size
- boundary swap (PBC ↔ OBC) → plateau invariance (fermion model)
- block-range stability → small change of slope a1
- size scan → mild shrink of RMS/SE with L

Outputs are written to ./gaussian_outputs/.

Usage examples:
  # Free fermion, L=200, blocks 20..100, chemical-potential shift δμ=0.01
  python gaussian_capacity_probe.py --model fermion1d --L 200 --ell-min 20 --ell-max 100 --delta-m 0.01

  # Harmonic chain (scalar), L=200, mass shift δm=0.01
  python gaussian_capacity_probe.py --model scalar1d --L 200 --ell-min 20 --ell-max 100 --delta-m 0.01

  # Referee pack (quick validations) on fermion model
  python gaussian_capacity_probe.py --model fermion1d --quick-validate

Author: (c) 2025
License: MIT
"""

from __future__ import annotations
import argparse
import json
import math
import os
from dataclasses import dataclass
from typing import Tuple, List, Dict

import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt


# ------------------------------ Utilities ------------------------------ #

EPS_EIG = 1e-12  # clipping for eigenvalues near 0 or 1
OUTDIR = "./gaussian_outputs"


def ensure_outdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def se_from_samples(x: np.ndarray) -> float:
    if x.size < 2:
        return 0.0
    return float(np.std(x, ddof=1) / math.sqrt(x.size))


def linfit_logx(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """
    Fit y ≈ a0 + a1 * log(x) using least squares.
    Returns (a0, a1).
    """
    X = np.column_stack([np.ones_like(x), np.log(x)])
    coef, _, _, _ = npl.lstsq(X, y, rcond=None)
    return float(coef[0]), float(coef[1])


# ------------------------------ Fermion 1D ------------------------------ #

def build_fermion_H(L: int, t: float = 1.0, mu: float = 0.0, pbc: bool = True) -> np.ndarray:
    """
    Tight-binding single-particle Hamiltonian:
      H = -t ∑ (|i><i+1| + h.c.) - mu ∑ |i><i|
    pbc: include wrap-around hopping.
    """
    H = np.zeros((L, L), dtype=float)
    for i in range(L - 1):
        H[i, i + 1] = H[i + 1, i] = -t
    if pbc:
        H[0, L - 1] = H[L - 1, 0] = -t
    H -= mu * np.eye(L)
    return H


def fermion_groundstate_C(H: np.ndarray) -> np.ndarray:
    """
    Build the real-space correlation matrix C = <c_i^\dagger c_j>
    by diagonalizing H and filling all negative-energy single-particle modes.
    """
    evals, evecs = npl.eigh(H)
    occ = evals < 0.0
    # C = U occ U^\dagger
    C = (evecs[:, occ] @ evecs[:, occ].T.conj()).real
    return C


def entropy_fermion(CA: np.ndarray) -> float:
    """Entanglement entropy S(C) = -Tr[C log C + (I-C) log(I-C)] for fermionic Gaussian state."""
    lam, _ = npl.eigh((CA + CA.T) * 0.5)  # ensure Hermitian numerical symmetrization
    lam = np.clip(lam, EPS_EIG, 1.0 - EPS_EIG)
    S = -np.sum(lam * np.log(lam) + (1.0 - lam) * np.log(1.0 - lam))
    return float(S)


def modular_kernel_fermion(CA0: np.ndarray) -> np.ndarray:
    """
    For fermions, the quadratic modular generator (one-body kernel) is:
      h0 = log[(I - C0)/C0]  in the one-body subspace of region A.
    We compute via eigen-decomposition of C0: C0 = V diag(λ) V†, then
      h0 = V diag(log((1-λ)/λ)) V†.
    """
    lam, V = npl.eigh((CA0 + CA0.T) * 0.5)
    lam = np.clip(lam, EPS_EIG, 1.0 - EPS_EIG)
    diag = np.log((1.0 - lam) / lam)
    h0 = (V @ np.diag(diag) @ V.T.conj()).real
    return h0


def run_fermion_chain(L: int,
                      ell_list: List[int],
                      delta_mu: float,
                      pbc: bool = True) -> Dict[str, np.ndarray]:
    """
    Compute δS, δ⟨K⟩ across ℓ for the fermion chain with chemical potential shift μ→μ+δμ.
    Baseline μ = 0 (half-filled symmetric TB chain).
    """
    H0 = build_fermion_H(L, t=1.0, mu=0.0, pbc=pbc)
    H1 = build_fermion_H(L, t=1.0, mu=delta_mu, pbc=pbc)

    C0 = fermion_groundstate_C(H0)
    C1 = fermion_groundstate_C(H1)

    dS_list, dK_list = [], []
    for ell in ell_list:
        A = slice(0, ell)
        C0A = C0[A, A]
        C1A = C1[A, A]

        S0 = entropy_fermion(C0A)
        S1 = entropy_fermion(C1A)
        dS = S1 - S0

        h0 = modular_kernel_fermion(C0A)
        dK = float(np.trace((C1A - C0A) @ h0))  # exact quadratic formula

        dS_list.append(dS)
        dK_list.append(dK)

    return {
        "ell": np.array(ell_list, dtype=float),
        "dS": np.array(dS_list, dtype=float),
        "dK": np.array(dK_list, dtype=float),
    }


# ------------------------------ Scalar 1D (harmonic chain) ------------------------------ #

def scalar_1d_correlators(L: int, m: float, kappa: float = 1.0, pbc: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build ground-state correlators for a harmonic chain (quadratic boson):
      H = 1/2 ∑ [ p_i^2 + m^2 x_i^2 + κ (x_i - x_{i+1})^2 ]
    Return (X, P) with X_{ij} = <x_i x_j>, P_{ij} = <p_i p_j>.
    Uses Fourier sum (PBC) or approximate cosine modes (OBC fallback).
    """
    X = np.zeros((L, L), dtype=float)
    P = np.zeros((L, L), dtype=float)

    if pbc:
        ks = 2.0 * np.pi * np.arange(L) / L
        omega = np.sqrt(m * m + 2.0 * kappa * (1.0 - np.cos(ks)))
        # avoid exact zero mode in massless limit
        omega = np.clip(omega, 1e-9, None)
        phase = np.exp(1j * ks)
        # Toeplitz structure: corr depends on (i-j) mod L
        for r in range(L):
            # X(r) = (1/2L) ∑_k (1/ω_k) e^{ik r}
            Xr = 0.5 * np.sum(np.cos(ks * r) / omega) / L
            # P(r) = (1/2L) ∑_k ω_k e^{ik r}
            Pr = 0.5 * np.sum(np.cos(ks * r) * omega) / L
            for i in range(L):
                j = (i + r) % L
                X[i, j] = Xr
                P[i, j] = Pr
    else:
        # OBC: use sine/cosine modes; good enough for small L in validations
        ns = np.arange(1, L + 1)
        ks = np.pi * ns / (L + 1)
        omega = np.sqrt(m * m + 2.0 * kappa * (1.0 - np.cos(ks)))
        omega = np.clip(omega, 1e-9, None)
        phi = np.sqrt(2.0 / (L + 1))  # normalization
        for i in range(L):
            for j in range(L):
                cosi = np.sin(ks * (i + 1))
                cosj = np.sin(ks * (j + 1))
                X[i, j] = 0.5 * np.sum((phi * phi) * cosi * cosj / omega)
                P[i, j] = 0.5 * np.sum((phi * phi) * cosi * cosj * omega)

    return X, P


def entropy_scalar_from_XP(XA: np.ndarray, PA: np.ndarray) -> float:
    """
    Entanglement entropy of a bosonic Gaussian state with no <xp> cross terms.
    Symplectic eigenvalues ν_j = sqrt(eig(XA @ PA)), with ν_j >= 1/2.
    S = ∑ [ (ν+1/2) ln(ν+1/2) - (ν-1/2) ln(ν-1/2) ].
    """
    # Numerical symmetrization
    XA = (XA + XA.T) * 0.5
    PA = (PA + PA.T) * 0.5
    M = XA @ PA
    # Guard against tiny negative due to rounding
    w, _ = npl.eigh((M + M.T) * 0.5)
    w = np.clip(w, 1e-12, None)
    nu = np.sqrt(w)
    # Ensure minimal uncertainty
    nu = np.clip(nu, 0.5 + 1e-12, None)

    S = np.sum((nu + 0.5) * np.log(nu + 0.5) - (nu - 0.5) * np.log(nu - 0.5))
    return float(S)


def run_scalar_chain(L: int,
                     ell_list: List[int],
                     m0: float,
                     delta_m: float,
                     pbc: bool = True) -> Dict[str, np.ndarray]:
    """
    Compute δS and (by first law) δ⟨K⟩ ≈ δS for harmonic chain under m→m+δm.
    We label δ⟨K⟩ as 'dK_firstlaw' in the JSON for clarity.
    """
    X0, P0 = scalar_1d_correlators(L, m0, kappa=1.0, pbc=pbc)
    X1, P1 = scalar_1d_correlators(L, m0 + delta_m, kappa=1.0, pbc=pbc)

    dS_list, dK_list = [], []
    for ell in ell_list:
        A = slice(0, ell)
        S0 = entropy_scalar_from_XP(X0[A, A], P0[A, A])
        S1 = entropy_scalar_from_XP(X1[A, A], P1[A, A])
        dS = S1 - S0
        dS_list.append(dS)
        dK_list.append(dS)  # first-law linear response for small deformation

    return {
        "ell": np.array(ell_list, dtype=float),
        "dS": np.array(dS_list, dtype=float),
        "dK": np.array(dK_list, dtype=float),  # equals dS here by construction
    }


# ------------------------------ Plotting ------------------------------ #

def plot_dK_vs_logl(ell: np.ndarray, dK: np.ndarray, a0: float, a1: float, outpng: str) -> None:
    plt.figure(figsize=(5.0, 4.0), dpi=150)
    x = np.log(ell)
    plt.scatter(x, dK, s=18, label=r"$\delta\langle K\rangle$")
    xx = np.linspace(x.min(), x.max(), 200)
    yy = a0 + a1 * xx
    plt.plot(xx, yy, lw=2, label=rf"fit: $a_0 + a_1 \log \ell$, $a_1={a1:.3e}$")
    plt.xlabel(r"$\log \ell$")
    plt.ylabel(r"$\delta\langle K\rangle$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()


def plot_residual_plateau(ell: np.ndarray, residual: np.ndarray, outpng: str) -> None:
    plt.figure(figsize=(5.0, 4.0), dpi=150)
    plt.scatter(ell, residual, s=18)
    plt.axhline(0.0, color="k", lw=1, alpha=0.5)
    plt.xlabel(r"block size $\ell$")
    plt.ylabel(r"residual after $[1,\log\ell]$ subtraction")
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()


def plot_firstlaw_check(dS: np.ndarray, dK: np.ndarray, rms: float, outpng: str) -> None:
    plt.figure(figsize=(5.0, 4.0), dpi=150)
    plt.scatter(dK, dS, s=20, label=r"$(\delta\langle K\rangle,\ \delta S)$")
    lo = min(dK.min(), dS.min())
    hi = max(dK.max(), dS.max())
    plt.plot([lo, hi], [lo, hi], "k--", lw=1, label="y=x")
    plt.xlabel(r"$\delta\langle K\rangle$")
    plt.ylabel(r"$\delta S$")
    plt.title(rf"First-law RMS: {rms:.3e}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()


# ------------------------------ Runner & Validations ------------------------------ #

@dataclass
class RunConfig:
    model: str
    L: int
    ell_min: int
    ell_max: int
    ell_step: int
    pbc: bool
    delta_m: float
    m0_scalar: float
    outdir: str


def run_core(cfg: RunConfig) -> Dict[str, float]:
    ensure_outdir(cfg.outdir)

    # Assemble block sizes
    ell_list = list(range(cfg.ell_min, cfg.ell_max + 1, cfg.ell_step))
    if len(ell_list) < 3:
        raise ValueError("Need at least 3 block sizes for a stable fit.")

    # Compute dS, dK
    if cfg.model == "fermion1d":
        out = run_fermion_chain(cfg.L, ell_list, delta_mu=cfg.delta_m, pbc=cfg.pbc)
        model_note = "fermion (exact modular kernel)"
    elif cfg.model == "scalar1d":
        out = run_scalar_chain(cfg.L, ell_list, m0=cfg.m0_scalar, delta_m=cfg.delta_m, pbc=cfg.pbc)
        model_note = "scalar (first-law δK=δS)"
    else:
        raise ValueError("Unsupported model. Use --model fermion1d or scalar1d.")

    ell = out["ell"]
    dS = out["dS"]
    dK = out["dK"]

    # Fit dK ≈ a0 + a1 log ℓ
    a0, a1 = linfit_logx(ell, dK)
    residual = dK - (a0 + a1 * np.log(ell))
    plateau_mean = float(np.mean(residual))
    plateau_se = se_from_samples(residual)

    # First-law RMS
    fl_rms = float(np.sqrt(np.mean((dS - dK) ** 2)))

    # Plots
    plot_dK_vs_logl(ell, dK, a0, a1, os.path.join(cfg.outdir, "dK_vs_logl.png"))
    plot_residual_plateau(ell, residual, os.path.join(cfg.outdir, "residual_after_subtraction.png"))
    plot_firstlaw_check(dS, dK, fl_rms, os.path.join(cfg.outdir, "first_law_check.png"))

    # Summary JSON
    summary = {
        "model": cfg.model,
        "model_note": model_note,
        "L": cfg.L,
        "pbc": cfg.pbc,
        "ell_min": cfg.ell_min,
        "ell_max": cfg.ell_max,
        "ell_step": cfg.ell_step,
        "delta_m": cfg.delta_m,
        "m0_scalar": cfg.m0_scalar if cfg.model == "scalar1d" else None,
        "fit": {"a0": a0, "a1": a1},
        "plateau": {"mean": plateau_mean, "se": plateau_se},
        "first_law_rms": fl_rms,
        "notes": "For scalar1d, δK equals δS by first-law (linear response). Fermion1d uses exact modular kernel.",
    }
    with open(os.path.join(cfg.outdir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    # Console print
    print(f"[GAUSS] {cfg.model} | L={cfg.L} | PBC={cfg.pbc}")
    print(f"  Plateau mean ± SE: {plateau_mean:+.3e} ± {plateau_se:.3e}")
    print(f"  First-law RMS(δS-δK): {fl_rms:.3e}")
    print(f"  Fit a1 (slope vs log ℓ): {a1:+.3e}")
    print(f"  Outputs: {cfg.outdir}/summary.json")
    print(f"           {cfg.outdir}/dK_vs_logl.png")
    print(f"           {cfg.outdir}/residual_after_subtraction.png")
    print(f"           {cfg.outdir}/first_law_check.png")

    return summary


def quick_validate(cfg: RunConfig) -> None:
    """
    Small, fast validation pack for referees:
      1) δg-scan (delta_m): measure slope of first-law RMS vs δg
      2) boundary swap (PBC↔OBC): plateau invariance (fermion only)
      3) block-range stability: change in fitted slope a1
      4) size scan: mild shrink of RMS/SE with L
    """
    ensure_outdir(cfg.outdir)

    # 1) δg-scan
    deltas = np.array([cfg.delta_m * x for x in [0.5, 1.0, 2.0]])
    rms_list = []
    for d in deltas:
        subcfg = RunConfig(cfg.model, cfg.L, cfg.ell_min, cfg.ell_max, cfg.ell_step, cfg.pbc, d, cfg.m0_scalar, cfg.outdir)
        res = run_core(subcfg)
        rms_list.append(res["first_law_rms"])
    rms_list = np.array(rms_list, dtype=float)
    A = np.column_stack([np.ones_like(deltas), deltas])
    coef, _, _, _ = npl.lstsq(A, rms_list, rcond=None)
    slope = float(coef[1])
    # crude R^2
    yhat = A @ coef
    SSres = float(np.sum((rms_list - yhat) ** 2))
    SStot = float(np.sum((rms_list - np.mean(rms_list)) ** 2))
    R2 = 1.0 - SSres / (SStot + 1e-30)

    np.savetxt(os.path.join(cfg.outdir, "quick_dg_scan.csv"),
               np.column_stack([deltas, rms_list]),
               delimiter=",", header="delta,rms_firstlaw", comments="")
    plt.figure(figsize=(4.5, 3.6), dpi=150)
    plt.plot(deltas, rms_list, "o-", label="RMS vs δ")
    plt.xlabel("δ (deformation size)")
    plt.ylabel("First-law RMS")
    plt.title(f"slope={slope:.3e}, R^2={R2:.3f}")
    plt.tight_layout()
    plt.savefig(os.path.join(cfg.outdir, "quick_dg_scan.png"))
    plt.close()

    # 2) boundary swap (fermion only)
    boundary_delta = None
    if cfg.model == "fermion1d":
        subP = RunConfig(cfg.model, cfg.L, cfg.ell_min, cfg.ell_max, cfg.ell_step, True, cfg.delta_m, cfg.m0_scalar, cfg.outdir)
        subO = RunConfig(cfg.model, cfg.L, cfg.ell_min, cfg.ell_max, cfg.ell_step, False, cfg.delta_m, cfg.m0_scalar, cfg.outdir)
        resP = run_core(subP)
        resO = run_core(subO)
        boundary_delta = abs(resP["plateau"]["mean"] - resO["plateau"]["mean"])
        with open(os.path.join(cfg.outdir, "quick_pbc_compare.json"), "w") as f:
            json.dump({"pbc_plateau": resP["plateau"], "obc_plateau": resO["plateau"]}, f, indent=2)

    # 3) block-range stability
    # widen (or shift) block window slightly
    ell_min2 = max(2, cfg.ell_min - max(1, cfg.ell_step))
    ell_max2 = cfg.ell_max + max(1, cfg.ell_step)
    subR = RunConfig(cfg.model, cfg.L, ell_min2, ell_max2, cfg.ell_step, cfg.pbc, cfg.delta_m, cfg.m0_scalar, cfg.outdir)
    res_base = run_core(cfg)
    res_wide = run_core(subR)
    delta_a1 = abs(res_base["fit"]["a1"] - res_wide["fit"]["a1"])
    with open(os.path.join(cfg.outdir, "quick_block_compare.json"), "w") as f:
        json.dump({"a1_base": res_base["fit"]["a1"], "a1_wide": res_wide["fit"]["a1"], "delta_a1": delta_a1}, f, indent=2)

    # 4) size scan
    Ls = [max(60, cfg.L - 40), cfg.L, cfg.L + 40]
    rows = []
    for Lval in Ls:
        subS = RunConfig(cfg.model, Lval, cfg.ell_min, cfg.ell_max, cfg.ell_step, cfg.pbc, cfg.delta_m, cfg.m0_scalar, cfg.outdir)
        resS = run_core(subS)
        rows.append([Lval, resS["first_law_rms"], resS["plateau"]["se"]])
    rows = np.array(rows, dtype=float)
    np.savetxt(os.path.join(cfg.outdir, "quick_size_scan.csv"),
               rows, delimiter=",", header="L,firstlaw_rms,plateau_se", comments="")
    plt.figure(figsize=(4.5, 3.6), dpi=150)
    plt.plot(rows[:, 0], rows[:, 1], "o-", label="RMS(δS-δK)")
    plt.plot(rows[:, 0], rows[:, 2], "s--", label="plateau SE")
    plt.xlabel("L")
    plt.ylabel("magnitude")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(cfg.outdir, "quick_size_scan.png"))
    plt.close()

    # Report
    rep = []
    rep.append(f"[VALIDATE] dg-scan: slope={slope:.3e}, R^2={R2:.3f} (expect ~linear)")
    if boundary_delta is not None:
        thr = max(1e-4, 5.0 * res_base["plateau"]["se"])
        verdict = "PASS" if boundary_delta < thr else "DRIFT"
        rep.append(f"[VALIDATE] boundary swap: |Δplateau|={boundary_delta:.2e} vs thr={thr:.2e} → {verdict}")
    rep.append(f"[VALIDATE] block-range: |Δa1|={delta_a1:.2e} (expect small)")
    rep.append("[VALIDATE] size-scan complete (see quick_size_scan.csv/png)")
    with open(os.path.join(cfg.outdir, "validation_report.txt"), "w") as f:
        f.write("\n".join(rep))
    for line in rep:
        print(line)


# ------------------------------ CLI ------------------------------ #

def parse_args() -> RunConfig:
    p = argparse.ArgumentParser(description="Gaussian substrate probe for modular-response plateau.")
    p.add_argument("--model", type=str, default="fermion1d",
                   choices=["fermion1d", "scalar1d"],
                   help="Choose substrate model.")
    p.add_argument("--L", type=int, default=200, help="System size (chain length).")
    p.add_argument("--ell-min", type=int, default=20, help="Minimum block size ℓ.")
    p.add_argument("--ell-max", type=int, default=100, help="Maximum block size ℓ.")
    p.add_argument("--ell-step", type=int, default=10, help="Step for block size ℓ.")
    p.add_argument("--pbc", action="store_true", help="Use periodic boundary conditions (default True).")
    p.add_argument("--obc", action="store_true", help="Use open boundary conditions (overrides --pbc).")
    p.add_argument("--delta-m", type=float, default=0.01,
                   help="Small deformation: fermion uses δμ (chem. potential shift); scalar uses δm (mass shift).")
    p.add_argument("--m0-scalar", type=float, default=0.10,
                   help="Baseline mass m0 for harmonic chain (avoid IR divergence).")
    p.add_argument("--outdir", type=str, default=OUTDIR, help="Output directory.")
    p.add_argument("--quick-validate", action="store_true", help="Run referee-friendly quick validations.")
    args = p.parse_args()

    # Boundary condition resolve
    pbc = True
    if args.obc:
        pbc = False
    elif args.pbc:
        pbc = True

    cfg = RunConfig(
        model=args.model,
        L=args.L,
        ell_min=args.ell_min,
        ell_max=args.ell_max,
        ell_step=args.ell_step,
        pbc=pbc,
        delta_m=args.delta_m,
        m0_scalar=args.m0_scalar,
        outdir=args.outdir,
    )
    return cfg


def main():
    cfg = parse_args()
    summary = run_core(cfg)
    if cfg.quick_validate:
        quick_validate(cfg)


if __name__ == "__main__":
    main()