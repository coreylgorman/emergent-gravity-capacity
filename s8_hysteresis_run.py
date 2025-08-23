#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
s8_hysteresis_run.py
First-principles S8 / growth analysis with entropic state-action (ΔS≥0).

Outputs:
  - s8_state_action_summary.json
  - fs8_comparison.png, E_of_z_check.png, gw_em_ratio.png
  - s8_residuals.png, s8_bestfit_lines.png
  - s8_p_sweep.png (kernel power sweep; optional fσ8 overlay from CSV if present)

No CLI needed; defaults are conservative and referee-ready.
"""

import json, math, os, sys
import numpy as np
import matplotlib.pyplot as plt

# Model label and summary JSON basename (referee-friendly)
MODEL_LABEL = "Entropic State-Action (ΔS≥0)"
SUMMARY_JSON_BASENAME = "s8_state_action_summary.json"

def ensure_rsds_csv(path="fs8_compilation.csv"):
    """
    If an observational RSD compilation is not present, write a small,
    widely used fallback set (6dFGS, SDSS-MGS, BOSS LOWZ/CMASS, WiggleZ, VIPERS).
    Replace with your curated file for paper use.
    """
    if os.path.exists(path):
        return path
    rows = [
        # label,   z,     fσ8,    err
        ("6dFGS",   0.067, 0.423, 0.055),
        ("SDSS-MGS",0.150, 0.490, 0.140),
        ("BOSS-LOWZ",0.320,0.384, 0.095),
        ("BOSS-CMASS",0.570,0.441, 0.043),
        ("WiggleZ", 0.440, 0.413, 0.080),
        ("WiggleZ", 0.600, 0.390, 0.063),
        ("WiggleZ", 0.730, 0.437, 0.072),
        ("VIPERS",  0.800, 0.470, 0.080),
    ]
    with open(path, "w") as f:
        f.write("z,fs8,err,label\n")
        for label, z, val, err in rows:
            f.write(f"{z:.3f},{val:.3f},{err:.3f},{label}\n")
    print(f"[info] Wrote a small fallback RSD compilation to '{path}'. Replace with your curated file when ready.")
    return path

# ------------------------
# Background cosmology (flat LCDM for H(a))
# ------------------------
Om0 = 0.315    # Planck-like reference (can read from your paper's constants if desired)
Ol0 = 1.0-Om0
H0  = 67.4     # km/s/Mpc Planck
h   = H0/100.0

def E_of_a(a):
    """ E(a) = H(a)/H0 """
    return np.sqrt(Om0*a**-3 + Ol0)

def dlnH_dla(a):
    """ d ln H / d ln a """
    Om_a = Om0*a**-3 / E_of_a(a)**2
    return -0.5*3*Om_a  # since H^2 ∝ Om0 a^-3 + Ol0

def Om_of_a(a):
    return Om0*a**-3 / E_of_a(a)**2

# ------------------------
# Linear growth ODE with Jordan-frame α_M(a) friction
# D'' + [ 2 + dlnH/dlna + α_M(a) ] D' + (3/2) μ(a) Ω_m(a) D = 0
# (Sign fix: the source term is +3/2 μ Ω_m D in d/d ln a form.)
# ------------------------
def solve_growth(a_grid, alphaM_of_a, mu_of_a=None, normalize=False):
    D   = np.zeros_like(a_grid)
    dD  = np.zeros_like(a_grid)
    # ICs in matter era: D ~ a, set at small a_ini
    D[0]  = a_grid[0]
    dD[0] = 1.0
    for i in range(len(a_grid)-1):
        a  = a_grid[i]
        ap = a_grid[i+1]
        la = math.log(a); lap = math.log(ap)
        dln_a = lap - la
        A = 2.0 + dlnH_dla(a) + alphaM_of_a(a)
        mu = 1.0 if mu_of_a is None else mu_of_a(a)
        # Correct source-term sign: D'' = -A D' + (3/2) μ Ω_m D
        B =  +1.5 * mu * Om_of_a(a)
        # 1st-order step for D' (Euler); small steps => stable enough for our grids
        dD[i+1] = dD[i] + (-A*dD[i] + B*D[i]) * dln_a
        D[i+1]  = D[i]  + dD[i] * dln_a
    if normalize:
        D /= (D[-1] if D[-1] != 0 else 1.0)
    return D

# ------------------------
# Entropic state-action (hysteresis)
# ε(a) = ε0 + c_log * ln[ 1 + J(a)/J_* ]
# with J(a) = ∫^{a} K(a,a') 𝔈(a') d ln a',  K ∝ (a'/a)^p,  𝔈 ∝ D(a')^2
# ------------------------
def cumulative_J(a_grid, D, p=5):
    # compute J(a_j) for each j with trapezoid in ln a'
    ln_a  = np.log(a_grid)
    D2    = (D**2)
    J = np.zeros_like(a_grid)
    for j in range(1, len(a_grid)):
        aj   = a_grid[j]
        K    = (a_grid[:j]/aj)**p
        integrand = K * D2[:j]
        J[j] = np.trapezoid(integrand, ln_a[:j])
    # Enforce monotonic non-decrease due to numerical noise
    J = np.maximum.accumulate(J)
    return J

def calibrate_normalization(J, a_grid, omega_from_pipeline, eps0=0.0):
    """
    Fix c_log and J_* from zero-mode Clausius balance INCLUDING irreversibility floor:
      - Naturalize J_* to the present cumulative exposure: J_* = J(a=1)
      - Enforce ∫ ε(a) d ln a = Ω_Λ  (unit-solid-angle Noether balance)
    This yields c_log uniquely; ε0 may be set by irreversibility floor (small, >0).
    """
    Jstar = J[-1] if J[-1] > 0 else 1.0
    ln_a  = np.log(a_grid)
    # raw integral with c_log=1 (variable part only)
    eps_raw = np.log(1.0 + J/Jstar)
    I_raw   = np.trapezoid(eps_raw, ln_a)
    # Total Clausius balance: ∫[eps0 + c_log * ln(1+J/J*)] d ln a = Ω_Λ
    ln_range = ln_a[-1] - ln_a[0]  # = ln(a_end/a_start)
    rhs = omega_from_pipeline - eps0*ln_range
    if I_raw <= 0 or rhs <= 0:
        c_log = 0.0
    else:
        c_log = rhs / I_raw
    return c_log, Jstar

def epsilon_of_a(J, a_grid, c_log, Jstar, eps0=0.0):
    return eps0 + c_log*np.log(1.0 + J/Jstar)

# ------------------------
# α_M(a) from the same Noether/segment normalization used for a0
# α_M(a) = κ · ξ · ε(a)  with κ≈2 fixed by operator counting (quadratic gradient),
# and ξ fixed by the same Noether normalization used for the 5/12 weak-field factor.
# For this runner we take a conservative, theory-locked ξ=2.5 (no fitting).
# ------------------------
def alphaM_of_a_factory(epsilon, kappa=2.0, xi=2.5):
    def alphaM(a):
        # nearest grid sample
        # (a_grid assumed monotonic; you can replace with interpolation if desired)
        return kappa*xi*epsilon(a)
    return alphaM

def mu_of_a_factory(epsilon, eta=5.0/12.0, lam=1.0):
    """
    Effective Poisson coupling μ(a)=G_eff/G_N driven by throttling.
    Use a smooth *saturating* law instead of a hard floor:
        μ(a) = 1 / (1 + η ε(a))               (default)
    which has the correct linearized limit μ ≈ 1 - η ε for small ε,
    stays positive, and asymptotes as ε grows.  If you prefer a Padé form,
    set `lam` accordingly and use the commented line below.
    """
    def mu(a):
        e = max(0.0, float(epsilon(a)))
        # Simple saturating map:
        return 1.0 / (1.0 + eta*e)
        # Padé alternative (commented):
        # return 1.0 - (eta*e)/(1.0 + lam*e)
    return mu

def rescale_alphaM_to_gw_bound(a_grid, alphaM_vec, z_max=1000.0, tol=1e-3):
    la = np.log(a_grid)
    a_min = 1.0/(1.0+z_max)
    i0 = np.searchsorted(a_grid, a_min)
    integ = np.trapezoid(alphaM_vec[i0:], la[i0:])
    if integ == 0.0:
        scale = 0.0
    else:
        # d_GW/d_EM ≈ exp(-0.5 ∫ α_M d ln a) → target integral ≲ 2*tol
        target = 2.0*tol
        scale = min(1.0, target/abs(integ)) * (1.0 if integ > 0 else -1.0)
    return alphaM_vec*scale, scale

def interp_1d(x, y):
    # simple linear interpolation closure over a monotonic grid
    def f(xq):
        if xq <= x[0]: return y[0]
        if xq >= x[-1]: return y[-1]
        i = np.searchsorted(x, xq) - 1
        t = (xq - x[i])/(x[i+1]-x[i])
        return y[i]*(1-t) + y[i+1]*t
    return f

# ------------------------
# Observables and plots
# ------------------------
def sigma8_today():  # just adopt the Planck-like σ8; we report S8 shifts relatively
    return 0.83

def S8_from_D(D, a_grid):
    # S8 ∝ σ8 * sqrt(Ω_m0/0.3); relative changes track D(a=1), which we've normalized to 1
    # So S8 baseline is ~0.83 for our H0, Om0 (we print absolute values for readability)
    return 0.83*np.sqrt(Om0/0.3)

def fs8_of_z(a_grid, D, sigma8_0):
    """fσ8(z) from the growth solution.
    σ8(z) = σ8,0 * D(z)/D(1); f = d ln D / d ln a."""
    la  = np.log(a_grid)
    dD  = np.gradient(D, la)
    f   = dD / D
    z   = 1.0/a_grid - 1.0
    sigma8_z = sigma8_0 * (D / D[-1])
    fs8 = f * sigma8_z
    return z[::-1], fs8[::-1]

def gw_em_ratio(alphaM_grid):
    """
    Very conservative proxy: d_GW/d_EM ≈ exp( -0.5 ∫ α_M d ln a ) → report max deviation over range.
    (No extra tensor sector; this friction-only mapping is tiny for our ε values.)
    """
    a = alphaM_grid['a']
    am= alphaM_grid['alphaM']
    la= np.log(a)
    integ = np.trapezoid(am, la)
    ratio = math.exp(-0.5*integ)
    return ratio

def run_and_report():
    # 1) Grid and LCDM background
    a_grid = np.linspace(1e-3, 1.0, 1200)

    # 2) Pull Ω_Λ from pipeline prediction (not observational):
    omega_from_pipeline = 0.683474127  # you can also read your results.json; kept inline for portability

    # 3) Growth with α_M=0 to source J(a)
    alphaM_zero = interp_1d(a_grid, np.zeros_like(a_grid))
    D_LCDM = solve_growth(a_grid, alphaM_zero, mu_of_a=None, normalize=False)

    # 4) Entropic state-action: J(a) with kernel p=5, calibrate c_log,J_*
    J = cumulative_J(a_grid, D_LCDM, p=5)
    eps0 = 0.01  # small irreversibility floor (irreversible hysteresis; >=0 by 2nd law)
    c_log, Jstar = calibrate_normalization(J, a_grid, omega_from_pipeline, eps0=eps0)
    eps_vec_raw = epsilon_of_a(J, a_grid, c_log, Jstar, eps0=eps0)
    # Enforce ΔS ≥ 0 (monotone ε) by clamping tiny numerical dips; report if any
    eps_vec = np.maximum.accumulate(eps_vec_raw)
    min_step = float(np.min(np.diff(eps_vec_raw)))
    if min_step < -1e-8:
        print(f"[note] ε(a) had small negative step(s) down to {min_step:.2e}; "
              "clamped to satisfy coarse-grained ΔS ≥ 0.")
    eps     = interp_1d(a_grid, eps_vec)

    # Diagnostics: extrema of ε and μ(a)
    eps_min, eps_max = float(eps_vec[0]), float(eps_vec[-1])
    mu_debug = np.array([mu_of_a_factory(lambda aa: eps(aa), eta=5.0/12.0)(ai) for ai in a_grid])
    mu_min, mu_max = float(mu_debug.min()), float(mu_debug.max())

    # 5) α_M(a) and μ(a)
    # Build an α_M(a) ∝ ε(a) from the same Noether mapping and *rescale* to satisfy GW bound.
    kappa, xi = 2.0, 2.5
    alphaM_vec_raw = kappa*xi*eps_vec
    alphaM_vec, alphaM_scale = rescale_alphaM_to_gw_bound(a_grid, alphaM_vec_raw, z_max=1000.0, tol=5e-3)  # 0.5% GW bound over full history
    alphaM = interp_1d(a_grid, alphaM_vec)
    # Weak-field normalization from the same Noether/5/12 factor used in a0 (η = 5/12 ≈ 0.4167), with smooth saturation
    mu_of_a = mu_of_a_factory(eps, eta=5.0/12.0)

    # 6) Growth with hysteresis α_M(a) and μ(a)
    D = solve_growth(a_grid, alphaM, mu_of_a=mu_of_a, normalize=False)

    # 7) Observables
    S8_LCDM = 0.83*np.sqrt(Om0/0.3)
    S8      = S8_LCDM * (D[-1]/D_LCDM[-1])
    z, fs8  = fs8_of_z(a_grid, D, sigma8_0=S8_LCDM*(D[-1]/D_LCDM[-1]))

    # 8) Diagnostics
    gw_ratio = gw_em_ratio({'a':a_grid, 'alphaM':alphaM_vec})
    s8_drop  = (S8/S8_LCDM - 1.0)*100.0

    # --- Robustness: S8 vs kernel power p ∈ {4,5,6} ---
    def compute_S8_for_p(p_val):
        Jp = cumulative_J(a_grid, D_LCDM, p=p_val)
        c_log_p, Jstar_p = calibrate_normalization(Jp, a_grid, omega_from_pipeline, eps0=eps0)
        eps_vec_p = epsilon_of_a(Jp, a_grid, c_log_p, Jstar_p, eps0=eps0)
        eps_vec_p = np.maximum.accumulate(eps_vec_p)
        # α_M and μ for this p
        alphaM_vec_raw_p = kappa*xi*eps_vec_p
        alphaM_vec_p, _ = rescale_alphaM_to_gw_bound(a_grid, alphaM_vec_raw_p, z_max=2.0, tol=5e-3)
        alphaM_p = interp_1d(a_grid, alphaM_vec_p)
        mu_p     = mu_of_a_factory(interp_1d(a_grid, eps_vec_p), eta=5.0/12.0)
        Dp = solve_growth(a_grid, alphaM_p, mu_of_a=mu_p, normalize=False)
        return float(S8_LCDM * (Dp[-1]/D_LCDM[-1]))
    p_vals = [4,5,6]
    S8_p   = {int(p): compute_S8_for_p(p) for p in p_vals}
    # Figure: S8 vs p
    plt.figure(figsize=(5.0,3.2))
    plt.axhspan(0.77,0.79, color="#3b82f6", alpha=0.15, label="S8 obs. band")
    xs = sorted(S8_p.keys())
    ys = [S8_p[k] for k in xs]
    plt.plot(xs, ys, marker="o", lw=2, label="model")
    plt.axhline(S8_LCDM, ls="--", color="k", lw=0.8, label="ΛCDM")
    plt.xlabel("kernel power p")
    plt.ylabel("S8")
    plt.xticks(xs)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig("s8_p_sweep.png", dpi=180)
    plt.close()

    # 8) Write summary JSON
    summary = {
        "model_label": MODEL_LABEL,
        "summary_file": SUMMARY_JSON_BASENAME,
        "Omega_m0": Om0,
        "Omega_Lambda_pipeline": omega_from_pipeline,
        "H0_km_s_Mpc": H0,
        "kernel_p": 5,
        "Jstar": float(Jstar),
        "c_log": float(c_log),
        "eps0": float(eps0),
        "kappa": 2.0,
        "xi_Noether": 2.5,
        "S8_today": float(S8),
        "S8_LambdaCDM": float(S8_LCDM),
        "S8_percent_change_vs_LCDM": float(s8_drop),
        "GW_EM_distance_ratio": float(gw_ratio),
        "max_abs_dGW_over_dEM_minus_1": float(abs(gw_ratio-1.0)),
        "S8_by_kernel_p": S8_p,
        "mu_model": "mu=1/(1+eta*epsilon), eta=5/12",
    }
    with open(SUMMARY_JSON_BASENAME,"w") as f:
        json.dump(summary, f, indent=2)

    # 9) Plots (same style as your referee set)
    # (a) S8 — band + lines (improved readability)
    plt.figure(figsize=(7.2, 4.2))
    plt.axhspan(0.77, 0.79, color="#60a5fa", alpha=0.18, label="S8 observational band")
    plt.axhline(S8_LCDM, color="0.20", ls="--", lw=1.4, label="ΛCDM baseline")
    plt.axhline(S8, color="C0", lw=2.6, label=MODEL_LABEL)

    ax = plt.gca()
    info_txt = (
        f"S8 = {S8:.3f}  (Δ = {s8_drop:.1f}% vs ΛCDM)\n"
        f"max|d_GW/d_EM−1| ≈ {abs(gw_ratio-1.0):.3e}\n"
        f"ε0 = {eps0:.3g},  c_log = {c_log:.3g}"
    )
    ax.text(
        0.02, 0.98, info_txt,
        transform=ax.transAxes,
        fontsize=12, ha="left", va="top",
        bbox=dict(facecolor="white", alpha=0.9, edgecolor="0.8")
    )

    plt.ylim(0.74, 0.86)
    plt.ylabel("S8", fontsize=12)
    plt.xticks([])
    plt.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=11)
    plt.tight_layout()
    plt.savefig("s8_bestfit_lines.png", dpi=240)
    plt.close()

    # (b) fσ8(z)
    plt.figure(figsize=(5.4,3.2))
    z_l, fs8_l = fs8_of_z(a_grid, D_LCDM, sigma8_0=S8_LCDM)
    plt.plot(z_l, fs8_l, 'k--', lw=1, label="ΛCDM")
    plt.plot(z,   fs8,   lw=2,  label=MODEL_LABEL)
    # Optional overlay if a CSV is available: columns z,fs8,err,label (label optional)
    csv_path = ensure_rsds_csv("fs8_compilation.csv")
    try:
        data = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=None, encoding=None)
        z_col   = np.array(data["z"],   dtype=float)
        fs8_col = np.array(data["fs8"], dtype=float)
        err_col = np.array(data["err"], dtype=float)
        labels  = data["label"] if "label" in data.dtype.names else np.array(["RSD"]*len(z_col))
        # Plot each label group separately for a clean legend.
        for lab in np.unique(labels):
            m = (labels == lab)
            plt.errorbar(z_col[m], fs8_col[m], yerr=err_col[m],
                         fmt="o", ms=4, lw=1, alpha=0.9, label=str(lab))
    except Exception as e:
        print(f"[warn] Could not parse '{csv_path}' for overlay: {e}")
    plt.gca().invert_xaxis()
    plt.xlabel("z")
    plt.ylabel("fσ8")
    plt.legend(frameon=False, loc="upper left")
    plt.tight_layout()
    plt.savefig("fs8_comparison.png", dpi=180)
    plt.close()

    # (c) E(z) check (unchanged here by construction)
    z_plot = np.linspace(0,2,200)
    plt.figure(figsize=(5.4,3.2))
    plt.plot(z_plot, [E_of_a(1/(1+zz)) for zz in z_plot], lw=2)
    plt.ylabel("E(z)")
    plt.xlabel("z")
    plt.tight_layout()
    plt.savefig("E_of_z_check.png", dpi=180)
    plt.close()

    # (d) GW/EM ratio
    plt.figure(figsize=(5.4,3.2))
    a = a_grid
    la= np.log(a_grid)
    # cumulative ratio from 1 to a
    cum = np.cumsum(alphaM_vec*np.gradient(la))
    ratio_curve = np.exp(-0.5*(cum - cum[-1]))  # normalize to 1 at a=1
    plt.plot(1/a - 1, ratio_curve, lw=2)
    ax = plt.gca()
    ax.axhline(1.0 + 5e-3, color="r", ls=":", lw=1)
    ax.axhline(1.0 - 5e-3, color="r", ls=":", lw=1)
    ax.text(0.03, 0.90, f"maxΔ≈{abs(gw_ratio-1.0):.3e}", transform=ax.transAxes,
            fontsize=9, ha="left", va="top")
    plt.gca().invert_xaxis()
    plt.xlabel("z")
    plt.ylabel("d_GW/d_EM (≈1; GR-like)")
    plt.tight_layout()
    plt.savefig("gw_em_ratio.png", dpi=180)
    plt.close()

    # (e) Residuals vs LCDM (toy)
    plt.figure(figsize=(5.4,3.2))
    diff = S8 - S8_LCDM
    plt.axhspan(0.75 - S8_LCDM, 0.79 - S8_LCDM, color="#3b82f6", alpha=0.15)
    plt.axhline(0.0, color="k", lw=0.8, ls="--")
    plt.bar([0],[diff], width=0.4, color="C0")
    plt.ylabel("S8 - S8(ΛCDM)")
    plt.xticks([])
    plt.tight_layout()
    plt.savefig("s8_residuals.png", dpi=180)
    plt.close()

    print(f"=== {MODEL_LABEL} run summary ===")
    print(f"S8 (state-action): {S8:.3f}   |  S8(ΛCDM): {S8_LCDM:.3f}   (Δ={s8_drop:.1f}% vs ΛCDM)")
    print(f"GW/EM distance ratio (0<z<~1000): {gw_ratio:.6f}   |  |Δ|≈{abs(gw_ratio-1.0):.3e}")
    print(f"Wrote: {SUMMARY_JSON_BASENAME} and figures.")
    print("S8 robustness vs kernel p:", ", ".join([f"p={k}: {v:.3f}" for k,v in S8_p.items()]))
    print(f"α_M overall scale to satisfy GW bound (±0.5%): {alphaM_scale:.3e}")
    print(f"ε(a) ∈ [{eps_min:.4f}, {eps_max:.4f}]  |  μ(a) ∈ [{mu_min:.4f}, {mu_max:.4f}]")

if __name__ == "__main__":
    run_and_report()