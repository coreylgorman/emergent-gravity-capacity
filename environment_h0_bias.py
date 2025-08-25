#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
environment_h0_bias.py — First-principles THEORY+ H0 adjustment with weak-field β, entropic mapping, and SALT-mapped SN physics.

DEFAULT (no CLI): THEORY+ with SN-cap = 0.05 mag, EM geometry GR-like (alpha_M=0), GW/EM=1.
Auto-discovers a CSV at ./data/host_catalog.csv (or similar) unless H0_HOST_CSV is set.
"""

from __future__ import annotations
import argparse, dataclasses, importlib.util, json, math, os, sys
from typing import Dict, List, Tuple, Optional
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

DEFAULT_HOST_CSV_CANDIDATES = (
    "./data/host_catalog.csv",
    "./data/hosts.csv",
    "./data/shoes_hosts.csv",
    "./data/mock_shoes_hosts.csv",
    "./mock_shoes_hosts.csv",
    "./hosts.csv",
)

@dataclasses.dataclass
class HostCategory:
    name: str
    g_over_a0: float

@dataclasses.dataclass
class MixDefinition:
    name: str
    fractions: Dict[str, float]

@dataclasses.dataclass
class AnalysisConfig:
    eta: float = 5.0/12.0
    eps0: float = 0.01
    kernel_p: int = 5

    gate_exponents: Tuple[int, ...] = (1,2,3,5)
    K_SN_grid: Tuple[float, ...] = tuple(np.linspace(-2.0, -0.6, 8))
    K_Ceph_grid: Tuple[float, ...] = tuple(np.linspace(-0.3, 0.3, 7))

    # Theory+ SALT / microphysics (first principles; no fitted amplitude)
    alpha_salt: float = 0.14
    beta_salt: float  = 3.1
    gamma_ni: float   = 0.6
    s_t: float        = 6.0
    c_t: float        = 0.02
    K_Ceph_th: float  = 0.0  # Cepheid rung coefficient (default 0; auto-clipped to same-host cap if nonzero)

    cap_Ceph_mag: float = 0.03
    cap_SN_mag: Optional[float] = None

    H0_Planck: float = 67.4
    H0_SHOES: float  = 73.0
    H0_TRGB: float   = 70.4

    outdir: str = "./outputs"

    # CSV ingestion
    host_csv: Optional[str] = None
    col_sample: str = "sample"
    col_weight: Optional[str] = None
    col_g_over_a0: Optional[str] = None
    col_vcirc_kms: Optional[str] = None
    col_radius_kpc: Optional[str] = None
    col_Menc_Msun: Optional[str] = None
    col_Sigma_Msunpc2: Optional[str] = None
    a0_m_s2: float = 1.2e-10

    categories: Tuple[HostCategory, ...] = (
        HostCategory("bulge", 10.0),
        HostCategory("arms", 3.0),
        HostCategory("outskirts", 0.5),
        HostCategory("intergroup", 0.1),
    )
    mixes: Tuple[MixDefinition, ...] = (
        MixDefinition("baseline_Cal", {"bulge":0.30,"arms":0.50,"outskirts":0.20,"intergroup":0.00}),
        MixDefinition("baseline_HF",  {"bulge":0.05,"arms":0.40,"outskirts":0.40,"intergroup":0.15}),
    )

def ensure_outdir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def F_gate(g_over_a0: float, n: int) -> float:
    x = float(g_over_a0)
    return 1.0 / (1.0 + x**n)

def mu_from_eps(eps: float, eta: float) -> float:
    return 1.0 / (1.0 + eta*eps)

def write_invariants(cfg: AnalysisConfig, outdir: str) -> None:
    inv = {
        "EM_geometry_GR_like": True,
        "alphaM_in_growth_solver": 0.0,
        "GW_over_EM_ratio": 1.0,
        "uses_beta_gating_only_in_hosts": True,
        "Planck_regime_beta_active": False,
        "SN_host_cap_mag": cfg.cap_SN_mag,
        "Cepheid_same_host_cap_mag": cfg.cap_Ceph_mag,
    }
    with open(os.path.join(outdir,"invariants.json"),"w") as f:
        json.dump(inv, f, indent=2)

def load_s8_run_module() -> object:
    try:
        import s8_run  # type: ignore
        return s8_run
    except Exception:
        pass
    for p in [os.path.join(os.getcwd(),"s8_hysteresis_run.py"), "/mnt/data/s8_hysteresis_run.py"]:
        if os.path.exists(p):
            spec = importlib.util.spec_from_file_location("s8_run", p)
            module = importlib.util.module_from_spec(spec)  # type: ignore
            assert spec and spec.loader
            spec.loader.exec_module(module)  # type: ignore
            sys.modules["s8_run"] = module
            return module
    raise ImportError("Could not import s8_run; place s8_hysteresis_run.py in CWD or PYTHONPATH.")

def compute_epsilon_today(s8_run: object, eps0: float, p: int) -> float:
    a = np.linspace(1e-3, 1.0, 2000)
    alphaM_zero = s8_run.interp_1d(a, np.zeros_like(a))  # α_M ≡ 0 → GR EM distances
    D = s8_run.solve_growth(a, alphaM_zero, mu_of_a=None, normalize=False)
    ln_a = np.log(a); D2 = D**2
    J = np.zeros_like(a)
    trapz_fn = getattr(np, "trapezoid", np.trapz)
    for j in range(1, len(a)):
        aj = a[j]
        K = (a[:j] / aj)**p
        J[j] = trapz_fn(K * D2[:j], ln_a[:j])
    J = np.maximum.accumulate(J)
    try:
        s8_run.np.trapezoid = trapz_fn  # type: ignore
    except Exception:
        pass
    c_log, Jstar = s8_run.calibrate_normalization(J, a, omega_from_pipeline=0.683474127, eps0=eps0)
    eps_arr = s8_run.epsilon_of_a(J, a, c_log, Jstar, eps0=eps0)
    return float(eps_arr[-1])

def epsilon_env(epsilon_today: float, eps0: float, g_over_a0: float, n: int) -> float:
    return eps0 + (epsilon_today - eps0) * F_gate(g_over_a0, n)

def lnmu_for_g(epsilon_today: float, eps0: float, eta: float, g_over_a0: float, n: int) -> float:
    eps_env = epsilon_env(epsilon_today, eps0, g_over_a0, n)
    mu = mu_from_eps(eps_env, eta)
    return float(np.log(mu))

def _normalize_sample_label(x: str) -> Optional[str]:
    if not isinstance(x, str): return None
    t = x.strip().lower()
    if t in {"cal","calibrator","calibrators","ceph","cepheid"}:
        return "Cal"
    if t in {"hf","hubbleflow","hubble_flow","hflow","hubble"}:
        return "HF"
    return None

def _compute_g_over_a0_from_row(row: pd.Series, cfg: AnalysisConfig) -> Optional[float]:
    a0 = cfg.a0_m_s2
    if cfg.col_g_over_a0 and (cfg.col_g_over_a0 in row) and pd.notnull(row[cfg.col_g_over_a0]):
        try:
            return float(row[cfg.col_g_over_a0])
        except Exception:
            pass
    G=6.67430e-11; MSUN=1.98847e30; KPC=3.085677581e19; KM=1000.0; PC=3.085677581e16
    if cfg.col_vcirc_kms and cfg.col_radius_kpc and (cfg.col_vcirc_kms in row) and (cfg.col_radius_kpc in row):
        if pd.notnull(row[cfg.col_vcirc_kms]) and pd.notnull(row[cfg.col_radius_kpc]):
            try:
                V = float(row[cfg.col_vcirc_kms])*KM; R = float(row[cfg.col_radius_kpc])*KPC
                return float((V*V/R)/a0)
            except Exception:
                pass
    if cfg.col_Menc_Msun and cfg.col_radius_kpc and (cfg.col_Menc_Msun in row) and (cfg.col_radius_kpc in row):
        if pd.notnull(row[cfg.col_Menc_Msun]) and pd.notnull(row[cfg.col_radius_kpc]):
            try:
                M = float(row[cfg.col_Menc_Msun])*MSUN; R = float(row[cfg.col_radius_kpc])*KPC
                return float((G*M/(R*R))/a0)
            except Exception:
                pass
    if cfg.col_Sigma_Msunpc2 and (cfg.col_Sigma_Msunpc2 in row) and pd.notnull(row[cfg.col_Sigma_Msunpc2]):
        try:
            Sigma = float(row[cfg.col_Sigma_Msunpc2])*MSUN/(PC*PC)
            return float((2.0*math.pi*G*Sigma)/a0)
        except Exception:
            pass
    return None

def load_host_csv(cfg: AnalysisConfig) -> pd.DataFrame:
    assert cfg.host_csv and os.path.exists(cfg.host_csv), "host_csv must exist."
    df = pd.read_csv(cfg.host_csv)
    if cfg.col_sample not in df.columns:
        raise ValueError(f"CSV missing sample column {cfg.col_sample!r}")
    df["group"] = df[cfg.col_sample].apply(_normalize_sample_label)
    df = df[df["group"].notnull()].copy()
    gvals = []
    for _, row in df.iterrows():
        g_over_a0 = _compute_g_over_a0_from_row(row, cfg)
        gvals.append(np.nan if g_over_a0 is None else g_over_a0)
    df["g_over_a0"] = gvals
    df = df[pd.notnull(df["g_over_a0"])].copy()
    if cfg.col_weight and (cfg.col_weight in df.columns):
        w = pd.to_numeric(df[cfg.col_weight], errors="coerce").fillna(0.0).to_numpy()
    else:
        w = np.ones(len(df), dtype=float)
    df["weight"] = w
    return df

def ladder_bias_delta_mu(K_SN: float, K_Ceph: float, lnmu_cal: float, lnmu_hf: float) -> float:
    return float(K_SN * (lnmu_hf - lnmu_cal) + K_Ceph * lnmu_cal)

def delta_H_over_H_from_mag(delta_mu_mag: float) -> float:
    return - (math.log(10.0) / 5.0) * float(delta_mu_mag)

def cepheid_cap_ok(K_Ceph: float, lnmu_cal: float, cap_mag: float) -> bool:
    return abs(K_Ceph * lnmu_cal) <= cap_mag + 1e-12

def sn_host_cap_ok(K_SN: float, lnmu_hf: float, lnmu_cal: float, cap_mag: Optional[float]) -> bool:
    if cap_mag is None: return True
    return abs(K_SN * (lnmu_hf - lnmu_cal)) <= cap_mag + 1e-12

def K_SN_theory(alpha_SN: float) -> float:
    return (2.5 / math.log(10.0)) * 1.5 * alpha_SN

def K_SN_eff_theoryplus(alpha_salt: float, beta_salt: float, gamma_ni: float, s_t: float, c_t: float) -> float:
    return 1.6286*(gamma_ni - 0.5) - 0.75*alpha_salt*s_t - beta_salt*c_t

def get_lnmu_cal_hf(cfg: AnalysisConfig, s8_run: object, n: int) -> Tuple[float, float]:
    epsilon_today = compute_epsilon_today(s8_run, cfg.eps0, cfg.kernel_p)
    if cfg.host_csv:
        df = load_host_csv(cfg)
        df["lnmu"] = [lnmu_for_g(epsilon_today, cfg.eps0, cfg.eta, float(x), n) for x in df["g_over_a0"]]
        lnmu_cal = float(np.average(df.loc[df["group"]=="Cal","lnmu"], weights=df.loc[df["group"]=="Cal","weight"]))
        lnmu_hf  = float(np.average(df.loc[df["group"]=="HF","lnmu"],  weights=df.loc[df["group"]=="HF","weight"]))
        return lnmu_cal, lnmu_hf
    # fallback toy-mix
    cat = {c.name: c for c in cfg.categories}
    mixes = {m.name: m for m in cfg.mixes}
    lnmu_cal = sum(mixes["baseline_Cal"].fractions[k]*lnmu_for_g(epsilon_today, cfg.eps0, cfg.eta, cat[k].g_over_a0, n)
                   for k in mixes["baseline_Cal"].fractions)
    lnmu_hf  = sum(mixes["baseline_HF"].fractions[k]*lnmu_for_g(epsilon_today, cfg.eps0, cfg.eta, cat[k].g_over_a0, n)
                   for k in mixes["baseline_HF"].fractions)
    return float(lnmu_cal), float(lnmu_hf)

def run_theoryplus(cfg: AnalysisConfig, s8_run: object) -> pd.DataFrame:
    n = cfg.kernel_p
    lnmu_cal, lnmu_hf = get_lnmu_cal_hf(cfg, s8_run, n)
    Delta_lnmu = lnmu_hf - lnmu_cal
    K_eff = K_SN_eff_theoryplus(cfg.alpha_salt, cfg.beta_salt, cfg.gamma_ni, cfg.s_t, cfg.c_t)

    # 1) Raw (SN only)
    dmu_raw_no_ceph = ladder_bias_delta_mu(K_eff, 0.0, lnmu_cal, lnmu_hf)
    dH_raw_no_ceph  = delta_H_over_H_from_mag(dmu_raw_no_ceph)
    H0_raw_no_ceph  = cfg.H0_SHOES * (1.0 + dH_raw_no_ceph)

    # 2) Cap (SN only)
    if (cfg.cap_SN_mag is not None) and (abs(Delta_lnmu) > 0):
        K_cap = math.copysign(min(abs(K_eff), cfg.cap_SN_mag/abs(Delta_lnmu)), K_eff)
    else:
        K_cap = K_eff
    dmu_cap_sn_only = ladder_bias_delta_mu(K_cap, 0.0, lnmu_cal, lnmu_hf)
    dH_cap_sn_only  = delta_H_over_H_from_mag(dmu_cap_sn_only)
    H0_cap_sn_only  = cfg.H0_SHOES * (1.0 + dH_cap_sn_only)

    # 3) Cap + Cepheid (Cepheid auto-clipped to same-host cap)
    K_Ceph_eff = cfg.K_Ceph_th or 0.0
    if lnmu_cal != 0.0 and (cfg.cap_Ceph_mag is not None):
        if abs(K_Ceph_eff * lnmu_cal) > cfg.cap_Ceph_mag:
            K_Ceph_eff = math.copysign(cfg.cap_Ceph_mag/abs(lnmu_cal), K_Ceph_eff if K_Ceph_eff != 0 else 1.0)
    dmu_cap_sn_plus_ceph = ladder_bias_delta_mu(K_cap, K_Ceph_eff, lnmu_cal, lnmu_hf)
    dH_cap_sn_plus_ceph  = delta_H_over_H_from_mag(dmu_cap_sn_plus_ceph)
    H0_cap_sn_plus_ceph  = cfg.H0_SHOES * (1.0 + dH_cap_sn_plus_ceph)

    return pd.DataFrame([{
        "mode":"theoryplus","n":n,
        "lnmu_cal":lnmu_cal,"lnmu_hf":lnmu_hf,"Delta_lnmu":Delta_lnmu,
        "K_SN_eff_raw":K_eff,"K_SN_eff_cap":K_cap,"K_Ceph_eff":K_Ceph_eff,
        "Delta_mu_raw_no_ceph":dmu_raw_no_ceph,
        "Delta_mu_cap_sn_only":dmu_cap_sn_only,
        "Delta_mu_cap_sn_plus_ceph":dmu_cap_sn_plus_ceph,
        "Delta_H_over_H_raw_no_ceph":dH_raw_no_ceph,
        "Delta_H_over_H_cap_sn_only":dH_cap_sn_only,
        "Delta_H_over_H_cap_sn_plus_ceph":dH_cap_sn_plus_ceph,
        "H0_raw_no_ceph":H0_raw_no_ceph,
        "H0_cap_sn_only":H0_cap_sn_only,
        "H0_cap_sn_plus_ceph":H0_cap_sn_plus_ceph,
    }])

def plot_H0_points(cfg: AnalysisConfig, pts: Dict[str, float], outpath: str, cap_text: str = "") -> None:
    labels = [
        "Planck (CMB)", "TRGB (CCHP)", "SH0ES (Cepheid/SN)",
        "Theory+ raw (SN only)", f"Theory+ cap (SN only){cap_text}",
        f"Theory+ cap+Ceph{cap_text}"
    ]
    xpos = [cfg.H0_Planck, cfg.H0_TRGB, cfg.H0_SHOES,
            pts["H0_raw_no_ceph"], pts["H0_cap_sn_only"], pts["H0_cap_sn_plus_ceph"]]
    ypos = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8.8,4.2), dpi=150)
    for x,y,m in zip(xpos, ypos, ["s","o","^","D","D","D"]):
        ax.scatter([x],[y], s=50, marker=m)
        ax.text(x+0.18, y+0.10, f"{x:.2f}", va="center")
    ax.set_xlabel("H₀  (km s⁻¹ Mpc⁻¹)"); ax.set_yticks(ypos); ax.set_yticklabels(labels)
    ax.set_title("H₀ comparison: Planck / TRGB / SH0ES vs Theory+ (first principles)")
    ax.grid(True, axis="x", linestyle=":", linewidth=0.8)
    fig.tight_layout(); fig.savefig(outpath, bbox_inches="tight"); plt.close(fig)

def run_default_theoryplus(cfg: AnalysisConfig, s8_run: object) -> None:
    cfg.cap_SN_mag = float(os.environ.get("H0_SN_CAP", "0.05"))
    if cfg.host_csv and os.path.exists(cfg.host_csv):
        try:
            _df = pd.read_csv(cfg.host_csv, nrows=5)
            if "g_over_a0" in _df.columns and (cfg.col_g_over_a0 is None):
                cfg.col_g_over_a0 = "g_over_a0"
            if "w" in _df.columns and (cfg.col_weight is None):
                cfg.col_weight = "w"
            if "sample" in _df.columns:
                cfg.col_sample = "sample"
            print(f"[DEFAULT] Using host CSV: {cfg.host_csv}")
        except Exception as _e:
            print(f"[DEFAULT] Warning: could not pre-scan host CSV header ({_e}).")
    else:
        print("[DEFAULT] No host CSV detected; using toy mixes. Set H0_HOST_CSV or place CSV under ./data/.")
    for attr, envname in [("alpha_salt","H0_ALPHA_SALT"),("beta_salt","H0_BETA_SALT"),("gamma_ni","H0_GAMMA_NI"),
                          ("s_t","H0_S_T"),("c_t","H0_C_T")]:
        if envname in os.environ:
            try: setattr(cfg, attr, float(os.environ[envname]))
            except ValueError: pass
    if "H0_K_CEPH" in os.environ:
        try: cfg.K_Ceph_th = float(os.environ["H0_K_CEPH"])
        except ValueError: pass

    ensure_outdir(cfg.outdir)
    df = run_theoryplus(cfg, s8_run)
    out_csv = os.path.join(cfg.outdir, "theoryplus_summary.csv")
    df.to_csv(out_csv, index=False)
    pts = {"H0_raw_no_ceph": float(df["H0_raw_no_ceph"].iloc[0]),
           "H0_cap_sn_only": float(df["H0_cap_sn_only"].iloc[0]),
           "H0_cap_sn_plus_ceph": float(df["H0_cap_sn_plus_ceph"].iloc[0])}
    fig = os.path.join(cfg.outdir, "H0_points_theoryplus.png")
    plot_H0_points(cfg, pts, fig, cap_text=f" (|Δm_SN|≤{cfg.cap_SN_mag:.02f} mag)")
    write_invariants(cfg, cfg.outdir)
    print("DEFAULT THEORY+ complete.")
    print(json.dumps(pts, indent=2))
    print("CSV:", out_csv)
    print("Plot:", fig)
    print("Invariants:", os.path.join(cfg.outdir, "invariants.json"))

def main(argv: Optional[List[str]] = None) -> None:
    if argv is None: argv = sys.argv[1:]
    if len(argv) == 0:
        host_csv_env = os.environ.get("H0_HOST_CSV", None)
        if not host_csv_env:
            for _cand in DEFAULT_HOST_CSV_CANDIDATES:
                if os.path.exists(_cand):
                    host_csv_env = _cand
                    break
        cfg = AnalysisConfig(outdir="./outputs_paper_ready", host_csv=host_csv_env)
        s8_run = load_s8_run_module(); run_default_theoryplus(cfg, s8_run); return

    parser = argparse.ArgumentParser(description="First-principles Theory+ H0 adjustment (scan | theory | theoryplus).")
    sub = parser.add_subparsers(dest="mode", required=True)
    def add_common(p):
        p.add_argument("--outdir", type=str, default="./outputs", help="Output directory.")
        p.add_argument("--cepheid-cap", type=float, default=0.03, help="Same-host Cepheid cap (mag).")
        p.add_argument("--planck", type=float, default=67.4); p.add_argument("--shoes", type=float, default=73.0)
        p.add_argument("--trgb", type=float, default=70.4)
        p.add_argument("--host-csv", type=str, default=None)
        p.add_argument("--col-sample", type=str, default="sample")
        p.add_argument("--col-weight", type=str, default=None)
        p.add_argument("--col-g-over-a0", type=str, default=None)
        p.add_argument("--col-vcirc-kms", type=str, default=None)
        p.add_argument("--col-radius-kpc", type=str, default=None)
        p.add_argument("--col-Menc-Msun", type=str, default=None)
        p.add_argument("--col-Sigma-Msunpc2", type=str, default=None)
        p.add_argument("--a0", type=float, default=1.2e-10)
        p.add_argument("--kernel-p", type=int, default=5)
        p.add_argument("--eta", type=float, default=5.0/12.0)
        p.add_argument("--eps0", type=float, default=0.01)

    p_tplus = sub.add_parser("theoryplus", help="First-principles Theory+ (Arnett + diffusion/opacity + SALT) → K_SN_eff")
    add_common(p_tplus)
    p_tplus.add_argument("--sn-cap", type=float, default=0.05)
    p_tplus.add_argument("--alpha-salt", type=float, default=0.14)
    p_tplus.add_argument("--beta-salt", type=float, default=3.1)
    p_tplus.add_argument("--gamma-ni", type=float, default=0.6)
    p_tplus.add_argument("--s-t", type=float, default=6.0)
    p_tplus.add_argument("--c-t", type=float, default=0.02)
    p_tplus.add_argument("--k-ceph", type=float, default=None)

    args = parser.parse_args(argv)
    ensure_outdir(args.outdir)
    cfg = AnalysisConfig(
        eta=args.eta, eps0=args.eps0, kernel_p=args.kernel_p, outdir=args.outdir,
        H0_Planck=args.planck, H0_SHOES=args.shoes, H0_TRGB=args.trgb,
        host_csv=args.host_csv, col_sample=args.col_sample, col_weight=args.col_weight,
        col_g_over_a0=args.col_g_over_a0, col_vcirc_kms=args.col_vcirc_kms,
        col_radius_kpc=args.col_radius_kpc, col_Menc_Msun=args.col_Menc_Msun,
        col_Sigma_Msunpc2=args.col_Sigma_Msunpc2, a0_m_s2=args.a0,
    )
    s8_run = load_s8_run_module()
    cfg.cap_Ceph_mag = args.cepheid_cap

    if args.mode == "theoryplus":
        cfg.cap_SN_mag = args.sn_cap
        cfg.alpha_salt = args.alpha_salt; cfg.beta_salt = args.beta_salt
        cfg.gamma_ni = args.gamma_ni; cfg.s_t = args.s_t; cfg.c_t = args.c_t
        if args.k_ceph is not None: cfg.K_Ceph_th = args.k_ceph
        df = run_theoryplus(cfg, s8_run)
        out_csv = os.path.join(cfg.outdir,"theoryplus_summary.csv"); df.to_csv(out_csv, index=False)
        pts = {"H0_raw_no_ceph": float(df["H0_raw_no_ceph"].iloc[0]),
               "H0_cap_sn_only": float(df["H0_cap_sn_only"].iloc[0]),
               "H0_cap_sn_plus_ceph": float(df["H0_cap_sn_plus_ceph"].iloc[0])}
        fig = os.path.join(cfg.outdir,"H0_points_theoryplus.png")
        plot_H0_points(cfg, pts, fig, cap_text=f" (|Δm_SN|≤{cfg.cap_SN_mag:.02f} mag)")
        write_invariants(cfg, cfg.outdir)
        print("THEORY+ results saved to:", out_csv); print(json.dumps(pts, indent=2))
        print("Plot:", fig); print("Invariants:", os.path.join(cfg.outdir,"invariants.json")); return

if __name__ == "__main__":
    main()
