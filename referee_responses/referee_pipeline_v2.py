"""
referee_pipeline_v2.py — helper utilities for scheme invariance and window scans.

This script does not recompute β from field theory; it uses the provided results.json and
derives cross-checks that are useful for referee exchanges.
"""
import json, math, sys, pathlib, argparse

def load_results(path="results.json"):
    with open(path, "r") as f:
        return json.load(f)

def cap_cos_from_cgeo(cgeo: float) -> float:
    # ΔΩ = 2π(1 - cosθ), c_geo = 4π / ΔΩ => c_geo = 2/(1 - cosθ)
    return 1.0 - 2.0/float(cgeo)

def cgeo_from_cos_theta(cos_theta: float) -> float:
    return 2.0/(1.0 - float(cos_theta))

def omega_lambda(beta: float, f: float, cgeo: float) -> float:
    return beta * f * cgeo

def a0_from_omega_lambda(omega: float, H0_km_s_Mpc: float) -> float:
    c = 299792458.0
    H0 = H0_km_s_Mpc * 1000.0 / 3.0856775814913673e22  # s^-1
    return (omega**2 / 2.0) * c * H0

def invariance_check(results_path="results.json", out_json="invariance_check.json"):
    data = load_results(results_path)
    beta = data["beta_outputs"]["beta"]
    f_A  = data["schemes"][0]["f_total"]
    omega_A = data["omega_lambda_by_scheme"]["A"]
    cgeo_A = omega_A / (beta * f_A)

    # Scheme B — fixed by convention
    f_B   = 3.125
    cgeo_B= 10.49

    omA = omega_lambda(beta, f_A, cgeo_A)
    omB = omega_lambda(beta, f_B, cgeo_B)
    a0A = a0_from_omega_lambda(omA, data["H0_km_s_Mpc"])

    out = {
        "beta": beta,
        "scheme_A": {"f": f_A, "cgeo": cgeo_A, "cos_theta": cap_cos_from_cgeo(cgeo_A)},
        "scheme_B": {"f": f_B, "cgeo": cgeo_B, "cos_theta": cap_cos_from_cgeo(cgeo_B)},
        "omega_lambda_A": omA,
        "omega_lambda_B": omB,
        "delta_frac": (omB - omA)/omA if omA else None,
        "a0_A_m_s2": a0A,
    }
    with open(out_json, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default="results.json")
    parser.add_argument("--out", default="invariance_check.json")
    args = parser.parse_args()
    invariance_check(args.results, args.out)
