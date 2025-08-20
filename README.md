# Emergent Gravity from Finite Information Capacity  
**Referee-Ready Reproducibility Pipeline**

---

## Overview

This repository provides a **referee-ready kit**: a self-contained, reproducible bundle of code, data, and the manuscript.  

It accompanies the paper:

> **Emergent State-Dependent Gravity from Local Information Capacity:  
> A Conditional Thermodynamic Derivation with Cosmological Predictions**  
> Corey Gorman, 2025 (preprint)

The goal is **transparency**: referees can (i) read the PDF, (ii) inspect the LaTeX source, and (iii) run the Python pipeline to reproduce every number quoted in the paper.

---

## What the Pipeline Does

The Python script `referee_pipeline.py` performs the following:

1. **QFT / MI Calculation of β**  
   - Uses the Casini–Huerta–Myers (CHM) modular Hamiltonian of a ball.  
   - Applies *mutual-information subtraction* (“moment-kill”) to eliminate UV area/contact divergences.  
   - Yields the dimensionless modular-response coefficient:  
     \[
     \beta = 2\pi\,C_T\,I_{00}, \qquad C_T=\tfrac{3}{\pi^4}.
     \]

2. **Geometric Normalization (f and c_geo)**  
   - Maps local Clausius flux to global FRW dynamics.  
   - Two bookkeeping conventions are supported (Scheme A and B).  
   - Only the **invariant product** β·f·c_geo is observable.

3. **Cosmology**  
   - Predicts the dark-energy fraction:  
     \[
     \Omega_\Lambda = \beta \, f \, c_{\rm geo}.
     \]  
   - Matches observation without introducing a bare cosmological constant.

4. **Phenomenology (MOND-like a₀)**  
   - Weak-field Clausius flux renormalization fixes the acceleration scale:  
     \[
     a_0 = \frac{\Omega_\Lambda^2}{2}\,c\,H_0,
     \]  
   - No new parameters required. Order of magnitude matches galactic rotation-curve data.

---

## Quick Start

Clone and run:

```bash
git clone https://github.com/<yourname>/emergent-gravity-capacity.git
cd emergent-gravity-capacity
pip install -r requirements.txt
python referee_pipeline.py# emergent-gravity-capacity
Emergent State-Dependent Gravity — Referee Pipeline
```
Default outputs (baseline)
	•	β ≈ 0.0209
	•	Scheme A: f = 0.8193, c_geo = 40
	•	Ω_Λ ≈ 0.685
	•	a₀ ≈ 1.6 × 10⁻¹⁰ m/s² (for H₀ = 70 km/s/Mpc)

Results are saved to results/run_<timestamp>/:
	•	results.json (machine-readable)
	•	summary.txt (human-readable log with equations)
	•	beta_sensitivity.csv (if sweep enabled)

Options
	•	Schemes
    python referee_pipeline.py --scheme both
    Shows both Scheme A and B → invariant product β f c_geo is the same.

	•	Hubble Constant
    python referee_pipeline.py --H0 74
    Recomputes a₀ for chosen H₀.

    •	Sensitivity Sweep
    python referee_pipeline.py --preset sweep
    Sweeps σ₁,σ₂,u_gap to confirm β plateau robustness. Results in beta_sensitivity.csv.

Report (PDF + LaTeX)

The full manuscript is included:
	•	report/emergentStateDepGravityFromLocalInfo.pdf
→ Referee-ready PDF.
	•	report/emergentStateDepGravityFromLocalInfo.tex
→ LaTeX source (compile with revtex4-2, natbib, AMS packages).

Referees are encouraged to check both the PDF narrative and the LaTeX derivations. Every number in the paper can be reproduced by running referee_pipeline.py.

⸻

Requirements
	•	Python 3.10+
	•	mpmath
	•	tqdm (optional)
Install with:
pip install -r requirements.txt

Sample Output

To help referees preview without running code, we include a curated sample run:
	•	results/x_summary.txt — Human-readable log
	•	results/x_results.json — Machine-readable values

These correspond to the baseline manuscript settings (Scheme A, H0=70).

License & Citation
	•	License: Apache 2.0 — attribution required (see NOTICE).
	•	Citation: see CITATION.cff.

If you use this code, results, or manuscript, please cite both the paper and this repository.

References
	•	Jacobson (1995), Thermodynamics of spacetime
	•	Casini–Huerta–Myers (2011), CHM modular Hamiltonian
	•	Iyer–Wald (1994), Noether charge entropy
	•	Osborn–Petkou (1994), Stress-tensor normalization

Note to Referees

This repository is intended as a transparent, reproducible submission package.
	•	The report/ folder contains the manuscript PDF and LaTeX source.
	•	The pipeline reproduces all quantitative results.
	•	The sensitivity sweep demonstrates robustness of β under parameter variations.

Only the invariant product β f c_geo is physical; internal bookkeeping differences (Scheme A vs B) are conventions.

We invite referees to review, run, and test. Clear falsifiers are described in the manuscript (GW/EM luminosity-distance ratio, bounds on \u02d9G/G).