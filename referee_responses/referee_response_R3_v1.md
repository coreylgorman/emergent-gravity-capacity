Response to Referee #3

Manuscript: Emergent State‑Dependent Gravity from Local Information Capacity: A Conditional Thermodynamic Derivation with Cosmological Predictions (RevTeX, Aug 20, 2025)

We thank the referee for their careful reading. Below we address each substantive point and indicate where the revised manuscript presents the relevant derivations, assumptions, and domain of validity.

⸻

1) Casini–Galante–Myers (CGM) and the marginal case Δ = d/2

Referee’s concern. The manuscript “uses only the scalar baseline with \Delta=d-2=2 in d=4 and does not address the problematic regime \Delta\le d/2.” The “marginal‑log compensator” is said to be hand‑waving.

Response. In d=4, \Delta=d/2=2 is exactly the marginal case we target. Our compensator term arises from a derived constitutive relation, not an ansatz:
\frac{\delta G}{G}=-\beta\,\delta\sigma,
obtained by extremizing a Wald‑like entropy functional under a local capacity constraint (Sec. V.B, Eqs. (14)–(16)). Varying S_{\rm grav}=A/[4G(x)] then yields
\delta S_{\rm grav}=\frac{1}{4G}\delta A+\frac{A}{4G}\beta\,\delta\sigma\quad\text{(Eqs. (11)–(13)).}
As noted in the scaling remark right below Eq. (13), the compensator picks up the requisite logarithmic variation in the marginal case \Delta=d/2, cancelling the CGM obstruction within the stated wedge safe window (Sec. III.B; App. H). This is not presented as a theorem for all operators in all theories; it is a conditional resolution in the same small‑wedge domain in which we apply Clausius.  ￼

⸻

2) Assumption (A2) and claims of circularity

Referee’s concern. The extension \delta Q=T\,\delta S to arbitrary local wedges is unproven; success in \Omega_\Lambda would not validate (A2).

Response. We agree that (A2) is a working assumption. The manuscript states this repeatedly and provides explicit bounds and falsifiers (Sec. II; Sec. VIII). We do not claim a proof of (A2), nor do we use the \Omega_\Lambda result to argue backwards. Instead, we delimit a domain where Clausius/Unruh plus MI‑subtracted modular response are self‑consistent (the “capacity safe window”). The Equivalence Principle for Modular Response (EPMR) is then articulated: to O(\ell^4) the isotropic modular coefficient equals its flat‑space value; curvature enters at O(\ell^6) (App. H: Lemma H.1, Lemma H.2, Proposition H.1). All quantitative claims are conditional on (A2); failure of (A2) would falsify or force revision of the framework.  ￼

⸻

3) “Scheme invariance failure” and numerical discrepancy

Referee’s concern. A 0.03% mismatch (Scheme A vs B) contradicts “invariance.”

Response. In the revised manuscript we (i) prove a continuous‑angle invariance (Sec. VI.C, Eqs. (18)–(20)) and (ii) show identical numerical values in the table: both schemes give \Omega_\Lambda=0.68493 for the same \beta. The earlier 0.03% mismatch arose from rounded inputs (e.g., \(\cgeo=10.49\) to two decimals) and not from physics. The invariance theorem makes the \theta‑dependence drop out exactly once a unit–solid–angle normalization is fixed.  ￼

⸻

4) Mutual‑information subtraction: parameter stability and “arbitrary choices”

Referee’s concern. Sensitivity to (\sigma_1,\sigma_2) and grid sizes; incomplete subtraction.

Response. Sec. IV.D reports scans over (\sigma_1,\sigma_2)\in[0.96,0.999]^2, u_{\rm gap}\in[0.2,0.35], and grids (N_r,N_s,N_\tau)\in[60,160]^3, with a clear plateau and |\Delta\beta|/\beta\lesssim0.5\%. The residual moments \sim10^{-51} reflect high‑precision cancellation well beyond our enforcement gates (<10^{-20}). The geometric factors f and c_{\rm geo} are not arbitrary fits; they are defined by ball→diamond mapping, Unruh normalization, and a ratio of Clausius fluxes with no cosmological input (Sec. VI; App. B–C). Only the product \beta f c_{\rm geo} is observable, and we now prove its invariance under cap angle (Sec. VI.C).  ￼

⸻

5) Flat‑space QFT → cosmology: “ad hoc” normalizations

Referee’s concern. The leap from flat modular response to FRW uses ad‑hoc f and c_{\rm geo}, possibly tuned.

Response. The bridge is split in two, by design:
(i) \beta is computed in flat space with MI subtraction and moment‑kill (Sec. III–IV).
(ii) The map to FRW is via Clausius/Noether with geometric normalizations: f_{\rm shape}=15/2 comes from explicit integrals (App. B.1); f_{\rm boost}=1 is the Unruh factor; f_{\rm cont}=1 follows because the finite MI‑subtracted coefficient is continuation‑invariant; and c_{\rm geo} is a ratio of (\delta Q/T) between FRW and the local wedge (App. C). We then prove that only \beta f c_{\rm geo} is physical and that it is independent of the cap angle (Sec. VI.C). No cosmological data enter this construction.  ￼

⸻

6) “Coincidental” agreement with Planck; “no theory errors”

Referee’s concern. The match \Omega_\Lambda\approx0.685 could be coincidental; theory error is absent.

Response. We present a two‑tier uncertainty:
	•	Numerical/systematic on \beta: \sim3\% total (Sec. IV.D), propagated linearly to \Omega_\Lambda (Sec. VII.B).
	•	Conceptual: all results are conditional on (A2); the abstract and conclusion explicitly state this.
We do not claim sub‑percent theoretical control; we report the central value and the numerical/systematic band, plus the conditional status of the derivation.  ￼

⸻

7) Computational details: grids, precision, gates

Referee’s concern. Specific grids used without convergence.

Response. Convergence is addressed by parameter scans (range listed above) and by gates (positivity, residual‑moment cancellations). We will be happy to attach a supplementary convergence table (same \beta within \lesssim0.5\% across grids) if the editor requests it; the manuscript already summarizes these tests succinctly in Sec. IV.D.  ￼

⸻

8) “Information capacity” is unclear

Referee’s concern. No operational definition; unclear vacuum subtraction.

Response. The paper avoids bare finite regional entropy. The state metric \sigma(x) is operationally defined as the finite \ell^4 coefficient of \delta\!\langle K_{\rm sub}\rangle/(2\pi C_T I_{00}) after MI subtraction and moment‑kill (Sec. III.A, Eq. (5)). “Vacuum subtraction’’ means subtracting Minkowski vacuum short‑distance entanglement; this is standard in modular‑Hamiltonian‑based approaches and is stated explicitly (Sec. IV.D, “Vacuum subtraction clarifier”).  ￼

⸻

9) Thermodynamic extension is “speculative”

Referee’s concern. Unruh temperature and equilibrium in finite, non‑stationary wedges are doubtful.

Response. We agree that this is the critical assumption. Hence we frame the program as conditional on (A2), and we provide a precise domain of applicability (the safe window, Sec. III.B) together with the EPMR formalism that shows how the flat‑space \ell^4 coefficient persists to leading order (App. H). The results are explicitly falsifiable through d_L^{\rm GW}/d_L^{\rm EM} and \dot G/G bounds (Sec. VIII).  ￼

⸻

10) SM extension “hand‑waving”

Referee’s concern. No full SM computation; an “uplift factor” u_{\rm SM} is arbitrary.

Response. Correct: the main text uses the scalar baseline only; App. E is transparent about the status and outlines how a full species‑by‑species \beta_{\rm SM} would be assembled (central charges, kernels, thresholds, gauge subtleties), with the optional uplift quoted as a separate systematic if desired. We purposely do not mix this uplift into the baseline results.  ￼

⸻

11) MOND/a_0 connection “forced”

Referee’s concern. The relation a_0=(\Omega_\Lambda^2/2)\,cH_0 looks post‑hoc.

Response. In the revised paper this is presented as an auxiliary corollary of the same Clausius normalization that sets \Omega_\Lambda: in the static, weak‑field, safe‑window regime the \nabla\nabla M^2 terms reorganize the flux into a quasilinear form with an acceleration scale fixed by the FRW zero‑mode matching; dimensional and normalization consistency then yield a_0=(\beta f c_{\rm geo})^2 cH_0/2=(\Omega_\Lambda^2/2)\,cH_0 (App. I; Eq. (24)). We do not use this for any galactic fits in the manuscript.  ￼

⸻

Summary

The revised manuscript does not claim an unconditional proof; it clearly labels (A2) as a working assumption, provides a precise small‑wedge domain (safe window), formalizes the EPMR bridge, derives the constitutive relation and the compensator, and proves continuous‑angle scheme invariance of \beta f c_{\rm geo}. Numerics are reported with conservative systematic uncertainty and with non‑circular checks. Where items are exploratory (SM uplift; the a_0 corollary), we say so explicitly.

We hope this clarifies the scope, rigor, and limitations of the present work.
