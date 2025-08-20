Manuscript: Emergent State-Dependent Gravity from Local Information Capacity: A Conditional Thermodynamic Derivation with Cosmological Predictions (RevTeX, dated August 20, 2025)
Thank you for the careful, critical review and the follow-up report. We appreciate the opportunity to clarify the logic and to improve the rigor and presentation. Below we respond point-by-point and list the precise changes made in the revised manuscript.

⸻

Executive Summary of Revisions
	1.	Cap angles and “tuning” → Removed as inputs. We now prove a continuous-angle scheme invariance: with a unit–solid–angle boundary normalization, the observable product \beta\,f(\theta)\,c_{\rm geo}(\theta) is independent of the cap half-angle \theta. The earlier angles (e.g., \cos\theta=19/20, \cos\theta\simeq 0.80934) are worked examples only. See Sec. VI.C; Appendix C gives the spherical-cap form c_{\rm geo}(\theta)=2/(1-\cos\theta).
	2.	Assumption (A2): Clausius for small, non-stationary wedges → Conditional theorems with explicit bounds. We added Lemma H.1 (first‑law domain), Lemma H.2 (moment‑kill + MI subtraction), and Proposition H.1 (curvature suppression to O(\ell^6) and EPMR: flat‑space modular kernel valid to O(\ell^4)). These statements make the scope of (A2) precise and testable rather than asserted. See Appendix H.
	3.	Flat‑QFT \rightarrow gravity bridge → Formalized as EPMR. We show that the \ell^4 coefficient defining \beta is a flat‑space object (Sec. III.A, Eq. (5)), while curvature enters only at O(\ell^6) after MI subtraction and moment‑kill (Appendix H). Gravity appears only after the Clausius/Noether mapping and geometric normalization.
	4.	CGM critique and constitutive relation → Derived, not posited. We added a variational capacity closure: extremizing a Wald‑like functional with a local capacity constraint yields \delta G/G=-\beta\,\delta \sigma. Substituting into \delta S_{\rm grav}=A(4G)^{-1} produces the compensator term (A/4G)\beta\,\delta\sigma that cures the marginal‑log obstruction in the stated window. See Sec. V.B (Eqs. (15)–(16)) and Eqs. (11)–(13).
	5.	Operational definition of “capacity” → No bare entropy needed. We define the state metric \sigma(x) as the finite \ell^4 coefficient of the MI‑subtracted, moment‑killed modular response (Sec. III.A, Eq. (5)). This is regulator‑independent at this order and directly measurable in principle via modular response—no standalone finite regional entropy is assumed.
	6.	Normalization procedures and circularity → Invariant product. Only the product \beta f c_{\rm geo} enters the observable \Omega_\Lambda. The revised Sec. VI.C shows \theta-invariance, removing any appearance of post‑hoc angle selection. Sec. VII.B retains a non‑circularity check: varying only \beta moves \Omega_\Lambda linearly.
	7.	Weak‑field scale a_0 → Derived from the same invariants. Appendix I derives
a_0=\frac{\Omega_\Lambda^2}{2}\,c\,H_0=\frac{(\beta f c_{\rm geo})^2}{2}\,c\,H_0,
as an auxiliary consequence of the same Clausius normalization (static, weak‑field, safe‑window regime). This is not used for galaxy fits here.
	8.	Uncertainties → Two‑tier accounting. We separate (i) numerical/systematic uncertainty on \beta (\sim 3\%) from (ii) conceptual scope set by (A2). The abstract and conclusions explicitly state the conditional status of all quantitative claims.

⸻

Point‑by‑Point Responses

1) “Circular reasoning—cap angles appear tuned”

Response. We prove that cap angles are not inputs: defining f_{\rm bdy}(\theta)=f_{\rm bdy}^{\rm unit}\,\Delta\Omega(\theta) with \Delta\Omega(\theta)=2\pi(1-\cos\theta) and c_{\rm geo}(\theta)=4\pi/\Delta\Omega(\theta)=2/(1-\cos\theta), the observable product
\beta f(\theta)c_{\rm geo}(\theta)=\beta\,f_{\rm shape}f_{\rm boost}f_{\rm cont}f_{\rm bdy}^{\rm unit}(4\pi)
is independent of \theta (Sec. VI.C, Eqs. (18)–(20)). The earlier angles in Schemes A/B are worked examples illustrating the same invariant. Appendix C gives the spherical‑cap form explicitly (Eq. (C2)).

2) “(A2) insufficiently justified; curvature order; first‑law domain”

Response. Appendix H upgrades the status of (A2) to conditional theorems with explicit bounds:
	•	Lemma H.1: first‑law domain for small diamonds in Hadamard states.
	•	Lemma H.2: MI subtraction + moment‑kill cancel contact and curvature–contact terms through O(\ell^2).
	•	Proposition H.1: surviving isotropic term is O(\ell^4) with the flat‑space coefficient; curvature enters at O(\ell^6).
We also collect safe‑window inequalities in Sec. III.B. The manuscript continues to label all results as conditional on (A2).

3) “Flat‑space modular Hamiltonians used in gravity without proof”

Response. Sec. III.A defines \sigma(x) via Eq. (5) (MI‑subtracted, moment‑killed response with leading \ell^4 scaling). EPMR (Appendix H) shows the \ell^4 coefficient is the flat‑space one, with curvature sensitivity at O(\ell^6). Thus \beta is computed purely in flat space; gravity appears only after the Clausius/Noether step and geometric normalization (Secs. VI–VII).

4) “CGM not resolved; \delta G/G=-\beta\delta\sigma assumed”

Response. We derive the constitutive relation via a variational capacity constraint (Sec. V.B, Eqs. (15)–(16)), replacing a bare ansatz. Varying S_{\rm grav}=A/[4G(x)] with running G(x) gives the compensator (A/4G)\beta\,\delta\sigma (Eqs. (11)–(13)), which cancels the marginal‑log obstruction within the safe window—precisely where CGM flagged issues.

5) “Information capacity ill‑defined; UV‑sensitive”

Response. We do not rely on a standalone finite entropy. The state metric is defined operationally by the finite \ell^4 coefficient in \delta\!\langle K_{\rm sub}\rangle/(2\pi C_T I_{00}) (Sec. III.A, Eq. (5)), after MI subtraction and moment‑kill. This removes UV/area and curvature–contact contributions by construction (Appendix A).

6) “Eq. (5) assumed; scaling needs proof”

Response. Appendix A (moment‑kill identities) explains how the chosen linear combination of balls cancels the r^0 and r^2 moments, isolating \ell^4. The curvature‑suppression order and the flat‑space nature of the \ell^4 coefficient are formalized in Appendix H (Proposition H.1).

7) “Normalization f, c_{\rm geo} arbitrary; circular construction”

Response. Only the product \beta f c_{\rm geo} is physical (Sec. VII). The new continuous‑\theta section (Sec. VI.C) shows this product is angle‑independent, eliminating any hidden tuning. Sec. VII.B includes a non‑circularity check: varying \beta alone translates one‑to‑one into \Omega_\Lambda.

8) “Physical mechanism; local‑to‑global connection”

Response. The mechanism is stated up front (Sec. I): finite capacity \Rightarrow minimal four‑geometry response \Rightarrow local time dilation \Rightarrow global emergent gravity. The FRW zero‑mode mapping is clean: \Omega_\Lambda=\beta f c_{\rm geo} (Sec. VII, Eq. (23)) with consistency under the Bianchi identity (Eq. (27)); the field equation (Eq. (28)) shows how running M^2 enters at background/linear order.

9) “Weak falsifiability; generic tests”

Response. We emphasize frame‑invariant signatures: d_L^{\rm GW}/d_L^{\rm EM} depends only on the integral of \alpha_M (Sec. VIII.1), and bounds on \dot G/G map directly to \alpha_M(0) (Sec. VIII.2). These tests directly falsify the (A2)‑based Clausius closure and the capacity mechanism.

10) “a_0 appears reverse‑engineered”

Response. Appendix I derives
a_0=\tfrac{\Omega_\Lambda^2}{2}\,c\,H_0
from the same invariant that fixes \Omega_\Lambda, using the static weak‑field limit of the constitutive equations and the identical Clausius normalization. We present a_0 as an auxiliary consequence (not used for any fits here).

⸻

Specific Changes in the Manuscript
	•	Sec. VI.C (Continuous‑angle normalization and scheme invariance): new derivation showing \beta f(\theta)c_{\rm geo}(\theta) is \theta-independent; cap angles demoted to examples.
	•	Sec. V.B (Variational capacity closure): derivation of \delta G/G=-\beta\,\delta\sigma and compensator structure (A/4G)\beta\,\delta\sigma.
	•	Appendix H (Small‑wedge Clausius domain and curvature suppression / EPMR): Lemma H.1, Lemma H.2, Proposition H.1.
	•	Appendix I (Weak‑field flux law; derivation of a_0=(\Omega_\Lambda^2/2)cH_0): normalization matched to the same Clausius invariant.
	•	Sec. VII (Cosmological constant sector): retained clean mapping \Omega_\Lambda=\beta f c_{\rm geo}; added explicit non‑circularity check.
	•	Sec. III (State metric and safe window): clarified the operational definition and the bounds.

All other content (numerical baselines, replication preset, and species bookkeeping) remains unchanged.

⸻

Closing

We hope these changes address the referee’s concerns: the cap‑angle issue is removed by construction; the Clausius step is made conditional with explicit domain and curvature order; the flat‑\beta \rightarrow gravity bridge is formalized (EPMR); the constitutive closure is derived variationally; and the weak‑field a_0 scale follows as a parameter‑free corollary. We keep the manuscript explicit about its conditional nature: if (A2) fails in the stated window, the framework is falsified or must be revised.

We thank the referee for the rigorous feedback and would be happy to consider further suggestions for clarity.