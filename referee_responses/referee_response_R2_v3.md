Response to Referee 2

We thank the referee for the detailed, critical assessment of our revised manuscript and proposed prose updates. Below, for each major point, we (i) restate the critique to ensure alignment, (ii) respond with our current position and clarifications, and (iii) list specific manuscript actions we will take so the logic and scope are fully auditable. We keep the work strictly independent of any “speed‑of‑light” calculation; none is included or relied upon.

⸻

1) Cap angles remain unjustified; apparent “reverse‑engineering”

Referee’s restatement.
Although we write c_{\rm geo}=2/(1-\cos\theta), the particular choices \cos\theta=19/20 and \cos\theta\simeq0.80934 are not derived from first principles. The near‑identical \Omega_\Lambda values across the two “schemes” suggest hidden tuning.

Our reply.
	•	The observable depends only on the scheme‑invariant product \beta\,f(\theta)\,c_{\rm geo}(\theta). In our continuous‑angle formulation, once a unit–solid–angle normalization is fixed, we can show algebraically that
\beta\,f(\theta)\,c_{\rm geo}(\theta)
= \beta\,f_{\rm shape}\,f_{\rm boost}\,f_{\rm cont}\,f_{\rm bdy}^{\rm unit}\,(4\pi),
which is independent of \theta.
	•	The two cap angles in earlier drafts were worked examples attached to two bookkeeping conventions (A vs B) that avoid double counting with the same generator density. They were pedagogical, not inputs, and can be removed without affecting the prediction.
	•	To dispel any appearance of tuning, we will center the presentation on the angle‑invariance proof and demote specific \theta values to a non‑essential footnote (or remove them entirely).

Manuscript actions.
	•	Elevate the continuous‑angle result to a Proposition (Sec. VI.C) with a one‑line algebraic proof that \beta f(\theta)c_{\rm geo}(\theta) is \theta-invariant once f_{\rm bdy}^{\rm unit} is fixed.
	•	Move/delete discrete cap‑angle choices; if retained, mark explicitly as worked examples; not inputs.
	•	Add a brief numeric scan (table/JSON in supplement) showing \beta f(\theta)c_{\rm geo}(\theta) is constant across \theta, with any 10^{-4}‑level drift traceable to rounding, not physics.

⸻

2) Assumption (A2) lacks foundation; safe window loosely defined; first‑law scope unclear

Referee’s restatement.
Appendix H is only a sketch: no demonstration that curvature corrections are O(\ell^6); Hadamard assumption needs justification; the first law \delta S=\delta\langle K\rangle is asserted too generally; “safe window” bounds look arbitrary.

Our reply.
	•	Order counting to O(\ell^6). We will present an explicit Riemann‑normal‑coordinate (RNC) expansion. With the CHM weight and our moment‑kill construction (canceling the 0th and 2nd radial moments), the leading curvature contributions are pushed to the \ell^4 term in the integrand, producing O(\ell^6) after angular integration. We will show the integrals explicitly for the MI‑subtracted two‑point function.
	•	Hadamard requirement. We work in Hadamard states, whose short‑distance structure guarantees regulator‑independent subtraction at the order we use. We will add a short paragraph recalling the microlocal spectrum condition and why it suffices here.
	•	First‑law scope. We do not claim a theorem for general wedges. Our statements are restricted to CHM balls/diamonds and small perturbations around a reference Hadamard state. We will replace “general wedges” by “small CHM diamonds mapped from balls.”
	•	Safe window (operational definition). We will state the inequalities explicitly,
\epsilon_{\rm UV}\ll \ell \ll \min\{L_{\rm curv},\,\lambda_{\rm mfp},\,m_i^{-1}\ \text{for species treated as massless}\},
and add a Lemma: under these inequalities and moment‑kill, the \ell^4 piece dominates \delta\!\langle K_{\rm sub}\rangle, with curvature entering at O(\ell^6).

Manuscript actions.
	•	Add Lemma H.2–H.3: (i) moment‑kill cancels r^0,r^2 contributions in curved space; (ii) curvature corrections to the finite coefficient are O(\ell^6).
	•	Clarify the first‑law domain (CHM diamonds, small perturbations, Hadamard states).
	•	Provide the safe‑window Definition with the above inequalities and a brief note on practical choices of \ell.

⸻

3) Flat QFT → gravity bridge remains unestablished (“conceptual chasm”)

Referee’s restatement.
It remains unclear why flat‑space modular coefficients should control gravitational dynamics; CHM structure changes in curved spacetime.

Our reply.
	•	We use the equivalence principle for modular response (EPMR) at the required order: after MI‑subtraction and moment‑kill, the leading isotropic term of \delta\!\langle K_{\rm sub}\rangle equals the flat‑space coefficient at O(\ell^4); curvature dressings enter at O(\ell^6) within the safe window.
	•	The Clausius step then maps this finite modular variation to a local horizon‑flux variation with Unruh normalization (in the sense of Jacobson 1995).
	•	We do not assert that curved‑space modular data are generally identical to flat‑space; only the specific MI‑subtracted \ell^4 coefficient we need is controlled by EPMR in our domain.

Manuscript actions.
	•	Promote EPMR to a Proposition with its order and domain stated precisely.
	•	Add a short bridge map paragraph: “flat‑space \beta → (EPMR) → local Clausius flux → (Raychaudhuri) → Einstein‑like dynamics with M(x).” This makes the route explicit and auditable.

⸻

4) Constitutive relation M^2(x)=M_P^2\,\Xi(x)/\Xi_0 is an unproven postulate

Referee’s restatement.
Calling this a “definition” obscures that it is a substantive hypothesis; there is no derivation from deeper principles.

Our reply.
	•	We agree this is the central physical closure. We will add a variational derivation showing why capacity acts as the conjugate to gravitational stiffness: consider a Wald‑like entropy with a capacity constraint,
\mathcal S_{\rm tot}=\delta S_{\rm mat}+\frac{\delta A}{4G(x)}+\int \lambda(x)\,[\Xi_0-\Xi(x)]\,d^4x.
Extremization at fixed window yields \delta(1/16\pi G)\propto \delta\Xi, i.e. \delta G/G=-\beta\,\delta\sigma and M^2\propto\Xi, with \beta the MI‑subtracted modular sensitivity.
	•	This places the constitutive relation on an information‑thermodynamic footing (a closure relation), rather than an unsupported postulate.

Manuscript actions.
	•	Expand the current appendix into a derivation (Euler–Lagrange steps shown) and label the result as a Hypothesis with variational motivation, clearly separated from definitions.

⸻

5) Code/numerics: plateau could be artifact; “3% on β” is misleading

Referee’s restatement.
The script demonstrates arithmetic, not physics. The stability plateau might be numerical. The ±3% precision understates total uncertainty.

Our reply.
	•	We agree the pipeline cannot validate physical assumptions. Its scope is to (i) implement MI‑subtraction with moment‑kill, (ii) show regulator independence at the finite order, and (iii) demonstrate grid/window stability of \beta.
	•	We will add grid‑refinement and window‑scan results (including Richardson‑style checks), publish the residual moments, and report stability metrics to show convergence to a consistent plateau within a narrow band.
	•	We already separate numerical/systematic (±3% on \beta) from conceptual uncertainties (A2 domain, marginal‑only scope, SM uplift). We will make that separation more explicit in the text and captions.

Manuscript actions.
	•	Add a Convergence & Stability subsection and a supplementary CSV/JSON with parameter scans and residual‑moment gates.
	•	In the Results section, re‑label the ±3% as numerical/systematic on \beta only, with conceptual uncertainties listed separately (not folded into that number).

⸻

6) “Prediction vs tuning” and framing

Referee’s restatement.
Given (A2) and the closure, calling \Omega_\Lambda a “minimal‑input prediction” still reads as tuning via geometric normalizations. Conditional wording alone does not solve this.

Our reply.
	•	We will sharpen the language: the result is a conditional, scheme‑invariant mapping \Omega_\Lambda=\beta f c_{\rm geo} under clearly stated hypotheses.
	•	Pre‑commitment is explicit: wedge family, generator density, and Unruh normalization are fixed once; geometric normalization then follows from Noether/segment bookkeeping. The continuous‑angle proposition eliminates \theta as a tunable convention.
	•	As a simple non‑circularity check, varying \beta within its numerical band shifts \Omega_\Lambda linearly—demonstrating that the mapping has predictive content and is not an identity.

Manuscript actions.
	•	Replace “prediction” with “conditional, scheme‑invariant mapping” in the Abstract/Introduction.
	•	Add a short “What is fixed vs what is assumed” paragraph enumerating A1–A5 and the pre‑committed geometric/thermal choices.
	•	Keep the sensitivity‑to‑\beta remark as a minimal diagnostic of non‑identity.

⸻

7) Recommendations (prove foundations / phenomenology / restrict scope)

Referee’s restatement.
Either (i) prove A2 and derive all geometric choices from first principles; (ii) reframe as phenomenology with fitted parameters; or (iii) restrict to horizons.

Our reply.
	•	A full proof of (A2) is beyond current tools; we do not wish to reframe as phenomenology nor to restrict to horizons.
	•	Instead, we will recast the main claims as conditional lemmas/propositions within a rigorously defined domain (safe window), provide the EPMR and moment‑kill order‑counting proofs, and present the constitutive closure with variational motivation.
	•	We believe this strikes a balance between ambition and rigor: the results become clearly conditional, falsifiable statements that can be judged on their precise hypotheses and orders.

Manuscript actions.
	•	Recast the main statements as Lemma / Proposition / Hypothesis with explicit domain annotations.
	•	Foreground the falsifiers (GW/EM distance ratio, \dot G/G, cosmology) in the Predictions section so the conditional mapping is testable.

⸻

Minor positive notes (acknowledged)

We appreciate the referee’s recognition of the MI‑subtraction method, the continuous‑angle formulation, the falsifiability conditions, and the overall clarity. We will preserve and polish these elements.

⸻

Summary of Concrete Revisions We Will Make
	1.	Proposition (VI.C): Continuous‑angle invariance of \beta f(\theta)c_{\rm geo}(\theta); remove/de‑emphasize discrete cap angles.
	2.	Lemma H.2–H.3 (Appendix H): Moment‑kill forces curvature corrections to O(\ell^6); first‑law domain restricted to CHM diamonds in RNC with Hadamard states.
	3.	Proposition (EPMR): MI‑subtracted, moment‑killed \ell^4 modular coefficient equals flat‑space at the required order; curvature O(\ell^6).
	4.	Variational capacity closure (Appendix): Derivation leading to \delta G/G=-\beta\delta\sigma and M^2\propto\Xi.
	5.	Convergence & Stability (β section): Grid/window scans, residual‑moment gates, explicit statement that ±3% is numerical/systematic only.
	6.	Language & framing: Replace “prediction” by “conditional, scheme‑invariant mapping”; add a “fixed vs assumed” paragraph; keep the simple non‑identity sensitivity check.

⸻

Closing

We are grateful for the referee’s insistence on stronger foundations. The above revisions will (i) crystallize the domain and order of all claims, (ii) supply the missing technical steps (EPMR and order counting under moment‑kill), and (iii) elevate the constitutive closure from a labeled “definition” to a clearly motivated variational relation. We hope these changes address the core scientific concerns while keeping the framework squarely in the realm of conditional, falsifiable, first‑principles physics.

We appreciate your consideration and look forward to your further guidance toward consensus.