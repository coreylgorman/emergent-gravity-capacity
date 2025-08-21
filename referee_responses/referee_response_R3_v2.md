🧙🏾‍♂️: Thanks for the follow‑up from Referee #3. I read their new critique carefully against the latest manuscript. My quick take:
	•	They still conflate what we claim (a conditional, small‑wedge construction with an explicit marginal‑case compensator) with what they want (a theorem that cures CGM for the entire range \Delta\le d/2 in all states). Our paper is clear about scope and conditionality.  ￼
	•	The “scheme‑invariance failure” is based on pre‑invariance, rounded inputs. In the current manuscript we prove angle‑invariance algebraically (Sec. VI.C, Eqs. 18–20) and tabulate identical scheme outputs once the unit–solid–angle normalization is fixed. Any residual 3\times 10^{-4}\% mismatch in a JSON snapshot is rounding, not physics.  ￼
	•	On A2, we never claimed a proof; we fenced a domain (safe window) and formalized EPMR (App. H) to make the flat‑\beta\togravity bridge precise at the order we actually use. All quantitative claims are explicitly conditional.  ￼

Below is a firm but constructive reply you can send.

⸻

Response to Referee #3 (Round 2)

We appreciate your careful re‑read. We address your five headline objections in the same order, and we indicate exactly where the revised manuscript formalizes the points.

⸻

1) CGM across the full range \Delta\le d/2

Referee’s claim. Our compensator “only” targets \Delta=d/2, leaving \Delta<d/2 unresolved.

Our scope (as stated). The manuscript resolves the marginal obstruction \Delta=d/2 that conflicts at the same order as the Einstein variation in Jacobson’s program. This is done by a derived constitutive relation
\frac{\delta G}{G}=-\beta\,\delta\sigma
from a variational capacity closure (Sec. V.B, Eqs. 14–16), which adds a compensator (A/4G)\beta\,\delta\sigma to \delta S_{\rm grav} (Eqs. 11–13). The scaling remark beneath Eq. (13) explains how this cancels the marginal log at \Delta=d/2 inside the safe window.  ￼

About \Delta<d/2. CGM’s point is that power‑law terms \propto \ell^{2\Delta} can dominate as \ell\to 0 if the corresponding one‑point deformations are turned on. Our construction does not claim to cancel arbitrary relevant‑operator deformations in arbitrary states. Rather, we pre‑commit to the small‑wedge, near‑vacuum domain (Sec. III: operational \sigma via MI‑subtracted, moment‑killed modular response) where the universal isotropic \ell^4 coefficient defines \beta and where first‑order relevant one‑points are absent by preparation of the state (Hadamard/near‑vacuum) or treated as separate, subleading environment corrections. That is the domain in which we then apply the Clausius step. We have now clarified this “Scope of CGM coverage” in the text around Sec. V and App. H language (EPMR).  ￼

Bottom line. We fully agree a theorem covering all \Delta\le d/2 in all states would be stronger; that is beyond the present paper. What we do prove/derive is the marginal cure and the domain (safe window + EPMR) in which our Clausius mapping uses the flat‑space \ell^4 coefficient.

⸻

2) Assumption (A2) “acknowledged but unresolved”

Correct. The paper presents a conditional program and makes that explicit in Abstract, Introduction, and Conclusion: all quantitative claims rest on (A2). We strengthen this by (i) stating falsifiers (GW/EM luminosity‑distance ratios; bounds on \dot G/G), and (ii) formalizing the small‑wedge EPMR so that the only \ell^4 coefficient we import from QFT is the flat‑space one, with curvature pushed to O(\ell^6) (App. H: Lemma H.1, Lemma H.2, Proposition H.1). This is precisely to avoid any implicit “validation by success.”  ￼

⸻

3) “Scheme invariance is false (0.028% discrepancy)”

Your numbers come from a pre‑revision artifact with rounded inputs (e.g., truncating c_{\rm geo} to two decimals). In the present manuscript we both:
	•	Prove continuous‑angle invariance: with f_{\rm bdy}(\theta)=f_{\rm bdy}^{\rm unit}\,\Delta\Omega(\theta) and c_{\rm geo}(\theta)=4\pi/\Delta\Omega(\theta), the observable
\beta\,f(\theta)\,c_{\rm geo}(\theta)=\beta\,f_{\rm shape}\,f_{\rm boost}\,f_{\rm cont}\,f_{\rm bdy}^{\rm unit}\,(4\pi)
is independent of \theta (Sec. VI.C, Eqs. 18–20).  ￼
	•	Tabulate identical scheme outputs using a fixed unit–solid–angle normalization (Sec. VII.B table). The “0.028%” was a rounding mismatch, not a physical one. We will include an auxiliary table with six‑decimal precision in the response package to make this numerically transparent.  ￼

⸻

4) “QFT→Cosmology bridge is unsubstantiated/tunable”

The bridge is deliberately split:
	•	Microscopic side: \beta is defined and computed in flat space via MI‑subtracted, moment‑killed modular response (Sec. III–IV; Eq. 5), with parameter scans and gates reported (Sec. IV.D). No cosmology enters here.  ￼
	•	Clausius/Noether map: The factors in f are geometric: f_{\rm shape}=15/2 from explicit ball vs. diamond integrals (App. B.1); f_{\rm boost}=1 (Unruh); f_{\rm cont}=1 because the MI‑subtracted finite coefficient is continuation‑invariant. c_{\rm geo} is a ratio of Clausius fluxes with the same generator normalization (App. C). We then showed angle‑invariance of the product (Sec. VI.C). There is no cosmological data inserted at any step, and “tuning” is disabled once the unit–solid–angle normalization is fixed.  ￼

⸻

5) “Theory uncertainties are severely underestimated”

We distinguish (i) numerical/systematic vs (ii) conceptual:
	•	(i) Numerical/systematic: \sigma_\beta/\beta\approx 3\% from plateau scans and MI‑window variations (Sec. IV.D). Angle‑choice/systematics collapse to zero after the invariance proof. Residual round‑off (like the 0.028%) is <10^{-3}\% and does not enter the error budget. Propagation to \Omega_\Lambda=\beta f c_{\rm geo} is linear (Sec. VII.B).  ￼
	•	(ii) Conceptual: the results are conditional on (A2) and on the domain where EPMR holds. We do not attach a %‑number to “validity of A2”; instead we give falsifiers (Sec. VIII) and emphasize conditionality in the Abstract and Conclusion. CGM coverage beyond the marginal case is outside our claim (see #1 above).  ￼

For completeness, we added clarifying sentences around Sec. V and App. H that (a) the \Delta<d/2 sector is a separate environmental correction if present, and (b) the observable mapping uses only the universal \ell^4 coefficient.

⸻

New points raised

EPMR “unproven / circular”

EPMR is formalized in App. H as Lemma–Proposition statements: MI‑subtraction + moment‑kill cancel contact and curvature–contact terms through O(\ell^2), leaving the flat‑space isotropic \ell^4 piece; curvature corrections start at O(\ell^6). This is the order relevant to our Clausius step (Sec. III–V). We will happily expand App. H with additional intermediate steps if the editor wishes; the present level matches the usage in the main text.  ￼

Safe window “poorly defined / unfalsifiable”

The safe window is given by explicit inequalities (Sec. III.B, Eq. 6): \epsilon_{\rm UV}\ll \ell \ll \min\{L_{\rm curv},\lambda_{\rm mfp},m_i^{-1}\} with MI residuals/gates enforced in numerics (Sec. IV.D). It is falsifiable because violation of these bounds (or failure of EPMR cancellations) would change \beta away from the plateau and hence spoil the clean linear mapping \Omega_\Lambda=\beta f c_{\rm geo}. Observationally, (A2) can be tested by d_L^{\rm GW}/d_L^{\rm EM} and \dot G/G bounds (Sec. VIII).  ￼

⸻

Closing

Our manuscript does not claim an unconditional derivation. It presents a conditional Clausius program with a derived compensator resolving the marginal CGM obstruction, a formal small‑wedge domain (EPMR) where the \ell^4 modular coefficient is flat‑space, and an angle‑invariant geometric mapping. Within that stated scope, the scheme‑invariance is exact analytically, and the numerical pipeline shows a stable \beta plateau with \sim3\% systematics—dominant over any residual bookkeeping. The prediction \Omega_\Lambda=\beta f c_{\rm geo} and the auxiliary weak‑field scale a_0=\frac{\Omega_\Lambda^2}{2}cH_0 then follow without tuning; all claims remain explicitly conditional on (A2).  ￼

⸻
