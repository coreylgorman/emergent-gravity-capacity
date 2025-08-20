
# Manuscript Revisions — Addendum v2 (for insertion)

## A. Small‑Diamond Response Lemma (formal version)
**Lemma.** Let the state be Hadamard in a neighborhood of point \(x\). In Riemann normal coordinates about \(x\),
for diamond size \(\ell\) satisfying \(\epsilon_{\rm UV}\ll\ell\ll L_{\rm curv}\), the MI‑subtracted, moment‑killed modular response obeys
\[\delta\!\langle K_{\rm sub}(\ell)\rangle = (2\pi\,C_T I_{00})\,\ell^4\,\delta\sigma(x) + \mathcal O(\ell^6\,\mathcal R)\, ,\]
where \(\mathcal R\) denotes curvature invariants with appropriate dimension. The coefficients of \(\mathcal O(\ell^0)\) and \(\mathcal O(\ell^2)\)
(including curvature‑contact terms) cancel by construction.

*Sketch.* Expand correlators and geometric factors in the RNC ball. The subtraction \(K(\ell)-aK(\sigma_1\ell)-bK(\sigma_2\ell)\)
annihilates the first two radial moments of any smooth integrand, eliminating contact structures and the curvature‑contact pieces
that share their moment profile. The first nonvanishing isotropic term is \(\propto \ell^4\).

## B. Safe‑window inequality and remainder control
Choose \(\ell\) such that
\[\epsilon_{\rm UV}\ll \ell \ll \min\{L_{\rm curv},\,\lambda_{\rm mfp},\,m_i^{-1}\}\,.\]
Then the remainder bound becomes \(|\mathcal O(\ell^6 \mathcal R)|\lesssim C\,(\ell/L_{\rm curv})^2\,\ell^4\,|\delta\sigma|\) for some \(C=O(1)\),
ensuring dominance of the \(\ell^4\) term. (This section is purely QFT/geometry; no Clausius input is used.)

## C. Geometric multiplicity \(c_{\rm geo}\) via the Noether‑Flux Equipartition Principle (NFEP)
**Definition (NFEP).** Fix the Unruh‑normalized \(\chi^a\) and generator weight \(\hat\rho_D(u)\) used in the modular calculation. Let
\[\Phi_\wedge(\theta_\star)\equiv\int_{\rm wedge(\theta_\star)}\frac{\delta Q}{T}\,,\qquad \Phi_{\mathbb S^2}\equiv\int_{\mathbb S^2}\frac{\delta Q}{T}\,.\]
Then \(c_{\rm geo}\) is the unique multiplicity such that disjoint spherical caps of half‑angle \(\theta_\star\) **flux‑equipartition** the sphere:
\[c_{\rm geo}\,\Phi_\wedge(\theta_\star)=\Phi_{\mathbb S^2}\,.\]
This prevents angular double counting by construction and fixes \(\theta_\star\) (hence \(c_{\rm geo}=2/(1-\cos\theta_\star)\)) once \(\hat\rho_D\) is fixed.
All prior numerical instances (e.g. \(c_{\rm geo}=40\), \(10.49\)) are now clearly labelled **historical conventions** and moved to examples.

## D. Hypothesis \(H_{\sigma\to G}\) and variational motivation
We elevate the constitutive relation to a hypothesis:
\[\boxed{~H_{\sigma\to G}:\quad \delta G/G = -\beta\,\delta\sigma(x).~}\]
Motivation: maximize the total wedge entropy \(S_{\rm tot}=S_{\rm mat}+S_{\rm grav}\) at fixed Clausius flux and Bianchi consistency.
With \(S_{\rm grav}=A/[4G(x)]\), \(\delta S_{\rm grav}=(1/4G)\delta A - (A/4G)\,\delta G/G\). Extremality with respect to local variations of
\(\delta\sigma\) delivers a linear closure between \(\delta G/G\) and \(\delta\sigma\) with proportionality \(\beta\). This anchors the compensator term
without postulating it independently.

## E. Theorem/Assumption labeling
Results are tagged as **Theorem** (purely QFT/geometric), **NFEP** (definition), or **Assumption** (A2), **Hypothesis** (H\(_{\sigma\to G}\)).
Every formula in the main text is annotated accordingly; all global claims depending on A2 or H\(_{\sigma\to G}\) are flagged **conditional**.

## F. Error‑budget demarcation
- **Tier‑N (Numerical/QFT):** affects \(\beta\); quoted as \(\approx 3\%\) from stability scans and gates.
- **Tier‑M (Model):** arises from A2, H\(_{\sigma\to G}\), and NFEP; quoted qualitatively as **dominant** for claims about \(\Omega_\Lambda\).
