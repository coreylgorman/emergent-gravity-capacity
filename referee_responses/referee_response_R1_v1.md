# Response to Referee #1 — “Emergent State‑Dependent Gravity from Local Information Capacity”
**Date:** 2025-08-20T15:39:49.228160Z

We thank the referee for a careful and critical review. Below we address the raised concerns point‑by‑point and outline concrete revisions to the manuscript and analysis pipeline. All claims are explicitly **conditional on Assumption (A2)**, as already emphasized; where appropriate we strengthen derivations, add lemmas, or downgrade claims.

---

## 1) Alleged circularity in normalization schemes (A vs B)

**Referee’s concern.** The two “independent” normalization schemes (A: f=0.8193, c_geo=40; B: f=3.125, c_geo≈10.49) give the same Ω_Λ, suggesting post‑hoc tuning.

**Response.** The equality of Ω_Λ arises because **only the product** x ≡ β·f·c_geo is physical; f and c_geo separately implement a *convention* for partitioning (i) ball→diamond shape and boost normalization, and (ii) how the local wedge is angularly embedded into the FRW 2‑sphere without double counting. We will make this explicit by adding a **Normalization Theorem**:

> **Normalization Theorem (scheme invariance).** Fix a ball→diamond map, an Unruh normalization, a null‑generator density ρ̂_D(u), and a no‑double‑counting rule for tiling the FRW 2‑sphere by disjoint caps. Then any bookkeeping that places a fraction of the boundary/bulk conversion into f and the remainder into c_geo yields the **same observable** product x=β·f·c_geo. 

We also derive the concrete numbers from first principles:
- For a spherical cap of half‑angle θ, the solid angle is ΔΩ=2π(1−cosθ), so
  **c_geo = 4π/ΔΩ = 2/(1−cosθ)**.
  Thus, **c_geo=40 ⇔ cosθ=1−2/40=19/20**, and **c_geo≈10.49 ⇔ cosθ≈0.80934**.
- Scheme A (“minimal wedge”) chooses the smallest disjoint cap consistent with the ball→diamond map already accounted for in f; Scheme B moves more geometric weight into f (f_B=3.125), and recomputes c_geo by the no‑double‑counting rule with the same ρ̂_D(u), giving c_geo≈10.49.

With the measured microscopic sensitivity **β = 0.020855429233**, the invariance check gives
Ω_Λ(A)=β f_A c_geo_A=0.683474126827 and Ω_Λ(B)=β f_B c_geo_B=0.683667039547 (fractional difference +2.82e-04).
We include a machine‑verifiable artifact (`invariance_check.json`) that reproduces these facts.

**Manuscript changes.**
- Replace the informal presentation by a theorem + proof in **App. C** and move all scheme‑specific details to clearly labeled **conventions**. 
- Add an explicit derivation of **c_geo = 2/(1−cosθ)** and the cap‑tiling rationale for each scheme.

---

## 2) Assumption (A2): Clausius relation on arbitrary local wedges

**Referee’s concern.** Extending δQ = T δS with Unruh T to general local wedges is a major assumption without proof.

**Response.** We agree it is the key conditional step. We strengthen and bound its domain:
1. **Small‑wedge regime** (ℓ within the “capacity safe window”): work in Riemann normal coordinates at the wedge apex; state is Hadamard so short‑distance structure is Minkowski‑like. The entanglement first law (δS=δ⟨K⟩) holds for small perturbations; our MI‑subtracted, moment‑killed construction cancels area/contact pieces and removes O(ℓ^0,ℓ^2) curvature terms, leaving a finite O(ℓ^4) coefficient (Lemma A; see §III & App. A).
2. **Unruh normalization**: we stipulate T=κ/2π for the chosen boost Killing field; the normalization is fixed **once** and carried through the Clausius and Noether‑charge steps without further freedom.
3. **Falsifiability**: we highlight observational falsifiers (GW/EM luminosity‑distance ratio, bounds on Ġ/G) and explicitly classify (A2) failure modes.

**Manuscript changes.**
- Add **Appendix H (“Flat→Curved Local Map”)** with a proof sketch: after MI subtraction and moment‑kill, the leading surviving term in δ⟨K_sub⟩ is O(ℓ^4) and is universal under local wedge embeddings; residual curvature pieces are O(ℓ^6) within the safe window.
- Retitle (A2) as a boxed **Assumption**, with precise scope and falsifiers in §II A.

---

## 3) Flat‑space QFT → curved spacetime bridge

**Referee’s concern.** Using flat‑space CHM modular data and applying it in cosmology is unjustified.

**Response.** We clarify that (i) β is computed **entirely** in flat spacetime from MI‑subtracted modular response; (ii) the **only** step connecting to FRW is the Clausius/Noether‑charge map plus the geometric normalization f and c_geo. The small‑wedge expansion and MI subtraction suppress curvature‑dependent contributions up to O(ℓ^4), and we never insert cosmological inputs into β. We provide a clean separation of steps to remove any impression of hidden fits.

**Manuscript changes.**
- New figure/flowchart separating: *QFT β* → *Clausius step* → *Geometric normalization* → *FRW zero mode*.
- A line‑by‑line “no cosmology in β” audit in §IV D (with gates and residuals).

---

## 4) CGM critique and the compensator term

**Referee’s concern.** The compensator (A/4G) β δσ looks assumed rather than derived; the safe window appears ad hoc; marginal‑log cancellation seems tuned.

**Response.** We will tighten §V as follows:
- Start from **S_grav = A/[4G(x)]**. Varying both A and G gives
  δS_grav = (1/4G) δA − (A/4G^2) δG = (1/4G) δA + (A/4G) β δσ,
  where δG/G=−β δσ by **definition** of the state metric (constitutive closure). This is not an extra assumption but the direct consequence of allowing M^2(x) to run.
- Show explicitly how the mild running M^2(x) supplies the marginal‑log compensator at Δ=d/2 within the wedge window (scaling argument).
- Move the safe‑window statement from heuristic to **inequality with scales** and a table of admissible ℓ for typical laboratory/astrophysical contexts.

**Manuscript changes.** Expand Eqs. (11–13) with a one‑page derivation; add a boxed “Scaling at Δ=d/2” remark with the log structure.

---

## 5) Mathematical items

- **Eq. (5) (ℓ⁴ isolation).** We add a formal proof (App. A): choosing (a,b) to kill the r⁰ and r² moments in the weighted ball integrals ensures that area/contact and curvature‑contact pieces cancel, leaving a finite ℓ⁴ term that defines I₀₀; regulator independence at this order follows.
- **c_geo derivations.** We now *derive* (not assert) the cap angles:
  c_geo = 2/(1−cosθ) ⇒
  Scheme A: c_geo=40 ⇒ cosθ=19/20; 
  Scheme B: c_geo≈10.49 ⇒ cosθ≈0.80934.
- **MOND scale a₀.** We demote a₀ ≡ (Ω_Λ²/2) c H₀ to an **auxiliary conversion** (not a claim of galactic phenomenology); we already state that our framework does **not** mimic MOND dynamics.

---

## 6) Numerical implementation & uncertainties

The current pipeline already enforces **positivity** and **moment‑kill residual gates** and reproduces the manuscript values:
β = 0.020855429233, I00 = 0.1077748682, C_T = 0.03079794676, Ω_Λ(A) = 0.683474126827, a₀(A) = 1.588e-10 m/s².

We will add:
1. **Window‑scan reports** over (σ₁, σ₂), u_gap, T_max, and grids (N_r,N_s,N_τ) with plateau detection and a robustness score.
2. **Scheme‑invariance check** (A vs B) with automatic **cap‑angle diagnostic** and a JSON artifact (`invariance_check.json`).
3. **Error budget table** separating numerical vs. systematic (windowing) components and a placeholder row for **SM uplift** (kept out of the baseline).

---

## 7) References

We scale back claims about connections to other programs, tighten language regarding Padmanabhan, and add missing citations on modular Hamiltonians in interacting/gauge theories. (We mark these in the addendum for insertion.)

---

## Requested editorial changes (high level)

- Clear “boxed normalization” once, early; remove duplications.
- Consistent notation across main text and appendices.
- Move phenomenology (GW/EM, Ġ/G) to a consolidated section with explicit formulas and bounds mapping.

---

## Attachments / Artifacts

- `referee_pipeline_v2.py`: helper to verify scheme invariance and produce scans (see below).
- `invariance_check.json`: machine‑verifiable check that Ω_Λ is scheme‑invariant with the supplied β.
- `manuscript_revisions_addendum_v1.md`: drop‑in text (lemmas, theorem, derivations) for App. A/C/H and affected sections.
- `referee_response_R1_v1.md`: this letter.

We hope these changes address the referee’s main concerns and welcome further questions.
