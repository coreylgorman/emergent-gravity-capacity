# Response to Referee #2 — Round 2
**Date:** 2025-08-20T16:03:04.021636Z

We appreciate the forthright critique. Below we **tighten the scientific scope**, replace ambiguous conventions by
**single, principled definitions**, and clearly separate *conditional* statements from *derived* statements.

## Executive adjustments (high level)
1. **Cap angles / c_geo** — We drop the dual‑scheme presentation as a “prediction.” We now adopt a **single
first‑principles definition** of the geometric multiplicity based on a *Noether‑Flux Equipartition Principle (NFEP)*
(see Addendum v2, §C). This replaces numerical choices like cosθ=19/20. All appearance of specific values
(c_geo=40, 10.49) in the main text are removed and moved to an **illustrative appendix** only.
2. **Assumption (A2)** — We sharpen its status and add a **controlled small‑diamond error bound**:
\[\delta\!\langle K_\text{{sub}}(\ell)\rangle = (2\pi C_T I_{{00}})\,\ell^4\,\delta\sigma(x) + \mathcal O(\ell^6 \, \mathcal R\, \Lambda_\text{{UV}}^0)\, ,\]
with a concrete safe‑window inequality and a statement of what is *proved* vs *assumed* (see Addendum v2, §A–B).
3. **Flat→Curved bridge** — We present a **theorem/assumption split**: what follows from QFT (first law, MI subtraction,
moment‑kill, RNC expansion) vs what is **postulated** (Clausius balance for general wedges). Any result depending
on (A2) is labeled **conditional**.
4. **Constitutive closure** — We relabel \(\delta G/G = -\beta\,\delta\sigma\) as **Hypothesis H\(_{{\sigma\to G}}\)** and provide a
**variational motivation** from maximizing wedge entropy subject to Bianchi consistency (Addendum v2, §D).

## Point‑by‑point replies

### 1) “Arbitrary cap angles” → **NFEP definition; single convention**
**Referee’s concern.** The specific angles were not first‑principles; near identity of schemes suggests tuning.

**Action & reply.** We now define **c\_geo** uniquely by NFEP:
> **NFEP.** Let \(\chi^a\) be the fixed Unruh‑normalized boost at the wedge apex and \(\hat\rho_D(u)\) the adopted null‑segment
weight (same as in the modular calculation). Define the **wedge Clausius flux** \(\Phi_\wedge \equiv \int_\wedge \delta Q/T\) and the
**FRW patch flux** \(\Phi_\text{{FRW}} \equiv \int_{\mathbb S^2}\delta Q/T\). Then \(c_\text{{geo}}\) is the unique multiplicity that makes a tiling by
**disjoint spherical caps of half‑angle \(\theta_\star\)** *flux‑equipartitioned*, i.e.
\[ c_\text{{geo}}\,\Phi_\wedge(\theta_\star;\,\hat\rho_D)\;=\;\Phi_\text{{FRW}}(\hat\rho_D) .\]
This prevents angular double counting **by construction** and leaves no free cap parameter once \(\hat\rho_D\) is fixed.
We move all numerics (including the previously quoted \(c_\text{{geo}}\approx 10.49\)) to an *illustrative* appendix and
retain **only** the NFEP definition in the main text. The former “Scheme A” is demoted to a historical convention.

**Outcome.** The near‑equality across schemes is now explained as a **redundant bookkeeping** that has been removed.
There is **no tuning** in the prediction once NFEP is adopted; all geometric freedom is fixed by \(\hat\rho_D\) and the
ball→diamond map already used in the modular calculation.

### 2) “(A2) insufficiently justified” → **Sharp theorem/assumption split + bound**
We add a **Small‑Diamond Response Lemma** (Addendum v2, §A): in Riemann normal coordinates at the apex,
for Hadamard states and \(\ell\) within the safe window \(\epsilon_\text{{UV}}\ll\ell\ll L_\text{{curv}}\),
MI subtraction plus moment‑kill cancels \(\mathcal O(\ell^0,\ell^2)\) *including curvature‑contact pieces*;
the leading surviving term is \(\mathcal O(\ell^4)\) and **state‑universal** at this order. The extension to a **Clausius balance**
for arbitrary wedges is retained as **Assumption (A2)**, now boxed and explicitly flagged wherever used.

### 3) “Flat‑space CHM → gravity” → **QFT core vs thermodynamic postulate**
We emphasize that **β is a flat‑space QFT object** (Eq. (3) in the paper) and provide a logic flow that cleanly separates
(i) QFT statements (first law, MI subtraction, moment‑kill, finite \(I_{{00}}\)), (ii) **postulated** Clausius balance (A2),
and (iii) geometric bookkeeping (NFEP‑defined \(c_\text{{geo}}\)). Any conclusion beyond (i) is marked **conditional**.
This removes the appearance of smuggling gravitational content into β.  (Paper §§III–VI updated accordingly.)

### 4) “CGM: assumed compensator; fine‑tuning” → **Hypothesis H\(_{{\sigma\to G}}\) + scaling check**
We relabel the closure \(\delta G/G = -\beta\,\delta\sigma\) as **Hypothesis H\(_{{\sigma\to G}}\)** and motivate it by an
**entropy‑maximization** argument at fixed wedge flux (Addendum v2, §D). The marginal‑log cancellation at \(\Delta=d/2\)
is shown to arise from a **slow running** of \(M^2\) (\(\delta\sigma\propto \log\ell\) within the safe window), removing the sense of
tuning. All formulas depending on this hypothesis are now marked as such.

### 5) “Pipeline helper is arithmetic only” → **Scope clarified**
Correct; `referee_pipeline_v2.py` is designed to verify **scheme invariance** and produce machine‑checkable artifacts,
not to compute β (which is a QFT calculation already documented in the manuscript and the original code). We now
state this explicitly in the README and add placeholders for **window scans** and **NFEP cap solving** hooks.

### 6) “3% β uncertainty is misleading” → **Two‑tier error budget**
We split uncertainties into:
- **Tier‑N (numerical/QFT):** integration grids, window scans, MI residual gates → \(\sigma_\beta/\beta\approx 3\%\) (unchanged).
- **Tier‑M (model):** dependence on (A2) and H\(_{{\sigma\to G}}\), NFEP assumption, FRW mapping → **quoted separately** and
**dominant** for any global claim. The abstract and conclusions are updated to reflect this hierarchy.

## Manuscript changes delivered
- **Main text**: remove dual‑scheme numerics as “predictions”; insert **NFEP definition** and theorem/assumption
annotations; adjust abstract/intro/conclusions to mark conditional scope.  
- **Appendix A (v2)**: formal moment‑kill identities; small‑diamond bound.  
- **Appendix C (v2)**: NFEP‑based c_geo definition; historical note on prior conventions; numerical values relegated to examples.  
- **Appendix H (v2)**: RNC expansion details; explicit safe‑window inequality.  
- **CGM section (Sec. V, v2)**: recast as depending on H\(_{{\sigma\to G}}\); expanded scaling at \(\Delta=d/2\).

We hope these changes address the core concerns by replacing conventions with a single physical definition (NFEP),
demarcating postulates from derivations, and clarifying the error hierarchy.
