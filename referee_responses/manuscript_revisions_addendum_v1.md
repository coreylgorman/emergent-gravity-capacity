
# Manuscript Revisions — Addendum v1

## A. Lemma A (Moment‑kill & ℓ⁴ isolation)

**Statement.** Let \(K_{B_\ell}\) be the CHM modular Hamiltonian with weight \(W_\ell(r)=(\ell^2-r^2)/(2\ell)\).
Define \(K_{\text{sub}}(\ell)=K(\ell)-aK(\sigma_1\ell)-bK(\sigma_2\ell)\), with \(a,b\) chosen so that the weighted
integrals of any smooth radial \(F(r)=F_0+F_2 r^2+O(r^4)\) satisfy
\[\int W_\ell F - a\int W_{\sigma_1\ell}F - b\int W_{\sigma_2\ell}F = O(\ell^6).\]
Then, to leading nontrivial order,
\[\delta\!\langle K_{\text{sub}}(\ell)\rangle = (2\pi\,C_T\,I_{00})\,\ell^4\,\delta\sigma(x)+O(\ell^6),\]
with the \(O(\ell^0,\ell^2)\) pieces cancelled and the \(\ell^4\) coefficient finite and regulator independent.

**Proof (sketch).** Expand correlators and curvature contributions in Riemann normal coordinates around the apex.
Area/contact terms scale as \(\ell^0\) and \(\ell^2\); curvature–contact terms share the same moment structure and are
canceled by the chosen \(a,b\). The surviving isotropic piece scales as \(\ell^4\) and defines \(I_{00}\).

## B. Normalization Theorem (Scheme invariance)

**Statement.** Fix: (i) ball→diamond map and Unruh normalization, (ii) generator density \(\hat\rho_D(u)\), and
(iii) a no‑double‑counting cap tiling of the FRW 2‑sphere. Then any bookkeeping split into \(f\) and \(c_{\rm geo}\)
preserves the product \(x=\beta f c_{\rm geo}\).

**Proof (sketch).** The Clausius/Noether charge for one diamond is linear in the bulk modular response and in the
boundary/horizon conversion. Changing the split multiplies the diamond contribution by a constant and divides the
cap multiplicity by the same constant, leaving the product invariant.

## C. Derivation of \(c_{\rm geo}\) and cap angles

For a spherical cap of half‑angle \(\theta_\star\), \(\Delta\Omega=2\pi(1-\cos\theta_\star)\). If the FRW 2‑sphere is tiled by
\(c_{\rm geo}\) disjoint caps of equal flux (no double counting), then
\[\boxed{~c_{\rm geo}=\frac{4\pi}{\Delta\Omega}=\frac{2}{1-\cos\theta_\star}~}.\]
Hence:
- Scheme A (minimal cap): \(c_{\rm geo}=40 \Rightarrow \cos\theta_\star=1-2/40=19/20\).
- Scheme B (equal‑flux cap under the adopted \(\hat\rho_D\)): numerical evaluation yields \(c_{\rm geo}\simeq10.49\)
  \((\cos\theta_\star\simeq0.80934)\).

## D. Appendix H (Flat→Curved local map; proof sketch)

Work in Riemann normal coordinates centered at the wedge apex. For \(\ell\) within the safe window,
the state is Hadamard and the local vacuum two‑point function admits the Minkowski leading term plus curvature
corrections. After MI subtraction and moment‑kill, all local terms up to \(O(\ell^2)\) cancel; the universal \(O(\ell^4)\)
coefficient multiplies \(\delta\sigma(x)\). Curvature contamination appears at \(O(\ell^6)\) and is subleading for the
windows we use.

## E. CGM resolution (expanded derivation)

Allow \(G\) to run via \(M^2(x)=(8\pi G)^{-1}\) and \(\delta G/G=-\beta\delta\sigma\). Varying
\(S_{\rm grav}=A/[4G(x)]\) gives
\[\delta S_{\rm grav}=\frac{1}{4G}\delta A+\frac{A}{4G}\beta\,\delta\sigma.\]
With Unruh \(T=\kappa/2\pi\) and the linearized Raychaudhuri relation this restores the Clausius balance at the same
order as the CHM response, including the marginal log at \(\Delta=d/2\).

## F. Phenomenology note

We keep \(\alpha_T=\alpha_B=0\) and only \(\alpha_M\) active. Observable, frame‑invariant signatures include
\(\,d_L^{\rm GW}/d_L^{\rm EM}\) and a direct mapping \(\alpha_M(0)=-\dot G/(H_0 G)\). We demote
\(a_0=(\Omega_\Lambda^2/2)cH_0\) to an auxiliary scale.
