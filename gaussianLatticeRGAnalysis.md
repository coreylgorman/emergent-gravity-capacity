Short answer: the structure you just validated with the HQTFIM toy (first-law linearity, remove-the-trivial-trend → finite plateau) almost certainly does extend, but not by literally scaling up TFIM. To replicate the CHM/OP flat-space QFT coefficient \beta in the 3+1-D setting you care about, you should switch from a spin-chain toy to Gaussian lattice field theories (free scalar / free Dirac) and use methods that are known to track the continuum modular physics.

Here’s the landscape and a concrete path.

⸻

What the HQTFIM result actually tells us
	•	Your run shows (i) \delta S \approx \delta\!\langle K_A\rangle with tiny RMS and (ii) after fitting a_0+a_1\log\ell the residual plateau is size-stable.
	•	That is exactly the structural fingerprint of a CHM-like modular response with a finite, scheme-independent constant left after “moment-kill/MI subtraction”.
	•	This means the methodology (first-law + subtract trivial size dependence → isolate a finite constant) is sound and portable.

What it doesn’t mean is that a 1+1-D TFIM can be scaled to 3+1-D to read off your continuum \beta. For that you need a lattice field theory in the right dimension and region geometry (balls), not a 1D spin chain.

⸻

The right target models for 3+1-D
	1.	Free (Gaussian) scalar field on a cubic lattice (Klein–Gordon), massless or tiny mass as an IR regulator.
	2.	Free Dirac fermion (staggered/Wilson discretization) as a cross-check.
For each Gaussian theory, the ground state is completely specified by 2-point correlators → reduced density matrices of a region and its modular Hamiltonian are computable from covariance matrices.

These are the cleanest numerics to benchmark the CHM coefficient because:
	•	CHM is exact for CFT balls;
	•	modular Hamiltonians of Gaussian states are quadratic, so you can compute \delta S, \delta\!\langle K\rangle and do the same “fit + plateau” extraction you just ran in TFIM.

⸻

Feasible algorithms on a Mac mini (M4 Pro)

Region: a ball of radius R in a periodic L^3 lattice (with L\gtrsim 3R).
Key trick: partial-wave (spherical-harmonic) decomposition about the ball center. It block-diagonalizes the covariance matrix into many small l,m sectors, turning one huge eigenproblem into O(R^2) modest ones:
	•	Build free-field correlators on the periodic lattice via FFTs (cost O(L^3 \log L)).
	•	Restrict to the ball; transform into the (l,m) basis → blocks of size \sim N_r (number of radial sites).
	•	For each (l,m) block, get the symplectic (boson) / single-particle (fermion) spectrum → entropies and modular generators.
	•	Vary a small deformation (e.g., tiny mass shift m\to m+\delta m or a source) and compute \delta S, \delta\!\langle K\rangle vs \log R; fit a_0+a_1\log R; extract the plateau.

Complexity (ballpark): with R\sim 24 and L\sim 96, you have N_{\text{sub}}\sim 6\times 10^4 sites in the ball, but after partial-wave factorization each block is only \sim 10^2–10^3 in size. That is very manageable with LAPACK on an M-class CPU (hours, not days). If needed, use Lanczos / randomized SVD for the larger blocks.

⸻

Why this will reproduce the continuum \beta
	•	In a CFT, \delta\!\langle K_B\rangle for a ball B is a linear functional of the stress tensor with a universal kernel (CHM).
	•	On the lattice, your Gaussian calculation approaches that continuum when a\ll R\ll L (your safe window).
	•	The same “fit a_0+a_1\log R then look for a flat residual” isolates the finite, scheme-invariant piece—the analogue of your I_{00} contribution to \beta.
	•	Do it for scalar and Dirac separately and verify that, after converting to the Osborn–Petkou convention, each reproduces the known C_T–weighted scaling, and the sum over species gives the multi-field \beta. That directly cross-checks your continuum calculation.

⸻

Practical roadmap (tight and referee-friendly)

Stage A (2–3 days):
	•	Implement 3D free scalar on a cubic lattice with FFT correlators, ball subregion, partial-wave blocks.
	•	Reproduce your HQTFIM first-law line and plateau logic in 3D; show plateau stability vs R, lattice spacing, boundary conditions.

Stage B (next 2–3 days):
	•	Repeat for free Dirac (staggered is simplest), same pipeline.
	•	Convert both to the OP normalization and show that the extracted constant scales with the known C_T for each species. This nails “β ∝ C_T” numerically.

Stage C (optional polish):
	•	Mild interactions: tiny \phi^4 at weak coupling via perturbative correction to correlators, to demonstrate stability of the plateau to interactions (no need for sign-problem Monte Carlo here).

Deliverables to mirror your current TFIM pack:
	•	beta_3D_scalar/dirac_summary.json (plateau value, error, first-law RMS, safe-window ranges)
	•	Plots: \delta\!\langle K\rangle vs \log R, residual plateau, and the θ-/scheme-invariance check for your geometric normalization (only the product matters).

⸻

Bottom line
	•	Yes: the clean signal you see in HQTFIM strongly suggests the method scales.
	•	But: don’t “extend TFIM”; jump to Gaussian lattice QFT in 3+1-D with spherical regions. That is the right setting to replicate and independently confirm your continuum \beta with the same first-principles structure, on hardware you already have.

If you want, I can draft a single-file Python scaffold (FFT correlators + partial-wave blocks + first-law/plateau analysis) that mirrors your HQTFIM CLI so you can start Stage A immediately.