# Session 28: three narrow historical corrections

This package proposes changes to exactly three existing targets. No repository edits have been made by this lane. `repair-proposals/manifest.json` records exact preimage and proposed SHA256 hashes; `.preimage` and `.proposed` files carry complete bytes, and `corrections.patch` is a unified diff. A publisher must compare the live source with the recorded preimage before applying; stop on mismatch rather than overwrite concurrent changes. Existing source, including leading shebangs or maintained banners, is preserved. The available originals of these three targets did not contain a maintained warning at byte zero.

## 1. THM-134: downgrade only the unsupported continuous claim

Target: `01-canon/theorems/THM-134-paley-local-max-parseval.md`.

Machine-readable status becomes:

    OPEN (continuous-extension claim; historical proof incomplete)

The warning quotes the exact historical status, `PROVED (p=7,11 — algebraic + computational)`, rather than silently erasing it. Discrete interpolation and permutation symmetry do not specify an intrinsic continuous Hamiltonian-count function on the whole Parseval simplex, and symmetry does not force every stationary point to be uniform. A finite-sample-vanishing symmetric polynomial can alter the uniform Hessian arbitrarily. The continuous-extension inference is therefore unproved. This is a proof gap, NOT a counterexample to the local maximum of an explicitly specified polynomial. Direct Hessian calculations for such a polynomial, and separately established discrete counts, survive within their own scopes. This warning was independently reviewed by the spectral audit lane.

## 2. Historical n−2 reflection: point to the existing canonical refutation

Target: `07-reflections/n-2-negative-eigenvalues.md`.

Add a prominent historical warning, preserving the old discussion. THM-300 already records the exact n=9 inertia (20 positive,0 zero,8 negative), refuting the claimed7 negative. The session independently reproduces it but does not claim a new refutation. Finite observations through n=8 remain valid; an all-order replacement is not certified.

## 3. Tournament L-function script: comments only, identical Python AST

Target: `04-computation/tournament_L_function.py`.

Add warning comments after the shebang, retaining the original docstring and computation. The finite reciprocal product is1/D(s), with D entire, so it has no finite zeros. For positive Q its factor equations describe poles with real part log(Q)/log(p), not nontrivial zeros. Paley Q=(p+1)/4, not sqrt(p); quadratic-residue orientations give tournaments for p=3 mod4, not p=1 mod4. Fourier-mode labels are not rational primes; no functional equation follows from the definition. The warning withdraws the unsupported RH/automorphy interpretations, not the explicitly defined numerical function.

`package_corrections.py` verifies equality of parsed Python ASTs, ignoring source locations. Thus the only script change is comments; computational behavior and docstring are unchanged.

## Suggested MISTAKES ledger entry (publisher assigns unique ID or uses dated bullet)

### 2026-10-22: Conflating sampled spectra, signed coordinates, and spectral interpretations

Three historical scope corrections were applied with preimage hashes:

- THM-134's prior continuous Parseval-simplex PROVED status was unsupported by discrete polynomial fitting and permutation symmetry. The status is now OPEN for the continuous-extension claim. Specified-polynomial Hessians and verified finite counts are not retracted. Symmetry makes the uniform point stationary under appropriate differentiability, but does not make it the only stationary point, nor determine Hessian sign.
- The n−2 negative-eigenvalue reflection omitted the existing THM-300 refutation at n=9. Its historical finite observations are retained with a warning. A fixed-path Hessian at a transitive corner must not be confused with the full-arc unbiased Hessian or a spectral interpolation Hessian.
- The exploratory tournament L-function script called poles of a finite reciprocal product zeros and inferred RH-style functional-equation/automorphy structure without a derivation. A comment-only, AST-preserving warning now states the correct pole geometry and Paley normalization.

Prevention: specify the coordinate space, extension, group action, positivity class, and exact analytic object before transferring spectral vocabulary. For continuous conclusions, either prove an identity on an open set or explicitly declare the chosen extension. For moment claims, test exact polynomial-square positivity. A familiar singularity or numerical resemblance does not supply a self-adjoint positive realization or a modular transformation law.
