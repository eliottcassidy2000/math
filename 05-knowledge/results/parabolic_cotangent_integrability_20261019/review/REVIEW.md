# Independent publication review: parabolic cotangent integrability

19 October 2026. **PASS, within the stated scope.** No substantive mathematical defect was found in the frozen main proof or its formal-unfolding extension. This is an independent proof audit with finite arithmetic controls, not a historical-priority assessment.

## Main results

1. **Formal invariant ring.** The cotangent lift is tangent to identity and its logarithm raises total-degree filtration. Thus its fixed ring equals the kernel of the Hamiltonian derivation. Expanding in momentum reduces the kernel equation to `(A_m/v^m)'=0` over the Laurent field. Characteristic zero gives precisely the stated ring of all formal compositions in `pv`; no polynomial-only restriction is hidden in the argument. Hamiltonian uniqueness is understood with zero additive constant.

2. **Meromorphic integral iff convergent iterative logarithm.** A convergent two-variable quotient embeds into Laurent series in `p` over the field of meromorphic `q`-germs. Its leading coefficient satisfies the required Julia equation. If that coefficient has exponent zero, it is constant and can be subtracted; nonconstancy then yields a nonzero exponent. A nonzero power of `v` is consequently meromorphic, and its positive formal order plus an analytic unit root proves convergence of `v`. This handles poles and nonlinear dependence on momentum, rather than assuming an integral linear in momentum. The converse and time-one embedding equivalences follow from the holomorphic flow.

3. **All-rational obstruction.** The leading Laurent exponent in momentum can be negative. The valuation at infinity then gives `s=m`, and the leading coefficient gives `d^m=1`. This is impossible for nonzero integer `m` in characteristic zero. The proof covers every nonlinear polynomial, without requiring a parabolic fixed point or a birational lift. The earlier repository result [quadratic_clock_return_20261012a.md](../../quadratic_clock_return_20261012a.md), Section 5, item 3, proves the polynomial-invariant version; the manuscript correctly credits it and identifies the extension to the full rational field.

4. **Finite-jet family and period-three residue.** The total-degree jet calculation includes the extra momentum factor correctly. Integrability differs despite arbitrarily long jet agreement, while the globally ramified lifts remain distinct from polynomial automorphisms. The period-three quotient calculation gives `b/a^2=-2/49` and residue `51/49` at all three orbit points. Time is explicitly measured per third return, so no factor of three is lost. The polynomial-completion lower-bound argument is also sound.

5. **Sectorial gluing.** The final text specifies canonical univalent attracting/repelling Fatou coordinates on sectors of opening greater than pi with their common formal asymptotic. This is important: arbitrary narrow-sector Abel solutions need not have the same uniqueness. The cotangent overlap formula and the condition `h'=1` are correct. Uniform sector asymptotics give the removable singularity after gluing. Constant logarithmic monodromy is compatible with gluing; it is not itself the obstruction.

## Analytic inputs checked against primary sources

[Aschenbrenner–Bergweiler, *Julia's equation and differential transcendence*](https://www.mat.univie.ac.at/~maschenbrenner/pdf/itlog-final.pdf), printed page 3, explicitly records the classical convergence classification for globally meromorphic parabolic maps, with only parabolic Mobius exceptions. The application to nonlinear polynomials is valid. Differential transcendence over entire coefficients is not being substituted for this divergence theorem. Page 4 separately treats the stronger convergent-coefficient result, including the quadratic cusp.

[Dudko–Sauzin, *The resurgent character of the Fatou coordinates of a simple parabolic germ*](https://arxiv.org/pdf/1307.8093), Definition 1, Theorem 2 and Section 2, supplies the sector domains, uniqueness up to constants, differentiated asymptotics and lifted horn maps used here. The manuscript's gluing conclusion is its own elementary consequence of this input and the displayed overlap map.

## Formal unfolding extension

The joint `(q,delta)`-adic construction is correctly distinguished from a point-adic logarithm at nonzero numerical `delta`. Divisibility by `q^2+delta^2`, formal Weierstrass division, evaluation on the two fixed sections, and formal linearization give the stated unique `A`, `B`, and `U`. The apparent singularities in `A,B` are removable and their formulas converge. No convergence of the whole time form or real transit-error estimate is inferred from finite or formal algebra.

The final extension also includes a sound further consequence: `U(q,0)` must diverge. Otherwise `(1+q)/q^2+U(q,0)` would be a meromorphic reciprocal of the cusp iterative logarithm, forcing that logarithm to converge. The numerical-delta primitive is correctly described as multivalued analytic, not an element of the joint formal power-series ring. Neither observation is needed for the main invariant theorems.

## Reproduction and finite scope

- The frozen upstream main checker passes **186** always-active exact checks in normal and optimized Python modes.
- The frozen upstream unfolding checker passes **15** always-active exact checks in both modes.
- `independent_controls.py` uses only Python rational arithmetic, imports no author code and no SymPy, and passes **129** checks in both modes. It constructs the cusp logarithm via an independent Julia-coefficient recurrence (through degree 20, comparing the published degree 14), reconstructs return coefficients from the square factorization, checks the finite-jet family for `N=2,...,40`, and inverts the scalar logarithm to check the discrete period and unfolding coefficients.
- Finite coefficients do not prove divergence, universal rational-invariant exclusion, or meromorphic nonintegrability. Those conclusions rely on the proofs and identified analytic input.

The reviewed source hashes are recorded in [reviewed-source-sha256.json](reviewed-source-sha256.json). [Publication overlays](../publication-overlays.json) document relative-link, reproduction-path and presentation-only updates; theorem statements and mathematical computations are unchanged.
