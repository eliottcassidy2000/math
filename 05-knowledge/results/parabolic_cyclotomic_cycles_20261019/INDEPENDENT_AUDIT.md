# Independent audit of the parabolic–cyclotomic supplement

Date: 2026-10-19. Verdict: **ACCEPT the corrected supplement.**

The rational neutral-cycle exclusion and the three-point cyclotomic affine bijection are valid. One misleading preservation statement was identified and corrected: the affine set map changes the cycle multiplier from 8 to 1. The final Section 6 now explicitly distinguishes permutation intertwining from derivative or multiplier preservation.

## Main claims

### Rational neutral cycles

For any exact rational cycle of length `n≥3`, the cited denominator theorem gives the common exact point denominator `b` with `4|b`. Since each numerator is coprime to `b`, its multiplier has

`v_2(λ)=n(1−v_2(b))<0`.

This immediately excludes both rational real-neutral possibilities `λ=1` and `λ=−1`, independently of any conjectural upper bound on the period. At an odd denominator prime, the multiplier valuation is likewise exactly `−n v_p(b)`, proving local repulsion there. The hypotheses are precisely those of Eliahou–Fares Proposition 2.1 and Theorem 2.12, verified directly in the preceding independent denominator audit; switching the sign convention for the parameter does not change denominator valuations.

For exact period two, subtraction forces `u+v=−1`. Multiplier one forces the double root `u=v=−1/2`, so this is not an exact two-cycle. Multiplier minus one forces a monic quadratic of discriminant two, precluding rational points. At period one, `λ=2x` gives exactly the two stated fixed points. Thus the claimed classification covers all rational periods, without assuming the Flynn–Poonen–Schaefer conjecture.

### Cyclotomic affine bijection

The cubic `g(t)=t³+t²−2t−1` is irreducible by the rational-root test and has discriminant 49. An independent resultant computation gives

`Res_z(1+z+⋯+z⁶, tz−z²−1)=g(t)²`.

This verifies its identification with the three real seventh-cyclotomic conjugates. The map `σ(t)=t²−2` has exact order three on those roots, and reduction modulo `g` gives `σ²(t)=−t²−t+1`.

For `L(t)=−t−1/2`, direct substitution proves

`f_(-7/4)∘L = L∘σ²`

on the specified roots. Since `L` is invertible, this is an actual bijection of the two three-point sets with the asserted inverse orientation. Its image polynomial is `8X³+4X²−18X−1`; it is irreducible by the affine substitution, and its roots form an exact three-cycle of multiplier one. The original cyclotomic cycle under `f_(-2)` has multiplier eight. Thus the corrected text properly avoids claiming a global dynamical conjugacy or multiplier preservation.

## Other symbolic and analytic statements

The independent verifier reconstructs the generic third dynatomic polynomial by dividing `f_c³(X)−X` by `f_c(X)−X`. It confirms its factorization into the two trace cubics after `c=−T²−T−2`, their square discriminant, and the multiplier formula by polynomial reduction. Eliminating the trace gives exactly the displayed quadratic in the multiplier and both stated square-root branches.

It also reconstructs the marked parameterization from the inverse parameter `t=x_0²+x_0+c`, rather than importing the primary verifier's expression for `c`. This independently checks both unfolding identities and the attracting rational control at `t=5/4`.

The cubic roots are simple and avoid `t=0,−1`; the coefficient of `g(t)` in the multiplier identity is nonzero at each root. Therefore the parameter map has quadratic ramification, while the multiplier crosses one to first order. Rational parameters on the attracting side approach the boundary, giving infinitely many distinct rational parameters locally and disproving any uniform positive real multiplier gap. Exact neutrality is nevertheless excluded by the denominator argument above.

The `c=−2` hostile is sound: the two trace cubics are irreducible by the rational-root test and have rational multipliers `8` and `−8`. Rational trace, multiplier, and square discriminant do not imply pointwise rational roots.

## Sources and independent computation

The current primary [Morton–Raianu paper](https://arxiv.org/pdf/2105.07435), arXiv v5 dated January 2022, was read directly. Its equations (3)–(4) give the marked parameterization and its inverse; the text also explicitly describes the seventh-cyclotomic norm for the parameter numerator. The supplement does not import its graph conjectures or claim priority for the cyclotomic arithmetic.

The independent script is `verify_bridge_independent.py`. It imports no primary script and uses generic dynatomic division, resultants, modular polynomial iteration, and reconstruction from the parameter inverse. Normal and optimized Python runs pass and produce byte-identical outputs:

- `independent_verification.json`
- `independent_verification_optimized.json`

The arbitrary-period neutral exclusion is a mathematical consequence of the cited denominator theorem; it is not presented as a finite computational check.

Reproduce with:

```sh
python verify_bridge_independent.py
python -O verify_bridge_independent.py
```

No checkout or Git writes were performed.
