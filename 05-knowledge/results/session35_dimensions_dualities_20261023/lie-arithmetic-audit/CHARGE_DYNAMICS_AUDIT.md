# Independent audit: elliptic rational charge dynamics

Verdict: **PASS**, no material errors. Audited the full `session35-charge-dynamics/ELLIPTIC_CHARGE_DYNAMICS.md` and producer `checks.py` / `checks.json`. No Git operations.

## Independent derivation and scope

The root t=(9+i√319)/20 has norm one, trace 9/10, and irreducible polynomial 10t²−9t+10. A root of unity would have integral rational trace, contradicting 9/10. Thus θ/π is irrational, all positive powers avoid ±1, and all positive-index charges are strictly between 1 and 25. The conjugate branch reverses θ and preserves every charge.

If θ/π were algebraic irrational, Gelfond–Schneider applied to base −1 and the logarithm iπ would make exp(iθ) transcendental. This contradicts the displayed quadratic equation. Consequently θ/π is transcendental; this does not use or assert that a product of transcendental numbers must be transcendental.

The integer recurrence A_(m+1)=9A_m−100A_(m−1), A_0=2, A_1=9 gives A_m odd and prime to 5 for every m≥1. For m≥2, the unreduced charge numerator 13·10^m−6A_m has 2-adic valuation exactly one and no factor 5. The separate m=1 value 38/5 has precisely the claimed denominator. Hence reduced denominators are 10^m/2, are strictly increasing, and certify distinctness. Negative powers give the same charges, so “all integer powers” must mean consecutive positive exponents for distinctness, as explicitly indexed in the note.

An independent useful version of the discriminant identity is: if t^m=a_m+b_m√(−319), then Δ(c_m)=−319(12b_m)². Here b_m is rational and nonzero by infinite order. This directly certifies the exact square class and quadratic field without relying on the Chebyshev recurrence.

The geometric-sum proof gives vanishing nonconstant Fourier averages along all consecutive positive integers. Uniform trigonometric approximation yields weak convergence on the circle. The cosine map c=13−12cosφ pushes that measure to density 1/[π√((c−1)(25−c))]. Endpoint singularities are integrable and nonatomic, so the interval-frequency and density conclusions are valid. No uniform Fourier rate is assumed.

Constant Laurent coefficients give mean 13, variance 72, centered even moments 36^k binomial(2k,k), odd moments zero, and fixed-lag covariance 36(t^r+t^−r)=36S_r. Lag-one correlation is 9/20. These are deterministic empirical limits, not independence claims.

The dyadic subsequence m=2^k is expressly excluded from the equidistribution theorem. Invariance of the arcsine measure under P_2 does not establish genericity of this seed. The unit-circle branch supplies neither positive-real mirror bases nor an intermediate-E7.5 VOA construction. These exclusions are essential and correctly stated.

## Exact independent checks

`charge_independent_checks.py` uses binomial expansions of (9+√(−319))^m, rather than the producer recurrence, through m=80. It checks norms, denominators, open interval, discriminant square class, distinctness, affine recurrence and doubling identity. Laurent-polynomial constant coefficients check centered moments through degree 20 and the symbolic fixed-lag coefficient 36(v+v^−1). Producer checks through m=150 also reran successfully. Finite tests do not certify the analytic limiting or transcendence theorems; those were proof-reviewed separately.

## Source spot-check

Freshly opened the author-hosted Richard Elman, *Lectures on Abstract Algebra*, [§74, Theorem 74.29 and Corollary 74.30](https://www.math.ucla.edu/~rse/algebra_book.pdf), printed p.442, PDF page 450 (zero-based index 449). The theorem requires an algebraic base different from 0 and 1 and algebraic irrational exponent. Its proof writes the base as exp(γ) for an arbitrary applicable logarithm, validating the value with γ=iπ. The source also explicitly treats base −1 in the nearby examples. No fresh verification of Gelfond's 1934 original paper is claimed or needed for this audit. The equidistribution argument is proved directly in the note rather than imported with unchecked hypotheses.

## Final synthesis scope pins

The separate ungraded completion result is over C: finite-dimensional quadratic Lie algebras containing the specified e7⋉H57 have minimum dimension 248, and equality forces simple E8 with an emergent compatible contact grading. The separately audited real statement concerns the split e7(7) and its absolutely irreducible real 56-dimensional module; equality then forces split e8(8). It is not a classification for arbitrary real forms. Neither statement produces a VOA.
