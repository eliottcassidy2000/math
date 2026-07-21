# Hostile audit of the Frobenius lowest-face proof of NC2

**codex-2026-07-21. Verdict: PASS.** I found no load-bearing gap in THM-2022.
This note records the checks independently, especially where earlier NC2
attempts failed.

## The proof invariant

For a fixed exact support write

```text
M_m(c)=sum_(|r|=m, q.r=0) Mult(m;r) A(r)! c^r.
```

Minimizing `a.x` over the balanced probability simplex exposes a face `F` with
`a_i-lambda*q_i=delta`. DvdK gives a nonzero complete face sum
`Q=CT(f_F^m0)`. At a good prime `p`, the desired invariant is not a unique
channel but the complete residue after common-factorial normalization:

```text
M_(p*m0)/(p*A0)! = Q^p mod pfrak.                         (A)
```

Identity (A) is exactly what survives arbitrary cross-atom cancellation.

## Audit checklist

1. **Algebraic descent is legitimate.** A complex torus zero makes the
   localization of the rational moment quotient nonzero. A maximal residue
   field is a finite extension of `Q` by Zariski's lemma, and localization keeps
   every exact-support coefficient nonzero.
2. **The LP face is rational and balanced.** Rational primal data give a
   rational dual pair `(lambda,delta)`. Complementary slackness puts an optimal
   balanced point on the equality face, so zero lies in its charge hull.
3. **No angular projection collision is hidden.** Equal charges on the equality
   face force equal `a`, hence equal `b`; exact monomial support then forces the
   indices equal. DvdK therefore applies to the actual face coefficients.
4. **The face height is integral where used.** Every balanced length-`m0` face
   allocation has `A=delta*m0=A0`; existence of a nonzero term in `Q` makes
   `A0` a nonnegative integer.
5. **A good prime exists.** Only finitely many rational primes ramify or support
   a nonzero valuation of one of the finitely many `c_i` or `Q`. Choose outside
   them and above `m0`.
6. **Kummer's layer is exact.** Since the base-`p` digits of `p*m0` are
   `(0,m0)` with `m0<p`, no units-column carry occurs iff every `r_i` has units
   digit zero, i.e. `r=p*s`. Every other allocation pays at least one valuation.
7. **The common normalization is integral.** The face inequality gives
   `A(r)>=p*A0`, so every channel factorial is divisible by `(p*A0)!` before
   any reduction. No choice `p>A0` and no local division by a nonunit is needed.
8. **Off-face dilation also vanishes.** If `r=p*s` is balanced but uses an
   off-face index, strict dual inequality gives the integer gap
   `A(s)>=A0+1`. The quotient `(p*A(s))!/(p*A0)!` contains the factor
   `p*(A0+1)`.
9. **The tied face cannot cancel.** On the face the factorial quotient is one.
   Multinomial Lucas and residue-field Frobenius turn the complete tied sum
   into `Q^p`, nonzero because `Q` was chosen a unit.
10. **Endpoint cases survive.** A neutral-only lowest face has `m0=1` and
   `Q=c_i`; `A0=0` is allowed. A many-circuit face is also allowed because (A)
   preserves the sum rather than separating its monomials.
11. **The final quantifiers are correct.** The argument excludes algebraic
    torus null points on each finite exact support; algebraic descent excludes
    complex ones. A nonzero support not strictly positive or strictly negative
    has zero in its one-dimensional charge hull. Strict one-sidedness gives
    moment nullity by charge, and THM-1540 then gives GMC(2).

## Assumptions challenged

- The minimum face need not have one or two vertices.
- The minimum valuation need not select one channel.
- Archimedean magnitude dominance is unnecessary.
- Primitive atoms need not vanish separately.

The faithful object is the entire face initial form. Kummer suppresses every
other layer; Frobenius preserves the face sum. This is precisely the
phase-sensitive certificate that the earlier tournament and radial-channel
programs were seeking.

Cross-links: THM-2022, THM-1630, THM-2019, THM-2020, MISTAKE-211,
MISTAKE-213, HYP-8765, HYP-8766.
