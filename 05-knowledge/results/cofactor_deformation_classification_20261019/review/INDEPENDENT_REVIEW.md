# Independent review of the cofactor-deformation results

19 October 2026. Review of the three supplied notes and their verifiers. Reviewed input hashes are preserved in `../SOURCE_PROVENANCE.json`; the published receipts record their corresponding published hashes.

## Verdict

**The three mathematical results pass this review.** I found no gap in the finite odd-part reduction, the complete eventual-periodicity proof for the fixed deformation parameter `t=2`, or the classification and sharp bounds for all two-cycles at arbitrary positive integer `t`.

In particular, the `t=2` argument is an all-exponent proof. Its conclusion does not depend on treating the finite exponent grid or a finite orbit census as exhaustive. The stated large-parameter regions, exact dyadic guards, and finitely many leftover exponent pairs cover every input. Well-founded descent then proves entry into the listed fixed points or six nontrivial cycles.

There are no required mathematical corrections. A checker-mode limitation and minor optional wording improvements are recorded below.

## 1. Domain and notation

Throughout, `t` is a fixed positive integer deformation parameter, and

    F_t(n) = n - n/c + 2tc - 1,
    c = the exact power of two dividing the positive even integer n.

The superscript in `F_t²` means composition of this same map twice. It does not change the parameter. The rule maps positive even integers to positive even integers.

The distinction from the aliquot map is essential. When `n=cp` and `p` is an odd prime, multiplicativity gives

    s(cp) = (2c-1)(p+1)-cp = (c-1)p+2c-1 = F_1(cp).

This identity is restricted to that input domain. It is not an equality between the full extended maps, and the prime-odd-part domain need not be preserved by iteration. For example, `F_1(18)=12`, whereas `s(18)=21`. Even on prime-odd-part inputs the `t=2` rule is different: `F_2(14)=14`, whereas `s(14)=10`.

The reviewed results therefore establish no general aliquot convergence theorem. They also establish no eventual-periodicity theorem for every `t>=3`. No literature-priority claim was checked or inferred.

## 2. Finite odd-part reduction

Let `p+1=dk` with `d` a power of two and `k` odd. Growth is precisely `p+1<2tc`. If `k>=2t+1`, this forces `d<c`. Since both are powers of two, `c/d` is even. Consequently

    F_t(n)/d = (c/d)(dk+2t-1)-k

is odd, proving the first image's exact core is `d`; this guard is not merely a divisibility assertion.

The displayed identity

    d(F_t²(n)-n)
      = (c+d)(2t(d-1)-p) + 2(t-1)d + p+1

is correct. From `p>= (2t+1)d-1` and the strict growth inequality, its right-hand side is strictly smaller than

    c(1-d)-d²-d < 0.

All signs and strictness are preserved. There is no hidden exponent bound. The residual possibilities are exactly the `t` odd coefficients `1,3,...,2t-1`, with two still-unbounded dyadic exponents. Fixed points are exactly `c(2tc-1)`, since the bracket is positive and odd.

## 3. Exhaustion of the t=2 proof

First separate immediate descent (`p+1>4c`) and fixed points (`p+1=4c`). Under growth, the previous theorem handles every odd coefficient `k>=5`. Thus only `k=1,3` remain.

For `d<c`, the common first and second iterate formulas are correct, and the exact first-image core is `d`.

### The k=3 lane

Growth forces `d<=c`.

- If `d<c`, write `e=c/d>=2`. The cases `d>=4,e>=4`; `e=2,d>=8`; and `d=2,e>=8` cover all but `(d,e)=(2,2),(2,4),(4,2)`. These give starts `20,40,88`, respectively, exactly as treated in the note.
- If `d=c`, the uniform six-step route applies for `c>=32`. The only smaller dyadic cores are `2,4,8,16`, giving starts `10,44,184,752`. The first two are in the exceptional set; the last two have the listed six-step descents.

### The k=1 lane

Growth allows exactly `d<c`, `d=c`, or `d=2c`.

- For `d=2c`, the large branch begins at `c=8`; the smaller starts are `6,28`.
- For `d=c`, it begins at `c=16`; the smaller starts are `2,12,56`.
- For `d<c`, put `e=c/d`. If `e<d`, the dyadic alternatives are `e<=d/4` or `e=d/2`. The first implies `d>=8`. The second is uniform for `d>=16`; its smaller starts are `24,224`.
- If `e=d`, the uniform descent begins at `d=16`; the smaller starts are `4,48,448`.
- If `e>d`, put `h=e/d>=2`. For `d>=16`, the three classes `h=2`, `h=4`, and `h>=8` cover all possibilities. For `d=8`, only `h=4` remains outside those uniform routes, giving `1792`. For `d=4`, `h=2` is already covered; `h=4,8,16,32` give `192,384,768,1536`, and `h>=64` has the dedicated five-step route. For `d=2`, `h=2,4,8` give `8,16,32`, and `h>=16` has the dedicated four-step route.

Thus the leftover finite starts are precisely those supplied in the exceptional set and twelve descent chains. No extra exponent region is being dismissed on numerical evidence.

## 4. Exact cores and strict inequalities

Every subsequent asserted core was checked as an exact valuation. The independent script gives **23 universal odd-quotient certificates**. It substitutes the stated minimum divisibility into each claimed quotient and verifies that the resulting polynomial has an odd constant term and even nonconstant coefficients. These certificates apply to arbitrary integer multipliers, and hence to all allowed dyadic exponents.

Representative quotient certificates include:

- `k=3,d=c`: at successive threshold substitutions, the quotients after division by the asserted cores are `3`, `28u²-1`, `84u²+3`, `126u²+3`, and `882u²+13`.
- `k=1,e=d/2,d=16u`: the third iterate divided by four is `512u³+6u-1`.
- `k=1,h=4,d=16u`: the third iterate divided by eight is `2048u³+128u²-26u+1`.
- `k=1,d=4,h=64u`: the relevant quotients are `504u+3` and `882u+13`.

The threshold restrictions matter. For example, at start `224` the third iterate has core `8`, not `4`; at `1792` it has core `32`, not `8`; and at `752` the fifth iterate has core `2`, not `4`. The proof correctly routes these smaller values to explicit chains rather than applying a failed guard.

The independent script also supplies **15 universal strict-negativity certificates**. After shifting each parameter by its stated lower threshold, the relevant descent polynomial has a negative constant term and no positive coefficient. This proves negativity on the whole nonnegative real orthant and, in particular, on every dyadic exponent value. For the `e<=d/4` branch it certifies the note's upper bound `-d²/4+5d/4+2` for `d>=8`; the reduction to that upper bound is valid because `3d+1>0` and `-3d/e<0`.

The certificates are printed in full in the receipts. They supplement the written inequalities; finite testing is not being substituted for the universal argument.

## 5. Finite exceptions, cycles, and the global conclusion

All twelve displayed descent chains are consecutive literal iterates, end below their starts, and have at most six steps. All six listed cycles close exactly, have distinct members within each cycle, and are pairwise disjoint. All thirteen exceptional starts enter one of these cycles.

For every other nonfixed input, the proof supplies a smaller positive even iterate within six steps. Strong induction on the starting integer is an equivalent formalization of the final descent argument: the smaller iterate already eventually reaches a listed terminal object, so the initial input does as well. This excludes divergence and any additional cycle. The infinitely many fixed points cause no problem for this induction.

An additional consequence is that the universal six-step bound cannot be replaced by five while keeping the same exceptional set. The displayed chain

    56 -> 80 -> 138 -> 76 -> 72 -> 94 -> 54

has its first strict descent at step six.

## 6. All-t two-cycle classification and sharpness

Solving the two actual affine cycle equations gives exactly the odd parts stated in the note. Their fractional terms sum to `4t-2`. The common denominator and both numerators are odd, so integrality forces the fractional terms to be odd integers. Positivity and exact dyadic cores then follow; there is no unproved formal-core realization assumption.

For `c=2^a<d=2^b`,

    gcd(c+d-1,2d-1) = gcd(2c-1,2d-1)
                   = 2^gcd(a+1,b+1)-1.

For the first equality, use `2(c+d-1)-(2d-1)=2c-1` and the oddness of `2d-1`. The second is the ordinary Euclidean algorithm on exponents for Mersenne numbers. Therefore the stated divisibility criterion is both necessary and sufficient, and the solution is unique for its two distinct cores. A same-core two-cycle is impossible because the second iterate of the common affine branch has the same unique fixed point.

The finite bounds are also correct. For an existing cycle, the integral odd quantity `j` obeys

    2t-1 < j <= 4t-3,
    j ≡ 2t-1 (mod c).

Hence `c<=2t-2`. Rearranging the equation for `j` gives a positive integer denominator in

    d = [j(c-1)+(2t-1)] / [(4t-2)-j],

and bounding that denominator below by one yields `d<=8(t-1)²`. This proves finiteness for each fixed parameter without an orbit census.

The family `t=c/2+1`, `d=2c²` has the stated positive odd parts and genuine two-cycle. It attains **both** bounds simultaneously, for infinitely many parameters. It does not produce infinitely many two-cycles for a single parameter. The specialization `t=2` gives exactly the two-cycle `40 <-> 66`, consistent with the global six-cycle list.

## 7. Reproduced and independent controls

The supplied `t2-verification.json` and `verify_two_cycles.out` were reproduced byte-for-byte. The supplied two-cycle verifier also reproduced its receipt under `python -O`.

The separate `independent_controls.py` imports neither author verifier and uses explicit `require` checks rather than Python assertions. Both ordinary and optimized full runs passed, and their JSON results agree. The runs include:

- 36 exact symbolic identities, 23 universal exact-core certificates, and 15 universal strict-negativity certificates.
- 524,288 `t=2` exponent-grid inputs, with both exponents through 256 and odd coefficients `1,3,5,7,9,11,17,257`.
- Literal testing of every positive even input through 2,000,000. The only nonfixed inputs with no descent within six steps are exactly the thirteen stated exceptions. There are nine fixed inputs in this interval.
- 193,929 qualifying growth inputs for the all-t reduction, across `t=1,...,64` and dyadic exponents through 40.
- Exact two-cycle classification for `t=1,...,1024`, producing 1,699 cycles in total; 479,232 independent gcd/divisibility/bound checks, including core ranges beyond the proved bounds.
- A direct two-cycle search from 100,000 even starts per parameter for `t=1,...,32`, agreeing exactly with the algebraic classification within the searched starting-value range.
- The sharp two-cycle family through exponent 256.

These finite counts are controls only. The universal identities, guards, inequalities, exponent partition, and induction are what establish the unbounded results.

### Checker-mode limitation

The author's `verify_t2_global.py` uses `assert`, so it must be run with assertions enabled. An optimized run must not be counted as verification. A deliberately corrupted cycle in an in-memory copy fails in ordinary mode but still produces `PASS` under `python -O`; this diagnostic is recorded in `original-assert-mutation-control.json`. This does not invalidate the authentic ordinary-mode receipt, and the note does not promise optimized-mode safety.

The independent checker's intentional negative control fails in both modes, as recorded in `independent-negative-controls.json`. A corrupted cycle is also rejected in both modes; see `independent-cycle-mutation-controls.json`. An AST check confirms that the independent checker contains no `assert` statements.

## 8. Optional wording improvements

These are editorial and do not affect the verdict:

1. The finite-reduction note still calls a bounded-block proof for `t=2` a next target. Once the three notes are published together, a cross-reference to the completed `t=2` note would remove that chronological ambiguity.
2. “The restriction of `F_1` to inputs with odd prime odd part agrees with the aliquot map” is a precise formulation of the source-domain statement. It should not be read as a claim that accidental pointwise equalities with the aliquot sum can never occur elsewhere. For example, `F_10(54)=s(54)=66`.
3. For the cycle equations in their literal displayed ordering, the usual coefficient matrix has determinant `-(c+d-1)`; reversing one row gives the positive determinant used in the text. Saying “nonzero determinant of magnitude `c+d-1`” would remove this harmless sign convention.

## Reproduction

From this review directory:

    python independent_controls.py --output independent-controls.normal.json
    python -O independent_controls.py --output independent-controls.optimized.json

The two receipt files should agree. Each of these must fail with a nonzero exit status:

    python independent_controls.py --negative-control
    python -O independent_controls.py --negative-control

The source author checks can be rerun with ordinary `python`; only the two-cycle author checker is safe to count under `python -O`.
