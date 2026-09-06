---
id: THM-3010
title: "Ballot-column Newton ratios, the bronze log-concavity threshold, and metallic maximal alternation"
status: CORRECTED / PROVED identities and recurrence; historical finite controls retained
source: klein-S428
depends_on:
  - THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer
  - THM-3004-circuit-sign-change-cluster-law-and-classifier-refutation
related:
  - THM-438-paley-cluster-integrals-are-catalan
  - 04-computation/three_strand_sequence.py (kind-pasteur-2026-03-14-S87)
  - THM-224-golden-exceptional-points
  - THM-3003-antipodal-circuit-rigidity-and-the-multipole-spread-criterion
script: 04-computation/gmc_ballot_column_newton_ratios_and_metallic_alternation_thm3010.py
output: 05-knowledge/results/gmc_ballot_column_newton_ratios_and_metallic_alternation_thm3010.out
script_sha256: 11b115d861eae1f31416e5fc8fe7089f023fc58c93412cdc4d79511a218ff703
output_sha256: 00cee818d95c7dd288e1d16bc9e501db0d26a7adf4aacb90e158e637163245b3
hash_basis: LF-normalized bytes
---

# THM-3010 -- ballot-column Newton ratios and metallic alternation

**CORRECTED, 2026-09-06 continuing8.** The displayed identities and metallic
recurrence survive. The former nonrational-Fuss assertion and the general
norm/alternation dichotomy are **REFUTED**. The ballot quadratic and bronze
classification now have all-integer proofs in the
[correction and proof](../../05-knowledge/results/continuing8_20260906_ballot_repair.md).
The source/output pins above retain the original historical computation; they
do not verify the corrected universal statements.

Notation is THM-3000 section 1, but the sequence `h_0=1, h_1, h_2, ...` is taken
directly from a classical integer family rather than from a polynomial's
coefficients; `R_k=h_k^2/(h_(k-1)h_(k+1))` and `c_k=log(R_k/R_(k-1))` as usual.

## 1. Ballot columns have exactly rational Newton ratios

**Theorem.**  For `h_k=binom(2k+a,k+b)` the Newton ratio is an exact rational
function of `k`,

    R_k = 1 - Q_(a,b)(k)/D_(a,b)(k),                       (1)

with a canonical monic quadratic before cancellation and a positive denominator.
Precisely, put `delta=a-2b` and assume `k>=max(1,1-b,1-a+b)`. Then

    Q = k^2+(a+1-delta^2)k+[a(a+2)-(2a+1)delta^2]/4,
    D = (k+b)(k+a-b)(2k+a+2)(2k+a+1)/2.

Literal factorial cancellation proves (1) for every integer `(a,b)`. A reduced
numerator can have smaller degree. The column is log-concave (`R_k>1`) exactly
where `Q(k)<0`. Computed rows:

| `h_k` | `1-R_k` |
|---|---|
| `binom(2k,k)` | `1/(k(2k+1))` |
| Catalan `binom(2k,k)/(k+1)` | `3/((k+1)(2k+1))` |
| `binom(2k,k-1)` `= 1,4,15,56,210` | `(k^2-3k-1)/((k-1)(k+1)^2(2k+1))` |
| `binom(2k+1,k-1)` `= 1,5,21,84,330` | `(k^2-7k-6)/((k-1)(k+1)(k+2)(2k+3))` |

The five-letter circuit words at the originally sampled indices `k=3..7` are
Catalan `+++++`, `binom(2k,k)` `+++++`, `binom(2k+1,k-1)` `-----`, and
`binom(2k,k-1)` `---++`. These are finite windows. The third column also changes
concavity later: `k^2-7k-6` is negative at 7 and positive at 8.

**Correction: every integer order has rational adjacent ratios.** For `p>=2`,

    binom(pk,k)/binom(p(k-1),k-1)
      = product_(i=1)^p[p(k-1)+i]
        / {k product_(i=1)^(p-1)[(p-1)(k-1)+i]}.

Dividing the column by `(p-1)k+1` retains rationality and gives its Fuss--Catalan
analogue. For `p=3` that adjacent ratio is `3(3k-1)(3k-2)/(2k(2k+1))`.
The quadratic numerator phenomenon does not extend automatically: the reduced
numerator of `1-R` for the `p=4` Fuss--Catalan column has degree four.

## 2. The unique metallic discriminant in the family is BRONZE

Call `Q` **metallic** when it is monic of the shape `k^2-nk-1`, the defining
equation of the `n`-th metallic ratio (`n=1` golden, `n=2` silver, `n=3` bronze),
whose two roots have product `-1`.

**Theorem, now global.** Over all integer `a,b`, the only columns with metallic
canonical monic quadratic are `(a,b)=(0,-1)` and its mirror `(0,+1)` -- exactly

    h_k = binom(2k,k-1) = 1, 4, 15, 56, 210, 792, ...

with `Q=k^2-3k-1`, the **bronze** equation.  Consequently that sequence is
log-**concave** for `k` below the bronze ratio and log-**convex** above it, with
the crossing at

    (3+sqrt13)/2 = 3.302775638...,

bracketed exactly between `k=3` and `k=4` in the computed pattern `>><<<<<`.
No golden discriminant occurs anywhere in the family. For the proof, the
constant coefficient condition gives `a(a+2)+4=(2a+1)delta^2`. With `D=2a+1`,

    4D delta^2 = D^2+2D+13, so D divides 13.

The negative divisors give `delta^2=-3`. The positive divisors `1,13` give
`delta^2=4`: `a=0` yields the bronze quadratic, whereas `a=6,b=2,4` yields
`k^2+3k-1`, outside `n>=1`. Discriminant 13 prevents cancellation against
the rational linear denominator factors, so the reduced degree-two version
has the same classification.


## 2a. Connection to the repo's three-strand sequence (kind-pasteur-2026-03-14-S87)

`04-computation/three_strand_sequence.py` already identified the sequence

    1, 1, 2, 3, 4, 6, 10, 15, 20, 35, 56, 70, 126, 210, 252, 462, 792, 924, ...

as the interleaving of three Pascal strands
`A: binom(2n+1,n)`, `B: binom(2n+2,n)`, `C: binom(2n,n)`.  Equivalently -- and
this is the compact description -- it is the **central band of Pascal's triangle
of width two**, i.e. the row-distinct values of `binom(n,k)` with `|n-2k|<=2`.
Checked exactly: rows `n=1..9` give precisely its first `13` terms.

Section 1 now gives every strand a Newton-circuit meaning, and the answer is a
clean trichotomy:

| strand | `h_k` | `1-R_k` | behaviour |
|---|---|---|---|
| A | `binom(2k+1,k)` | `1/D` | log-convex throughout, no change of sense |
| C | `binom(2k,k)` | `1/(k(2k+1))` | log-convex throughout, no change of sense |
| B | `binom(2k,k-1) = 1,4,15,56,210` | `(k^2-3k-1)/D` | **flips at the bronze ratio** |

So of the three strands, two have constant discriminant `Q=1` and never change
sense, while the third -- and only the third -- carries a metallic discriminant
and changes from log-concave to log-convex at `(3+sqrt13)/2`.  That is the
Newton-circuit content of kind-pasteur's decomposition, and it identifies strand
`B` as the distinguished one.

## 3. Metallic recurrences attain the MAXIMAL circuit alternation

Let `x^2-nx-1=0` with `n>=1`, and `a_k=n a_(k-1)+a_(k-2)`, `a_0=0`, `a_1=1`
(`n=1` is Fibonacci).  The two roots have product `-1`, and the Simson/Catalan
identity

    a_(k-1)a_(k+1)-a_k^2 = (-1)^k                          (2)

is precisely the **norm form** of that quadratic order.  Therefore

    R_k = a_k^2/(a_k^2+(-1)^k)                             (3)

alternates strictly about `1` for `k>=2`; the circuit alternates for `k>=3`, and the
sign-change count is the **maximum** a circuit of that length admits.  Verified
exactly for `n=1,2,3,4,5`: `9` of `9` changes in each case, with (2) returning
`+1,-1,+1,-1,...` throughout.

**Relation to THM-3004.** The separated double-cluster polynomial family and
this recurrence both attain maximal circuit alternation in their respective
types. This is a shared conclusion, not an identification of recurrence
characteristic roots with coefficient-polynomial root parameters. The golden
ratio is the `n=1` member of the recurrence family.


## 3a. Repaired two-element distinction; the general dichotomy is false

Metallic quadratics `x^2-nx-1` have root **product `-1`**, so the pair is
`{lam, -1/lam}`: closed under `r -> -1/r`, an *anti*-reciprocal map, and
containing a negative element.  THM-3003 section 1's rigidity requires
`{r} = {mu/r}` with `mu = e_d^(2/d) > 0` -- norm `+1` and positive roots.

The particular positive two-element reciprocal pair cannot equal the metallic
characteristic pair: their products have opposite signs. This repairs the
containment claim attributed to HYP-9070. It does **not** separate the possible
circuit words of the two different constructions.

An exact hostile to the former general claim is

    N(n)=(n+1)^2(n+3)^2(n+9)^2.

Its root parameters are fixed by `r -> 9/r`. With `h_k=e_k/binom(6,k)`, its
five Newton ratios are

    (65/57, 4693/4005, 71289/61009, 4693/4005, 65/57).

Its four circuit signs are `+,-,+,-`: antipalindromic **and** maximally
alternating. Thus the general implication from positive reciprocal symmetry
to nonalternation, and the claimed converse from maximal alternation to a
metallic norm, are refuted. Other JC or Euclidean-algorithm assertions are not
dependencies of this repair.

## 4. Boundaries

- Sections 1 and 2 now have all-integer proofs; the original scanned rectangle
  remains historical finite evidence, not the source of their quantifiers.
- Section 3 is proved from (2), which is classical; the maximality statement is
  about the `h`-sequence, and these `h` need not come from a real-rooted
  polynomial.  Nothing here claims a real-rooted realization.
- Section 3a preserves only the stated two-element product distinction. Its
  former general norm/circuit dichotomy is refuted by the exact six-root hostile.
- Nothing here bears on no-return for the first-gap family, GMC(2), or ULC, and
  nothing here is a bridge to JC(2) -- MISTAKE-237 retracted one such bridge and
  HYP-9070 is explicitly a stratification, not a reduction.

## 5. Reproduction

    python3 04-computation/gmc_ballot_column_newton_ratios_and_metallic_alternation_thm3010.py

The historical producer reported four parts `True`, but did not test its false
prose implications. Use the linked continuing8 correction's standalone source,
raw outputs and exact certificate for the repaired statements.
