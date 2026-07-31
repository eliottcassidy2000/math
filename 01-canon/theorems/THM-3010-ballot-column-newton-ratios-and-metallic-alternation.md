---
id: THM-3010
title: "Ballot-column Newton ratios, the bronze log-concavity threshold, and metallic maximal alternation"
status: PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT
source: klein-S428
depends_on:
  - THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer
  - THM-3004-circuit-sign-change-cluster-law-and-classifier-refutation
related:
  - THM-438-paley-cluster-integrals-are-catalan
  - THM-224-golden-exceptional-points
  - THM-3003-antipodal-circuit-rigidity-and-the-multipole-spread-criterion
script: 04-computation/gmc_ballot_column_newton_ratios_and_metallic_alternation_thm3010.py
output: 05-knowledge/results/gmc_ballot_column_newton_ratios_and_metallic_alternation_thm3010.out
script_sha256: c49c17605939b804333b3bb16881892b97d6d8149ddcd4a6d730ff7252d3a079
output_sha256: 2818ed7d24d1a4c1c1b6d65282cf51b972eaaa77d23f15736216b4a07120dcc3
hash_basis: LF-normalized bytes
---

# THM-3010 -- ballot-column Newton ratios and metallic alternation

**PROVED + VERIFIED-EXACT.**

Notation is THM-3000 section 1, but the sequence `h_0=1, h_1, h_2, ...` is taken
directly from a classical integer family rather than from a polynomial's
coefficients; `R_k=h_k^2/(h_(k-1)h_(k+1))` and `c_k=log(R_k/R_(k-1))` as usual.

## 1. Ballot columns have exactly rational Newton ratios

**Theorem.**  For `h_k=binom(2k+a,k+b)` the Newton ratio is an exact rational
function of `k`,

    R_k = 1 - Q_(a,b)(k)/D_(a,b)(k),                       (1)

with `Q` monic in `k` of degree at most `2` and `D>0` on the range.  Hence the
column is log-concave (`R_k>1`) exactly where `Q(k)<0`.  Computed rows:

| `h_k` | `1-R_k` |
|---|---|
| `binom(2k,k)` | `1/(k(2k+1))` |
| Catalan `binom(2k,k)/(k+1)` | `3/((k+1)(2k+1))` |
| `binom(2k,k-1)` `= 1,4,15,56,210` | `(k^2-3k-1)/((k-1)(k+1)^2(2k+1))` |
| `binom(2k+1,k-1)` `= 1,5,21,84,330` | `(k^2-7k-6)/((k-1)(k+1)(k+2)(2k+3))` |

Circuits, computed exactly: Catalan `+++++`, `binom(2k,k)` `+++++`,
`binom(2k+1,k-1)` `-----`, and `binom(2k,k-1)` `---++` -- the only one of the
four with a sign change.

**`p=2` is forced.**  `binom(2k,k)/binom(2k-2,k-1)=2(2k-1)/k` is rational in `k`,
but `binom(pk,k)/binom(p(k-1),k-1)` is not a rational function of `k` for `p>=3`.
So there is no Fuss--Catalan analogue of this table; the hypergeometric-term
condition singles out the `p=2` ballot triangle.

## 2. The unique metallic discriminant in the family is BRONZE

Call `Q` **metallic** when it is monic of the shape `k^2-nk-1`, the defining
equation of the `n`-th metallic ratio (`n=1` golden, `n=2` silver, `n=3` bronze),
whose two roots have product `-1`.

**Theorem.**  Over `a in [-2,4]`, `b in [-3,2]`, the only columns with metallic
discriminant are `(a,b)=(0,-1)` and its mirror `(0,+1)` -- that is, exactly

    h_k = binom(2k,k-1) = 1, 4, 15, 56, 210, 792, ...

with `Q=k^2-3k-1`, the **bronze** equation.  Consequently that sequence is
log-**concave** for `k` below the bronze ratio and log-**convex** above it, with
the crossing at

    (3+sqrt13)/2 = 3.302775638...,

bracketed exactly between `k=3` and `k=4` in the computed pattern `>><<<<<`.
No golden discriminant occurs anywhere in the family.

## 3. Metallic recurrences attain the MAXIMAL circuit alternation

Let `x^2-nx-1=0` with `n>=1`, and `a_k=n a_(k-1)+a_(k-2)`, `a_0=0`, `a_1=1`
(`n=1` is Fibonacci).  The two roots have product `-1`, and the Simson/Catalan
identity

    a_(k-1)a_(k+1)-a_k^2 = (-1)^k                          (2)

is precisely the **norm form** of that quadratic order.  Therefore

    R_k = a_k^2/(a_k^2+(-1)^k)                             (3)

alternates strictly about `1`, the circuit alternates at **every** index, and the
sign-change count is the **maximum** a circuit of that length admits.  Verified
exactly for `n=1,2,3,4,5`: `9` of `9` changes in each case, with (2) returning
`+1,-1,+1,-1,...` throughout.

**Relation to THM-3004.**  THM-3004 proved the maximal count `d-3` is attained by
the real-rooted family `prod_(j)(n+T^j)^2`.  Section 3 is the **h-sequence face of
the same phenomenon**: both extremes are "two geometric terms whose ratio is
negative" -- geometrically as interleaved bands whose dominance alternates, and
algebraically as a quadratic order of norm `-1`.  The golden ratio is the `n=1`
member, so `phi` is not incidental here: it is the smallest metallic parameter,
hence the extreme case of maximal Newton-circuit alternation.

## 4. Boundaries

- Section 1 is an identity for each displayed row; the "degree `<=2`" statement is
  verified over the scanned rectangle, not proved for all `(a,b)`.
- Section 2's uniqueness is over `a in [-2,4]`, `b in [-3,2]` only.  It is a
  FINITE-EXACT statement about that window, not a theorem for all `(a,b)`.
- Section 3 is proved from (2), which is classical; the maximality statement is
  about the `h`-sequence, and these `h` need not come from a real-rooted
  polynomial.  Nothing here claims a real-rooted realization.
- Nothing here bears on no-return for the first-gap family, GMC(2), or ULC.

## 5. Reproduction

    python3 04-computation/gmc_ballot_column_newton_ratios_and_metallic_alternation_thm3010.py

Three parts, all reporting `True`.
