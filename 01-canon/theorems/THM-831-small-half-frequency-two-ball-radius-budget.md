---
id: THM-831
title: Small half-frequency folded targets are the exact no-switch radius frontier
status: PROVED (all primitive ratios, sharp switch obstruction, and gcd quotient) + FINITE-EXACT (518-pair formula audit and 6,160 packet replay)
source: codex-2026-07-15-S10 continuation
depends_on: [THM-774, THM-824]
related: [THM-821, HYP-6820]
verification:
  - 04-computation/lrc13_small_half_frequency_two_ball_radius_budget_codex_S10.py
  - 05-knowledge/results/lrc13_small_half_frequency_two_ball_radius_budget_codex_S10.out
---

# THM-831 — the exact no-switch frontier of the folded two-sheet target

Write `T=R/Z`, and let `|| ||` denote circular distance to zero.  Suppose

```text
A>B>=1,       g=gcd(A,B),       alpha=A/g,       beta=B/g,
gcd(alpha,beta)=1,              alpha-beta odd.              (1)
```

Put

```text
r=alpha-beta,       s=alpha+beta,       epsilon=2/13.         (2)
```

Then `r,s` are coprime and odd.  The folded target is

```text
H_(A,B)={t:||At||+||Bt||>=11/13}.                           (3)
```

The theorem has three parts.

1. Every component of (3) has an exact Bezout-offset centre and radius.
2. For a primitive target (`g=1`), symmetric Minkowski containment has a
   no-switch radius factorization **exactly** when `4<=alpha<=9`.
3. A common gcd creates raw deck switches, but the quotient-scaled radius
   factorization remains exact in the same range.

Thus THM-824's `(A,B)=(9,4)`, equivalently odd exceptions `(13,5)`, is one
row of a sharp sixteen-row frontier rather than an isolated coincidence.

## 1. Exact component formula

Set

```text
n=floor(2alpha/13+1/2),       u=r^(-1) mod s.              (4)
```

For `j=1,...,n`, let

```text
q_j=2j-1,       p_j in {1,...,s-1},       p_j=u q_j mod s,
a_j=min(4r,4alpha-13q_j),
d_j=min(0,4beta-13q_j),
h_j=a_j/(26rs),
c_j=p_j/s+d_j/(26rs) mod 1.                                (5)
```

### Theorem 1 — the Bezout-offset ball decomposition

The primitive target is the disjoint union

```text
H_(alpha,beta)
 = disjoint_union_(j=1)^n (B_T(c_j,h_j) union B_T(-c_j,h_j)). (6)
```

In particular it is empty for `alpha<=3`.  The original target is the
`g`-fold pullback

```text
H_(A,B)
 = disjoint_union_(k=0)^(g-1) disjoint_union_(j=1)^n
     (B_T((c_j+k)/g,h_j/g) union B_T((-c_j+k)/g,h_j/g)).    (7)
```

The primitive component radii are all equal if and only if

```text
n=1       or       4beta>=13(2n-1).                        (8)
```

### Proof

THM-774's folded-diamond identity says that `t` lies in (3), in the primitive
case, precisely when there are integers `p,q` of opposite parity with

```text
|st-p|<=epsilon,       |rt-q|<=epsilon.                    (9)
```

Put `m=rp-sq`.  Since `r,s` are odd, `m` is odd exactly when `p,q` have
opposite parity.  The two teeth in (9) meet exactly when

```text
|p/s-q/r|=|m|/(rs)<=2/(13s)+2/(13r),
```

or

```text
13|m|<=4alpha.                                            (10)
```

The positive odd solutions of (10) are exactly `q_j=2j-1`, `j=1,...,n`.
For each one, `rp=m mod s` gives the unique `p_j` in (5), and the negative
offset is its reflection.

The `s`-tooth has radius `2/(13s)`, the wider `r`-tooth has radius
`2/(13r)`, and their centre offset is `q_j/(rs)`.  The first tooth is wholly
contained in the second exactly when

```text
13q_j<=4beta.                                             (11)
```

In that case the intersection has centre `p_j/s` and radius `2/(13s)`.
Otherwise its midpoint is shifted left by `(13q_j-4beta)/(26rs)` and its
radius is `(4alpha-13q_j)/(26rs)`.  These are exactly (5).  Distinct offsets
give distinct `s`-teeth, so the intervals are disjoint.  This proves (6).
Multiplication by `g` has the `g` inverse branches displayed in (7).

If `n=1`, equality of the two reflected radii is automatic.  For `n>=2`, all
radii equal exactly when the last offset is still in the containment regime
(11); otherwise the partial-intersection radii decrease strictly with `q_j`.
This proves (8). ∎

## 2. The weighted centre-switch object

The binary tournament on components is not the proof-bearing combinatorial
object.  For any disjoint circular balls

```text
K=disjoint_union_(i in I) B_T(c_i,h_i),                    (12)
```

define the ordered, `j,k`-symmetric three-hyperedge weight

```text
omega(i|j,k)
 =||2c_i-c_j-c_k||-(2h_i+h_j+h_k).                        (13)
```

### Lemma 2 — sharp no-switch criterion

The following are equivalent.

1. Whenever `t,t+z,t-z` all lie in `K`, their three component labels are the
   same.
2. `omega(i|j,k)>0` for every nonconstant label triple `(i,j,k)`.

When they hold, define the clearance of `t in K` by

```text
clear_K(t)=h_i-||t-c_i||
```

for its unique component `i`.  For every nonempty compact `E` and compact
`R=-R` containing zero,

```text
E+R subset K
iff E subset K and rho_0(R)<=min_(t in E) clear_K(t).       (14)
```

For equal radii `h_i=h`, (14) is the radius formula

```text
E+R subset K
iff rho_C(E)+rho_0(R)<=h,                                 (15)
```

and condition 2 becomes

```text
sigma(C):=min_(i,j,k not constant)||2c_i-c_j-c_k||>4h.    (16)
```

### Proof

Give `t,t+z,t-z` component labels `i,j,k`.  The circular triangle inequality
gives

```text
||2c_i-c_j-c_k||<=2h_i+h_j+h_k.                           (17)
```

Thus condition 2 forces constant labels.  Conversely, if one inequality in
(17) holds, choose the balanced representative of
`2c_i-c_j-c_k`.  The interval of values of `2e_i-e_j-e_k` with
`|e_l|<=h_l` is exactly

```text
[-(2h_i+h_j+h_k), 2h_i+h_j+h_k].
```

It therefore contains the negative centre defect.  Solving

```text
t=c_i+e_i,
t+z=c_j+e_j,
t-z=c_k+e_k
```

constructs a switched symmetric triple.  This proves the equivalence,
including equality and endpoint cases.

Under no-switch, every `z in R` keeps `t,t+z,t-z` in the same lifted interval.
The exact real identity

```text
max(|e+z|,|e-z|)=|e|+|z|
```

then gives `||z||<=clear_K(t)`.  Taking compact extrema proves necessity in
(14), and the triangle inequality proves sufficiency.  Formula (15) is the
equal-radius specialization. ∎

## 3. The exact primitive frontier

For `4<=alpha<=9`, (4) has `n=1`.  There are precisely sixteen admissible
coprime opposite-parity pairs.  Write the two components as
`B(c,h),B(-c,h)` with `0<c<=1/2`, and put

```text
Delta=min(||2c||,||4c||),       slack=Delta-4h.            (18)
```

The complete exact table is:

| `alpha` | `beta` | `c` | `h` | `Delta` | `slack` |
|---:|---:|---:|---:|---:|---:|
| 4 | 1 | 49/130 | 1/130 | 16/65 | 14/65 |
| 4 | 3 | 25/182 | 3/182 | 25/91 | 19/91 |
| 5 | 2 | 23/78 | 1/78 | 7/39 | 5/39 |
| 5 | 4 | 1/9 | 2/117 | 2/9 | 2/13 |
| 6 | 1 | 381/910 | 11/910 | 74/455 | 4/35 |
| 6 | 5 | 1/11 | 2/143 | 2/11 | 18/143 |
| 7 | 2 | 17/78 | 1/78 | 5/39 | 1/13 |
| 7 | 4 | 4/11 | 2/143 | 3/11 | 31/143 |
| 7 | 6 | 1/13 | 2/169 | 2/13 | 18/169 |
| 8 | 1 | 719/1638 | 19/1638 | 100/819 | 62/819 |
| 8 | 3 | 261/1430 | 19/1430 | 193/715 | 31/143 |
| 8 | 5 | 4/13 | 2/169 | 3/13 | 31/169 |
| 8 | 7 | 1/15 | 2/195 | 2/15 | 6/65 |
| 9 | 2 | 551/2002 | 23/2002 | 101/1001 | 5/91 |
| 9 | 4 | 5/13 | 2/169 | 3/13 | 31/169 |
| 9 | 8 | 1/17 | 2/221 | 2/17 | 18/221 |

Every slack is positive, so Lemma 2 proves the symmetric compact-set radius
factorization on all sixteen rows.

This frontier is sharp in a stronger sense than equal radii.  If
`alpha>=10`, the `q=1` and `q=3` components both occur.  Abbreviate

```text
d_q=min(0,4beta-13q),       a_q=min(4r,4alpha-13q).        (19)
```

For the label triple `(-1|-3,+1)`, the main `u q/s` terms cancel and give

```text
centre defect = ||c_3-3c_1||=|d_3-3d_1|/(26rs),
radius tax    = 3h_1+h_3=(3a_1+a_3)/(26rs).               (20)
```

Always

```text
|d_3-3d_1|<=3a_1+a_3.                                    (21)
```

Indeed:

- if `beta<=3`, the two sides reduce to `8beta` and `16alpha-78`;
- if `4<=beta<=9`, they are `39-4beta` and
  `12r+min(4r,4alpha-39)`; if `r>=2` the latter is at least 24 while the
  former is at most 23, and the sole `r=1` edge is immediate;
- if `beta>=10`, the left side is zero.

Thus Lemma 2 constructs a switched symmetric triple for every admissible
`beta` as soon as `alpha>=10`.  At the first boundary `alpha=10`, the four
components `q=+/-1,+/-3` already have negative switch slack for all
`beta=1,3,7,9`.

We have proved the exact primitive statement:

```text
H_(alpha,beta) has the symmetric no-switch radius factorization
iff 4<=alpha<=9.                                          (22)
```

Targets with `alpha<=3` are empty; (22) concerns the nonempty regime.

## 4. Common gcd: raw switches, exact quotient radii

If `g>1`, (7) contains, for any primitive centre `c`, the deck progression

```text
(c+k-1)/g,       (c+k)/g,       (c+k+1)/g.                (23)
```

Its centre defect is zero.  For `g=2` the two outer labels coincide modulo the
deck but differ from the middle label; for `g>2` all three can be distinct.
Hence the raw component no-switch property fails for every `g>1`.

The quotient statement remains exact.  For `4<=alpha<=9`, let
`C={c,-c}` and `h` be the corresponding primitive table row.  For nonempty
compact `E` and compact symmetric `R` containing zero,

```text
E+R subset H_(A,B)
iff max_(e in E) d_C(ge)+max_(z in R)||gz||<=h.            (24)
```

This follows by applying multiplication by `g` to the Minkowski sum and then
using (15) on the primitive target.  The quotient preserves folded membership
and addition while deliberately collapsing raw sheet identity.

When the row is in the full-`s`-tooth regime `beta>=4`, one may rewrite (24)
in phase form only with the target-cell guard retained:

```text
E+R subset H_(A,B)
iff gE subset H_(alpha,beta) and
    max_(e in E)||sge||+s max_(z in R)||gz||<=2/13.        (25)
```

Without the first conjunct, points near the wrong `s`-tooth give the same
false rewrite exposed by the corrected THM-824 guard.

## 5. Tournament Analysis and preservation boundary

The proof-bearing discrete object is the weighted ordered three-hypergraph
`omega(i|j,k)`, not a binary tournament.  In the full-tooth regime its centres
are the `u`-gauge image of the truncated odd comb

```text
u{+/-1,+/-3,...,+/-(2n-1)} subset Z/sZ.
```

The fatal relation `2u=3u+(-u)` is invariant under the unit gauge and is an
additive-energy witness.  A pairwise tournament cannot retain which ordered
pair supplies the right side of that relation.

For planning telemetry only, the verifier ranks the sixteen frontier ratios
by raw no-switch slack and by slack per target radius.  The two transitive
tournaments flip nine edges; each has score histogram `0,...,15`, singleton
SCCs, no directed cycle, and one Hamiltonian path.  The observable preserves a
priority order for proof search and destroys centres, endpoints, owner
transport, the three-hyperedge witness, and the LRC predicate.  Vertices are
ratio types, not runners.

The challenged assumption is that common dilation merely copies a valid raw
component quotient.  It does not: the deck progression (23) creates switches.
Only the quotient map `m_g`, which preserves the folded predicate and
Minkowski addition, licenses the scaled radii in (24).

## 6. Exact verification and scope

The Fraction-exact verifier:

1. independently reconstructs all sixteen target pairs by an affine sweep and
   by the Bezout tooth intersections;
2. checks the table, the positive second-difference gaps, and 6,160 exact
   symmetric compact packets with zero radius/direct mismatches;
3. matches formula (6) on all 518 primitive opposite-parity pairs through
   `alpha=50`, including the two empty rows;
4. verifies (8) on all 516 nonempty rows and the weighted switch frontier on
   16 positive rows through `alpha=9` versus 500 nonpositive rows from
   `alpha=10` through 50; and
5. prints the four exact negative `alpha=10` switch witnesses and Tournament
   Analysis telemetry.

The frozen digests are

```text
source  de22dfe33ce641eff186d7fb2aec5638a911ebeb1e4da4ca54f33ef2aeb3c63f,
output  8ec9a1091494c03ae1862fbbac1de2797881b9f7cbd016b50bfd2edab129fe68,
small-frontier certificate
        a39c46d37d72888b88827e3da9ef00f273a8b0a53606e1b9f916e05c3d81c7a6,
formula audit
        8cf16040f4adcf2663c867ef378bfc1d979c87ffd82b7bacea4bec7f679db05c.
```

The universal sharpness beyond 50 is the algebraic proof (19)--(21), not an
extrapolation from the finite audit.  The theorem classifies the geometry and
the legal symmetric-return quotient.  It does not prove that a two-sheet LRC
core realizes or violates the corresponding radius budget, does not factor an
individual nonsymmetric satellite, and does not close the `n=12` sporadic
branch.
