---
id: THM-3428
title: "Rough full-order half-twist rank-seven classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For every odd Q>=3 with
  spf(Q)>7, at most seven literal half-twist blocks all having quotient order
  Q cover the Q sheets iff Q is 11,13,23,or29; the exact ranks are 6,7,6,7.
  Reflection accounting, the THM-3426 odd-block clique, and an un-cancelled
  shortest-gap relation prove the Q>=512 half; an exact normalized bank closes
  all 116 smaller rough moduli, independently reproduced by a coefficient-level
  exhaustive bank.
  Mixed lower quotient orders, arbitrary common time, and LRC(14) remain open.
source: root-rough-maximal-order-rank-seven-2026-08-15
audit: independent reflection-invoice, odd-deletion, shortest-gap, threshold, and scope reconstruction; normal/optimized/stored analytic replay; independent exhaustive finite-boundary replay with coefficient multiplicities retained
depends_on:
  - THM-3426-rough-composite-odd-interval-collision-and-dyadic-clique-law
related:
  - THM-3421-prime-half-twist-rank-seven-classification
  - THM-3425-half-twist-rank-six-primitive-breaker-profile-closure
script: 04-computation/rough_maximal_order_half_twist_rank7_thm3428.py
output: 05-knowledge/results/rough_maximal_order_half_twist_rank7_thm3428.out
script_sha256: db668302a1a99eec477531b5684e1a892f15f5897fd883a2b961d034e35720cf
output_sha256: c27a498e4ee5cd41c15e10e335638e3649dfdc02af655c8098c526fa1fad6c5d
semantic_sha256: 4e81e4be4ef811b910a12537432dae02202d8e3913ea8b28633a1bc71615d3f9
independent_finite_boundary_script: 04-computation/rough_maximal_order_half_twist_rank7_finite_boundary_20260815.py
independent_finite_boundary_output: 05-knowledge/results/rough_maximal_order_half_twist_rank7_finite_boundary_20260815.out
independent_finite_boundary_script_sha256: ca3ea3d64703a88987c20d1ad148695543db6f64bb6abcd9619905ee79969fe7
independent_finite_boundary_output_sha256: d1268dc5ec58ffba6c1c32f2bd4e59790996cb07eba1593b64c66cb6b4528961
independent_finite_boundary_semantic_sha256: b61db8412c1a45419cb4be595f75c802461472f24525b017b3d50d01af9f35ef
hash_basis: LF-normalized bytes
---

# THM-3428 -- rough full-order half-twist rank-seven classification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement

For odd `Q` and a residue `r` modulo `2Q`, put

```text
B_(Q,r)={ell in Z/QZ: ||r(2ell+1)/(2Q)||<1/14}.          (1)
```

The quotient order of `(1)` is `Q/gcd(Q,r)`.  Assume

```text
Q>=3 is odd,        spf(Q)>7.                            (2)
```

Then a family `R` satisfying

```text
|R|<=7,
gcd(Q,r)=1 for every r in R,
union_(r in R) B_(Q,r)=Z/QZ.                             (3)
```

exists if and only if

```text
Q in {11,13,23,29}.                                    (3a)
```

At those four moduli the exact minimum cardinalities are respectively
`6,7,6,7`.

Every block in `(3)` is required to have the full quotient order `Q`.  No
primitivity hypothesis is needed beyond that pointwise condition.  The four
positive cases regress THM-3421.  The analytic exclusion below closes every
`Q>=512`; the finite-exact boundary closes all smaller rough moduli and finds
no composite positive.

## 2. Inheritance and connection ledger

THM-3421 closes prime rank seven by a scalar reflection invoice plus a sharp
odd-interval packing theorem.  THM-3426 proves that the same odd packing has
clique number three over every rough modulus in `(2)`, even though the odd
interval itself may contain nonunits.  The remaining residue-one even block
does not require a quotient law: its original short-gap relation already
exhibits an intersection.

| field | exact connection |
|---|---|
| source | full-quotient-order half-twist masks `(1)` |
| target | reflection-paired overlap debt plus odd/even multiplier collision graphs |
| map | parity of `r`, the fixed sheet of `ell -> -1-ell`, and multiplicative normalization |
| preserved | literal union, every overlap multiplicity, strict endpoints, and full quotient order |
| destroyed | owner identity inside a parity bank and all lower quotient-order arms |
| required sidecar | THM-3426 for odd blocks; the un-cancelled centered-interval gap for even blocks |
| cheapest tests | exact small bank through `Q<512`; rough controls `Q=517,533`; lower-order scope hostile `Q=513` |

The corrected near miss is scalar capacity: at seven blocks it leaves four
residue classes exactly at nonnegative mass.  The least-used coordinate is
the location of the reflection-fixed sheet.  The proof spends its overlap
before taking any density quotient.

## 3. Block sizes and the reflection invoice

Assume first that `Q>=512`.  It is enough to exclude a seven-block multiset:
any shorter cover can be
padded by full-order blocks without losing its union.  Write

```text
Q=14k+s,       s in {1,3,5,9,11,13},                    (4)
```

where `s=7` is absent because `spf(Q)>7`.  The reflection

```text
sigma(ell)=-1-ell                                       (5)
```

has the unique fixed sheet `ell_0=(Q-1)/2`.  Every mask `(1)` is
`sigma`-invariant.  If `r` is even, `(1)` contains `ell_0` and has size

```text
2k+1.                                                    (6)
```

If `r` is odd, it avoids `ell_0` and has size

```text
2k       for s in {1,3,5},
2k+2     for s in {9,11,13}.                             (7)
```

Let `e` be the number of even blocks.  Coverage of the fixed sheet gives
`e>=1`.  If `Omega` is total block mass minus `Q`, the fixed sheet contributes
`e-1`, and every other contribution occurs in a reflection pair.  Hence

```text
Omega=(e-1)+2B,        B>=0 an integer.                  (8)
```

Equations `(6)--(8)` give

```text
s in {1,3,5}:  Omega=e-s;
s in {9,11,13}: Omega=14-s-e,                            (9)
```

and therefore

```text
s=3,5:   impossible;
s=1:     B=0;
s=9:     1<=e<=3, B=3-e;
s=11:    1<=e<=2, B=2-e;
s=13:    e=1, B=0.                                      (10)
```

The first line is already a capacity contradiction with the fixed-sheet
invoice; no incidence estimate is used there.

## 4. Spending overlap among odd blocks

At an off-fixed reflection orbit met by `t` odd blocks, delete `t-1` of those
blocks.  Performing this for every orbit deletes at most the odd/odd portion
of `B`, hence at most `B` blocks total, and leaves a pairwise-disjoint odd
subfamily.  Thus the `7-e` odd blocks contain at least

```text
7-e-B                                                   (11)
```

pairwise-disjoint members.

For `s=9,11,13`, formula `(10)` makes `(11)` respectively `4,5,6`.
In each case the odd block is a unit dilate of the maximal odd interval from
THM-3426 with `h=7` and

```text
(s,c)=(9,2),(11,4),(13,6).                              (12)
```

Condition `(2)` implies both of that theorem's thresholds.  Its sharp clique
number is three, contradicting `(11)`.

## 5. The residue-one even collision

Suppose `s=1`.  In the odd-sheet coordinate modulo `Q`, an even block is a
unit dilate of

```text
E={-k,...,-1,0,1,...,k}.                                (13)
```

The zero is the fixed sheet.  For any unit `lambda modulo Q`, place

```text
0,lambda,2lambda,...,k lambda
```

on the circular `Q`-gon.  A shortest gap gives

```text
b lambda=+-A (mod Q),
1<=b<=k,        1<=A<=floor(Q/(k+1))<=13.               (14)
```

Since `Q>=512` gives `k>=13`, both `b` and `A` already lie in
`E minus {0}`.  Equation `(14)` is therefore an off-fixed intersection of
the two even dilates.  Crucially, no gcd is cancelled and neither endpoint
must be a unit; this lemma holds over every odd modulus in its range.

But `B=0` in `(10)`, so two even blocks are impossible.  Hence `e=1`, and
the six odd blocks are pairwise disjoint.  Here the THM-3426 parameters are
`h=7,c=8`; condition `(2)` is exactly large enough for its first threshold
and exceeds its cubic threshold.  Clique number three again gives a
contradiction.  This closes the last residue class and proves that `(3)` is
impossible for every `Q>=512` satisfying `(2)`.

## 6. Exact finite boundary and companion

For `Q<512`, the companion makes the finite universe literal:

```text
Q odd, 3<=Q<512, spf(Q)>7, gcd(Q,r)=1.                 (15)
```

There are exactly `116` such moduli.  The reflection-fixed sheet belongs only
to blocks with even `r`, so every cover contains an even block `B_(Q,2a)`.
Because `a` is a unit modulo `Q`, it has an odd unit lift modulo `2Q`;
multiplying the odd sheet coordinates by that lift permutes the full block
family and sends `B_(Q,2a)` to `B_(Q,2)`.  The search therefore fixes
`B_(Q,2)` without loss.

For each modulus the companion collapses only literally equal masks, branches
on an uncovered sheet with the fewest available masks, uses the sum of the
largest remaining gains as a rigorous pruning bound, and memoizes the exact
pair `(covered mask, remaining depth)`.  Iterative depths through seven prove
the minimum, not merely existence.  Across `4,735` visited states it finds

```text
112 negative moduli;
(Q,minimum rank)=(11,6),(13,7),(23,6),(29,7).           (16)
```

It directly replays witnesses

```text
Q=11: (1,2,3,5,7,9);
Q=13: (1,2,3,5,7,9,11);
Q=23: (1,4,5,7,9,11);
Q=29: (1,5,7,8,12,13,22).                              (17)
```

Mask collapse is harmless here: two coefficients defining the same literal
set are interchangeable in a union cover, and repeating either cannot reduce
the minimum rank.

The independent coefficient-level boundary bank retains distinct coefficient
choices even when their masks coincide.  It checks the same `116` moduli and
all feasible owner counts through seven via `249` scalar profiles, `26,333`
states, and `646,035` branches.  Its support and minimum ranks agree exactly
with `(16)`, and no one of the `23` composite moduli is positive.

The standard-library companion directly reconstructs the literal masks on
the first two rough composite controls in every possible residue class:

```text
s=1:  533,589;       s=3: 689,703;       s=5: 551,649;
s=9:  527,583;       s=11:529,697;       s=13:517,559.  (18)
```

For all twelve moduli it checks the sizes `(6)--(7)`, fixed-sheet ownership,
reflection, and the exact normalized odd packing number three.  In the two
residue-one controls it also checks every normalized even block meets the
base block away from the fixed sheet.  It verifies `(14)` for `1,988` values
`13<=k<=2000`.

Two scope controls are frozen separately:

- the order-29 family in `(17)` is the largest positive finite boundary;
- at `Q=513`, the five-block family `(1,57,285,342,399)` covers with quotient
  orders `(513,9,9,3,9)`, showing why a theorem that simply drops the
  maximal-order and roughness gates would be false.

Normal and optimized transcripts are byte-identical for both companions.  The
bounded banks prove exactly the finite side `(15)`; the arbitrary large range
is supplied by Sections 3--5, not inferred from either bank.

## 7. Boundaries and non-consequences

- `512` is only the seam between the analytic proof and the exact finite bank;
  it is not asserted to be a natural sharp constant.
- The even collision `(14)` is stronger than a field ratio-set statement and
  deliberately retains nonunit endpoints.
- THM-3426's least-prime-factor gate is used only for the odd block.  Small
  prime factors create the breaker arms isolated by THM-3425.
- The pointwise full-order hypothesis in `(3)` is stronger than joint period
  `Q`.  Families mixing proper quotient orders are not classified here.
- The result concerns the half-twist Boolean layer of zero complete mode
  cochains.  It gives no arbitrary common-time, physical-runner, or LRC(14)
  conclusion.
