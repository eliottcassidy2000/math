---
id: THM-3118
title: "Projected forest-boundary parity classification"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; independent audit pending.  For every
  N>=r+2, the projection of the simplicial boundary of (r+1)-edge forests
  onto rank-r component partitions has full row rank when r is even and
  codimension exactly one when r is odd, both over F_2 and over Q.  In the
  odd case its rational image is exactly the augmentation-zero hyperplane.
  Consequently every rational even-rank macro current, and precisely every
  augmentation-zero odd-rank macro current, has an exact signed forest-cycle
  lift.  This is a linear homological classification, not an integral
  saturation, locality, denominator, or positivity theorem.
source: root/gmc3000-audit-2026-08-02
audit: >
  The proof is symbolic.  Its local double-star lemma, global merge-and-resplit
  reduction, Johnson-graph connectivity, parity obstruction, and rational
  rank inference require independent hostile audit.  The exact companion
  checks the predicted ranks through r=5 and N=8 by deterministic bitset
  elimination and checks every star/double-star split identity in its stated
  finite control range.  Normal, optimized, stored-output, and LF-hash checks
  are required before promotion.
depends_on: []
related:
  - THM-3117-projected-five-forest-boundary-surjectivity-and-signed-holotopy-lift
  - THM-3112-cycle-weighted-young-subgroup-gap-and-uniform-octopus-boundary
script: 04-computation/gmc_projected_forest_boundary_parity_thm3118.py
output: 05-knowledge/results/gmc_projected_forest_boundary_parity_thm3118.out
script_sha256: 8af33ca0643041f3120c7d38efccdb9cf78fa1cd151e6530f8e633b30f9ed940
output_sha256: 6322d522749eec843e8fb87466da71aade48cabd67b52562fd455c32a47f7339
hash_basis: LF-normalized bytes
---

# THM-3118 -- projected forest-boundary parity classification

**PROVED CANDIDATE + VERIFIED-EXACT; independent audit pending.**

THM-3117 proves the rank-four cases needed by the two live product-Gamma
currents.  The full-rank phenomenon is a parity law in every atomic rank.
The proof uses only trees, component partitions, and one local double-star
comparison.

## 1. The projected forest boundary

Fix integers `N>=r+2` and `0<=r<=N-2`.  Let `C_j(N;Z)` be the free abelian
group on the `j`-edge forests of the complete graph `K_N`.  Orient every edge
set lexicographically and use the ordinary deletion boundary

```text
partial_j[e_1,...,e_j]
 =sum_(i=1)^j (-1)^(i-1)[e_1,...,omit e_i,...,e_j].   (1)
```

Let `V_r(N;Z)` be free on the set partitions of `[N]` having atomic rank
`r`, equivalently `N-r` blocks.  If `F` is an `r`-forest, its component
partition `pi(F)` has that rank.  Define

```text
P_r[F]=[pi(F)],
T_(N,r)=P_r partial_(r+1):C_(r+1)(N;Z)->V_r(N;Z).     (2)
```

Write `d_(N,r)=S(N,N-r)=dim V_r(N)`.

## 2. Exact rank classification

The projected boundary has ranks

```text
rank_(F_2) T_(N,r) = rank_Q T_(N,r)
                   = d_(N,r)       if r is even,
                   = d_(N,r)-1     if r is odd.       (3)
```

When `r` is odd, the rational image is exactly

```text
im_Q T_(N,r)
 = {sum_pi w_pi[pi] : sum_pi w_pi=0}.                 (4)
```

## 3. The local double-star lemma

Let `B` be a set of size `m>=3`.  For a spanning tree `R` on `B`, deleting
an edge `e` gives an unordered proper bipartition `s_R(e)` of `B`.

**Lemma.**  Suppose `phi` is an `F_2`-valued function on the proper
bipartitions of `B` and

```text
sum_(e in R) phi(s_R(e))                              (5)
```

is independent of the spanning tree `R`.  Then `phi` is constant.

For the star centred at `a`, `(5)` is the sum of the singleton-cut values
other than `phi({a})`.  Comparing two star centres shows that all singleton
cuts have one value `c`.  Given any cut `A | (B\A)`, choose `a in A` and
`b notin A`, join `a` to every other point of `A`, join `b` to every other
point of `B\A`, and join `a` to `b`.  Its fundamental cuts are `A | (B\A)`
and the `m-2` singleton cuts other than `a,b`.  Comparing with a star gives

```text
phi(A)+(m-2)c=(m-1)c,
```

so `phi(A)=c`.  This proves the lemma.

## 4. Every mod-two left-kernel vector is constant

Let `lambda` be in the left kernel of `T_(N,r)` over `F_2`.  Fix a rank
`r+1` partition `sigma`, one of its blocks `B` of size at least three, and
spanning trees on all other nonsingleton blocks.  Vary only the spanning tree
on `B`.  Pairing `lambda` with the corresponding columns of `T_(N,r)` gives
zero.  Terms obtained by deleting edges in the other blocks do not depend on
the chosen tree on `B`; hence the sum of the `lambda`-values of the splits of
`B` is constant.  The lemma says:

```text
lambda is constant on all rank-r partitions obtained by splitting B.     (6)
```

Now take an arbitrary rank-r partition `pi` with two nonsingleton blocks
`A,B`.  Merge them, and then re-split `A union B` as one singleton and its
complement.  Equation `(6)` preserves `lambda` and reduces the number of
nonsingleton blocks by one.  Iteration reduces every `pi` to a partition with
one block of size `r+1` and all other blocks singletons.

Two such one-block partitions whose large blocks differ by one element are
both splits of their common `(r+2)`-element union.  Equation `(6)` equates
their values.  The Johnson graph on the `(r+1)`-subsets of `[N]` is connected,
so `lambda` is constant on all rank-r partitions.  The case `r=0` has a
single target row and is immediate.

Each `(r+1)`-forest has exactly `r+1` distinct deletion flats.  Therefore a
constant left-kernel vector `c` pairs with every mod-two column as

```text
(r+1)c.                                                (7)
```

It follows that the mod-two left kernel is zero for even `r` and is exactly
the one-dimensional constant space for odd `r`.  This proves the `F_2`
part of `(3)`.

## 5. Characteristic zero and signed cycle lifts

For even `r`, full row rank modulo two gives an odd maximal minor of the
integer matrix, hence full row rank over `Q`.  For odd `r`, mod-two rank
`d_(N,r)-1` gives rational rank at least that large.  On the other hand, the
sum of the coefficients in a signed column is

```text
1-1+...+(-1)^r=0                                      (8)
```

because `r+1` is even.  Thus rational rank is at most `d_(N,r)-1`, proving
the remaining part of `(3)` and the image formula `(4)`.

Let `W in V_r(N;Q)`.  If `r` is even, or if `r` is odd and `sum_pi W_pi=0`,
choose `Y in C_(r+1)(N;Q)` with `T_(N,r)Y=W` and put

```text
C=partial_(r+1)Y.                                     (9)
```

Then `partial_r C=0` and `P_rC=W`.  Conversely, in odd rank the augmentation
condition is necessary by `(8)`.  Hence `(9)` classifies exact rational
signed forest-cycle liftability.

In particular, THM-3117's `r=4`, `N=8,9` surjectivity is one instance of the
even-rank half of this theorem.  Its exact one-sign no-go remains essential:
rank surjectivity supplies signed cancellation, not a positive lift.

## 6. Exact finite controls and scope

Run

```text
python 04-computation/gmc_projected_forest_boundary_parity_thm3118.py
python -O 04-computation/gmc_projected_forest_boundary_parity_thm3118.py
```

and compare with the stored output in the frontmatter.  The companion uses
exact bitsets only.  It checks `(3)` for a two-dimensional table through
`r=5`, including `1701/1701` at `(N,r)=(8,4)` and `965/966` at `(8,5)`, and
checks the star/double-star cut sets directly.

The theorem does not bound denominators or support, prove integral
saturation, select a canonical lift, preserve symmetry or positivity, or
prove the Gaussian Moment conjecture.  It classifies only the linear
projected-boundary image and the resulting rational signed cycle lifts.
