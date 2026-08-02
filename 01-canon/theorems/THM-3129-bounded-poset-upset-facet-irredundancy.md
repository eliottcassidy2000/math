---
id: THM-3129
title: "Bounded-poset upset-facet irredundancy"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.  In a finite
  connected poset, the extreme rays of the isotone cone modulo constants are
  exactly indicators of upsets whose Hasse subgraph and complementary Hasse
  subgraph are connected.  Hence in every bounded poset every nontrivial
  upset inequality is a distinct facet of the upward Hasse-boundary cone.
  For integer partitions, no universal coarsening-filter inequality in the
  THM-3127 dual description can be discarded.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
script: 04-computation/gmc_bounded_poset_upset_facets_thm3129.py
output: 05-knowledge/results/gmc_bounded_poset_upset_facets_thm3129.out
script_sha256: a52b93b38d870bce13fa9328f6a8fd880e71b4603ea5da706ff7dc11ee995c5b
output_sha256: b5b05962c01c6803ef581e37133dfcf05f0c736681c8b0d4270593a0de0b3a2b
hash_basis: LF-normalized bytes
---

# THM-3129 -- bounded-poset upset-facet irredundancy

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-3127 turns positive refinement transport into one inequality for every
coarsening upset.  A natural hope is that only principal upsets, connected
cuts, or some much smaller filter basis is needed.  The exact convex geometry
says otherwise: on a bounded poset every nontrivial upset is already an
extreme dual ray and hence defines a facet.  The filter family is large
because the refinement cone genuinely has that many faces, not because the
duality proof was wasteful.

## 1. The order cone and the Hasse cone

Let `P` be a finite poset whose undirected Hasse graph `Gamma(P)` is
connected.  Put

```text
I(P)={f:P->R : u<=v implies f(u)<=f(v)},                       (1)
```

and quotient by the line of constant functions.  The resulting cone
`Ibar(P)=I(P)/R1` is pointed.  Its dual is the upward Hasse-boundary cone

```text
C(P)=cone{e_v-e_u : u lessdot v}
     subset {c:sum_x c_x=0}.                                  (2)
```

The underlying graph is connected, so `(2)` spans the whole zero-mass
hyperplane and the two cones have dimension `|P|-1` in their respective
quotients.

For `U subset P`, let `Gamma[U]` denote the induced undirected Hasse
subgraph.  Call a nonempty proper upset `U` a **poset bond** if both

```text
Gamma[U] and Gamma[P\U] are connected.                         (3)
```

## 2. Extreme-ray classification

The extreme rays of `Ibar(P)` are exactly

```text
R_+[1_U]                                                       (4)
```

as `U` ranges over the poset bonds.

### Every extreme ray has two levels

Let an isotone `f` have distinct values `a_1<...<a_k`.  Its exact layer-cake
decomposition is

```text
f=a_1 1_P+sum_(j=2)^k(a_j-a_(j-1))1_{f>=a_j}.                 (5)
```

Modulo constants, every summand is the indicator of a nonempty proper upset.
If `k>=3`, `(5)` contains at least two distinct nested indicators with
positive coefficients, and their quotient classes are not proportional.
Thus an extreme ray has exactly two levels and is represented by `[1_U]` for
one nontrivial upset `U`.

### Disconnected sides split the ray

If `U` has connected components `U_1,...,U_s` with `s>=2`, then each `U_i`
is itself an upset.  Indeed, a saturated chain from `x in U_i` to any
`y>=x` remains in `U` and joins `x` to `y`.  Hence

```text
[1_U]=sum_i[1_(U_i)],                                         (6)
```

a nontrivial positive decomposition.

If the downset `D=P\U` has components `D_1,...,D_s`, each `P\D_i` is an
upset and

```text
sum_i 1_(P\D_i)=s 1_P-1_D,
[1_U]=sum_i[1_(P\D_i)]                                       (7)
```

modulo constants.  Thus disconnectedness on either side prevents extremality.

### Connected sides force extremality

Conversely suppose `(3)` holds and

```text
[1_U]=[f]+[g],                  f,g in I(P).                   (8)
```

Absorb a constant into `f` so that `(8)` holds pointwise.  Along every Hasse
edge internal to `U` or to `P\U`, the left side has increment zero.  Both
increments on the right are nonnegative, so each is zero.  Connectivity in
`(3)` makes `f` and `g` separately constant on each of the two sides.  Their
crossing increments are nonnegative and sum to one.  Therefore, modulo
constants,

```text
[f]=lambda[1_U],       [g]=(1-lambda)[1_U],    0<=lambda<=1.  (9)
```

This proves `(4)`.

## 3. Facets and bounded posets

Duality between the full-dimensional cones `(1)--(2)` turns every extreme
ray `(4)` into the facet inequality

```text
c(U)=sum_(u in U)c_u>=0.                                     (10)
```

Now suppose `P` has a unique minimum `0hat` and maximum `1hat`.  Every
nonempty upset contains `1hat`, and every vertex in it is joined to `1hat`
by a saturated chain lying inside the upset.  Every proper upset has a
nonempty downset complement containing `0hat`, with the analogous downward
chains.  Hence **every** nontrivial upset satisfies `(3)`.  Consequently:

```text
P bounded  ==>  every nontrivial upset inequality is a facet. (11)
```

In particular none of these inequalities is a nonnegative consequence of
the others.  Removing even one enlarges the universal Hasse-boundary cone.

## 4. Integer-partition refinement

The refinement poset `P_N` of integer partitions is bounded by

```text
0hat=(1^N),                 1hat=(N).                          (12)
```

Thus every nonempty proper coarsening upset is a distinct facet in the
THM-3127 Strassen/filter description.  The first branching case `N=4` has
five such facets.  One is the nonprincipal upset

```text
{(4),(3,1),(2,2)},                                          (13)
```

which supplies THM-3127's sharp hostile to principal-filter sufficiency.

The exact facet census through degree eleven is

```text
N       partitions    Hasse edges    upsets    nontrivial facets
2            2              1            3              1
3            3              2            4              2
4            5              5            7              5
5            7              9           10              8
6           11             17           27             25
7           15             28           47             45
8           22             47          168            166
9           30             73          573            571
10          42            114         3588           3586
11          56            170        19542          19540.    (14)
```

Thus a putative all-degree normalized filter proof faces `19,540` genuinely
distinct universal facets already at degree eleven.  Symmetry or special
identities of the product-Gamma response may still prove many facets at once,
but no facet can be deleted from the description of the ambient refinement
cone.

## 5. Exact companion and scope

Run

```text
python 04-computation/gmc_bounded_poset_upset_facets_thm3129.py
python -O 04-computation/gmc_bounded_poset_upset_facets_thm3129.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_bounded_poset_upset_facets_thm3129.out.
```

The companion constructs every partition Hasse graph and every upset in
`(14)`, verifies connectedness of both induced sides for all `23,949`
nontrivial facets, and hashes every exact upset bank.  It also checks the two
minimal nonbounded controls: a disconnected upset splits into its components,
while a disconnected complement splits modulo the constant line.  Every
load-bearing check remains active under optimization.

This theorem classifies the universal facet family.  It does not say that a
particular signed response visits every facet, or that every facet needs a
separate analytic argument.  It proves no new product-Gamma sign, arbitrary
support flag, Gaussian moment statement, LRC(14), JC(2), or DC(2).  Its sharp
message is a stopping theorem: principal filters and generic cut pruning
cannot replace the full coarsening-upset dual.

**QED (candidate pending independent audit).**
