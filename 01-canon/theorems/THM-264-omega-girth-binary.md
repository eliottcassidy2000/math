# THM-264: Finite-Exact Binary Girth of Tournament Conflict Graphs

**Status:** FINITE-EXACT for `3 <= n <= 6`; original global statement
**SUPERSEDED** by MISTAKE-507
**Filed by:** kind-pasteur-2026-03-21-S18b; scope repaired 2026-08-25

## Definitions

For a tournament `T`, let `Omega(T)` have one vertex for each distinct
directed odd cycle of `T`, taken up to cyclic rotation; distinct cycles with
the same vertex support remain distinct vertices.  Two cycle-vertices are
adjacent when their supports intersect.  Write
`alpha_1(T)=|V(Omega(T))|`, and use `girth=infinity` for an acyclic graph.

## Finite-exact statement

For every labelled tournament on `n` vertices with `3 <= n <= 6`,

- `girth(Omega(T))` is either `3` or `infinity`;
- `girth(Omega(T))=3` exactly when `alpha_1(T)>=3`; and
- `girth(Omega(T))=infinity` exactly when `alpha_1(T)<=2`.

In the same finite universe, the complement of `Omega(T)` is acyclic.

This is an exhaustive finite result, not a theorem for arbitrary order.

## Exact census

The computation exhausts all `2^binom(n,2)` labelled tournaments for
`n=3,4,5,6`, namely `8,64,1024,32768` objects.

- At `n=5`, 544 tournaments have conflict girth 3 and 480 have conflict
  girth infinity.  The former have `alpha_1 in {4,5,6,7}`; the latter have
  `alpha_1 in {0,1,2}`.
- At `n=6`, 28,848 tournaments have conflict girth 3 and 3,920 have conflict
  girth infinity.  The former have `alpha_1>=4`; the latter have
  `alpha_1<=2`.  The value `alpha_1=3` does not occur.

The structural observations in the original proof explain the small-order
census: every odd-cycle support has size at least three, so any two such
supports intersect when `n<=5` (and a five-cycle meets every nonempty
support); at order six the tournament constraints force a conflicting triple
once three or more odd cycles exist.  They do not extend the quantifier
beyond the computed range.

## Global failure boundary

The former implication

```text
alpha_1(T)>=3  =>  girth(Omega(T))=3
```

is false already at `n=7`.  THM-343 records ten labelled tournaments with
`alpha_1=3` whose conflict graph is `K_1 disjoint_union K_2`, hence has
girth infinity.  Thus the threshold criterion and the slogan “three cycles
force a conflict triangle” are **REFUTED** globally.

Whether every tournament conflict graph nevertheless has girth in
`{3,infinity}` is **OPEN** here: the order-seven witness refutes the proposed
criterion but does not create an intermediate cycle.

The finite complement statement does fail globally.  On nine vertices, take
three disjoint directed triangles `A,B,C` and orient every cross-block edge
according to the transitive order `A->B->C`, also `A->C`.  The strongly
connected components are exactly the triangles, so the only directed odd
cycles are those three disjoint cycles.  Hence

```text
Omega(T)=edgeless graph on three vertices,
complement(Omega(T))=K_3.                                  (1)
```

Thus global complement acyclicity is **REFUTED** (with a witness of order
nine); the finite `n<=6` statement remains exact.

Consequently, the old use of THM-264 as a universal exclusion of Petersen,
Moore, cage, or other finite-girth graphs is also superseded.  The exclusion
is proved only for conflict graphs arising in the finite universe
`3 <= n <= 6`.

## Reproduction

```bash
python3 04-computation/omega_girth_fixed_s18b.py \
  > /tmp/omega_girth_fixed_s18b.out
cmp /tmp/omega_girth_fixed_s18b.out \
  05-knowledge/results/omega_girth_fixed_s18b.out
```

The script includes positive controls `K_3,K_5`, hostile controls `C_5` and
Petersen (both girth five), the acyclic control `P_3`, and the order-nine
three-triangle complement hostile `(1)` before the tournament census.  Its
cycle generator fixes the least vertex as the first entry, which
removes cyclic rotations without collapsing multiple directed cycles on one
support.

## Correction lineage and related results

- MISTAKE-507 -- global-scope, Lie-carrier, flow, and arithmetic-analogy audit
- THM-343-H7-impossible -- exact order-seven hostile `K_1 disjoint_union K_2`
- THM-261-petersen-root-orthogonality -- exact Petersen/Kneser support carrier

## Source

`04-computation/omega_girth_fixed_s18b.py` and
`05-knowledge/results/omega_girth_fixed_s18b.out`
