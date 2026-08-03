---
id: THM-3307
title: "Aggregate origin-affine criterion and six minimal visible kernel circuits"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY CROSS-CHECKED.  In the finite static
  decorated relation-walk universe of THM-3287/3288, over THM-3277's selected
  eleven-edge tree and its eleven optional full-core edges, observer-visible
  adjacency zero modes admit an exact nine-scalar aggregate reduction.  With
  the edge indicators defined below, visibility is equivalent to `w=1` or to
  `w=0, u=1, 1-s=r(1-v), q-p+v != 0`.  This criterion gives exactly six
  inclusion-minimal optional-edge supports: one persistent singleton and five
  triples.  It also recovers their primitive kernel directions and exact
  masses `3/5`, four copies of `12/35`, and `3/7`.  This is a theorem about
  finite static symmetric compatibility graphs, not a tournament,
  chronological response system, or an FC/GMC/LRC mechanism.
audit: >
  The structural companion replays the pinned selector-paired half-transfer
  source, reconstructs the exact `14 x 12` selector block, verifies all twelve
  column equations and an explicit row-span certificate making the full-side
  observer invisible, checks aggregate lifting formulas, and derives the six
  supports from a written `r=0/r=1` split.  A 64-case reduced Boolean table is
  used only as a hostile check.  Exact Gram projections are computed only for
  the six derived supports.  The previous independent audit is hash-pinned but
  not executed; it separately reconstructed and exhausted all `2^11=2048`
  graphs from the older relation source and obtained the same six labelled
  primitive supports and masses.  Normal replay byte-matches the stored output;
  the structural source has zero assertion nodes and zero floating literals.
source: root/creative-synthesis-recover/2026-08-03
depends_on:
  - THM-3287-weighted-backbone-dominance-witness-section-and-selector-cut
  - THM-3288-maximizing-witness-lifted-walk-rational-series
  - THM-3305-single-edge-visible-zero-mode-and-rank-two-resolvent-update
related:
  - THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary
  - THM-3282-complete-boolean-seam-decoder-width-spectrum-and-maximal-cut-class
script: 04-computation/gmc_visible_kernel_aggregate_circuit_reduction_20260803.py
output: 05-knowledge/results/gmc_visible_kernel_aggregate_circuit_reduction_20260803.out
script_sha256: 46a168223455bc83d2c8031ef073eaffb81ef5dd3fd4faa6f8860e9953687f39
output_sha256: 83849f50ba9a2751b14db7de5800db1c475621db565d911aa5007023e9ad8cf2
independent_crosscheck_script: 04-computation/gmc_multi_edge_visible_kernel_bordered_resolvent_independent_audit_20260803.py
independent_crosscheck_output: 05-knowledge/results/gmc_multi_edge_visible_kernel_bordered_resolvent_independent_audit_20260803.out
independent_crosscheck_script_sha256: 75a7afacc5056976efff5686b3afe44534e8c64f1e2a9fee18d26b1c70f03e52
independent_crosscheck_output_sha256: ce7f9239b9ba68813ceb919662458b1c2ec32dea98d41507edab341a5417f651
hash_basis: LF-normalized bytes
---

# THM-3307 -- aggregate origin-affine criterion and six minimal visible kernel circuits

**PROVED + VERIFIED-EXACT + INDEPENDENTLY CROSS-CHECKED.**

This theorem replaces a `2^11` graph census by a short exact criterion.  The
criterion explains why the six visible circuits occur and separates the one
persistent circuit from the five fragile circuits.

## 1. Fixed universe and edge indicators

Let `T_*` be THM-3277's selected eleven-edge row tree and let `C` be
THM-3287's twenty-two-edge response core.  For a set `E` of optional row
edges, let `A_E` be the symmetric adjacency matrix of the decorated static
graph `Gamma(T_* union E)`, and let `1_E` be its active-node observer.  Define

```text
mu_0(E)=1_E^T P_ker(A_E) 1_E.
```

The eleven optional edges, in the fixed order inherited from `C`, receive
binary indicators as follows:

| edge | indicator | decorated addition |
|---|---:|---|
| `(2,10)` | `p` | `(2,Q+4)--(10,Q-1)` |
| `(3,22)` | `e` | `(22,Q+4)--(3,Q-1)` |
| `(7,18)` | `w` | `(18,Q+1)--(7,Q-7)` |
| `(10,16)` | `h` | `(16,Q+1)--(10,Q-8)` |
| `(10,22)` | `q` | `(22,Q+4)--(10,Q-1)` |
| `(11,13)` | `j` | `(11,Q+1)--(13,Q-8)` |
| `(13,18)` | `r` | all four `(18,Q+k)--(13,Q-1)`, `k=1,2,4,5` |
| `(13,22)` | `s` | `(22,Q+4)--(13,Q-1)` |
| `(16,21)` | `k` | `(16,Q+1)--(21,Q-8)` |
| `(17,22)` | `u` | `(22,Q+4)--(17,Q-1)` |
| `(19,22)` | `v` | `(22,Q+4)--(19,Q-1)` |

The selector bipartition of THM-3278/3287 splits the active full-core universe
into fourteen small-side and twelve full-side decorated nodes.  In this order,

```text
A_E = [ 0    N_E  ] .                                           (1)
      [ N_E^T 0   ]
```

All fourteen small-side nodes are already active in `Gamma(T_*)`.  The three
future full-side nodes `(10,Q-8)`, `(13,Q-8)`, `(21,Q-8)` have active-observer
coefficients `h,j,k`, respectively.

## 2. Exact aggregate reduction

Write a small-side coefficient vector in the labelled order

```text
x;
a1,a2,a4,a5;
b;
c1,c2,c4,c5;
d1,d2,d4,d5,
```

corresponding to rows `2,11,16,18,22`, and put

```text
A=a1+a2+a4+a5,   C=c1+c2+c4+c5,   D=d1+d2+d4+d5.             (2)
```

The twelve equations `N_E^T xi=0` reduce exactly to

```text
b=0,                    j a1=0,
d4=-x,                  a4=e x,
w c1+D=0,               A+(p-q)x=0,
(1-s)x+rC=0,            (1-u)x=0,
C+(1-v)x=0.                                                    (3)
```

The small-side observer on such a vector is

```text
S(xi)=x+A+C+D.                                                 (4)
```

No full-side kernel vector can supply a missing observer component.  Indeed,
if `E_(row,state)` denotes the corresponding row equation in `N_E z=0`, then
the active full-side observer is the following explicit row combination:

```text
Q-8 part     = E_(16,Q+1)+E_(11,Q+1)-E_(11,Q+2),
non-Q-8 part = E_(11,Q+4)+E_(22,Q+1)+E_(2,Q+4)
                -p E_(11,Q+2).                                (5)
```

Thus the full-side observer belongs to the row space of `N_E` and annihilates
`ker(N_E)`.

Conversely, every aggregate solution of `(3)` lifts to node coefficients.
One may take `a1=a5=c4=c5=d2=d5=0`,

```text
a4=e x,       a2=A-e x,       b=0,       d4=-x,
```

and then use either

```text
w=1: c1=-D, c2=C+D, d1=D+x,
w=0: D=0, c1=0, c2=C, d1=x.                                  (6)
```

The displayed lifts verify that the aggregate system loses no visible
solution.

## 3. Visibility criterion

If `w=1`, choose `x=A=C=0` and `D != 0` in `(6)`.  This gives a kernel vector
with observer sum `D`, for every choice of the other ten indicators.

Suppose `w=0`.  Equation `(3)` gives `D=0`.  If `x=0`, then `(3)` also gives
`A=C=0`, so `(4)` vanishes.  A visible solution therefore requires `x != 0`.
Dividing the remaining equations by `x` gives

```text
u=1,             1-s=r(1-v),             q-p+v != 0,          (7)
```

where the final inequality is exactly

```text
S(xi)/x=q-p+v.                                                 (8)
```

Equations `(6)` show sufficiency.  We have therefore proved the exact
classification

```text
mu_0(E)>0  iff
    w=1
    or [w=0, u=1, 1-s=r(1-v), q-p+v != 0].                    (9)
```

The four indicators `e,h,j,k` can change active nodes, kernel dimension,
projection mass, or recurrence order, but they cannot create first
visibility.  This is a structural irrelevance statement only for the Boolean
predicate `(9)`.

## 4. Six and only six minimal circuits

The minimal supports follow without a graph-lattice enumeration.

- If `r=0`, equation `(7)` forces `s=1`.  With `s,u` fixed, the smallest way
  to make `q-p+v` nonzero is to choose exactly one of `p,q,v`.
- If `r=1` and `v=0`, equation `(7)` forces `s=0`; exactly one of `p,q` is
  then required.
- If `r=1` and `v=1`, equation `(7)` forces `s=1`; the support already
  contains the visible proper subset `{s,u,v}` and is not minimal.

Together with the `w` singleton, this gives exactly

```text
{(7,18)},
{(2,10),(13,18),(17,22)},
{(2,10),(13,22),(17,22)},
{(10,22),(13,18),(17,22)},
{(10,22),(13,22),(17,22)},
{(13,22),(17,22),(19,22)}.                                   (10)
```

The earlier independent `2^11` audit obtains precisely the same list by a
different route.

## 5. Primitive directions and exact masses

Let `e_(r,Q+k)` denote the corresponding decorated coordinate, and abbreviate

```text
L11=e_(11,Q+1)+e_(11,Q+2)+e_(11,Q+5),
L18=e_(18,Q+1)+e_(18,Q+2)+e_(18,Q+4)+e_(18,Q+5),
L22=e_(22,Q+1)+e_(22,Q+2)+e_(22,Q+5).                         (11)
```

Exact projection onto `ker(N_E^T)` gives the following primitive direction
for each circuit.  Every other kernel direction is observer-orthogonal, so
the displayed mass is `(1^T z)^2/(z^T z)`.

| circuit type | primitive direction `z` | `(1^Tz,z^Tz)` | `mu_0` |
|---|---|---:|---:|
| `{w}` | `-3e_(18,Q+1)+(L18-e_(18,Q+1))+L22` | `(3,15)` | `3/5` |
| `{p,r,u}` or `{p,s,u}` | `-12e_(2,Q+4)+4L11+3L18-4L22+12e_(22,Q+4)` | `(12,420)` | `12/35` |
| `{q,r,u}` or `{q,s,u}` | `12e_(2,Q+4)+4L11-3L18+4L22-12e_(22,Q+4)` | `(12,420)` | `12/35` |
| `{s,u,v}` | `3e_(2,Q+4)+L22-3e_(22,Q+4)` | `(3,21)` | `3/7` |

The singleton direction is annihilated by every one of the ten remaining
optional-edge adjacency deltas.  Hence it remains a fixed kernel direction
on every superset containing `(7,18)`, and

```text
mu_0(E) >= 3/5 whenever (7,18) is in E.                       (12)
```

The five triples are fragile.  For example, starting from `{p,r,u}`, adding
`q`, `s`, or `v` violates `(7)` or `(8)` and removes visibility unless another
new circuit is simultaneously supplied.  This explains the previously
observed nonmonotonicity at the level of the aggregate equations.

## 6. Origin-affine interpretation and boundary

The small-side visibility condition is

```text
N_E^T xi=0,             1^T xi != 0.                          (13)
```

After normalizing `1^T xi=1`, equation `(13)` says that the origin lies in the
affine hull of the active small-row neighbourhood vectors of `N_E`.  The six
sets in `(10)` are therefore minimal *edge-parameter supports* that enable an
origin-affine circuit.  This is not a claim that the optional row edges are
the elements of an ordinary linear matroid.

The reduction proves a finite structural classification and supplies a cheap
compiler for visibility.  It does not turn the symmetric graph into a
tournament, recover chronological response dynamics, or prove an asymptotic
closed-form complexity bound.  Those losses remain explicit sidecars.
