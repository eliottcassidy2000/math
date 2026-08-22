# Multi-edge visible-kernel arity and a bordered resolvent compiler

**Status:** **FINITE-EXACT** procedural synthesis, with no theorem ID.  The
primary computation is
[the multi-edge scout](../04-computation/gmc_multi_edge_visible_kernel_helly_bordered_resolvent_scout_20260803.py)
and its
[stored output](../05-knowledge/results/gmc_multi_edge_visible_kernel_helly_bordered_resolvent_scout_20260803.out).
It replays the proved single-edge result
[THM-3305](../01-canon/theorems/THM-3305-single-edge-visible-zero-mode-and-rank-two-resolvent-update.md)
by pinned hash.

This note concerns static symmetric relation-walk graphs.  It proves no
tournament, chronology, response composition, `GMC`, `FC`, or `LRC` claim.
The phrase *visible-kernel Helly arity* below names an exact minimal-subset
census; it does not invoke or prove an abstract Helly theorem.

## Inheritance pass

The closest proved mechanism is THM-3305: among the eleven missing row
relations over THM-3277's selected tree, the single addition `(7,18)` uniquely
creates an all-ones-visible adjacency zero mode.  Its primitive incidence
direction is a local three-versus-six relation, its zero mass is `3/5`, and
its scalar walk series is an exact rank-two Woodbury update.

The canonical hostile example remains the relation path `2-17-16`: these
graphs encode symmetric static compatibility, not a chronological response
law.  The corrected near miss is to retain only a recurrence denominator.
A visible zero eigenvalue supplies no nontrivial resolvent pole in `x`; it
instead changes the degree-zero spectral mass and the finite recurrence head.
The least-used sidecars restored here are therefore:

```text
active-node births; the kernel projection itself; the edge subset supporting
it; and the local update endpoints before scalar recurrence elimination.
```

The live board was:

```text
edge subset; incidence circuit; all-ones visibility; persistence under the
next edge; active-node birth; local resolvent border; finite recurrence head.
```

## Exact closure under all edge subsets

The computation first exhausts all `55` unordered two-edge additions and then
all `2^11=2048` subsets of the missing core edges.  For every subset it rebuilds
the variable active-node graph exactly and records

```text
(subset size, active nodes, directed arrows, kernel dimension,
 all-ones zero mass).
```

The serialized lattice has `37,894` bytes and SHA-256

```text
53f740b8fd3e8b1baa88a460e210662d51b6aa69ff8a157fcdfc9fd5d8510ad0.
```

The visibility census by subset size is:

| added edges | subsets | zero mass | nonzero mass | distinct nonzero masses |
|---:|---:|---:|---:|---:|
| 0 | 1 | 1 | 0 | 0 |
| 1 | 11 | 10 | 1 | 1 |
| 2 | 55 | 45 | 10 | 1 |
| 3 | 165 | 115 | 50 | 3 |
| 4 | 330 | 188 | 142 | 10 |
| 5 | 462 | 212 | 250 | 20 |
| 6 | 462 | 169 | 293 | 25 |
| 7 | 330 | 91 | 239 | 23 |
| 8 | 165 | 29 | 136 | 18 |
| 9 | 55 | 4 | 51 | 13 |
| 10 | 11 | 0 | 11 | 7 |
| 11 | 1 | 0 | 1 | 1 |

All ten visible pairs merely inherit `(7,18)`.  There is no new
inclusion-minimal two-edge circuit.  At arity three, however, five new minimal
visible sets appear:

```text
{(2,10),(13,18),(17,22)}       mass 12/35
{(2,10),(13,22),(17,22)}       mass 12/35
{(10,22),(13,18),(17,22)}      mass 12/35
{(10,22),(13,22),(17,22)}      mass 12/35
{(13,22),(17,22),(19,22)}      mass  3/7.
```

Together with the singleton `{(7,18)}`, these are exactly the six
inclusion-minimal visible sets.  Thus the literal minimality profile is

```text
arity one: 1,  arity two: 0,  arity three: 5.                 (1)
```

Pairwise invisibility is therefore not enough to certify joint invisibility.
The missing coordinate is an edge-set circuit, not an individual edge score.

## Persistent and fragile circuits are different objects

For a graph `A` and all-ones vector `1`, write

```text
mu(A)=1^T P_ker(A) 1.
```

The singleton circuit has a primitive vector `v` with

```text
1^T v=3,   v^T v=15,   mu=3/5.                              (2)
```

Every one of the ten remaining edge deltas annihilates this same `v`.
Consequently `v` remains in the kernel of every graph above `{(7,18)}` in
the subset lattice, and orthogonal projection gives the exact lower bound

```text
mu(A_later) >= (1^T v)^2/(v^T v)=3/5.                       (3)
```

The census verifies `(3)` on all `1024` supersets.  This is the persistent
circuit.

The five arity-three circuits are fragile.  Four have primitive statistics
`(1^T v,v^T v)=(12,420)` and one has `(3,21)`.  For the first triple, adding
one more edge changes the mass as follows:

| fourth edge | zero mass |
|---|---:|
| `(3,22)` | `12/59` |
| `(7,18)` | `14/15` |
| `(10,16)` | `12/35` |
| `(10,22)` | `0` |
| `(11,13)` | `12/37` |
| `(13,22)` | `0` |
| `(16,21)` | `12/35` |
| `(19,22)` | `0` |

So visible-kernel mass is not monotone under edge addition.  An added edge
may preserve the circuit, rotate the full kernel projection, replace it by a
different visible circuit, or annihilate visibility entirely.  In particular,
adding `(7,18)` breaks this triple's primitive direction but introduces the
persistent singleton direction; the resulting mass `14/15` must not be
misread as preservation of the old circuit.

This yields a more useful sidecar than a Boolean visible/invisible flag:

```text
(minimal support, primitive direction, compatible next edges,
 breaking next edges, observer mass).                         (4)
```

## The bordered closed-form compiler

The full-core graph has `26` possible active vertices.  Embed the selected
tree adjacency as a fixed `26 x 26` matrix `B`, padding future vertices as
isolates.  Let `w_0` mark the old active vertices and `w` the active vertices
after an update.  If `b` vertices are born, then

```text
w^T (I-xB)^(-1) w = w_0^T (I-xB)^(-1) w_0 + b.              (5)
```

For the new undirected decorated edges, choose endpoint columns `U` and a
block-diagonal swap matrix `C`, so

```text
Delta=U C U^T,   C^(-1)=C.
```

The scalar generating series therefore satisfies the exact identity

```text
G_(B+Delta,w)(x)-G_(B,w_0)(x)
  = b + x s(x)^T (C-xK(x))^(-1) s(x),                       (6)

s(x)=U^T(I-xB)^(-1)w,
K(x)=U^T(I-xB)^(-1)U.
```

Equation `(6)` separates three things that a raw recurrence hides: observer
births, the fixed base propagation, and a small local inverse.

It can be compiled coefficientwise without a symbolic `26 x 26` inverse.  If

```text
s(x)=sum s_n x^n,  K(x)=sum K_n x^n,
(C-xK(x))^(-1)=sum Q_n x^n,
```

then

```text
s_n=U^T B^n w,      K_n=U^T B^n U,
Q_0=C,
Q_n=C sum_(i=1)^n K_(i-1) Q_(n-i),                          (7)
```

and the correction to coefficient `n>=1` is

```text
sum_(i+j+k=n-1) s_i^T Q_j s_k.                              (8)
```

Across all `55` two-edge updates, the local dimension is only `4` or `10`:

```text
(new decorated edges, local dimension, births): count
(2,4,0):21  (2,4,1):21  (2,4,2):3  (5,10,0):7  (5,10,1):3.
```

The first `20` coefficients from `(5)--(8)` agree exactly with direct matrix
powers for every pair.  This is an efficient exact recurrence compiler in a
fixed finite universe.  It is not yet a bit-complexity theorem, and it does
not claim that the local dimension stays bounded in an asymptotic family.

## Cross-frontier synthesis without type collapse

The arity profile `(1)` is procedurally analogous to the factorial
three-support-face obstruction: every lower-order view can look compatible
while a three-way closure exposes a new object.  The preserved predicate is
only this:

```text
test minimal joint supports, not just individual or pairwise shadows.
```

The source and target objects are different.  Here a circuit is a literal
linear dependence in a symmetric adjacency kernel.  In the factorial lane,
the obstruction is a typed ancestry/reset condition under row-labelled
successors.  There is no map identifying their states, and no implication in
either direction.

The Jacobian comparison is likewise procedural.  A resultant is a scalar
shadow of a coefficient ideal; a recurrence is a scalar shadow of a bordered
resolvent.  In both cases retaining the local presentation distinguishes a
persistent component from a fragile specialization.  But an exceptional
divisor and an adjacency-kernel circuit are not the same mathematical object.

Tournament language should not be forced here.  One can form a directed
compatibility relation “edge `e` breaks circuit `v`,” but many pairs tie and
the underlying graph remains symmetric.  The correct next object is a
bipartite circuit--edge incidence structure, not a cosmetic tournament.

## What changed on the live board

1. **Edge subset.**  The first genuinely new visibility after the singleton
   occurs at arity three, not two.
2. **Incidence circuit.**  Minimal support and primitive direction are the
   correct pre-scalar data.
3. **Visibility.**  Zero mass is nonmonotone and cannot label an edge once and
   for all.
4. **Persistence.**  Delta-annihilation of the primitive direction gives an
   exact upper-interval certificate and lower bound.
5. **Births.**  Variable active-node boundaries contribute the explicit
   constant `b` in `(6)`.
6. **Closed form.**  The sequence update reduces to base powers plus a local
   formal inverse of dimension at most ten on every tested pair.
7. **Finite head.**  The visible zero mode lives in numerator/initial-data
   information even though zero supplies no nonconstant resolvent pole.

## Next decisive tests

The cheapest positive continuation is not another scalar recurrence fit.
It is to derive the six minimal circuits directly from the decorated
incidence columns, explaining why the four `12/35` directions occur in a
two-by-two switch family and why the `3/7` direction has smaller support.
That would replace the `2048`-subset discovery census by a structural circuit
enumerator.

Then extend `(6)--(8)` from pairs to all subsets, quotient first by any
persistent circuit directions, and record determinant/numerator updates of
the local border.  The hostile checks are essential: preserve active-node
births, verify direct coefficients, and test whether fragile circuit factors
cancel after a breaking edge.  Only after that closure should one ask for an
asymptotic closed-form or complexity theorem.
