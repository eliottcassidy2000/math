# Five-point interval walls, ternary prefixes, and what the preorder forgets

**Reflection / mechanism note, 2026-08-15.**  The truth source is
[THM-3462](../01-canon/theorems/THM-3462-five-point-line-metric-interval-sum-preorder-atlas.md),
which is **PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED**.  This note
records connections and boundaries; it is not an additional theorem.

## The object that finally stayed faithful

For five ordered line points, the natural variables are not the five points
and not a tournament on five vertices.  They are the four positive gaps and
the ten nonempty contiguous interval sums.  Their `45` pair comparisons split
into:

```text
25 strict containments + 20 mixed comparisons -> 15 equality walls.
```

The faithful finite state is the sign covector of those fifteen walls.  It
keeps every distance comparison and every tie, while forgetting all
magnitudes inside a cell.  This choice makes three ideas line up cleanly:

1. line geometry supplies interval indicators;
2. the wall arrangement supplies ternary signs `-,0,+`; and
3. the distance preorder supplies the consequence object.

The exact atlas has `477` states.  A bounded scan sees this number quickly but
does not explain why it is complete.  The decisive change was to certify
nonempty sign cones through pointed-polyhedron vertices and, independently,
through Fourier--Motzkin projection.

## Active concept board after the proof

| object / representation | preserved predicate | lost coordinate | cheapest next test |
|---|---|---|---|
| fifteen-wall covector | complete labelled distance preorder | magnitudes within a cell | derive an `n=6` wall generator before counting cells |
| ternary prefix feasibility tree | exact extendability of partial signs | geometric witness unless reattached | compare wall orders by live-prefix width |
| strict comparison tournament | edge order inside one chamber | all ties and metric values | measure which tournament statistics are constant on reversal classes |
| Fibonacci two-wall lane | `c=a+b`, `d=b+c`, and edge ties | Cassini/Pell unit and recurrence origin | classify generalized-Fibonacci rays inside `C054` |
| primitive representative profile | arithmetic complexity of each cell | real-cell volume and adjacency | correlate minimum total with face dimension, then attack hostiles |

This board separates exact carriers from attractive but lossy shadows.

## Connection 1: the four-point atlas extends, but its tiny chamber grammar does not

The source is THM-3457's three-gap/five-wall arrangement; the target is
THM-3462's four-gap/fifteen-wall arrangement.  The map adjoins one positive gap
and all intervals ending at the new endpoint.  It preserves the rule

```text
distance equality <=> equality of disjoint contiguous gap blocks after overlap cancellation.
```

What is destroyed by forgetting the new endpoint is every comparison using
`34,24,14,04`.  Conversely, the old six-edge preorder is the restriction of
the new ten-edge preorder, but many of the `477` states restrict to the same
one of the old `25` states.  The needed sidecar is therefore the new wall-sign
fiber, not a scalar count of added ties.

This is a real extension mechanism, but the growth `5 -> 15` walls and
`25 -> 477` cells warns against extrapolating a simple recurrence from the
first two atlas sizes.

## Connection 2: a ternary tree is present, but feasibility is the grammar

Each wall contributes one of three symbols.  In the frozen lexicographic wall
order, exact recursive feasibility has level sizes

```text
1,3,9,15,19,25,75,125,175,213,253,283,339,389,427,477.
```

Thus the atlas is literally the leaf set of a pruned ternary prefix tree.  The
tree is not the full `3^15` language, and its nodes are not Berggren branches,
Rule 30 states, or LRC owners.  The source-to-target contract is:

```text
partial ternary word
    -> margin-one linear feasibility problem
    -> exact extendability of a wall covector.
```

The map preserves existence of a positive gap realization.  It destroys the
realization itself; an exact Cramer witness restores it.  Reordering the walls
changes prefix counts and tree shape but not the final `477` leaves.  That
makes maximum live-prefix width a possible algorithmic invariant of an
elimination order, not an invariant of the line metric.

This also clarifies the relation to
[THM-3459](../01-canon/theorems/THM-3459-rule30-ternary-intersection-factorial-truth-lift-and-keller-boundaries.md):
both use ternary symbols as a defect/feasibility compiler, and neither obtains
orientation merely from the ternary alphabet.  THM-3459's symbols separate
owners in mask-intersection tests; THM-3462's symbols record signs of linear
forms.  Shared arity is not a carrier equivalence.

## Connection 3: XOR recovers `K5` incidence, not the distance order

Encode edge `ij` by `e_i+e_j` in `F_2^5`.  For two distinct `K5` edges,

```text
wt((e_i+e_j) XOR (e_k+e_l)) = 2  iff they share one endpoint,
                               4  iff they are disjoint.
```

So the ten comparison vertices form the weight-two layer of the Boolean
five-cube, and XOR recovers the line graph `L(K5)`.  This preserves incidence
and disjointness, including which mixed comparisons cancel to disjoint gap
blocks.  It loses every numerical inequality and supplies no direction.  The
fifteen wall signs are the missing orientation sidecar.

This is the exact survivor behind “XOR and the earlier concepts are a
tournament”: XOR supplies a symmetric incidence geometry; a strict metric
chamber separately supplies a transitive `T10`; a boundary supplies a
preorder.  There is no intrinsic tournament of size five carrying all three.

## Connection 4: Fibonacci occupies a codimension-two lane, not a family classifier

Four consecutive Fibonacci gaps satisfy

```text
c=a+b,  d=b+c,  a<b.
```

For stable indices they occupy `C054`, a projective one-dimensional stratum
with two tied edge pairs.  The boundary `(1,1,2,3)` adds `a=b` and lands at the
vertex `C057`.  But `(1,3,4,7)` has the same walls and the same preorder while
its Cassini form is `5`, not `+/-1`.

The map from Fibonacci windows to the atlas preserves the two recurrence
equalities and their induced ties.  It destroys the marked origin and Pell
unit.  Restoring those two sidecars is exactly what distinguishes the
Fibonacci ray from all positive generalized-Fibonacci seeds on the same lane.

## Connection 5: every atlas subfamily can index a harmonic subseries, but that is not yet structure

After assigning stable IDs `C000,...,C476`, any selected atlas family is a
finite subset of the natural numbers and hence a finite subseries of the
harmonic series.  Repeating or extending those IDs can manufacture infinite
subsets as well.  This observation preserves membership only.  It depends on
the arbitrary ID order and forgets wall adjacency, reversal, face dimension,
and metric meaning.

A meaningful harmonic transfer would need a canonical infinite tower—perhaps
the `n`-point atlases together with a specified endpoint-deletion map—and a
weight invariant under relabelling.  The cheapest falsification test for any
proposed density is to permute the stable IDs: if the harmonic claim changes,
it is about the encoding rather than the geometry.  THM-3462 alone supplies no
such density theorem.

## Why the two sharp bounds are conceptually useful

The total-23 hostile `(6,5,4,8)` and the height-10 hostile `(2,3,4,10)` are not
large random witnesses.  Each is forced by a short inequality chain:

```text
c<b<a<d, d<b+c, a+b<c+d       -> total >=23,
a<b<c<a+b, a+b+c<d             -> max gap >=10.
```

They show that representative complexity is controlled by narrow integer
crevasses between interval-sum walls.  This resembles an atom floor more than
a chamber count: the obstruction is a chain of strict integer separations.
For a future six-point atlas, searching such crevasses before attempting the
full census may expose lower bounds and hostile controls cheaply.

## Honest frontier

The package closes the five-point positional preorder problem.  It leaves open
the combinatorics and arithmetic growth for six or more points, adjacency and
volume data for the `477` cells, and any canonical deletion/extension
recurrence between atlas sizes.

It does not reconstruct points from an unlabelled distance multiset.  It does
not turn XOR into orientation, a ternary feasibility tree into a Berggren or
Rule 30 conjugacy, or an arbitrary subset of IDs into a harmonic invariant.
It supplies no LRC physical clock/current and no Keller or Jacobian map.
