# Five private interiors cannot force the four paired palettes

**Status: PROVED.** Exact geometric argument below; finite independent controls
in `check_five_interiors.py`. This strengthens the four-interior result in
`../root_ports_closed_carriers_20261028/GEOMETRIC_SELECTORS.md`.

Let eight distinct points of the Euclidean plane be partitioned into four
disjoint named pairs. Prescribe a proper colouring of the unit-distance graph
on these ports using at most six colours, constant on each pair. Add at most
five other distinct points, the private interiors. Then the prescribed port
colouring extends to a proper six-colouring of the full unit-distance graph.
Extra unit edges are retained throughout. There are no coincidences, reused
ports, assumed palettes, or abstract edges standing in for distances.

## The exact critical quotient

Merge port pairs with the same prescribed colour into colour groups. There
are at most four groups. Add all edges between these groups to preserve their
distinct colours, and join a group to an interior whenever at least one port
in it has that unit contact. The resulting graph has at most nine vertices.
A six-colouring of it extends the prescription after a permutation of colour
names. Conversely an extension gives such a colouring.

Every graph of chromatic number at least seven has an induced seven-critical
subgraph: delete vertices until chromatic number first becomes seven, then
minimize. Its minimum degree is at least six.

**Lemma.** A seven-critical graph on at most nine vertices is either `K7` or
`K4 join C5`.

For seven vertices this is immediate. For eight vertices its complement is a
matching with isolates; any matching edge violates vertex-criticality. For
nine vertices its complement has maximum degree two, hence is a disjoint
union of paths and cycles. Colouring the original graph is clique-covering
the complement. The saving from the singleton cover is the matching number,
plus one additional unit for each triangle component. Vertex-criticality
requires this saving to remain unchanged when **any** vertex is deleted.
A nontrivial path, even cycle, or triangle fails this test. The surviving
components are isolates and odd cycles of length at least five. Total saving
is `9-7=2`, so there is exactly one `C5` and four isolates. This proves the
lemma without a general critical-graph classification theorem.

A `K7` would contain at least three interior vertices. Four interior vertices
cannot form a unit clique in the plane. With exactly three, they form an
equilateral unit triangle and the other four vertices must be four separate
port pairs. Each pair must touch all three triangle vertices, so at least one
member touches two. There are exactly three external double-contact points
of a unit triangle; four disjoint pairs cannot use them. Thus `K7` is excluded.

In `K4 join C5` there must be four different boundary groups, each consisting
of one original pair, and all five interiors. The boundary groups form a
four-clique. Write `k` for how many belong to the universal `K4`. A cycle
clique has size at most two, so `k=2,3,4`.

* `k=2`: the interiors include two universal vertices and a three-vertex
  path in the cycle. They contain a unit `K4`, impossible.
* `k=3`: the interiors form `K1 join P4`. The four path vertices lie around
  the hub as four consecutive vertices of a regular unit hexagon. An
  external point cannot touch three of these five interior points: three
  outer points have only the hub as a unit-circle centre; a triple consisting
  of the hub and two outer points either has no unit-circle centre or has one
  of the existing outer points as its centre. Each universal boundary pair
  must cover all five interiors, requiring a new triple-contact point. This
  is impossible even for one such pair.
* `k=4`: all interiors form an **induced unit pentagon**. Each of the four
  boundary pairs must cover its five vertices. The next lemma supplies the
  obstruction.

## The induced-pentagon contact lemma

**Lemma.** An induced unit `C5` has no external point touching four vertices,
and has at most three external points touching three vertices.

For four contacts, the four selected cycle vertices induce `P4`. Around the
external point they must be four consecutive vertices of a unit hexagon:
consecutive angular steps are `±60 degrees`, and changing direction repeats
a vertex. Their endpoints are antipodal. The fifth cycle vertex, adjacent to
both endpoints, is forced to equal their midpoint, the purported external
point. This contradicts distinctness.

For three contacts there are two combinatorial types of selected triples.
If they induce an edge and an isolate, translate their common centre to zero
and call their unit vectors `a,b,c`, with `|a-b|=1`. In cyclic order the whole
pentagon must be

```
a, b, b+c, c, a+c.
```

Indeed the two common unit neighbours of `b,c` are `0,b+c`; the first is
excluded, and likewise for `a,c`. Degeneracies only force repeated vertices.
But `|(b+c)-(a+c)|=|b-a|=1` is a diagonal of the pentagon. This contradicts
inducedness. This is the missing local geometric constraint that a generic
critical-graph model does not see.

Thus every external triple-contact point belongs to one of the five
consecutive triples. Such a triple has at most one unit-circle centre, and
its middle vertex has turning angle `±60 degrees`. Four different centres
would impose this condition at four of the five vertices. All five edge
directions would then lie in the same six triangular-lattice directions.
Distinct unit lattice edges have no crossings, so this is a simple lattice
pentagon. Inducedness forbids turns of `±120 degrees` (they create a unit
diagonal); reversal repeats a vertex. Each remaining turn is `0` or
`±60 degrees`. Five such turns cannot sum to the required `±360 degrees`.
Contradiction. There are at most three centres.

Finally each pair covering five cycle vertices needs a member with at least
three contacts. Four disjoint pairs require four distinct such points,
contradicting the lemma. This closes the `k=4` case and the extension theorem.

## Boundary and next question

This is a lower bound on this **specific paired-port selector architecture**,
not a lower or upper bound on the chromatic number of the plane. A selector
obstructing one of these assignments needs at least six private interiors.
The contact lemma does not extend unchanged to six-cycles: a regular unit
hexagon has a common unit-contact centre. That is a useful next geometric
object, but one centre still cannot be reused by distinct ports.

`check_five_interiors.py` independently enumerates all possible maximum-
degree-two complement types in the critical quotient, computes the fan's
circle centres over rational triangular coordinates, and checks every rooted
unit triangular-lattice pentagon. These checks support the proof, not an
unbounded search claim.
