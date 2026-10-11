# Boundary palette relations: exact compression for geometric gadget assembly

2026-10-11. **PROVED:** the elementary interface and saturation statements
below. **FINITE-EXACT:** the declared small-graph controls. **OPEN:** a supplied
six-chromatic unit-distance graph. This packet does not prove a new plane
chromatic bound or universal Collatz coverage.

## Inheritance and precise target

The [finite algebraic witness mechanism](finite_algebraic_witness_20261007.md)
already converts an accepted unrestricted plane lower bound into a terminating
search for a finite algebraic obstruction. The new
[critical witness packet](plane_algebraic_witness_20261011.md) retains the
geometry, noncolourability refutation, colouring, and vertex-deletion witnesses.
Here the target is smaller: compress a graph fragment by exactly what another
fragment can observe at its named boundary, and retain a lift to a full colouring.

The closest mechanism is palette permutation symmetry. The hostile is two
two-edge paths whose boundaries each use two colours but whose union is a
triangle. The corrected near miss is replacing an equality partition by its
number of blocks. The least-used sidecar is the complete relation of possible
boundary partitions, including evidence that omitted states are impossible.
The live board is **named port / equality partition / witness lift / complete
relation / geometric realization / terminal**.

## 1. The interface and its exact operations

Let F be a finite simple graph, B=(b_1,...,b_b) an ordered list of distinct
vertices, and k>=1 a fixed palette. Define

\[
 \Sigma_k(F,B)=\{\pi(c|_B): c\text{ is a proper }k\text{-colouring of }F\},
\]

where pi remembers which named ports have equal colours and forgets colour
names. Encode a partition as its restricted growth string, for example
(red,blue,red) becomes (0,1,0). Retain one full proper colouring for every
listed partition. There are at most sum_{j=0}^k S(b,j) possible states, with
the empty partition included when b=0. This is a bound on the interface,
not on the cost of discovering or certifying it.

**Gluing.** If F and H intersect exactly in their same ordered boundary B,
and their union has precisely their declared edges, then

\[
 \Sigma_k(F\cup_B H,B)=\Sigma_k(F,B)\cap\Sigma_k(H,B).       \tag{1}
\]

A colouring of the union restricts to both sides. Conversely, witnesses of
the same partition induce a bijection between their used boundary colours.
Extend this to a permutation of the full k-colour palette and apply it to
one witness. The witnesses now agree at every boundary vertex and glue.
This also constructs a full colouring lift for each surviving state. QED.

**Forgetting.** Removing a port from B, while leaving its graph vertex in
F, restricts and recanonicalizes each partition. This is existential
projection; duplicate projected states coalesce, retaining any one lift.
Later contexts may attach only through retained ports. Adding an edge to a
forgotten vertex is not an operation of this interface.

**Adding a boundary edge.** Adding b_i b_j filters out the states with
pi_i=pi_j. Their existing lifts remain valid for every retained state.

These operations can be composed recursively. Their proof is induction on
the operation DAG. A complete table can therefore be replaced by a checked
construction from complete leaf tables, with matching port names at every
node. The present prototype recomputes completeness directly on small
graphs; it is not a production implementation of that proof DAG.

**Completeness is separate from positive lifts.** A valid colouring for
every listed state does not prove that every possible state was listed.
Dropping a state can manufacture a false empty intersection. The checker
rejects missing states by recomputing the whole small relation.

## 2. A minimal hostile, and the geometric sidecar

On boundary (0,1,2), take F with edges 02,12 and H with edges 01,12. For k=2,

\[
 \Sigma_2(F,B)=\{(0,0,1)\},\qquad
 \Sigma_2(H,B)=\{(0,1,0)\}.
\]

Both have the same colour-count profile {2}, but their union is K_3 and
their relations have empty intersection. The three ports can be placed at
the corners of an equilateral unit triangle, so the loss is relevant even
for geometrically realizable subgraphs. The fragments here are subgraphs,
not induced graphs, of those three points.

For unit-distance assembly, retain exact coordinates or an exact realization
certificate and allowed isometries. Ports must occupy the same points;
distinct vertices must remain distinct. If the claim concerns the induced
unit-distance graph, account for every cross-unit edge. Additional unit
edges do not invalidate a lower bound obtained from a noncolourable subgraph,
but the purported subgraph still needs distinct, valid unit coordinates.
Equation (1) by itself does not certify geometric realizability.

A table is also palette-specific. With empty boundary, K_3 and the Moser
spindle both have empty two-colour tables. Their three-colour tables differ.
Thus one fixed-k table cannot in general answer the next-palette question.
No minimal-state theorem for geometrically restricted contexts is claimed.

## 3. A circle boundary that forces an entire palette

Let k>=2 and G be a k-vertex-critical graph: chi(G)=k and chi(G-v)=k-1 for every
vertex v. Set F=G-v and B=N_G(v). Then

\[
 \varnothing\ne\Sigma_{k-1}(F,B)
 \subseteq\{\pi: |\pi|=k-1\}.                         \tag{2}
\]

Otherwise a proper (k-1)-colouring of F leaves a colour unused on B and
extends to v, contradicting chi(G)=k. The star centred at v accepts exactly
the boundary partitions with at most k-2 blocks in this palette. Gluing it
to F therefore has empty relation. This gives a concrete construction target:
force the full palette on a circle, then add its centre.

For a unit-distance realization, every point of B lies on the unit circle
about v. Its local unit edges change angle by exactly +/-pi/3. Each connected
component is a subgraph of a six-cycle, hence is bipartite. Consequently, for
k>=4, the full-palette requirement in (2) must be imposed by the rest of F; it is not
explained by the local edges among circle points. This is a precise local
versus global distinction, without identifying the circle with a simplex.
The restriction k>=4 matters: in an equilateral triangle, the single edge
between the two neighbours already forces the full two-colour palette.

At palette five, saturation of a five-port boundary has only one equality
partition; six ports have S(6,5)=15; seven have S(7,5)=140. At palette six,
six ports have one and seven have 21. These counts classify equality data,
not tournament isomorphism classes. The port geometry and lift stay attached.

To seek a six-chromatic graph, construct a geometrically realizable F with
a circle boundary, a nonempty five-colour relation containing only five-block
partitions, and room for its distinct centre. The added centre then proves
non-five-colourability. An explicit six-colouring, or critical extraction,
supplies exact chromatic number six. In the alternative plane-seven case,
the corresponding target is a six-colour relation saturated on a circle.
Neither existence statement supplies a small boundary-degree bound.

### An equivalent finite extension target

For a fixed k>=1, the plane is k-colourable iff the following property holds:
for every finite k-colourable point set F and every p outside F, there is a
proper k-colouring of F that uses at most k-1 colours on
B={x in F: distance(x,p)=1}.

The forward implication restricts a global colouring, retaining the colour
of p as an unused boundary colour. Conversely, induct on the size of a
finite point set: colour its deletion F, use the property to recolour F
if necessary, and extend to p. Every finite point set is then k-colourable;
graph compactness supplies a global colouring. This is an exact reformulation,
not a shortcut proof of its universal premise. Recolouring F is essential:
the property does not require extending every existing colouring.

For k=6 this locates the case split at one interface. If the plane needs
six colours, every finite six-colourable fragment releases at least one
colour at each proposed centre. If it needs seven, deleting the centre of
a seven-critical witness gives a six-colourable fragment whose entire
boundary relation is saturated. The latter is a finite refutation; a finite
collection of successful extension tests cannot certify the universal case.

## 4. Exact controls and what is compressed

The [script](../../04-computation/experiments/boundary_palette_compression_20261011.py)
and [saved output](boundary_palette_compression_20261011.out) check:

- all 76 labelled simple graphs on 0..4 vertices, all boundary subsets, and
  palettes 1..3: 3,297 complete tables, compared with independent literal
  k^n enumeration rather than the producer's restricted-growth enumeration;
- every boundary-forgetting and boundary-edge operation in that universe;
- 192 pairs of two-port fragments with disjoint interiors, including full
  witness reconstruction after palette alignment;
- all seven circle centres of the exact algebraic Moser spindle: nonempty
  saturated three-colour deletion relations, bipartite local neighbourhoods,
  impossible three-colour star gluing, and possible four-colour gluing;
- the named-partition and palette-change hostiles, omitted states, forgotten
  ports, undeclared interior sharing, reordered ports, and malformed types.

The Moser profile (centre, degree, number of allowed saturated partitions) is
[(0,4,4),(1,3,1),(2,3,1),(3,3,1),(4,3,1),(5,3,1),(6,3,1)]. Its algebraic
coordinates use the inherited primitive element sqrt(3)+sqrt(11), isolated
in (5,6) with polynomial t^4-28t^2+64. This is a positive exact control for
the interface and saturation mechanism, not a six-chromatic witness.

The map is F -> Sigma_k(F,B) with a witness sidecar. It preserves extension
through declared ports under the specified operations. It destroys interior
geometry, colour names and multiplicities of lifts. Exact realization,
complete-table evidence, and source-specific reconstruction must be retained
when a later operation needs them. This is the useful connection to the
[Collatz proof-DAG routing packet](multiplication_certificate_routing_20261011.md):
both compress future-observable information while retaining a route back to
the literal witness. It does not turn graph colouring into Collatz descent.

Reproduce with:

```text
python -B -X utf8 04-computation/experiments/boundary_palette_compression_20261011.py
python -B -O -X utf8 04-computation/experiments/boundary_palette_compression_20261011.py
```

There are 14,884 always-active exact checks. No assertion disappears under
Python optimization. The cached producer validates public inputs before cache
lookup, retaining a hostile for bool/int equality aliases; see the corresponding
2026-10-11 entry in the [mistakes ledger](../../01-canon/MISTAKES.md).

An independent peer audit checked the interface, geometric hypotheses, and
optimized replay. Its small-k hostile produced the explicit k>=4 restriction
above and is retained in the final exact controls.
