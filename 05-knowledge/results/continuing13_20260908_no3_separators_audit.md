# Independent referee: cell-separator repair and a cycle-excess corollary

**Verdict: PASS, without mathematical repair.** The separator theorem,
weighted conditional compiler, connected geometric family, and stated
failure boundaries are valid. The complete finite replay agrees with the
producer. Section 5 adds a proved sufficient cycle-excess bound; it changes
no frozen primary claim and introduces no probabilistic or extremal result.

## 1. Exact target and independent routes

The target is continuing13_20260908_no3_cell_separators. Its frozen source
SHA256 is 76027e9d78bde4d140b921025115d3f8e6c2708a62f66a3a7ceae73fa4721c5e,
and its certificate SHA256 is
418f449018069ae13fbc7c5c121a3a80a518a952667104b012511b89cb23c271.
The producer has 45383 always-active gates. No mathematical producer is
imported or executed by this referee.

The referee reconstructs ordinary biconnected blocks by **recursive vertex
deletion**, splitting at an articulation until no articulation remains.
This is independent of the producer's Tarjan lowlink route. It then merges
blocks through line vertices, and verifies every original edge, whole
line ownership, separator membership, and the resulting forest.

For each conditional bag state, the first local route enumerates every
literal retained subset against the original owned-line constraints and
the actual continuation costs. A second route branches on its free local
high-incidence cells and enumerates integral residual graph edge sets with
their exact weights and capacities. These agree in every exercised parent
state. The assembled message optimum is separately compared with direct
whole-board deletion-subset enumeration. No triple-hitting recursion is
used by this whole-instance reference.

## 2. Analytic audit of ownership, messages, and compiler

The block-cut forest establishes the claimed bags. Contracting connected
groups through articulation *line* vertices assigns every original line,
with all of its cells and its full constraint, to exactly one bag. The
remaining intersections are articulation cells, and the bag-cell incidence
graph remains a forest. Empty and disconnected systems and bridge blocks
are included. Removing all articulation cells before forming bags is not
an equivalent construction; the actual cyclic hostile correctly retains
the lost cycle.

For a cell separator p, the message

    C_p(s)=1-s+sum_child_bags C_B(s)

pays its own deletion exactly once. A bag excludes the cost of its fixed
parent cell, uses child messages as unary costs, and keeps every owned
line's original retained-cardinality constraint <=2. A leaf-to-root
induction now partitions both variables and constraints exactly. Both
parent states remain feasible, since one may retain that parent and delete
every other descendant cell. There is no hidden infinity case in the
reward subtraction.

Writing baseline=sum c_p(0) and reward w_p=c_p(0)-c_p(1) is exact even
when rewards are negative. A fixed retained parent subtracts one unit of
capacity on **every** owned line through it, even if it has three or more
local incidences. It is not re-enumerated as a free high-incidence cell.
After fixing the remaining local high cells, every free cell is exactly
one ordinary edge or an edge with a private capacity-one leaf. Integral
weighted b-matching therefore preserves every feasible retained set and
its exact cost. Upper capacity constraints allow negative optional edges
to be omitted. Unweighted matching or bipartite flow would not preserve
this problem.

There are at most two parent states per nonroot bag and at most 2^h_B
free high-cell assignments in each state. Thus the theorem's bound
2 sum_B 2^h_B is correct. It counts graph oracle calls, with polynomial
message/decomposition overhead; the finite verifier is not itself a
polynomial-time matching implementation or a timing claim.

For three distinct directions every global exceptional cell has exactly
three incidences. An articulation cell splits these among at least two
bags, leaving at most two in each. For four or more directions this
conclusion fails, exactly as the retained fifteen-point example shows.
The general local-incidence definition h_B, rather than blanket removal
of articulation cells, is therefore essential.

## 3. Complete finite replay

The complete geometric bank is generated independently as all row-pair
products, followed only by exact column-degree-two tests:

| Square dimension | Complete simple two-regular boards | With global exceptional cells | With negative continuation reward |
| --- | ---: | ---: | ---: |
| 2 | 1 | 0 | 0 |
| 3 | 6 | 0 | 0 |
| 4 | 90 | 0 | 0 |
| 5 | 2040 | 16 | 12 |

All 2137 boards have h_B=0 for slopes 1,-1,2, and every localized optimum
equals direct original-cell deletion. Every aggregate count and graph-call
total agrees with the producer. No higher-incidence board was filtered out.

The referee also independently enumerates the whole declared abstract
linear triple-incidence bank on zero through four labelled lines, with
counts 1,1,2,9,96 (109 systems total). Shared cells have size two or three;
each pair of lines occurs together at most once, each line has shared
degree at most three, and private cells fill it to a triple. Every system
and every bag rooting is checked. These are abstract controls, not an
assertion that every one is geometrically realizable in the chosen slopes.

Additional actual-coordinate checks include all eight connected chains,
every inherited named board in three and four directions, the joint-choice
theta, the cyclic naive quotient, the four-direction boundary, and the
new saturated core. All root choices are exercised on these controls.

* The theta's forced **deletion** states 00,01,10,11 have exact costs
  3,4,4,2. Its two exceptional cells remain in one bag.
* The naive twelve-point construction has a genuine cycle in its
  delete-articulations-first quotient, but the correct quotient is a
  forest with four bags and exact cost three.
* The fifteen-point four-direction board has two articulation cells
  that each retain three incidences in the same bag. Its local counts
  are 0,0,2; its cost is two; all rerootings exercise four fixed-high
  parent states in total. Its negative reward -2 is retained.
* The six-by-six board (01,01,23,45,35,24) has one bag, h_B=1, and cost
  four, including its actual four-cell overfull line. Together with the
  complete n<=5 census, dimension six is minimal for positive h_B in
  the stated simple two-regular class and three fixed directions.
  No complete dimension-six enumeration is inferred.

## 4. The connected unbounded construction is analytically paid

For the centres (2i,2i+2 floor(i/2)), the labels in the three slopes are
2 floor(i/2), -2 ceil(i/2), and 4i+2 floor(i/2). Only adjacent centres
share a selected line, alternating slopes one and two. This gives exactly
2k+1 selected lines and a connected chain of k degree-three centres.

The private-point filling excludes only finitely many integer parameters
at each step: occupied rows and columns, other selected-line labels, and
previously used labels in the other two directions. Distinct nonzero
slopes make each such equality exclude at most one integer. Installing
all centres first and retaining all exclusions ensures the claimed
singleton nonselected lines and unique occupied rows and columns for
every k, not only the eight tested values.

The completed graph has 4k+3 cells, 2k+1 lines, and 6k+3 incidences,
so it is a connected tree. Its bags each own one complete line and
h_B=0. Deleting the k centres attains repair k; the k pairwise disjoint
slope-minus-one triples force at least k deletions. The uniform bound is
4k+2 graph calls and the declared rooting uses 4k+1. For k=1 the central
message is (1,2), giving reward -1 and showing why weighting is necessary.
These are finite integer boards with empty rows and columns in their
bounding squares, not saturated two-regular constructions or random-board
frequency estimates.

## 5. Additional proved corollary: cycle excess bounds the local parameter

Let K range over the ordinary **nonbridge biconnected blocks** inside a
merged bag B. Write beta(K)=|E(K)|-|V(K)|+1 for its cycle rank. Then

    h_B <= sum_(K inside B, nonbridge) (2 beta(K)-2).   (1)

First, any one cell belongs to at most one ordinary block inside a given
merged bag. If two blocks incident to that articulation cell were merged
through line articulation vertices, the connecting path avoiding the
cell, together with their two cell incidences, would make a cycle in the
block-cut forest. This is impossible. Hence a cell's local bag degree is
its degree within one ordinary block. A cell of local degree at least
three cannot belong to a bridge block.

Every vertex in a nonbridge biconnected block has degree at least two.
The number of its high-degree cell vertices is bounded by
sum_cell_vertices(deg-2), which is at most

    sum_all_vertices(deg-2)=2|E|-2|V|=2 beta(K)-2.

Summing proves (1). In particular, if the original bipartite incidence
graph is a **cactus** (every nonbridge block is a single cycle), every
h_B is zero for *any finite direction set*. Such instances require no
exceptional-cell branching; exact weighted integral graph optimization
inside the bags remains necessary.

This is sufficient, not an equivalence. The abstract linear triple system
obtained by subdividing every edge of K4 has four line vertices, six cells,
cycle rank three and cycle excess four, but every cell has degree two and
h_B=0. The referee checks this hostile. The inherited projected triangle
also reminds us that a cyclic bag with no exceptional cell can still
require nonbipartite integral matching. Equation (1) is a coarse structural
bound, not a replacement for the exact local-incidence parameter.

The referee checks the one-block-per-cell fact and (1) on every bag in
its complete finite universes and rerooted controls. The proof above,
rather than those tests, supplies the all-graph assertion.

## 6. Scope and frozen reproduction

The theorem is an exact original-cell repair reduction for arbitrary
finite boards and finite direction sets. If an adaptive process produces
a final no-three-in-line set A from an original board B, then B intersect A
satisfies every retained selected-direction constraint. Consequently
|B minus A| is at least the exact repair cost. This preserves original-cell
loss without claiming a bound on final cardinality under arbitrary
insertions. No new mean, concentration coefficient, saturated asymptotic
family, or extremal no-three-in-line constant is inferred.

The independent source passes 85491 always-active gates. Normal and
optimized runs produce identical raw LF output and certificate bytes.
Reproduce after relocation with

    python 04-computation/continuing13_20260908_no3_separators_audit.py
    python -O 04-computation/continuing13_20260908_no3_separators_audit.py

The runner accepts --producer for an external frozen certificate directory.
It writes its own certificate next to the report after repository
relocation or beside the source in the external packet. No maintained
repository file or frozen producer was edited by this independent audit.
