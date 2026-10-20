# Exact palette forcing outside the triangular lattice

**PROVED, with FINITE-EXACT geometric certificates:** a 16-point unit-distance
graph with nine properly precoloured ports forces a distinguished interior
vertex into the three colours used on those ports, complementary to the
allowed palette of the other six interior vertices. Adding two private ports
gives an 18-point graph whose *entire proper six-colour boundary relation is
full*. Adding a third gives a 19-point graph with an explicit nonextending
proper boundary colouring. All three graphs have chromatic number four.

The positive result is a boundary relation, not a six-chromatic construction.
The remaining operation is to realize the input palette conditions using
unit-distance geometry. The general list-colouring mechanism is classical;
the concrete port arrangement and its exact contact audit are the contribution
of this packet. No minimality among all possible forcing gadgets is claimed.

## 1. Inheritance and the live concept board

The closest proved mechanism is the wheel-and-three-wings relation in
[adaptive_budgets_smooth_interfaces_20261020, §3](../adaptive_budgets_smooth_interfaces_20261020/README.md):
its single interior vertex has `6 - number of distinct rim colours`
extensions. Its blocked states are precisely rainbow rims. The corrected
near miss is that this graph, and unit-edge gluings remaining in the same
triangular lattice, are still three-colourable.

The canonical loss example is the triangle-fibre strip in
[prime_batch_triangle_20261020, §6](../prime_batch_triangle_20261020/README.md):
the three-colour orientation obstruction has already vanished at four
colours after two links. It does not supply a six-colour palette synchronizer.

The less-used relevant sidecar is degree choosability. The repository already
uses the full Gallai-tree obstruction in
[THM-4661, odd-degree carriers](../../../01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md).
The short diamond proof below needs none of that theorem's field hypotheses.
Conversely, its carrier-specific size bounds cannot be transferred to the
new geometry. Likewise,
[THM-4614, three-colourable contexts in Q(zeta55)](../../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md)
has a fixed port set and field; it is not a universal assertion about every
three-colourable exterior and every choice of ports.

The board is: full boundary lists; degree-list reducibility; rotated spindle
geometry; port covers by matchings; and the common-anchor obstruction. The
anchor is six-colour boundary forcing, the niche is an explicit list kernel,
and the wildcard is turning three perfect matchings into exact geometric
palette incidences. A connection must retain the actual boundary graph and
the colours at its named vertices.

## 2. A boundary extension certificate

Let G be a finite simple graph, B a set of boundary vertices, H=G-B its
interior, and c a proper q-colouring of the **actual induced graph G[B]**.
For v in H put

    L(v) = {0,...,q-1} minus {c(w): w in B and vw in E(G)}.

The boundary colouring extends if and only if H is L-colourable. The boundary
edges are not discarded by this equivalence: their properness is a premise.
If every interior vertex has total degree at most q, then

    |L(v)| >= q - d_B(v) >= d_H(v).                         (1)

**Diamond certificate.** If H is connected and contains an induced diamond
K4 minus one edge, every list assignment satisfying (1) is colourable.

Proof. Colour vertices outside the diamond in decreasing distance from it.
Each has a neighbour of smaller distance still uncoloured, so at most
d_H(v)-1 colours have been forbidden by already coloured neighbours.
At the end, the residual diamond lists have sizes at least 2,2,3,3, with
the two size-two positions being the nonadjacent tips. Restrict to those
sizes. If the tip lists intersect, give the tips one common colour; both
adjacent middle vertices retain at least two colours and can differ.
If the tip lists are disjoint, their union has size four. After assigning
tip colours a,b, failure at the middle edge can occur only if both middle
lists were the same three-set {a,b,c}. If the middle lists are different,
this never happens. If they agree on a three-set S, select a tip colour
outside S, possible because the tip-list union has size four. Both middle
lists then retain at least two colours. This completes the colouring. QED.

This is a special case of the classical degree-choosability theory already
inherited through THM-4661, with a self-contained proof adapted to this
boundary problem. It gives a useful hostile test: with this interior, a
blocked q-colour boundary must give some interior vertex total degree >q.

For comparison, two adjacent interior vertices with available sets A,B
have exactly `|A||B|-|A intersect B|` extensions. If both sets are nonempty,
failure is exactly the equal-singleton case. Coupled lists matter even when
each vertex separately has an available colour. This elementary identity
does not assert any new geometric realization.

## 3. Seven core vertices and nine palette ports

Use the real coordinate field K=Q(sqrt(3),sqrt(11)), and the following
distinct points:

| i | p_i |
|---|---|
| 0 | (0,0) |
| 1 | (1,0) |
| 2 | (1/2,sqrt(3)/2) |
| 3 | (3/2,sqrt(3)/2) |
| 4 | (5/6,sqrt(11)/6) |
| 5 | ((5-sqrt(33))/12,(5sqrt(3)+sqrt(11))/12) |
| 6 | ((15-sqrt(33))/12,(5sqrt(3)+3sqrt(11))/12) |

Their unit-distance graph is the Moser spindle, with edge set

    01,02,04,05,12,13,23,36,45,46,56.                      (2)

The two diamonds share p0; the second is rotated by the complex unit
`(5+i sqrt(11))/6`. In any three-colouring, the first diamond forces
c(0)=c(3), the second forces c(0)=c(6), and edge 36 makes this impossible.
This classical spindle argument is inherited, not claimed as new; see
de Grey's [2018 paper, §1](https://arxiv.org/html/1804.02385v3).

If k is a common unit neighbour of i,j, the point `p_i+p_j-p_k` is the
reflection of p_k through the midpoint of p_i,p_j and is another common
unit neighbour. Add nine such points, numbered in row order:

| Boundary vertices | Fixed colour | Triples (i,j,k), defining p_i+p_j-p_k |
|---|---:|---|
| 7,8,9 | 3 | (1,3,2), (2,4,0), (5,6,4) |
| 10,11,12 | 4 | (1,4,0), (2,6,3), (3,5,6) |
| 13,14,15 | 5 | (1,5,0), (2,3,1), (4,6,5) |

Each row's first two entries form a perfect matching of {1,...,6}. The
exact pairwise-distance audit proves all nine points are distinct, outside
the core, and adjacent to exactly their displayed two core vertices. No
new point sees p0. There are exactly two boundary edges, **8--10 and
10--13**; they join different fixed colours. Thus this boundary colouring
is proper, but the boundary is not an independent set. The induced
16-point unit-distance graph has 31 edges: 11 core, 18 core-port, 2 port-port.

**Palette gate.** Under the displayed nine-port colouring, the six noncentral
core vertices have lists {0,1,2}. The output p0 can take exactly {3,4,5}.

If p0 were in {0,1,2}, (2) would be three-coloured, impossible. For each
choice of p0 outside that palette, triangles 123 and 456 each receive
{0,1,2}, with the sole additional restriction c(3)!=c(6). There are
`3*2*2*2=24` assignments. The exact output multiplicities at colours
0,...,5 are therefore `(0,0,0,24,24,24)`.

The map is now explicit: three geometrically realized matching covers
carry three named excluded colours to six interior vertices. The preserved
predicate is this extension relation. Forgetting the port colours destroys
it; identifying equal-coloured ports changes the geometry and is not a
permitted quotient without a separate realization proof.

## 4. The sharp degree boundary: eighteen versus nineteen points

Add three more points on the unit circle around p0:

    p16=(0,1), p17=(-3/5,4/5), p18=(-4/5,3/5).

They come from the rational parametrization
`((1-t^2)/(1+t^2),2t/(1+t^2))`, at t=1,2,3. Exact checking of every point
pair proves that their only new edges are 0--16, 0--17, 0--18. They do
not coincide with earlier points or create extra port-port contacts.

**Full eighteen-point relation.** Delete any one of these three private
central ports, and regard the remaining eleven ports as boundary. Every
proper six-colouring of the induced boundary graph extends to the seven
core vertices. Indeed, all seven interior vertices now have total degree
six. The spindle is connected and contains an induced diamond, so (1)
and the diamond certificate apply. This is an all-boundary-colourings
proof, not a conclusion extrapolated from the displayed colouring.

**Blocked nineteen-point assignment.** Keep all three private ports and
give them colours 3,4,5 respectively, in addition to the nine colours
above. The entire twelve-port colouring is proper. Every core list is
now {0,1,2}, so no extension exists. The interior degree vector is
`(7,6,6,6,6,6,6)`: one vertex has crossed the degree certificate's bound.

Every one of the twelve port constraints is necessary for this particular
blocked witness: deleting any one port makes its restricted boundary
colouring extend. The checker records a colouring for every deletion.
The universal full-relation assertion above is made only for deletion of
one of the three **central** ports. The other nine deletion checks are
positive controls for the one displayed assignment.

The full nineteen-point graph is four-colourable: one verified colour
vector in vertex order 0,...,18 is

    (0,1,2,0,1,2,3,2,0,0,2,0,1,0,1,0,1,1,1).

Together with the spindle obstruction, this proves its chromatic number
is exactly four, and the same holds for its 16- and 18-point subgraphs.
A proper but nonextending six-colour boundary assignment does not imply
chromatic number six. Conversely, proving a whole graph not six-colourable
would force at least seven colours; a six-colour lower bound requires
failure of five-colourability.

## 5. A small exact geometric minimum, with its quantifiers retained

For this fixed embedded core, no point outside the seven core points is
at unit distance from three core vertices. The certificate examines all
35 triples. None is collinear. For 25 triples the circumradius is not
one. For the remaining ten, the unit circumcentre is already a core
vertex. A common unit neighbour of a noncollinear triple must be that
unique circumcentre, so this is exhaustive over the entire real plane,
not merely over K-rational candidate points.

Consequently any external port can directly forbid its one fixed colour
at at most two core vertices. To forbid the same three named colours
at all six noncentral targets requires at least three ports per colour,
hence nine. To forbid them at all seven targets requires at least four
ports per colour, hence twelve. Our arrangements attain both bounds while
respecting their induced boundary edges.

This is a minimum only for **direct neighbour exclusion of three fixed
colours at these targets in this embedding**. It is not a lower bound for
arbitrary internal palette gadgets, different core embeddings, indirect
forcing, or all unit-distance graphs.

## 6. The next forcing problem and the obstruction to the easy shortcut

The gate now has an exact useful interface, but its nine input colours
are prescribed externally. Replacing each colour class of ports by one
common anchor vertex is not a geometric solution. Two distinct fixed
anchors have at most two common unit neighbours, whereas all six target
vertices would have to be common neighbours of every pair of the three
anchors. Thus the naive three-anchor replacement is impossible in the
plane. Allowing different coordinates or extra edges does not evade that
two-circle intersection bound while retaining the six distinct targets.

The strongest next target is a unit-distance **palette synchronizer**:
an actual geometric assembly whose boundary relation supplies these
three monochromatic triples of distinct colours, or an alternate input
relation strong enough to force the same output. Its test must retain
named ports, actual unit contacts, all permissible input palettes, and a
proper global colouring as a hostile control. A gluing based only on
abstract equalities silently removes the central difficulty.

This is consistent with the existing THM-4661 equality-tester warning:
forcing equality in every proper q-colouring while remaining q-colourable
already requires chromatic number q. The present graph evades no such
barrier because its equality pattern is a boundary condition.

## 7. Reproduction and validity scope

Run from the repository root:

    python 05-knowledge/results/kernels_games_detours_20261020f/forcing_check.py
    python -O 05-knowledge/results/kernels_games_detours_20261020f/forcing_check.py

The checker uses rational coefficients in the exact basis
`1,sqrt(3),sqrt(11),sqrt(33)`. Every one of the 171 point pairs is checked;
there is no floating tolerance and no inherited geometric filter. It
independently enumerates the output relation, finds and validates the
global four-colouring and twelve deletion witnesses, and examines all
35 core triples. It also checks 100,585 minimum-size diamond list
assignments over palettes of sizes 3 through 6 and all 4,096 pairs of
available six-colour lists for the two-hub identity. These finite list
checks supplement the proofs, rather than replace their quantifiers.

`forcing_check.json` stores exact coordinates, the complete edge list,
colouring witnesses, and the circumcentre certificate. `forcing_check.out`
records identical normal and optimized-Python stdout. All validity guards
use explicit exceptions and remain active under `-O`.

The separate parent audit `independent_audit.py` reconstructs every point
and all 171 distances with SymPy 1.14.0, without importing the custom field
arithmetic above. It agrees on the complete contact graph, interior degree
vector, four-colouring, and independently enumerated output multiplicities.
