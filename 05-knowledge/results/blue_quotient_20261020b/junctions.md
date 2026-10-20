# Seven-port cubic junctions: the missing coordinates and a complete repair

Status: PROVED for the explicit null vectors, Kempe inequality and tree
pairings; FINITE-EXACT for the rank-84 lower bound. Their combination gives
an exact, reproducible seven-port reconstruction theorem. No general
Petersen-minor theorem is claimed.

## 1. Object and inheritance

The inherited `../boundary_support_20261020_incoming.md` repairs the
parity-impossible equal-colour triple by a genuine cubic vertex: its
three ports have three distinct colours. Together with two free wires,
this is a seven-pole J with 54 supported boundary words, each with one
extension. There are binom(7,3)*3=105 labelled choices of this closure.
The companion `snarks.md` now realizes 4J as a **connected** 17-vertex
seven-pole with an explicit rooted minor back to J.

These are distinct objects from the inherited 105 perfect matchings on
**eight** ports, partitioned into 7+42+56 by the F_2^3 direction atlas.
The eight-port measurement matrix has rank 91 on 274 parity orbits;
the seven-port matrix here has rank 84 on 91 parity orbits.
Numerical equality of the number of closures does not identify the carriers.

Let f_P(x) count proper three-edge-colourings extending the specified
seven-port word x. Every internal vertex is cubic; free wires are allowed.
Global colour permutations preserve f_P. The parity lemma requires all
three colour multiplicities to be odd. The 546 admissible words split into
91 free S_3 orbits:

- 21 orbits of multiplicity pattern (5,1,1), indexed by the two singleton
  port labels;
- 70 orbits of pattern (3,3,1), indexed by a singleton port and an
  unordered partition of the other six ports into two triples.

Use one coordinate per orbit and the ordinary inner product on these
91 coordinates. Actual closure colouring counts are six times these
inner products. This factor matters for the constants below. If gluing
two free wires creates a vertex-free closed wire circle, retain that
circle with three colourings, as prescribed by tensor contraction;
do not erase it with weight one. For poles without isolated free wires
the closures have the ordinary cubic-graph interpretation directly.

## 2. Seven explicit missing directions

For each port i define w_i on the 91 orbits by

    w_i(x) = +1, if i is a singleton and x has pattern (5,1,1);
             -1, if i is a singleton and x has pattern (3,3,1);
              0, otherwise.

Every one-junction/two-wire kernel J satisfies

    <J,w_i> = 0.                                                (1)

If i is on a wire it cannot be a singleton. Otherwise fix its unique
colour at the central vertex. Each wire uses one of the other two
colours; their four assignments have signs +,-,-,+. This proves (1).

The seven w_i are independent: on an orbit of pattern (3,3,1) with
singleton i, only w_i is nonzero. Thus the 105 by 91 junction matrix
A has rational rank at most 84.

The exact program `junctions.py` constructs all rows and columns and
finds rank 84 modulo 101. A nonzero minor modulo a prime is nonzero
over Q, so the upper bound is attained:

    rank_Q A=84,      ker_Q A=span{w_0,...,w_6}.                 (2)

This lower-bound step is FINITE-EXACT, with a reproduction command and
pivot columns in `junctions.json`; it does not use floating eigenvalues.
As a further control, a column of pattern (5,1,1) is hit by 15 closures,
and one of pattern (3,3,1) by 9.

An equivalent integer description is useful. Write A_ij for the
(5,1,1) coordinate with singleton ports i,j, and B_i;U,V for the
(3,3,1) coordinate. Every vector invisible to all 105 closures has
the form

    A_ij = c_i+c_j,       B_i;U,V = -c_i.                       (3)

Equation (3) describes the ambient linear blind space. It does **not**
assert that an arbitrary vector in it is the difference of two
realizable multipole extension tensors.

## 3. What the missing coordinates mean for an actual cubic pole

There is a sign constraint stronger than arbitrary ambient linear algebra:

    <f_P,w_i> <= 0.                                            (4)

More precisely, let N_i count full colourings of P in which the colour
at port i is unique on the boundary, and every one of the three paths
formed by the other two colours has differently coloured boundary ends.
Then

    6 <f_P,w_i> = -N_i.                                        (5)

Proof. Fix a colouring whose boundary colour at i is unique. The other
two colours form cycles and exactly three paths joining the remaining
six ports. Let r count paths with differently coloured ends. Both
boundary colour multiplicities are odd, so r is odd, hence 1 or 3.
Switching the colours along any of these paths preserves the proper
colouring and the unique colour at i. The three independent switches
partition these colourings into orbits of size eight.

If r=1, switching the two paths with equally coloured ends gives signs
+,-,-,+, cancelling their contribution to w_i. If r=3, all boundary
words have multiplicities (3,3,1) and contribute -1. Sum over all such
switch orbits and all three choices of the unique colour. Finally pass
from words to the free six-element global colour orbits. This gives (5).
Internal two-colour cycles do not affect the argument; switching only
the three boundary paths already gives the required partition.

Thus the missing coordinates count a concrete Kempe-path phenomenon.
The disconnected junction J and its connected realization 4J have
N_i=0 for every i. They occupy the equality boundary of these seven
inequalities. This also explains why a cyclic realization can preserve
J without becoming a five-vertex cubic tree.

## 4. Seven actual tree closures complete the measurements

Take the path on five internal vertices 0-1-2-3-4. Add ports so that
every vertex is cubic: their incident vertices, in label order, are

    (0,0,1,2,3,4,4).

For each i, relabel the sole port at the middle vertex 2 as port i;
fix any bijection of the other six labels. Call its extension kernel K_i.
It has 96 boundary words, with one extension each. Then

    <K_i,w_j> = -4 delta_ij.                                   (6)

Proof. If the unique boundary colour occurs at the middle vertex,
the internal edges of that colour are forced to be 0-1 and 3-4.
The other two colours form three boundary paths, with respectively
1,3,1 internal vertices. All three paths have differently coloured
ends. For each of three unique-colour choices there are 2^3 colourings,
so N_i=24 and (5) gives -4.

If the singleton port is at an adjacent vertex 1 or 3, the other four
vertices admit no internal matching of its colour, so there is no such
colouring. If it is at an endpoint vertex, the two-colour paths have
2,2,1 internal vertices; only one has differently coloured ends.
The switch cancellation in Section 3 applies. This proves (6) at
every port, independently of the remaining labels.

**Reconstruction theorem.** The 105 junction closures and these seven
tree closures determine the entire extension tensor of any cubic
seven-pole.

Indeed, if two invariant parity tensors give the same 105 junction
counts, their difference is sum_i c_i w_i by (2). Equality of the
seven tree counts gives -4c_i=0 by (6). Thus the difference is zero.
Equivalently the 112 by 91 augmented matrix has rank 91.

This determines exact extension multiplicities, and hence support.
The tests can therefore certify context-independent colour-equivalence.
They do **not** certify that one pole is a rooted minor of another, or
that a completed reduction remains simple and bridgeless. Those need
the disjoint paths and graph-class checks in `snarks.md`.

## 5. Map, loss, tests, and the Pythagorean interface

Source: a cubic seven-pole with labelled ports.
Target: its 91-coordinate extension vector, then 112 physical closure counts.
Map: enumerate extensions and pair with the displayed closure kernels.
Preserved predicate: equality of all boundary extension counts.
Destroyed information: internal graph structure and minor branch sets.
Needed sidecar: actual rooted minor paths plus the completed graph class.
Cheapest decisive tests: (1), (6), the rank-84 minor modulo 101, and
comparison of internal-edge enumeration with cut-XOR constraints.

For f,g on this finite boundary space, put M=||f||^2||g||^2 and
C=<f,g>. The inherited sum/difference construction on f tensor g
and g tensor f gives orthogonal vectors with squared lengths

    2(M+C^2),             2(M-C^2).                             (7)

For nonzero nonnegative f,g, the nondegenerate isosceles case is exactly
C=0, namely disjoint boundary supports. Complete measurements can
recover C; the 105 junction tests alone cannot recover an arbitrary
tensor because of (2). Neither (7) nor equal measured coordinates
supplies a graph-minor map.

The alternative flower strip sequence
(4^h-1,2^(h+1),4^h+1) remains a primitive Pythagorean family through
(3,4,5), but becomes flatter, not isosceles. No identification of
complete graphs, bipartite graphs and snarks with all primitive triples
has been established.

## Reproduction and exact boundary

From the repository root:

    python3 05-knowledge/results/blue_quotient_20261020b/junctions.py

The independent controls enumerate all internal tree-edge colours and
all 3^7 boundary words, comparing the counts to the cut-XOR criterion.
The receipt is `junctions.json`. A separate general edge-colouring
engine in `../../../04-computation/blue_quotient_20261020b_snarks.py`
also counts actual Kempe path components and independently checks (5):
the connected 17-vertex pole has all N_i=0; the path tree has N_i=24
only at its middle port; the Petersen-minus-{0,1,3} seven-pole has
N_i=24 at two ports and zero at the other five.
The exploratory `junction_probe.py`
also reports orders 3,5,9, but its order-nine modular rank is not used
as an all-order theorem or as an exact rational rank claim.
