# Graph interfaces, right triangles, and a Petersen/K5 support identity

**Status: PROVED (elementary statements below) + FINITE-EXACT (the named
graph census).** This is a research note, not a reserved canon theorem. No
claim of literature priority is made. There is no proof of Kuratowski's
theorem or a general Petersen-minor/flow theorem here.

## Inheritance and live concepts

The closest proved mechanism is
[THM-4116, boundary-state gluing](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md):
the number of three-edge-colorings is the dot product of ordered boundary
extension vectors. Its Petersen adjacent-vertex interface is the hostile;
its K4 interface is the positive control. The corrected near miss is
[MISTAKE-507](../../01-canon/MISTAKES.md): three-edge-coloring of a cubic graph
means an F2^2 flow, not an F3 flow, and numerical graph/arithmetic labels do
not supply a transport map. The underused sidecar is the actual boundary
gluing permutation, before forgetting color and edge addresses.

[THM-261, Petersen root orthogonality](../../01-canon/theorems/THM-261-petersen-root-orthogonality.md)
provides the exact pair-indexed carrier E(K5). It identifies Petersen with
the disjointness graph on those ten pairs, not with a Pythagorean triple.
[THM-3339, Fibonacci/Berggren current](../../01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md)
shows why three matching channels still require their signs and operations.
The prior [Kuratowski/Tait note](procgen_kuratowski_20260925_tait_kempe_triples.md)
is useful provenance, but its historical source claims do not replace the
elementary proofs here or a current primary-source audit.

The concept board was: (1) Pythagorean polarization, (2) flow group, (3)
boundary support, (4) gluing permutation, (5) maximum-cut defect, (6) oddness.
The anchor is a rigorous isosceles interpretation, the niche is the C5
five-pole, and the wildcard is whether cut frustration has an arithmetic
limit. The five-pole produced structure; normalized cut frustration fails
as a universal snark-angle coordinate.

## 1. An exact right triangle from a coloring interface

Let G be a finite simple cubic graph. Fix a nontrivial vertex cut X and one
common ordering of its boundary edges. For a boundary word sigma in
{0,1,2}^delta(X), let f(sigma) and g(sigma) count proper three-edge-colorings
on the two sides, including the prescribed boundary colors. Constraints
are imposed only at vertices of the relevant side. Write

```text
C = #Col_3(G),   U = <f,f>,   V = <g,g>.
```

THM-4116's restriction-and-gluing bijection gives C=<f,g>. Set

```text
L_- = ||f-g||,   L_+ = ||f+g||,   H = sqrt(2(U+V)).
```

**Proposition 1 (PROVED).**

```text
L_-^2 + L_+^2 = H^2,
L_+^2 - L_-^2 = 4 C.
```

Consequently, when both side vectors are nonzero, this is an isosceles
right triangle exactly when G is not three-edge-colorable. It may be
degenerate in the colorable case f=g. A snark lands exactly on the real
isosceles boundary, not merely near it.

*Proof.* Expand the two squared norms and use C=<f,g>. All coordinates are
nonnegative, so C=0 is precisely disjoint support. The nonzero-side condition
ensures that the boundary certificate describes two locally feasible pieces
which cannot be glued, rather than an already impossible side. QED.

This is a real Euclidean triangle. Its legs generally contain square roots.
No nonzero integer Pythagorean triangle is exactly isosceles: a=b>0 would
force c/a=sqrt(2). Therefore an integer-triple encoding must represent this
boundary through an infinite sequence or a limit, and must state the rule
choosing that sequence.

**Guaranteed nonvacuity (PROVED).** Every nonempty simple cubic graph has a
cut with f,g both nonzero. Take any cut locally maximal under single-vertex
switches. Every vertex has at most one neighbor on its own side; otherwise
switching it increases the cut. Each induced side is therefore a matching
plus isolated vertices. Color every internal matching edge, then assign
the remaining colors freely and distinctly to each vertex's boundary
semiedges. This gives at least one side coloring on both sides.

The map in Proposition 1 has source a graph with a marked cut and gluing,
target a real right triangle, and preserved predicate C=0. It loses the
support addresses, extension counts, graph reconstruction, and minor
operations. The full vectors and gluing permutation are the restoration
sidecar. It is a correct diagnostic, not a cheaper coloring algorithm or a
proof that a particular minor must exist.

## 2. The Petersen boundary really lives on E(K5)

Give each vertex of a 5-cycle one dangling boundary edge. Identify the five
boundary positions with {0,1,2,3,4}.

**Proposition 2 (PROVED).** The boundary extension vector of this five-pole
is a 0/1 vector with exactly 30 nonzero coordinates. Its support consists
of the words with color multiplicities (3,1,1) whose two singleton-colored
positions are adjacent in the 5-cycle. Every such word has exactly one
extension.

*Proof.* For each color c, counting incidences gives 5=2m_c+k_c, where m_c
is its number of internal edges and k_c its number of boundary appearances.
Thus all three k_c are positive odd integers summing to five, hence are
3,1,1. Let a be the color appearing three times. Exactly one internal edge
has color a. Its endpoints cannot have boundary color a, so they must be
the two singleton positions; in particular those positions are adjacent.
Conversely, color their joining edge a. Along the remaining four-edge path
the other two colors alternate uniquely, giving the desired extension.
There are 5 choices for the singleton pair, 3 for a, and 2 for the ordering
of the other colors: 30. QED.

Now form a cubic graph by taking two 5-cycles C and D and joining their
vertices by a perfect matching. Use the matching to identify both vertex
sets with one five-element set. The two cycles are then edge subsets of
the same K5.

**Corollary 2.1 (PROVED).**

```text
#Col_3(G(C,D)) = 6 |E(C) intersect E(D)|.
```

Both pole norms squared are 30. Their support intersection has six words
for each common singleton pair; restriction and gluing prove the formula.

For Petersen, the outer cycle C and inner pentagram D are complementary
5-cycles in K5. Their supports are disjoint because

```text
E(K5) = E(C) disjoint_union E(D).
```

This gives a complete elementary proof that Petersen is not
three-edge-colorable. The vectors f,g themselves are orthogonal and have
equal norm sqrt(30), so the triangle (||f||,||g||,||f-g||) has sides
(sqrt(30),sqrt(30),sqrt(60)). Proposition 1 instead gives
(sqrt(60),sqrt(60),sqrt(120)). Both are exactly isosceles right triangles.

For the pentagonal prism D=C. Here f=g, C=30, and Proposition 1 gives
(0,sqrt(120),sqrt(120)). The local pieces and their individual norms are
identical in both examples. The obstruction is the matching identification
of their addresses. Forgetting that permutation destroys the answer.

There are exactly 12 undirected Hamilton cycles on five fixed labels. With
C fixed, exact enumeration gives the overlap histogram

| Common edge count | Number of labelled D | #Col_3(G(C,D)) |
|---:|---:|---:|
| 0 | 1 | 0 |
| 2 | 5 | 12 |
| 3 | 5 | 18 |
| 5 | 1 | 30 |

The complement is the unique uncolorable matching gluing within this
specific labelled two-pentagon family. This is not a classification of
all snarks. The actual bridge to K5 is the pair-support carrier; it does
not identify K5's own edge-coloring problem with the Petersen problem.

### The three graphs also have literal contraction/subdivision maps

In the same Petersen presentation, contract the five matching spokes
{i,i+5}. The quotient edges are exactly the outer cycle and its pentagram
complement, hence the quotient is K5. The branch sets are the five displayed
pairs. The retained two-color decomposition of E(K5) and the way each
degree-four vertex splits into two degree-three vertices reconstruct the
original Petersen graph; the undecorated K5 quotient forgets that pairing.

There is also an explicit K3,3 subdivision. With the script's standard
Petersen labels, delete vertex 0 and use shores {2,8,9} and {3,6,7}. Its nine
internally disjoint paths are

```text
2-3,     2-1-6,   2-7,
8-3,     8-6,     8-5-7,
9-4-3,   9-6,     9-7.
```

Thus Petersen contains both Kuratowski minors by explicit certificates.
Conversely, one obtains Petersen from this K3,3 by subdividing the three
edges of a perfect matching and joining their subdivision vertices to one
new vertex. These graph operations give a precise relation among the
user's three graphs. They supply no correspondence with Euclid parameters
by themselves, and coloring of a minor does not lift without checking the
discarded local constraints.

## 3. Small exact examples and the graph-only normalization

The session's companion inheritance note derives the cut-independent
quantity

```text
s(G) = |E(G)| - MaxCut(G),
q(G) = C(G)^2 * 3^s(G) / 6^|V(G)|.
```

For any maximum cut, each side is a matching plus isolated vertices. An
isolated vertex has six extensions; an internal matching edge has twelve.
Each extension has a unique boundary word. Consequently UV=6^|V|/3^s, and
q=C^2/(UV) is the squared cosine between the two side vectors. This also
proves independence of the maximum cut chosen. The tensor pair f tensor g,
g tensor f has equal norm squared UV and inner product C^2, giving a
graph-invariant right triangle with normalized squared legs (1-q)/2 and
(1+q)/2. Its exact isosceles endpoint is q=0.

| Graph | n | C(G) | s(G) | q(G) | Oddness |
|---|---:|---:|---:|---:|---:|
| K4 | 4 | 6 | 2 | 1/4 | 0 |
| K3,3 | 6 | 12 | 0 | 1/324 | 0 |
| triangular prism | 6 | 6 | 2 | 1/144 | 0 |
| cube | 8 | 24 | 0 | 1/2916 | 0 |
| pentagonal prism | 10 | 30 | 2 | 25/186624 | 0 |
| Petersen | 10 | 0 | 3 | 0 | 2 |
| Heawood | 14 | 48 | 0 | 1/34012224 | 0 |
| dodecahedron | 20 | 60 | 6 | 25/34828517376 | 0 |
| flower J5 | 20 | 0 | 3 | 0 | 2 |

Oddness here is the minimum number of odd cycles in a spanning 2-factor.
The table is FINITE-EXACT. The script enumerates every perfect matching,
then counts the even components of its complementary 2-factor. If there
are r even components and no odd ones, that perfect matching contributes
2^r colorings with its edges assigned the distinguished first color.
An independent recursive edge-coloring count agrees through ten vertices;
the full count also agrees with independently generated boundary-vector
dot products for every named graph. A separate inherited perfect-matching
engine checks the four boundary-orbit multiplicities.

The formula distinguishes exact zero, but small positive q does not
identify a near-snark. The cube and Heawood graph are already close to the
isosceles endpoint in this normalization despite being bipartite. The
companion note's connected even-prism family gives an explicit colorable
sequence q -> 0. Normalized overlap carries a graph-size dilution effect.

## 4. Cut frustration goes toward the bipartite boundary in flower snarks

Use the standard flower graph on vertices a_i,b_i,c_i,d_i, 0<=i<k. Each a_i
is joined to b_i,c_i,d_i. The b vertices form a k-cycle, and the c and d
vertices form the 2k-cycle

```text
c_0,...,c_(k-1),d_0,...,d_(k-1),c_0.
```

**Proposition 3 (PROVED).** For every odd k>=3,

```text
MaxCut(J_k)=6k-3,
s(J_k)/|E(J_k)| = 1/(2k).
```

*Proof of the upper construction.* Alternate the colors of the centers
with i and give every leaf the opposite color to its center. All star
edges and all nonwrapping horizontal edges cross the cut. Precisely the
three wrapping edges fail to cross, so MaxCut>=6k-3.

*Proof of the lower obstruction.* Consider an arbitrary two-coloring. If
some module i has all three leaves the same color, follow the three
edge-disjoint horizontal walks of length k from this module around the
flower: the b-cycle returns b_i to itself, and the other two paths run
between c_i and d_i. All three have odd length and equally colored ends, so each
contains an uncut edge. If no module has three equally colored leaves,
every star contributes at least one uncut edge, for at least k>=3 uncut
edges. Thus every cut leaves at least three edges uncut. QED.

The exact transfer program checks k=5,7,9,11,21 and independently agrees
with exhaustive cut enumeration for J5. Therefore the usual odd flower
snarks approach bipartite graphs in the normalized maximum-cut sense.
They do not approach a universal positive interior angle in that metric.
This is a failure of a proposed scalar coordinate, not of the exact
support-compatibility interpretation in Proposition 1.

A smaller hostile is Petersen versus the dodecahedron: both have maximum
cut fraction 4/5, but the former has C=0 while the latter has C=60 and is
planar. The lost coordinate is extension compatibility.

## 5. What the elementary invariants do and do not separate

All cubic graphs have |E|/|V|=3/2. Feeding that ratio to Euclid's formula
produces (5,12,13) for every cubic graph, including K3,3, Petersen, the cube,
and all snarks. An order/size ratio cannot assign those families distinct
branches.

The corrected flow dictionary is structural. A nowhere-zero F3 flow on a
cubic graph can be oriented with every edge value 1. At a cubic vertex,
conservation then forces all three edges to point in or all three to point
out; the graph is bipartite. Conversely a bipartition gives this flow.
A nowhere-zero F2^2 flow is equivalent to proper three-edge-coloring:
three nonzero elements sum to zero exactly when they are all distinct.
The smallest separator is K4: three-edge-colorable but not bipartite.
K5 is Eulerian and has a nowhere-zero 2-flow. These are different flow
strata, not a uniform (3,4,5) labeling of the user's three graphs.

The exact adjacency characteristic polynomials are

```text
K5:       (x-4)(x+1)^4,
K3,3:     x^4(x-3)(x+3),
Petersen: (x-3)(x-1)^5(x+2)^4.
```

These are checks, not a proof of a graph-family/Pythagorean-tree
isomorphism. The companion arithmetic lane examines a different Seidel
spectral assignment. Any claim transferring Kuratowski or a Tutte theorem
still owes an operation-compatible graph map and the relevant preserved
predicate; Proposition 1 supplies only a coloring predicate.

## Reproduction and scope

```powershell
python 04-computation/graph_pythagorean_20261018_graph_probe.py
```

Script: [graph probe](../../04-computation/graph_pythagorean_20261018_graph_probe.py).
Saved output: [graph probe output](../../04-computation/graph_pythagorean_20261018_graph_probe.out).

The universe is the nine named cubic graphs in the table, all twelve
labelled C5 matching gluings, the three unlabelled simple cubic graphs
through seven vertices supplied by NetworkX's complete small-graph atlas,
and the five specified flower transfer lengths. There are no hidden
filters. Positive controls include K4, K3,3, the cube and prisms; hostile
controls include Petersen, J5, the Petersen/dodecahedron cut collision,
and the prism/Petersen equal-pole-norm collision. Maximum cuts through
twenty vertices are exhausted modulo global cut complementation. All
reported decisions use integers or rational numbers. Characteristic
polynomials are computed symbolically. No floating-point threshold decides
colorability, planarity, cut optimality, oddness, or a Pythagorean identity.
