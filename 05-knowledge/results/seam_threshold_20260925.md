# The 15/16 seam and a genuine triangular-ten halving construction

**Date:** 2026-09-25. **Status:** PROVED elementary graph and decoder
statements; FINITE-EXACT enumerations with explicit universes. This note
does not establish an equivalence of square-sum Hamiltonicity, planarity,
or Collatz convergence. No novelty or canon promotion is claimed.

## Inheritance and board

Closest mechanism: the exact rotation and K_(3,3) subdivision in
[decoder minors](decoder_minors_20260925.md), proving that the literal
`+2`/doubling graph G_N is planar iff N<=15. The endpoint-aware treatment
in [prime/square decoding](decoder_prime_square_20260925.md) prevents
unconditional degree-two forcing. Canonical hostile: Q_15 is Hamiltonian
while G_15 is still planar. Corrected near miss: confusing a triangular
number of edges with the same number of vertices. Least-used sidecar:
the support pair indexing each edge-object, recovered from
[THM-261, Petersen-root-orthogonality](../../01-canon/theorems/THM-261-petersen-root-orthogonality.md).
That canon already identifies disjoint pairs with the Petersen graph;
we retain the support-index interpretation and do not identify Lie algebras.

Anchor: explain the actual two thresholds. Niche: a support-preserving
map from triangular ten to regular five. Wildcard: whether three cyclic
triples can be actual edge-disjoint triangles of K5.

| Concept | Exact operation | Cheap control / lost information |
|---|---|---|
| Q_N square-sum graph | attach the new label and its square partners | Q14 has three forced endpoints |
| G_N +2/doubling graph | attach a leaf or a two-edge ear | G16 ear completes K_(3,3) |
| Triangular ten | ten support pairs from five vertices | no intrinsic orientation of 45 pairs follows |
| L(K5) matching | pair incident edge-objects and retain their center | 144 choices, rather than a unique half |
| Regular-five decoder | two outgoing arcs per owner | preserve original support/owner to recover direction |

## 1. The nearby thresholds are different, but the local operation is shared

Use Q_N for the simple graph on 1,...,N with adjacency iff the sum is a
square; use G_N for adjacency `{x,x+2}` or `{x,2x}`. Then:

| N | Hamiltonian path in Q_N | Q_N planar | G_N planar |
|---|---|---|---|
| 14 | no | yes | yes |
| 15 | yes, unique up to reversal | yes | yes |
| 16 | yes, unique up to reversal | yes | no |

The first nontrivial positive square-sum order is **15**, not16. The
one-vertex path is vacuous; Q_N has no Hamiltonian path for N=2,...,14.
The latter bounded statement is independently checked by unpruned DFS,
alongside the explicit local proof below. Thus "the first15 are
impossible" is not the literal standard indexing.

### Q14 to Q15: an ear repairs a three-leaf obstruction

Q14 is exactly the tree with three arms meeting at3:

```
3--1--8
3--6--10
3--13--12--4--5--11--14--2--7--9.
```

It cannot have a Hamiltonian path because its three leaves8,9,10 would
all have to be endpoints. Adding15 supplies the two-edge ear
`1--15--10`, with sums16 and25. It absorbs leaf10 and creates the
five-cycle `1--3--6--10--15--1`. The remaining leaves8 and9 force the
endpoints; vertices15,10,6 are then certified interior and force the
long side of that cycle. The only possible Hamiltonian path is

```
8,1,15,10,6,3,13,12,4,5,11,14,2,7,9.                    (T1)
```

It uses every edge except1--3. Adding16 supplies only9--16, which extends
the endpoint and preserves uniqueness. All square partners were checked
directly; no argument forces both incident edges at an unclassified
degree-two endpoint.

### G15 to G16: an ear completes a nonplanarity obstruction

Adding15 to G14 supplies just13--15, a leaf. Adding16 to G15 supplies
the ear `8--16--14`: a doubling edge and an additive-step-two edge.
That ear is precisely the missing path in the already recorded
K_(3,3) subdivision, with shores `{5,8,12}` and `{6,10,14}`. The other
eight paths lie in G15. Thus this attachment destroys planarity.

There is a genuine shared operation: introduce a degree-two vertex
between two specified older vertices. Suppressing that vertex recovers
a marked added edge. In Q15 the marked edge is1--10, whose sum11 is
not square; in G16 it is8--14, which is neither a doubling nor a
step-two edge. Retaining the two edge types/labels is necessary to
reverse either suppression in its arithmetic category.

The source/target/predicate distinction matters: an ear in a tree
creates a unicyclic planar graph and here repairs Hamiltonicity; an ear
between the particular G15 attachments completes a forbidden planar
configuration. This is a structural analogy with an explicit common
operation, not a preserved-predicate map between the two graph families.
Their cyclomatic ranks change from0 to1 and from5 to6 respectively.

## 2. What T4=10 counts, and why order10 arose earlier

With `T_m=m(m+1)/2`, the identity `T4=10` counts the edges of K5, or
the independent arc-orientation slots of a five-vertex tournament.
A four-vertex tournament has T3=6 edges; a ten-vertex tournament has
T9=45 edges. These are different objects.

The earlier ten-vertex obstruction came from `3A+1` at the smallest
nontrivial regular odd core, A=3. Thus `Q[C3,C3,C3,1]` has10 vertices.
Regular odd order1 is trivial; order3 first supplies an internal cycle
and no pair module. This explains that particular ten without assuming
that the ten vertices are edges of a five-vertex graph. Already A=5
would give16 vertices, not a triangular number.

There is nevertheless a useful exact **edge-object** construction.
Take the ten unordered support pairs of `{1,2,3,4,5}` as vertices.
Join two when they share one endpoint. This is L(K5), with30 edges and
degree6. Its complement joins disjoint supports: the Petersen graph,
with15 edges and degree3, as in THM-261. The split

```
45 = 30 incident-support pairs + 15 disjoint-support pairs
```

is an exact organization of the45 pair slots on ten objects. It supplies
an incidence relation, not an orientation of every pair. No cosmetic
tournament is imposed on missing or tied comparisons.

## 3. A marked 10-to-5 decoder from regular tournaments

Let H be a regular tournament on the five original vertices. Each
vertex v has exactly two outgoing arcs. Pair their two underlying
support edges and call the pair P_v. Then:

1. The five P_v partition the ten support edges, because each edge has
   exactly one tail.
2. Each P_v is a connected two-vertex set in L(K5), centered at v.
3. Every two distinct sets P_u,P_v have an edge between them: if u->v,
   the edge-object `{u,v}` lies in P_u and meets each edge-object in P_v
   at v.

Thus the five pairs are an explicit branch-set model contracting L(K5)
to K5. More strongly, the original orientation decodes exactly from
the support and center: the support edge `{u,v}` belongs to P_u iff
the original arc was u->v. This is a genuine reversible representation
of a regular five-tournament by a marked pairing of ten edge-objects.

It is ordinary connected-pair contraction, not uniform module halving:
L(K5) has no two-vertex module. For example, incident supports `{1,2}`
and `{1,3}` are distinguished by `{2,4}`; relabeling proves the general
incident case. Disjoint supports are distinguished by a support joining
one endpoint to the fifth vertex. Therefore the pair labels and centers
must be retained rather than inferred from twin neighborhoods.

### Exact choice count: 144 pairings, 24 regular decoders

Every perfect matching in L(K5) pairs incident support edges. Orient
each support edge away from its pair's common endpoint. All outgoing
degrees are even, since every pair contributes two outgoing arcs.
Conversely, any orientation of K5 with even outgoing degrees, together
with a pairing of the outgoing arcs at each vertex, gives such a
perfect matching. These operations are inverse.

The possible score multisets are `(2,2,2,2,2)` and `(0,2,2,2,4)`.
Indeed every score is0,2,4, their sum is10, and at most one vertex can
have score0 or score4. The second family is necessarily a source over
a cyclic triple over a sink. It has `5*4*2=40` labeled orientations.
The source's four outgoing arcs have three pairings, while the other
nonempty outgoing stars have unique pairings.

There are24 regular labeled five-tournaments. An elementary count fixes
vertex1, chooses its two outneighbors in six ways, and orients the edge
inside its outneighbor pair and the edge inside its inneighbor pair in
two ways each. Regularity then uniquely forces the four cross edges.
Every outgoing star has size two and hence has a unique pairing.
Consequently

```
# perfect matchings of L(K5) = 24 + 3*40 = 144.           (T2)
```

Exactly the24 matchings with five distinct centers are the regular
10-to-5 decoders. The other120 use four distinct centers, one twice.
This establishes the construction and its lack of automatic uniqueness.
It does not transfer the integer's numerical Collatz guard or supply
an inverse-Collatz route from root4.

### A cheap hostile to a stronger triangular identification

The ten edges of K5 cannot partition into three edge-disjoint triangles
and one leftover edge. The union of triangles has even degree at every
vertex, whereas deleting one edge from K5 leaves its two endpoints with
odd degree3. Therefore the three cyclic blocks in the ten-vertex
substitution cannot literally be three edge-disjoint triangular circuits
of K5 under an incidence-preserving identification. Edge-objects and
cycles require separate sidecars.

## 4. Reproduction, controls, and next question

Run:

```
python 04-computation/experiments/seam_threshold_20260925.py
python -O 04-computation/experiments/seam_threshold_20260925.py
```

See [script](../../04-computation/experiments/seam_threshold_20260925.py)
and [output](seam_threshold_20260925.out). The core uses exact integers
and explicit exceptions, unaffected by optimization. Universes: every
Q_N for N1..16, local G14..G16, all1024 labeled orientations of K5,
every perfect matching of L(K5), and every triple of K5 triangles.
No random samples or undeclared filters. Independent optional NetworkX
planarity checks confirm all six local flags; the earlier rotation and
subdivision remain the planarity proof.

The positive controls are the unique Q15/Q16 paths and all24 regular
branch-set decoders. Hostiles are Q14's three leaves, the non-module
pairings, and the impossible three-triangle packing. These refine the
board: triangular ten carries a useful support/owner structure; the
ear operation explains both nearby thresholds but preserves different
target predicates. The next concrete question is whether a selected
edge-owner decoder can carry an arithmetic guard without storing an
already-known Collatz certificate. That additional transfer remains OPEN.
