# Square-sum graph: paths, separators, centers and distances

Let G_n be the simple graph on {1,...,n}, with distinct vertices adjacent exactly when their sum is a square. No loops are allowed. The user's three named sets have union S={1,...,12,16,17,18}; distinguish the induced graph G18[S] from G18 and from G12.

## The three displayed sets

G18[S] has exactly three connected components:

* A={10,6,3,1,8,17}, the path 10–6–3–1–8–17 (6 vertices, 5 edges).
* B={12,4,5,11}, the path 12–4–5–11 (4 vertices, 3 edges). Its central edge is 4–5.
* C={2,7,9,18,16}, the tree with edges 2–7, 7–9, 7–18, 9–16 (5 vertices, 4 edges). It is not a path: vertex 7 has degree 3.

Indeed deleting 13,14,15 from G18 gives precisely these components. In G18, 13 links A to B, 14 links B to C, and 15 adds edges to 1 and 10. The unique cycle is 1–3–6–10–15–1; all other 13 edges are bridges. The graph is connected and has three leaves 16,17,18, so it has no Hamilton path.

## The center graph and the 3–4–5 relationship

G12 has components P5=10–6–3–1–8, P4=12–4–5–11, P3=2–7–9. Their respective graph centers (minimum eccentricity within each component) are {3}, {4,5}, {7}. The graph induced on their union {3,4,5,7} has only edge 4–5. It is isomorphic to G4, whose only edge is 1–3: for example 1↦4, 3↦5, 2↦3, 4↦7. There are exactly four such isomorphisms, arising from independent endpoint and isolate swaps; none preserves numerical order. In this example the center-retaining operation is stable after one application because all its components are isolated vertices or edges.

The relationship is not uniform in n. G10 has the same center graph as G12. G11 has center set {3,5,7}, inducing three isolated vertices, so it is not isomorphic to G3. G13 has centers {7,12,13}, inducing one edge, isomorphic to G3. G14 has centers {4,5}, inducing K2, not isomorphic to G2 (two isolated vertices). These are finite exact checks, not a general classification.

The finite range n=14,...,22 exhibits a further change. The graph centers are {4,5} at n=14,15,17,18,19,21; {5} at n=16; {4} at n=20; and {14,22} at n=22, where the radius is 5. The edge 4–5 stays a bridge through G21 and first ceases to be a bridge in G22: the cycle 4–5–11–14–22–3–13–12–4 witnesses this. These are exact finite observations and do not assert an extension past n=22.

G18 itself has graph-center vertices {4,5}, with eccentricity 7. Indeed its 16–17 distance is 13, so its radius is at least 7; the only vertices within 7 of both are 4 and 5, and both have eccentricity 7.

## Hamilton paths through order 18

Among positive orders n≤18, precisely n=1,15,16,17 have Hamilton paths, counting the one-vertex path. For n=15,16,17 the path is unique up to reversal. Explicitly, the G15 path is
8–1–15–10–6–3–13–12–4–5–11–14–2–7–9.
At n=16 append 16, and at n=17 also prepend 17. The G18 graph has three leaves and no Hamilton path.

Proof: for 2≤n≤12 the vertices 1 and 2 lie in distinct components of G12 and hence of G_n; at n=13 the C component is still separate. G14 is a tree with leaves 8,9,10, so no Hamilton path. At n=15 the new vertex joins 1 to 10, creating the 5-cycle 1–3–6–10–15–1, and removes leaf 10. The remaining leaves 8 and 9 must be path endpoints. The cycle lies between its attachments at 1 and 3; to include all cycle vertices, the path must traverse the four-edge arc 1–15–10–6–3, not the edge 1–3. Everything else is forced. Vertex 16 extends the 9 endpoint, vertex 17 extends the 8 endpoint, and vertex 18 branches at 7, giving three leaves. Consequently in each of G15, G16 and G17, edge 4–5 is forced in every Hamilton path because it is a bridge; edge 1–3 is present but appears in no Hamilton path because the full cycle must be traversed via the other arc; and 2–4 is absent from the graph. This trichotomy is non-vacuous because those graphs do have Hamilton paths. G18 has no Hamilton path at all. Since G14 is a tree, the G15 5-cycle is also the first odd cycle (indeed, the first cycle) in this nested family.

## Three-way decompositions

Exactly two 3-vertex deletion sets yield three equal 5-vertex components:

* Delete {3,6,14}: {1,8,10,15,17}, {2,7,9,16,18}, {4,5,11,12,13}.
* Delete {3,14,17}: {1,6,8,10,15}, {2,7,9,16,18}, {4,5,11,12,13}.

This is an exhaustive finite classification, not just two examples; see check.py and the independent audit. Three deletions are minimal because the number of surviving vertices must be divisible by 3, and the k=0 graph is connected.

The degree-sensitive component identity certifying this enumeration is
c(G18−D) = sum(deg(v) for v in D) − |E(D)| − |D| + β,
where β is 1 if D avoids the unique cycle, else 0. It follows from c=|V|−|E|+β for the surviving subgraph. For triples, the candidate signatures (degree sum, induced-edge count, β) and their counts when c=3 are (5,0,1):75, (6,0,0):162, (6,1,1):43, (7,1,0):63, (7,2,1):2, (8,2,0):4. Checking their actual component sizes leaves exactly the two listed cuts.

For edge deletion, two cuts are the minimum for three components. Since G18 is unicyclic, a two-edge deletion produces three components iff both removed edges are bridges. Thus the 13 bridges give exactly 78 three-component decompositions. Equal 6/6/6 is impossible: either edge in such a cut would individually separate a 6-vertex component from the other 12, but the only bridge with a 6/12 split is 11–14. The optimally balanced profile 5/6/7 is attained by exactly two bridge pairs:

* {3–13, 2–14}: components {1,3,6,8,10,15,17}, {2,7,9,16,18}, {4,5,11,12,13,14}.
* {3–13, 11–14}: components {1,3,6,8,10,15,17}, {2,7,9,14,16,18}, {4,5,11,12,13}.

## Which edge is central?

Under standard unnormalized shortest-path edge betweenness, the bridge 4–5 has value 10×8=80 in G18. The bridge 4–12 has value 9×9=81 and is the unique global maximizer; 12–13 also has value 80. These statements differ from 4–5 being the central edge of component B, and from {4,5} being the vertex center of G18.

Proof of the global bound: for any edge uv of a connected n-vertex graph, partition vertices into those closer to u, closer to v, and equidistant. A shortest path traversing uv has its two endpoints in the first two classes, so edge betweenness is at most the product of the two class sizes, at most floor(n²/4). Equality requires a balanced partition with no ties and every cross-class shortest path using uv. If a different cross-class edge existed, its endpoints would have their own one-edge shortest path, contradicting equality; thus uv is a balanced bridge. Conversely a balanced bridge attains the bound. For G18, 4–12 is the unique bridge splitting 9 and 9.

## The missing 2–4 edge and infinite-graph distance

1–3 and 4–5 are square-sum edges; 2–4 is not. Stronger: each of the pairs {1,3}, {2,4}, and {4,5} has no common neighbor in the entire positive-integer square-sum graph. A common neighbor would make two positive squares differ by 2 or 1. A difference of squares cannot be 2 modulo 4; the only squares differing by 1 are 0 and 1, which cannot be the sums of positive vertices. In particular the two present edges lie in no triangle.

Nevertheless dist(2,4)=3 in the infinite graph, with witness 2–23–77–4. This is the first such 3-edge path as n grows: G14 already contains 2–14–11–5–4, and distance is 4 for every 14≤n≤76 and 3 for every n≥77. To prove the threshold, any three-edge path 2–a–b–4 has a=s²−2, b=t²−4, and s²+t²−6=u². If a,b≤76, then 3≤t≤8 and (u−s)(u+s)=t²−6 is one of 3,10,19,30,43,58. The three even values are 2 modulo 4; factorization of 3 forces s=1, while 19 and 43 force s=9 and s=21, giving a=79 or 439. All violate positivity or the bound. For t=9, 75=5·15 gives u=10,s=5,a=23,b=77.

## Reproduction and limits

Run `python check.py`, `python -O check.py`, `python center_check.py`, `python -O center_check.py`, `python hamilton_check.py`, and `python -O hamilton_check.py`; matching normal and optimized outputs are retained. The verifier uses explicit checks rather than optimization-sensitive assert statements and rejects malformed vertex/edge inputs. The center and separator conclusions concern the definitions above. The square-sum graph's first triangle threshold is 30 (6,19,30) and first 4-cycle threshold is 46 (1,3,46,35), independently exhausted in check.py; no claim of absence beyond those thresholds is made. A separate [independent verifier](independent_check.py) checks the Hamilton and center claims; the [core audit](../core_audit/AUDIT.md) checks additional graph claims.
