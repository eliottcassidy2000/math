# Independent square-sum graph audit

Use the simple graph G_n with vertices 1,...,n and an edge between distinct x,y when x+y is a perfect square. Results below are scoped to this convention; the three sets in the prompt have 6, 4, and 5 elements, not 3, 4, and 5. The interpretation as 3/4/5 *component sizes* is supported by G_12, but the prompt does not specify n and does not assert a graph triangle.

## Structural decomposition

G_12 has exactly three components, all paths:
- 10–6–3–1–8 (5 vertices),
- 12–4–5–11 (4 vertices),
- 2–7–9 (3 vertices).

The union of the user's three sets is S={1,...,12,16,17,18}. The induced graph on S has no cross-edges and exactly the original three groups as components. Vertices 17, 16, 18 respectively add 8–17, 9–16, 7–18. Thus the first component remains a six-vertex path, the second stays a four-vertex path, and the third is a five-vertex fork with degree-3 vertex 7, not a path.

The decomposition is real and exact but is not the unique 4/5/6 component split of G_18 after deleting three vertices: exhaustive enumeration gives 20 such deletion sets. Deleting {13,14,15} is natural because it preserves the entire G_12 and attaches the above three new vertices.

## Three highlighted pairs

Pairs 1–3 and 4–5 are edges since their sums are 4 and 9. Pair 2–4 is a nonedge since its sum is 6. Yet all three pairs have a property independent of n: they have no common neighbor. If p,q share neighbor z, then q−p=(z+q)−(z+p) is the difference of two positive squares. Differences 1 and 2 cannot be such differences (difference 2 is impossible modulo 4). In particular, neither edge 1–3 nor 4–5 lies in a triangle, even though triangles occur elsewhere for large n. More generally any two vertices whose difference is 2 modulo 4 have no common neighbor.

The gap between 2 and 4 has a sharp distance profile: they are disconnected in G_n up to n=13, have distance 4 for 14≤n≤76, and distance 3 for n≥77. For the upper bounds use 2–14–11–5–4 and 2–23–77–4. Disconnection in G_13 and no shorter G_76 path are computationally verified; monotonicity supplies the stated ranges, while the common-neighbor obstruction excludes distance 2 for every n.

## 4–5 centrality depends on graph and metric

In G_18 this edge is a bridge. Its deletion produces sets of sizes 10 and 8, hence it carries 80 unordered pairs' unique cross-cut routes. Edge 4–12 has even greater ordinary unnormalized shortest-path edge betweenness, 81 (a 9/9 bridge), and is the unique global maximum. Thus 4–5 is very central in this explicit sense, but not most central.

It is a bridge in G_15 through G_21, but ceases to be a bridge starting at G_22. In G_21−{4–5} the two vertex sets are {1,3,4,6,8,10,12,13,15,17,19,21} and {2,5,7,9,11,14,16,18,20}. Vertex 22 connects to 3 and 14 across this cut, yielding the alternate path 4–12–13–3–22–14–11–5. The same alternate path remains in all G_n with n≥22.

## Triangles

For x<y<z positive, the three edge sums are increasing squares A=x+y < B=x+z < C=y+z. Conversely three increasing positive squares A<B<C give a triangle exactly when A+B>C and A+B+C is even, with x=(A+B−C)/2, y=(A+C−B)/2, z=(B+C−A)/2. For G_30, C≤59 and possible squares are 4,9,16,25,36,49. Exactly one triple satisfies the constraints: (25,36,49), giving (x,y,z)=(6,19,30). Consequently G_29 has no triangles and G_30 has exactly one. An explicit infinite family follows by multiplying all its vertex labels by k² for k≥1.

## Additional audited 3-way splits of G_18

Exactly two 3-vertex deletions yield three five-vertex components: {3,6,14} and {3,14,17}. Two edge deletions cannot yield 6/6/6. The two 5/6/7 edge splits are exactly {2–14,3–13} and {3–13,11–14}. These calculations treat 3-way decomposition as three connected components after explicit deletion, not as a triangle or an unspecified covering.

## Reproduction

From a Python 3 environment run:

    python check.py
    python -O check.py

The script uses only the standard library and explicit exceptions, not assertions; both modes pass. It independently builds each graph from square sums and exhaustively checks relevant triples and cuts, computes shortest paths and Brandes edge betweenness, and validates the claims above.
