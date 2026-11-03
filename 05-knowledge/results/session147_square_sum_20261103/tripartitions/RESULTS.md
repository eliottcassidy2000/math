# Square-sum tripartitions: exact graph and structural extensions

Working definition: the simple graph on a set V of positive integers has an edge {x,y}, x≠y, iff x+y is a perfect square. The specified set V is exactly {1,2,3,4,5,6,7,8,9,10,11,12,16,17,18}. It is not the prefix graph Q_15: 13, 14, 15 are absent and 16, 17, 18 are present.

## The specified 3-way decomposition

Its exact edges, grouped as the user specified, are:

* A={10,6,3,1,8,17}: 10–6–3–1–8–17. The successive sums are 16,9,4,9,25. Six vertices, five edges.
* B={12,4,5,11}: 12–4–5–11. Sums 16,9,16. Four vertices, three edges.
* C={2,7,9,18,16}: 2–7, 7–9, 7–18, 9–16. Sums 9,16,25,25. Five vertices, four edges. This is a tree, not a path: vertex 7 has degree 3.

No edge joins different groups: these are exactly the three connected components. This is obtained by direct enumeration of pair sums 4,9,16,25 (36 exceeds 17+18). Thus the groups have edge counts 5,3,4, a permutation of 3,4,5; their vertex counts are 6,4,5. The Pythagorean identity 3²+4²=5² is true but does not alone force the observed three components or their sizes.

The minimum number of vertex-disjoint paths covering V is exactly four. The first two components each require one path and the third requires two because its degree-3 vertex cannot lie in a single path covering it. Three minimum covers are obtained, and only these, by removing respectively 7–2, 7–9 or 7–18 from C, leaving A and B intact. Removing 9–16 leaves the degree-3 branch.

## A canonical arithmetic 3-way partition, not a proper 3-coloring

For any positive vertex set, partition V=E∪O1∪O3 into even, odd ≡1 (mod 4), and odd ≡3 (mod 4). O1 and O3 are independent; odd–odd edges cross between them. The even vertices split further into E0 (≡0 mod 4) and E2 (≡2 mod 4), with no edges between them. In the prefix graph Q_n:

* E0 induces precisely 4·Q_floor(n/4): 4u+4v is a square iff u+v is a square.
* E2 consists of 4u+2 with 0≤u≤floor((n−2)/4), and distinct such vertices are adjacent iff u+v+1 is a square.
* For odd r and even s modulo 8, a mixed-parity edge must satisfy (r,s) in {(1,0),(3,6),(5,4),(7,2)}. This is necessary, not sufficient: the actual sum must be square.

Reason: even squares are 0 or 4 mod 8, odd squares are 1 mod 8. The map u↦2u is NOT a copy of Q_floor(n/2): 1–3 is an edge, but 2–6 is not. In particular, one cannot infer that E itself is independent, so this 3-part partition is not a proper 3-coloring. Stronger: no coloring of all positive integers that depends only on x mod m can be proper. Choose k with 2k²>m, x=2k²−m, y=2k²+m; they are distinct positive numbers congruent mod m and x+y=(2k)².

The user examples fit: 4–5 and 1–3 are edges, while 2–4 is not (sum 6≡2 mod 4). More strongly, 2 and 4 have no common neighbor in the whole positive graph: their two putative square sums would differ by 2, impossible modulo 4. Their distance is exactly 3, witnessed by 2–23–77–4 with sums 25,100,81.

## Difference-of-squares endpoint test

For x<z put d=z−x. Common neighbors y correspond bijectively to factorizations d=uv, with 0<u<v, u≡v (mod 2), via a=(v−u)/2, b=(v+u)/2, y=a²−x, subject to y>0 and y distinct from x,z. For Q_n also require y≤n. Proof: x+y=a², z+y=b², so d=(b−a)(b+a); the converse is immediate. In particular d≡2 (mod 4) excludes a common neighbor. The converse without the positivity/distinctness filters is false.

## Complete triangle parameterization and the 3–4–5 boundary

Every triangle on positive distinct vertices is uniquely determined by distinct positive integers a<b<c satisfying:

1. a²+b²>c², and
2. a²+b²+c² is even (equivalently an even number of a,b,c are odd).

Its ordered vertices are

x=(a²+b²−c²)/2, y=(a²+c²−b²)/2, z=(b²+c²−a²)/2.

The displayed inequalities yield 0<x<y<z, and x+y=a², x+z=b², y+z=c². Conversely square sums on a triangle are distinct because the vertices are distinct, and solving their three linear equations gives these formulas. For a triangle contained in Q_n, additionally z≤n.

If n≤30, the largest square sum is at most 2n−1≤59, so c≤7. Direct inspection of root triples 1≤a<b<c≤7 under the two conditions leaves only (5,6,7), yielding (6,19,30). Hence Q_29 is triangle-free and Q_30 has precisely this triangle. This does not say smaller Q_n are bipartite: Q_15 already has a 5-cycle (prior THM-4505).

There is an exact link to the 3–4–5 right triangle. Allowing zero, roots 3,4,5 give vertices (0,9,16) and square-sum edges 0–9, 0–16, 9–16. Zero is NOT in the original graph; the last edge is in C. More generally, roots (2t+1,2t+2,2t+3) give vertices (2(t²−1), 2t²+4t+3, 2t²+8t+6). For integer t=1 this is (0,9,16); for every t≥2 all vertices are positive, starting with (6,19,30). The formula proves an infinite triangle family; it does not claim all triangles belong to it.

## General three-block path concatenation

Let P1,P2,P3 be three pairwise vertex-disjoint, fixed paths, each with at least two vertices. For each permutation π and choices ε_i∈{0,1} of orientations, concatenate the three oriented paths. This is a path precisely when both junction endpoint pairs are edges. Consequently the exact number of directed concatenations is the sum over the 3!·2³=48 orientation/permutation choices of the product of those two edge indicators. The number of undirected paths is half this count. This is an endpoint compatibility criterion, not a claim that all paths on their union respect those blocks. For the user's exact components there are no cross edges, so no such concatenation exists; moreover C is not itself a path.

## Prior work and audit

THM-4505 at pinned commit f139791 gives the unique Hamiltonian path in Q_15 and says Q_15 has an unused edge {1,3}; its domain is 1..15, not the present V. THM-4851 summarizes Q_n cycle results and THM-4850 concerns different paired Hamiltonian paths. The elementary observations above are proved here but novelty relative to the whole literature/repository is not asserted.

Run `python audit.py` and `python -O audit.py`. Explicit checks, not assert statements, compare both sides of the triangle formula to Q_200, validate the exact user graph and all its minimum path covers, check parity and scaling identities, the common-neighbor factorization for x<z≤60 with y≤200, and the boundary/family formulas. Both modes passed.
