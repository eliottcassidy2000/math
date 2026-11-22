# Planar unit-distance geometry of literal tournament seams

A **unit-height seam** is a pair of tournaments differing by reversal of just one arc with absolute Hamilton-path count difference two. This packet gives a faithful Euclidean soft amplifier inside the literal seam graph, and two ways in which the geometric analogy necessarily stops: hard color relations and geometric encoding of the Hamilton height. It does not improve a bound for the chromatic number of the plane.

## 1. An explicit seam soft amplifier with planar ports at distance two

Fix \(k\ge2\) consecutive disjoint ordered triples in an ordinal sum. In any triple \(a,b,c\), use three local states:
\[
 T_0=(a,b,c),\qquad T_2=(b,c,a),\qquad
 T_1=(a\to b\to c\to a).
\]
The first two are transitive orders; the middle state is a directed triangle. The edges \(T_0\leftrightarrow T_1\leftrightarrow T_2\) are literal single-arc reversals, and \(T_0,T_2\) differ at two arcs. Allow any tuple \(t\in\{0,1,2\}^k\) with at most one coordinate 1. The ordinal product rule for Hamilton counts gives height 1 at the all-even tuples and height 3 at the tuples with one 1. A flip anywhere else cannot join two distinct selected vertices: a comparison in two different coordinates needs at least two flips, and a change between 0 and 2 at one coordinate also needs two. Hence the induced labelled seam graph is exactly the 1-subdivision of the \(k\)-cube.

This graph has a faithful rational plane unit-distance representation. Choose rational unit vectors \(u_1,\dots,u_k\) and send \(t\) to \(\sum_i t_i u_i\). They can be chosen constructively so that among all coefficient vectors \(c_i\in\{-2,-1,0,1,2\}\), no nonzero sum vanishes, and the only unit sums use a single coefficient \(1\) or \(-1\). Given vectors through \(j-1\), for each nonzero old sum \(v\) and new coefficient \(c_j\ne0\), the equation \(v+c_ju_j=0\) excludes at most one direction and \(\|v+c_j u_j\|=1\) at most two (after squaring it is a nonconstant linear functional on the circle). At most \(12(5^{j-1}-1)\) directions are excluded. Test the injective rational parametrization \(((1-s^2)/(1+s^2),2s/(1+s^2))\) at integer parameters from 0 through \(12\cdot5^{j-1}\). The old sum zero is harmless, as by induction it had all old coefficients zero. This is an exact finite algorithm, not a floating-point drawing.

Take ports differing in coordinate \(a\) alone, from 0 to 2. Their Euclidean distance is exactly two. For each other coordinate \(b\), the route in the binary cube which changes \(b,a,b\) yields a length-six route in the subdivided cube. These routes are internally disjoint. Omitting the direct length-two route leaves their induced union: each selected midpoint has just its two selected binary endpoints as neighbors. Thus this is simultaneously an induced rational unit-distance graph with ports at distance two and an induced labelled graph of actual tournament unit seams. The \(3k\)-vertex tournaments here are ordinal sums, not strong. Relabelling quotients can collapse the labelled routes and are not part of this positive result.

For proper three-colorings, a length-six path has 22 possible extensions when port colors agree and 21 when they differ (the transfer matrix is \(J-I\), whose sixth power is \(I+21J\)). The parallel routes independently multiply these numbers. The ratio is consequently
\[
 \left(\frac{22}{21}\right)^{k-1}.
\]
This gives unbounded soft equality bias at literal tournament seams. It is an analogue of the repository's distance-two soft gadgets, but here the physical flips and graph vertices are fully specified.

More generally with three colors and a path of length \(L\), the unequal-port extension count is the exact Jacobsthal number
\[
 J_L=\frac{2^L-(-1)^L}{3},
\]
and the equal-port count is \(J_L+(-1)^L\), directly from \((J-I)^L\) (where \(J\) in the matrix expression is the all-ones matrix). There is no identification with a Hamilton path count: these are boundary-color extension counts on a seam subgraph. Any prescribed finite list of lengths \(L_i=4d_i+2\) with integers \(d_i\ge1\) occurs at the same kind of distance-two ports. In the binary cube turn on \(d_i\) fresh auxiliary coordinates, flip the port coordinate once, then turn those auxiliaries back off. The path is simple, and using disjoint auxiliary sets for the different arms gives internal disjointness. Subdivision gives \(L_i\) unit edges. Since selected midpoints have exactly their selected binary endpoints as neighbors, these are induced parallel arms and their Jacobsthal boundary multiplicities multiply.

### The ambient low-height fibre is a parity incidence graph

More generally the full labelled subgraph on heights 1 and 3 is the incidence of all transitive vertex orders with cyclic rotations of their consecutive triples. A tournament of height 3 has one strongly connected directed triangle and all its other SCCs singleton: SCC counts multiply, and a tournament strong component of order \(m\) has at least \(m\) paths by cutting a spanning cycle. Starting from a transitive order, reversing a consecutive arc leaves it transitive; reversing an arc spanning exactly three positions produces the ordinal directed triangle; reversing a longer arc creates a strong SCC on the entire interval, hence no height 3. Each directed triangle has exactly three transitive restorations, the three cyclic rotations of its order. Thus two-edge moves on the height-one part act by neighboring 3-cycles on permutations. They generate the alternating group (induct on the last position and conjugate its adjacent 3-cycle by the alternating group on earlier positions). Hence for \(n\ge3\) this incidence graph has exactly the two parity components of the transitive orders. The soft amplifier above is a coordinate subgraph built by taking independent disjoint local rotations in this ambient fibre.

## 2. An exact Kempe connection formula, and why hard gadgets do not transport

Any unit-height seam graph, including its isomorphism-class quotient, is bipartite by \(H\bmod4\). For a finite bipartite graph and distinct ports \(x,y\), fix their colors at \((1,1)\) or \((1,2)\), with any fixed palette size at least 3, and write \(N_=\), \(N_{\ne}\) for the extension counts. On every proper coloring where the ports are separated in the subgraph on colors 1 and 2, swap 1 and 2 in the component of \(y\). This is an involution exchanging the two boundary prescriptions. It leaves unmatched exactly those colorings in which the two ports share a Kempe component. In that component all paths alternate colors, so the possible unmatched sign is determined solely by the ports' bipartition parity. Therefore \(N_=-N_{\ne}\) is the signed number of Kempe-connected colorings with the favored prescription: positive for same-side connected ports and negative for opposite-side connected ports. Strictness follows using a two-coloring of their connected component. Different connected components give zero difference.

In particular no bipartite graph can be a hard three-color equality gadget at distinct ports. Opposite sides are distinguished already by a two-coloring; same-side ports can be distinguished by recoloring one port with color 3. Likewise two nonadjacent ports cannot be forced unequal: same-side ports agree in a two-coloring, and if on opposite sides recolor both with color 3. An adjacent pair does force inequality for the trivial edge reason. In a literal seam graph the sign of a nonzero soft bias is thus prescribed by \(H\bmod4\), not a freely switchable gadget parameter.

## 3. A strong seven-vertex obstruction to a planar Hamilton-height coordinate

Even allowing extra unit contacts, one cannot in general embed a unit-height seam graph injectively in the plane with the Hamilton count encoded by an affine linear coordinate. In a circle of radius one about the image of the parent, at most two points have any given projected height. The following strong seven-vertex tournament \(H=71\) has three pairwise nonisomorphic strong one-flip children all of count 69. All masks and indices are zero-based.

| Role or reversed pair | Outgoing masks of the displayed representative |
|---|---|
| Center | \((4,9,66,37,15,23,59)\) |
| \((1,2)\) | \((2,64,11,35,15,23,61)\) |
| \((2,3)\) | \((8,5,33,70,15,27,55)\) |
| \((5,6)\) | \((4,9,34,69,15,27,55)\) |

The child representatives are obtained from the literal flipped masks by the respective transpositions \((1,2),(2,3),(5,6)\). The first child has a different score multiset from the other two. The latter have a unique outdegree-one vertex, whose sole out-neighbor has outdegree 3 and 2 respectively. Thus this is a genuine obstruction even after isomorphism quotienting. The companion replay checks all counts, incidences and strongness.

The full labelled orientation cube itself does have faithful rational Euclidean embeddings by the unit-vector construction above using only coefficients \(0,1\). Conversely in an injective planar realization of all cube edges as unit distances, any 4-face is a rhombus: the other two points are distinct intersections of the circles about opposite corners, hence are symmetric about their midpoint. Opposite edge vectors are equal and each cube coordinate has one fixed unit vector everywhere; the placement is therefore additive. A linear Hamilton-height coordinate would then make \(H\) an affine function of the bits. But \(H(S)=H(S^\complement)\) by reversing every Hamilton path in the converse tournament, so such an affine function would have all coefficients zero. On \(n\ge3\), a transitive tournament and a single directed triangle in an ordinal sum give heights 1 and 3, contradiction. The positive Euclidean construction is a transport of the adjacency relation, not a projection of Hamilton height.

## 4. Inherited context and scope

The comparison is with the [induced Euclidean distance-two amplifiers](../session254_four_frontiers_20261111/hadwiger_nelson/PROOF.md) and the [hard equality gadget](../session269_centres_ports_switches_20261111/hadwiger/PROOF.md). The hard relation does not transfer to unit-height seams; the soft relation does, with an exact color-count mechanism and an explicit family of planar distance-two realizations. The Jacobsthal link is an exact boundary-color count, not a claim about arithmetic Collatz dynamics, Hamilton height recurrences, or plane realizability of all tournament flip quotients.

The hard port obstruction has an exact chromatic formulation. For a \(k\)-colorable abstract graph \(G\) and nonadjacent ports \(u,v\), hard equality is equivalent to \(\chi(G+uv)=k+1\). Hard inequality is equivalent to \(\chi(G/(u=v))=k+1\), where the two ports are identified. Adding the edge encodes unequal colors; identifying the ports encodes equal colors. In either operation an additional color suffices, which makes the threshold exact. If \(\chi(G)\le k-1\), a spare color permits either relation: recolor one port for inequality or both nonadjacent ports for equality. Thus a bipartite literal ±2 graph cannot directly force a hard equality for any palette \(k\ge3\), regardless of how strongly its extension counts are biased. Neither the edge addition nor the identification is presumed to preserve planar unit-distance realizability.
