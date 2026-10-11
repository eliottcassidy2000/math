# Graph unions: mixed palettes, exact colour contexts, and the missing geometry

**Status: PROVED elementary extensions and interfaces; CITED graph-existence input; FINITE-EXACT controls.** The supplied paper is used for its stated graph results, not presented as a new theorem of this workspace. No unit-distance realization of its extremal graphs, universal plane six-colouring, or six-versus-seven decision is obtained.

## 1. Source and inheritance

The attached five-page paper is **Wouter Cames van Batenburg, Maria Chudnovsky, Linda Cook, James Davies, Seokbeom Kim and Sang-il Oum, _On the chromatic number of the union of comparability graphs_, August 13, 2026**. The inspected source is [the supplied PDF](C:/Users/Eliott/Downloads/on_the_chromatic_number_of_the_union.pdf); the [primary arXiv record](https://arxiv.org/abs/2606.09415) identifies the paper. Earlier search metadata lists five authors; the attachment itself has the six-author heading above.

Page 1 defines an **edge union** and states Theorem 1: for every positive $k,d$, some union of $d$ comparability graphs has clique number $k$ and chromatic number $k^d$. The exponent was checked visually on pages 1–2. Theorem 2 on page 2 makes the same value equal to the fractional chromatic number and inverse independence ratio. Page 2 constructs the comparability layers by first differing coordinates of a product colouring of a high-girth graph. Page 3 distinguishes the explicit construction for Theorem 1 from the requested constructive version of Theorem 4; pages 4–5 give the probabilistic appendix. These are source statements, not new computational claims here.

The closest inherited interfaces are:

- [Marked unit-triangular carriers](unit_rhombus_interfaces_20261011c.md): a **vertex cover with retained geometric addresses** gives two disjoint three-colour palettes.
- [Complete boundary-colouring relations](boundary_palette_compression_20261011.md): gluing needs compatible boundary partitions, not merely the two chromatic numbers.
- [High-girth plane crowding](plane_girth_crowding_20261011b.md): an abstract graph operation owes a separate geometric realization; fixed-resolution geometry adds restrictions absent from abstract graphs.

The hostile is replacing an edge-union product bound by the vertex-cover sum bound. The corrected near miss is interpreting a quotient colouring obstruction as an obstruction to recolouring the original graph. The least-used marker is the **actual joint colour state of each original vertex**, together with all original edges. The live board is edge ownership / palette products / overlap relations / first-difference paths / geometric realization.

## 2. Three different union operations

Let all graphs be finite and simple.

**Edge union.** If $G=G_1\cup\cdots\cup G_d$ and each $G_i$ has a proper colouring $c_i$ with $k_i$ colours, add isolated vertices to make their labelled vertex sets agree. Then

\[
v\longmapsto(c_1(v),\ldots,c_d(v))
\quad\text{gives}\quad
\chi(G)\le\prod_i k_i.                                      \tag{1}
\]

Every edge belongs to a component in which the corresponding coordinate changes. No agreement of the individual colour names on overlaps is required. The bound is sharp even without the paper: on the product set $\prod_i[k_i]$, let $G_i$ be the complete multipartite graph joining exactly pairs with different $i$-th coordinates. Their union is complete. In particular the $2\times3$ example needs six colours, exceeding $2+3=5$; the $3\times3$ example needs nine, exceeding six. Each component is a comparability graph, oriented by increasing coordinate. Their union clique is the full product, so these small examples do not reproduce the paper's stronger fixed-clique assertion.

**Vertex cover.** If sets $V_i$ cover $V(G)$ and the full induced graph $G[V_i]$ has a supplied $k_i$-colouring, assign every vertex to one of its carriers and use the tagged colour $(i,c_i(v))$. This gives

\[
\chi(G)\le\sum_i k_i.                                      \tag{2}
\]

Cross-carrier edges are handled by disjoint palettes. This is the operation in the inherited two-unit-triangular-carrier theorem. Coordinates, scale, carrier membership and the formula $a-b\bmod3$ authenticate that theorem's hypotheses. An edge-layer decomposition does not supply a vertex-carrier cover.

**Same-palette gluing.** If two fragments meet on labelled boundary $B$ and their interiors are disjoint, a common $k$-colouring exists exactly when their complete boundary partition relations have a common member, as proved in the inherited boundary note. Merely knowing both fragments are bipartite is insufficient: paths of lengths two and three with just their endpoints identified give $C_5$. The first path forces equal endpoint colours, the second forces unequal ones. This is an overlap constraint, not a failure of (1).

These operations also apply as certificates for plane unit-distance graphs, but only when the input graph contains every required unit edge. Under the accepted lower bound of six and the usual upper bound of seven, a $2\times3$ edge-colouring certificate or a $3+3$ vertex-carrier certificate would settle six on its authenticated source. Neither finite certificates nor the abstract union theorem supply those hypotheses for the whole plane.

## 3. Mixed palettes with a sharper girth requirement

Here is an elementary extension of the paper's page-2 construction. No literature-priority claim is made.

Let $d\ge2$, let $k_1,\ldots,k_d\ge2$, and give a graph $G_0$ a proper product colouring

\[
\phi:V(G_0)\longrightarrow\prod_{r=1}^d[k_r].
\]

Assign each original edge to its first differing coordinate $r$ and orient it in the increasing direction of that coordinate. In each layer take the transitive closure; write its underlying comparability graph as $H_r$, and put $G=\bigcup_rH_r$. An edge of $H_r$ is certified by an original $r$-monotone path. All earlier coordinates are constant along the path, its $r$-th coordinate increases, and its length is at most $k_r-1$.

Thus every closure edge still has its **unique first-difference owner**. The $r$-th coordinate is a proper $k_r$-colouring of $H_r$, and $G_0\subseteq G$.

**Mixed-palette theorem.** Put

\[
g_*=1+\max_{r<s}\{\,2(k_r-1)+(k_s-1)\,\}.                   \tag{3}
\]

If $G_0$ has girth at least $g_*$, then

\[
\omega(G)\le \max_r k_r,\qquad
\chi(G_0)\le\chi(G)\le\prod_r k_r.                          \tag{4}
\]

If in addition $\chi(G_0)=\prod_rk_r$, then

\[
\chi(G)=\prod_rk_r,\qquad \omega(G)=\max_rk_r.              \tag{5}
\]

**Proof.** For any triple of distinct product labels, either all three first-difference owners agree, or they have the form $r,r,s$ with $r<s$. Indeed, at the first coordinate where the three labels are not all equal, either all three values differ, or two agree and only their mutual pair first differs later.

A triangle of $G$ whose owners are not all equal therefore has two original paths of owner $r$ and one of owner $s$. The $s$-path is edge-disjoint from both $r$-paths, because an original edge has only one owner. The union of the two $r$-paths connects the endpoints of the $s$-path. It follows that their total union contains a cycle, even if the paths meet at additional vertices. Its length is at most $2(k_r-1)+(k_s-1)<g_*$, a contradiction.

Consequently every triangle of $G$ is contained in one layer. Every clique of size at least three is then entirely in one layer: fix an edge, and use the triangles through it to identify the owners of all other edges. A clique in a comparability graph is a chain, so its size is at most that layer's $k_r$. Edges and isolated vertices satisfy the same bound because $k_r\ge2$. This proves the clique bound and (4).

If the source chromatic number is the full product, every $H_r$ must have chromatic number $k_r$, since a smaller value in any coordinate would give a smaller product colouring of $G$. For a finite poset its comparability graph can be coloured by longest-chain level and has a clique of maximum-chain size. Thus each $H_r$ has a chain of size $k_r$, giving the reverse clique inequality in (5). ∎

For nondecreasing palettes, (3) becomes

\[
g_*=k_d+2k_{d-1}-2.
\]

In particular **palettes $(2,3)$ require only girth five**; the order $(3,2)$ gives the sufficient threshold six. Equal palettes recover the paper's $3k-2$. The one-layer case is separate: there is no mixed triangle, and every such transitive closure has clique number at most $k_1$ without a girth hypothesis.

The classical high-girth prescribed-chromatic existence theorem quoted as Theorem 3 on page 2 supplies a suitable finite $G_0$ with chromatic number $\prod_i k_i$. Therefore (5) yields abstract mixed-palette extremal graphs. In particular there exists a union of a bipartite comparability graph and a height-three comparability graph with $\chi=6$ and $\omega=3$. This is an existence consequence of the cited input plus the proved extension, not an explicit six-chromatic unit-distance graph.

### The geometry destroyed by closure

Transitive completion is not a unit-distance-preserving operation. The collinear points $0,1,2$ have the unit path $0-1-2$. With increasing one-coordinate labels, its closure adds $02$, whose retained Euclidean length is two. The receipt stores the original two-edge path; it does not relabel its endpoints as unit distance.

Thus the paper's abstract $\omega=3,\chi=9$ example cannot be transferred to the plane's unit graph merely because both are geometric-looking graph constructions. Likewise, the mixed $\chi=6$ example still needs exact coordinates and every added edge checked if proposed as a plane witness. Clique bound alone does not supply those coordinates.

## 4. Exact compression relative to supplied colourings

For supplied component colourings in (1), let $s(v)$ be the joint colour tuple. Let $S$ be the set of tuples that actually occur, and define a graph $Q$ on $S$ by

\[
st\in E(Q)
\quad\Longleftrightarrow\quad
\exists uv\in E(G)\text{ with }s(u)=s,\ s(v)=t.              \tag{6}
\]

There are no loops, by (1). Keep the original vertex-to-state map and the original source graph; these determine all edges of $Q$.

**Fixed-state compression theorem.** For any map $f:S\to[q]$,

\[
f\circ s\text{ is proper on }G
\quad\Longleftrightarrow\quad
f\text{ is proper on }Q.                                   \tag{7}
\]

Both directions follow by reading the witnessing original edge in (6). Thus $\chi(Q)$ is exactly the least palette size obtainable by recolouring **only as a function of the retained joint tuple**. In particular two states can be identified while leaving every other state distinct iff they are nonadjacent in $Q$.

**Saturation corollary.** If $\chi(G)=P=\prod_i k_i$, then every collection of proper component $k_i$-colourings uses all $P$ joint tuples and has $Q=K_P$. Indeed $G\to Q$ gives $P\le\chi(Q)\le|S|\le P$, and a graph on $P$ vertices has chromatic number $P$ only if complete. In the paper's extremal examples this can happen while $\omega(G)=k\ll P$: every pair of states has some conflict, but its witnessing edges need not share representatives that form one clique.

The converse conclusion about the original graph is false. On $C_8$, list the joint states in cyclic order as

\[
0,1,0,2,1,3,2,3,\quad
0=00,\ 1=01,\ 2=10,\ 3=11.
\]

Put edges that change the first bit in the first bipartite layer, and all remaining edges in the second. These are proper binary component colourings. Every pair of the four states occurs on an edge, so $Q=K_4$, but the source $C_8$ is bipartite. Four colours are necessary only for functions of these particular states; a fresh source recolouring uses two. This is the sharp information-loss boundary for this interface.

There is also a tight source control: the 11-vertex Mycielski graph of $C_5$ has $\omega=2,\chi=4$, and a product $2\times2$ colouring. Its height-two closures add no edges, and its state-conflict graph is $K_4$. The implementation independently rejects every three-colouring. The C8 example separates accidental codec saturation from saturation forced by the source chromatic number.

This interface differs from the inherited **complete boundary relation**: here the local colourings are fixed and the decoder is restricted to their tuple fibres. Boundary tables instead retain every admissible partition and permit contextual recolouring. A negative result in one fixed-state quotient does not reject all those alternatives.

## 5. Receipt implementation and controls

The [script](../../04-computation/experiments/graph_union_context_20261011d.py) is self-contained. **compile_union** accepts the exact source graph, ordered palette sizes and one product label per source vertex. Every added edge carries a monotone path made entirely from original edges. **verify_packet** checks all paths, source-edge inclusion, transitive completeness and palette propriety; thus omitted closure edges and forged extra edges are both rejected. The girth theorem is stated only after checking its separate hypothesis.

**joint_quotient** retains the source graph and projection. **lift** reconstructs the complete conflict graph before consuming a proposed state colouring, then authenticates the resulting source colouring. It does not claim that every source colouring factors through that projection. Endpoints, palette entries and colours must be exact integers; Python boolean/float aliases are not accepted. Graph endpoints are checked before conversion to a set, so an invalid duplicate cannot disappear through equality aliasing.

The complete finite universe is all graphs with compatible $(2,3)$ product labels on $0,\ldots,4$ vertices: **49,885 packets**, of which **34,555** satisfy the sufficient girth condition. Additional controls include all **3,120** proper product colourings of $C_5$, the exact Mycielski control, every recolouring of the four C8 states using palettes $1,\ldots,4$, complete multipartite product extremals, overlap-parity and metric hostiles, and typed malformed packets. These do not enumerate all graphs of arbitrary order.

~~~text
python -B -X utf8 04-computation/experiments/graph_union_context_20261011d.py
python -B -O -X utf8 04-computation/experiments/graph_union_context_20261011d.py
~~~

Normal, optimized and [saved output](graph_union_context_20261011d.out) agree at **140,950 exact checks**, with LF SHA256 **f389d7fe95161f6259aa22317b36d6406a6673b32552f53b4ecd37fdd6b61e9c**. The mathematical extension is the mixed-palette theorem and its marked path compiler; the source theorem, product colouring and basic graph quotient are inherited mechanisms with explicit provenance. The actionable plane target remains producing valid source-wide decompositions or compatible geometric boundary relations, not merely small local chromatic numbers.
