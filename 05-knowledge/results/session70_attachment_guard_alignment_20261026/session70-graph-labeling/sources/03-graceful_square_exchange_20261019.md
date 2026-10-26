# Graceful and square-sum paths: which exchanges actually work?

**Status:** PROVED (elementary exchange criteria and obstructions, with no external-priority claim); FINITE-EXACT (fully enumerated reconfiguration graphs and the sharp `28 -> 29` surgery cost); OPEN (two explicitly scoped connectivity conjectures). The Graceful Tree Conjecture is not proved. General square-sum existence is not being proposed as open: the repository already cites Gerbicz's construction for every `n >= 25`, and square-sum cycles for every `n >= 32`.

Companion: [main checker](../../04-computation/experiments/graceful_square_exchange_20261019.py); [deterministic output](graceful_square_exchange_20261019.out); [independent checker](../../04-computation/experiments/graceful_square_exchange_20261019_audit.py); [audit output](graceful_square_exchange_20261019_audit.out); [order-14 extension](graceful_square_exchange_20261019.path14.out). The standard run checks graceful paths through 13 vertices, gracefully labeled trees through 11, and square-sum paths through 36. A separately recorded `--graceful-max 14` extension is specified below. All path counts identify reversal, **not** label complementation.

## 1. Inheritance and the actual missing operations

The inherited mechanisms are:

- **THM-2761, graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge:** sum/difference collision factors; the bipartite sign gauge does not construct a graceful labeling.
- **THM-2787, signed-path-sum-weyl-orbit-and-gap-tail-leaf-insertion:** graceful labelings are bounded, nonzero signed-path integrations; a marked leaf can be inserted across a rank gap precisely when its crossing edges form a difference suffix. Its all-tree marked extension remains OPEN.
- **THM-2793, oriented-ramp-reference-and-rooted-graceful-gap-reconstruction:** edge sums are the missing first-moment sidecar beside absolute differences; no existence follows from reconstruction.
- **THM-2798, star-graceful-defect-minimax-barrier-for-bounded-window-moves:** a fixed-width vertex-order grammar incurs an unbounded collision barrier on stars.
- **THM-4505, small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet**, and `crossroads_family_20260926_squares.md`: the rigid first square-sum path, the first square-sum four-cycle/two-edge cycle switch at 46, and the difference between a selected path and its ambient graph.
- `collatz_mod6_20260922_w6_square_sum_hamiltonicity.md`: exact path counts, Gerbicz's actual construction, and the endpoint-aware correction to degree-two forcing. `decoder_prime_square_20260925.md` repairs the former endpoint-blind argument.

The closest proved operation is the gap-tail insertion. The canonical hostile is the star's fixed-width barrier; the corrected near miss is forcing both edges at a degree-two **possible endpoint**. The neglected coordinate here is whether an operation keeps one endpoint, both endpoints, neither endpoint, or the underlying tree shape.

Concept board: (i) signed-path integration, (ii) fixed-endpoint versus moving-endpoint exchanges, (iii) a fundamental cycle and its edge colors, (iv) interval-incidence bases, (v) square-sum insertion seams. This note replaces a generic sum/difference analogy with exact operations and sharply different failure boundaries.

### Literature discipline

Endpoint rotations are classical Pósa rotations, not a new operation. A primary exposition is Michael Krivelevich, *Long paths and Hamiltonicity in random graphs*, section 4.2, https://www.math.tau.ac.il/~krivelev/HamGnd.pdf (the indexed definition was checked; direct PDF retrieval timed out during this session). Graceful permutations and their enumeration are established subjects: Michal Adamaszek, *Efficient enumeration of graceful permutations*, https://arxiv.org/abs/math/0608513. The Walecki construction is classical. Rainbow spanning-tree reconfiguration is also an existing subject: Shun-ichi Maezawa and Yutaro Yamaguchi's 2024 lecture abstract, https://www.en.sdu.edu.cn/info/1158/8226.htm, states general nonconnectedness and a sufficient condition for connectivity. The abstract does not state that condition; its applicability to the interval-difference specialization here has not been checked. GT1 is open in this packet, pending that literature comparison, not a claim of a newly introduced general framework. No novelty claim is made for the elementary exchange lemmas, the involution, or the existence of graceful path labelings. The contributions here are the combined operation audit, explicit obstructions, exact reconfiguration censuses, and precise next hypotheses.

## 2. One-edge graceful exchanges: midpoint rotations plus a double-end splice

Represent a graceful path by a permutation `p=(p_0,...,p_(n-1))` of `0,...,n-1` whose consecutive absolute differences are `1,...,n-1`. Paths are identified up to reversal. Numeric labels remain fixed; changing the path's edges is equivalent to changing the graceful permutation of the abstract path.

### 2.1 Midpoint criterion and frozen long edges (PROVED)

A Pósa rotation with fixed first endpoint, pivot `v=p_i`, removed neighbor `w=p_(i+1)`, and moving endpoint `x=p_(n-1)` produces

`(p_0,...,p_i,p_(n-1),p_(n-2),...,p_(i+1))`.

Only `{v,w}` is removed and `{v,x}` added. The new path is graceful iff

`|v-x|=|v-w|`, equivalently `x+w=2v`.

The other solution `x=w` is excluded by a nontrivial rotation. Thus every rotation reflects the endpoint through its pivot and preserves the parity of the moving endpoint. Moreover, if `d=|v-x|`, the distinct labels `v-d,v,v+d` lie in `[0,n-1]`, so `2d <= n-1`. Consequently every edge longer than half the label range is literally fixed by every sequence of graceful midpoint rotations.

This is a barrier even to an arbitrarily long reversed suffix. It is different from THM-2798's bounded consecutive support. It is **not** an invariant of arbitrary one-edge exchanges: the next move can change both endpoints.

### 2.2 A universal double-end splice (PROVED)

For `n>=3`, close a graceful path's endpoints `a,b` by adding `{a,b}`. Let `d=|a-b|`. There is exactly one old edge `e_d` of length `d`. Delete `e_d` from the newly formed spanning cycle.

The result is another graceful Hamiltonian path: its vertex set is unchanged, deleting one cycle edge leaves a spanning path, and one length-`d` edge replaces another. Its endpoints are the ends of `e_d`, hence its endpoint gap is again `d`.

Applying the move again restores `e_d` and deletes `{a,b}`. The map is therefore a fixed-point-free involution on unoriented graceful paths, preserving the endpoint-gap fiber. It cannot be fixed because the old endpoints were not adjacent in a path with at least three vertices.

The involution alone is a perfect matching, not a connected search algorithm. It also pairs each endpoint-gap fiber into even cardinality, but that parity is **not a new counting conclusion**: ordinary label complementation gives another fixed-point-free pairing modulo reversal. Indeed a complement-reversal fixed path with `n>=3` would have equal differences at mirrored edge positions, contradicting gracefulness. The useful content of the splice is its local edge-exchange realization.

### 2.3 Exhaustion of all one-edge path exchanges (PROVED)

Remove an edge from a path, yielding two path segments. Each has two ends. Every replacement edge making a spanning path must join one end of each segment. Apart from restoring the original edge, the possibilities are a rotation retaining one original endpoint, or joining the two original endpoints (the double-end splice).

Thus **every** one-edge graceful path exchange is a midpoint rotation or the double-end splice. Degenerate choices at an original endpoint merely coincide; they do not introduce another operation. The exact `E1` census below is therefore also the exact connectivity census of this complete one-edge grammar.

### 2.4 Walecki's exact mod-three rotation boundary (PROVED)

The only graceful path having endpoint 0 is

`W_n=(0,n-1,1,n-2,2,n-3,...)`.

Proof: the largest difference forces the first edge `0--(n-1)`. Delete 0; the remaining path uses labels `1,...,n-1` and differences `1,...,n-2`, with an extreme endpoint. Translate and complement to apply induction. This gives the stated zigzag uniquely.

A nontrivial rotation retaining endpoint 0 is therefore impossible. For a rotation moving endpoint 0 and retaining the other endpoint, the pivot must be a lower zigzag value `v`; its predecessor in the original order is `n-v`. The midpoint equation is `2v=n-v`. Hence

`W_n has a graceful midpoint rotation iff 3 divides n`,

and then it has exactly one: remove `{n/3,2n/3}` and add `{0,n/3}`.

For every `n>=2` with `3` not dividing `n`, Walecki is isolated under midpoint rotations **even if arbitrary endpoint-preserving graceful rearrangements are added**: any such rearrangement would retain endpoint 0, where uniqueness forbids change. The double-end splice escapes this obstruction for `n>=3`. This identifies a missing endpoint operation, rather than only reporting disconnected finite search graphs.

## 3. Why fixed-endpoint graceful switches need three edges

### 3.1 No nontrivial two-edge switch (PROVED)

For four distinct labels `a<b<c<d`, the three perfect matchings have length multisets

- `{b-a,d-c}`;
- `{c-a,d-b}`;
- `{d-a,c-b}`.

The first has strictly smaller total than either other matching; the latter two have equal total but different maxima (`d-a` is strictly largest). All three multisets are therefore different.

Two Hamiltonian paths with the same unordered endpoint set have the same degree at every vertex. If exactly two edges are replaced, their symmetric difference must be an alternating four-cycle, hence two distinct matchings on four distinct labels. Those matchings cannot have the same length multiset. Thus distinct graceful paths with fixed endpoints cannot differ by a two-edge switch, at **any** size.

Compare the square-sum result: fixed-endpoint square-sum two-edge switches need a square-sum four-cycle, first possible at `n=46`. The graceful matching obstruction never disappears with height. Allowing endpoints to change is essential; two-edge exchanges with changed endpoints need not be alternating four-cycles.

### 3.2 Sharp three-edge positive examples (PROVED by displayed certificates)

Already on five labels,

`P=(1,4,0,2,3)`, `P'=(1,2,4,0,3)`

are graceful with endpoints `{1,3}`. Remove `{1,4},{0,2},{2,3}` and add `{1,2},{2,4},{0,3}`; the common edge is `{0,4}`. The changed length multiset is `{1,2,3}` on both sides. The length-3 edge changes, escaping the midpoint long-edge invariant. The census through four vertices shows five is the smallest vertex count for such a graceful path trade.

There is also a clean disjoint-matching trade on six labels:

`(2,0,5,1,4,3) -> (2,1,4,0,5,3)`.

The exchanged matchings are `{(0,2),(1,5),(3,4)}` and `{(0,4),(1,2),(3,5)}`, each of lengths `{1,2,4}`. Their union is an alternating six-cycle; common edges have lengths 3 and 5.

A constructive family explains this example. For six distinct labels arranged in complementary pairs `u,L-u`, `a,L-a`, `b,L-b`, replace a contiguous segment

`(u,a,L-a,b,L-b,L-u)`

by

`(u,b,L-b,a,L-a,L-u)`.

Endpoints and vertex set agree. Both segment difference multisets are

`{|u-a|, |L-2a|, |L-a-b|, |L-2b|, |u-b|}`.

Hence the move preserves any ambient graceful path. More generally, reverse a segment and complement each of its labels by `x -> L-x`, when its label set is closed under this reflection and its endpoints sum to `L`. This leaves the endpoint labels and the segment difference multiset unchanged. These operations need not connect all labelings.

## 4. Two connectivity hypotheses with explicit finite support

### 4.1 Three-edge graceful-path exchange conjecture (OPEN)

Let `E_r(n)` have as vertices all graceful paths on the **fixed** label set `0,...,n-1`, modulo reversal. Join two when at most `r` old edges are removed and the same number added. Endpoints may change.

**Conjecture GP3:** `E_3(n)` is connected for every `n>=2`.

Complete enumeration establishes this through `n=13` in the standard companion, and through `n=14` in the separately recorded extended run (4,956 paths). Radius two already fails at `n=8` with two components of ten paths, and fails again at 11,12,13,14. Thus replacing three by two is false.

| n | unoriented paths | midpoint components | E1 components | E2 components | E3 components |
|---|---:|---:|---:|---:|---:|
| 2 | 1 | 1 | 1 | 1 | 1 |
| 3 | 2 | 1 | 1 | 1 | 1 |
| 4 | 2 | 2 | 1 | 1 | 1 |
| 5 | 4 | 4 | 2 | 1 | 1 |
| 6 | 12 | 4 | 2 | 1 | 1 |
| 7 | 16 | 12 | 6 | 1 | 1 |
| 8 | 20 | 20 | 10 | 2 | 1 |
| 9 | 60 | 40 | 14 | 1 | 1 |
| 10 | 148 | 102 | 42 | 1 | 1 |
| 11 | 324 | 214 | 92 | 6 | 1 |
| 12 | 664 | 520 | 220 | 2 | 1 |
| 13 | 1600 | 1132 | 448 | 2 | 1 |
| 14 | 4956 | 3364 | 1286 | 2 | 1 |

The table does not imply a bounded graph diameter, an efficient algorithm, or any statement about other tree shapes. In particular, midpoint rotations plus **fixed-endpoint** three-edge trades still fail to connect all paths at several small sizes; GP3 deliberately permits both endpoints to change.

### 4.2 Graceful-tree one-color exchange conjecture (OPEN)

Color every edge `{a,b}` of `K_n` by `|a-b|`, with labels `0,...,n-1`. A gracefully labeled tree is exactly a rainbow spanning tree with colors `1,...,n-1`.

For any nonedge `{u,v}` of such a tree, let `d=|u-v|` and let `e_d` be its unique old color-`d` edge. Replacing `e_d` by `{u,v}` is another gracefully labeled tree **iff `e_d` lies on the original tree path from u to v**. This is exactly the fundamental-cycle exchange condition, and its path version is the double-end splice.

**Conjecture GT1:** The graph of all these gracefully labeled trees, with one-color exchanges as edges, is connected for every `n>=2`.

The underlying tree shape is allowed to change. This is not the Graceful Tree Conjecture: it asserts connectivity only among the tree labelings that already exist, without asserting that every unlabeled tree shape occurs.

The complete finite counts are:

| n | gracefully labeled trees | exchange components |
|---|---:|---:|
| 2 | 1 | 1 |
| 3 | 2 | 1 |
| 4 | 4 | 1 |
| 5 | 12 | 1 |
| 6 | 40 | 1 |
| 7 | 164 | 1 |
| 8 | 752 | 1 |
| 9 | 4020 | 1 |
| 10 | 23576 | 1 |
| 11 | 155632 | 1 |

**Exact interval-matrix reformulation (PROVED).** Encode a color-d edge as `{a_d,a_d+d}`, where `0<=a_d<=n-1-d`. Relative to the oriented unit path-incidence basis, its column is the interval indicator `1_[a_d,a_d+d)` in `Z^(n-1)`. Selecting one interval of each length `1,...,n-1` gives a square matrix. The selected graph is a spanning tree iff this matrix is nonsingular; its determinant is then `+1` or `-1`, since a tree's reduced incidence matrix is unimodular. Otherwise its determinant is zero. Thus GT1 is exactly connectivity of the nonsingular one-interval-per-length matrices under replacement of one column's start. A move may change that start by more than one; no unit-slide connectivity claim is made.

The interval/difference structure matters. In the proper three-color one-factorization of `K_4`, the only rainbow spanning trees are the four stars. Each is isolated under one-color exchange: the other edge of the same color is disjoint from the removed star edge and does not reconnect the isolated leaf. Therefore arbitrary rainbow spanning-tree systems need not satisfy GT1, even though every color class is a forest.

**Acyclic upper-tail completion (PROVED).** Any acyclic selection of one edge of each length `k,k+1,...,n-1` extends to a gracefully labeled tree. At the next stage `d=k-1`, the current forest has `d+1` components. The graph consisting of every length-d edge on `0,...,n-1` has exactly d components, the residue chains modulo d. If all length-d edges were internal to the forest's components, each residue chain would lie inside one forest component, forcing at most d forest components, a contradiction. A length-d edge therefore crosses the forest cut and can be added without creating a cycle. Continue down to d=1. This is a useful construction of a graceful tree containing any prescribed acyclic upper tail; it preserves neither a requested tree shape nor a previously chosen lower tail, and does not prove GT1.

**Shortest-edge buffer (PROVED).** Every gracefully labeled tree on `n>=3` admits another length-one exchange. Remove its length-one edge. The two components have 0 and n-1 on the same side, because their length-(n-1) edge remains. Reading component membership along the ordered labels 0,1,...,n-1 therefore crosses the cut a positive even number of times. Each crossing adjacent pair is a possible length-one replacement. There are at least two choices, and all reconnect the same two components. Consequently every nonempty fiber obtained by fixing all edges of lengths 2,...,n-1 is an even clique of size at least two. In particular GT1 has no isolated vertices. This is an intrinsic interval-order feature, but does not prove that the cliques connect to one another. The audit checks these fibers through n9.

A plausible descent is already false. For the code `(a_1,a_2,a_3,a_4)=(2,0,1,0)`, the edges are `{2,3},{0,2},{1,4},{0,4}`. Every attempt to reduce any positive `a_d` alone destroys the tree property. This is a non-star local minimum of `sum a_d`; the star at 0 is not reachable by a monotone coordinate descent from it. The exact counts of such decreasing-coordinate local minima are 2 at n5, 4 at n6, 10 at n7, 27 at n8, 87 at n9, 312 at n10, and 1207 at n11. Connectivity needs reversible detours or another potential.

## 5. Square-sum extensions: rank gaps fail, rotations are insufficient, three cuts succeed

Let `Q_n` join distinct positive labels `x,y<=n` when `x+y` is a square.

### 5.1 Why the graceful rank-gap operation does not transfer (PROVED)

Inserting a vacant integer t into `1,...,n` shifts each old label at least t upward by one. Every old edge incident with a shifted label increases its sum by one or two. A positive square edge sum is at least 4, and the next square is at least five larger. Such an edge is no longer square-sum.

Thus on a nontrivial path **no nonempty rank-gap shift can preserve all old square edges**. Only the outside gap `t=n+1`, with no old labels shifted, survives. In contrast, THM-2787's difference suffix is specifically designed to absorb one-unit increments. This is a target-alphabet spacing obstruction, not a failed sign gauge.

### 5.2 Exact two-cut insertion criterion (PROVED)

For any graph defined by an allowed edge relation, take a path containing, in this order, two disjoint cut edges `{a,b}` and `{c,d}`:

`(...,a,b,...,c,d,...)`.

Reverse the middle segment and insert a new vertex z:

`(...,a,c,...,b,z,d,...)`.

This is a path on the enlarged vertex set with the same endpoints iff its new edges `{a,c}`, `{b,z}`, `{z,d}` are allowed. For square sums the criterion is exactly

`a+c, b+z, z+d are squares`.

The symmetric version inserts z between a and c and uses the new old-vertex edge `{b,d}`. Conversely, every endpoint-preserving insertion with exactly two removed old edges has one of these forms: the cuts leave three old path segments; their outer endpoint segments are fixed, so the middle segment must reverse, with z at one of its two joins. If it did not reverse, one removed old edge would be restored. Adjacent cut edges give a one-vertex middle segment and cannot realize two genuinely changed old edges this way.

The criterion is not asserted for operations changing the old endpoints.

### 5.3 The first global endpoint-rotation obstruction in the viable tail: 28 -> 29

At `n=28`, all 52 Hamiltonian paths fall into seven Pósa-rotation components, of sizes `2,2,2,3,7,9,27`. None has an endpoint in `{7,20}`, the neighbors of 29 in `Q_28`. Yet `Q_29` has 19 Hamiltonian paths.

This blocks **every** sequence of valid endpoint rotations followed by appending 29, not merely one seed or a bounded number of rotations.

There is an elementary endpoint-aware explanation. Vertex 18 is a leaf, so it must be one endpoint; its neighbor 7 cannot be the other endpoint in a spanning path of more than two vertices. Suppose 20 were the other endpoint. Vertex 16 is then a certified internal degree-two vertex, forcing `20--16--9`; hence `{20,5}` cannot be used. Internal vertex 5 must use `{5,4}` and `{5,11}`. Internal vertex 25 forces `{25,11}` and `{25,24}`, so `{11,14}` cannot be used. Internal vertex 14 then forces `{14,2}` and `{14,22}`. Internal vertex 27 forces `{27,9}` and `{27,22}`. Together with the leaf edge `{18,7}`, the saturated vertex 9 excludes `{9,7}`; vertex 7 must use `{7,2}`. The forced edges already form the complete endpoint-to-endpoint path

`20--16--9--27--22--14--2--7--18`,

omitting 19 vertices. This is impossible. Every use of degree two here is guarded by the fixed endpoint set `{18,20}`.

### 5.4 A sharp three-edge surgery (FINITE-EXACT minimum; explicit positive certificate)

The following source and target have the same endpoints 1 and 18:

`P28=(1,24,25,11,14,22,27,9,16,20,5,4,12,13,3,6,19,17,8,28,21,15,10,26,23,2,7,18)`

`P29=(1,24,25,11,5,4,12,13,3,6,19,17,8,28,21,15,10,26,23,2,14,22,27,9,16,20,29,7,18)`.

Remove `{2,7},{5,20},{11,14}`. Add `{2,14},{5,11},{7,29},{20,29}`. Their respective sums are `16,16,36,49`, all squares. Every other edge is retained. The vertices are exactly `1..28` and `1..29`.

The independent-source/target comparison of all `52*19=988` pairs proves that three old edges is the **minimum possible**, even allowing the old endpoints to change. Across all 52 source paths, the smallest distance to any target has histogram

`3:21, 4:11, 5:15, 6:5`.

Thus no source path admits a zero-, one-, or two-edge surgery. A second audit independently deletes every set of at most two old edges from each source, and tries every orientation and ordering of the resulting segments together with singleton 29: none can be rejoined by square edges. The negative three-edge threshold here rests on exhaustive exact enumeration; the endpoint-rotation obstruction above has a separate elementary proof. No unbounded-n minimum-surgery claim is made.

### 5.5 Rotation components are not an existence test

| n | paths | rotation components | components with an endpoint adjacent to n+1 |
|---|---:|---:|---:|
| 25 | 10 | 2 | 1 |
| 26 | 12 | 3 | 2 |
| 27 | 35 | 4 | 2 |
| 28 | 52 | 7 | 0 |
| 29 | 19 | 3 | 1 |
| 30 | 20 | 3 | 3 |
| 31 | 349 | 1 | 1 |
| 32 | 392 | 1 | 1 |
| 33 | 669 | 4 | 3 |
| 34 | 4041 | 5 | 5 |
| 35 | 17175 | 1 | 1 |
| 36 | 12960 | 2 | 1 |

An isolated rotation component, or one that cannot append the next integer, is an obstruction to that search grammar only. The larger graph may have solutions in other components or require a different insertion seam. Gerbicz's multi-block construction already supplies a global existence route; this note measures a smaller operation's actual boundary.

## 6. Connection ledger

| Source -> target | Exact map / preserved predicate | Lost information / required sidecar |
|---|---|---|
| Graceful path -> rainbow Hamiltonian path | Color an edge by its absolute numeric difference | Keeping one endpoint, both endpoints, or neither changes the move graph |
| Rainbow tree -> fundamental-cycle exchange | Insert a color-d nonedge, remove its unique old color-d edge | Tree iff the removed edge lies on the fundamental cycle; shape need not survive |
| Difference-colored tree -> interval matrix | Express incidence columns in the unit-path basis | Length spectrum retained; nonsingularity is essential; columnwise monotonicity fails |
| Graceful endpoint rotation -> midpoint reflection | `x+w=2v` | Long edges frozen; double-end splice evades that invariant |
| Square-sum path -> local insertion seam | Reverse the middle segment and check new edge sums | Dense coverage, distinct labels, endpoint constraints, and connectedness remain explicit |
| Graceful rank-gap insertion -> square-sum insertion | **No predicate-preserving transport** | Square spacing kills every shifted old edge |

## 7. Reproduction, independence, and stopping boundary

From the directory containing the two checker scripts, run the standard companion normally and with optimization:

```
python3 graceful_square_exchange_20261019.py > graceful_square_exchange_20261019.out
python3 -O graceful_square_exchange_20261019.py > graceful_square_exchange_20261019.optimized.out
cmp graceful_square_exchange_20261019.out graceful_square_exchange_20261019.optimized.out
```

Optional path-only extension (stored as `graceful_square_exchange_20261019.path14.out`):

```
python3 graceful_square_exchange_20261019.py --graceful-max 14 --tree-max 2 --square-max 14
```

Independent surgery/shortest-edge audit:

```
python3 graceful_square_exchange_20261019_audit.py
python3 -O graceful_square_exchange_20261019_audit.py
```

The path enumeration is complete DFS over all permutations, pruning only repeated vertex labels or repeated absolute differences. It is independently compared with literal permutation enumeration through n8. Trees are enumerated by choosing one edge of each length, pruning only cycles; independent Prüfer-word enumeration agrees through n7. The square path DFS fixes a leaf as first endpoint whenever one exists; this loses no unoriented paths. Its only other pruning checks are zero residual degree, too many residual degree-one vertices, and disconnected residual graph. The counts agree with the repository's separately audited square-sum census. No incomplete search or timeout is interpreted as a proof.

Exchange components are computed by an independent common-edge signature: two paths differing in at most r edges share a signature obtained by deleting r edges from each. Thus the E1/E2/E3 census does not rely on the classification of rotations and splices that it corroborates. Square rotation components use explicit sequence transforms; the insertion-distance minimum compares literal edge sets of independently generated source/target paths. A SHA-256 of the complete Q28/Q29 path bank is printed.

GP3 and GT1 remain OPEN beyond the stated finite universes. Neither implies the Graceful Tree Conjecture as stated. The results do not identify a mechanism for Collatz, Petersen Hamiltonicity, planarity, LRC, or any other prize problem. The most useful next step is a proof or hostile for the **interval-basis exchange graph**, or a structural classification of the first radius-two graceful-path components; increasing a bounded census alone is not a substitute.
