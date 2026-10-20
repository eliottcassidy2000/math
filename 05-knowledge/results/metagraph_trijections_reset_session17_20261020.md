# Session 17: exact correspondences, fibre boundaries, and reset certificates

20 October 2026. Proof-first continuation. This document distinguishes inherited theorems, new constructions within this research session, finite checks, and open questions. No literature-priority claim is made. Companion packets contain full proofs, executable certificates, source provenance, and independent audits.

## 1. The useful interpretation of “beneath the mathematics”

An invariant is meaningful only after specifying which operations preserve it. Tournament relabelling, converse, cut switching, changing an observer order, changing a distinguished Hamiltonian path, and deleting vertices are different operations. They cannot be interchanged simply because they produce equal-sized sets.

The organizing principle of this session is an **operation-preserving correspondence**: state the source and target objects, the actual reversible map, the marking it needs, the operation it transports, and the information it discards. A failed naturality claim is itself useful when its smallest obstruction is explicit.

For an ordered labelled vertex set, a tournament is exactly the simple graph of its backwards arcs. Changing the observer order acts by symmetric difference with its inversion graph. This is an exact labelled correspondence, but not a relabelling-equivariant unmarked identification of all tournaments with all simple graphs. A transposition fixes the empty graph and fixes no tournament under its natural vertex action. At order four there are four unlabelled tournament types but eleven unlabelled simple graphs.

## 2. Why the three four-vertex shapes really differ

The twenty labelled three-edge graphs on four vertices split into twelve paths, four stars, and four triangles with an isolated vertex. The condition “one fewer edge than vertices” does not imply a tree. In general the cycle rank is m−n+c; when m=n−1 it equals c−1. Disconnection can therefore pay for cycles. The triangle-plus-isolate has one of each; the path and star are connected trees.

This difference controls tournament switching charts exactly. Prescribing arcs on a spanning constraint graph F with m edges and c components imposes one consistency condition per independent cycle. A consistent switching class contains exactly 2^(c−1) realizations. Consequently a prescribed spanning tree gives one representative of every cut-switching class, while a cyclic/disconnected constraint may miss classes and represent the surviving ones more than once.

The graph complement also matters: a free-coordinate star leaves a triangle-plus-isolate constraint, while a free-coordinate triangle-plus-isolate leaves a star constraint. There is no license to treat the three eight-state coordinate cubes as interchangeable quotient charts.

See [tree gauges and trijections](tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md).

## 3. An all-order aggregate, with a strict lifting boundary

Aggregate all labelled fixed-path complement incidences by cut-switching classes. The resulting weighted graph is a Cayley graph on edge bits modulo cuts. Every class has weighted row sum n!, because every ordered Hamiltonian word has exactly one tree-gauge representative in each class. Its generators are the cut classes of complements of spanning paths.

For n≥3 its support has no loops. Its complete component classification is:

- even n: connected and nonbipartite;
- n≡3 mod4: connected and bipartite;
- n≡1 mod4: exactly two components, each nonbipartite.

For n≥8 the different undirected paths give different generators, so the support degree is n!/2 and every positive edge has weight two. Order four gives the eight-vertex cross-polytope graph, with edge weight four.

The proof uses cut parity and the span of differences of path generators. These are all-order arguments, not extrapolation from computed pictures.

**Essential limit:** this is an aggregate of a partition, not an equitable quotient or an automorphism quotient. Consecutive aggregate edges can require different representatives in the same switching class. Aggregate paths, squares, and connectivity must not be presented as legal paths, detours, or connectivity in the original or black metagraph.

### 3.1 The exact order-five model corrects a tempting extrapolation

Each order-five component is the Hamming-distance-four graph on the 32 even six-bit vectors, not the usual halved six-cube (whose edges have distance two). Its six coordinates are the six complementary pairs of Hamiltonian five-cycles, equivalently the six Sylow-5 subgroups of S5. The exceptional doubly transitive six-point action explains why the fifteen generators become all weight-four vectors.

The explicit invertible linear map proves the graph isomorphism. Its spectrum is 15 once, 3 ten times, −1 fifteen times, and −5 six times; edge weight eight multiplies these eigenvalues by eight. Its local graph is KG(6,2), and its largest clique has four vertices. The halved six-cube instead has a six-clique and a different spectrum. Equal vertex and degree counts did not determine the shape.

See [the distance-four model](tournament_graph_distance_four_20261020/DISTANCE_FOUR_SIX_CUBE.md).

### 3.2 Fourier modes and bridges, without enormous enumeration

The dual characters of the cut quotient are exactly Eulerian edge sets Z, meaning even degree at each vertex. This degree condition is not the different automorphism-sign notion of “even graphs” in the earlier tournament-equivalence literature. The weighted eigenvalue at Z is (−1)^|Z| times the signed sum of all ordered Hamiltonian vertex words, with each word weighted by (−1) to the number of its path edges in Z.

Expanding that signed sum gives an exact sum over linear forests F contained in Z: each contributes (−2)^|F| times 2 to the number of its nontrivial components, times (n−|F|)!. Contracting the forest's path components proves the formula.

For one triangle Z, the weighted eigenvalue simplifies to

−(n−2)!(n−4)(n−9).

Thus order nine has at least 168 zero modes on the full aggregate, or at least 84 on either parity component. The complement characters explain the factor two. No enumeration of its 2^28 states is needed, and these are lower bounds, not the complete order-nine spectrum.

Every support edge of the cut aggregate lies on a commuting-generator square for n≥4, so that aggregate has no bridges. In contrast, a transitive tournament is a degree-one vertex in the original labelled path-complement graph, giving an actual bridge at every n≥3. This supplies a concrete failure of lifting aggregate cycles to a specified tournament representative. At n=3 the aggregate support is K2 with one bridge, but its weight six represents parallel incidences, none individually a multigraph bridge.

The [Fourier/square addendum](tournament_cut_fourier_squares_20261020/FOURIER_MODES_AND_COMMUTING_SQUARES.md) supplies the precise small-order conventions and executable certificates.

## 4. Balanced deletion has a literal path-level explanation

In the two-return tournament T_(a,b), the complete Hamiltonian-path decoder is a product of two pointed Boolean choice sets. Deleting arbitrary permitted interior vertices is restriction of those choices. Finite coordinates have powers-of-two fibres; exceptional infinity coordinates have one lift. This proves deletion is lawful for every path in this particular family, rather than assuming a false general deletion theorem.

Balanced central growth gives the inherited 4,2,1 fibres. The extension proves arbitrary rectangular thinning, all-depth correction terms, and the exact propagation of the full endpoint-count matrix. Blue counts require retaining the blue endpoint submatrix separately: uncoloured endpoints do not determine whether an individual path is blue.

A sharper boundary invariant is new to this packet: for every k≥2, the full (2k+1)-square endpoint matrix E_k has determinant and permanent both equal to one. Its support has a unique perfect matching, whose row-to-column permutation is a (2k+1)-cycle and whose selected entries all equal one. Thus a large rank-one interior block coexists with a full-rank integral matrix; the small boundary entries cannot be dropped.

A five-vertex counterexample proves even homogeneous transitive twins do not make arbitrary Hamiltonian-path deletion safe. The exceptional family works because its actual decoder survives restriction.

See [balanced fibres and endpoint transport](balanced_fibre_endpoint_extensions_20261020/BALANCED_DELETION_EXTENSIONS.md).

## 5. The arithmetic bridge preserves fibres, not numerical resemblance

With F counting proper nontrivial divisors, S their squarefree members, and U their prime members, the inherited classification is exactly N=p, p³, or p²qr, with distinct primes p,q,r and N≥2. U is a subset of S, so F=S+U is not a partition into the original two subsets.

For p²qr, an explicit involution pairs p with p², q with p²q, and r with p²r, and fixes the other four divisors. Its seven quotient fibres have three doubletons: 10=7+3.

On the graph side, first quotient the twenty three-edge K4 graphs by ordinary graph complement. On those ten pairs, tetrahedral opposite-edge replacement has four fixed points and three transpositions. A marked tetrahedron orientation gives an explicit reversible correspondence with four corners and six edges, transporting the arithmetic involution. The graph-side correspondence is A4-equivariant; odd relabellings introduce edge opposition. A transposition fixes two objects on one side and four on the other, proving that unmarked S4-equivariance is impossible.

Thus there is a genuine fibre-preserving bridge, with a necessary orientation marking. It preserves neither multiplication nor divisibility nor divisor complement, and supplies no Collatz conclusion.

See [prime balance and tetrahedral pairs](session17_prime_shapes_20261020/PRIME_BALANCE_TETRAHEDRAL_BRIDGE.md).

## 6. Exactly what 625 and 576 count

Let S=T_(3,3), a seven-vertex tournament with 25 Hamiltonian paths. The ordered join S→S has exactly 625 paths, one for each ordered pair of component paths. Excluding the standard path in each component removes a row and a column of the 25-by-25 grid:

625−25−25+1=576.

Deleting the literal diagonal instead leaves 600. These are different predicates. In raw component coordinates the blue set is the twisted diagonal {(P,r(P))}, where r is the source anti-reversal, not {(P,P)}. The literal diagonal contains only five blue and twenty black join paths, so deleting it leaves twenty blue and 580 black paths. Conjugating the second coordinate by r turns the blue set into an ordinary diagonal. The actual join has 25 blue and 600 black paths in total; its 576-path interior has 24 blue and 552 black paths.

A chosen five-state marking identifies HP(S) with F5², sending the standard path to zero and choice-complement to vector negation. The 576 pairs are precisely two-column matrices with both columns nonzero: 96 singular and 480 invertible. Independently quotienting each column by sign gives 144 ordered pairs on twelve nonzero sign classes, split as 12 diagonal, 12 antipodal, 60 adjacent, and 60 distance-two pairs.

The relation det(v,w)²=1 on these twelve classes is explicitly an icosahedron. An explicit state table identifies it with one previously proved regular-five two-complement return component. Vector doubling becomes tournament converse; the actual balanced-path anti-reversal becomes converse composed with a specified even relabelling. The induced SL2(5)/{±I} action agrees with the sixty even relabellings, not just with their cardinality.

The marking remains indispensable. The equally sign-compatible relabelling x↦x³ changes an edge to a nonedge. Consequently the determinant adjacency is not canonically determined by the path involution alone, and its transported action is not an automorphism group of the rigid seven-vertex source tournament.

See [the marked F5 and icosahedral correspondence](session17_prime_shapes_20261020/F5_PATH_PAIR_ICOSAHEDRON.md).

### 6.1 Extending the field gives surfaces, then singular gluings

The determinant graph is the classical Platonic graph over a finite field, not a newly discovered graph family. Fill its triangular faces, but no higher-dimensional cliques. For an odd prime p, this gives a connected orientable closed surface with

V=(p²−1)/2, E=p(p²−1)/4, F=p(p²−1)/6,
g=1+(p²−1)(p−6)/24.

The first cases are the tetrahedron, the icosahedron, and a genus-three surface at p=7. No component in this family has genus one, so this construction does not supply a torus or a torus-knot interpretation.

For q=p^r with r>1, a vertex link splits into q/p disjoint p-cycles. The complex is therefore not a surface. Splitting each vertex by those link components produces q(q²−1)/(p(p²−1)) copies of the prime-field surface, indexed by the left-coset space SL2(Fq)/SL2(Fp), not a quotient group. Each edge determines its unique prime-field plane, which proves that sheets share vertices but never edges.

At q=9 the complex has 40 vertices and normalizes into thirty tetrahedron boundaries, three through each original vertex. At q=25 the proved formula gives 130 icosahedra, five through each original vertex. The latter is a formula consequence, not a separate exhaustive computation. Substituting the singular complex's Euler characteristic into the surface genus formula would be invalid.

See [finite-field surfaces and normalization](session17_prime_surfaces_20261020/FINITE_FIELD_SURFACES.md), including the primary-literature attribution and exact characteristic-two and clique-complex boundaries.

### 6.2 Vertex identifications create a different kind of topology

Write m=q/p, c=q(q²−1)/(p(p²−1)), and V=(q²−1)/2. The incidence graph between the c normalized surfaces and the V original vertices has mV edges and is connected. Its cycle rank is b=(m−1)V−c+1.

Replacing each vertex identification by a contractible star, and then moving attachment points within each connected sheet, proves that the complex is homotopy equivalent to a wedge of c genus-g_p surfaces and b circles. Therefore its integral H1 is free of rank 2c g_p+b and H2 is free of rank c. This is a homotopy statement, not a homeomorphism turning singular points into manifold points.

In particular, the F9 complex has H1=Z^51 and H2=Z^30; the F25 complex has H1=Z^1119 and H2=Z^130. These integer and torsion statements come from the decomposition proof, not only from modular boundary-matrix experiments.

See the separately audited [homotopy supplement](session17_prime_surfaces_20261020/independent_audit/HOMOTOPY_SUPPLEMENT.md).

There is also an embedding distinction. For prime p, the graph's minimum orientable genus is exactly g_p: the Euler lower bound is attained by the triangular surface. For a proper extension q>p, no triangular embedding exists, because a vertex rotation would have to be a Hamilton cycle in its disconnected neighbor-link graph. Consequently its graph genus is at least 2+(q²−1)(q−6)/24. This gives at least 12 for F9 and at least 496 for F25; neither is asserted to be an attained minimum. The singular complex's formal Euler “genus” is not the graph's genus.

See the separate [graph-genus proof](session17_prime_surfaces_20261020/independent_audit/GRAPH_GENUS_SUPPLEMENT.md).

## 7. Reset-debt closure is valid, but coverage is a separate theorem

For the forward estimate x_out≤(x_in+δ)/λ, with λ>0 and δ≥0, chronological composition is

(λ1,δ1) followed by (λ2,δ2) gives (λ1λ2, δ1+λ1δ2).

Hence δ≤λ−1 survives arbitrary finite composition, with weighted additive slack. More generally a positive physical-height barrier B is certified by δ≤B(λ−1). A macro may pay for a growing internal reset, so forbidding every local increase is unnecessarily strong.

The concrete guarded macro with odd valuations (1,2,3) sends x to (27x+23)/64 on x≡11 mod128. At x=11 it gives 11→17→13→5. The first reset block, with valuations (1,2), grows from 11 to 13, although its individual reset step 17→13 decreases. The complete macro pays that block’s debt. This is verified on its entire arithmetic progression, not only at the displayed example.

Exact valuation cylinders retain their odd-terminal bit, prefix legality, inverse-image restriction, and ROOT stop. Arbitrarily deep finite cylinders do not establish an infinite trajectory from one positive integer. Repeating this macro has a unique 2-adic source 23/37 and no positive integer in the infinite intersection.

See [guarded reset composition](guarded_reset_debt_20261020/RESET_DEBT_COMPOSITION.md).

### 7.1 A finite-state feasibility theorem and a universal-coverage obstruction

For a finite directed block graph, edge u→v with λ>0 and δ≥0 admits positive state barriers b_v satisfying δ≤λ b_v−b_u exactly when every directed cycle has multiplier product at least one and every product-one cycle has zero debt on every edge. The proof is constructive via logarithms, and entirely rational for rational edge data. It handles growing individual edges when the complete cycle can pay their debt.

But in the physical-height model, a finite-state graph admitting such positive physical-height barriers, with finitely many fixed nonempty valuation words, cannot cover every unresolved odd trajectory. If it has V≥1 states and a word-length bound L≥1, every source congruent to −1 modulo 2^(VL+1) has its first VL valuations equal to one. Any V legal block edges would repeat a state and create a cycle with multiplier below one, contradicting the barrier criterion. ROOT is not reached during that growing prefix.

This is a precisely scoped impossibility theorem, not a Collatz counterexample. It identifies what a viable universal route must change: permit unbounded-length parameterized rules, unbounded state, or a resource other than these finitely many physical-height barriers. Finite banks can still prove useful subclasses and arbitrarily deep finite certificates.

See the [finite-state barrier and coverage theorem](guarded_reset_debt_20261020/FINITE_STATE_BARRIERS_ADDENDUM.md).

## 8. Open questions and next proof obligations

- Nothing here proves universal Collatz entry into a certified bank, repairs the separate native signed-row 7X inequality, or proves global convergence.
- It is not proved that {7,21} are the only missing positive odd raw tournament Hamiltonian-path counts; even counts are already excluded by the parity theorem. The earlier density obstruction outside the intercyclic subclass remains a reason not to confuse that subclass with all tournaments.
- Original and black metagraph connectivity, general bridge classification, and full wiring realizability remain distinct from the solved cut aggregate.
- Recent incoming endpoint-three work proves three-black detours for its stated endpoint layer and n≥11. It does not solve all layers.
- A canonical unmarked replacement for the necessary field/orientation markings would need new structure, not a different choice of labels.

## 9. Incoming work: retain its maps, repair its small-order domain

Concurrent work at the inspected commit `b0c9a1ee240e53779f0b22cbb41aada9de0e200e` in [the three marked carriers](https://github.com/eliottcassidy2000/math/blob/b0c9a1ee240e53779f0b22cbb41aada9de0e200e/05-knowledge/results/trijection_four_nodes_three_edges_20261020.md) (file SHA-256 `68d73fe7b0f17abce4bedeaa349a31b2e93a371be0db697ceb05d81a27209dc7`) distinguishes line, cycle-plus-solo, and star markings through different orbit sets. This is different from our cut-switching aggregate and must not be merged into it by matching counts.

At that inspected commit, its stated n≥3 domain for all three carriers needs a small-order restriction: the cycle-plus-solo construction requires n≥4 in simple tournaments. At n=3 its C_(n−1) would be a directed two-cycle, requiring both directions of one edge, which a tournament forbids. The proposed encoding also sends the two cycle positions to the same undirected pair. The line and star cases remain meaningful at n=3. We therefore do not inherit an all-three n≥3 theorem from that source. This is a domain correction, not a rejection of the distinct n≥4 marked constructions.

## 10. Verified publication ledger

Each entry below was published as a coherent commit and checked against its frozen payload hashes. The final entry corrects the historical domain error discussed in §9; that section deliberately keeps its inspected-source citation pinned. Packet checks and document checks are not a claim that an unrelated whole-repository test suite passed.

- Tree gauges, marked trijections, and cut-Cayley aggregates: [`ae23c720da4da91798625fdead911d519844df32`](https://github.com/eliottcassidy2000/math/commit/ae23c720da4da91798625fdead911d519844df32); 14 files, 104,358 bytes.
- Guarded reset-debt composition: [`1a32bd3b9dac0a5acd4d3b9527635226932e683d`](https://github.com/eliottcassidy2000/math/commit/1a32bd3b9dac0a5acd4d3b9527635226932e683d); 15 files, 80,921 bytes.
- Prime-balance fibres and marked tetrahedral/icosahedral models: [`b1ea84c714e5da417a5b4a0356ec92d12666e7cd`](https://github.com/eliottcassidy2000/math/commit/b1ea84c714e5da417a5b4a0356ec92d12666e7cd); 24 files, 173,298 bytes.
- Balanced deletion fibres and endpoint transport: [`80fe6ac6771dfe7e9d0747b11d568647219b06d3`](https://github.com/eliottcassidy2000/math/commit/80fe6ac6771dfe7e9d0747b11d568647219b06d3); 19 files, 132,911 bytes.
- The order-five distance-four six-cube model: [`ae588ff1fba93f5462bc3f35d78852ac055b6119`](https://github.com/eliottcassidy2000/math/commit/ae588ff1fba93f5462bc3f35d78852ac055b6119); 12 files, 39,558 bytes.
- Finite-state reset barriers and coverage obstruction: [`885d370fb328359cca48387722f874f5e39af248`](https://github.com/eliottcassidy2000/math/commit/885d370fb328359cca48387722f874f5e39af248); 12 files, 56,594 bytes.
- Finite-field surfaces, normalization, and homotopy: [`8f5f6ecfbf22a1b6cdf42d68f6a9de07f4f6ff4d`](https://github.com/eliottcassidy2000/math/commit/8f5f6ecfbf22a1b6cdf42d68f6a9de07f4f6ff4d); 11 files, 60,801 bytes.
- Minimum orientable graph-genus boundary: [`6d9b3fbd5371b549f12420f194bc6b5433b0143d`](https://github.com/eliottcassidy2000/math/commit/6d9b3fbd5371b549f12420f194bc6b5433b0143d); 1 file, 8,180 bytes.
- Fourier eigenvalues and commuting-square boundaries: [`77a66bbcfe9241b25eb686da9aac716830b80d56`](https://github.com/eliottcassidy2000/math/commit/77a66bbcfe9241b25eb686da9aac716830b80d56); 13 files, 51,236 bytes.
- Incoming marked-carrier domain correction: all three cases at n≥4, line/star retained at n=3: [`718c3ab5a741f5906952abca2826be3a0b3d1229`](https://github.com/eliottcassidy2000/math/commit/718c3ab5a741f5906952abca2826be3a0b3d1229); 1 file, 18,014 bytes.
