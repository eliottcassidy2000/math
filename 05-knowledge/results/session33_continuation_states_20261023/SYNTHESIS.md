# Session33 synthesis: compositional information, lawful operations, and realizability

Started 2026-10-23 02:12 UTC. This packet distinguishes proved statements, inherited results, exact finite certificates, sampled negatives and open problems. It makes no claim that all derived statements are new to the literature.

## The answer to the underlying-structure question

The strongest common explanation is **compositional sufficiency**: retain exactly the information needed by the operations you still intend to perform. Then ask a second, separate question: can that compressed state actually be realized by nonnegative counts, legitimate source owners, graph incidences, integer clocks or geometric objects?

There is no evidence here for a universal number, entropy or scalar invariant underneath all the topics. Several precise theorems instead explain both successes and failures:

1. A scalar correct now can fail after a small future operation.
2. The minimal sufficient state depends on which future operations are permitted.
3. Linear rank measures signed linear information, not automatically the number of positive or physical states.
4. A quotient may preserve every local adjacency while losing global closure or geometry.
5. Same marginal counts do not certify simultaneous compatibility.

These are not merely research slogans in this session: each is instantiated by an exact theorem, native example or explicit hostile below.

## 1. A sharp experimental budget for tournament factors

Let F_S(c)=c! times the number of spanning covers of a tournament S by c unordered directed paths. When S is a uniform module inside a larger tournament, a Hamiltonian path can enter it repeatedly. Every exit-return requires at least two exterior vertices, because a single exterior vertex cannot both beat and lose to the whole module.

Consequently all contexts with at most m exterior vertices observe exactly

F_S(1),...,F_S(floor(m/2)+1).

The reverse direction is constructive: cyclic three-block probes yield a triangular measurement system with nonzero diagonal (j!)². Thus this is an iff continuation theorem, not just an upper bound on what a particular method detects. Odd budgets2j+1 buy nothing beyond2j.

A concrete first separator is smallest at order5: two factors both have H=15 but F(2)=78 versus90. Two exterior vertices separate them into123 versus135 Hamiltonian paths; no one-exterior-vertex context can do so. At a fixed known order, inherited parity identities reduce the number of coordinates needed; we do not claim every formal coordinate is independently realizable.

The top response coefficient at budget2r is r!per(B), where B records the directions between outgoing and incoming exterior shores. Thus the deepest accessible response exists iff Hall's matching inequalities hold. Two explicit seven-vertex quotient tournaments have the same labelled scores and triangle census, but their top coefficients are0 and6. The missing information is a two-row common-neighbor deficit, invisible in the margins.

See [main theorem](tournament/BOUNDARY_OBSERVABILITY.md) and [triangle/Hall addendum](tournament/ADDENDUM.md). The deterministic20,000-draw order7 probe found no deeper response twins; this remains a sampled negative, not a nonexistence theorem.

## 2. The Fano numbers become an exact operation model

Recovered older work identifies the128 line orientations of a labelled Fano plane with a binary seven-cube. The16 zero-syndrome orientations have H=189; the112 others have H=171. These base facts are inherited and independently rechecked, not newly discovered here.

The new continuation analysis makes their roles precise:

- Current H has only two values.
- Specified line reversals require all eight syndrome states. Two nonzero syndromes with equal H are separated by one specified reversal.
- Uniform random line reversals admit a two-state Markov quotient instead. The consumer changed, so the correct state changed.
- The syndrome map Q7→K8 is locally bijective, but quotient triangles lift to open paths with a nonzero codeword displacement.

Character orthogonality shows that H has exactly seven nonconstant Walsh terms, all degree4. Every conditional average fixing at most three individual line choices is the same693/4; some four-choice averages differ. Its independent-bias polynomial has zero gradient, Hessian and third derivative at the unbiased point, but a nonzero quartic saddle. This is a concrete limitation of low-order diagnostics, not an identification with the separate rank88 Gaussian problem.

The larger operation family matters immediately. Allowing every directed-triangle reversal reaches the third regular order7 class, H=175, outside the fixed Fano chart. The complete2640-labelled-tournament audit gives class-neighbor matrix

[[0,14,0],[2,9,3],[0,7,7]]

in the order189,171,175. Its stationary masses are(1,7,3)/11 and mean Paley-class return is11, versus8 in the restricted Fano walk. Those different numbers come from different exact processes, not contradictory numerology.

See [Fano spectrum](FANO_CONTINUATION_SPECTRUM.md) and [ambient reversal walk](REGULAR7_TRIANGLE_CONTINUATIONS.md). Neither operation is the fixed-path chord complement of the merged metagraph.

## 3. A genuine Hadwiger–Nelson guardrail from the same cube

The128 cube vertices have an exact rational induced unit-distance realization in the plane: sum subsets of the seven unit vectors parametrized by slopes1,...,7. All8128 pair distances were checked without rounding. This graph is two-colorable, yet the syndrome quotient is K8.

More generally, the inherited generic-angle product construction embeds every finite cube as an induced plane unit-distance graph. The Hamming syndrome produces covers Q_(2^r−1)→K_(2^r). Therefore a plane unit-distance graph can have a locally bijective quotient with arbitrarily high chromatic number.

This proves that even local-star-exact quotienting is not a valid coloring lower-bound transfer. A homomorphism pulls quotient colorings back and gives the opposite inequality. The quotient's odd cycles fail to lift closed; a lower-bound-preserving section is absent. The same missing deck-coordinate that mattered in the Fano operation model now explains the geometric failure.

See [exact geometry and obstruction](PLANE_UNIT_DISTANCE_QUOTIENT_OBSTRUCTION.md). This is not a new numerical bound for the chromatic number of the plane.

## 4. Direct progress on the actual tournament metagraph

For the actual fixed-path complement tau, which satisfies tau²=id, and merged node output o, the minimal deterministic continuation state is simply

(o(x),o(tau x)).

It remembers the current node and its actual complement partner. If color is an output, retain color too: order6 has explicit same-pair blue and black presentations. This is a native metagraph statement, not a substituted auxiliary graph.

The minimal signed-linear observer dimension is2r-c-b, where r is the number of node fibers, c the number of support components and b the number of bipartite components. This uses the inherited defect-rank theorem. The new exact memory calculation determines every projected return kernel from the normalized incidence matrix K. For fiber-uniform initial impulses these response matrices alternate I,K,I,K,..., while K^t describes the different experiment of re-randomizing the tiling inside its current fiber after each step. Arbitrary nonuniform initial sources additionally require the initial sidecar eta=y1−Ky0; K alone does not determine that forcing.

At order5 the actual figures are:

-64 raw tilings;
-10 merged nodes;
-38 deterministic/positive continuation states;
-19 signed-linear observer dimensions.

The positive minimum assumes every pure tiling can be initialized and exact node outputs are required. Positivity forces distinct current/next pairs into disjoint hidden-state supports. A19-dimensional signed realization is therefore not a19-state stochastic tiling model.

An additional operation boundary is sharp at every order: if every labelled free-arc flip is allowed, no nontrivial continuation quotient remains. The transitive node has a unique fixed-path presentation; translate any chosen tiling to that anchor and any distinct tiling becomes distinguishable. The signed-linear closure also contains all translated delta functions. This is a state-model obstruction, NOT an algorithmic lower bound against faster static A000568 enumeration.

See [native involution memory](metagraph-memory/NOTE.md) and [full-flip boundary](metagraph-memory/ADDENDUM_FULL_FLIPS.md).

## 5. Exact compression and exact refusal to compress

The square-sum graph at25 supplies both sides. Against its specified surrounding graph, its three blocks'5/101/149 endpoint-labelled states collapse to3/5/5 response classes including dead classes, with rational ranks2/4/4. Simultaneous weighted contraction still gives all ten Hamiltonian paths. Endpoint-polynomial weights retain their full endpoint spectrum.

Against unrestricted continuations, distinct perfect boundary pairings cannot be merged. Small matching matrices also exhibit different ranks over Q and F2. A parity representation can erase a positive count; it is not automatically a counting representation.

Electrical elimination is an especially useful hostile: a four-vertex path and star have the same two-terminal Dirichlet response and the same eliminated determinant, but Hamiltonian path counts1 and0. Their Gaussian/quadratic boundary statistics are genuinely equal; the omitted vertex-degree and connectivity obligations are not. No scalar resistance correction repairs the lost property.

See [port connection algebra](ports/PORT_CONNECTION_ALGEBRA.md).

## 6. Credit coverage: an exact realization test, not invented extra budget

At each actual arrival cut, an independent source bank is described by total capacities for each nonempty exact eligibility subset of obligations. This mask profile is necessary and sufficient for equality of responses to all nonnegative real demand vectors, or to all nonnegative integral demand vectors when capacities are integral. Empty-mask capacity is unobservable. Boolean inversion recovers the profile from the Hall rank function; a restricted demand alphabet need not determine it.

Nonnegative recovered coefficients characterize realizable independent banks. The perfectly valid shared-gate rank r(S)=min(|S|,2) on three labels has a coefficient−1, proving that it cannot be represented as independent donor masks. A missing shared bottleneck cannot be repaired merely by renaming credits.

This gives a concrete interface test for future source-preserving certificates. It does not improve the unrestricted scalar annulus estimate: when every credit can pay every obligation, the result reduces to the existing max(p−b,0). Restrictions on legal use may make a certificate more meaningful while weakening its numeric bound. The residual R', genuine actual-clock supply and uniform153N coverage remain open.

See [clock and credit proof](clocks/PROOF.md).

## What these connections really buy

They replace vague parallels by specified maps and decisive tests. The most useful new terminology here is operational:

- A **continuation budget** is the permitted size or length of a test context.
- **Response twins** agree on all contexts within that budget.
- A **memory defect** measures the extra information needed after a chosen operation.
- A **realizability obstruction** says a formal response cannot come from the demanded positive, integral, geometric or source-owned carrier.

These are defined tools, not a claim of a new universal foundation of mathematics. Their value is that they tell us which additional invariant to seek and which apparently powerful shortcut cannot possibly work.

All specific open targets and stopping conditions are in [NEXT_PROOF_TARGETS.md](NEXT_PROOF_TARGETS.md); the [concept board](CONCEPT_BOARD.md) records the tested cross-topic maps. Independent reports accompany the packet. The separate Session32 attachment-derived publication remains excluded.
