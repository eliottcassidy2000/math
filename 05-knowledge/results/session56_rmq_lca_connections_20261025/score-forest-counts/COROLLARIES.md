# Parity and weighted complement matching of labelled score fibres

Status: PROVED consequences of the classical score/forest theorem and the cycle-fibre connectivity theorem. This supplement leaves REPORT.md, check.py, and the core certificate unchanged. No tournament-isomorphism or A000568 identification is claimed.

## 1. Forest-count parity

Let G be a finite loopless graph, with its vertices and edges distinguished; the simple-graph case suffices here. Reversing all edges of an orientation sends its labelled outdegree vector r to

    iota(r)=deg_G−r.

This is an involution on the set of realizable labelled score vectors. Its only possible fixed vector is r=deg_G/2. If any vertex degree is odd, this vector is not integral, so there are no fixed vectors. If all vertex degrees are even, orient an Euler circuit in each nontrivial connected component; isolated vertices contribute zero. This realizes deg_G/2 and proves that exactly one score vector is fixed. Pairing the remaining score vectors therefore proves

    number of nonempty labelled score fibres is odd
    iff every vertex degree of G is even.

By REPORT.md's score/forest identity, equivalently

    number of spanning forests of G is odd
    iff every vertex degree of G is even.

“Spanning forests” includes isolated vertices and the empty edge subset. A disconnected graph with every vertex degree even is allowed; no connectedness hypothesis is implicit in the Euler-circuit construction. If G has no edges, its unique orientation, score vector, and forest give the expected odd count one.

For G_n=K_n−P_n with n≥3, endpoint degrees are n−2 and internal degrees are n−3. When n is odd the endpoints have odd degree; when n is even the internal vertices have odd degree. Thus the score-fibre/forest count F_n is even for every n≥3. The n=1,2 free graphs have no edges and F_n=1.

## 2. Keep complement edges when contracting cycle-connected fibres

Start with the labelled orientation state space of G. Add legal directed cycle-reversal edges; its connected components are exactly the labelled score fibres. Separately add an undirected native complement edge between D and its complete reversal D-bar. If G has at least one edge, this involution has no fixed orientation, so these native edges form a perfect matching on orientation states, with 2^(|E(G)|−1) edges.

Now contract each cycle-connected score fibre, but KEEP the native complement edges. This is a quotient multigraph construction: it does not identify complement-paired fibres. Its vertices are the realizable score vectors, and complement edges can only join r to deg_G−r.

Let

    w_r = [x^r] product_{uv in E(G)} (x_u+x_v).

Each monomial choice records the tail of every edge, so w_r is exactly the number of orientations in score fibre r. Complete reversal is a bijection between fibres r and deg_G−r, hence

    w_r=w_(deg_G−r).

For r different from deg_G−r, the corresponding two quotient vertices are joined by exactly w_r parallel complement edges. This is w_r, not 2w_r: counting each orientation on one chosen side counts every complement pair once.

When every vertex degree is even and |E(G)|>0, the one central fibre r=deg_G/2 is self-paired. Its w_r states are paired without fixed orientations, so w_r is even and it contributes exactly w_r/2 internal complement pairs, which become loops after contraction. All other fibres are paired as above. If |E(G)|=0, there is instead one fixed orientation: treating complement adjacency as a loop records that fixed point, while a convention requiring two distinct states gives no native edge. The formula 2^(|E|−1) and the w_r/2 pairing formula do not apply to this zero-edge case.

## 3. The fixed-path specialization

Keep the directed path 0→1→...→n−1 locked, and reverse only the residual edges of G_n=K_n−P_n. In this path-gauged labelled state space, directed residual triangles and quadrilaterals connect exactly the labelled residual-score fibres, by the separate cycle-move theorem. Full tournament scores are residual scores plus the same locked contribution, so they induce the same partition.

For n≥3, there is no fixed score vector. After contracting each cycle-connected fibre and retaining complement edges, the underlying simple skeleton is therefore a perfect matching with

    F_n/2 coarse edges.

Its parallel-edge multiplicity on the pair {r,deg_Gn−r} is w_r. The total number of retained native complement edges, before or after contraction with multiplicity, is

    2^(|E(G_n)|−1),  where |E(G_n)|=(n−1)(n−2)/2.

Thus F_n/2 counts pairs of score fibres, NOT the individual native complement lines between orientation states. In particular, for n=3,...,7 the coarse-pair counts are 1,4,26,239,2853, whereas the native-edge totals are 1,4,32,512,16384. The two happen to agree at n=3,4 because these free graphs are forests and all their score fibres are singletons.

This is NOT the original metagraph obtained by merging tournament-isomorphism classes. Labelled score vectors need not descend to those classes; colours, class multiplicities, and loops in that different quotient are not classified here. In particular the absence of loops in this path-gauged score skeleton does not rule out loops in the original isomorphism-merged metagraph.

## 4. Terminology guard

Throughout the parity theorem, the condition is explicitly “every vertex degree is even,” sometimes called the Eulerian degree condition. Do not confuse it with the different notion of “even graph” in Royle et al., “Tournaments and Even Graphs are Equinumerous,” arXiv:2204.01947v2. That paper's condition concerns the parity of edge-orientation reversals induced by automorphisms, not parity of vertex degrees. Its equinumerosity result cannot be substituted into this forest-parity theorem.

Primary source: https://arxiv.org/abs/2204.01947v2

For a small warning example, there are three unlabelled simple graphs on four vertices with every vertex degree even (empty graph, triangle plus an isolated vertex, and 4-cycle), but four unlabelled tournaments on four vertices. This does not contradict either statement: their graph predicates differ.

## 5. Verification boundaries

The independent audit in [the corollary audit](../audit-forest-mirror/corollary/REPORT.md) (check.py, REPORT.md, and byte-identical passing normal.json/optimized.json) verified all 1,099 labelled simple graphs on n=1,...,5 vertices and the path-complement cases n=3,...,7. The supplied corollary_check.py is a separate reproducible standard-library checker for those finite claims: it compares direct forest parity against the degree condition, enumerates score fibres, and verifies complement pairing, multiplicity, central-fibre parity, and path totals. Its checks use explicit exceptions, not removable Python assertions.

The all-size claims above follow from the proofs, not from these finite tests. This supplement was prepared without rerunning the frozen core checker or expanding the research scope.
