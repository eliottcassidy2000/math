# A prime even-order obstruction to the endpoint-parity criterion

The assertion that every even-order tournament has a vertex with an odd sum of its Hamilton-start and Hamilton-end counts is false. This is an obstruction to that sufficient criterion, not a full named-response collision.

## 1. Exact eight-name counterexample

Use physical names0 through7 and the out-neighbor bit rows

    Q=(238,4,240,118,195,146,34,74).

For example bit j in row i means i->j. All diagonal bits vanish and every unordered pair has exactly one direction. Exact enumeration of all8!=40320 vertex permutations gives

    S=(87,13,13,42,31,19,18,26),
    E=(13,65,29,22,13,29,52,26),
    H(Q)=249.

Every coordinate S_i+E_i is even: the start and end parity masks are both55. The verifier independently reproduces these vectors by subset dynamic programming and literal permutation enumeration.

Every proper vertex subset of size at least two is checked against every exterior vertex. None is a module, so Q is prime. In particular this example is not caused by an ordered join or a proper homogeneous block.

A complete endpoint-parity census of all labeled tournaments at orders2,4,6 has no balanced example. Those finite controls locate order8 as the smallest positive even order at which this criterion can fail; they are not used to infer any larger-order behavior.

## 2. What this does and does not obstruct

The all-order endpoint transfer in `../parity_profile_transfer/PROOF.md` uses a nonzero bit S_i+E_i as a sufficient separator for arbitrary child substitutions. This Q has no such anchor. Moreover that packet's full modulo-four law shows that when every child is nontrivial, EVERY child-uniform quotient-converse query for Q[A_i] is zero modulo4. This is an unbounded statement of coarse modular blindness only.

Neither exact coarse equality nor complete named-response equality follows. The verifier retains the complete named response bank for the displayed eight-name Q and its converse, including every bidegree-at-most-(2,2) coordinate, zero coefficient, cubic and quartic term, and full-mask checks. It explicitly finds a nonzero physical response difference.

It also checks an exact shared-cubic separator for the16-name pair inflation Q[TT2,...,TT2], with pair interiors fixed. The inherited eight-pair quotient theorem supplies a nonzero coefficient modulo8 when its shared-cubic matrix is nonzero. Here that matrix is recomputed directly from even-subset endpoint parities, and the selected coefficient is recovered by the signed sum of eight actual physical queries on each side. Thus the full response is visibly different even though the entire coarse modulo-four endpoint test is blind.

## 3. Exact controls

The census uses integer bit-vector recurrences, not random sampling. The counterexample is separately checked by all permutations. The physical Hamilton DP uses unsigned64-bit integers only for augmented orders at most18;18! is below its range. Every physical bank value is an entire Hamilton count, not a truncated series. Source and optimized modes must reproduce exactly.

The result rules out one proposed all-order route to even-quotient separation. It does not weaken the previously proved eight-pair separation theorem, assert a response collision, or classify all endpoint-balanced quotients.
