# Sparse graph carriers that preserve the entire tournament automorphism group

**Status: PROVED elementary constructions and asymptotic bounds; FINITE-EXACT verification.** No novelty claim. Graphs are finite simple undirected graphs; tournaments are nonempty. This strengthens the class-faithful padded codec in `GRAPH_CARRIERS.md`: class reconstruction and exact automorphism preservation are different requirements.

## 1. A rigid orientation codec

For each tournament vertex v make core vertices c_v,t_v,h_v, edges c_vt_v,c_vh_v, and pendant paths with the following lengths:

- at c_v: one path of length 1;
- at t_v: paths of lengths 1 and 2;
- at h_v: paths of lengths 1, 2 and 3.

Add the undirected edge t_u h_v exactly when u→v. Call this graph B(T).

It has 13n vertices and 12n+binom(n,2) edges. It is connected. Its core vertices are exactly the vertices of degree at least three. Every other vertex belongs to a pendant path. Counting these paths intrinsically identifies c,t,h by the signatures (1),(1,2),(1,2,3); each c has precisely its matching t and h among its core neighbors. Thus T is reconstructed without a chosen vertex order. Every graph isomorphism induces a tournament isomorphism; conversely every tournament isomorphism extends uniquely, because paths of distinct lengths have no interchangeable arms and each path fixes its attachment end. In particular,

    Aut(B(T)) ≅ Aut(T),

with the isomorphism given by the actual action on centers, not just equal group orders.

This removes the independent leaf permutations in the earlier 9n-vertex graph. The retained port matching is essential; independently quotienting two unaligned graph layers does not provide this property.

## 2. Balance the global edge deficit without adding symmetry

Let

    k = binom(n,2) − n + 1 = (n−1)(n−2)/2.

For n=1,2, k=0 and B(T) already has E=V−1. For n≥3, B(T) has cycle rank k and contains a cycle. Append k pairwise nonisomorphic rigid trees. Their components cannot mix with B(T), because they are acyclic. Each contributes one to V−E, and no nontrivial automorphism. Therefore the resulting graph G(T) satisfies

    |E(G(T))| = |V(G(T))|−1,
    G(T)≅G(U) iff T≅U,
    Aut(G(T)) ≅ Aut(T).

When k>0 the unique cyclic connected component recovers B(T); for k=0 there is only that component. Hence decoding is intrinsic. The choice of padding depends only on n, not a labelled presentation of T. No tournament canonicalization oracle is used.

### Explicit logarithmic-size rigid padding bank

For k≥1 choose L=max(1,ceil(log2 k)) and k distinct binary words of length L. For each word construct a tree with a root and four branches. Three are simple pendant paths of lengths 1,2,3. The fourth is a backbone with L bit vertices followed by four fixed terminal vertices. Attach one pendant leaf at bit position i exactly when its bit is 1.

The root is the unique degree-four vertex. Its fourth branch has height at least L+4, so it is distinguishable from the three short arms. Along the backbone every continuing suffix has height at least four, so it cannot swap with an optional one-edge leaf. Consequently each tree is rigid and its length-L word is intrinsically recoverable. Distinct words of this fixed length give nonisomorphic trees.

Each tree has 11+L+popcount(word) vertices, at most 11+2L. Therefore

    |V(G(T))| ≤ 13n + k(11+2L) = O(n² log n).

The carrier is also sign-even in the automorphism-sign sense of `GRAPH_CARRIERS.md`: the edge-reversal sign is a homomorphism from its odd-order automorphism group to {+1,−1}, hence trivial. This is not the assertion that every vertex degree is even.

This is one graph, not a necessary pair or triple. It is not claimed to be a minimum-size general graph encoding.

## 3. Why symmetry-preserving forest padding costs more than isolates

The simpler class-faithful construction appends k identical isolated vertices: that adds k vertices but also a factor S_k to its automorphism group. That is harmless to reconstruction of isomorphism classes, but not to exact symmetry.

**Lemma.** A finite forest has odd-order automorphism group if and only if all its tree components are rigid and pairwise nonisomorphic. In that case its automorphism group is trivial.

Proof. Equal components admit an order-two swap. A nonrigid finite tree has an automorphism of order two: its center is one vertex or one edge; a center-edge exchange gives such an involution, while in the center-fixed case rooted automorphisms arise by permuting isomorphic branches, or recursively from a nonrigid branch. Any nontrivial branch permutation supplies a transposition. Induction proves the claim. Conversely distinct rigid components leave no choices.

Every finite tournament automorphism group has odd order: an even-order group has an involution, and a nonidentity involution on vertices swaps a pair, reversing the arc between them. Thus any forest padding preserving the tournament's automorphism group must use pairwise distinct rigid trees.

There are at most 4^r unlabelled tree types on r vertices, since every unlabelled tree admits a rooted plane presentation, and rooted plane trees are counted by the Catalan number C_(r−1)<4^r. Hence the number of tree types of order at most r is at most (4^(r+1)−4)/3. Given k distinct trees, choose r of order (log_4 k)/2. For sufficiently large k fewer than k/2 can have at most r vertices, so the total number of vertices is Ω(k log k). The explicit bank above gives O(k log k).

Therefore the padding overhead Θ(k log k) is optimal **within the operation of appending a forest to this fixed cyclic core while preserving odd-order symmetry**. This is not a lower bound for every possible faithful encoding of a tournament, nor for arbitrary changes to the core, nor for class-only reconstruction.

## 4. Three exact fidelity levels

1. A global equation E=V−1 does not force a forest when disconnectedness is allowed: cycle rank β equals number of components minus one.
2. Dense-core plus identical isolates can retain the source's isomorphism class with O(n²) total vertices, while adding symmetry.
3. Dense rigid core plus distinct rigid trees retains its actual automorphism group, with O(n² log n) total vertices and optimal forest-padding overhead for that core.

This is a concrete instance of the session's loss ledger: one must name whether the target is a scalar statistic, an isomorphism class, or the full symmetry action. A representation can solve one task and fail another.

## Reproduction

Run `python3 check_sparse_symmetry.py` and `python3 -O check_sparse_symmetry.py`. The standard-library checker directly counts graph automorphisms by color refinement followed by adjacency-preserving backtracking; it does not assume the theorem's group formula. It checks all 75 labelled tournaments through order four and 24 deterministic sampled order-five tournaments, before and after vertex relabelling; it compares graph automorphism orders to direct tournament permutation counts. It also checks all 126 bit trees of lengths 1 through 6 and the hostile double-copy forest, whose automorphism order is 2.

An independent auditor additionally used recursive tree canonization on 1,023 bit trees of lengths 0 through 9, including the optional L=0 four-arm tree. The public construction fixes L≥1 and only requires distinctness at a single chosen length.
