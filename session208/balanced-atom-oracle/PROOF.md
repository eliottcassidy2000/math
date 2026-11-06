# Balanced atoms and an exact arc-flip oracle

Status: PROVED local counting identity and construction bound; VERIFIED explicit witness. Search is bounded and nonexhaustive. No global minimum or optimal asymptotic rate is asserted.

## 1. A balanced strong atom

There is a strongly connected tournament on13 vertices with exactly759375=15^5 directed Hamiltonian vertex-word paths. Its adjacency code is241005792788926351623915: bit k, starting at zero, belongs to the kth lexicographic pair i<j, and equals one exactly when i->j. The code is stored as an arbitrary-precision decimal integer.

Independent reconstruction, all-source reachability and a suffix subset recurrence verify this count and strong connectivity. Dominance joins of q copies give count15^(5q) on13q vertices. For k=5q+s,0<=s<5, add s order-five H15 atoms to obtain

 m(15^k)<=13q+5s,   lim_k m(15^k)/k<=13/5.

This improves the preceding8/3 construction. It does not claim m(15^5)=13 or optimality of the asymptotic rate. The subadditivity limit and dominance-join rule are inherited from the preceding audited packets.

## 2. Exact local oracle from the missing endpoint data

Let T be a tournament on V. For every nonempty A subset V and u in A, let E(A,u) count directed Hamiltonian paths of T[A] ending at u; let B(A,u) count those beginning at u. The usual endpoint subset recurrences compute both tables. For distinct u,v define

 D_uv = sum_{A subset V, u in A, v notin A} E(A,u) B(V minus A,v).

This expression is defined regardless of the orientation of the pair u,v. It counts vertex words whose prefix on A is a directed path ending at u, whose suffix on V minus A is a directed path beginning at v, and whose adjacency u,v is left unconstrained. In particular, if u->v, D_uv counts exactly the Hamiltonian paths of T using that arc consecutively.

Theorem. If u->v in T and T' reverses only that arc, then

 H(T')-H(T) = D_vu-D_uv.

Proof. Paths that do not use the two vertices consecutively are unchanged by the reversal. The paths destroyed are exactly the D_uv words. The paths created are exactly the D_vu words. The endpoint factors in the latter expression can be computed in T rather than T': u and v belong to different subsets in every summand, so neither induced subgraph contains the reversed edge. Cutting at the consecutive pair gives a bijection in each case, with no multiplicity: the prefix vertex set is uniquely determined. QED.

Both endpoint tables and all binom(n,2) differences can be computed using O(n^2 2^n) integer arithmetic operations and O(n 2^n) stored integers. Once computed, every one-arc neighbor's count is known exactly. Computing each neighbor independently by subset DP would instead require O(n^4 2^n) operations. Integer bit lengths are O(n log n), since all counts are bounded by n!; these are arithmetic-operation bounds, not unit-cost bit-complexity claims.

A further check is

 sum_{u->v in T} D_uv = (n-1) H(T),

because each Hamiltonian path contains exactly n-1 consecutive arcs. This checks the normalization and is not an estimate.

## 3. What the search does, and what it does not prove

The source `target_search.cpp` uses the exact difference table to select an arc reversal bringing H closer to a specified target. At a local minimum it makes a pseudorandom move or restarts. Every chosen difference is checked against a fresh exact DP count after the move. Seed20813 and a budget of300 restarts found the order13 witness at zero-based restart70, step3. The witness proof uses its exact verified code, not a convergence assertion about this search.

The orientation-hypercube search reverses ONE arc. It is not the repository's native path-complement metagraph move, which retains one Hamiltonian path and reverses all other arcs. No degree, fibre, or native-metagraph claim is inferred from this search graph.

The structural bridge is that Session207's induced-subset endpoint data, which restores information lost by the scalar H, also provides an exact local improvement oracle. The attached paper's mixture idea enters after discovery: atom exponent vectors can be combined by rational nonnegative weights, then denominators cleared by dominance joins. This balanced atom already lies on the target exponent ray, so a one-type mixture suffices. Nothing here uses or certifies the paper's main logarithmic spanning-tree theorem.

Failure to hit a target in a finite search is not a spectrum exclusion. A rate13/5 upper bound does not rule out larger efficient atoms, a lower global rate, or a zero limiting linear rate.
