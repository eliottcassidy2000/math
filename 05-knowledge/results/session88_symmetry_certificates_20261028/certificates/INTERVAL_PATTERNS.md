# Interval-column recovery from maximum spanning trees

**Status: CLASSICAL/CITED for Fulkerson–Gross (1965), Theorem 2.1, and the Hunter union bound; PROVED self-contained proofs and algorithmic specializations; FINITE-EXACT controls in `check_certificates.py`.** In particular, the arbitrary-set Gram rigidity in §6 is a classical theorem, independently rederived here before its exact primary-source identification—not a new theorem. No claim of new maximum-spanning-tree theory or literature priority. The concrete repository improvement is a quadratic-arithmetic specified/partial-pattern query and output-sensitive polynomial recovery of the **entire column multiset**, in place of THM-2789's displayed exponential Boolean inversion. It does not reconstruct gap order or prove graceful extensions exist.

## 1. The old source and new question

**THM-2789**, `01-canon/theorems/THM-2789-interval-gram-tomography-and-graceful-gap-tail-quadratic-detector.md`, considers intervals `I_1,...,I_m` on `X={1,...,N}` and their Gram matrix

    G_ij=|I_i intersect I_j|.

The quantitative interval Helly identity is

    |intersection_{i in A} I_i| = min_{i,j in A} G_ij       (A nonempty).

Together with `h(empty)=N`, Boolean Möbius inversion recovers the count of every exact incidence column. The source explicitly warns that the quadratic size of `G` does **not** make its displayed exponential reconstruction a quadratic-time algorithm.

THM-2789 already gives exact graceful suffix-column counts, not just a generic information-theoretic sufficiency assertion. Its §3 formula (13) is a Boolean sum over supersets of the suffix; it explicitly calls the data quadratic-*size*. The improvement here is computational, not discovery that Gram determines those counts. Often the actual question asks for only one pattern `A`; there is no need to rebuild the entire column distribution first.

## 2. Union identity for intervals

For intervals `K_1,...,K_k` on a finite line, possibly empty, form the complete graph on `[k]` with nonnegative weights

    w_ij=|K_i intersect K_j|.

Let `W` be its maximum spanning-tree weight, with `W=0` for `k<=1`. Zero-weight edges may connect disjoint components. Then

    |union_i K_i| = sum_i |K_i| - W.                      (1)

**Universal inequality.** For any finite sets and any tree `T`, at a point belonging to `r>0` sets, the induced subgraph of `T` has at most `r-1` edges. Summing over points gives

    sum_{ij in T}|K_i intersect K_j|
       <= sum_i|K_i| - |union_i K_i|.

Hence every tree gives an upper bound on the union. This is the finite counting form of the classical spanning-tree/Hunter union bound.

**Equality for intervals.** Ignore empty intervals initially, and order nonempty intervals by increasing left endpoint (ties arbitrary). For each interval after the first, choose as its parent a previous interval with maximum right endpoint. The portion of the current interval already covered by the preceding union is exactly its intersection with this parent. This is true even across a gap: both are empty. Adding the current interval therefore increases union size by its length minus its parent overlap. Parent indices decrease, so these edges form a tree; attach empty intervals arbitrarily with zero-weight edges. Summing the increments proves equality for this tree. The universal inequality makes it maximum. ∎

Equivalently there exists a tree on the intervals such that the indices of intervals containing each point induce a connected subtree. This is an incidence-structure statement, not a claim that arbitrary pairwise-intersecting sets satisfy (1).

## 3. Exact-pattern formula using only `G,N`

For a desired index set `A subset [m]`, write

    n(A)=#{x in X : {i:x in I_i}=A}.

Let `B=[m]\A` and `J=intersection_{a in A}I_a`, taking `J=X` if `A` is empty. Define

    L = min_{a,b in A} G_ab,           A nonempty;
    L = N,                            A empty.

For each `i in B`, define

    k_i = min(L, min_{a in A}G_ai),    A nonempty;
    k_i = G_ii,                       A empty.

Then `k_i=|J intersect I_i|`. For `i,j in B`, define

    w_ij = min(k_i,k_j,G_ij).

The quantitative Helly identity gives `w_ij=|J intersect I_i intersect I_j|`. Thus these are precisely the pairwise overlaps of the intervals `K_i=J intersect I_i`. Applying (1),

    n(A) = L - sum_{i in B} k_i + MaxST(w).               (2)

All boundary conventions are explicit:

* `A=empty` uses ambient length `N`, which Gram alone does not supply;
* `A=[m]` has `B=empty`, so the answer is `L`;
* empty `J` gives `L=k_i=w_ij=0`, hence zero;
* repeated or empty intervals are allowed;
* `m=0` gives `n(empty)=N`;
* zero-weight edges remain available in the complete graph.

Computing `L,k,w` and dense Prim uses `O(m^2)` comparisons/arithmetic operations and `O(m^2)` stored input; a lazy weight accessor can use `O(m)` extra working memory. The code materializes `w` for clarity. Counts have `O(log(N+1))` bits, and bit-complexity must include that arithmetic. This is genuinely polynomial in the explicit dense Gram input size, **not** in the number of all `2^m` patterns output by a full reconstruction. Querying all suffix patterns separately is at most `O(m^3)` arithmetic operations; no optimized all-suffix bound is claimed.

## 4. Partial patterns and polynomial full-multiset recovery

The same formula works when `A` is the set of **required** rows and `B` is any disjoint set of **forbidden** rows; rows outside `A union B` are unconstrained. Keep the same `J,L,k_i,w_ij`, with indices only in this `B`. It counts exactly

    #{x : x in every I_a (a in A), x in no I_b (b in B)}.

The proof is unchanged: this is `|J \ union_{b in B}(J intersect I_b)|`. Thus every partial binary pattern has an `O(m^2)` oracle.

### At most `2m+1` realized patterns

Along the integer line each nonempty interval changes membership at its left endpoint and just after its right endpoint. There are at most `2m` change positions, hence at most `2m+1` constant runs and at most

    s <= min(N,2m+1)

**distinct** realized column patterns. Repeated patterns on separated runs are counted together, and the empty pattern includes all uncovered positions. For `N=0`, `s=0`; empty/repeated intervals cannot increase the bound. The discussion is on the declared integer line: it does not silently switch to real-measure conventions with exceptional endpoint atoms.

### Positive-prefix trie

Fix the row order `1,...,m`. Start with the unconstrained pattern of count `N`. At a positive prefix whose first `d` memberships are fixed, query its child requiring row `d+1`; the other child's count is the parent's count minus this answer. Recurse only into children with positive count. At depth `m`, emit the named pattern and its multiplicity.

Every realized full pattern has one prefix per depth, so at each depth there are at most `s` positive prefixes. There are at most `m s` internal positive nodes and hence at most `m s` oracle calls. Therefore the entire incidence-column multiset is recovered using

    O(s m^3) arithmetic operations,  hence O(m^4),

with `O(s m)` output bits for pattern names plus `s` multiplicities of `O(log(N+1))` bits. This is a deliberately simple bound, not claimed optimal. The implementation can keep a depth-first stack and emit leaves; no list of all `2^m` candidate patterns is ever built. For the degenerate empty family `m=0`, the output is just the empty pattern of count `N` if `N>0`.

This strengthens the computational recovery result, not its information content: THM-2789 already proved that the multiset is determined. The new procedure still cannot order those columns along the line. Enumerating all unobserved patterns explicitly would of course still require exponential output.

## 5. A tree is also a one-sided certificate

For any supplied spanning tree `T` on `B`, put

    c_T = L - sum_i k_i + sum_{ij in T} w_ij.

The universal inequality proves `c_T<=n(A)`. Therefore `c_T>0` is a sufficient existence witness: a verifier need only check the tree, its local edge weights, and the integer inequality. It need not verify that `T` is maximum. With literal intervals, a specified point is an even shorter positive witness; the tree is useful when the available input is Gram data.

For an exact count, or an absence claim `n(A)=0`, maximum-tree optimality matters. It can be established by recomputing dense Prim in `O(m^2)`, or by a tree certificate satisfying the maximum-tree cycle condition: every non-tree edge has weight at most the minimum edge weight on its tree path. This is a local optimality certificate with all-pairs checks, not permission to infer absence from a low-scoring arbitrary tree.

The phrase “short certificate” excludes the already supplied dense Gram input and the interval-realizability premise. The next lemma shows that an explicit interval realization of `G` can certify this premise even if the original sets were not supplied in interval order. Finding that realization is a separate task; PSD and elementary pairwise bounds alone are insufficient.

## 6. A short interval-realization witness certifies every set realization

**Classical Fulkerson–Gross rigidity, specialized notation.** Suppose `G` has some realization by intervals on an `N`-point line. Then every family of arbitrary subsets of any `N`-point universe with Gram matrix `G` has exactly the same full column multiset. In particular the formulas and trie above apply to that family even when its named ground-set order does not make its sets intervals.

This is exactly **Fulkerson–Gross (1965), Theorem 2.1**, [primary paper, p. 837](https://msp.org/pjm/1965/15-3/pjm-v15-n3-p11-p.pdf), checked 2026-10-28. Their incidence matrices have **points as rows and sets as columns**, the transpose of this note's convention. Their equality `A^T A=B^T B`, together with consecutive ones and equal row count, forces `B=PA` for a permutation matrix `P`. Thus their row-permutation conclusion is precisely our equality of point-column multisets. The following independently rederived zero-slack proof is included to expose the mechanism used by our certificates; it does not establish novelty.

**Self-contained proof for triples.** If a pair overlap is zero, the triple intersection is zero in every realization. Otherwise three realizing intervals pairwise overlap. Choose intervals supplying the leftmost and rightmost union endpoints; their union contains a remaining interval because those two overlap. If a single interval supplies both extremes, choose it and any other, with the same conclusion for the remaining one. Thus for one index `i`,

    I_i subset I_j union I_k,
    |I_i intersect I_j intersect I_k| = G_ij+G_ik-G_ii.

For intervals the triple size is also `min(G_ij,G_ik,G_jk)`. Any arbitrary sets with the same Gram satisfy the universal lower bound `triple >= G_ij+G_ik-G_ii` and upper bound `triple <= min pair overlap`; these bounds now coincide. Hence their triple intersection is forced to the interval value.

**Larger intersections.** For `|A|>=2`, choose distinct `a,b in A` minimizing the pair overlap over `A`. The minimum over pairs including diagonals is attained by a cross pair too, because `G_ab<=G_aa` and `G_ab<=G_bb`. Every triple `{a,b,c}` has intersection size `G_ab`, so `S_a intersect S_b` is contained in each `S_c` in the finite-set sense. Therefore the full intersection over `A` has size `G_ab`. Singletons are the diagonal values. Every higher intersection is thus forced, and Boolean Möbius inversion proves that all nonempty-pattern multiplicities are identical; the common ambient size `N` also fixes the empty pattern. This inversion is a proof of uniqueness, not the implementation. ∎

**Zero-slack mechanism.** The forcing is not a mysterious recovery of joint data. The two bounds above are nonnegative set identities:

    triple - (G_ij+G_ik-G_ii) = |S_i \ (S_j union S_k)| >= 0,
    G_ab - |S_a intersect S_b intersect S_c|
        = |(S_a intersect S_b) \ S_c| >= 0.

An interval realization shows that the available pair data make the lower and upper endpoints equal. Their nonnegative slacks must therefore vanish in **every** realization. The forbidden atom regions disappear, and a minimum pair becomes the full-intersection witness. In atom-coordinate language this is a zero-slack face certificate in a nonnegative cone. It provides the missing joint information without reconstructing the original point order. The inspected THM-2789 states recovery inside the interval category; the broader arbitrary-realization uniqueness is already Fulkerson–Gross's classical Theorem 2.1. Our contribution here is recovering that connection and applying its rigidity to short interval-realization certificates and the explicit partial-pattern/trie pipeline, not claiming a new binary-factorization uniqueness result.

An explicit certificate is merely one list of `m` interval endpoint pairs, with empties allowed. After range checks, the verifier computes each pairwise overlap from endpoints and compares it to `G_ij`, taking `O(m^2)` integer operations. Its length is `O(m log(N+1))` bits. It need not give a point-by-point bijection from the original set family to this realization, because the lemma already certifies the decision-relevant column multiset. It does **not** identify which original point carries which pattern.

The finite control inspects all `65,536` families of four arbitrary subsets of four points. There are `3,462` interval-realizable Gram matrices in that universe; `56,536` arbitrary families match one, and all have the predicted unique column multiset. These are ordered families, with no quotient or omitted empty sets. Independent audit separately verified the elementary proof; primary-source checking subsequently identified the theorem exactly as Fulkerson–Gross (1965), Theorem 2.1.

## 7. Sharp hostiles and information not restored

1. **Arbitrary sets/circular arcs.** Take `{1,2}`, `{2,3}`, `{1,3}`. Diagonal overlaps are two, off-diagonal overlaps one, but the triple intersection is empty. The interval formula incorrectly reports `n({1,2,3})=1` and `n(empty)=-1`. This Gram matrix has no interval realization. An interval-realizability promise or checked interval witness is load-bearing; PSD and pairwise entry bounds alone do not suffice.
2. **Address loss remains.** Reflecting every interval on `1,...,N` preserves `G` and every column count while moving the column addresses. Thus (2) does not recover the position of a useful graceful gap or the side of a marked vertex.
3. **Pattern names remain.** Permuting row names permutes which specified pattern is queried. An anonymous spectrum of overlap values is not enough.
4. **A recovered multiset is not all future operations.** Even after recovering every realized pattern and multiplicity, this algorithm does not recreate a faithful ordered incidence matrix for future labelled operations. THM-2793's additional ramp and root-zero data remain meaningful.

## 8. Exact transfer to the graceful gap query

For a graceful `m`-edge tree, row `d` is the interval of cuts crossed by its unique edge of difference `d`. An unmarked suffix-gap candidate is a column with exact support `A={k,...,m}`. Formula (2) counts exactly those candidates from the inherited interval Gram, without the `2^(k-1)` Möbius sum. This is the map and preserved predicate.

What is destroyed is still gap order/address and the position of the proposed attachment vertex. The restoration sidecar is the oriented ramp moment vector plus root anchor in **THM-2793**, `01-canon/theorems/THM-2793-oriented-ramp-reference-and-rooted-graceful-gap-reconstruction.md`. The cheapest decisive test is to compare (2) against literal column enumeration; the accompanying checker does that exhaustively on its declared small universe, including empties and repetitions.

## Literature scope

The union-tree inequality predates this work: David Hunter, [*An upper bound for the probability of a union*](https://www.cambridge.org/core/journals/journal-of-applied-probability/article/an-upper-bound-for-the-probability-of-a-union/092D711504BA968EF0D1D903A2685D60), Journal of Applied Probability 13 (1976), 597–603, DOI `10.2307/3212481`. Publisher metadata/abstract checked 2026-10-28; 2016 is the online republication date, not the original year. The self-contained proof above is the mathematical dependency; a focused search identified the classical bound, not novelty of the interval equality. The repository-specific result is the explicit `G -> clipped overlaps -> one MaxST -> partial-pattern count -> positive-prefix trie` implementation and its application to THM-2789. No exhaustive literature priority search is claimed.

**Classical rigidity source:** D. R. Fulkerson and O. A. Gross, [*Incidence matrices and interval graphs*](https://msp.org/pjm/1965/15-3/pjm-v15-n3-p11-p.pdf), Pacific Journal of Mathematics 15(3) (1965), Theorem 2.1, p. 837. This supplies exact priority for §6, not merely an analogy. The independent proof, finite controls, endpoint-witness verification, and MaxST/trie specialization are retained as an explicit implementation and application of classical structure. No global novelty claim is made for the algorithmic pipeline either.
