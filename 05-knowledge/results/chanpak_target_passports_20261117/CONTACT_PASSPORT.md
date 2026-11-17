# Target-first one-exit cores: a contact ledger, exact small ranks, and a growth barrier

A prescribed prime cannot in general be placed by choosing contacts on a deletion of bounded Hamiltonian count. This is an all-order restriction on the one-exit geometry, even allowing arbitrary growth in the number of vertices. It also yields the first feasible deletion counts for an exact normalized ±1 offset. It does not restrict alternative geometries or arbitrary intercepts.

## Two tests on the actual contacts

Let R be a deletion tournament, g=H(R), and pick Q⊆V(R). Add a vertex q beating precisely Q. Now attach a nonempty module B by R→B→q. Write (P,h)=(P₂(B),H(B)). The [two-run law](PLUS_RAY.md#an-arbitrary-module-plus-reader) gives H(output)=WP+Kh, where W is the weighted starting cover count excluding the singleton q, and K is the normalized constant coordinate.

They can be read directly from R without assembling the output. For every ordered two-slot cover (X,Z) of R, allow X empty but require Z nonempty. Put

    W = # such covers with start(Z) in Q,
    U = # such covers with X nonempty, end(X) in Q, and start(Z) in Q.
    K = g−U.

Indeed q can be prepended to Z in every cover counted by W. This joins to a full path X,q,Z exactly when X is empty or its last vertex beats q, i.e. is not in Q. Therefore W−U equals the number L of full paths of R+q not ending in q, and K=L+g−W. In particular K≤g. The coefficient W is linear in the contact bits and U is a nonnegative quadratic in them. The inverse at prescribed W,K is therefore a physical weighted contact problem, not only a congruence. Rédei parity applied with B a singleton shows K is odd, so U is even.

For later use, the ledger is compositional. For an ordinal join R₁→R₂ use subscripts for each block, let W_s,W_e count tagged starts/ends of one nonempty path with the other optional, and let S,E count full starts/ends at tagged vertices. Then

    g = g₁g₂,
    W_s = W_s,₁ P₂(R₂)+g₁ W_s,₂,
    U = g₂U₁+g₁U₂+E₁S₂+W_s,₁ W_e,₂.

For U, if both marked endpoints are in one block the other block is forced to one path. If end(X) is in R₁ and start(Z) in R₂, both blocks are forced to full paths in different slots. In the reverse ordering of marked endpoints the two splitting choices are independent. These different terms explain why matching only W does not preserve an exact correction.

## Trimming a bounded deletion to dyadic slopes

Take the strongly connected components of R in ordinal order. If Q is nonempty, mark the interval from the first component meeting Q to the last one. Any singleton components strictly between these two give at least 2^s distinct covers counted by U, where s is their number. To see this, choose Z to start at a tagged vertex of the first component and X to end at a tagged vertex of the last one, use full paths within those components, and distribute the s singletons freely between X and Z. Assign the rest of the exterior before the first selected component to X, and the exterior after the last to Z. Full paths with the required initial and terminal vertices exist in each strong component by rotating a Hamiltonian cycle; the singleton case is immediate. Thus if the endpoints lie in different components, s≤log₂ U.

Singletons outside the marked interval are untagged. Delete all of them, obtaining R*. The value g does not change. In the formula for W_s, the wholly untagged prefix is forced to X, hence contributes only its Hamiltonian count; the wholly untagged suffix contributes P₂ and each removed singleton there contributes a factor 2. The quadratic U similarly only sees full counts of the untagged prefix and suffix. Thus trimming preserves K and changes W only by a power of two. If Q is empty the slope is zero and cannot give a nontrivial ray.

The extent of the retained graph is bounded by its path count. For a nontrivial strong component with count g_i and size n_i, [Busch's lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) gives n_i≤3 log₃ g_i: it is exact for a triangle; for n_i≥4 use n_i≤4(n_i−1)/3 and 4/log₂5 < 3/log₂3. Consequently

    |R*| ≤ 3 log₃ g + 2 + log₂ max(1,U).

This formula also covers the case of one selected component. The [odd-pack expansion](https://arxiv.org/abs/2307.05569) gives P₂(R*)≤2^|R*| g: for each partition into two slots, disjoint odd-cycle packs in the two induced sets occur among the packs of the whole, so the product of the two Hamiltonian counts is ≤g, and there are 2^|R*| partitions. Therefore the **odd part of the slope** satisfies

    oddpart(W) ≤ 4 max(1,g−K) · g^(1+3 log₃ 2)       (W>0).

In particular for fixed g,K the possible positive slopes belong to finitely many dyadic rays, even if the graph order and all contacts are chosen from the target. For any bounded bank g≤G and K=±7 there are only O_G(log X) possible outputs of exact unit offset up to X. The offset being exactly ±7 means Kh=±7; since h is odd and cannot be 7, it forces h=1, so B is transitive.

For a target 7p with p≡3 or 5 mod 8, either unit-offset sign requires `oddpart(W)=7 oddpart(p∓1)≥7(p−1)/4`. Consequently a witness in this architecture must have

    g ≥ c·p^[1/(2+3 log₃ 2)]   (exponent about 0.257).

This is a lower bound on the **deletion count**, not just the number of vertices being used. It does not assume prime coverage or claim a similar bound for an arbitrary exact representation 7p=WP+Kh with a different intercept. The two-axis goal remains the [Chan–Pak spectrum conjecture](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf).

## The first deletion counts that pass the unit-offset test

The same trimming makes all orders finite when g is small. For g≤15, Busch bounds each nontrivial strong component by order six. The products g=9 and g=15 also allow respectively two triangles, and a triangle and a strong count-five core in either order. Enumerate those ordered strong sequences and insert at most 2+floor(log₂(g−K)) singleton vertices, testing every contact set. This suffices because larger singleton pads outside the selected interval do not change the odd part or K. Zero or negative g−K is handled as in the proof.

| g | odd parts of W divisible by 7, K=−7 | same, K=+7 |
|---:|:---|:---|
| 1,3,5 | none | impossible by K≤g |
| 9,11 | none | none |
| 13 | 49, 119 | none |
| 15 | none | 21 |

The g=7 row is forbidden. Thus **g=13 and g=15 are the sharp first deletion counts** for the negative and positive signs respectively. The literal first cores give all-order recipes (outgoing masks for R and mask Q):

| R | Q | g | (W,K) | output for B=TT_m |
|---|---:|---:|---|---|
| (44,53,56,50,33,0) | 36 | 13 | (49,−7) | 7(7·2^m−1) |
| (62,24,42,48,36,2) | 3 | 13 | (119,−7) | 7(17·2^m−1) |
| (6,12,24,17,3) | 5 | 15 | (42,7) | 7(6·2^m+1) |

The last odd part is 21, realized here with one extra factor of two. In each example a strong output is checked at m=1; homogeneously inflating that one block vertex preserves strength. The replay certifies the finite small-g table by subset-DP endpoints, and independently compares the contact ledger with the rooted read on all labelled cores through order four. The growth bound does not use the finite census.

Replay: `node replay_contact_passport.js`.
