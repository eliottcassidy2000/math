# Session230: OCF is exactly the attachment cut kernel

Status: **PROVED**, identities below, conditional only on the established OCF identity. **FINITE-EXACT**, independent checks in `check.py` / `RESULTS.json`. No claim of all-odd realization, all-order first-creator genealogy, or new priority for OCF or the exclusion of 7 and 21.

## Inheritance and live concept board

Prior-art correction: THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature already proves the opposite Boolean curvatures and every mixed support difference below. THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension already supplies the split-path cut-field representation; THM-4114 section 6 represents the full response as a nonnegative s-t cut with two auxiliary terminals. These are inherited, not new. The refinement here is the explicit OCF even-path kernel on the original vertices, representing the excess divided by two, the quantitative bound (3), and its first-creator/min-cut applications. THM-4118-ear-response-lattice-and-stateful-unit-component-intervals forbids inferring solid cut spectra from a gcd condition.

Closest mechanisms: THM-4942, `05-knowledge/results/first_creator_generation_20261108/PROOF.md`, section 1 (split-path quadratic) and section 5 (typed analogy); THM-1900-insertion-response-calculus, neutrality; THM-029-alpha1-gap and THM-343-H7-impossible, repaired local strong-union proof; THM-1370-h-spectrum-omits-7-21-all-n, the trusted all-order exclusion. Session224 `session224_first_occurrence_20261107/README.md` gives the order/value ledger, not all-order realization. The historical body of THM-115-h21-permanent-gap contains superseded intersection assertions; we use only its OCF-component perspective and the later exclusion.

Canonical hostile: pairwise-intersecting cycles need not have a common vertex. Corrected near miss: deleting a vertex from a Hamiltonian path need not leave a Hamiltonian path. Least-used sidecar here: the even-cardinality endpoint path sum, weighted by the path count on the omitted vertices.

Board: odd-cycle polymers; attachment Boolean quadratic; induced-subset Boolean polynomial; condensation cuts; first-creator threshold. Anchor: turn section 5's analogy into an identity. Niche: positive coefficients for all support derivatives. Wildcard: distinguish the opposite curvature of the two Boolean spaces.

## 1. Exact positive cut kernel (PROVED; OCF refinement of existing cut calculus)

Write H(empty)=1. A directed path is an ordered sequence of distinct vertices. For distinct a,b of a tournament T on V, define

    K_T(a,b) = sum_{P: a-to-b directed path, |V(P)| even} H(T - V(P)).

Set K_T(a,a)=0. Add v beating exactly A subset V, and put y_a=[a in A]. Then

    H(T+v_A) = H(T) + 2 cut_K(A),
    cut_K(A) = sum_{a in A, b outside A} K_T(a,b).                 (1)

Proof. OCF writes H as the sum of 2^|F| over vertex-disjoint families F of directed odd cycles, including the empty family. Partition families in T+v by whether v is covered. Those not covering v contribute H(T). In any other family exactly one cycle C contains v. Removing C leaves an arbitrary odd-cycle family on T-V(C), contributing 2 H(T-V(C)). Cut C immediately after v: C is uniquely v,a,...,b,v, where the a-to-b path P has even cardinality, a in A, b outside A. Summation proves (1). No insertion of v into an old Hamiltonian path is assumed. QED.

This is not merely interpreting a graph analogy: the odd cycles in the exclusion proofs compute every coefficient in the generation operation.

### Coefficient bridge to the split-path calculus

Use THM-4942's convention x_a=[a -> v]=1-y_a, and its L_a, M_ab, W[a][b]. Uniqueness of multilinear polynomials on the Boolean cube gives

    L_a = 2 sum_b K_T(b,a),
    M_ab = W[a][b]+W[b][a] = 2(K_T(a,b)+K_T(b,a)).               (2)

In particular every L_a and every M_ab is even. This is stronger than deducing evenness just from odd-valued attachments: (2) supplies positive combinatorial witnesses for the halves.

A single arc a->b contributes the two-vertex path (a,b), so

    K_T(a,b)+K_T(b,a) >= H(T-{a,b}) >= 1.                       (3)

Therefore the mixed Boolean derivative of the extension function in coordinates a,b is the strictly negative constant

    -2(K_T(a,b)+K_T(b,a)) <= -2 H(T-{a,b}).                     (4)

Thus the attachment function is strictly submodular in every coordinate pair. This is an exact weighted directed-cut energy, not an analogy to one. The kernel entries themselves need not be symmetric; T=C3 already has K equal to its directed adjacency matrix.

### Neutrality, quantitatively

Equation (1) is nonnegative term by term. If there is an arc from A to its complement then its two-vertex contribution is positive. Conversely if there is no such arc then no directed path can leave A, hence every cut entry is zero. Thus neutral attachments are precisely terminal unions of SCCs, recovering THM-1900(B) without the invalid deletion-of-path premise.

For strong T of order n>=3, every nonempty proper A has positive cut weight. The smallest strictly increasing extension is exactly

    H(T)+2 min_{nonempty proper A} cut_K(A).                    (5)

It is a minimum directed cut problem with integer nonnegative capacities. Once K is computed, the minimum can be obtained by n(n-1) ordered source/sink minimum-cut computations: every nontrivial cut separates some ordered pair, and each pair-minimum is attained by a nontrivial cut; no claim is made that constructing K is polynomial-time.

## 2. A rigorous first-creator certificate from the same kernel

Let S(n) be the attained spectrum at order n and M_n=max S(n). For n>=3, every extension of a strong n-vertex T with a nonempty proper beaten set is a first creator at order n+1 if

    H(T)+2 mincut(K_T) > M_n.                                  (6)

Proof: (1) bounds every such extension strictly above every earlier count; padding makes S(n) contain every smaller-order count. QED.

This is sufficient, not necessary: an extension can land in an interior hole of S(n). The exact fertility condition is

    H(T)+2 cut_K(A) not in S(n).

Thus scalar first-count labels require the cut kernel (or an equivalent attachment sidecar) to control generation. This criterion does not say every strong extension is fertile, or that every first creator has a first-creator deletion.

## 3. How 7 and 21 constrain the same energy

Use the all-order exclusion THM-1370, not any conjecture about other odd numbers. For every realized base T,

    cut_K(A) != (7-H(T))/2 and (21-H(T))/2                       (7)

whenever the displayed target is nonnegative. For instance H(T)=5 forbids cut level 1; H(T)=19 forbids level 1 for the other reason. The consequence is exact but is not an independent proof of the omissions. Their local OCF obstruction (too few cycles once connected unions close under induced cycles) is a restriction on which positive cut kernels are realizable by tournaments. An arbitrary nonnegative directed cut kernel need not satisfy this restriction: on two vertices put K_01=3, K_10=0 and baseline 1, obtaining a formal output 7. Such an energy is not a tournament OCF kernel.

For a transitive base in the order 1->2->...->n the formula is explicit:

    K(i,j)=0 for i>j;
    K(i,i+1)=1;
    K(i,j)=2^(j-i-2) for j>=i+2.

Indeed paths use increasing subsequences, all residual tournaments have H=1, and even path cardinality means an even-size subset of the j-i-1 interior positions. Consequently all one-apex extensions of transitive tournaments are generated by a single upper-triangular binary-weight cut model.

**Direct noncircular subclass exclusion.** Cut levels 3 and 10 cannot occur in this model, independently of THM-1370. Encode A by a binary word, strip its initial zeros and terminal ones (which participate in no forward cut arc). A nonzero cut has a residual word starting 1 and ending 0, of length m>=2. If m>=7 the endpoint pair alone contributes 2^(m-3)>=16. Therefore a cut equal to 3 or 10 would have m<=6. Exhausting the at most 16 interior words at each remaining length yields:

- m=2: {1};
- m=3: {2};
- m=4: {4,5};
- m=5: {7,8,9,11};
- m=6: {12,13,15,16,18,19,22,23,25}.

Neither 3 nor 10 appears. This is a complete elementary proof of the omissions for all one-apex transitive tournaments, equivalently all tournaments whose cyclic-triangle clutter has a transversal of size at most one. It is not the general theorem. Thus the same clutter blocker condition provides both the generative construction (choose the apex and a cut word) and the exclusion certificate (the short residual word table). The finite table is verified in `check_transitive.py`; its maximum universe is 16 words, not a sample or bounded-order extrapolation. For T=C3 the kernel is the cyclic adjacency matrix, every nontrivial cut has weight 1, and every nonneutral attachment has H=5, another direct base control.

## 4. Inherited opposite Boolean geometry (THM-4114): support is absolutely monotone

Fix a tournament T. For each vertex-disjoint family F of odd directed cycles let U(F) be its vertex union. Define the squarefree polynomial

    Z_T(z) = sum_F 2^|F| product_{u in U(F)} z_u.

For a Boolean vector 1_S, Z_T(1_S)=H(T[S]) by OCF. All coefficients are nonnegative integers. Hence for disjoint B,R subset V,

    sum_{Q subset R} (-1)^(|R|-|Q|) H(T[B union Q])
      = sum_{F: R subset U(F) subset B union R} 2^|F| >= 0.      (8)

Proof: each monomial survives inclusion-exclusion iff it covers every marked vertex R and avoids vertices outside B union R. QED. Thus support addition is supermodular, and in fact every mixed support derivative is nonnegative, while attachment-orientation variation is strictly submodular by (4). Confusing these two Boolean spaces would reverse the inequality.

For R={v} and B=V-{v}, (8) specializes to the useful exact marginal

    H(T)-H(T-v) = 2 sum_{odd cycles C containing v} H(T-V(C)).   (9)

For larger R, cycles need not share one common point or be a single cycle: an entire disjoint family can jointly cover R. For example in C3 ordinal-sum C3, choose one marked vertex from each component. The mixed derivative is 4 even though no cycle meets both marked vertices. This is a hostile control against replacing the family-cover condition in (8) by a single-cycle intersection condition.

A set D hits all directed odd cycles iff T-D is transitive, since every nontransitive tournament contains a triangle. Hence the minimal supports of nonconstant monomials of Z_T are precisely cyclic-triangle vertex sets, and their blockers are exactly minimal feedback vertex sets. Full Z_T additionally remembers longer directed cycles and disjoint-cycle multiplicities, which a bare triangle clutter discards. This provides the precise transversal bridge: the clutter governs the zero-cycle boundary H=1; the weighted polynomial governs all count increments and the forbidden targets.

## 5. Loss ledger and scope

Source: fully oriented tournament and its odd-cycle families. Map: weighted even-path kernel K_T. Preserved target: the entire one-vertex attachment function, exact increments, neutral cuts, minimum positive increment, and fertility when S(n) is supplied. Lost information: a priori the vertex orientation itself, higher-vertex interactions, and module path-cover data; injectivity of K is not asserted. Needed sidecars: S(n) for first-creator decisions, labeled endpoints for attachment, full support polynomial or higher-marked data for multi-vertex operations. Cheapest decisive test: compare (1) against direct Hamiltonian-path DP on every attachment of every labeled small base; compare (2) separately against split-path DP.

No modular substitution formula is inferred from H or K alone. No Chan-Pak construction is assumed to generate every realized value or every tournament. The refinement is the explicit OCF even-path identification between two existing calculi and its scoped min-cut/creator consequences. The induced-support derivative theorem and opposite-curvature distinction are inherited from THM-4114, reproduced to make the bridge self-contained.

## 6. Reproduction

Run `python3 check.py` and `python3 check_transitive.py`. Universe: every labeled base tournament of orders 1 through 5 and every one-vertex attachment, so every labeled tournament of orders 2 through 6 appears. Independent path: direct subset DP counts Hamiltonian paths, while K sums explicit even-cardinality paths. Also checks split-path coefficients, strict pair inequalities, neutrality, forbidden outputs, and every disjoint B,R mixed support derivative. The script searches a hostile pair with equal base order and H but different attachment value profiles; details are in RESULTS.json. This is a small exact audit, not a replacement for the all-order proofs.

Audit summary: 1,099 bases; 33,866 attachments; 220,387 mixed support derivatives; zero validation failures. The equal-H hostile has n=5, H=3, bit codes 2 and 5 (lexicographic vertex pairs, bit 1 means smaller endpoint beats larger). Its extension-value sets differ, as recorded in RESULTS.json.

Execution-hardening audit: all validations in both scripts use explicit `require` calls raising RuntimeError, never Python assert. Normal and `python3 -O` output must each byte-match the corresponding frozen JSON.
