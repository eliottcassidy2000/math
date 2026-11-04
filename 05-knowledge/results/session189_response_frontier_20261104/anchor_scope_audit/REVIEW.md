# Separate scope and battery review of incoming THM-4864

## Scope of this review

This note independently challenges the boundary, exposed-end hypothesis, finite query bound, and source-relative integration of the incoming anchor proof at commit 97632252. It does not repeat the independent weighted-cycle enumeration and physical control banks in ../session189-incoming-anchor-audit. The incoming exact proof has Git blob 51390c8fe6423e7bd5a7326495ae81e01121c3aa; SOURCE.json pins the retained bytes and upstream URL.

The anchor proof is mathematically sound with an explicit quotient-order restriction m≥2 in Theorem B and a clarified branch-specific battery. Generic quotient-converse separation for arbitrary children is not established. The statements below use the incoming proved cycle contraction and weighted E1 identity, separately audited in the companion packet.

## 1. The order-two boundary is valid; order one is not

Let m be the number of children of the quotient Q. Splitting one distinguished child into its exposed endpoint u and the remaining block beta gives a weighted kernel of order n=m+1. The weighted E1 identity requires n≥3, so m≥2 suffices, including the two-child case.

Theorem B must state m≥2. At m=1, take A'=A=TT3. Q=Q^op is a singleton and W'=W on every physical name, so no query separates. The proof's step n=m+1≥3 silently excludes exactly this boundary. Restricting the original-module application to a highest prime node makes the boundary harmless there, but does not fix the standalone theorem's wording.

For m=2 every split kernel has three vertices. If the other child is a singleton, its weighted E1 functional works. Otherwise the hard-case proof applies unchanged: for k≥3 the two-vertex quotient's endpoint parity at the distinguished child is one, supplying the coarse fallback; for k=2 the mixed-coefficient argument is valid. No order-two exception remains after adding m≥2.

## 2. The actual hypothesis is a common physical exposed end

The transitivity assumption can be strengthened as follows.

**Common-end lemma.** Let Q have m≥2 children. Compare Q[A_i] with Q^op[A'_i] on the same child sets, with |A_i|=|A'_i|. Suppose a distinguished child has at least two names and a physical name u is a source or sink in BOTH A_j and A'_j. No transitivity is required of either remaining child interior. Then an actual two-port query separates modulo four, uniform on every other child and on beta=A_j−{u}.

**Proof.** Because u is an end in each child, beta is a module on both sides. Its internal tournament can be arbitrary; the port-cycle contraction uses only its nonzero cardinality, assigning weight [|beta|=1]. The other children similarly contract using their cardinalities, and the exposed name u has weight one. The source split kernel is Q[u→beta] after choosing the source-end orientation; the sink case follows by whole converse and port converse-complement, preserving query cardinality and uniformity.

If u changes from source to sink, the target weighted kernel is the full converse. If u remains a source, the kernels differ by quotient converse while retaining u→beta. These are precisely the two cases of the incoming proof. Its swap argument is still valid when u and beta have different weights: only cycles using the flipped u-beta edge survive, and every such cycle contains both names, so swapping the two preserves total weight. Unchanged cycles cancel before this swap. The endpoint-parity fallback uses only the quotient, not the discarded beta interior. Thus every proof step survives verbatim. QED.

This is a pairwise common-end condition. It does not say that an arbitrary target preserves a source endpoint, or that every nontransitive child can be separated by this argument. Under a full-response collision, it implies that no physical child vertex may remain an exposed end on both sides.

## 3. A complete source-fixed battery within the claimed bound

For a kernel on n vertices, each coefficient in L_v is a quadratic [a_i c_v] or a cubic [a_i a_v c_v]. Its Boolean finite-difference cube uses c=0 or c=e_v and a=0, singletons, and selected pairs e_i+e_v. Their union has at most4n distinct actual queries.

There are two source-determined branches.

1. If another child is a singleton z, use the L_z battery. The reversed-end target has contrast equal, up to a constant, to the full-converse contrast and L_z=1. The unchanged-end target agrees with it on tau-orbit sums, and L_z is tau-invariant because z sees u and beta identically. Thus the SAME L_z battery detects both target end directions. At most4n queries suffice.
2. If no other child is a singleton, use the union of the coefficient cubes for L_u, L_beta, and the extra cubic [a_u a_beta c_beta], taking u→beta. Let d be the number of outside common outneighbors of u,beta. The exact union contains12 internal queries,7 for each outside outneighbor, and3 for each outside inneighbor. Thus its cardinality is

       12+7d+3(n-2-d)=3n+6+4d≤7n−2<9n.

   The two coarse fallback queries (a,c)=(0,0) and (e_u+e_beta,0) are already in the twelve internal queries. The same union covers the reversed-end target through L_u and the unchanged-end target through the stated hard-case functionals.

Consequently the claimed9(m+1) bound is safe. The incoming battery sentence should explicitly mention the singleton-z branch rather than naming only L_u,L_beta and the extra cubic. The above construction depends on the source and the fixed physical endpoint, not on the target's choice between its two endpoint directions.

verify.py independently enumerates the actual finite-difference query supports for kernel orders3 through64 and every possible common outside outdegree. It checks the exact7/3/12 formula, the coarse fallback inclusion, and the4n one-functional bound. These finite checks illustrate the elementary union count, not the weighted parity theorem.

## 4. Highest-node collision integration

Assume an original-source-module XOR target U has the same full named response as T. Then H(T)=H(U), and the reviewed hereditary equality theorem supplies equal complete cover profiles on every fixed named subset. Highest changed nodes of the ORIGINAL strong-module tree form a disjoint antichain, and their boundaries are unchanged.

At a highest changed prime node M, its quotient becomes its converse, and its original physical child sets remain modules in the comparison. The anchor theorem itself does not require their profiles to agree. It separates when there is a singleton child or a common exposed child endpoint. Give both ports the original incidences of M on all exterior names. The reviewed modulo-four ambient lift preserves the nonzero gap. Every other highest changed node is a module in that same actual augmentation and has equal complete profile, so restore it exactly, one disjoint node at a time. This yields the original target U at the same physical query and a contradiction.

For a transitive source child A, hereditary equality forces the target child to be transitive because its Hamilton count remains one. For |A|=2 every target order is A or A^op. For |A|=3 the two source ends and two target ends must intersect. For |A|≥5 the separately reviewed distance theorem restricts the target physical order to A or A^op. The common-end anchor therefore excludes all these cases. At |A|=4, absence of a common endpoint leaves exactly the orders whose two ends are the two original middle names. This note does not eliminate those crossed orders.

The frontier is only a necessary condition on highest changed original nodes. It does not assert a collision in the remaining cases, pairwise injectivity among arbitrary outputs, or validity for sequences using newly created modules.

## 5. A source-relative root-family battery needs the scalar step

For a fixed source core in any fixed ambient, a target-independent root-family battery can combine the existing ingredients. Include the core-zero query first, lifting both ports with the core's exterior incidences. Original-module monotonicity makes all relevant profiles coordinatewise nondecreasing, and the occupied-run Hamilton coefficient is positive. If the core Hamilton count increases, that actual query separates. Otherwise hereditary equality gives complete child-profile equality and keeps every transitive source child transitive.

Only AFTER that equality reduction may the distance battery for a child of order at least five be invoked to restrict its target order to source or converse. The anchor then handles both. A two-name child needs one source-end anchor. A three-name child needs the union of BOTH source-end anchor batteries, because which source end remains a target end can depend on the target. This union is still fixed from the source. A singleton child has its direct coarse battery.

It would be false to assume that every original-module XOR target of a transitive child stays transitive without the equality reduction. In TT3, toggling both original intervals {0,1} and {1,2} creates a directed triangle. The scalar step detects that strict increase instead.

All of this restricts operations to the original core's source modules. For collisions with other changed highest nodes, use the exact disjoint restoration in Section4, not a claim that the root-family battery automatically absorbs arbitrary exterior changes.
