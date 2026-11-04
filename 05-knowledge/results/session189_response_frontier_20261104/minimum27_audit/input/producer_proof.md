# The 27-name cyclic-child boundary has no original-module response collision

**Finite-exact result and source-relative integration, pending independent response/scope review.** The quotient census completeness and parity filter have a separate independent exhaustive certificate. This packet does not claim unrestricted two-port response injectivity, closure of every order-nine quotient with arbitrary children, or an arbitrary-ambient lift of its new modulo-eight scalar gaps.

## 1. Physical statement and prior inputs

For a named tournament T, let F_T(a,c) be the number of Hamilton paths after adjoining physical ports x->y, with a_v=[v->x] and c_v=[v->y]. Let Q be any prime tournament on nine names. On the 27 old names 0,...,26, give child i the three names 3i,3i+1,3i+2 and the directed cycle

    3i -> 3i+1 -> 3i+2 -> 3i.

Between children use Q. Denote the resulting source by T=Q[C3,...,C3]. The quotient-only target is U=Q^op[C3,...,C3], with all nine physical interiors unchanged.

**Finite theorem.** F_T differs from F_U for every such Q. After transporting Q to one of the certified isomorphism/converse representatives, a fixed battery of at most 20 child-uniform actual queries separates them. Nine of those are endpoint queries. The remaining 12-query battery has one query in common with those nine, and suffices for every prime representative with all endpoint parities zero.

**Source-relative consequence.** Combining this finite theorem with the separately approved inputs below, a nontrivial full named-response collision obtained by simultaneous XOR of modules of the original source requires at least 28 old names. This integration is subject to the separate review identified above.

The exact prior inputs, pinned in DEPENDENCIES.json, are:

1. The approved Session189 highest-node frontier: in such a collision every highest changed original strong quotient is prime of order at least nine. Its proof includes hereditary complete-profile equality, exclusion of linear nodes, and restoration of other highest changed nodes at the same actual query. See frontier_nine/PROOF.md and nine_integration_audit/REVIEW.md.
2. The approved transitive-child synthesis in tt4_cross_audit/REVIEW.md Section 6: every child of a surviving highest changed prime node is nontransitive. It combines the corrected incoming common-end theorem, the transitive distance result, the general conversed-child lemma, and the crossed-TT4 positive quartics. Its original-module and same-query exterior qualifications are essential.
3. The all-order endpoint parity theorem recorded and reviewed in prime8/PROOF.md Section 2 and prime8_audit/REVIEW.md Section 1. For p_i=S_i(Q)+E_i(Q) modulo two nonzero, the actual child-uniform query (a,c)=(e_i,0) has source-minus-target contrast 2 modulo four. It applies to arbitrary nonempty children with matching profiles; no order-eight bound is used in this particular lemma.

Incoming E4 already establishes local nontrivial module-converse detectability. The new finite work here concerns the quotient-only converse with unchanged cyclic interiors, not another proof of that local theorem. No weighted-deletion reconstruction claim is used.

## 2. Complete quotient census, independently certified

The compact census.tsv retains every node ID and its exact 36-bit adjacency code from the reproduced incoming order-nine census. Decode bits in descending order along (i,j)=(0,1),(0,2),(1,2),...,(7,8); a one means i->j. This gives exactly 97,144 representatives modulo isomorphism and converse.

Completeness is not assumed from a reported total or a hash. The separate c3_census_audit independently canonically labels every representative and its converse using only sorted outdegree cells and all permutations inside each cell. Every tournament isomorphism preserves outdegrees. Therefore minimizing the full row-major edge code over these permutations gives an exact isomorphism key, and the number of permutations giving the minimum is exactly the automorphism-group order.

The audit checks that each pair of converse keys belongs to exactly one merged representative. It obtains 191,536 distinct ordinary isomorphism classes, of which 2,752 are self-converse. It then checks the exact disjoint labeled-orbit sum

    sum over ordinary classes of 9!/|Aut(Q)| = 68,719,476,736 = 2^36.

The orbits are disjoint subsets of all labeled nine-name tournaments and have the full cardinality of that universe. This proves exhaustive coverage without trusting the incoming canonicalizer or inductive generator. The independent audit examines 33,867,002 degree-cell permutations at order nine and retains every canonical key and automorphism count. It also reconciles the source generator's ordinary-extension construction, but that is not needed for the orbit-mass completeness certificate.

The producer tests all proper nonsingleton physical subsets for the module property. Independently, the auditor starts with every pair and repeatedly adds any outside name distinguishing the current set; a proper fixed closure is a module, and any proper module contains a seed pair whose closure stays inside it. Both methods give exactly 69,445 prime merged representatives.

For every prime representative, the producer's predecessor subset recurrence computes all nine start and nine end counts. The auditor independently enumerates every directed Hamilton path by literal depth-first traversal. All 1,250,010 prime endpoint coordinates agree. Exactly 67,097 prime representatives have some p_i nonzero and are covered by input 3; the remaining 2,348 have p=0. filter.tsv retains the full module/endpoint data. endpoint_witnesses.tsv records the chosen p_i, its actual expanded physical mask, and the certified contrast 2 modulo four for every excluded prime row.

## 3. Exact C3 run-word model modulo eight

Let q_h(A) count ordered covers of A by h nonempty directed paths. For a directed triangle the complete profile is

    q_1(C3)=3, q_2(C3)=6, q_3(C3)=6, q_h(C3)=0 for h>=4.

The first value counts its three directed Hamilton paths. A two-cover consists of one directed edge and the remaining singleton in either order, giving six. Three singleton paths have six orders.

Take an actual query uniform on every child. Decompose any augmented physical Hamilton path into maximal consecutive runs in its children, together with the single occurrences of x and y. Consecutive equal child labels are forbidden; all cross-run arcs are determined by the augmented quotient. If child i occurs k_i times, the physical fillings of its ordered runs are counted by q_(k_i)(C3). This is a bijection, including all named internal edges and all port placements.

Consequently its contribution after dividing out the common odd factor 3^9 is one for a child used once and two for a child used twice or three times. A word with at least three repeated children is divisible by eight and contributes zero modulo eight. No degree truncation of F is made. In particular, words containing triple runs, mixed double/triple runs, and all cubic/quartic response effects are retained whenever their weight can be nonzero modulo eight.

The producer's finite-state recurrence uses nine visit counts in {0,1,2,3} and two port counts in {0,1}. Keep exactly those states with at most two old children whose visit count is at least two. There are 94,208 such states, including the empty state and all four port-occupancy patterns. Initialize every single-symbol word. Append a symbol only along an actual augmented-quotient arc and below its capacity. The transition weight is two on an old child's second visit and one otherwise. Sum states in which every old child and both ports have appeared, over all final symbols; multiply by 3^9 modulo eight.

Every transition strictly increases the packed state code, so the dynamic program is acyclic. Summing complete states of every allowed length includes each admissible run word once. Missing states have weight divisible by eight, as proved above. This establishes the exact physical response residue, not an endpoint-product approximation or an extrapolation from small children.

For coarse nine-bit masks a,c, the augmented quotient rows are explicitly

    r_i = Q_i | (a_i << 9) | (c_i << 10),
    r_x = (511 xor a) | 1024,   r_y = 511 xor c.

For the target replace Q_i by 511 xor Q_i xor 2^i, and keep the same a,c. Physical masks replace each selected child bit 2^i by 7*2^(3i). Ports are the fixed names x=27 and y=28, with x->y on both sides.

## 4. The complete fixed battery and finite result

On the certified labels the 12 coarse (a,c) queries, in order, are

    (0,1), (1,0), (1,2), (1,3), (3,1), (3,3),
    (3,5), (0,7), (0,18), (40,256), (4,3), (130,128).

battery.json records both coarse and expanded physical masks. residues.tsv retains all 56,352 residues: source and target at each of these 12 queries for all 2,348 p-zero prime representatives. Every representative has a nonzero difference. Every nonzero difference in this bank is 4 modulo eight.

The first separating query histogram is

    index 0: 1159; 1: 630; 2: 296; 3: 133; 4: 68;
    5: 0; 6: 30; 7: 15; 8: 7; 9: 4; 10: 4; 11: 2.

The last two rows first separated at index 11 are node 4615, with residues (7,3), and node 78886, with residues (3,7). physical_witnesses.json retains every one of the 2,348 selected queries, both complete 29-name augmented outgoing-row lists, and both residues. It is not merely a table of unlabeled classes or a digest of the calculation.

Adding all nine endpoint queries (2^i,0) gives 20 distinct actual coarse queries, since (1,0) already appears above. This union covers all 69,445 prime merged representatives. The 12-query claim alone is restricted to the p-zero subdomain.

An arbitrary labeled prime Q is isomorphic either to a retained representative R or to R^op. Apply the same physical child permutation to source and target. In the converse case their roles exchange, preserving nonvanishing. Within each cyclic child, simultaneously choose a relabeling to the displayed three-cycle; every directed triangle is isomorphic to it. Because the target interiors are unchanged, this is the same relabeling on both sides. Thus no independent source/target relabeling is used.

This is a bounded exact theorem about C3 children at order nine. It makes no statement about arbitrary nine-quotient child profiles. It also does not assert that a scalar four-modulo-eight gap survives arbitrary ambient substitution: no such lift is needed below.

## 5. All original-module families at the 27-name boundary

Suppose T and a simultaneous XOR target U from modules of the original T satisfy F_T=F_U but T!=U, and suppose |T|<=27. Input 1 supplies a highest changed original strong node M with prime quotient of order k>=9. Input 2 says every source child of M is nontransitive, hence has at least three physical names. Therefore

    27 >= |T| >= |M| >= 3k >= 27.

Every inequality is equality. In particular |T|=27, M is the whole source, k=9, and all nine children are directed triangles. There is no exterior ambient and no other disjoint highest changed node. This equality argument is why the new scalar finite gaps need no all-cover or ambient certificate.

A prime quotient has no proper original modules meeting more than one child. A C3 child has no nontrivial proper modules. Ignoring singleton selections, which change no edge, a family with the root selected oddly is therefore the root together with a subset of the nine whole children. Its target quotient is Q^op. Each target triangle is either unchanged or its literal converse, with the same complete cover profile in either case.

All 512 possible vectors of target triangle orientations have the same complete child profile (3,6,6), child by child. At every child-uniform physical query the occupied-run bijection of Section 3 therefore gives an exact integer equality between the actual family target and the quotient-only target with unchanged interiors. This is not an independent relabeling of either response: the physical names and port masks stay fixed, while the exact ordered-cover factors agree.

Consequently the same 20 coarse queries from Section 4 separate the source from every root-containing original family. The scalar residue for a chosen query is independent of which target triangles have been conversed. The quotient-only target is the special family selecting the root and all nine children, whose internal triangle toggles cancel child by child.

Thus no nontrivial original-module full-response collision exists on at most 27 old names, conditional only on the explicitly approved inherited inputs and the present finite certificate. Equivalently any still-possible such collision needs at least 28 old names. This does not exclude collisions between arbitrary tournaments or arbitrary pairs of family outputs, and does not authorize inversion operations on modules created after earlier operations.

For a fixed 27-name source Q[C3^9], the 20-query battery is thus source-dependent but family-uniform over root-containing original families. In the p-zero subdomain its 12-query part suffices. No additional exposed-child battery or nonuniform attachment is needed for this final finite family theorem. The general conversed-child lemma remains an inherited ingredient of the already approved transitive-child exclusion, not a new ingredient here.

## 6. Same-query restoration and the exact novelty boundary

The previous 27-name argument leaves no outside names. For completeness, its use of the inherited all-transitive-child exclusion is not a tacit assumption that other highest changes cannot cancel a local gap. That exclusion is already reviewed with the following exact restoration: full equality implies hereditary complete-profile equality; other highest changed original nodes are disjoint; both ports imitate the chosen core at every exterior name, so they are uniform on those other nodes; equal-profile replacement preserves the actual queried count at that same query. No modular gap is added to an unrelated scalar gap.

At the final new step there are no such nodes to restore and no virtual contracted exterior. Contracting C3 interiors is only the exact uniform-port run-word bijection of Section 3. The physical quotient and both ports are restored explicitly by the masks and 29-name row lists in the packet.

The contribution is complete finite closure of the minimal cyclic-child boundary, together with this source-relative size integration. It is distinct from incoming local-converse E4, weighted deletion decks, unweighted deck collisions, or a generic prime-nine/nontransitive-child reconstruction theorem.

## 7. Reproduction and review boundaries

Run python3 verify.py and python3 -O verify.py. The script compiles the retained C++17 sources with the already installed compiler, recomputes the complete producer module/endpoint filter and the full 12-query bank, checks every physical row/mask, reconciles every first witness, and replays the retained prior exact-profile controls. Explicit exceptions remain active under optimized Python. No software installation or network access is needed.

The census auditor supplies its own exact compact-census replay and separate prime/endpoint calculation. The response/scope auditor uses a different physical-subset recurrence modulo the independent rotations of each C3, with no run-weight normalization or at-most-two-repeats truncation. That audit recounts every selected witness, not necessarily every redundant entry of the full 12-query bank. Its status and scope must be read from its separately frozen review, rather than inferred from this producer document.

The bank records actual values of the complete response at selected physical masks. It is not a full 27-name Boolean polynomial bank, and no claim is made to retain every individual named coefficient for all 69,445 source types. All cubic and quartic contributions to every retained actual query are included exactly by the run-word bijection.
