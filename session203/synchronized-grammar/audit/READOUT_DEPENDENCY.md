# Freely prescribed additive differences from endpoint-sensitive contexts

**PROVED, independently audited in the stated scope.** This is a constructive connection between inherited rooted amalgamation and a binary endpoint source. It does not hold a base count fixed, supply an operation on arbitrary Hamilton counts, or prove an eventual threshold. No new claim of priority is made for the rooted-amalgamation specialization.

## 1. A physical correction that is not the total module count

For a rooted nonempty tournament (U,v), let h=H(U), s=start_v(U), and L count ordered spanning two-path covers (P,Q), where P is nonempty and ends at v, while Q may be empty. Build E(U) by adding a with v->a and a->U-v. The elementary rooted-amalgamation formula gives

    H(E(U)) = h+L-s.

Let P(U) be the ordered join of a directed triangle before E(U). Thus H(P(U))=3(h+L-s). This parent is weak, since it is an ordered join of two nonempty tournaments. An equivalent rooted core has outgoing rows (2,0,11,19,7), root0, with every other core vertex dominating U-v.

Let C(U) use the rooted core with rows (2,8,3,5), root0, and again every nonroot core vertex dominating U-v. Preserve U literally. This four-vertex core has (h_G,e_G,R_G,d_G)=(5,2,6,3). The inherited one-sided formula

    H(G*U)=d_G h+(R_G-d_G)L+(e_G-d_G)s
             +(h_G-e_G-R_G+d_G)H(U-v)

therefore gives

    H(C(U))=3h+3L-s=H(P(U))+2s.                 (1)

This is a literal pair of tournaments, not permission to add an integer to an arbitrary tournament count. The parent has |U|+4 vertices and the child |U|+3, so this is a non-induced replacement construction, not a birth extension. C(U) is strong whenever U is strong: its core is strong and U shares the root; both induced strong subgraphs reach one another through that root.

The two fixed cores have equal (d,R,h-e)=(3,6,3), while their e coordinates are0 and2. Their deletion two-bin values differ:24 and12. Consequently this scalar equality is not a full-profile rebase and cannot simply be iterated as a closed endpoint fiber.

## 2. The correction can be prescribed as any positive even integer

For m>=1 form W_m from the transitive chain p_0,...,p_(m-1) and two additional vertices x,y, with

    chain -> x -> y -> chain.

The Hamilton-path start counts are

    start_(p_i)(W_m)=2^(m-i-1),
    start_x(W_m)=start_y(W_m)=1.

Indeed a path starting p_i has x,y consecutively. Every p_j with j<i lies after y, and every p_j with j>i independently chooses its increasing block before x or after y. This gives exactly2^(m-i-1). The two remaining initial vertices each have one forced path. Hence H(W_m)=2^m+1.

Given any integer s>=1, take m=bit_length(s). Add a root v that loses to x,y and beats p_i precisely when bit m-i-1 of s is1. Its Hamilton-start count is the sum of the start counts of its out-neighbors, and is therefore exactlys. Call this actual rooted tournament U_s. W_m is strong, and v has both an out-neighbor and an in-neighbor in W_m, so U_s is strong. Its order is m+3.

Applying (1) gives explicit odd-count tournaments P_s,C_s satisfying

    H(C_s)-H(P_s)=2s,
    |P_s|=bit_length(s)+7,
    |C_s|=bit_length(s)+6.

Thus the endpoint response is freely prescribed with logarithmic order. The parent count remains a correlated function of the actual U_s. These pairs do not yield an interval of possible child counts around any fixed parent.

## 3. Exact controls and remaining positive step

`free_additive.cpp` independently counts all labelled root-zero inputs of orders1 through4,75 inputs, by vertex-subset path DP. It computes L by actual two-path partitions, builds both literal cores, verifies the triangle-join parent and the+2s identity, and constructs all freely specified s from1 through127. It does not use the core tensor bank. Optimized and debug outputs agree.

    g++ -std=c++17 -O3 free_additive.cpp -o free_additive
    ./free_additive > FREE_ADDITIVE_NORMAL.json
    g++ -std=c++17 -O0 -D_GLIBCXX_ASSERTIONS free_additive.cpp -o free_additive_debug
    ./free_additive_debug > FREE_ADDITIVE_DEBUG.json
    cmp FREE_ADDITIVE_NORMAL.json FREE_ADDITIVE_DEBUG.json

A useful next target is an explicitly realizable two-interface family that can change s while holding the relevant base h+L-s fixed, or that provides a returning update of the entire required boundary profile. The freely prescribed difference by itself supplies neither. Endpoint-unit and uniform-incidence obstructions still apply in their stated classes; this multiplier-one pair does not contradict them.

Inherited source: session197/path-polynomials/ROOTED_ORDERED_AMALGAMATION.md, equations1-4. The exact same+2s context pair was concurrently recovered by the endpoint and two-port workers. This proof records the free-endpoint construction and the precise remaining coupling, rather than claiming a new tensor-search result. The exploratory negative companion from the source note is excluded from this release.
