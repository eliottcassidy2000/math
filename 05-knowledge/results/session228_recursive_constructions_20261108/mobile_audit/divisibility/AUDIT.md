# Independent audit: sharp triangle and universal SCC divisibility

Verdict: **PASS with no corrections required** on frozen `inputs/DIVISIBILITY.md` and `inputs/GENERAL_DIVISIBILITY.md`. Earlier audit packets remain unchanged. This packet includes an independent standard-library verifier, normal and optimized receipts, frozen inputs and hashes.

## Endpoint-block mechanism

For any ordered decomposition of the root-deletion graph into nonempty tournament blocks, deleting the root from a Hamilton path yields an ordered pair(P,R), with P ending at a permitted inneighbor when nonempty and R starting at a permitted outneighbor when nonempty. Fix the first occupied R block i and the last occupied P block j, using missing markers for empty paths. This partitions paths into disjoint cases; it is not an inclusion-exclusion argument.

The at-most-two boundary blocks have arbitrary local endpoint restrictions. Every other block contributes H or Q, where Q counts two ordered paths covering the block and permits either path to be empty. For i<j, the strictly intermediate blocks split freely, earlier blocks lie wholly in P, later blocks wholly in R. For i=j all other blocks are assigned wholly. For j<i an intervening block makes the case impossible; if there is no gap, P and R occupy a prefix and suffix. A missing path leaves at most one constrained block. With both missing the deletion graph is empty. These cover every case, including empty local path pieces inside the i<j interval.

All local nonempty pieces concatenate uniquely because cross-block edges point forward. Thus the product factorization is valid: a Hamilton path cannot leave a block and later return to it. Boundary constraints cannot secretly impose new conditions on a third block.

## General arithmetic conclusion

Writing g_i=gcd(H_i,Q_i), each case count is divisible by the product of all g_i outside its boundary set C. Therefore the gcd of those products, allowing every |C|≤2, divides the total D. This gcd equals

    K = (product g_i) / lcm { product_(i in C) g_i : |C|≤2 }.

Every boundary product divides the full product, so the denominator divides it as well and K is an integer. Prime by prime, its valuation is the sum of all valuations minus the two largest, padding by zeros. Distinct indices matter: the same block cannot be removed twice from the product. The frozen producer explicitly handles this correctly.

For zero or one block, and also for two blocks, K=1. The empty-root-deletion case has D=1. No singular convention or false nontrivial bound is introduced at these boundaries.

The SCC condensation of a tournament is itself a transitive tournament. Thus the ordered-block theorem applies canonically to every rooted tournament using its deletion SCCs. Strongness of individual blocks is not needed for the proof, only for canonical selection. This is genuinely a universal necessary condition, not a theorem only about the triangle construction class.

The prime-count corollary follows directly: if D=q is prime, for p≠q at most two deletion SCCs can have p dividing g_i. For p=q there can be at most three, and if there are exactly three, the smallest of their positive q-valuations is1. Larger retained sums contradict v_q(D)=1. This is a restriction on each root-deletion decomposition, not a sufficient condition for D=q.

The gcd cannot be replaced indiscriminately by H. The strong four-vertex tournament has H5 and Q32, hence g1. Independently computed cyclic regular tournaments give order5 signature(H,Q,g)=(15,120,15), and order7 signature(175,1988,7). These heterogeneous signatures also show the obstruction is not merely a restatement of the triangle valuation.

## Triangle bound and sharpness

A triangle has(H,Q,g)=(3,12,3), and a singleton has(1,2,1), giving3^max(t−2,0)|D. The two constrained blocks can use up two triangle factors, exactly matching the stated exponent.

The two-triangle adjacent hostile has D41. Prepending t−2 triangle components to the entire rooted graph makes an ordinal sum, so D becomes3^(t−2)41. Its deletion still has t triangle components. This proves sharpness for every t≥2, with no implicit strongness claim for these sharpness witnesses. The separated family has D=3^(t−1) times a factor congruent to−t modulo3, so its valuation is exactly t−1 whenever3 does not divide t. This latter claim inherits MULTI's t≥1 scope. It does not predict the exact valuation when3 divides t.

## Independent verification

`check.py` imports only the earlier independent audit module `../verify.py`. For each chosen block it directly enumerates all local directed paths, computes H and Q, and builds its root-neighborhood six-state matrix from local path pairs. It tests the divisor as the gcd of complementary products, an independent formulation of the producer's lcm quotient. Literal controls build full adjacency matrices and count paths by subset DP.

Reproduction from this directory:

    python3 check.py > normal.json
    python3 -O check.py > optimized.json
    cmp normal.json optimized.json
    sha256sum -c SHA256SUMS

Both runs PASS and outputs are byte-identical. Checks use explicit errors, not assertions.

- 6,000 deterministic seeded heterogeneous block words, lengths1–8, including arbitrary root masks for all blocks of order≤5 and four selected order7 masks. 3,299 have a nontrivial predicted divisor.
- Five independently computed block signatures: singleton, C3, strong4, cyclic5 and cyclic7.
- 100 literal profile controls with at most10 deletion vertices.
- Three additional nontrivial literal controls: three triangles yield D417,K3; triangle+cyclic5+cyclic5 yield D37305,K3; three cyclic5 blocks yield D239085,K15. Exact root masks are retained in the receipt.
- Exhaustive9,330 words of lengths1–5 over the two singleton root types and four triangle cardinality types; every triangle divisibility bound passes.
- Explicit sharpness values3^(t−2)41 for every t2–40.

The finite checks corroborate the endpoint-block proof; they do not establish universality by themselves. No literature-priority claim was audited or asserted. Passing the universal divisibility condition does not establish physical membership, a prescribed Hamilton count, or Chan–Pak completeness.
