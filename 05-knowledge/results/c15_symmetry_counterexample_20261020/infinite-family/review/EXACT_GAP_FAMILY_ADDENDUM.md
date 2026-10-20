# Exact deficit one in the infinite Paley-block family

**Verdict: PASS.** This strengthens the strict-inequality review in `INFINITE_FAMILY_REVIEW.md`. Under precisely the same hypotheses there is the exact formula

    tau_E(H)=tau_T(H)-1.

The extra argument is an intrinsic recovery of the two large independent blocks, up to sign-even permutations inside false-twin classes. All root-only permutations are already allowed in the normalizer, so they cause no additional sign escape.

## 1. N-evenness forces independent large blocks

The odd group H_p on a prime p=3 mod4 is transitive on unordered pairs, so an H-fixed graph restricts to either empty or complete on each large block. Inversion x->-x belongs to AGL(1,p), the normalizer. On K_p it swaps (p-1)/2 pairs, an odd number, and has orientation sign -1. Since external adjacency is uniform on the block, this automorphism extends with the same sign. Therefore an N-even H-fixed graph must be independent inside both A and B.

It remains to analyze graphs built from blocks of sizes a,b and k singleton roots, with arbitrary uniform crossing edges and arbitrary edges among the roots. This is a complete reduction, not an extra filter in the count.

## 2. Recover the large roles from false twins

False-twin classes are defined by equality of open neighborhoods. They are intrinsic independent sets. Any permutation inside a false-twin class is sign-even, because its transposition generators exchange nonadjacent twins.

Both A and B are contained in false-twin classes. There are two cases.

### The large blocks merge

If A and B belong to one false-twin class, it has size at least a+b. Every other false-twin class consists solely of roots and has size at most k<a. Thus the unique large class is intrinsically recognized. A full isomorphism sends it to the corresponding large class in the target graph.

### The large blocks remain separate

Write their class sizes as

    w_A=a+r_A,   w_B=b+r_B,   r_A+r_B<=k=b-a.

Every root-only class has size at most k<a. Also

    a<=w_A<=b<=w_B.

Therefore the A- and B-containing classes are individually recognized by their sizes unless w_A=w_B=b. Equality forces r_A=k and r_B=0: all roots are false twins of A, while B forms the other class. The whole graph then consists of exactly two independent classes of size b. If they were nonadjacent, they would be one false-twin class, contrary to the case assumption. Thus the sole equal-role exception is K_b,b.

This analysis is invariant under isomorphism. The merged and separated cases cannot be interchanged: they have respectively one or two false-twin classes of size at least a. For separated unequal roles, the two such classes have distinct sizes. Root-only class permutations are unrestricted but never confuse a root-only class with either large role.

## 3. Adjust an isomorphism into the normalizer

Let X,Y be two reduced H-fixed graphs that are not the K_b,b exception, and let g:X->Y be any full isomorphism. The role recovery above sends the large class or classes to the corresponding ones in Y. Permutations internal to Y's false-twin classes can therefore adjust g so that it sends each vertex of A to itself and each vertex of B to itself. In the merged case the same adjustment fixes both disjoint subsets within the single large cell; their images are disjoint and have the right sizes.

After adjustment, the map fixes A and B pointwise and permutes only F. It therefore belongs to S_F<=N. Hence full isomorphism does not merge distinct N-orbits of reduced graphs.

If g was an automorphism of X, every correction used above is an internal false-twin automorphism, hence sign-even. The adjusted automorphism has the same sign as g and lies in N. Consequently an N-even reduced graph other than K_b,b is globally even.

For the exception, there is only one H-fixed presentation: B is one whole side, and A union F is the other. Indeed B has b vertices and cannot share a bipartition side with A. Thus K_b,b contributes exactly one N-orbit and one full graph type. Its side swap is globally odd, while N preserves each side and is even, as proved in the strict review.

These facts prove both needed equalities:

- every N-even reduced graph except K_b,b is globally even;
- each globally even H-symmetric graph type corresponds to exactly one N-orbit.

The inherited restricted twisted Burnside theorem therefore gives

    tau_E(H)=#(N-even H-fixed graph orbits)-1
            =tau_T(H)-1.

No assumption about true-twin root classes is needed separately. If a permutation of such roots were odd, it would already belong to S_F and violate N-evenness. All root-only class permutations are handled by the same adjustment argument.

## 4. Independent finite control at order 22

`exact_gap_family_check.py` does not import the producer's checker. For a=7,b=11,k=4 it enumerates all 2^15=32768 normalized quotient graphs, quotients by the 24 root permutations, then independently computes each full weighted false-twin canonical form and every weighted-cell automorphism sign. Tournament quotient orbits are enumerated independently with the affine orientation action of those root permutations.

The exact results are:

    normalized presentations: 32768
    N-orbits of reduced graphs: 1992
    full types of reduced graphs: 1992
    N-even graph types: 1408
    globally sign-even graph types: 1407
    tournament types: 1408
    sole normalizer-even/global-odd exception: K11,11.

The counts 1992 concern the reduced independent-large-block graphs, not all H-fixed graphs including automatically odd complete internal blocks. The checker confirms no hidden graph fusion, and it identifies the exceptional graph by its two false-twin classes of sizes 11,11 and adjacency matrix [[0,1],[1,0]]. Its quotient mask is 481 under the documented six-node pair order.

    python3 -B exact_gap_family_check.py
    python3 -B -O exact_gap_family_check.py

Output and exact JSON agree byte-for-byte under optimization. The earlier n=10 cyclic example remains a separate smaller construction with deficit one. No minimal-order assertion is made.
