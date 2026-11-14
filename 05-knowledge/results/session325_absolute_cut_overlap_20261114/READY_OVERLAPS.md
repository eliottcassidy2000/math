# Repeatably ready overlaps and adjacent strong corridors

Use the ordered [Busch seam](../session324_busch_inverse_and_gap_structure_20261113/BUSCH_INVERSE.md). For a strong tournament A rooted at v, write H for its Hamilton count, h=H(A−v), and L for the number of ordered two-path covers of A−v with second slot empty or starting at an out-neighbor of v. A right near-transitive block R_j (j≥1), glued at its first vertex, adds j+1 vertices and acts at the **same** old root by

    (H,h,L) → (H+2^j L, h, 2^(j+1)L−h).

The new tournament is strong. The whole update, not only its first coordinate, is essential below.

## Two different strong cores synchronize at the same absolute target

The following rooted strong sources are neither isomorphic nor converse, as a direct check of the six-vertex permutations verifies:

| outgoing masks | root | H | h | L |
|---|---:|---:|---:|---:|
| (2,48,35,7,13,25) | 0 | 33 | 15 | 36 |
| (48,9,3,5,14,30) | 4 | 33 | 15 | 72 |

Use j=2 on the first and j=1 on the second. Both resulting rooted states are exactly (177,15,273). Their orders differ by one, and continue to do so under any common later attachment word. Thus arbitrary future seam choices preserve equality at every step despite different starting cores, not merely along one chosen sequence. More generally j versus j−1 synchronizes the full states for every j≥2, giving the common absolute ray 33+36·2^j. The coalescence is possible because the two source L values are in dyadic ratio, the counts H agree, and the shared h makes the post-step subtraction agree as well.

## Three consecutive odds persist at every later step

Three literal order-seven sources, all at root 5, have the following profiles:

| outgoing masks | H | h | L |
|---|---:|---:|---:|
| (72,17,67,38,13,23,58) | 135 | 45 | 303 |
| (72,5,73,50,7,23,58) | 137 | 45 | 303 |
| (66,24,67,37,13,23,58) | 139 | 45 | 303 |

At every depth, apply any common word of right near-transitive attachments, always at the old root. The h and L responses stay identical on all three tracks. The outputs are therefore three consecutive odd counts at the same order, indefinitely, with both gaps exactly two. Unlike earlier adjacent pairs with different slot responses, this triple matches the **entire active seam state other than the intended count gap**. An order-five pair already has this property: rows (8,17,3,6,13) rooted at 3 and (16,9,3,5,14) rooted at 0 give (H,h,L)=(13,3,15),(15,3,15).

The order-seven triple even reaches every allowed residue simultaneously modulo 2^a 3^b 5^c using these same physical rooted moves. Choose its initial member by the target residue modulo 3 if b>0. Throughout the moves, h=45, L stays odd, v₃(L)=1, and L remains a 5-adic unit. First fix successive binary digits: if bit j≥1 is wrong, one R_j changes precisely that bit and no lower bits, because its increment is 2^jL and L is odd. Use j only through a−1. Next choose a common j≥a such that 2^(j+1)≡1 modulo the odd part 3^b5^c. Repeating this same attachment t times gives, modulo the odd part,

    L_t=L−th,    H_t−H=[tL−h·binom(t,2)]/2.

Modulo 5 the polynomial has nonzero linear coefficient L/2 and all quadratic corrections divisible by 5. Modulo 3, divide the difference by 3 first: L/3 is a unit and h/3 is divisible by 3. At each prime, shifting t by a digit times a prime power changes the next digit by a unit times that digit. Thus all necessary higher digits lift uniquely after selecting the starting residue, and CRT chooses one nonnegative t. These repetitions leave the prescribed binary residue fixed. All three tracks continue to realize consecutive odd counts. This is modular coverage, not equality to an arbitrary prescribed absolute count.

## Exact word algebra, disjoint size shells, and an asynchronous collision

For a word of attachments put r_i=j_i+1≥2, s_i=r_1+…+r_i, and R=s_k, the total number of new vertices. Iteration gives

    L_k=2^R L_0−h Σ_i 2^(R−s_i),
    2(H_k−H_0)=L_0 Σ_i 2^s_i−h Σ_(i<j) 2^(s_j−s_i).

The last identity uses 2(H_i−H_(i−1))=L_i+h. The physical histories are exactly the nonadjacent binary breakpoint sets {s_i}: the first breakpoint is at least 2 and all successive differences at least 2.

For one fixed initial root, the size shells are disjoint. Indeed L_0>h (the second slot can start at an out-neighbor), and summing the separated powers gives

    (L_0−h/3)2^R ≤ 2(H_k−H_0) ≤ (4L_0/3)2^R.

The upper bound for R is strictly below the lower bound for R+1 because L_0>h. Thus an exact target forces the total added order even when the history itself is nonunique. The number of compositions of R into parts at least two is Fibonacci-sized, O(φ^R), where φ is the golden ratio. Each fixed source therefore gives at most O(X^(log₂φ)) counts through X, even allowing every history. This proves the permanently ready triple does not cover the absolute spectrum merely by using all its words.

Within a fixed shell there is genuine exact overlap. Compare one block of size R and two blocks of sizes s,R−s, each at least two. The difference in twice the output count is 2^s L_0−h2^(R−s); hence equality occurs **if and only if** L_0/h=2^(R−2s). Take the strong order-four source (8,1,3,6) rooted at 1, with (H,h,L)=(5,3,6). For every s≥2 set R=2s+1. Both the one-block history and the split (s,s+1) give

    H=5+6·4^s,   order 2s+5.

Their L values differ. They are also distinct up to converse as unrooted tournaments: the one-block graph's outdegrees are 1 through R+2 with extra copies of 1 and R+2, a symmetric multiset. In the split graph the first added block's degrees are shifted by s+1, the last block has degrees at most s, and old vertices have degrees at least R−1. Hence it omits degree s+1, so neither it nor its converse has the one-block score multiset. This is a same-source, same-size, exact infinite intersection of two construction histories.

All finite seed profiles, the first coalescence, a common word on the triple, and the first asynchronous collisions are checked by literal subset DP in `replay_ready_overlaps.py`. The infinite statements follow from the seam identity and the exact arithmetic above, not from finite agreement.
