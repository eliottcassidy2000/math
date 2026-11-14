# Two freely inflated modules preserve the 67/69 gap

Let R and R' have outgoing rows `(2,12,17,37,75,23,47)` and
`(2,12,17,37,43,71,31)` respectively, indexed from vertex zero.
They differ by reversing the directed triangle on {4,5,6}. Mark vertices zero and one. The original pair has H=67 and H=69.

## Independent substitutions

Replace these two vertices by transitive modules of independently arbitrary orders n,m≥1. Project each Hamilton path to a word in the seven original names, retaining the order of distinct runs. Each unmarked vertex occurs once, repeated symbols can occur only for the two marks, and each consecutive pair is an oriented edge. Conversely a word with i and j runs of the two marks lifts in A_i(n) A_j(m) ways, where A_l(x)=l! S(x,l) allocates x labeled vertices among l ordered nonempty runs. Each run is then written in unique transitive order. The two marked vertices induce a single arrow, so the projected word cannot cycle using only marks, and enumeration terminates. Valid finite words are enumerated by the recurrence in replay_two_scale.py: choose any allowed next vertex, remove it from the remaining mask precisely when it is unmarked, and increment its run coordinate otherwise.

The exact nonzero table for R is

    C = [67  81  24]
        [32  86  42]
        [ 0  22  22].

All other coefficients vanish. The table for R' is C plus 2 in the first row and first column, with no other changes. Since A_1(n)=A_1(m)=1, the two absolute Hamilton path counts are exactly B(n,m), B(n,m)+2, where

    B(n,m) = (1, 2^n−2, 3^n−3·2^n+3) · C · (1, 2^m−2, 3^m−3·2^m+3)^T.

This stays physically exact even when both modules support three runs. For m=1 it is 32·2^n+{3,5}. For m=2 it is 44·3^n+72·2^n+{−47,−45}, so the independent length exposes a genuine second exponential scale, not a translation given by adding inert vertices. More generally, with arbitrary nonempty replacement tournaments A,B the same projected-word proof uses their counts of ordered full path covers into i and j nonempty paths in place of A_i(n),A_j(m), giving **difference 2H(A)H(B)**. That is a multiplicative boundary law, not a unit gap for arbitrary modules.

Both parameter directions are strictly increasing for transitive modules, because the allocations are nondecreasing and C_(2,1)=32,C_(1,2)=81 are positive. Therefore the exact inverse within the two-scale family scans only O(log N) values of n and binary-searches m after doubling to bracket N. It returns literal witness rows when equality holds. The O(log² X) parameter pairs giving counts below X cannot prove universal placement.

## Marking the last clone of vertex one yields a seam-ready pair

In the rooted port state (H,e,o), e and o sum the Hamilton counts of the complement over rooted paths of even and odd order. Mark the last clone u of vertex one in the transitive second module. Treat u as a separate singleton, with all other clones of vertex one beating u and with the same exterior arcs as vertex one. A rooted slot table records pairs (P,Q) where P is a path starting at u and Q is a possibly empty complementary path, with total vertex cover. Its key is (parity of the number of fixed vertices of P, runs of module zero in P, runs of the remaining module one in P, runs of module zero in Q, runs of the remaining module one in Q).

Compute this exact table by extending P through an available arc and choosing its stopping time. For the remaining fixed vertices count all paths Q by the same recurrence as above, now allowing the empty path precisely when no residual fixed vertices remain. The entire signed table difference consists of +2 at each key

    (1,0,0,0,0), (1,0,0,0,1), (1,0,0,1,0), (1,0,0,1,1).

All clones must then lie in Q; exactly one key has the appropriate zero-run or one-run allocation for each of the two available clone sizes (n,m−1). That allocation has weight one. The surviving P has odd order. Thus the **full** port displacement is

    Δ(H,e,o)=(2,0,2)

for all independent n,m≥1. In particular the ordinary singleton root-change seam J maps these to another absolute two-scale window with counts B(n,m)+2e and B(n,m)+2e+2.

There is an all-order limit on iterating this exact gap without changing the source arrangement. At a physical marked-block seam (THM-5070) the count row is (eta,2q_o,2q_e), where eta is H of the reader and q denotes complementary path sums ending at its entry w. The displacement stays 2 only if eta=1 and q_e=0. The reader must then be transitive (a nontransitive tournament has an odd directed triangle and hence count at least 3 by odd-cycle packing), with entry at its source (if an arc enters w then there is a two-vertex ending path, so q_e>0). Its count increase over the old source is exactly 2e.

If the new output port is neither entry nor sink, it cannot reach the reader entry, hence cannot return into the old source. A rooted path ending after choosing an odd number of later transitive vertices is even, and its complement retains the entry. Its complement count reads the same first-row state increment +2. Thus Δe>0. If the new port is the entry and not the sink, the rooted path case that uses only the entry in the block before traversing an old odd path contributes a positive term from Δo to Δe. The same holds if the reader is a singleton (its entry is not an actual sink of the composed graph, since it beats the old root). Any later seam with Δe>0 makes the count gap strictly exceed 2, because q_o≥1 and all entries are nonnegative. If instead the output is the actual sink of the composed graph, its even response e is exactly zero and the state difference returns to (2,0,2). Further exact-gap seams starting there do not raise the baseline. Thus these filters allow at most one baseline-raising exact-gap marked-block seam.

## A contrasting signed root at vertex six

Mark the fixed vertex six instead. Then the displacement is

    Δ(H,e,o) = (2,−t,−t), where t=1+(2^n−2)2^(m−1).

To verify, form the same rooted slot table using two module orders n,m and singleton root six. If a,b are the run counts for module zero in P and Q, summing over distributions of its clones gives A_(a+b)(n); inserting the parity sign for the number in P instead gives `sum_(i=0)^n (−1)^i binom(n,i) A_a(i) A_b(n−i)`. Its exponential generating function is (exp(−z)−1)^a(exp(z)−1)^b, by the ordered-bin rule. Independently apply the same calculation to the second module. After contracting, the total-response difference for n,m>0 is −2−2A_2(n)−A_2(n)A_2(m)=−2t. For the alternating response e−o, the two signed exponentials leave only the coefficient −1 at base pair (0,0), hence zero when n,m>0. The replay certifies the complete integer contractions and reproduces two independent literal root probes at each tested pair.

For n≥2 no product of physical marked-block seams, even with arbitrary block orders, can take this displacement into a mixed-sign phase direction. Indeed every matrix row (a,b,c) of THM-5070 satisfies b+c≥a. For the count row, every Hamilton path of the reader splits at the entry w to an ending path and its complement, so q_o+q_e≥eta. For a rooted row, any term of a comes from a path P' wholly in the reader. If P' avoids w, the same ending-path split in the complement gives at least 2H(remainder) in b+c. If P' contains w, split it at w and insert a formal old even path there (the e coefficient). Its tail is precisely Q2 in the seam theorem, the complement weight is unchanged, and the selected row parity is unchanged. This is injective because the split at w is recoverable. Nonnegativity then implies each row gives 2a−t(b+c)≤0 for t≥3, and nonpositive vectors stay nonpositive. For n=1, t=1: the exact small alphabet of 264 order-≤5 base-2 operators has only the two further mixed states (0,−1,1) and (−2,1,−1), with no exits to any other mixed direction. The replay checks this finite transition closure. No all-order reader claim is made for t=1.

## Reproduction
`python replay_two_scale.py 262147` checks the complete finite tensor differences, sixteen literal paired module sizes, thirty-two literal seam counts, and returns an exact 19-vertex witness at 262147. `base2_order5_matrices.txt` is the deduplicated base-2 physical alphabet from THM-5070's published generator with maximum block order five.

## A separate four-odd bank obtained by matching the two-module interface

The tensor method can find more than two offsets when the marked core changes.
Each of these cores has C12=C21=12, C22=6, and all other multi-run coefficients zero:

| C11 | two marks, in order | outgoing rows |
|---:|---|---|
| 23 | 5,4 | (8,5,9,50,7,23,63) |
| 25 | 5,2 | (2,24,3,21,37,15,63) |
| 27 | 3,2 | (2,16,3,7,45,79,31) |
| 29 | 2,0 | (2,8,3,21,39,79,31) |

The readout simplifies to H=6·2^(n+m)+C11−24. Thus
`6·2^k+{-1,1,3,5}` is an exact four-odd window for every k≥2 (use n=1,m=k−1). Unlike the 67/69 tensor, this chart happens to depend only on the sum n+m. It still gives an unbounded exact family, including the finite eleven-block census misses 197 and 3073. Twelve direct inflations provide an independent check.

A parallel matched tensor, with rows (h,54,12),(54,45,12),(12,12,4), occurs for:

| h | marks | outgoing rows |
|---:|---|---|
| 69 | 4,1 | (32,5,9,19,7,94,31) |
| 71 | 3,1 | (16,5,33,7,78,27,47) |
| 73 | 3,1 | (16,5,33,7,46,75,31) |

It gives H=4·3^(n+m)+9·2^(n+m)+h−72, hence the three-odd window
`4·3^k+9·2^k+{-3,-1,1}` for all k≥2. This closes the full mixed-base family through the finite eleven-block misses 181, 469 and 1259. The mixed-base term is an exact two-scale composition of path-allocation weights, not an extrapolation.
