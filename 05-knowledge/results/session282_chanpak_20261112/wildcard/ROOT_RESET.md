# A parametric escape from all powers of 13

Status: PROVED, independently reviewed; exact finite-state proof with freshly replayed certificate. This addendum does not modify the frozen main proof.

Let A be the six-vertex rooted tournament in PROOF.md, with deletion D=(18,8,3,5,14), root outgoing cut9. Form A^k by physical rooted fusion, and let r be vertex1 inside the final A block (global label 1+5(k-1)). Write V_r for the ending mass at r. Then, for every k>=4,

    V_r(A^k) = 2*3^(k-4) modulo13.

Consequently the physical triangle extension at r has count

    H(A^k)+2V_r(A^k) = 4*3^(k-4) modulo13,

which is never zero. Since 13^(k-2) divides H(A^k), one actual reroot-and-triangle operation can remove arbitrarily many powers of13. The old two-jump obstruction does NOT imply a bounded loss under arbitrary rerooting. This is not a universal target selector or a renewable fixed receiver.

## Exact three-strand compiler

Retain the original common root z. By definition V_r counts a Hamiltonian path ending at r together with a disjoint complementary path (possibly empty). Classify the position of z:

1. z lies in the path ending r. Write that path P,z,Q, with Q ending r; call the complementary path R.
2. z lies in the complementary path. Write that path P,z,Q, and let R be the path ending r.

P and Q can be empty except that Q is nonempty in case1. R is nonempty in case2. A nonempty P must end at an inneighbor of z; a nonempty Q must begin at an outneighbor of z. Removing z leaves the ordinal chain D^k. Each of P,Q,R traverses its occupied blocks in increasing order, so its local pieces independently concatenate. This is a bijection, not a heuristic path splitting.

The state records P's last-vertex condition (empty, good, bad), Q's first-vertex condition (empty, good, bad), and whether R is nonempty. Thus there are18 states. For each block distribute its five vertices among the three strands. For each subset count local Hamiltonian paths, distinguishing the relevant first/last endpoint. The update replaces P's status on a nonempty local P, preserves Q's first nonempty status, and ORs R's occupancy. This is the explicit integer linear operator F in reroot_transfer.py. Its 243 assignments per block are exhaustive. In the final block force r into Q or R and count paths ending at r; adding those two final linear functionals yields V_r. Empty intermediate strands are allowed exactly as above.

Let s0 be the empty state vector, and sm=F^m s0. The checked certificate is the full18-state identity

    F s3 = 3 s3 modulo13.

Linearity proves sm=3^(m-3)s3 modulo13 for every m>=3. The final functional L satisfies L(s3)=2 modulo13. This proves the asserted formula with m=k-1, without extrapolation from a finite count census. REROOT_RESULTS.json records the complete s3 vector. The standalone checker constructs F from literal tournament masks and tests the full vector identity, including zero states, not merely the projected V sequence.

As controls the same compiler is compared with literal subset-DP V_r on A,A²,A³ (orders6,11,16). Normal and optimized Python replay outputs agree. Ordered three-path-cover count of D is471≡3 modulo13, whereas its one- and two-path totals13 and104 vanish modulo13. The third-strand channel explains the eigenvalue and why the original two-jump estimate does not extend to arbitrary rerooting.

## What is needed for iteration

The output count and V_r are not enough to iterate a receiver-preserving construction. One needs the actual rerooted deletion count h_r, the full packet (H,h,Q,U,V,d,e), and an identified root to which the next operation applies. One then needs a uniform rule returning these data to a reusable family while controlling the absolute target and selector range. The current theorem proves an all-order arithmetic escape, not such a return map. Enlarging the compiler to carry those marked endpoint coordinates is the next precise obligation. A single scalar residue or an unlimited supply of escaped counts does not meet it.
