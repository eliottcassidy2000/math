# Candidate H is a sieve coordinate, not a scalar cost

Take an odd target N coprime to 99. In the [literal S325 bank](../session325_absolute_cut_overlap_20261114/FULL_99_DIGIT.md), c is a digit in {1,3,…,197} with known phase r_c∈[0,4]. Using any nontrivial source T of count H>1 asks for

    N=99·2^r_c Q(T)+cH.

Reduction modulo 99 forces both c and H to be units. For a fixed proposed H, the unique odd c modulo 198 is already forced by cH≡N (mod 198). If one instead wants to sift many H, use the stronger physical gates:

    N≡cH (mod 99·2^(r_c+2)),
    N−cH is not divisible by 3·99·2^(r_c+2).

The first follows because every source of order at least two has Q divisible by four, and the second because Q is a unit modulo three whenever H is (Q≡(−1)^order H mod 3 by the odd-cycle packing formula). Impose positivity separately. These are necessary arithmetic gates, not sufficient for physical feasibility.

Take the common period 99·2^6=6336 for the first gate. There are 60 unit digits c; their phase histogram at r=0,1,2,3,4 is (20,18,10,9,3). A digit yields 2^(4−r_c) residues for H modulo 6336, because c is invertible in the local modulus. The classes belonging to different digits are disjoint already modulo 198. The total is 525 of the 1920 unit classes, fraction 35/128. Adding the mod-three gate multiplies the period by three (to 19008) and retains two of three lifts per class; each lift varies the resulting raw Q through all three residues. The linked replay confirms 1050 allowable classes among 5760 units, fraction 35/192, independent of the unit target N. This version pays the actual digit-specific phase instead of requiring five padding vertices uniformly.

The exact local example H=23 exposes why minimizing Q can be wrong even at fixed H and order. The six-vertex masks (2,12,9,17,39,15) have (H,Q)=(23,272), and the active one-apex word 11000 has the same H and order with Q=296. At N=30017, the smaller Q cannot work under any digit or nonnegative phase: modulo 23 the phased value 272·2^r would have to equal 20, whereas phase zero has residue 272≡19 and every positive phase makes 99·272·2^r>N. The larger Q=296 has residue 20 and furnishes the exact digit 31 at phase zero.

This is the precise least-resistance distinction: minimize over feasible routes of the full pair, not the first available scalar. The sieve moves the expensive check to a smaller pool; the arbitrary-Q decoder supplies literal sources when its grammar succeeds. A proof of coverage still requires that one valid H class contain an admissible absolute profile for each remaining target.
