# Ordinal-sum ear composition and a reverse-connector obstruction

Status: exact all-order identity, an explicit family, and a limited obstruction to one proposed gadget. It does not characterize tournament ear carriers or settle the path-count spectrum.

Let T=A→B be the ordinal sum of two nonempty tournaments: every A vertex dominates every B vertex. Put a=|A|, b=|B|, and H=H_A H_B. For S⊆A⊔B let S_A=S∩A and S_B=S∩B. Let E_A^+ (respectively E_A^−) count Hamiltonian paths of A ending inside (outside) S_A. Let P_B^+ (respectively P_B^−) count Hamiltonian paths of B starting inside (outside) S_B. Let X_A count ordered decompositions A=A1⊔A2 into directed paths, A1 possibly empty and A2 nonempty, with the first vertex of A2 inside S_A. Let Y_B count ordered decompositions B=B1⊔B2 into directed paths, B1 nonempty and B2 possibly empty, with the last vertex of B1 outside S_B. Empty paths count once.

Then the exact child count is

 F_T(S) = F_A(S_A) H_B + H_A F_B(S_B)
          − E_A^− P_B^− − E_A^+ P_B^+ − E_A^− P_B^+
          + X_A Y_B.                                             (1)

Proof. Every child path has at most one B→A transition, since the only way to make such a transition is through x. If it has one, it uniquely has word form A1 B1 x A2 B2, where B1 and A2 are nonempty and A1,B2 may be empty. The two A paths and two B paths are independent and determine x-adjacencies, giving X_A Y_B. Otherwise the A/B word is A* B*, with x located in the A portion or B portion. The two respective raw counts are F_A H_B and H_A F_B. In the former, an invalid A→B boundary occurs exactly when x is last in the A+x path and the first B vertex is outside S_B; subtract E_A^− P_B^−. In the latter, an invalid boundary occurs exactly when x is first in the B+x path and the last A vertex lies inside S_A; subtract E_A^+ P_B^+. The corrected classes overlap exactly in paths A x B, counted E_A^− P_B^+. Removing that overlap gives (1). Note that neither F_A H_B nor H_A F_B by itself may be included uncorrected: the adjacency to x at the block boundary matters.

## Full reverse connector

Take S=A, so B→x→A. Then the three boundary corrections vanish and F_A(A)=H_A, F_B(∅)=H_B. Counting marked breaks in permutations, X_A=a H_A+B_A and Y_B=b H_B+B_B, where B_U denotes the number of permutations of U with exactly one backward adjacent arc. Therefore

 F_T(A)=2H_A H_B + (a H_A+B_A)(b H_B+B_B) ≥ (ab+2)H.  (2)

The marked-break equality holds because a valid permutation contributes a−1 nonempty breaks and the empty first part once, while a permutation with exactly one backward adjacency contributes its unique bad break. The argument for B is symmetric. Consequently, in a proposed bank trying to realize 3H plus an offset bounded independently of H, the fully reversed connector cannot work along families with ab≥2 and unbounded H: its excess above 3H is at least (ab−1)H. This excludes only this specific orientation and block construction, not all ears of decomposable tournaments or other gadgets.

When A and B are transitive, H_A=H_B=1 and B_U=2^|U|−|U|−1 (the Eulerian one-descent count). Thus F_T(A)=2+(2^a−1)(2^b−1). Equivalently, the factor counts follow directly because each nonempty subset determines the second path uniquely, so X_A=2^a−1, and each nonempty subset determines the first B path uniquely. In particular a=b=1 gives 3, while (a,b)=(1,2) gives 5.

## Audit

`python3 check.py` independently enumerates permutations for Hamiltonian paths and decompositions and uses subset dynamic programming on literal child tournaments. It exhaustively checks (1) for all labeled blocks of sizes (1,1), (1,2), (2,2), (2,3), (3,3), and (1,4), and verifies (2) in the full-reverse case. It is a finite audit, not the proof. `out.txt` records the result.
