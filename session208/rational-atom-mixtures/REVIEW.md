# Independent spectrum audit

Status: PASS, 2026-11-06.

Reviewed PROOF.md and verify.py. Independently enumerated all vertex permutations of all fifteen retained tournament witnesses, confirming every claimed Hamiltonian-path count, and checked each rational dual inequality as the exact integer inequality 16a+24b<=11c. Replayed the supplied verifier, which checks strong connectivity and subset-DP counts as well as the residue DP. Confirmed the inherited census is byte-identical to both session207-size-spectrum/census.json and session207/inverse-prefix661/census.json. Exhaustiveness remains inherited from that prior census, not newly exhaustively enumerated here.

The 375^3*405^2=15^11 identity and forty-vertex dominance-join witness are correct. The dual vector (16/11,24/11) and primal weights (3/11,2/11) certify finite-dictionary optimum 40/11; zero slacks occur only at 375 and 405, proving uniqueness of the rational mixture. Support size two means two atom types, while the integer construction has five SCCs.

The general support-compression proof, bounded additive integrality gap, exact denominator multiples and infimum-over-all-finite-dictionaries characterization are valid. The numerical lower bound applies only to the fixed strong-order-at-most-eight dictionary. Global conclusions are upper bounds, not global optimality of forty vertices or global optimality of the rate 40/11. No use of the attached preprint's new logarithmic theorem is necessary.

Independent checker and output: independent_check.py and independent.out. This review does not claim a direct 2^40 path-DP computation: the forty-vertex count follows from the proved SCC product rule and independently checked eight-vertex factors.
