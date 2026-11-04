# Independent audit of modulo-eight module transfer

**PASS. No mathematical correction found.** The even-order unit multiplier, the physical-port lift, and the odd-order Hamilton-cycle correction are all valid. This review uses the primary odd-cycle expansion already checked in arXiv:2307.05569v1, Theorem 1.39, which also implies the cited modulo-four refinement.

## Profile congruences and boundaries

For an unordered k-partition, the product of part Hamilton counts is one plus twice the sum of odd directed-cycle counts inside individual parts, modulo four. A fixed odd cycle contracts to one atom, leaving exactly S(N-|gamma|+1,k) partitions. This justifies the displayed partition congruence without assuming path-cover log-concavity.

For k=2, the relevant Stirling number is odd precisely when the contracted ground set has at least two atoms. The excluded case is a spanning odd Hamilton cycle, explaining the subtraction of z. Multiplying the modulo-four count by 2 gives the stated modulo-eight q2 formula. For k=3, the contracted size has N's parity. S(t,3) is even for even t, odd for odd t>=3, and zero at t=1. Multiplication by 6 preserves the needed modulus because an omitted multiple of four becomes a multiple of 24. The same spanning-cycle exclusion produces the odd-order correction. For k>=4, the ordered-slot factorial is divisible by eight. The conventions at N=1,2,3 check directly.

An equivalent verification from the full odd-cycle-cover expansion is q_k=k! sum_C 2^c S(ell,k), where ell is the number of permutation cycles including singleton cycles. At even N, ell is even and at least two. Every nonidentity term of q2-2H vanishes modulo eight, leaving 2^N-4; every nonidentity q3 term vanishes modulo eight because S(ell,3) is even. This independently confirms the even formulas and their exact scope.

## Arbitrary substitution and actual ports

For same-size even modules, the higher differences reduce to Delta q2=2 Delta H and Delta q_k=0 for k>=3, modulo eight. The exact nonnegative maximal-run expansion therefore gives multiplier c1+2c2. The identity c1+2c2=H(Q[TT2]) is exact, since TT2 has profile (1,2). This count is odd, so multiplication preserves and reflects a zero modulo-eight gap. It need not preserve a residue two numerically; it always preserves a residue four.

For a physical local two-port query, both new ports receive the old module's same incidence to each exterior name. Consequently old module plus ports is a common module in both actual extended tournaments, with the same rooted exterior quotient. Its order is old order plus two, so evenness is preserved. The theorem acts on the literal query counts and does not drop any cubic or quartic information. It lifts the odd-pair quotient theorem through arbitrary exterior contexts without changing that theorem's local hypotheses.

For odd module order, the correction term -4(c2+c3)Delta z is necessary. The first multiplier remains odd; an even Delta z restores unit transfer. The retained nine-name port-core controls have odd cycle differences -13 and -15, with local gaps -68 and -84 and actual ambient gaps -600 and -720. The latter are zero modulo eight but not zero as integers. These are modular obstructions, not exact response collisions.

## Independent execution

The standalone audit imports no producer helper. It independently computes all complete profiles through order five (1,099 named sources), counts 1,033 odd-order Hamilton-cycle controls by literal permutations with the starting vertex fixed to remove rotations, and recomputes all 570 declared even-module substitutions. It also recomputes the complete profiles and actual strong-context counts for all four physical port-core controls.

Both unchanged producer runs agree with their saved receipt. Both independent normal/optimized outputs agree. `receipt.json` pins the source proof, source verifier and both transitive helper dependencies, including `response/verify.py`. Producer stdout is represented by its SHA256 rather than duplicated, since the source packet already retains both full receipts. All checks use explicit exceptions.

## Separate modulo-sixteen sharpness control

The added `SHARPNESS.md` boundary is also correct. Independent endpoint/profile dynamic programming reproduces the complete six-name profiles (11,178,828,1752,1800,720) and (3,90,612,1608,1800,720), and literal ambient counts 211 and 99. The core gap eight is nonzero modulo sixteen, while the ambient gap 112 vanishes modulo sixteen and remains nonzero as an integer. This refutes a universal even-order odd-unit transfer modulo sixteen, not a full-response or exact scalar collision. `sharpness_audit.py` replays the producer's independent permutation/cut proof controls normally and optimized; `sharpness.normal.json` and `sharpness.optimized.json` agree and bind all dependencies.
