# The pairing carrier must change at orders divisible by four

Status: the all-order arithmetic branch integrality is PROVED in `BRANCH_INTEGRALITY.md`; the direct matching-flip rule already fails at n = 4 by independent FINITE-EXACT enumeration; a fixed-point-free involution on all unadorned even-Aut types is impossible at n = 8 using the inherited FINITE-EXACT census. No claim refutes every possible decorated or weighted combinatorial pairing.

## Smallest hostile to the exact Theorem B operation

Theorem B in THM-4531 uses the operation

    q -> q + M_g,

where g is an involution in the automorphism group of q and M_g is its perfect matching. At n = 4 there are two switching-class isomorphism types. Their automorphism groups have orders 4 and 12, with odd-element proportions 1/4 and 3/4. For both types, and for every involution in their automorphism group, q + M_g is isomorphic to q. Thus the matching-flip operation is the identity on the two types, despite N_lev(4) = 1.

Independent exact encodings: normalize the tournament so that vertex 0 is a source, then encode the three remaining edges (1,2), (1,3), (2,3) in low-bit-first order. The two type representatives are masks 0 and 2. For mask 0, take g = (0 2)(1 3); the flip produces mask 7, in the same type, and h = (1 3) transports mask 0 to mask 7 while commuting with g. For mask 2, take g = (0 1)(2 3); the flip produces mask 5, in the same type, and h = (2 3) transports mask 2 to mask 5 while commuting with g. All three involutions for mask 2 give the same type.

The positive control n = 6 reproduces the existing Theorem B: the four even-Aut types have rho = 1/2, and the matching flip pairs canonical masks 0 with 40, and 4 with 41. The entire n = 4 and n = 6 universes were enumerated directly over the labelled switching quotient and all vertex permutations, independently of nauty and the inherited scripts. Both ordinary and Python -O runs give identical receipts.

## Stronger obstruction to any unadorned fixed-point-free involution

The inherited complete n = 8 census in section 5 of

`05-knowledge/results/procgen_tbij_20261001_natural_bijections.md`

(and the audited summary in THM-4531) gives the following 27 types with even-order automorphism group:

- 19 with Sylow 2-subgroup C_2, all rho = 1/2;
- 2 with Sylow C_4, both rho = 1/4;
- 2 with Sylow V_4, having rho = 1/4 and 3/4;
- 2 with Sylow C_8, both rho = 1/8;
- 2 with Sylow D_8, having rho = 1/8 and 5/8.

Their total number is odd, so no fixed-point-free involution on this entire unadorned set exists, even if it mixes Sylow types. Furthermore

    N_lev(8) = 19(1/2) + 3(3/4) + 1/4 + 3(7/8) + 3/8 = 15,

whereas half the number of these types is 27/2. The new `check_matching_pairing.py` verifies this arithmetic explicitly, but the n = 8 group census itself remains inherited, not newly enumerated here.

The total branch integrality is now proved uniformly. An interpretation by a pairing would need a different carrier, such as suitably decorated objects or a more general grouping with weights. The statement “by a combinatorial pairing” in HYP-9172(b) should not be marked PROVED by the arithmetic argument, nor declared universally false merely because these specific carriers fail.
