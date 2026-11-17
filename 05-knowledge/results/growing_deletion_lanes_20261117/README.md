# Change the deletion core, not just the tape: five ready phases and a growing sixth

**Result.** Larger rooted deletion cores give permanent unit switches at 7p in five prime classes modulo 18 while keeping the previously defined normalized ready state on one branch. A separate nonopaque construction grows the deletion count without limit and supplies the sixth phase, p≡1, at the literal pair **19535/19537**, where 19537=7·2791. This last construction preserves an exact unit gap but does not satisfy readiness. No universal exact-target theorem follows.

## 1. A same-count re-encoding that actually matters

Use the [previous ready-core compression](../ready_core_growth_and_mod18_20261117/README.md). At a marked root let (H,g,E,O) denote its count, its deletion count, and its weighted even/odd start data. A normalized root has O−E=g; write d=H−2E+g. A ready member has d=0. A companion has d=±2 with the same (g,E,O); a common binary tail at the retained root always preserves the count difference.

The count 91=7·13 is ready at g=5 but **cannot** have a matched unit there. Moving to g=9 gives ready count 91 with a matched 89 companion. This is a literal example of the [landing and re-encoding proposal](../chanpak_selection_retrospective_20261117/README.md): holding the scalar fixed, change its deletion core to obtain a continuation the lower state forbids. It is not a claim that one can get between the two realizations by a single edit.

A complete classification in this architecture at g≤73 gives the following table. For each listed class the displayed g is the first that admits a ready count with a matched unit and an actual prime example in that class. At p≡1 there is no such normalized ready pair anywhere in this deletion range.

| p mod 18 | g | prime p | ready 7p | matched neighbor |
|---:|---:|---:|---:|---:|
| 5 | 9 | 113 | 791 | 789 |
| 7 | 11 | 547 | 3829 | 3831 |
| 11 | 5 | 173 | 1211 | 1213 |
| 13 | 9 | 13 | 91 | 89 |
| 17 | 37 | 197 | 1379 | 1377 |
| 1 | none through 73 | | | |

**How completeness is checked.** For normalized d∈{0,−2,+2}, the positive even-path count K satisfies K=g−d/2. After harmless sources and tail powers of two are removed, retain at most floor(log₂(g+1))+4 singleton components. The remaining strong components have counts multiplying to g. The [strong floor](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) forces each component to have at most eight vertices for g≤73 (order nine has minimum 75). Enumerate all ordered factorizations and all locations of the retained singleton components, choose all root masks, and compare slopes via the odd part A of 2E. Near-sinks double the slope, so any two seeds with the same (g,A) can be aligned at the greater starting exponent. Testing the six exponent phases modulo 126 then reads the quotient p modulo 18. The g=1 and g=3 exclusions are analytic in the earlier packet, and 7 and 21 cannot occur as deletion counts. The lone would-be small lane A=g=5 always has count divisible by 5 with count >35, and is excluded from the prime test.

**The speedup from the older material.** Instead of rebuilding a subset-path matrix on the whole deletion for each ordinal arrangement, compute each strong component's weighted start/end menus and even/odd optional path-cover totals once. A path crossing components uses a start menu at its first component, an end menu at its last, optional menus between, and ordinary Hamiltonian counts outside. Convolve in the two parity coordinates. Use K's positivity to prune masks. This is the same weighted-path arithmetic as the old [target-pruned continuation search](../fresh_unit_prime_gluing_20261116/target_permuted/README.md), adapted to an exhaustive bank of deletion cores, not a search over arbitrary full tournaments at large H. The direct older g5 traversal gives the same rays; phases.py independently verifies literal subset counts for all five witnesses. Integer safety on the claimed range follows since H(S)H(D−S)≤H(D) by order-preserving path merge, and the deletion contains at most 21 vertices: the optional cover total is bounded by 73·2^21, so even twice that fits the integer implementation.

## 2. Growing the core opaquely: when it works and when it cannot

The old [terminal-block transfer](../fresh_unit_prime_gluing_20261116/target_permuted/README.md) gives an exact comparison. Append a nonempty block M beaten by all old nonroots and beating the root. Put u=H(M), R=P_M(2) (two ordered path slots, each allowed empty). At a normalized source the transfer is

    g' = ug,    E' = R E,    O' = R E + ug,
    H' = uH + 2E(R−u),       d' = u d.

Indeed the optional path slots split equally by parity. Weighted root paths other than the root alone operate on those slots; positive even paths through a new return contribute the indicated term. This both explains the mobility and identifies the obstruction: d=0 grows without difficulty, but the ±2 defect scales by u. Even choosing a different opaque module for each branch cannot preserve a nonzero unit on a companion whose module has u>1.

One genuinely useful old/new combination is the four-vertex strong tournament (2,4,9,3), for which (u,R)=(5,32). It has the **same two-cover weight as a five-vertex transitive block**, but different Hamiltonian count. In a fixed number of such slots, swapping one transitive tile for this one multiplies g by 5 while keeping 2E fixed; on the ready locus H+g is then unchanged. This grows deletion complexity without adding a third run: removing the sole return root leaves the appended ordinal region, and each side of the root visits it at most once.

A scope check on the old A=19 lane: with k of the H5 modules and any added neutral powers of two, its ready read is H=19·2^s−5^(k+1). To divide H by 7, k must be even and s≡k mod 3. To avoid a factor of 3, s must be even. Checking k mod 6 now gives H≡5 mod 18 in every case, exactly the old p≡11 lane. So this elegant growth changes the core but neither preserves the unit nor unlocks a new phase on that particular ray. Other source slopes can change phases, but the unit penalty remains.

## 3. Nonopaque growth that really preserves the unit

Take upstream A with a marked vertex v, deletion count 1 and a normalized ready profile. Use the complete g=1 rays from the previous packet: h=H(A)=a·2^t−1, a=3 or 5, t≥1. At v its even/odd weighted start counts are (h+1)/2 and (h+3)/2. Build A→B, then add one new global root r. In A the new root beats **just v**; in B its outcontact set is U. Change the physical B tile between these rows and masks:

| B rows | U |
|---|---:|
| 126,60,120,112,96,0,34 | 9 |
| 62,124,120,112,32,0,49 | 70 |

In both B tiles, H(B)=9, R_B=P_B(2)=336. The weighted even/odd start counts from U are (100,100). Count weighted ends in the incontact complement of U by parity, denoting these by T_e,T_o, and let J_B count weighted even paths starting in U and ending outside it. The triples (T_e,T_o,J_B) are (143,143,91) and (142,144,91). These are independent literal subset calculations.

By decomposing a root-start path, or an old even path closing a cycle through r, at the single A→B boundary, the new profiles are

    H = 486h+581 or 486h+583,
    g = 9h,    E = 268h+336,    O = E+g.

Indeed all old paths starting at v contribute their total even-plus-odd start weight h+2 times the optional B read R_B; paths starting inside B see multiplier h. Parity balances, so E=(h+2)R_B/2+100h in either variant. The full count is

    H = 9h + 9(h+1) + 2h J_B + (h+1)T_e + (h+3)T_o.

The constant difference is two while all coefficients of h agree. This uses opposite endpoint imbalances in **different** B tiles, rather than reading the same module blindly and multiplying its defect.

Both profiles permit every common binary sink/near-sink continuation at r, preserving +2. The tail can be visited at most twice. Even the growing upstream transitive tape has at most two visits forced by the ports: it can return from the exterior only by consuming v (if it exits through r, the only entrance back into A from r is v). Thus growing the deletion core here does not create the old third freely assignable run.

At a=5,t=3 the two counts are **19535,19537** and g=351. The upper endpoint is 7×2791 with prime 2791≡1 mod 18, the missing phase of the bounded ready bank. Generally its upper endpoint is 486a·2^t+97; once divisible by 7 its quotient is 1 mod 18. Neither endpoint is normalized ready for the previously specified defect: d=−41h−91 or −41h−89. This is a landing with a permanent unit repair, not a proof that its exact count has a ready lift in an unknown larger core. The exact target inverse in growth.py combines the upstream choice with the binary split and reports both hits and misses. The 19537 hit does not extrapolate to all primes in that phase.

**A precise next tile test.** Use instead the near-transitive Fermat upstream with h=2^m+1, g_v=1, weighted starts ((h−1)/2,(h+1)/2). Let T'_e,T'_o be weighted end counts *in* U (not the complementary menu above), K_B count even paths with both endpoints in U, and S the common even/odd start sum in U. For both output branches to keep a fixed defect as h varies, the exact equation is 2K_B+T'_e+T'_o=2H(B), and then d=−(T'_o−T'_e). This follows by applying the positive K test to the A→B decomposition. Since each specified vertex's unsigned end menu is at least H(B) (cut each full path at that vertex), and every two vertices of U already contribute a positive even path, the equation allows **only a singleton U**. A ready branch and a −2 companion would therefore need two B tiles with the same H(B) and the same combined read R_B+2S, with a singleton of total end weight 2H(B) in each, but signed end weights 0 and 2 respectively. Starting from the earlier endpoint-menu matches produced the growing repair; this stricter singleton condition is the concrete next test for turning it into a growing *ready* repair.

## Reproduce

Convert the order-eight class codes of the [archived n7/n8 data](https://github.com/eliottcassidy2000/math/tree/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities) to rows with the [previous ledger command](../chanpak_selection_retrospective_20261117/README.md#reproduction-of-the-diagnostic). Then run:

    g++ -O3 -std=c++17 ordinal.cpp -o /tmp/core_ordinal
    /tmp/core_ordinal n7_representatives.txt n8_rows.txt 5 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 > cores.out
    python phases.py cores.out
    python growth.py

The class-representative completeness is inherited. The phase witnesses and growing pair replay use direct literal counts independent of the ordinal parity compiler.
