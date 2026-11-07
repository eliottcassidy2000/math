---
id: HYP-9494
title: "Root-change words and Carmichael: (A) the K = 0 word grammar has (c + o(1)) sqrt(X) singleton values, c about 0.50 (4,083 of the 8,191 chamber words below 2^26); (B) K = 1 has infinitely many singletons but O(log X) of them; (C) for K = 2 no value N > 1 has exactly one word (for larger K singletons can sit at values of exact depth K: M_7(1897) = 1, the only word-singleton above 1 with K >= 2 below 10^9); (D) two root changes from the rooted edge cover every odd N > 3,780,967, and every odd N > 2647 is a strong C = 1 Hamiltonian-path count of root-change depth at most 7, the maximum 7 attained only at 1877 and 1897; (E) the mean two-root-change multiplicity has a 3-adic local factor rising with v3(N) towards about 3"
status: "(C) FINITE-EXACT through 10^9 for K = 2 (pre-registered P3 CONFIRMED, 6d0e88f0e8); the first draft's form for all K >= 2 was REFUTED by M_7(1897) = 1 (audit M4, MISTAKE-733). (D) FINITE-EXACT through 10^9 (P1 CONFIRMED: F_<=2 misses exactly 140 odd values below 10^9, largest 3,780,967; depth table exact below 10^9; the C = 1 spectrum misses exactly Session197's 29 values below 10^9). (A), (B) FINITE-EXACT counts to 2^26, OPEN beyond. (E) EXPLORATORY (exact 16-bit census to 2^24). P2 and P4 (window minima) FAILED and are recorded. AUDIT: independent, PASS WITH CORRECTIONS (every number reproduced with independent code, including the full 10^9 census; MISTAKE-733)."
session: kind-pasteur-2026-11-07-S14
date: 2026-11-07
related: [THM-4904, HYP-9480, THM-4888, THM-4884, HYP-9477, session200/carmichael-fibers/PROOF.md, session199/physical-adjoint/PROOF.md, session199/hidden-structure/PROOF.md, session198/c1-monotone/PROOF.md, session197/lowcycle/C1_COMPLETE.md, 05-knowledge/results/session218_smith_repunit_totient_20261107/README.md]
script: 05-knowledge/results/root_change_carmichael_chamber_20261107/rc_census.c; 05-knowledge/results/root_change_carmichael_chamber_20261107/rc_depth.c; 05-knowledge/results/root_change_carmichael_chamber_20261107/evaluate_prereg.py; 05-knowledge/results/root_change_carmichael_chamber_20261107/literal_witnesses.py; 05-knowledge/results/root_change_carmichael_chamber_20261107/audit/census_audit.c
---

# HYP-9494 - the Carmichael hierarchy of root-change words, two root changes to 10^9, and root-change depth

Packet: [README](../results/root_change_carmichael_chamber_20261107/README.md),
[PROOF](../results/root_change_carmichael_chamber_20261107/PROOF.md) (sections 5-7),
[AUDIT](../results/root_change_carmichael_chamber_20261107/audit/AUDIT.md). The proved structure is THM-4904.

## (A)-(C) The singleton hierarchy

M_K(N) is the number of O/Z/J words with at most K root changes and value N (THM-4904). Singleton values N > 1, per
dyadic window:
- K = 0: 1, 1, 0, 1, 2, 0, 1, 2, 3, 5, 8, 9, 10, 19, 30, 41, 61, 62, 98, 147, 234, 297, 418, 597, 854, 1183 (to 2^26,
  the first entry being N = 1). The ratio per doubling tends to sqrt(2). The total above 1 is 4,083 = 0.498 sqrt(X)
  below 2^26, out of 8,191 chamber words. That share, 0.4985, matches the root's Poisson miss rate of HYP-9480 (A).
- K = 1: about 3 per window to 2^26, with no growth; the first are 49, 161, 2457, 12553, 16777, ...
- K = 2: none above 1 through 10^9.
- K >= 3: a singleton must have root-change depth exactly K, because then M_K = B_K has no doubled term. Below 10^9
  the only one is M_7(1897) = 1, with word ZZJJJJJJJOOZ, while M_7(1877) = 4.

Heuristic: a chamber word is alone with probability exp(-lambda_K(N)). Here lambda_0 is bounded, lambda_1 is about
(1/2) ln N and lambda_2 is about c ln^2 N.

Proof routes:
- (A) needs an anti-concentration estimate for the root tail on the chamber.
- (C) needs either a second involution on the K = 2 chamber or an effective lower bound on M_2 outside a finite set.
  THM-4904 (3) already confines a K = 2 singleton to B J C C^dagger J B^dagger, with B ending in Z and C empty or
  beginning with O.

## (D) Two root changes and root-change depth

- F_<=2 misses exactly 140 odd N below 10^9 (largest 3,780,967), extending HYP-9480 (D) from 10^8.
- F_<=3 misses 44, the largest 27,781.
- With seven root changes (counted in source words, without the closing J) the misses below 10^9 are exactly
  Session197's 29 C = 1 omissions, and the depth table is:

| Root-change depth | Values below 10^9 |
| --- | --- |
| 3 | 96 values |
| 4 | 373, 1811, 2149, 2999, 4411, 4547, 7129, 17383, 19951, 22091, 27781 |
| 5 | 667, 3523 |
| 6 | none |
| 7 | 1877, 1897 (literal 14-vertex uniquely Hamiltonian witnesses) |

Session199 hidden-structure (section 6) has the analogous minimal-J-cost table for the balanced sub-readout to
1,000,001.

Conjecture: two root changes suffice for every odd N > 3,780,967, so the global maximal depth is 7. An effective N_0
for F_<=2, together with this census below N_0, would prove that every odd N > 2647 is the Hamiltonian-path count of a
uniquely Hamiltonian strong tournament. With the known complete spectrum below 2647, where every odd value except 7
and 21 is attained by some tournament, that would prove the Chan-Pak conjecture.

## (E) The 3-adic local factor (exploratory)

The mean of M_2 on [2^23, 2^24) by v3(N) = 0, ..., 5, >= 6 is 205.8, 332.0, 457.3, 510.7, 556.8, 614.2, 618.7. Among
3 ∤ N, the class N = 1 (mod 3) is lowest (199.96 against 211.73).

This is consistent with the misses and near-misses lying in N = 1 (mod 3): 14 of the 15 misses above 10^5, and 104 of
the 130 values with at most 12 words in (2^22, 2^28). No near-miss is divisible by 3. A 6% gap in the class means does
not by itself explain a 104:26 split.

The pre-registered window-minimum predictions P2 and P4 failed: extremes of a heavy-tailed statistic are noisy.

Mechanism, an illustration rather than a derivation: seeds with 3 | gcd(R, d) are confined mod 3 (Session206). The
dilated root JJ(t) = (5,2,6,3) gives H = 3 H_t(m) + 2 - 2(m mod 2), which is never 1 (mod 3).

## Tests

- Extend the F_<=2 census past 10^9, by segments (the audit's `census_audit.c` already segments).
- Derive the 3-adic factor from Session206's weighted capacities restricted to seeds with 3 | gcd(R, d).
- Prove that the K = 2 chamber has no singleton above a computable bound.
