---
id: THM-4850
title: "Ninja (weak B) pairs of the square-sum graph exist exactly for n >= 42: the gate is the unique Hamiltonian path of Q_41 from 1,3 to 10; session124's 49 digit operators work on the new pairs 42..52, so the base-49 ninja numeration n = 49^k b + sum d_i 49^i (b in [42, 2081], d_i in [24, 72], the owner's 49, 24, 72) starts at 42"
status: PROVED (exhaustive enumeration for n <= 53 + cited session124 certified base 53..2620 and recursion) + FINITE-EXACT (539 operator replays verified by session124's own validator); the reading 'ninja pair = weak B pair' is INTERPRETATION of the owner's term
session: opus-2026-10-31-S69
date: 2026-10-31
note: 05-knowledge/results/ninja_pairs_square_sum_number_system_20261031.md
script: 04-computation/experiments/ninja_pairs_square_sum_20261031.c, ninja_square_sum_number_system_20261031.py (+ .out), ninja_numeration_replay_20261031.py
related: [THM-4851, HYP-9447, THM-4505, THM-4510]
external:
  - "session124_recursive_square_prefix_20261031 (weak B pairs, 98 paired digit operators, recursion n = 49m + d, d in [24,72]); square_sum_B_base_interval_20261019 (certified base 53..2620)"
---

# THM-4850 - the ninja numbers start at 42

**Definition** (session124's weak B pair; the owner's "ninja pair"). At n it consists of:
- a Hamiltonian path S of the square-sum graph Q_n (Shorter) and L of Q_(n+1) (Longer);
- both start 1, 3;
- equal positional parity on every common vertex;
- the even-order member ends at 8 and the odd-order member at 10.

**(1) Threshold.** A ninja pair exists iff n >= 42.
- For n <= 53, by exhaustive enumeration: none for n <= 41; for 42 <= n <= 53 the numbers of compatible Shorter paths
  are 10, 4, 45, 3, 5, 11, 859, 108, 153, 230, 10710, 32242.
- For 53 <= n <= 2620, session124's certified base.
- For n >= 2621, session124's recursion B(49m + d) <- B(m).

**(2) Gate.** Odd-order Hamiltonian paths 1,3 -> 10 exist, among odd orders <= 41, only at order 41, where there is
exactly one:
1,3,6,30,19,17,32,4,21,15,34,2,7,18,31,33,16,9,40,41,8,28,36,13,23,26,38,11,5,20,29,35,14,22,27,37,12,24,25,39,10.
It is compatible with no partner at n = 40 or n = 41. Even-order paths 1,3 -> 8 first exist at order 34. So the
threshold 42 is set by the odd member.

**(3) Numeration.**
- Every n >= 42 has exactly one expansion n = 49^k b + sum_{i<k} d_i 49^i with b in [42, 2081] and d_i in [24, 72].
- Adding a node is the base-49 odometer, with 72 + 1 carrying to 24.
- The expansion is constructive: session124's paired operators B_p{m mod 2}_d{d}, applied to one ninja pair at each
  m = 42..52 and every d in [24, 72], give 539 outputs at n = 2082..2620. All pass session124's own `check_pair`.
- Session124's thresholds T_l = (107 * 49^l - 1)/2 are the numerals [53; 24, ..., 24].

**(4) Why consecutive orders.** The operator v -> 49v + (-1)^j c preserves square sums only because 49 is a square.
It draws low digits from the Shorter and high digits from the Longer (switch at 49 + (-1)^(n+1) c <= d). Together the
two consecutive orders cover a full residue system mod 49.

**Novelty.** The threshold 42, the gate, and the start of the numeration at 42 (the base block extended from 53 down
to 42 and replayed) are new. The definitions, the operators and the large-n existence are session124's.
