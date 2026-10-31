---
id: THM-4850
title: "Weak B pairs of the square-sum graph (session124's normalisation of Gerbicz's nice pairs; the owner's 'ninja pairs') exist exactly for n >= 42, sharply: none for n <= 41, the gate being the unique Q_41 path 1,3 -> 10, parity-incompatible with all partners at n = 40, 41; Gerbicz-normalised nice pairs start at 41 and free nice pairs at {15,16,25,26,27} u [29, oo); the base-49 numeration (radix 49, digits 24..72) runs from the base block [42, 2081] (prior work: session124 threshold42)"
status: PROVED (exhaustive enumeration for n <= 53 + cited session124 base 53..2620, threshold42 bases 42..52 and recursion) + FINITE-EXACT (counts; 539 operator replays verified by session124's own validator) + CITED; the reading 'ninja pair' = Gerbicz's 'nice pair' (weak B normalisation) is INTERPRETATION; independently audited SOUND WITH CORRECTIONS (MISTAKE-694)
session: opus-2026-10-31-S69
date: 2026-10-31
note: 05-knowledge/results/ninja_pairs_square_sum_number_system_20261031.md
script: 04-computation/experiments/ninja_pairs_square_sum_20261031.c, ninja_square_sum_number_system_20261031.py (+ .out), ninja_numeration_replay_20261031.py (+ .out)
audit: 04-computation/experiments/ninja_square_sum_20261031_audit_S69/AUDIT.md
related: [THM-4851, HYP-9447, THM-4505, THM-4510]
external:
  - "R. Gerbicz, mersenneforum 'The Square-Sum problem' post #22 (2018-01-21): 'nice pairs', the recursion 49n + res, res = 24..72, base 41..2032"
  - "session124_square_prefix_threshold42_20261031 (weak B bases 42..52, operators valid at m = 42, base block [42, 2081], S_l = (85*49^l - 1)/2, terminal cycles for N >= 2082); session124_recursive_square_prefix_20261031 (98 paired operators); square_sum_B_base_interval_20261019 (nice-pair base 53..2620); session124_gerbicz_reconstruction_20261031"
---

# THM-4850 - weak B pairs start sharply at 42

**Definition** (session124's weak B pair). At n it consists of:
- a Hamiltonian path S of the square-sum graph Q_n (Shorter) and L of Q_(n+1) (Longer);
- both start 1, 3;
- equal positional parity on every common vertex;
- the even-order member ends at 8 and the odd-order member at 10.

Dropping the prefix and endpoint conditions gives Gerbicz's **nice pair**. The owner's "ninja pair" most plausibly
echoes that term (INTERPRETATION).

**(1) Threshold, sharp.** A weak B pair exists iff n >= 42.
- For n <= 53, by exhaustive enumeration:
  - none for n <= 41;
  - for 42 <= n <= 53, the compatible Shorter paths number 10, 4, 45, 3, 5, 11, 859, 108, 153, 230, 10710, 32242,
    and the path pairs 20, 13, 221, 3, 11, 17, 3384, 276, 1090, 1351, 106235, 470586.
- Larger n: session124's certified bases (53..2620) and recursion; equivalently, the threshold42 base block [42, 2081]
  with recursion for n >= 2082.

**(2) Gate.** Among odd orders <= 41, odd-order Hamiltonian paths 1,3 -> 10 exist only at order 41, where there is
exactly one:
1,3,6,30,19,17,32,4,21,15,34,2,7,18,31,33,16,9,40,41,8,28,36,13,23,26,38,11,5,20,29,35,14,22,27,37,12,24,25,39,10.
- Odd-member scarcity forces n >= 40.
- At n = 40 and n = 41 the gate is parity-incompatible with all 11 and all 139 partners respectively.
- So the threshold is 42.

**(3) Normalisations** (audit, exhaustive). The pairs exist exactly for:
- free nice pairs: n in {15, 16, 25, 26, 27} u [29, oo);
- Gerbicz's normalisation (start 1; odd order ends 3, even order ends 8): n >= 41, Gerbicz's base 41..2032;
- weak B pairs: n >= 42.

**(4) Numeration.** Every n >= 42 has exactly one expansion n = 49^k b + sum_{i<k} d_i 49^i with b in [42, 2081] and
d_i in [24, 72]. The digit set is a complete residue system mod 49, which gives uniqueness.
- Session124's 98 paired operators (49 per input parity) make it constructive. Applied to new pairs at m = 42..52, all
  539 outputs pass session124's `check_pair`. This re-confirms threshold42.
- Threshold42's S_l = (85 * 49^l - 1)/2 are the numerals [42; 24, ..., 24].
- An output of order 49m + e takes exactly e - 24 of its 49 offset blocks from the Longer (those c with
  49 + (-1)^(m+1) c <= e) and the rest from the Shorter.
- The odometer n -> n+1 is index arithmetic. B(n) and B(n+1) are not related by inserting a vertex.

**Novelty.**
- Weak B pairs for every n >= 42, the bases 42..52, the base block [42, 2081], the validity of the operators at m = 42
  and S_l are prior work (session124_square_prefix_threshold42_20261031, independently audited there). The
  construction itself is Gerbicz's.
- New here:
  - the exhaustive non-existence for n <= 41, so 42 is sharp;
  - the uniqueness of the gate path;
  - the exact counts for n <= 53;
  - the normalisation table.
