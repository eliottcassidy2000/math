---
id: THM-4851
title: "The prehistory of the square-sum numbers: Q_32 and Q_33 each have a unique Hamiltonian cycle, the 33-cycle being the 32-cycle with one segment reversed and 33 inserted (both close 8-1-15); terminal-form cycles (8-1-3 for even N, 15-1-3 for odd N) exist exactly for N in {34,...,38} u [40, oo) (the owner's 35 is the first odd one; 39 is the only later hole); from 32 to 42 adding a node always deletes at least two old edges"
status: PROVED (exhaustive to N = 62; verified search 63..2620, pre-registered to 1000 as HYP-9447; cited session124 terminal theorem for N >= 2621; THM-4850 for even N >= 42) + CITED (OEIS A071983, A071984; Gerbicz) + FINITE-EXACT
session: opus-2026-10-31-S69
date: 2026-10-31
note: 05-knowledge/results/ninja_pairs_square_sum_number_system_20261031.md
script: 04-computation/experiments/ninja_pairs_square_sum_20261031.c, terminal_cycles_large_20261031.c (+ .out, _window.out), ninja_square_sum_number_system_20261031.py (+ .out)
related: [THM-4850, HYP-9447, THM-4505, THM-4510]
---

# THM-4851 - the prehistory of the square-sum numbers

**(a) Paths and cycles** (CITED, reproduced).
- Q_N has a Hamiltonian path iff N in {1, 15, 16, 17, 23} u [25, oo).
- Q_N has a Hamiltonian cycle iff N >= 32 (Gerbicz).
- The counts for N <= 42 equal OEIS A071983 and A071984.

**(b) The two unique cycles.** Q_32 has exactly one Hamiltonian cycle:
1,8,28,21,4,32,17,19,30,6,3,13,12,24,25,11,5,31,18,7,29,20,16,9,27,22,14,2,23,26,10,15.
- Q_33 has exactly one: the same with 31,18,7,29,20 reversed and 33 inserted. The edges 5-31 and 20-16 are deleted
  and 5-20, 31-33, 33-16 are added.
- Both contain 8-1-15.

**(c) Terminal spectrum.** Q_N has a Hamiltonian cycle through 8-1-3 (N even) or 15-1-3 (N odd) iff
N in {34, 35, 36, 37, 38} u [40, oo).
- N <= 62: exhaustive; the absences at 32, 33 and 39 come from complete searches. Q_39 has 50 cycles, none terminal.
- 63 <= N <= 2620: every N has a cycle found by search and verified independently.
- N >= 2621: session124's terminal-cycle theorem.
- Every even N >= 42 also follows from THM-4850, since the even ninja member closes by 8 + 1 = 9.

**(d) Insertion costs.** For 32 <= N <= 41, every Hamiltonian cycle of Q_(N+1), minus the vertex N+1, misses at least
two edges of every Hamiltonian cycle of Q_N. The minima over cycle pairs are 2, 3, 2, 3, 3, 4, 4, 3, 2, 2.

**Novelty.**
- (a) is classical.
- New to the repository:
  - the reversal-insertion description of the 32 -> 33 step;
  - the exact terminal spectrum and the hole at 39;
  - the insertion-cost table.
