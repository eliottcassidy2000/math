---
id: HYP-9447
title: "Terminal-form square-sum cycles (through 8-1-3 for even N, 15-1-3 for odd N) exist for every N >= 40; N = 39 is the last hole after the first terminal forms at 34/35; pre-registered search for 63 <= N <= 1000"
status: OPEN (PRE-REGISTERED; the confirmatory run follows this commit)
session: opus-2026-10-31-S69
date: 2026-10-31
related: [THM-4850, THM-4851]
note: 05-knowledge/results/ninja_pairs_square_sum_number_system_20261031.md
script: 04-computation/experiments/terminal_cycles_large_20261031.c
---

# HYP-9447 - terminal-form cycles from 40 on

## Statement

Call a Hamiltonian cycle of the square-sum graph Q_N **terminal** if it contains 8-1-3 when N is even and 15-1-3 when
N is odd. These are session124's closing forms: 1+3 = 4, 1+8 = 9, 1+15 = 16.

**Conjecture:** a terminal cycle exists for every N >= 40.

## Known before this pre-registration

- **Exhaustive (this session, `ninja_pairs_square_sum_20261031.c` and `terminal_cycles_large_20261031.c`):** terminal
  cycles exist for N in {34, 35, 36, 37, 38} and for every N in [40, 62]. They do not exist for N <= 33, nor for
  **N = 39**: Q_39 has 50 Hamiltonian cycles, none through 8-1-3 or 15-1-3; the complete search terminates.
- **Even N >= 42 is PROVED.** Every ninja (weak B) pair at n >= 42 (THM-4850) has an even-order member that starts 1,3
  and ends at 8, and closing it by 8+1 = 9 gives a terminal cycle.
- **Odd N >= 2621** is session124's terminal-cycle theorem.
- **Open window:** odd N in [63, 2619].

## Pre-registered predictions (confirmatory run, budget 10^9 search nodes per N)

Command: `terminal_cycles_large 63 1000 1 1000000000`.

- **(P1)** Every odd N in [63, 1000] is FOUND, with a verified terminal cycle through 15-1-3.
- **(P2)** Every even N in [64, 1000] is FOUND. This is a consistency check: it is implied by THM-4850 together with the
  certified base.
- **(P3)** No N in [63, 1000] returns EXHAUSTED (a complete search with no terminal cycle). So 39 remains the last hole.

## Reporting rules

- An UNDECIDED (budget exhausted) N is reported as such. P1 or P2 then counts as NOT CONFIRMED for that N, and is not
  refuted.
- An EXHAUSTED N refutes P3 and the conjecture.
- Every FOUND cycle is re-verified by the program's independent checker: vertex set, square sums and the closing edge.
