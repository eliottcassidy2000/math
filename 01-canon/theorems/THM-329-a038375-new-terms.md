---
id: THM-329
title: New Terms of A038375 (Max Hamiltonian Paths in Tournament)
status: LOWER BOUNDS (a(12) strongly believed exact; a(13) not confirmed)
session: opus-2026-05-27-S6
tags: [A038375, hamiltonian-paths, tournaments, computational]
related: [THM-326, THM-327, THM-328]
---

## Statement

A038375 = maximum number of Hamiltonian paths in a tournament on $n$ nodes.

**New terms (beyond OEIS):**

$$a(12) \geq 531205$$

$$a(13) \geq 3719831$$

## Evidence for a(12) = 531205

- 8 independent trials (90s each, distinct seeds): 4/8 achieve exactly 531205.
- All 8 trials converge to either 531175 or 531205 — two local optima.
- Multiple DISTINCT tournaments achieve 531205 (different adjacency lists), confirming it's not a single artifact.
- The search that finds 531175 appears to be a sub-optimal basin; 531205 is the larger basin.
- No trial has found anything > 531205 despite hundreds of restarts.

## Warm starts for n=12

The Paley tournament is undefined for n=12 (12 ≢ 3 mod 4). Circulant warm starts find ~530000; the best is reached via random restarts typically at restart #400-500.

## Evidence for a(13)

Multiple independent trials (3min to 10min, 8 distinct seeds):
- Values found: 3711175, 3711303, 3711611, 3714875, 3715419, 3716335, 3717455, 3717753, 3719435, **3719831** (best)
- Two independent seeds (0x33333333 and 0x77777777) find 3719831.
- Search has NOT converged as cleanly as n=12; true maximum unclear.
- Best tournament: adj = 3666 5352 6411 3009 686 3213 6708 5701 243 3366 4444 5522 825

Note: $n=13$ is prime but $13 \equiv 1 \pmod{4}$, so the Paley tournament is undefined (would create symmetric arcs, not a tournament).

## Lower bound for a(14)

A 3-minute search (seed 0x11111111) found $a(14) \geq 24786711$.
Best tournament: adj = 254 508 504 2032 3040 8128 16256 16128 15873 15367 14359 12303 8223 63

This is likely improvable with more search time.

## Context

Known values (OEIS A038375):
$a(1..11) = 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095$

For prime $p \equiv 3 \pmod{4}$: the Paley tournament $QR_p$ achieves $a(p)$.
- $p=3$: $a(3)=3$, $QR_3$ achieves it (immediately on warm start)
- $p=7$: $a(7)=189$, $QR_7$ achieves it (immediately)
- $p=11$: $a(11)=95095$, $QR_{11}$ achieves it (immediately)
- Next: $p=19$ ($19 \equiv 3 \pmod{4}$). Prediction: $a(19) = H(QR_{19})$.

## Growth rate

Ratios $a(n+1)/a(n)$:
$3, 5/3, 3, 3, 4.2, 3.5, 5.1, 4.7, 6.0, {\bf 5.6}, {\bf 7.0}$ (bold = our estimates)

The ratio oscillates but trends upward. The pattern for Paley positions (n=3,7,11): the next Paley prime is $p=19$, so $a(19)$ will likely show a large ratio.

## Algorithm

Local search (hill climbing) via bitmask DP:
- HP count via $O(2^n \cdot n^2)$ DP
- Warm starts: Paley (if $n \equiv 3 \pmod{4}$ prime), two circulant families, greedy extension from $(n-1)$-optimal
- Random restarts with best-so-far perturbation
- Validated: all known values $a(1..11)$ confirmed exactly

## Files

- Solver: `04-computation/a038375_solver.c`
- Results: `05-knowledge/results/a038375.out`
