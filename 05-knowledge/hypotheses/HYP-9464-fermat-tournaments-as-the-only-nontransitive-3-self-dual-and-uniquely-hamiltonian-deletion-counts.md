---
id: HYP-9464
title: "Coincidences with a mechanism. (A) The {-3}-self-dual tournaments of order 10 are exactly the transitive TT_10 and the almost transitive AT_10 = NT_10, whose path count 257 = 2^8 + 1 is the Fermat prime F_3; in general the non-transitive decomposable {-3}-self-dual tournament of order 2^i + 2 >= 10 is the one with H = F_i (Achour-Boudabbous-Boussairi, with THM-4749 (4)). (B) Every sole-exit-free tournament with exactly 3 Hamiltonian cycles (n >= 6) has exactly 3 vertices whose deletion is uniquely Hamiltonian; with exactly 4 cycles (n >= 7), exactly 2. (C) The few-cycle Fibonacci table is the rooted uniquely Hamiltonian count two orders down: u_3(n) = R(n - 2), R(m) = F(2m - 5) (THM-4868 (F)). Pre-registered at order 10"
status: PRE-REGISTERED (kind-pasteur-2026-11-04-S10); the order-10 runs are made after this commit
session: kind-pasteur-2026-11-04-S10
date: 2026-11-04
related: [THM-4868, THM-4869, THM-4870, THM-4749, THM-4866, HYP-9459, HYP-9461, HYP-9462]
note: 05-knowledge/results/fermat_self_duality_fewcycle_coincidences_20261104.md
script: 04-computation/experiments/selfdual3_fewcycle_20261104.c, 04-computation/experiments/uh_deletion_count_20261104.py
---

# HYP-9464 - Fermat tournaments, {-3}-self-duality and uniquely Hamiltonian deletions

## Inheritance

- **Almost transitive tournaments.** AT_n is the transitive tournament with the arc between its two extreme vertices
  reversed (Busch's "nearly transitive", Moon's NT_n). It is uniquely Hamiltonian, with a non-critical run of length
  n - 2 and 3 cut arcs (THM-4868 (C), (D); opus-S75, MISTAKE-709). Its path count is H(AT_n) = 2^(n-2) + 1
  (THM-1370, THM-4749 (4)), so H(AT_(2^i + 2)) = F_i.
- **The self-duality theorem** (Achour-Boudabbous-Boussairi, arXiv:1204.2513; CITED, abstract and statement read by
  S75). A decomposable tournament with at least 9 vertices is {-3}-self dual iff it is transitive or almost
  transitive. Here {-3}-self dual means that T - X is self-converse for every 3-set X.
- **Few-cycle counts.** THM-4869 (opus-S75, audited) gives
  u_k(n) = sum_d c_k(d) R(n - d + 1), with R(m) = F(2m - 5) the rooted (uniquely Hamiltonian, sole exit) count of
  THM-4868 (F). Its corollaries are u_3 = F(2n - 9) and u_4 = 2F(2n - 11).

## Exploratory evidence (FINITE-EXACT, before this commit; `selfdual3_fewcycle_20261104_explore.out`)

**{-3}-self-dual classes** at orders 6, 7, 8, 9: 56 (all), 10, 18, 2.
- At order 9 the two classes are TT_9 and AT_9 (H = 129, hc = 1). This checks the cited theorem at its threshold,
  with no indecomposable exception.
- Every module of every decomposable {-3}-self-dual class at orders 7-9 is transitive.
- At order 8 there are 16 decomposable exceptions besides TT_8 and AT_8, together with one non-self-converse
  indecomposable pair.

**Uniquely Hamiltonian deletions** in sole-exit-free tournaments:

| hc | orders | number of v with T - v uniquely Hamiltonian |
|---|---|---|
| 3 | 6-9 (2 + 5 + 13 + 34 = 54 classes) | always exactly 3 |
| 4 | 7-9 (4 + 10 + 26 classes) | always exactly 2 (at n = 6 the values are 2, 3, 4) |
| 5 | 7-9 | 2 or 3 |
| 2 | 5 (RT_5) | all 5 |

**(C) is algebra.** It follows from THM-4869's generating function: U_3(x) = x^2 (R(x) - x^2 - x^3). Each term
checks: R(4..8) = 2, 5, 13, 34, 89 = u_3(6..10).

## Pre-registered predictions (order 10)

- **(P1)** `selfdual3_fewcycle_20261104 selfdual3 10` finds exactly two {-3}-self-dual classes of order 10:
  - TT_10 (H = 1);
  - AT_10 (H = 257 = F_3, hc = 1).

  Both are decomposable, and there is no indecomposable one. The "no indecomposable" part is a genuine extrapolation;
  the decomposable part is the cited theorem.
- **(P2)** `hcdump 10 3 2` lists the 89 sole-exit-free classes with exactly 3 Hamiltonian cycles. For every one of them,
  exactly 3 vertices have T - v uniquely Hamiltonian.
- **(P3)** `hcdump 10 4 2` lists the 68 classes with exactly 4 cycles. For every one of them, exactly 2 vertices do.

Every outcome will be reported, with the class if a prediction fails.

## Results

(Filled after the order-10 runs.)
