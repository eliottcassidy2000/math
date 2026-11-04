---
id: HYP-9464
title: "Coincidences with a mechanism. (A) The {-3}-self-dual tournaments of order 10 are exactly the transitive TT_10 and the almost transitive AT_10 = NT_10, whose path count 257 = 2^8 + 1 is the Fermat prime F_3; in general the non-transitive decomposable {-3}-self-dual tournament of order 2^i + 2 >= 10 is the one with H = F_i (Achour-Boudabbous-Boussairi, with THM-4749 (4)). (B) Every sole-exit-free tournament with exactly 3 Hamiltonian cycles (n >= 6) has exactly 3 vertices whose deletion is uniquely Hamiltonian; with exactly 4 cycles (n >= 7), exactly 2. (C) The few-cycle Fibonacci table is the rooted uniquely Hamiltonian count two orders down: u_3(n) = R(n - 2), R(m) = F(2m - 5) (THM-4868 (F)). Pre-registered at order 10"
status: CONFIRMED at order 10 (pre-registered in 3fc61b533c, all three predictions held, no failing class); audited PASS WITH CORRECTIONS (MISTAKE-712). (A) FINITE-EXACT at orders 9-11 (11 by the audit) and CITED for decomposable n >= 9, indecomposable n >= 12 OPEN; (B) FINITE-EXACT at orders 6-14 (11-14 by the audit), general n OPEN; (C) PROVED
script_confirm: 04-computation/experiments/selfdual3_fewcycle_20261104_confirm.out, 04-computation/experiments/selfdual3_networkx_check_20261104.py, 04-computation/experiments/at_mirror_fixed_paths_20261104.py, 04-computation/experiments/selfdual3_fewcycle_20261104_records/
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
- At order 8 there are 16 decomposable classes, so 14 exceptions besides TT_8 and AT_8, together with one
  non-self-converse indecomposable pair (16 + 2 = 18). *Corrected after the order-10 runs:* the pre-registered text
  said "16 decomposable exceptions besides TT_8 and AT_8", but the record shows 16 decomposable classes in all. No
  prediction depends on this count.

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

The runs were made after the pre-registration commit 3fc61b533c. The record is
`selfdual3_fewcycle_20261104_confirm.out`. All three predictions are **CONFIRMED**, with no failing class.

- **(P1) CONFIRMED.** `selfdual3 10` checked all 9,733,056 classes of order 10 (38 s). Exactly two are
  {-3}-self-dual:
  - TT_10 (H = 1, hc = 0);
  - one class with scores 1 1 2 3 4 5 6 7 8 8, H = 257, hc = 1.

  Both are decomposable and self-converse; none is indecomposable. Independent networkx code
  (`selfdual3_networkx_check_20261104.py`) identifies the two classes as TT_10 and AT_10. It also confirms that AT_n
  and TT_n are {-3}-self-dual for n = 9, 10, 11 and that three random order-10 tournaments are not.
- **(P2) CONFIRMED.** `hcdump 10 3 2` lists 89 classes (= u_3(10) = F(11)). Each has exactly 3 vertices v with
  T - v uniquely Hamiltonian.
- **(P3) CONFIRMED.** `hcdump 10 4 2` lists 68 classes (= u_4(10) = 2F(9)). Each has exactly 2.

**Typing after the runs.**

| part | status |
|---|---|
| (A) at orders 9 and 10 | FINITE-EXACT (exhaustive): the {-3}-self-dual classes are exactly TT_n and AT_n |
| (A) at order 11 | FINITE-EXACT by the audit's complete pruned search (validated at orders 7, 8): only TT_11 and AT_11 |
| (A) the decomposable case, every n >= 9 | CITED (Achour-Boudabbous-Boussairi) |
| (A) no indecomposable {-3}-self-dual tournament for n >= 12 | OPEN |
| (A) H(AT_n) = 2^(n-2) + 1, hence F_i at order 2^i + 2 | PROVED (THM-4749 (4); re-derived in the note) |
| (B) hc = 3, exactly 3 uniquely Hamiltonian deletions | FINITE-EXACT, orders 6-14 (11-14 by the audit's `aud_fewcycle.c`); OPEN in general |
| (B) hc = 4, exactly 2 | FINITE-EXACT, orders 7-14 (11-14 by the audit); OPEN in general |
| (C) u_3(n) = R(n - 2) | PROVED (algebra from THM-4869); no bijection is known |

**The mirror law.** The note proves it for AT_n (`at_mirror_fixed_paths_20261104.py`). The converse involution of
AT_n fixes exactly H(AT_(k+2)) Hamiltonian paths when n = 2k + 2, and only one path when n is odd. So the Fermat
recursion F_i = (F_(i-1) - 1)^2 + 1 is the count of a fixed-point set: at order 2^i + 2 the mirror-fixed paths number
F_(i-1). The involution and the parity congruence are THM-582 and THM-647 (with THM-024 and THM-583). The AT_n
fixed-point count is the new part.

**Audit observation (not pre-registered).** For hc = 5 at orders 7-14 the deletion counts split term by term along
S75's u_5 = L(2n - 12) + F(2n - 13):
- L(2n - 12) classes have exactly 2 uniquely Hamiltonian deletions;
- F(2n - 13) classes have exactly 3.

This is a candidate extension of (B), to be pre-registered before any confirmatory run.

**Audit.** The independent adversarial audit gave PASS WITH CORRECTIONS (1 major, 11 minor; MISTAKE-712). Its report
is `04-computation/experiments/fermat_self_duality_20261104_audit_S10/AUDIT.md`. Two corrections concern this file:
- the order-8 count in the exploratory section, corrected above;
- the networkx check, now committed.
