---
id: HYP-9465
title: "Pre-registered checks for THM-4872: (P1) the sole-exit-free three-cycle tournaments of orders 13 and 14 split by automorphism group as F(2n - 10) with Aut = Z_3 (987, 2,584) and F(2n - 11) asymmetric (610, 1,597), and the sole-exit-free four- and five-cycle tournaments are all asymmetric; (P2) the unique-start rooted strong tournaments of order 14 number F(24) = 46,368; (P3) at order 14 at least one unit-start +-2 port pair exists (F = F', h' = h + 2)"
status: PRE-REGISTERED (opus-2026-11-04-S76); the order 13-14 runs are made after this commit
session: opus-2026-11-04-S76
date: 2026-11-04
related: [THM-4872, THM-4869, THM-4868, HYP-9461, HYP-9464]
script: 04-computation/experiments/fewcycle_symmetry_and_unit_pairs_20261104.c, 04-computation/experiments/unit_start_pm2_pairs_grammar_20261104.c
---

# HYP-9465 - unit-start +-2 port pairs and the symmetry split

## Evidence before this commit (FINITE-EXACT)

**The symmetry split** (`fewcycle_symmetry_and_unit_pairs_20261104.out`). At orders 6-12, the sole-exit-free classes
with exactly 3 Hamiltonian cycles have |Aut| = 1 for 1, 2, 5, 13, 34, 89, 233 classes and |Aut| = 3 for 1, 3, 8, 21,
55, 144, 377 classes. Every sole-exit-free class with hc = 4 or 5 at orders 7-12 is asymmetric.

**Why (proof sketch, THM-4869).** Aut(T) fixes the unique terminal core D, its tail u and head l. It also fixes the
uniquely Hamiltonian rooted prefix pointwise (rigidity, THM-4868 (E)) and the transitive chain. So Aut(T) = Aut(W, A)
of the block:
- Z_3 for (C3, all three vertices);
- trivial for the W4 blocks and the uniquely Hamiltonian 5-blocks.

The C3-block cores have generating function x^5/(1 - 2x). Multiplied by the prefix factor x(1 - 2x)/(1 - 3x + x^2),
this gives x^6/(1 - 3x + x^2), i.e. F(2n - 10).

**The unit-start +-2 pairs** (`plus_minus_two_port_routes_20261104.out`, `unit_start_pm2_pairs_20261104.out`,
`unit_start_pm2_pairs_grammar_20261104` run to order 13):
- At order 7, exactly one pair: (h, s, q, F) = (29, 1, 15, 157) and (31, 1, 9, 157).
- Pairs at orders 7-13: 1, 2, 5, 1, 6, 10, 4.

## Predictions

- **(P1)** Orders 13, 14: hc = 3 sole-exit-free classes with |Aut| = 3: 987, 2,584; asymmetric: 610, 1,597. All 1,220
  and 1,076 (order 13) and 3,194 and 2,817 (order 14) classes with hc = 4 and 5 are asymmetric.
- **(P2)** Unique-start rooted strong classes of order 14: 46,368 (FUNNELS F2), every one with s = 1.
- **(P3)** At least one unit-start +-2 port pair at order 14.

## Results

(Filled after the runs.)
