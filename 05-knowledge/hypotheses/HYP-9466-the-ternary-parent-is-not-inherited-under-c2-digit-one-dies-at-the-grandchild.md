---
id: HYP-9466
title: "Pre-registered order-7 checks of digit heredity under C2. For every class B of order 7 and every mark a: the child C2(B, a) has sigma_x = 0 and sigma_y = h(B) - e_a(B); exactly 6 children carry a digit-1 mark and exactly 4 a digit-2 mark; no grandchild C2(C2(B, a), u) carries a digit-1 mark and exactly 3 carry a digit-2 mark"
status: PRE-REGISTERED (opus-2026-11-04-S77); the order-7 run is made after this commit
session: opus-2026-11-04-S77
date: 2026-11-04
related: [THM-4873, THM-4862, THM-4865, HYP-9454]
script: 04-computation/experiments/c2_digit_heredity_20261104.c
---

# HYP-9466 - digit 1 is not inherited under C2

**Setting.** For a tournament B and a mark a, the C2 child adds x with x -> a only, and y beaten by all of B with
y -> x. The law is H = h(B) + 2 B_a(B) (multiport_defect_grammar law (4); klein THM-4865 at any vertex). The digit of a
mark u is sigma_u = B_u - H.

**Evidence before this commit** (`c2_digit_heredity_20261104`, all classes of order 1-6):
- sigma_x(child) = 0 and sigma_y(child) = h - e_a hold in every case.
- Children carrying a digit-1 mark: 1, 2, 3, 4, 5 at |B| = 2..6.
- Children carrying a digit-2 mark: 1, 4, 4, 4, 4.
- Grandchildren (all marks at both levels): 3, 8, 30, 96, 420, 2,688. None carries a digit-1 mark; exactly 3 carry a
  digit-2 mark at every order.

## Predictions for |B| = 7 (456 classes)

- **(P1)** sigma_x = 0 and sigma_y = h - e_a for all 3,192 pairs (B, a).
- **(P2)** Exactly 6 children carry a digit-1 mark, and exactly 4 a digit-2 mark.
- **(P3)** Of the 28,728 grandchildren, none carries a digit-1 mark, and exactly 3 carry a digit-2 mark.

## Results

(Filled after the run.)
