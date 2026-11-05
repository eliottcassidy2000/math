---
id: HYP-9466
title: "Pre-registered order-7 checks of digit heredity under C2. For every class B of order 7 and every mark a: the child C2(B, a) has sigma_x = 0 and sigma_y = h(B) - e_a(B); exactly 6 children carry a digit-1 mark and exactly 4 a digit-2 mark; no grandchild C2(C2(B, a), u) carries a digit-1 mark and exactly 3 carry a digit-2 mark"
status: CONFIRMED (P1-P3 all hold at |B| = 7; run after the pre-registration commit d96e81539f); recorded as THM-4873 (iii)
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

Output: `c2_digit_heredity_20261104.out`.
- **(P1) holds.** sigma_x = 0 and sigma_y = h - e_a, with 0 failures over 3,192 pairs (B, a).
- **(P2) holds.** 6 children carry a digit-1 mark and 4 a digit-2 mark.
- **(P3) holds.** None of the 28,728 grandchildren carries a digit-1 mark; 3 carry a digit-2 mark. Each of the three
  has two digit-2 marks, the native y2 and the old mark u (S77 audit M2). The smallest positive digit exceeds 20 in
  28,463 grandchildren.
- **The whole statement is PROVED at every order** (THM-4873 (ii), the S77 audit's proof).

The native-mark part is PROVED by the end-count lemma (THM-4873 (i), (ii)).
