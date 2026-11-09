# Pre-registration: ternary coverage with the complete order-10 strong census

opus-2026-11-09-S87. Written and pushed before `ternary_order10.c` was run.

## Objects

- **Ternary count.** H_3(T) = I(Omega(T), 3), the number of 3-coloured odd-cycle packings (THM-4940).
- **Spec_3.** The set of values H_3(T). It is multiplicative over order joins.
- **Known so far.** THM-4940 (9) and HYP-9530 (A): every value congruent to 1 (mod 3) up to 5398 except 10 and 31 is
  attained, already at order <= 9.
- **The threshold m*.** Let m* be the least alpha-sum m >= 11 (m != 3, 10) such that 1 + 3m is NOT in the
  multiplicative closure of the pieces below. The previous pieces gave m* = 1800.
- **Pieces of the new run:**
  - (i) every strong class of order <= 10, with its exact H_3 (order 10 is new);
  - (ii) every strong class with alpha_1 <= 80 (`few_cycle_strong_80.txt`, S86);
  - (iii) the strong one-apex tournaments and Session 16's two intercyclic rays.

## Calibration (in-sample, from S86's runs)

| strong pieces up to order | m* |
|---|---|
| 8 | 370 |
| 9 | 1800 |

The ratio between the two rows is about 4.9. The largest H_3 of an order-10 tournament is at most 3.375 * 15745 < 53,200,
since alpha_k = 0 for k >= 4 at order 10 and max H(10) = 15745. Hence single order-10 pieces reach alpha-sums of at
most about 17,700.

## Predictions

- **Q1.** m* >= 5000. Every value congruent to 1 (mod 3) up to 15,001, except 10 and 31, is attained.
- **Q2.** m* <= 20000. The construction universe of order <= 10 pieces runs out before alpha-sum 20,000.
- **Q3.** Every irreducible element of the monoid M_3 = {v = 1 (mod 3)} minus {10, 31} below 1 + 3m* is the H_3 of a
  STRONG tournament of order <= 10, or of a strong piece in (ii) or (iii). This holds by the definition of m* if the
  closure is formed correctly; it is listed as a consistency check of the lane reduction.

A failure of Q1 or Q2 says nothing about the truth of ternary Chan-Pak; it only measures the construction universe.

Outcome: to be appended below without editing the predictions.

## Outcome (appended after the run; the predictions above are unedited)

`ternary_order10.c` gave 10,070 distinct strong order-10 H_3 values (max 31,549), and `ternary_coverage10.py` formed
the closure (coverage10.out). The independent audit reproduced both, the order-10 set by a different method
(audit/ext10.out, audit/closure_check.out).

- **Q1 CONFIRMED.** m* = 9679 >= 5000. Every value = 1 (mod 3) up to 29,035, except 10 and 31, is attained.
- **Q2 CONFIRMED.** m* = 9679 <= 20000. The first uncovered value is 29,038 = 2 * 14,519, an irreducible of lane Q.
- **Q3 CONFIRMED.** All 4,031 irreducibles of M3 below 29,038 are piece values. As stated in the prediction, this
  holds automatically once the closure is formed correctly: an irreducible covered by the closure must be a single
  piece value.

Afterwards the audit observed that the strong classes of order <= 10 alone already give the same m*. Pieces (ii) and
(iii) contribute nothing below 29,038.
