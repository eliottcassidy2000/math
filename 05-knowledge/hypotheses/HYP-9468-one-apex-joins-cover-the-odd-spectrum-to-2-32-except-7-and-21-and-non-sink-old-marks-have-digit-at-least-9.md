---
id: HYP-9468
title: "Pre-registered out-of-sample checks. (P1) Every odd N <= 2^32 other than 7 and 21 is the Hamiltonian path count of a one-arc join of a one-apex tournament with a one-apex tournament or with a tournament of order <= 7. (P2) The one-apex x one-apex joins alone miss a vanishing proportion of odd numbers (growth ratio and density windows below). (P3) In C2 children of all order-8 bases, every old vertex other than the mark and the sink of B has digit >= 9, and the old-mark checks of THM-4875 (A) all pass"
status: CONFIRMED (P1, P2 and P3 all hold out of sample, run after the pre-registration commit 19562af223; the independent audit reproduced P1 and P2 to 2^32 and P3 at order 8); P3's all-order conjecture PROVED by the audit (THM-4875 (A3)); recorded as THM-4875 (D) and (A3)
session: opus-2026-11-04-S78
date: 2026-11-04
related: [THM-4875, THM-4873, THM-4866, THM-4865, HYP-9454]
script: 04-computation/experiments/chanpak_one_apex_join_cover_20261104.c; 04-computation/experiments/old_mark_digit_bound_20261104.c
---

# HYP-9468 - one-apex joins and the odd spectrum; the non-sink old digit floor

## Setting

**One-arc joins** (session29 ADDENDUM sections 2 and 5; THM-4866 (B)). For rooted tournaments (A, a) and (B, b), let
J = A => B with b -> a reversed. Write x = (h, e, R) = (h(A), e_a(A), R_a(A)), where R_a(A) is the root response
sum_{X subset V(A) - a} h(A[X]) s_a(A - X), and y = the same triple for (B^op, b). Then

    H(J) = x_h y_h - x_e y_e + x_R y_R.

**One-apex tournaments** (session25; Session191 one_apex/PROOF.md). For a binary word w with value v, let T(w) be the
transitive chain 1 -> 2 -> ... -> n together with x, where x -> i iff the i-th bit is 1. The rooted triple at x is

    (g(v), e(v), 1 + v),   g(0) = 1,  g(2v + 1) = g(v),  g(2v) = g(v) + v + (v mod 2),   e(v) = [v = 0 or v even].

This was checked by brute force for every word of length <= 7. The value R_x = 1 + v is the inherited binary state A + 1.
The join of two one-apex tournaments therefore has

    H = (1 + u)(1 + v) + g(u) g(v) - e(u) e(v),

checked by brute force on all 961 pairs of words of length <= 4.

**Families.**
- F0: the plain counts, namely every h of order <= 8 and every g(v).
- F1: one-apex x one-apex joins.
- F2: one-apex x bank joins, g(v) y_h - e(v) y_e + (1 + v) y_R, where y runs over the 2,732 distinct rooted triples of
  order <= 7. The bank is closed under the converse; its count matches THM-4866 (D).

## Exploratory evidence before this commit (`join_cover.c`, `bilinear_cover.c`)

- On [1, 2^22]: F0 leaves 1,315,007 odd numbers, F0 + F1 leaves 27,142, and F0 + F1 + F2 leaves exactly {7, 21}.
- F1 alone misses, in [2^k, 2^(k+1)) for k = 17..21: 1,634, 2,553, 3,906, 6,209 and 10,051.
  - Ratios: 1.56, 1.53, 1.59, 1.62.
  - Density in [2^21, 2^22): 0.0096.
- A greedy sub-bank needs 111 triples on top of F1 to close [23, 2^20]. The coverage is redundancy across many
  triples, not a small explicit family.
- Old marks, orders |B| <= 7:
  - every check of THM-4875 (A) passes;
  - the minimum old digit is 4, attained once per order;
  - the minimum digit at a non-sink old vertex is 9 at every order 2-7;
  - the minimum mark digit at out-degree 1 is 2.

## Predictions (out of sample)

- **(P1)** F0 + F1 + F2 attains every odd N in [1, 2^32] except 7 and 21. Any other miss refutes P1; the misses are to
  be recorded.
- **(P2)** F1 alone has miss counts m_k in [2^k, 2^(k+1)) for k = 22..31 with every ratio m_(k+1) / m_k in
  [1.45, 1.75]. Its miss density in [2^31, 2^32) lies in [0.0005, 0.003].
- **(P3)** For the 6,880 classes of order 8 and all 55,040 marks:
  - every check of THM-4875 (A) passes;
  - the minimum old digit is 4, attained by exactly one (class, mark, vertex) triple;
  - the minimum digit at a non-sink old vertex is 9;
  - the minimum mark digit at out-degree 1 is 2.
  - Conjecture for all orders: every old vertex w != a of C2(B, a) that is not the sink of B has digit >= 9.

## Results

**(P1) CONFIRMED.** `chanpak_one_apex_join_cover_20261104.out` (LOGX = 32): the odd misses of F0 + F1 + F2 in
[1, 2^32] are exactly 7 and 21.
- The suffix law for g was checked on 2,000,000 random v (0 failures).
- The run used the order-<=7 bank for F0 as well, so the confirmed statement is slightly stronger than the "order
  <= 8" written above.

**(P2) CONFIRMED.**
- F1-only misses m_k for k = 22..31: 15,827, 25,815, 40,446, 63,546, 99,877, 154,867, 241,881, 378,695, 593,645 and
  934,107.
- Ratios: 1.575, 1.631, 1.567, 1.571, 1.572, 1.551, 1.562, 1.566, 1.568, 1.574 (all in [1.45, 1.75]).
- Density in [2^31, 2^32): 0.000870 (in [0.0005, 0.003]).

**(P3) CONFIRMED.** `old_mark_digit_bound_20261104.out`, |B| = 8: 6,880 classes and 55,040 marks.
- 0 failures in every check.
- Exactly one digit-4 old vertex.
- The minimum non-sink old digit is 9, and the minimum mark digit at out-degree 1 is 2.

## After the independent audit (PASS WITH CORRECTIONS; MISTAKE-716)

Audit report: `04-computation/experiments/old_mark_one_apex_join_20261104_audit_S78/AUDIT.md`.

- **Reproduction.** With its own g, bank and sieve, the audit's 2^32 run misses exactly {7, 21}, and its 22 F1-only
  block counts equal the ones above. Its own class generation reproduces P3 at order 8. File times show both
  pre-registered outputs postdate 19562af223.
- **P2 indexing.** Read literally, P2 asked for m_(k+1)/m_k with k = 22..31, which needs m_32 (outside the run). The
  ratios listed above are m_k/m_(k-1) for k = 22..31. All nine ratios inside the literal window (k = 22..30) lie in
  [1.45, 1.75], as does m_22/m_21.
- **P3's all-order conjecture is PROVED** (audit U1; THM-4875 (A3)). Every non-sink old vertex w != a of C2(B, a) has
  digit >= 9, with equality iff B is transitive, a is its sink and w its penultimate vertex. The author-side labelled
  check `nonsink_digit_nine_check_20261104.py` (orders 2-6, 972,230 cases) finds exactly n! digit-9 cases per order.
- **F0 is redundant** (audit U3). g(2t) = F1(0, t) and g(2t + 1) = g(t), as also in Session 193 one_arc_joins/PROOF.md
  section 3. F1 + F2 alone miss only 7 and 21 up to 2^28 (and so up to 2^32).
- **Credits.**
  - The one-apex triple is proved in session179 spectrum/REPORT.md (1).
  - g is session124's f.
  - THM-4876 (6) is the same family.
  - F1 lies inside the two-feedback-vertex family, whose completeness is the OPEN conjecture of
    square_sum_spectrum_overnight_20261030/README.md. Session191 two_apex covers that family to 230,519.
  - P2's density statement is about a subfamily of it.
