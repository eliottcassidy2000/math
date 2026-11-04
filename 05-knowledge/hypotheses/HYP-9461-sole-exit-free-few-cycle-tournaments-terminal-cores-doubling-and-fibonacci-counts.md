---
id: HYP-9461
title: "Sole-exit-free tournaments with few Hamiltonian cycles factor through a terminal core: for n >= 7, every strong tournament with hc in {3, 4, 5} and no vertex of out-degree 1 is a one-arc join with a unique minimal one-arc-leaving set D (tail u), e_u(D) = hc and a uniquely Hamiltonian rooted prefix; the cores double (c_3(d) = 3 * 2^(d-6), c_4(d) = 2^(d-5), c_5(d) = 5 * 2^(d-7)), which gives klein's F(2n - 9) and 2F(2n - 11) (HYP-9459 (2), (3)) and predicts u_5(n) = 4F(2n - 12) - 3F(2n - 14). Pre-registered at orders 10-14"
status: CONFIRMED (P1-P10 all hold at orders 10-14, blind at 11-14 since the order-10 totals were klein's published values; run after the pre-registration commit c355721ce6) and PROVED for all n as THM-4869 (opus-2026-11-04-S75)
session: opus-2026-11-04-S75
date: 2026-11-04
related: [THM-4869, HYP-9459, THM-4868, THM-4867, HYP-9458, HYP-9460, THM-4862]
script: 04-computation/experiments/fewcycle_sole_exit_free_20261104.c
---

# HYP-9461 - sole-exit-free few-cycle tournaments and their terminal cores

## Notation

- hc(T) is the number of Hamiltonian cycles.
- A **sole exit** is a vertex of out-degree 1.
- t_k(n) counts the classes of order n with hc = k; u_k(n) counts those with no sole exit. This is klein's notation in
  HYP-9459, whose exhaustive table at orders 5-10 the program reproduces at 5-9; order 10 is the first run below.
- A **one-arc-leaving set** D is a non-empty proper vertex set with exactly one arc u -> l leaving it.
  - e_u(D) is the number of Hamiltonian paths of T[D] ending at u.
  - s_l(T - D) is the number of Hamiltonian paths of T - D starting at l.
- **Product law** (elementary): hc(T) = s_l(T - D) * e_u(D). Every Hamiltonian cycle traverses D as one block ending
  at u, followed by T - D as one block starting at l.
- R(m) is the number of rooted pairs (Q, z), with Q uniquely Hamiltonian of order m and z a sole exit of Q:
  - R(m) = F(2m - 5) for m >= 3 (THM-4868 (F));
  - R(2) = 1 (degenerate prefix: T - D is the single vertex l).
- A **core of type k and order d** is a pair (D, u) such that:
  - D is strong of order d, and e_u(D) = k;
  - no vertex of D other than u has out-degree 1 in D;
  - no non-empty Z inside D - u has exactly one arc into D - Z (minimality).

  c_k(d) is the number of cores up to isomorphism.

## Evidence before this commit (orders 3-9, all FINITE-EXACT)

The program generates every strong tournament with hc <= 5. The generation is complete because hc(T - v) <= hc(T) for
a non-critical v, and Moon 1966 gives such a v.
- **klein's table.** The program reproduces klein's t_k(n) and u_k(n) for k <= 5 at orders 5-9 exactly.
- **Decomposition.**
  - At orders 7-9, every sole-exit-free class with hc <= 5 has lambda = 1.
  - It has exactly one minimal one-arc-leaving set D, with s_l(T - D) = 1 and e_u(D) = hc. The product law holds on
    every class.
- **No case with e_u(D) = 1.** A strong core with e_u(D) = 1 would be uniquely Hamiltonian with u as its only sole exit.
  THM-4868 (C) and (G) exclude this: u would have to be y-- with y- a sole entry, which makes u critical.
- **Observed core counts:**
  - c_3: 2, 3, 6, 12 at d = 5-8;
  - c_4: 2, 2, 4, 8 at d = 5-8;
  - c_5: 0, 4, 5, 10 at d = 5-8.
- **2-strong classes.** 2-strong classes with hc <= 5 occur only at orders 5 (RT_5, hc 2) and 6 (hc 4: one class; hc 5:
  four classes).
- **Cut vertex.** No 2-arc-strong class with a cut vertex has hc <= 5.

**Algebra.** With R(x) = x^2 (1 - 2x)/(1 - 3x + x^2) and C_k(x) = sum c_k(d) x^d:
- the decomposition gives U_k(x) = C_k(x) R(x)/x;
- doubling cores make C_k(x)(1 - 2x) a polynomial, which cancels the factor (1 - 2x) of R.

The resulting series:
- U_3 = x^6 (2 - x)/(1 - 3x + x^2), i.e. F(2n - 9);
- U_4 = 2x^6 (1 - x)/(1 - 3x + x^2), i.e. 2F(2n - 11), apart from the one 2-strong class at order 6;
- U_5 = x^7 (4 - 3x)/(1 - 3x + x^2).

## Pre-registered predictions

Run: `fewcycle 14 5`. Orders 10-14 are run after this commit; order 14 only if it finishes in reasonable time.

**Sanity checks.**
- **(P1)** t_1(n) = F(2n - 6): 377, 987, 2,584, 6,765, 17,711 at n = 10-14.
- **(P2)** t_2(n) = (n - 2) F(2n - 9) - L(2n - 10): 589, 1,775, 5,257, 15,360, 44,394 at n = 10-14.
- **(P3)** At n = 10, t_3, t_4, t_5 = 1,193, 1,049, 1,500 (klein's exhaustive values).

**Predictions.**
- **(P4)** u_1(n) = u_2(n) = 0 at n = 10-14 (HYP-9459 (1)).
- **(P5)** u_3(n) = F(2n - 9): 89, 233, 610, 1,597, 4,181 at n = 10-14 (HYP-9459 (2)).
- **(P6)** u_4(n) = 2F(2n - 11): 68, 178, 466, 1,220, 3,194 at n = 10-14 (HYP-9459 (3)).
- **(P7)** u_5(n) = 4F(2n - 12) - 3F(2n - 14) = L(2n - 12) + F(2n - 13): 60, 157, 411, 1,076, 2,817 at n = 10-14
  (new; klein's order-10 value is 60).
- **(P8) Decomposition.** At n = 10-14, every class counted in u_3, u_4 or u_5 has lambda = 1 and exactly one minimal
  one-arc-leaving set D, with s_l(T - D) = 1, e_u(D) = hc, and 0 product-law failures.
- **(P9) Core table.** At order n, the number of classes of type k whose core has order d is c_k(d) R(n - d + 1). The
  cores are:
  - c_3(5) = 2, and c_3(d) = 3 * 2^(d - 6) for d >= 6;
  - c_4(5) = 2, and c_4(d) = 2^(d - 5) for d >= 6;
  - c_5(5) = 0, c_5(6) = 4, and c_5(d) = 5 * 2^(d - 7) for d >= 7.

  The new core orders d = 9-13 appear as the d = n - 1 entries:
  - c_3 = 24, 48, 96, 192, 384;
  - c_4 = 16, 32, 64, 128, 256;
  - c_5 = 20, 40, 80, 160, 320.
- **(P10)** No 2-strong class and no 2-arc-strong class with a cut vertex has hc <= 5 at n = 10-14.

Every outcome is reported. A failure of (P8) or (P9) refutes the decomposition picture as stated, and is reported with
the class.

## Results

Run: `fewcycle 14 5`, 27 s. Output: `04-computation/experiments/fewcycle_sole_exit_free_20261104_14_5.out`.

| | n = 10 | 11 | 12 | 13 | 14 | prediction |
|---|---|---|---|---|---|---|
| t_1 | 377 | 987 | 2,584 | 6,765 | 17,711 | P1 holds |
| t_2 | 589 | 1,775 | 5,257 | 15,360 | 44,394 | P2 holds |
| t_3, t_4, t_5 | 1,193, 1,049, 1,500 | 3,574, 3,390, 4,664 | 10,537, 10,705, 14,142 | 30,676, 33,194, 42,082 | 88,400, 101,420, 123,414 | P3 holds (n = 10) |
| u_1, u_2 | 0, 0 | 0, 0 | 0, 0 | 0, 0 | 0, 0 | P4 holds |
| u_3 | 89 | 233 | 610 | 1,597 | 4,181 | P5 holds |
| u_4 | 68 | 178 | 466 | 1,220 | 3,194 | P6 holds |
| u_5 | 60 | 157 | 411 | 1,076 | 2,817 | P7 holds |

**(P8) holds.** Every class counted in u_3, u_4 or u_5 has lambda = 1 and exactly one minimal one-arc-leaving set,
with s_l = 1 and e_u = hc. There are 0 product-law failures, on 217, 568, 1,487, 3,893 and 10,192 classes at orders
10-14.

**(P9) holds.** Every core-table entry equals c_k(d) R(n - d + 1). The new cores, at d = 9-13, are:
- c_3 = 24, 48, 96, 192, 384;
- c_4 = 16, 32, 64, 128, 256;
- c_5 = 20, 40, 80, 160, 320.

**(P10) holds.** There is no 2-strong class and no 2-arc-strong class with a cut vertex with hc <= 5.

**After the run, the statements were proved for all n (THM-4869).**
- u_2(n) = 0 and u_3(n) = F(2n - 9) for n >= 6; u_4(n) = 2F(2n - 11) and the u_5 law for n >= 7.
- The proof goes through the unique terminal core, the block list (C3 and W4 for k = 3), and an induction excluding
  2-strong tournaments with hc <= 5 at n >= 7 (hc <= 3 at n = 6).
- The doubling of the cores is the free choice of the arcs between u and the chain vertices w_2, ..., w_p.
