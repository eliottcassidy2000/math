---
id: HYP-9478
title: "Pre-registered out-of-sample tests for the representation mass of one-arc join families: (P1) the one-apex x bank family F2 has mean multiplicity exactly M_2 = sum_y (2/y_h) ln(1 + y_h/y_R) = 48.1976 in the dyadic blocks k = 28, 29; (P2) the one-apex x one-apex family F1 gains (pi^2/12) ln 2 = 0.5701 per doubling at k = 28, 29; (P3) F2 alone, with the fixed order-<=7 bank, attains every odd N <= 2^30 except 7 and 21; (P4) its lower tail keeps thinning; (P5) F1 u F2 attains every odd N <= 2^35 except 7 and 21; (P6) F2 alone attains every odd N <= 2^35 except 7 and 21; (P7) F1-only miss ratios stay in [1.50, 1.65] for k = 32..34"
status: PRE-REGISTERED (committed before any run beyond 2^28)
session: opus-2026-11-06-S82
date: 2026-11-06
related: [THM-4885, THM-4875, HYP-9468]
script: 04-computation/experiments/chanpak_f1f2_cover34_20261106.c; 04-computation/experiments/chanpak_join_multiplicity_20261106.c
---

# HYP-9478 - representation mass of the one-apex join families; F2 alone; coverage to 2^35

## Setting

Notation of THM-4875 (D) (audited; recertified by Session 203):
- g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2), and e(v) = [v = 0 or v even];
- F1 = {(1 + u)(1 + v) + g(u)g(v) - e(u)e(v) : 0 <= u <= v};
- F2 = {y_h g(v) - y_e e(v) + y_R (1 + v) : v >= 0, y in the 2,732 distinct rooted triples of order <= 7}.

For odd N, r1(N) and r2(N) count representations: pairs u <= v for F1, and pairs (y, v) for F2.

**Mass formulas** (derived before this commit; proof to be written in THM-4885). For a fixed triple (a, b, c) with
a, c >= 1:

    #{v >= 0 : a g(v) - b e(v) + c (1 + v) <= X} = (X / a) ln(1 + a / c) + o(X).

Why: g(v) = v rho(v) + O(log v), where rho(v) is the complement of the bits below the top bit, read from the least
significant end as a binary fraction. rho is uniform and asymptotically independent of v / 2^n within a dyadic block.

Consequences:
- The mean of r2 over odd N tends to M_2 = sum_y (2 / y_h) ln(1 + y_h / y_R) = 48.19761 for the order-<=7 bank.
- The mean of r1 over the odd N in a block [2^k, 2^(k+1)) grows by (pi^2/12) ln 2 = 0.57009 per doubling, since
  integral_0^1 ln(1 + t)/t dt = pi^2/12.

## In-sample evidence (before this commit; k <= 27)

- **Mean r2 by block.** 48.1960 (k = 16); 48.1977 (k = 19); 48.1976 for k = 24..27.
- **Mean r1 increments.** 0.56 to 0.58 per doubling for k = 18..27.
- **F2 alone.** It misses only 7 and 21 up to 2^28. The block minimum of r2 is 5 to 8 for k = 18..27.
- **Lower tail** of r2, as the fraction of odd N in block k:

  | k | r2 <= 10 | r2 <= 12 | r2 <= 20 |
  |---|---|---|---|
  | 22 | 1.00e-5 | 4.15e-5 | 4.53e-3 |
  | 27 | 4.87e-6 | 2.85e-5 | 3.42e-3 |

- **F1-only misses** match THM-4875 (D) exactly (65, 103, ..., 63,546 for k = 10..25).

## Predictions (out of sample)

- **(P1)** In the blocks k = 28 and k = 29, the mean of r2 lies in [48.1966, 48.1986].
- **(P2)** In the blocks k = 28 and k = 29, the mean of r1 exceeds the previous block's mean by an amount in
  [0.555, 0.585].
- **(P3)** F2 alone attains every odd N <= 2^30 except 7 and 21. The minimum of r2 in each of the blocks k = 28, 29 is
  at least 4.
- **(P4)** In block k = 29:
  - the fraction with r2 <= 12 is below 2.85e-5;
  - the fraction with r2 <= 20 is below 3.42e-3.
- **(P5)** F1 u F2 attains every odd N <= 2^35 = 34,359,738,368 except 7 and 21. That is eight times the current
  certified range of 2^32.
- **(P6)** F2 alone attains every odd N <= 2^35 except 7 and 21.
  - Confidence is about 75%. The extrapolated lower tail puts P(r2 = 0) between 1e-14 and 1e-12 per number.
  - Any miss refutes P6, but not P5. It would mean that a fixed finite bank stops working.
- **(P7)** The F1-only miss counts m_k in (2^k, 2^(k+1)) for k = 32, 33, 34 satisfy m_k / m_(k-1) in [1.50, 1.65]
  (m_31 = 934,107 is the value recorded in THM-4875 (D)).

If the 2^35 run cannot finish, P5 to P7 are reported as far as the computation reached and marked UNDECIDED beyond it,
not as confirmed.

## Results

(To be filled after the runs.)
