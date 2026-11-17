# Pre-registration ledger — klein-2026-11-16-tooth

Predictions in [prediction_tooth.txt](prediction_tooth.txt), two batches, each written before the runs it covers.

| # | prediction | outcome |
|---|---|---|
| R-1 | quadratic extension lemma, 0 mismatches on orders ≤ 7; at most 60 distinct child counts per order-7 class | lemma **confirmed** (62,422 extensions); the bound **refuted** (one class produces 97) |
| R-2 | `𝔅(1), 𝔅(3), 𝔅(5)` have odd parts {1}, {3}, {1, 3, 5, 7, 9}; 7 and 21 are not `h + 2^k B` nor dead-end appends | **confirmed** |
| R-3 | mean distinct parents ≥ 0.8 n at orders 8, 9; the most fertile class is not transitive for n−1 ≥ 5; `M_{n−1}` has fewer children than classes for n ≥ 5 | **confirmed** (0.915, 0.961; 21 < 24, 48 < 58, 106 < 121, 227 < 248; 9 < 12, …) |
| R-4 | tree branching at 8 → 9 has maximum > 100; 10–50 % of order-8 classes childless | **confirmed** (227; 43 %) |
| R-5 | genealogy out-degree mode at 1 or n−2; ≥ 60 % of first creators at orders 6–9 have a sole-exit/sole-entry genealogy edge | **refuted** (mode at the middle degree; 82 %, 72 %, 55 %, 38.5 %) |
| R-6 | < 15 % of first creators at orders 6–9 have a digit-0 vertex; the 63 pair has none | order 6 **refuted** (21 %), orders 7–9 **confirmed** (7.3 %, 1.7 %, 0.37 %); 63 clause **refuted** (one of the pair has one) |
| R-7 | ≥ 50 % of first creators at orders 5–9 have a first-creator tree parent | **refuted** at orders 5–7 (17 %, 27 %, 44 %), **confirmed** at 8, 9 (61 %, 73 %) |

Refuted clauses are recorded as refuted; none was revised after the fact.
