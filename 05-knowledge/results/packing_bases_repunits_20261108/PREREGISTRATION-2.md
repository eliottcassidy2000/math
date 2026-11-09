# Pre-registration 2: coverage of large alpha-sums by explicit products

opus-2026-11-08-S86. Written and pushed before `coverage_products.py` was written or run. The budget-80 tables of
PREREGISTRATION.md (P1-P5) are complete: they prove gaps. This addendum tests only attainment, which needs no
completeness. A value is covered if an explicit tournament realises it.

## Construction universe

H_b is multiplicative over order joins (THM-4940 (1)), so every product of strong pieces is realised by their join.
The strong pieces used are:
- (i) every strong class with alpha_1 <= 80 (`few_cycle_strong_80.txt`, with its full packing vector);
- (ii) the strong one-apex tournaments (intercyclic, alpha = (a, 0)) for every reduced word of length <= 23;
- (iii) the two intercyclic rays of Session 16, alpha_1 = 3 * 2^s and 9 * 2^s - 2.

## Predictions

The range is alpha-sums 80 < m <= 100,000.

- **P6 (ternary coverage).** In base 3, every alpha-sum m with 80 < m <= 100,000 is covered by a product of pieces
  (i)-(iii).
  - Together with P1 this would give: base 3 misses no value 1 + 3m with m <= 100,000 other than 10 and 31.
  - Falsified if any m in the range is uncovered.
  - An uncovered m would NOT prove a gap; it would only show the construction universe is too small.
- **P7 (binary control).** In base 2, every alpha-sum m with 80 < m <= 100,000 is covered, i.e. every odd H in
  (161, 200001]. This is consistent with the known all-odd coverage to 2^32 (THM-4875 (D)). It checks the pipeline.
- **P8 (bases 4-8).** In each base b in {4, ..., 8}, the products cover every alpha-sum m with 80 < m <= 20,000.
  - Basis: the exact tables showed at most one new gap per base on (44, 80], and the number of products grows quickly.
  - Falsified for a base if some m in the range is uncovered.
  - This would be evidence for HYP-9530 (B), not a proof.

Outcome: to be appended below without editing the predictions.
