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

## Outcome (appended after the run; the predictions above are unedited)

`coverage_products.py` writes `coverage_products.out`, in about 5 s. Within m <= 80 the product closure reproduces
the exact gap tables in every base 2-8, which checks the pipeline.

- **P7 (binary control) REFUTED.** 18,038 alpha-sums in (80, 100000] are uncovered, starting with m = 179 (H = 359).
  Every one of them is an H-value (all-odd coverage to 2^32, THM-4875 (D)). So the construction universe, not the
  spectrum, is the limit.
- **P6 REFUTED.** 24,059 uncovered alpha-sums in base 3, starting with m = 340. No gap is implied.
- **P8 REFUTED.** In bases 4-8 there are 5,182 to 6,561 uncovered alpha-sums in (80, 20000]. No gap is implied.
- **Why.** Mid-range values need single strong pieces that have disjoint cycles and more than 80 odd cycles. Neither
  the alpha_1 <= 80 inventory nor the intercyclic families supply them, and products of small pieces grow too fast.
- **The analogue that should have been run first.** Pieces with alpha_1 <= 20 (plus the intercyclic families) leave
  14 alpha-sums in (20, 80] uncovered, in both bases 2 and 3. Pieces with alpha_1 <= 40 leave one, m = 80, in each
  of bases 2 and 3. That
  would have predicted the failure. The prediction's basis ("the number of products grows quickly") was not checked
  against it. This is recorded with the session's MISTAKE entry.
- **Status of the ternary law.** HYP-9530 (A) rests on the exact tables to alpha-sum 80 (P1) and nothing more.

**Addendum (after 1d57e6b23c; not part of the pre-registered outcome).**
- An exploratory, non-pre-registered constructive check, `ternary_attainment.py`, attains every value = 1 (mod 3) up
  to 5398 except 10 and 31. It does so already with tournaments of order <= 9, and the audit reproduced it
  independently.
- With THM-4940 (3) and (4) this makes Spec_3 ∩ [1, 5398] exact (THM-4940 (9)).
- Beyond 5398, HYP-9530 (A) remains OPEN, and P6 stands REFUTED as registered.
- The audit follow-up found a truncation bug in the first version of `ternary_attainment.py`: inventory vectors were
  cut to alpha_1..alpha_4. It created 23 fake piece values, all above the reported ranges. After the fix, the base-3
  uncovered count to m = 10,000 is 1701, not 1698; the thresholds and witnesses are unchanged.
