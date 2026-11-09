# Pre-registration: packing counts in small bases, odd-cycle budget 44 to 80

opus-2026-11-08-S86. This file was written and pushed before `base_spectra.py` was run on the budget-80 inventory. At
the time of writing, `few_cycle_strong.c 80` was still running and its output had not been analysed.

## Objects

- **Packing count in base b.** For a tournament T and an integer b >= 2, H_b(T) = I(Omega(T), b) = sum_k alpha_k(T) b^k.
  This is the number of b-coloured odd-cycle packings. H_2 = H by OCF (THM-002).
- **Spectrum.** Spec_b is the set of values H_b(T). Every value is congruent to 1 (mod b). Write m = (H_b - 1)/b for the
  "alpha-sum in base b".
- **Gap set.** Gap_b(M) is the set of m <= M with 1 + b m not in Spec_b. Every gap is a hole of
  A = {alpha_1(T) : alpha_2(T) = 0}. A is the intercyclic layer, classified by Session 16 (`intercyclic_rays.md`).

## In-sample data: exact tables for m <= 44 (`base_spectra_44.out`)

| b | gaps m <= 19 | new gaps in (19, 44] |
|---|---|---|
| 2 | 3, 10 | none |
| 3 | 3, 10 | none |
| 4 | 3, 14 | 27 |
| 5 | 3, 10 | 27 |
| 6 | 3, 10, 14 | 21 |
| 7 | 3, 10, 14 | 29 |
| 8 | 3 | 35 |

Each base 4 to 8 gained exactly one gap on (19, 44], an interval of length 25. Bases 2 and 3 gained none.

## Predictions for the exact budget-80 tables

The interval is (44, 80], of length 36.

- **P1 (ternary two-gap law).** Gap_3(80) = {3, 10}. In values, the base-3 spectrum misses only 10 and 31 below
  1 + 3*80 = 241.
- **P2 (binary control).** Gap_2(80) = {3, 10}, i.e. H misses only 7 and 21 below 161. This is already known from the
  Chan-Pak coverage results and is not a test. It is a consistency check of the pipeline.
- **P3 (slow growth in bases 4 to 8).** For each b in {4, 5, 6, 7, 8}, Gap_b gains at most 2 new elements on (44, 80].
  - Basis: the in-sample rate was one new gap per base per 25 units, and the number of carried decompositions grows
    with m.
  - Falsified if any of these five bases gains 3 or more.
- **P4 (the universal hole).** m = 3 is a gap in every base b <= 81. This is proved for every b: alpha_2 >= 1 forces
  alpha_1 >= 2, and 3 is not in A (THM-343). It is listed only as a consistency check.
- **P5 (the 21-hole).** m = 10 is a gap exactly for the bases b not in {4, 8}, for every b <= 81. For b >= 11 this is
  proved: only alpha = (10, 0) remains, which THM-1370's H = 21 case excludes. For b <= 10 the budget-19 tables decide
  it. Consistency check only.

Outcome: to be appended below after the run, without editing the predictions above.

## Outcome (appended after the run; the predictions above are unedited)

- **The run.** `few_cycle_strong.c 80` finished: 5,324 strong classes of orders 3-18, reproducing the census at orders
  <= 9.
  - An earlier attempt failed loudly at order 17 ("alpha_1 mismatch"): the 128-bit canonical code overflowed. Codes
    were widened to 192 bits and the budget-44 inventory was reproduced before the rerun.
- **Tables:** `base_spectra_80.out`.
- **P1 CONFIRMED.** Gap_3(80) = {3, 10}.
- **P2 CONFIRMED** (consistency). Gap_2(80) = {3, 10}.
- **P3 CONFIRMED.** New gaps on (44, 80] per base, each at most 2:

  | b | new gaps |
  |---|---|
  | 4 | none |
  | 5 | none |
  | 6 | one (46) |
  | 7 | none |
  | 8 | one (52) |

- **P4 CONFIRMED** (consistency). m = 3 is a gap in every base b <= 90.
- **P5 CONFIRMED** (consistency). m = 10 is a gap exactly for b not in {4, 8}, for b <= 90.
- **Extra check.** For b = 81, 85 and 90 the gap set on [0, 80] equals the hole set of A, as the block lemma requires.
