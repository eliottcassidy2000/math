# Pre-registration ledger — klein-2026-11-16-seam

Predictions in [prediction_seam.txt](prediction_seam.txt), two batches, each written before the runs it covers.

| # | prediction | outcome |
|---|---|---|
| S-1 | contact seams (bases ≤ 8) miss 6, 8, 20, 22 and at most 10 other even centres in [4, 2000] | hole clause **confirmed**; count **refuted** (16 others) |
| S-2 | arc-flip lemma, 0 mismatches; `Σ M(a,b) = (n−1)h`; flip seams (orders ≤ 8) cover [4, 1000] but ≤ 20; flip digits cover [−30, 30] | lemma and identity **confirmed**; coverage **refuted** (flips at order ≤ 8 reach only centres ≤ 660); digits **confirmed** ([−153, 153]) |
| S-3 | the +2 window is a contact seam for ≥ 80 % of odd counts H ≤ 300 | **confirmed** (141 of 148) |
| S-4 | ≥ 20 % of first-witness seams change the number of strong components | **refuted** (2 of 980) |
| S-5 | contact digits cover [0, 100]; > 5 % of contact flips have digit 1 | range **confirmed** ([0, 399] at ≤ 8); share **refuted** (0.74 %) |
| S-6 | order 10: contacts in [9.4e7, 9.7e7], mean parents ≥ 9.7; max children ≥ 490; `M_9` not max but largest branch; childless 45–65 % | **confirmed** except the childless share (41.0 %, **refuted**) |
| S-7 | bases ≤ 9 miss only 6, 8, 20, 22 below 2000; ≤ 40 other misses below 10,000 | first clause **refuted** (13 others persist, all below 263); second **confirmed** |
| S-8 | no missed centre below 10,000 is `≡ 2 (mod 8)` | **confirmed** |
| S-9 | of the seven window failures at bases ≤ 8, all but H = 1 are recovered at order 9 | **refuted** (only 103 recovered; 11, 17, 19, 27, 35 remain) |
| S-10 | heaviest-parent tree: > 80 % of first creators hang on first creators; max branching at 8 → 9 below 150 | share **confirmed** (100 % but 63); branching **refuted** (247) |

Refuted clauses are recorded as refuted; none was revised after the fact.
