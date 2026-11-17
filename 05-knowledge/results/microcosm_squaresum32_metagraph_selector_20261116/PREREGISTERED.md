# Pre-registration ledger — klein-2026-11-16-micro

Predictions in [prediction_micro.txt](prediction_micro.txt), written before the runs.

| # | prediction | outcome |
|---|---|---|
| M-1 | block sizes 3, 1, 6, 1, 21 triangular; least orders 1, 3, 4 / 5, 5, 5, 5, 6, 6 / 6–8 with 63 at 8; 7·3^k a count for k ≥ 2 (63, 189, 567, 1701) | **confirmed** |
| M-2 | near-join misses in `[2^19, 2^21)` at least 1.5× the average in the classes 7, 21 (mod 64); 7 (mod 8) more than 1 (mod 8) | **refuted** (1.01, 0.99; 1.00 vs 0.98: flat) |
| M-3 | bases ≤ 10 reach `X_10` with 3,000,000 ≤ X_10 < 2^25 | **confirmed** (X_10 = 19,190,825, first miss 19,190,827) |
| M-4 | square-sum paths for all n in 25–34; no cycle for n ≤ 31; 1–3 cycles at 32; cycles at 33, 34 | **confirmed** (one cycle at 32; 1 and 11 at 33, 34) |
| M-5 | greedy cover of the 1,200 new order-9 counts by < 60 order-8 images; best class ≥ 25 %; a count-661 class in the top five | first **confirmed** (31); second **refuted** (18.8 %); third **refuted** (rank ≈ 4,843) |
| M-6 | every new order-9 count is in the image of some first creator of order 8 | **confirmed** |

Refuted clauses are recorded as refuted; none was revised after the fact.
