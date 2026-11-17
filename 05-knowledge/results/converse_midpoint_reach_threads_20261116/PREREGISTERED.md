# Pre-registration ledger — klein-2026-11-16-conv2

Predictions in [prediction_conv2.txt](prediction_conv2.txt), three batches, each written before the runs it covers.

| # | prediction | outcome |
|---|---|---|
| Q-1 | bases ≤ 9 reach `X_9 ∈ [100000, 500000]` with the seven exceptions; dead-end appends same exceptions, `X_9 ≤ X_9^ext < 2X_9` | **confirmed** (`X_9 = 398,129`; the dead-end appends add nothing, `X_9^ext = X_9`) |
| Q-2 | midpoint values above the reach are missed at ≥ half the odd gap density | **confirmed** (ratios 1.01, 0.99, 0.87, 0.61 by window) |
| Q-3 | same-shape form `h″ + 2^k P″` for < 10 % of strong `D ≤ 7`; never for `C_3`, strong4 | percentage **refuted** (70 of 397 = 17.6 %); `C_3`, strong4 **confirmed** (and proved at every order) |
| Q-4 | at least half but not all of the 62 strong-source pairs have a base induced inside `S` | **confirmed** (38 of 62) |
| Q-5 | `SC(n) = 1,1,2,2,6,8,44,72,736`; ratio maximal at `n = 4, 5` (both 2); strong uniquely maximal at 5 | values **refuted** (A002785 is 2, 2, 8, 12, 88, 176, 2752 for `n = 3..9`); peak location **refuted** as stated (maximal at `n = 5` alone, ratio 4); strong **refuted** (infinite through 5) |
| Q-6 | twin pair within 100 of `10^k` for ≤ 5 of `k ≤ 15`; `(7, 9)` only at `k = 4, 9` | count **refuted** (9 of 15); offsets **confirmed** |
| Q-7 | `21 | 2^{p−1} − 1 ⟺ p ≡ 1 (mod 6)` for primes `p ≤ 100` | **confirmed** (theorem) |
| Q-8 | `63·2^{s+2} − 11` through 64,501 near-join values from bases ≤ 9; 6015/6017/6019 from bases ≤ 8; no fixed `(h, B)` line | **confirmed** (all from bases ≤ 8, through 129,013) |
| Q-9 | one-vertex-core values contain 241 but not 493, 997, 2005 | **confirmed** |
| Q-10 | base `D + w` rewrites ≥ 80 % but not all of the 43 small sources | **confirmed** (36 of 43) |
| Q-11 | with `w` a sole exit, fewer than half | **confirmed** (13 of 43) |
| Q-12 | halved partner for fewer than half of 397; none for `C_3`, strong4; `strong4 → C_3` partners `(15, 216)` | count **confirmed** (53 of 397); `C_3` clause **refuted** (`C_3 → M_1` is a partner — the `k = 2` near-join in disguise); strong4 **confirmed**; `(15, 216)` **confirmed** (least partner order 6) |

Refuted clauses are recorded as refuted; none was revised after the fact.
