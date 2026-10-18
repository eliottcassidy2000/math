# Pre-registration (S44, second round): q-stability and p-saturation of transfer orders

Written and committed BEFORE running `vir_transfer_orders_20261018.py` on the models listed under "Runs". The first
round (`PREREGISTRATION_orders_67_57.md`, commit 94fbe3b30b) PASSED; see HYP-9315.

**Notation.** For the minimal model `Vir_(p,q)` with `N = (p-1)(q-1)/2`:
- `T_(p,q)(L, i)` is the multiset of transfer orders of `ker π` at length `N + i` and Li degree `L`. This is the dict
  `{order: multiplicity}` printed by `vir_transfer_orders_20261018.py` at weight `D = 2(N+i) + L`.
- `T^(k)` is its order-`k` multiplicity.
- `W(N) = {(L, i) : L <= N + i - 1}` is the conservative stable window. It is read off from the existing data: (3,5)
  leaves the common table at (7,1), (4,5) at (7,1) and (3,7) at (9,2), all outside `W(N)`. Nothing outside `W(N)` is
  predicted.

**Observation that motivates the test (already seen, NOT a prediction).**
- `T_(5,6) = T_(5,7)` at every bidegree with `L + 2i <= 10`.
- `T_(6,7)` equals both, except for one extra order-4 class at `(L, 0)` for `L >= 9`.

**Predictions.**
- **P4 (q-stability).** For fixed `p`, `T_(p,q)(L,i) = T_(p,q')(L,i)` at every computed `(L, i)` lying in
  `W(N_(p,q))` ∩ `W(N_(p,q'))`. This is tested on all pairs within each family:
  - `{(3,8), (3,10), (3,11), (3,13)}`;
  - `{(4,7), (4,9), (4,11)}`;
  - `{(5,6), (5,7), (5,8), (5,9)}`.
- **P5 (HYP-9315(b), new models).** The length-`N` row is `T(L, 0) = {1: 1, 2: 1, ..., m: 1}` with
  `m = min(p-2, floor((L-1)/2))`, for every computed `L >= 3`.
- **P6 (saturation in p).** If `k <= p - 2` and `p < p'`, then `T^(k)_p(L,i) = T^(k)_(p')(L,i)` at every common
  bidegree inside both stable windows.
- **P7 (monotonicity in p).** If `p < p'`, then `T^(k)_p(L,i) <= T^(k)_(p')(L,i)` for every `k` at every common bidegree
  inside both stable windows.
- **P8.** No mismatches: the recorded orders sum to `dim ker π` at every bidegree.

**Runs** (single prime `2^31 - 1`; window `L + 2i <= D - 2N`):

| model | `N` | `D` |
|---|---|---|
| (3,11) | 10 | 32 |
| (3,13) | 12 | 36 |
| (4,9) | 12 | 36 |
| (5,7) | 12 | 36 |
| (5,8) | 14 | 40 |
| (6,7) | 15 | 41 |
| (4,11) | 15 | 42 |
| (5,9) | 16 | 42 |

The comparison data already on disk:
- (3,8) and (3,10) to `D = 30`;
- (4,7) to `D = 30`;
- (5,6) to `D = 30`;
- (5,7) to `D = 34`.

**Side remark (not a prediction).** Of these models, (3,8), (3,10), (3,13), (4,7), (4,9), (4,11), (5,6), (5,8) and
(6,7) have `pq` even. They are nodes of the Berggren tree of primitive Pythagorean triples under
`(p, q) -> (q^2 - p^2, 2pq, q^2 + p^2)`; see the S44 note.

## Outcome (recorded after the runs): ALL PASSED

- The pre-registration was committed in dcbb7ffa8a at 14:29:53 and pushed. The runs started at 14:30. All eight
  outputs were written between 14:30:11 and 14:31:24.
- The evaluation is mechanical: `evaluate_stability_prereg_20261018.py`, with outputs `evaluate_stability_prereg.out`
  (informative counts) and `evaluate_stability_prereg_full.out` (every pair).

**Results.**
- **P8 PASSED.** There are 0 MISMATCH lines in all 12 outputs.
- **P5 PASSED.** All eight new models have exactly the predicted length-N row.
- **P4 PASSED.** There are 0 failures in 14 pairs. Counting only bidegrees where at least one side is nonempty:
  - (5,7) vs (5,8): 19 nonempty, 13 with order >= 2;
  - (4,9) vs (4,11): 19 nonempty, 13 with order >= 2;
  - (3,11) vs (3,13): 17 nonempty, 5 with order >= 2.
- **P6 PASSED.** There are 0 failures. The tests with content are:
  - order 2 between p = 4 and p = 5, 6: 7-14 bidegrees per pair;
  - orders 2-3 between p = 5 and p = 6: 8-11 bidegrees per pair.
- **P7 PASSED.** There are 0 failures.

**The stable table.** Wherever both are computed, every p-family is identical. Examples for orders 1 and 2:
- `(6,1) {1:3, 2:1}`, `(7,1) {1:4, 2:3}`, `(8,1) {1:5, 2:4}` (order 3 added for p >= 5), `(8,2) {1:4, 2:1}`.

The order-(p-1) part is a proper subset of the saturated count. For example, (9,1) has `{3:1}` for p = 4 against
`{3:3}` for p = 5, 6, and (10,1) has `{3:2}` against `{3:5}`.

## Outcome AMENDED after the independent audit (s44_audit_20261018, M3): P4 and P6 REFUTED as stated

- The evaluation above compared only bidegrees computed in both tables. The (5,6) data stopped at weight 30
  (m = L + 2i <= 10).
- The audit ran (5,6) to weight 32, at two primes. S44 reproduced it: `orders_56_D32.out`.
- That run contradicts P4 and P6 inside the stated window `W(N)`:
  - (L,i) = (9,1): (5,6) gives {1:6, 2:6, 3:2}; (5,7) and (6,7) give {1:6, 2:6, 3:3}.
  - (L,i) = (10,1): (5,6) gives {1:7, 2:7, 3:4}; (5,7) gives {1:7, 2:7, 3:5}.
- The break is forced. The Verma map of (5,6) has its second relation at relative level `2q - 1 = 11`, and
  `dim ker π(2N+m) = p(m-3) - p(m-2p+1) - p(m-2q+1)`.
- So "ALL PASSED" was an artifact of the data range, and `W(N)` was the wrong window. It was read off the data, not
  derived from the BGG thresholds (MISTAKE-640).
- P5, P7 and P8 remain consistent with all data.

**Corrected statements.** These are post hoc, NOT pre-registered.
- (U) All models agree for `m < min(2p-1, 2p'-1)`.
- (Q) For the same `p`, models agree for `m < 2 min(q,q') - 1`.
- (S) Orders `k <= p_min - 2` agree for `m < 2 q_min - 1`.
- `evaluate_corrected_windows_20261018.py` finds 0 disagreements in all 16 S44 models: 388 + 452 bidegree comparisons
  and 2401 order counts. Every same-`p` pair first disagrees exactly at `m = 2 q_min - 1`. The audit's 45 pairs agree.
