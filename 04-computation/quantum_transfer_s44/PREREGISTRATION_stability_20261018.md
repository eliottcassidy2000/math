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
