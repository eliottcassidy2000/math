# Pre-registered predictions (klein-2026-11-16-nearjoin)

Written to `prediction_nearjoin.txt` at 2026-11-16 15:07 MST, before any run of this session.

| # | prediction | verdict |
|---|---|---|
| P-1 | `B_z(T′) = H(T′)` for every near-join of a base of order ≤ 6, `k ≤ 3`; `z` carries the full count when the base is strong | **confirmed**: 1,260 cases, no failure; the 69 cases where `z`'s component is smaller all have a non-strong base |
| P-2 | `σ_u` = number of backward-junction two-path covers; `σ_z = e_{t_0}(S − z)` for sole exits | **confirmed**: 3,610 and 380 checks, no failure |
| P-3 | fewer than half of the strong digit-0 pairs `(S, u)` have `H(S) = h(S − u) + 2^k B_t(S − u)` | **confirmed**: 23 of 55 (orders 4–7) |
| P-4 | 11, 13, 47, 67 are never near-join values from any base | **confirmed** (THM-5263 (B): the strong census of orders ≤ 9 has no digit-0 vertex at these primes, and the near-join's `z` always is one) |

Not predicted: the out-degree-1 law for digit-0 vertices of strong tournaments (found in the data), the dead-end-append
coverage, and the reader computations of READERS.md.
