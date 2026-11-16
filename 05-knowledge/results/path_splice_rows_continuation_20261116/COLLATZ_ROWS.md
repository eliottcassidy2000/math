# 3q + 2, 3q + 4, the two Collatz maps, and the odd numbers in rows of doubling chains

Session klein-2026-11-16-splice. Script `scripts/rows_collatz.py`, receipt `rows_collatz.out`.

## 1. The two hosts 3q + 2 and 3q + 4 of the digit programme

The ternary-digit programme (fresh_unit_prime_gluing_20261116: triple_splice, quadruple_splice, ternary_deposit) needs,
from a source of count `q`, hosts of counts `3q`, `3q + 2`, `3q + 4` — the three odd numbers in `[3q, 3q + 4]`, one in
each residue class modulo 3 — with a common continuation. In the lineage vocabulary these are the C2 children
`3q + 2σ` with digits `σ = 0, 1, 2` (THM-4862, THM-4865: `3q + 2σ = q + 2B_a` with `σ = B_a − q`, so the three hosts
are the height-two values at a vertex whose root response is `B_a = q, q + 1, q + 2`).

**The two famous maps are the digits 0 and −2.** In count coordinates `h = 2x + 1` (CONDITIONAL.md of
wheels_spine14_shield_20261116) the Collatz odd step `x ↦ 3x + 1` is `h ↦ 3h` and the `3x − 1` odd step is
`h ↦ 3h − 4`, both on `h ≡ 3 (mod 4)`; the even step `x ↦ x/2` is `h ↦ (h + 1)/2` on `h ≡ 1 (mod 4)` for both maps.
So, as maps on the source count `q`:

| map on `q` | digit `σ` in `3q + 2σ` | in half-counts `x = (q − 1)/2` | row of the image (§2) |
|---|---|---|---|
| `3q − 4` (the `3x − 1` step) | `−2` | `x ↦ 3x − 1` | — |
| `3q` (the Collatz step; digit 0) | `0` | `x ↦ 3x + 1` | `L_oddpart(3q − 1)` |
| `3q + 2` (digit 1) | `1` | `x ↦ 3x + 2` | `L_oddpart(3q + 1)` = the Syracuse image of `q` |
| `3q + 4` (digit 2) | `2` | `x ↦ 3x + 3 = 3(x + 1)` | `L_3·oddpart(q + 1)` |

Equivalently `3q + 2 = 3(q + 1) − 1` and `3q + 4 = 3(q + 1) + 1`: the two hosts are the `3y − 1` and `3y + 1` images
of the even number `y = q + 1`, i.e. of the half-count's successor; clone_clocks_negative_collatz_20261113 recorded the
same identity as `3H + 2 = T(2H + 1) = 2T(H) + 1` with the three digits `3r − 1, 3r, 3r + 1` in `r = (H + 1)/2`.
Nothing here makes the digit programme depend on either conjecture (CONDITIONAL.md §1 proved that Collatz cannot bear
on Chan–Pak); it says exactly which arithmetic step each host is.

## 2. Odd numbers in rows of chains, each growing by doubling

Every odd `N > 1` is `N = 1 + c·2^s` with `c` odd and `s ≥ 1` in exactly one way; call `L_c = {1 + c 2^s : s ≥ 1}` the
row of `c` (THM-5126 (B)'s 2-adic lines). The row starts at its head `1 + 2c ≡ 3 (mod 4)` and grows by doubling
`N ↦ 2N − 1` (the dyadic lift of the signed descent, THM-5262; `N − 1` doubles). The number 1 is on no row.

- **Collatz in count coordinates is the Syracuse map on row indices.** `(h + 1)/2` walks down a row to its head
  `1 + 2c`; the step `3h` at the head lands on `3 + 6c = 1 + 2(3c + 1)`, i.e. on the row of `oddpart(3c + 1)`, the
  Syracuse image of `c`. So the Collatz dynamics on the odd numbers, read in count coordinates, is: descend the row,
  apply Syracuse to its index. The trivial cycle is the row `L_1 = {3, 5, 9, 17, …}` (fixed index 1).
- **The forbidden pair.** `7 = 1 + 3·2 ∈ L_3` (the head) and `21 = 1 + 5·4 ∈ L_5` (second element), and Syracuse sends
  `3 → 5 → 1`: the Collatz orbit of 7 in count coordinates is `7 → 21 → 11 → 33 → 17 → 9 → 5 → 3` (CONDITIONAL.md), i.e.
  head of `L_3` → `L_5` (descend to its head 11) → `L_1`. The pair `{7, 21} = {m, 3m}` for the head `m = 7`.
- **The seven near-join exceptions** `{1, 7, 11, 13, 21, 47, 67}` of HYP-9453 (also the Fermat-block non-cuts, and the
  `hc ≤ 2` gaps 21, 67 of HYP-9472) sit as follows: `{7, 13} ⊂ L_3` (positions 1, 2), `{11, 21} ⊂ L_5` (positions 1, 2),
  `47` = head of `L_23`, `67` = head of `L_33`, `1` on no row. Pre-registered P-5 **confirmed**. So the exceptional set is
  the first two rungs of the two rows `L_3`, `L_5` joined by one Syracuse step, plus two further heads; every later
  rung of every row is realised (THM-4875 (D): all odd numbers to `2^32` except 7, 21 by one-apex joins; THM-5126 (B):
  the 65 rows with `c ≤ 129` entirely, by explicit rays).
- **Where the three children of `q` land.** `3q + 2 = 1 + (3q + 1)` lies on the row of the Syracuse image of `q`, at
  height `v_2(3q + 1)`; `3q + 4 = 1 + 3(q + 1)` on the row `L_3·oddpart(q+1)`; `3q = 1 + (3q − 1)` on the row of the
  `3x − 1`-Syracuse image `oddpart(3q − 1)`. Thus the digit-1 child of every odd `q` is exactly one Syracuse step away
  on the row grid: the `+2` digit is the Collatz move between rows, the doubling along a row is the `×2` recursion.

## 3. The two recursions in one object

`THM-4865`: for a rooted tournament `(T, t)` the tower `T_k` has `h(T_1) = h − e + B`, `h(T_2) = h + 2B`, then tribonacci;
HYP-9453's near-join is `h + 2^k B` with the same root response `B`. The row `L_c` through a realised count is the
near-join row of `(h, B) = (1, c)` — a transitive base — while the C2 child `3h + 2σ` is `h + 2B` at `B = h + σ`. So the
`×2` recursion (rows) and the `+2` recursion (digits) are the `k`-direction and the `B`-direction of one array
`h + 2^k B`; the Collatz step between rows is the shear `B ↦ B + 1` at `k = 1`, and the whole programme's question is
whether every odd `N ≥ 69` is reached by some `(h, B, k)` (HYP-9453) — the rows are proved full for `c ≤ 129`
(THM-5126) and the digits exact up to the continuation problem of CONTINUATION.md.
