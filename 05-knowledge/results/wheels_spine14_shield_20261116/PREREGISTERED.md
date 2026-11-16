# Pre-registered predictions (klein-2026-11-16-wheels2)

Written to `prediction_wheels2.txt` at 2026-11-15 23:43 MST. At that moment no spine-14 slice had printed anything; the
`address_moves.py` run for `n = 20` had already written its output file in the background but had not been read, and the
`n = 22` run was still computing; the factorisations of `Φ_608(2)` and `Φ_736(2)` were still running.

| # | prediction | verdict |
|---|---|---|
| P-S1 | spine 14 (second apex beating ≤ 3 spine vertices) attains every odd `N < 2^18` except 7 and 21 | **confirmed**, and improved: every odd `N < 2^19` except 7, 21 (SPINE14.md) |
| P-S2 | the first miss above 21 lies in `[500,000, 560,000]` (doubling 77,533 → 117,345 → 256,521) | **refuted**: the first miss above 21 is 729,421 (ratio 2.84 over 256,521, not a doubling) |
| P-S3 | below `2^18` fewer than 0.1 % of attained values need overlap exactly 3, and `A → B` suffices for all | **confirmed**: 70 of 131,070 values (0.053 %) need overlap 3; every value below `2^18` is attained with `A → B` (below `2^19`, 18 values need `B → A`) |
| P-A1 | at `n = 20, 22` the fraction of addresses with an exact +2 partner stays in [2.0 %, 2.6 %], and single-bit +2 partners exist at the critical bit | **first clause confirmed** (2.2 % at both widths); **second clause refuted**: the single-bit partners at the critical bit realise `+4` at `n = 20` (512 addresses) and `+8` at `n = 22` (1,024), not `+2` |
| P-A2 | the fraction of addresses with an interior +2 edit (avoiding the top and bottom `n//5` bits) stays below 0.5 % | **confirmed**: 0.13 % (`n = 20`), 0.14 % (`n = 22`) |
| P-C1 (added 00:12, before `collatz_tree.py` ran) | in count coordinates the odd `h ≤ 2^20` whose Collatz orbit meets 7 or 21 form a few 2-adic rays, all one-apex counts | **refuted**: 373,393 of the 500,000 odd `h ≤ 10^6` meet 7 or 21 (every orbit through 5 passes 10 = h⁻¹(21)), and 187,200 of them are not one-apex counts (CONDITIONAL.md §1) |
| P-F1 | `ν(608) < 8` and `ν(736) < 12`, so `a(19) = a(23) = 6` | **first half confirmed**: `ν(608) = 6`, `c(5, 19) = 17 < 19`, `a(19) = 6`; `ν(736) = 3` (**confirmed**: `Φ_736(2)` = three primes of order 736, factordb factorisation verified locally), `c(5, 23) = 14 < 23`, `a(23) = 6` (period 1472) |

Not predicted (proved or computed before or independently of these runs): Theorems A–D of WHEELS2.md, the values
`a(11) = a(13) = 5`, `a(17) = 6`, the parity-sheet verdicts for `N ≤ 150`, the Erdős identification, and the address census
at `n ≤ 18`.
