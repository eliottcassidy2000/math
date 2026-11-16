# Pre-registered predictions (klein-2026-11-16-splice)

Written to `prediction_splice.txt` at 2026-11-16 09:30 MST, after the order-4..6 path-splice check had passed and
before the order-7 check, the continuation searches and the row census ran.

| # | prediction | verdict |
|---|---|---|
| P-1 | path-splice bound at order 7: 0 violations, shuffle count never below the bound | **confirmed** (4,143,216 instances; SPLICE.md §4) |
| P-2 | equality `H(U) = H(A)H(B)` in more than 15 % of the `|W| ≤ 1` instances at order 7 | **confirmed**: 479,878 of 2,419,536 = 19.8 % |
| P-3 | weakening the continuation to `(X, g)` lets a source of order ≤ 8 carry the three digits `3H, 3H+2, 3H+4` with common `(3X, 3g)` | **refuted as stated**: no rooted tournament of order ≤ 6 (two added vertices, exhaustive) and none of the three literal sources (orders 6–8) has the three digits under tripled `(X, g)`; every tripled-`(X, g)` hit was already a tripled-`(E, O, g)` hit. **Confirmed in the modified form** of CONTINUATION.md §3: with the root made non-terminal (a sink appended) the parity split vanishes and the three digits of the count-5 sources (15, 17, 19) share the complete continuation `(40, 40, 9)` |
| P-4 | from such a triple a second fresh step with a common continuation exists at order ≤ 10 | **refuted** at two added vertices from each of the three children (CONTINUATION.md §4): the counts 3H′ + d′ occur but never with a common `(X, g)`, let alone `(E, O, g)` |
| P-5 | rows of the seven exceptions: `{7, 13} ⊂ L_3`, `{11, 21} ⊂ L_5`, `47 ∈ L_23`, `67 ∈ L_33`, 1 on no row; Syracuse `3 → 5 → 1` | **confirmed** (COLLATZ_ROWS.md §2, receipt `rows_collatz.out`) |

Not predicted: the theorem itself (proved before the file was written), the `D + g ≥ 0` law (noticed in the data), the
order-7 one-vertex scan.
