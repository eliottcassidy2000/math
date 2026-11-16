# Pre-registered predictions (klein-2026-11-16-force)

Written to `prediction_force.txt` at 2026-11-16 13:50 MST (P-A–P-E) and before the strong-realiser census, the
gluing search and the rounding check (P-F–P-H).

| # | prediction | verdict |
|---|---|---|
| P-A | fewer than 10 counts `h ≤ 200` have no digit-1 realiser at order ≤ 8, all below 70; same for digit 2 | **refuted**: 47 counts lack a digit-1 realiser at order ≤ 8 (3, 15, 25, 27, 35, 41, 43, 53, 55, 57, 67, 75, 83, 89, 93, 99, …) and 43 lack digit 2 — digit vertices are scarce; the all-order criterion of THM-5260 (C) replaces the guess |
| P-B | the counts `h ≤ 100` with no charge-0 root at order ≤ 8 include 11, 29, 35, 47, 65; 13 and 41 have charge-0 roots | **mostly confirmed, one item refuted**: at orders ≤ 9 the counts without a charge-0 root include 11, 29, 35, 47 but **65 has one**; 13 and 41 have charge-0 roots |
| P-C | ready same-charge port moves at order 7 (≤ 2 added vertices) only at charges 0 (from `H = 1`) and 4, from fewer than 10 bases; at order 8 (1 added vertex) only from `H = 1`; none iterates | **refuted**: order 7 has ready same-charge moves at charges −2, 0, 4, 6 (bases of counts 15, 27, 71, 81 besides the trivial 1), order 8 with one vertex at charges −2, 0, 4, 6, 8, 16 including a charge-0 move 255 → 511 (PORTS78.md) |
| P-D | `B_t − e_t = 2E_t` for every rooted tournament of order ≤ 8 | **confirmed** (orders 4–9, no mismatch) |
| P-E | the Tribonacci row counts are `round(c τ^m)` for `3 ≤ m ≤ 40` | **confirmed** (`c = 0.435616389351097…`, holds from `m = 1`) |
| P-F | at most 6 primes `p < 125` have no strong realiser vertex with `2^k B_t = p + 1`, and 41 is among them | **refuted on the count, confirmed on 41**: ten primes (3, 41, 43, 53, 67, 83, 89, 101, 103, 113) have no tower-forced 3p + 2 at any order |
| P-G | the counts `h < 125` with no rooted doubling child at any order are 11, 19, 29, 31, 35 plus at most four more | **refuted**: twenty counts below 125 lack a rooted doubling child at any order (the five predicted plus 39, 43, 47, 53, 55, 61, 63, 67, 71, 73, 77, 85, 87, 101, 113) |
| P-H | an exact `H(U) = H(A) + H(B) − 1` gluing on a shared arc exists for pairs of order ≤ 4 beyond two triangles, and every such pair attains the path-splice bound | **weakly confirmed / refuted in spirit**: all additive cases up to (4,5) are two 3-cycles on the arc plus sources or sinks, every one attaining the bound; strong pairs of orders (5,5), (5,6) give none |
