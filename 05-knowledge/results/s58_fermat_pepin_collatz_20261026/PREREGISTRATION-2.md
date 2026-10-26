# S58 pre-registration 2: the algebraic-factorisation correction, out of sample

opus-2026-10-26-S58. Committed before `collatz_line_primes_20261026.py test2` was run.

## Outcome of pre-registration 1

| id | prediction | observed | outcome |
|---|---|---|---|
| L-P1 | Fermat-type lines `m = 8..1200`, window `[2583.3, 2830.9]` | 2548 | **FAIL** (z = -3.2) |
| L-P2 | Mersenne-type lines `p = 8..400`, window `[1012.3, 1162.9]` | 1034 | pass (z = -1.8) |
| L-P3 | genuine Fermat lines `m = 2^e`, `e = 3..12`, window `[11.8, 34.0]` | 17 | pass (z = -1.3) |

All three came in low.

## Explanation (post hoc)

Two-parameter values carry algebraic factorisations that a small-prime sieve cannot see.
- `3^j 2^b + 1 = w^r + 1` whenever an odd `r` divides `gcd(j, b) = gcd(j, m)`.
- `3^j 2^b + 1 = 4u^4 + 1` (Aurifeuille) whenever `j = 0 (mod 4)` and `m = 2 (mod 4)`.
- `3^k 2^c - 1 = w^r - 1` whenever `r = gcd(k, p) > 1`.

Removing these values gives corrected expectations:
- **L-P1:** 2619.67 ± 48.69, so the observed 2548 has z = -1.47;
- **L-P2:** 1019.22 ± 29.16, so the observed 1034 has z = +0.51;
- **L-P3:** 22.91, unchanged; no algebraic value occurs on a genuine Fermat line.

Genuine Fermat lines (`m = 2^e`) and `F_m` itself are immune: `gcd(j, 2^e)` has no odd factor, and `2^e` is not
`2 (mod 4)` for `e >= 2`.

## Out-of-sample test of the corrected model (new ranges; same sieve, `Q = 2^16`; values `>= 2^32`)

| id | family | corrected E (sd) | window (±2.5 sd) | uncorrected E |
|---|---|---|---|---|
| L-P4 | Fermat-type lines `m = 1201..1500` | 671.98 (25.61) | `[608.0, 736.0]` | 697.27 |
| L-P5 | Mersenne-type lines `p = 401..520` | 324.06 (17.55) | `[280.2, 367.9]` | 348.07 |

The two models differ by about one sd in each family, so these windows test the calibration of the corrected model;
they cannot by themselves decide between corrected and uncorrected. The discriminating evidence is the combined
post-hoc fit above, recorded as such.
