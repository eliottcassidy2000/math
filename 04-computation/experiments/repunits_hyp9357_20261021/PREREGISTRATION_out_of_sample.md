# Pre-registration (opus-2026-10-21-S51): the repunit rarity out of sample

This file is committed before `inv.py 2100` finishes. Nothing in the range `1100 < d < 2100` has been looked at.

## Training data

* `S(d)` is the number of strong tournament classes, of all orders, with exactly `d` Hamiltonian paths. The data are complete for `d < 1100` (an ear census by the audit's method from `audit_S50/b3_ear_census_625.py`, with the cap raised to 1100).
* The fitted model on `301 ≤ d ≤ 1098` (399 values, residual s.d. 0.100) is

  `log S(d) = cubic(log d) + 0.239[3|d] + 0.088[9|d] + 0.065[5|d] − 0.002[d prime] + 0.028[blocked] − 0.027[d ≡ 7 mod 8]`.

  The standard errors are 0.014, 0.018, 0.014, 0.014, 0.023 and 0.012.
* "Blocked" means: composite, with every factorization `a·b` (`a, b ≥ 3`) using 7 or 21.

## Working hypotheses

The repunit rarity is a coincidence of separate arithmetic biases of the class-level spectrum:

* **(i) A multiplicative / `C_3` bias.** Multiples of 3, 9 and 5 get extra classes. The cheapest pieces are `H(C_3) = 3`, `H(R_4) = 5` and `H(IC(2,2)) = 9`. Primality itself is neutral.
* **(ii) Gap shadows.** 63 = 3·21 = 9·7 cannot use its divisibility by 3 or 9, because 7 and 21 are not values.
* **(iii) A small-order 2-adic bias.** `H ≡ 7 (mod 8)` is depleted in the second octave of the order 6–8 spectra, and this decays with the order.

For odd `k`, `2^k − 1 ≡ 1 (mod 3)` while `2^k + 1 ≡ 0 (mod 3)`, which is where the comparison "31 vs 33, 127 vs 129, 511 vs 513" comes from.

## Predictions for `1101 ≤ d ≤ 2099`

These are evaluated by `S(d)` from `inv.py 2100`, which is complete for `d ≤ 2100`.

* **(R1)** The refitted coefficient of `[3|d]` lies in `[0.15, 0.40]`, and is the largest arithmetic coefficient.
* **(R2)** `|β(prime)| < 0.05`: there is no primality effect.
* **(R3)** `β(d ≡ 7 mod 8) ∈ [−0.05, +0.02]`: the 2-adic effect is small at this size.
* **(R4)** `S(2047) < S(2049)`. Here 2049 = 3·683, while 2047 = 23·89 ≡ 1 (mod 3).
* **(R5)** The base-4 repunit 1365 = 3·5·7·13 lies above the cubic-only smooth fit (positive residual).
* **(R6)** At the top-order level, the `[3|d]` coefficient restricted to order-14 classes (`f(14) = 1125`) is at least 0.4. Near the bottom of high orders the spectrum is chain-dominated, and chains carry factors 3, 5, 9.
* **(R7)** The sign rule at the class-count level, `N(n, h)` (all classes on `n` vertices, from the strong inventory below 1100 by sequence convolution; not yet computed):
  * for `n = 10, 11, 12` and `k ∈ {5, 6, 7, 9}`: `N(n, 2^k − 1) < N(n, 2^k + 1)`;
  * for `k ∈ {4, 8, 10}`: `N(n, 2^k − 1) > N(n, 2^k + 1)`.

Any failure is recorded as a failure.
