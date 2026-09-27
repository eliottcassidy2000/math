---
id: THM-4512
title: "Coefficient-descent classes are certified up to at most one member: for a Syracuse valuation word w of length j with first coefficient descent at j (3^j < 2^A, A = v_1+...+v_j), every odd n in the residue class of w modulo 2^A with n > N(w) = S_j/(2^A - 3^j) satisfies U^j(n) < n; N(w) < 2^A for every j <= 5000 (and for all j given an effective irrationality measure for log_2 3), so the only possibly uncertified member of a class is its representative. No class with j <= 14 has an uncertified member except the class of n = 1, and sigma(n) = sigma_inf(n) for every odd 3 <= n <= 10^7."
status: >
  PROVED (elementary: U^j(n) = (3^j n + S_j)/2^A with S_j = sum_(t<j) 3^(j-1-t) 2^(A_t) > 0,
  S_j <= 2^A ((3/2)^j - 1), Terras's inequality sigma >= sigma_inf) + FINITE-EXACT
  (2^A - 3^j > (3/2)^j - 1 at the minimal A for j <= 5000, worst ratio 0.507 at
  (j, A) = (5, 8); 606746 first-descent classes with j <= 14 and their representatives
  rho_w = -S_j 3^-j mod 2^A checked, only rho = 1 uncertified; sigma = sigma_inf for all
  odd 3 <= n <= 10^7, maximal sigma 155). Terras's conjecture sigma = sigma_inf for all
  n >= 2 remains open beyond these ranges; this theorem sharpens 'finitely many exceptions
  per class' to 'at most the representative'. No literature-priority claim.
source: opus-2026-09-26 session gilbreath6-collatz-precision-20260926
depends_on: [Terras 1976 (coefficient stopping time; the set of n with stopping time k is a union of residue classes mod 2^k), Lagarias 1985 survey section on stopping times, THM-4495 (exact order of the no-descent counts), reset_20260926_swaplift.md (the lane's certificate bank, for comparison)]
verification: 04-computation/experiments/collatz_precision_residual_20260926.py -> .out (residual densities D(k), thresholds N(w), exception enumeration for j <= 14); collatz_coefficient_stopping_20260926.py -> .out (sigma versus sigma_inf to 10^7; the gap inequality to j = 5000)
---

# THM-4512 -- coefficient-descent classes have at most one uncertified member

## Setting

`U(n) = (3n + 1)/2^v` on odd `n`. For a valuation word `w = (v_1, ..., v_j)`
put `A_t = v_1 + ... + v_t`, `A = A_j`. The odd `n` whose first `j`
valuations are `w` form one residue class modulo `2^A` (Terras), with
representative `rho_w = -S_j 3^(-j) mod 2^A`, where

```text
U^j(n) = (3^j n + S_j) / 2^A,      S_j = sum_(t=0)^(j-1) 3^(j-1-t) 2^(A_t)  > 0   (A_0 = 0).
```

*Coefficient descent* at step `j`: `3^j < 2^A`. *Actual descent*: `U^j(n) < n`.
`sigma_inf(n)` and `sigma(n)` are the first `j` with each property.

## Statement

1. Actual descent at `j` implies coefficient descent at `j` (`S_j > 0`), so
   `sigma(n) >= sigma_inf(n)` (Terras).
2. If `3^j < 2^A`, then `U^j(n) < n` for every `n > N(w) := S_j/(2^A - 3^j)`.
3. `S_j <= 2^A ((3/2)^j - 1)`, hence `N(w)/2^A <= ((3/2)^j - 1)/(2^A - 3^j)`;
   the right side is `< 1` for every `j <= 5000` at the minimal admissible `A`
   (hence at every `A`), with maximum `0.507` at `(j, A) = (5, 8)`, so each
   coefficient-descent class contains at most one member not certified by
   its coefficient descent, namely its representative `rho_w`, and only if
   `rho_w <= N(w)`. For `j` beyond `5000` the same follows from any effective
   irrationality measure for `log_2 3` (`2^A - 3^j >= 3^j j^(-mu)` beats
   `(3/2)^j` once `2^j > j^mu`).
4. For `j <= 14` the only class with an uncertified member is the class of
   `n = 1` (word `(2)`, `N = 1`, `U(1) = 1`). For every odd `3 <= n <= 10^7`,
   `sigma(n) = sigma_inf(n)` (maximal value `155`).

## Proof

`U^j(n) < n` iff `3^j n + S_j < 2^A n` iff `(2^A - 3^j) n > S_j`; if
`2^A <= 3^j` this fails, giving 1; if `2^A > 3^j` it holds exactly for
`n > N(w)`, giving 2. Since every valuation is `>= 1`, `A_t <= A - (j - t)`,
so `S_j <= 2^A sum_(m=1)^j 3^(m-1) 2^(-m) = 2^A ((3/2)^j - 1)`, giving 3; the
inequality `2^A - 3^j > (3/2)^j - 1` was checked exactly for `j <= 5000` with
`A = ceil(j log_2 3)` (the minimal admissible `A`; larger `A` only increase the
gap). For 4, the classes with `j <= 14` were enumerated for all valuations up
to `40` beyond the critical one; beyond that range `S_j/2^A <= 2^(-v) *
2((3/2)^j - 1) < 1 <= rho_w`, so no exception is possible; the direct
comparison to `10^7` is the script's part 1. ∎

## Remarks

* The obstruction the lane calls "insufficient precision" is therefore not
  the small members of certified classes (there is at most one, and none
  below `10^7` except `n = 1`): it is the words *without* coefficient descent,
  i.e. the no-descent set, whose density among odd integers is `D(k)` after
  `k` Syracuse steps (`0.0027` at `k = 41`; `0.00062` at `k = 60`) and whose
  exact order in `T`-coding is THM-4495.
* Terras's conjecture (`sigma = sigma_inf` for all `n >= 2`) would follow from
  `rho_w > N(w)` for every first-descent word; the theorem reduces it to a
  statement about representatives only.
