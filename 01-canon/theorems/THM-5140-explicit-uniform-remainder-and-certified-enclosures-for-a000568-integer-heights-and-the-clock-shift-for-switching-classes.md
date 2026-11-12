---
id: THM-5140
title: "The tournament census at giant orders. (B1) Explicit remainder of the level expansion of A000568, uniform in n and in the truncation order: with f_n = n! a(n)/2^C(n,2) = sum_k L_k(n), 0 <= f_n - sum_{k <= K} L_k(n) <= B_K(n) = sum_{k = K+1}^{floor(n/3)} C(n,2)^(2k) 2^(3k^2+k-2kn)/(2k)! + n! 2^(-ceil(n(n-1)/3)), log-convex in k, of the true exponential order 4^(-(K+1)n). Certified enclosures of a(n) by ball arithmetic for n up to 10^1000. (B2) a(n) is the expected cycle index of S_n at random powers of two 2^(H_p), H_p = sum_{d | p} J_d with independent integer-valued discrete Gaussians J_d (P(J_1 = j) ~ 2^(-j(j+1)/2), E[2^(M J_1)] = 2^C(M,2)); this field and the Gaussian field of THM-4845 have the same moments and mutually singular laws. (B3) Switching classes: A049313(n) = N_odd(n) + N_lev(n), N_odd(n) n!/2^C(n-1,2) = sum_k P_k(n) 4^(-k(n-1)) with the tournament polynomials P_k (clock n - 1), N_lev = 0 for odd n and at most 2^(n^2/4) for even n; first correction (64/3)(n)_3 4^(-n)"
status: >
  PROVED (elementary): (B1) from the bit-loss lemma kept by level (D >= 2kn - 3k^2 - k for defect 2k and 3k <= n,
  D >= n(n-1)/3 otherwise) and the count e_d(1, ..., n-1) <= C(n,2)^d/d! of permutations of defect d; (B2) by
  completing the square in the theta sums, with the exponent formula of THM-516; (B3) from the fixed-class count
  2^(o - l + 1) of THM-479 (Babai-Cameron) and a count of pair orbits. Source note, sections 3.2, 3.4, 3.6.
  Machine checks: the level sums add up to f_n exactly for 3 <= n <= 60; the bound in all 395 cases 6 <= n <= 60,
  n = 80, 100, 150, 200, 250, K <= 6; the enclosure with one level contains the exact a(n) and returns its 31
  leading digits for n = 50, 100, 150, 200, 250; the moments of the integer heights to 40 digits and a(n) for
  n <= 15; (B3) exactly for 3 <= n <= 19 against the OEIS values of A049313, N_odd = A049313 at every odd n <= 19.
  Enclosures (Arb): n = 300, 10^3, 10^4, 10^5, 10^6, 10^9, 10^12, 10^100, 10^1000 (exponent and 31 leading digits).
  KNOWN, credited: the level expansion, its polynomials P_k and the pole structure (Session98,
  ALL_POLES_CORE_EXPANSION.md, priority; THM-4845 (D)-(E); correction terms: Moon 1968, Wright 1970); the bit-loss
  lemma by support and the bound 0 <= f_n - 1 <= q^3/(1-q), q = n 2^(-(n-1)/3) (Session98,
  LIMITS_AND_NORMALIZATION.md section 2); the support envelope U_M(n), fixed-order asymptotics and 16-digit
  forecasts to n = 10^1000 (Session271, which lists a growing-order theorem and certified digits as open); the
  Gaussian representation (THM-4845 (A)); the theta shape of the exponent (THM-516, HYP-2538); the cycle gas and
  its Mayer expansion (April reflections, Session249), the factor 64 (Session271); A049313 = N_odd + N_lev with
  the fixed-class counts (THM-479); that exact integers cannot be obtained by rounding a truncation (Session98,
  Session270); log-normal moment indeterminacy (Stieltjes).
  NEW to the repository (after the audit's corrections): the explicit level-wise bound B_K(n), of the true
  exponential order, and its log-convexity (a bound by support valid for all n and all truncations is
  Session271's); outward-rounded, interval-certified digits; the application of the classical log-normal
  indeterminacy to the cycle index of A000568 (integer heights); the clock-shift rewriting of N_odd and the
  first correction (64/3)(n)_3 4^(-n) of A049313 with its remainder (N_lev = 0 at odd n and its values are
  THM-479; the bound on N_lev is immediate from it; the leading term is presumably folklore).
  NOT CLAIMED: any new asymptotic order (the orders are those of Session98 / THM-4845); exact integers.
  AUDIT: PASS WITH CORRECTIONS (05-knowledge/results/switching_quadratic_census_gauge_20261111/audit/AUDIT.md):
  (B1) checked on all 8,656 odd cycle types with 3 <= n <= 40 and in 437 cases of the bound; the enclosures
  reproduced by another code path (mpmath) at every n of the table, the ball computation reviewed; (B2), (B3)
  proofs correct, N_odd reproduced at odd n <= 13. Corrections of credit applied: Session271 already has a
  remainder bound by support valid for all n and all truncations (only the level-wise form and outward-rounded
  digits were open); N_lev = 0 at odd n and the values of N_lev are THM-479 part 3; moment indeterminacy with
  discrete theta solutions is classical. MISTAKE-766.
session: mac-mini-2026-11-11-gauge
source: 05-knowledge/results/switching_quadratic_census_gauge_20261111/README.md
related:
  - THM-4845
  - THM-514
  - THM-516
  - THM-479
  - THM-5139
  - HYP-2538
---

# THM-5140 - the census at giant orders: uniform remainder, certified enclosures, integer heights, switching classes

Source: [the note](../../05-knowledge/results/switching_quadratic_census_gauge_20261111/README.md), section 3.
Independent audit: [AUDIT.md](../../05-knowledge/results/switching_quadratic_census_gauge_20261111/audit/AUDIT.md), PASS WITH
CORRECTIONS (no mathematical error; three novelty claims narrowed).

## Statement

`a(n) = A000568(n)`, `f_n = n! a(n)/2^C(n,2) = sum_k L_k(n)`, `L_k(n) = P_k(n) 4^(-kn)` (level `k` = cycle types with odd
parts at least 3 and defect `2k`).

- **(B1)** `0 <= f_n - sum_{k <= K} L_k(n) <= B_K(n)` for all `n >= 3`, `K >= 0`, with
  `B_K(n) = sum_{k = K+1}^{floor(n/3)} t_k(n) + n! 2^(-ceil(n(n-1)/3))`, `t_k(n) = C(n,2)^(2k) 2^(3k^2+k-2kn)/(2k)!`;
  `t_k` is log-convex in `k`, so `B_K(n) <= (n/3) max(t_{K+1}, t_{floor(n/3)}) + n! 2^(-ceil(n(n-1)/3))`.
  Consequently `a(n) = (2^C(n,2)/n!)(1 + sum_{k <= K} L_k(n) + theta B_K(n))`, `0 <= theta <= 1`; table of certified
  values in the note, section 3.3.
- **(B2)** With independent integer-valued `J_d` (`d` odd), `P(J_1 = j)` proportional to `2^(-j(j+1)/2)`,
  `P(J_d = j)` proportional to `2^(-j^2/(2 phi(d)))` for `d >= 3`, and `H_p = sum_{d | p} J_d`:
  `E[2^(M J_1)] = 2^C(M,2)`, `E[2^(M J_d)] = 2^(phi(d) M^2/2)`, and `a(n) = E[Z_n(2^(H_p))]` where
  `sum_n Z_n x^n = exp(sum_{p odd} W_p x^p/p)`. The Gaussian field of THM-4845 has the same moments.
- **(B3)** `N_odd(n) = (2^C(n-1,2)/n!) sum_k P_k(n) 4^(-k(n-1))`; `N_lev(n) = 0` (`n` odd), `0 <= N_lev(n) <= 2^(n^2/4)`
  (`n` even); `A049313(n) n!/2^C(n-1,2) = 1 + (64/3)(n)_3 4^(-n) + ...`.

## Proof

Sections 3.2, 3.4, 3.6 of the note.
