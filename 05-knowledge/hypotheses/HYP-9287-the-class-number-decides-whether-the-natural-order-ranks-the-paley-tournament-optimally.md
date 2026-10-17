---
id: HYP-9287
title: "The class number decides whether the natural order ranks the Paley tournament optimally: the natural order of P_p (excess p h(-p), THM-4666) is a maximum acyclic ordering exactly when L(1, chi_(-p)) = pi h(-p)/sqrt(p) exceeds a slowly rising threshold (about 2 near p = 50, between 2.25 and 2.52 for 131 < p <= 500); in particular it is never optimal when h(-p) = 1 and p > 7"
status: >
  OPEN (CONJECTURE). Pre-registered test PASSED (paley_ranking_preregistration_20261017.md, pushed before the run):
  P1 nine of the ten primes with L >= 2.2 in (131, 500] unbeaten (only 439, L = 2.249, beaten; criterion: at most one);
  P2 all twenty primes with L <= 1.7 beaten by certified orders. FINITE-EXACT: optimality at p = 3, 7, 23 and
  non-optimality at 11, 19, 31 (exact MAS), 43 (cited MAS 543). "Unbeaten" is local-search evidence, not a proof of
  optimality; every "beaten" is a certified witness order.
session: opus-2026-10-17-S40
source: 05-knowledge/results/binary_information_class_groups_largest_anomalies_20261017.md
scripts:
  - 04-computation/experiments/paley_ranking_local_search_20261017.py (runs '131 30 150' and pre-registered '132:500 12 80'; + .out, .json)
  - 04-computation/experiments/paley_mas_exact_20261017.py (exact MAS, p <= 31)
  - 04-computation/experiments/paley_legendre_walk_positivity_20261017.py (L+ against L(1, chi))
related: [THM-4666, HYP-9286, THM-871, THM-4644, HYP-9277]
---

# HYP-9287 — class numbers and optimal rankings of Paley tournaments

## Statement

For primes `p ≡ 3 (mod 4)` write `L = L(1, χ_(-p)) = π h(-p)/sqrt(p)`. There is a threshold `τ(p)`, slowly increasing,
such that the natural order `0 < 1 < ... < p-1` is a maximum acyclic ordering of `P_p` iff `L >= τ(p)` (up to a thin
boundary layer). Necessary conditions are proved: the natural order must be insertion-locally optimal, i.e. `p ∈ L+`
(all partial sums of `(./p)` nonnegative, THM-4666 (5)); hence the set of such `p` has relative density 0
(Baker-Montgomery).

## Evidence

| p | h | L | natural excess | best found | natural optimal? |
|---|---|---|---|---|---|
| 3, 7 | 1 | 1.81, 1.19 | 1, 7 | 1, 7 (exact) | yes |
| 11, 19 | 1 | 0.95, 0.72 | 11, 19 | 15, 43 (exact) | no |
| 23 | 3 | 1.97 | 69 | 69 (exact) | yes |
| 31 | 3 | 1.69 | 93 | 105 (exact) | no |
| 43 | 1 | 0.48 | 43 | 183 (cited exact) | no |
| 47, 71 | 5, 7 | 2.29, 2.61 | 235, 497 | 235, 497 | unbeaten |
| 151, 199, 271, 439 | 7, 9, 11, 15 | 1.79, 2.00, 2.10, 2.25 | | beaten | no |
| 167, 191, 239, 263, 311, 359, 383, 431, 479 | 11-25 | 2.52-3.59 | | | unbeaten |

All other primes `p ≡ 3 (mod 4)` up to 500 (all with `L <= 1.77`) are beaten. The best excess found grows like
`0.75-0.96 p^(3/2)` at the beaten primes near 400-500, while the natural excess is `L p^(3/2)/π`; the observed threshold
matches `L/π ≈ 0.75-0.8`, i.e. the natural order wins only when it beats the generic optimum.

## Reading

The natural order is the arithmetic (affine) ranking: all affine orders have excess `±p h(-p)`. Its excess is the
area under the Legendre walk, and its local optimality is the positivity of that walk, the condition that also excludes
a real zero of `L(s, χ_(-p))`. The conjecture says the arithmetic ranking is globally optimal only when the class number
is large; for class number one it is never optimal (`p > 7`). Large `h(-p)` = small primes split = small steps go
forward. This is a tournament form of the `L(1, χ)` landscape; the conjecture would make `MAS(P_p) - p(p-1)/4`
equal to `p h(-p)/2` on a density-zero but infinite set (if `L(1, χ_(-p))` is unbounded along `p ≡ 3 (mod 4)`).

## Tests

* Exact MAS for `p = 47` or `71` (a symmetric branch and bound, or the ILP of Chindelevitch-Harutyunyan, who certified
  `p = 43`) would turn "unbeaten" into a proof or a refutation at the two smallest high-`L` primes.
* The threshold law: fit `τ(p)` on `500 < p <= 2000` (pre-register before running), in particular at the primes with
  `2.2 <= L <= 2.8`.
* Whether the best orders at beaten primes are perturbations of affine orders (structure of the witnesses in the JSON).
