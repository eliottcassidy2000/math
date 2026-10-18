# Pre-registration: the class number decides whether the natural order ranks the Paley tournament optimally (HYP-9287)

2026-10-17, opus-2026-10-17-S40. Written and pushed **before** any run on primes above 131.

## Background (already observed, not part of the test)

For p = 3 (mod 4) prime, the natural order 0 < 1 < ... < p-1 of the Paley tournament P_p has forward-minus-backward
excess p h(-p) (THM-4666 in preparation; Dirichlet's class number formula). Write L = L(1, chi_(-p)) = pi h(-p) / sqrt(p).

* Exact maximum acyclic subgraph (subset DP, brute-force cross-check at 11): the natural order is optimal at p = 3, 7, 23
  and not optimal at p = 11, 19, 31 (MAS 35, 107, 285 against 33, 95, 279); at p = 43 the certified MAS 543
  (Chindelevitch-Harutyunyan, arXiv:2607.26690, Thm 7.1) exceeds the natural 473.
* Iterated insertion local search (paley_ranking_local_search_20261017.py, 30 restarts x 150 kicks), p <= 131: the
  natural order was not beaten at p = 3, 7, 23, 47, 71 (L = 0.60, 1.19, 1.97, 2.29, 2.61; the value 1.81 printed for p = 3 omitted the unit factor w = 6, corrected after the audit) and was beaten at every other
  prime (L <= 1.77). The search reproduced every exact optimum (11, 19, 23, 31, 43).

## Predictions for 131 < p <= 500 (run with 12 restarts x 80 kicks; any found order is printed)

* **P1 (high L).** At each of the ten primes with L >= 2.2, namely 167, 191, 239, 263, 311, 359, 383, 431, 439, 479,
  the search does **not** find an order beating the natural one. (Evidence only: "not beaten" is not a proof of
  optimality.) Passes if at most one of the ten is beaten.
* **P2 (low L).** At each of the twenty primes with L <= 1.7, namely 139, 163, 179, 211, 223, 227, 251, 283, 307, 331,
  347, 367, 379, 419, 443, 463, 467, 487, 491, 499, the search **does** beat the natural order (a certified witness
  order). Passes if all twenty are beaten.
* **Exploratory, no prediction:** 151 (L = 1.79), 199 (L = 2.004), 271 (L = 2.099).

## What a pass or a failure means

A pass supports the hypothesis that optimality of the natural order is governed by L(1, chi_(-p)), with a threshold
near 1.8-2.0 at these sizes (it may drift with p). A P1 failure (two or more high-L primes beaten) refutes the
threshold form. A P2 failure (a low-L prime not beaten) would show either a weaker search or a non-monotone law; the
order found or not found is recorded either way. Nothing here bears on MAS values beyond the orders exhibited.

## Outcome (appended after the run, 2026-10-17)

Run: `python paley_ranking_local_search_20261017.py 132:500 12 80`
(`paley_ranking_local_search_p132_500_20261017.out` / `.json`).

* **P1 PASSED.** Unbeaten at 167, 191, 239, 263, 311, 359, 383, 431, 479; beaten at 439 (L = 2.249, the smallest L of
  the ten): one of ten, within the stated criterion.
* **P2 PASSED.** All twenty low-L primes beaten, each by a certified witness order (in the JSON).
* **Exploratory.** 151 (L = 1.79), 199 (L = 2.004) and 271 (L = 2.099) are beaten. On 131 < p <= 500 the largest
  beaten L is 2.249 and the smallest unbeaten L is 2.518: the threshold drifts upward from the small-p values (23 with
  L = 1.965 and 47 with L = 2.29 are optimal or unbeaten). HYP-9287 records the drifting-threshold form.

## Second pre-registration: 500 < p <= 1000 (written and pushed before this run)

Run: `python paley_ranking_local_search_20261017.py 501:1000 6 40` (fewer restarts and kicks than the first test, so
"unbeaten" is weaker evidence here; every "beaten" remains a certified witness order).

* **Q1 (high L).** At each of the seven primes with L >= 2.8, namely 503, 599, 647, 719, 839, 887, 911, the natural
  order is **not** beaten. Passes if at most one of the seven is beaten.
* **Q2 (low L).** At each of the 27 primes with L <= 2.0 (523, 547, 563, 571, 587, 607, 619, 631, 643, 659, 683, 691,
  727, 739, 751, 787, 811, 823, 827, 859, 883, 907, 919, 947, 967, 971, 991), the natural order **is** beaten. Passes
  if all 27 are beaten.
* **Exploratory, no prediction:** 743 (L = 2.42), 863 (L = 2.246), 983 (L = 2.705).

## Audit addendum to the first test (2026-10-17, after the independent S40 audit)

Certified rotation orders (rank `x` by `frac(alpha x)`) beat the natural order at 167 (1861), 263 (4079) and 383 (7519),
and at 439 with excess 9313. So P1 has 4 of 10 high-L primes beaten and, by the failure rule above, the threshold form
is **refuted**. The PASS recorded in the Outcome section only says that the pre-registered search did not find these
orders: it was too weak to support "unbeaten" (MISTAKE-626). P2 stands.

## Outcome of the second test (appended after the run and the audit)

Run: `python paley_ranking_local_search_20261017.py 501:1000 6 40`
(`paley_ranking_local_search_p501_1000_20261017.out` / `.json`).

* **Q2 PASSED.** All 27 primes with L <= 2.0 are beaten by certified witness orders.
* **Q1 FAILED.** As executed, the local search left all seven high-L primes unbeaten. But certified rotation orders beat
  503 (11657 against 10563), 599 (15507 against 14975), 647 (17661 against 14881) and 887 (28235 against 25723) (audit),
  and the exact Farey sweep (THM-4666 1c) also beats 911 (30327 against 28241). That is five of seven, against at most
  one allowed. The sweep shows that 719 and 839 keep the natural order as their best rotation order.
* **Exploratory.** 743 and 863 are beaten by the local search; 983 is beaten by rotation orders (audit: 34225 against
  26541; S40 scan: 32103).

HYP-9287 is rewritten accordingly: refuted in its threshold form, with a revised conjecture.
