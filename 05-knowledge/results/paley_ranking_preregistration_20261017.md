# Pre-registration: the class number decides whether the natural order ranks the Paley tournament optimally (HYP-9287)

2026-10-17, opus-2026-10-17-S40. Written and pushed **before** any run on primes above 131.

## Background (already observed, not part of the test)

For p = 3 (mod 4) prime, the natural order 0 < 1 < ... < p-1 of the Paley tournament P_p has forward-minus-backward
excess p h(-p) (THM-4666 in preparation; Dirichlet's class number formula). Write L = L(1, chi_(-p)) = pi h(-p) / sqrt(p).

* Exact maximum acyclic subgraph (subset DP, brute-force cross-check at 11): the natural order is optimal at p = 3, 7, 23
  and not optimal at p = 11, 19, 31 (MAS 35, 107, 285 against 33, 95, 279); at p = 43 the certified MAS 543
  (Chindelevitch-Harutyunyan, arXiv:2607.26690, Thm 7.1) exceeds the natural 473.
* Iterated insertion local search (paley_ranking_local_search_20261017.py, 30 restarts x 150 kicks), p <= 131: the
  natural order was not beaten at p = 3, 7, 23, 47, 71 (L = 1.81, 1.19, 1.97, 2.29, 2.61) and was beaten at every other
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
