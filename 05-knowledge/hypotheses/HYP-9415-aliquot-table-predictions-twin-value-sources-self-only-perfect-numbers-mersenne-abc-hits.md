---
id: HYP-9415
title: "After THM-4805: (1) twin collision values p + 1 in (10^4, 10^5] have exactly the two twin sources with share in [0.40, 0.60] and mean source count in [2.5, 3.2]; (2) 28 is the only perfect number that is its own sole aliquot source (33550336 has a further source); (3) abc hits with c <= 10^5: share with an odd perfect power in [0.65, 0.85], share whose even term is a power of two in [0.30, 0.50], share with a = 1 in [0.15, 0.30]; (4) every value v with s^-1(v) = {v} is 28; (5) no prime k is a Mersenne abc hit, i.e. every 2^p - 1 with p prime is squarefree (classical OPEN)."
status: "NUMERICAL: (1)-(3) pre-registered in the note before the C++ run; outcomes recorded in the note section 8a. (4) OPEN (FINITE-EXACT for the perfect numbers to 8128). (5) OPEN classical (equivalent to the squarefreeness of Mersenne numbers with prime exponent; Wieferich primes 1093, 3511 do not supply a counterexample since 364 and 1755 are composite)."
session: opus-2026-10-26-S60
source: 05-knowledge/results/aliquot_table_untouchable_primes_chains_seven_twentyone_20261026.md
related:
  - 01-canon/theorems/THM-4805-aliquot-sources-of-2-5-7-21-twin-prime-collisions-the-powers-and-doublings-chain-structure-and-abc-hits-on-the-row-of-one.md
  - 05-knowledge/results/s57_untouchables_root27_20261026/README.md
  - 05-knowledge/hypotheses/HYP-9400-untouchable-census-predictions.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
---

# HYP-9415: after the aliquot table

1. **Twin-value sources.** For twin primes `p - 2, p`, the value `p + 1` has the two sources `2(p - 2)` and `p^2` (THM-4805 B) and possibly others of shapes `2qr`, `2^a qr`, `2q^2 r`, ... Prediction (pre-registered): in `(10^4, 10^5]` the share of twin values with exactly the two twin sources lies in `[0.40, 0.60]` and the mean number of sources in `[2.5, 3.2]` (`[0.56, 2.68]` below `10^4`).
2. **Self-only perfect numbers.** `s^-1(28) = {28}` while 6, 496 and 8128 have companions. Prediction (pre-registered): `33550336` has at least one source other than itself. Conjecture (4): 28 is the only `v` with `s^-1(v) = {v}`.
3. **abc shares.** Pre-registered for `c <= 10^5`: odd perfect power present in `[0.65, 0.85]` of the hits (0.78 at `10^4`), power-of-two even term in `[0.30, 0.50]` (0.40), `a = 1` in `[0.15, 0.30]` (0.25).
4. See 2.
5. **Mersenne abc hits.** `(1, 2^k - 1, 2^k)` is an abc hit iff `2^k - 1` is not squarefree (THM-4805 D). The hit set is the union of the multiples of `ord_(p^2)(2)` over odd primes `p`; a prime `k` would be a hit iff `p^2 | 2^k - 1` for a prime `p` with `ord_(p^2)(2) = k`, which requires `p` Wieferich with `ord_p(2) = k` prime. Conjecture (classical): there is none.

Falsification: (1)-(3) by the printed C++ output; (4) by any perfect number with a single source beyond 8128, or any non-perfect fixed point; (5) by a Wieferich prime `p` whose `ord_p(2)` is prime.
