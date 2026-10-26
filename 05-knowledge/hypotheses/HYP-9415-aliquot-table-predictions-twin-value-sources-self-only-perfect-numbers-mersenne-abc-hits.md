---
id: HYP-9415
title: "After THM-4805: (1) twin collision values P + 1 in (10^4, 10^5] have exactly the two twin sources with share in [0.40, 0.60] and mean source count in [2.5, 3.2] (pre-registered; CONFIRMED 0.4848, 2.816); (2) 33550336 is not its own sole source (pre-registered; CONFIRMED, companion 44655764) and 28 is the only value v with s^-1(v) = {v}; (3) abc hits with c <= 10^5: odd-perfect-power share in [0.65, 0.85] (CONFIRMED 0.7709), power-of-two even-term share in [0.30, 0.50] (REFUTED: 0.2840), a = 1 share in [0.15, 0.30] (CONFIRMED 0.1504); replacement prediction for c <= 10^6: power-of-two share in [0.18, 0.28], a = 1 share in [0.09, 0.15]; (4) every value v with s^-1(v) = {v} is 28; (5) no prime k is a Mersenne abc hit (every 2^p - 1 with p prime is squarefree; classical OPEN)."
status: "NUMERICAL: (1) CONFIRMED, (2) first clause CONFIRMED, (3) two ranges CONFIRMED and one REFUTED, all pre-registered in the note before the C++ run (aliquot_sources_preregistered_20261026.out); the 10^6 replacement in (3) is untested. (4) OPEN (FINITE-EXACT for the perfect numbers to 33550336: only 28 is self-only). (5) OPEN classical (a prime hit k would need a Wieferich prime p with ord_p(2) = k prime; 1093 and 3511 have orders 364 and 1755)."
session: opus-2026-10-26-S60
source: 05-knowledge/results/aliquot_table_untouchable_primes_chains_seven_twentyone_20261026.md
related:
  - 01-canon/theorems/THM-4805-aliquot-sources-of-2-5-7-21-twin-prime-collisions-the-powers-and-doublings-chain-structure-and-abc-hits-on-the-row-of-one.md
  - 05-knowledge/results/session69_prime_ports_power_aliquot_20261026/README.md
  - 05-knowledge/results/aliquot_power_124_20261026/README.md
  - 05-knowledge/results/s57_untouchables_root27_20261026/README.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
---

# HYP-9415: after the aliquot table

1. **Twin-value sources.** For twin primes `P - 2, P` the value `P + 1` has the two sources `2(P - 2)` and `P^2` (session69's twin identity; THM-4805 B' (i) says these are the only row-type sources) and possibly others of shapes `2qr`, `2^a qr`, `2q^2 r`, .... Pre-registered for `(10^4, 10^5]`: exactly-two share in `[0.40, 0.60]`, mean in `[2.5, 3.2]`. Observed `0.4848`, `2.816` (max 11). Next: the share at `10^6`, predicted in `[0.42, 0.52]`.
2. **Self-only perfect numbers.** `s^-1(28) = {28}` while 6, 496, 8128 and 33550336 have companions (25; 652; 10682, 10828; 44655764 = 4 x 1669 x 6689). Pre-registered: 33550336 is not self-only; confirmed. Conjecture (4): 28 is the only `v` with `s^-1(v) = {v}`. Falsified by a perfect number with no companion, or any non-perfect fixed point.
3. **abc shares.** Pre-registered for `c <= 10^5`: odd perfect power present in `[0.65, 0.85]` (observed 0.7709), power-of-two even term in `[0.30, 0.50]` (observed 0.2840, **refuted**; the share was 0.40 at `10^4`), `a = 1` in `[0.15, 0.30]` (observed 0.1504). Replacement, untested: at `c <= 10^6` the power-of-two share lies in `[0.18, 0.28]` and the `a = 1` share in `[0.09, 0.15]`.
4. See 2.
5. **Mersenne abc hits.** `(1, 2^k - 1, 2^k)` is an abc hit iff `2^k - 1` is not squarefree (THM-4805 D'). The hit set is the union of the multiples of `ord_(p^2)(2)` over odd primes `p`; a prime `k` is a hit iff some prime `p` has `p^2 | 2^k - 1`, i.e. `p` is Wieferich with `ord_p(2) = k`. Conjecture (classical, equivalent to the squarefreeness of all Mersenne numbers with prime exponent): there is none.
