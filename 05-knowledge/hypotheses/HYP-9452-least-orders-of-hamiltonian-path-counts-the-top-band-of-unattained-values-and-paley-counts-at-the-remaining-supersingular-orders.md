---
id: HYP-9452
title: "After THM-4861: (1) for n >= 8 every odd h <= 0.85 * max H(n) other than 7 and 21 is attained at order n (the unattained odd values below the maximum form 7, 21 and one top band: from 611 of 661 at n = 8, from 2,883 of 3,357 at n = 9); (2) hence the least order n(h) of a count h equals the least n with max H(n) >= h, except in the top band (true for the five Fermat primes); (3) the Paley counts H(P_q) are known for q <= 31 (and GF(27)) but not at the supersingular orders 47, 59, 71 nor at the non-supersingular 43 (HYP-9357), 67, 79, 83, and no formula is known (CONJ-002 is refuted), so the supersingular/non-supersingular contrast at Paley orders cannot be tested"
status: >
  OPEN (klein-2026-11-02). (1)-(2) are read off the census of orders 6-9 (THM-4861 (b); the order-9 gap list is already in tournament_pseudoprimes_20261021/analyze_census.out) and A038375; the band's
  lower edge 611/661 = 0.92 and 2883/3357 = 0.86 suggests a constant fraction, but order 7 has an isolated gap at 63
  (0.33), so the claim is stated from n = 8. (3) is a statement of reach: q = 31 was done in paley_h31_exact_rate_and_apery_wheel_20261014.md; q = 43 needs
  2^42 * 43 cells by the plain DP; an orbit-counting method under the affine group of order q(q-1)/2 (HYP-9273 B: the
  orbit count is 1 mod 8) would reduce it.
session: klein-2026-11-02 (codex/session-ramsey-20261030)
source: 05-knowledge/results/fermat_monster_hamiltonian_counts_prime_classes_20261102/README.md
related:
  - 01-canon/theorems/THM-4861-the-five-fermat-primes-are-hamiltonian-path-counts-of-tournaments-with-least-orders-3-4-6-8-11-a-vertex-transitive-tournament-has-a-prime-count-only-for-c3-and-fermat-primes-are-free-vertices-of-the-prime-square-sum-graph.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
  - 01-canon/theorems/CONJ-002-paley-ham-paths.md
---

# HYP-9452 — least orders of Hamiltonian path counts, the top band, and the unreachable Paley orders

## 1. The spectrum at a fixed order is an initial segment up to a top band (OPEN)

THM-1370 proves that 7 and 21 are never attained and that the spectrum is monotone in the order; that every other odd number is attained is the Chan–Pak conjecture (OPEN, MISTAKE-559). The census of THM-4861 (b) shows the finer shape at a fixed order: at `n = 8` the unattained odd values below the maximum 661 are 7, 21 and `611, 615, 617, 619, 623, 625, 635, 647, 655`; at `n = 9` they are 7, 21 and a band from 2,883 to 3,355 with few exceptions. **Claim:** for `n ≥ 8` every odd `h ≤ 0.85 · max H(n)` other than 7 and 21 is attained at order `n`. (Order 7 violates such a bound with its isolated gap at 63; order 6 has the gaps 35, 39 of 45.)

## 2. Least orders (OPEN as a general law; PROVED for the Fermat primes)

Define `n(h)` as the least order with `H = h`. By monotonicity `n(h) ≥ min{n : max H(n) ≥ h}`; the claim is equality unless `h` lies in the top band of that order. For the Fermat primes equality holds: `n(3, 5, 17, 257, 65537) = 3, 4, 6, 8, 11` (THM-4861 (a)). A test: the least order of 611 (the first gap at order 8) is 9, and `max H(8) = 661 ≥ 611` — so the band exceptions are real and the law must carry the band.

## 3. Paley counts at the remaining supersingular orders (reach, not a conjecture)

`H(P_q)` is known for `q = 3, 7, 11, 19, 23` (THM-4747, CONJ-002, THM-336; recomputed in THM-4861) and for `q = 31` (`paley_h31_exact_rate_and_apery_wheel_20261014.md`). The supersingular Paley orders `47, 59, 71` and the non-supersingular `43, 67, 79, 83` need `2^{q−1} · q` 64-bit cells by the plain subset DP. CONJ-002's formula `|Aut| · 3^{(q−3)/2}` is refuted, and the quotients `H/|Aut| = 1, 9, 1729, 6857869865, 62293308207033` show no pattern beyond oddness. Until `H(P_31)` or an orbit-counting method exists, the question "do supersingular Paley orders behave differently" has no data.
