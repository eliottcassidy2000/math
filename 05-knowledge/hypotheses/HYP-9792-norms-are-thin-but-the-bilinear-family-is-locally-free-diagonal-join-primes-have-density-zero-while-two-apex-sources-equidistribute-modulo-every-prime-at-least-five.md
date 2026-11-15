---
id: HYP-9792
title: "Norms are thin, the bilinear family is locally free: (i) the primes that are diagonal-join counts N(beta_A) = (H + P2)^2 + 2H^2 have relative density zero among the primes = 3 (mod 8), by an archimedean angle obstruction on top of the mod-8 reciprocity obstruction; (ii) over the explicit two-apex sources of order m + 2, the pair (H, Q) mod l tends to the uniform distribution on F_l^2 for every prime l >= 5, geometrically in m. So a prime-count route through Z[sqrt -2] must be bilinear (Re(beta_A conj beta_B), a Type II form) rather than a norm route"
status: >
  OPEN (CONJECTURE), mac-mini-2026-11-15-synthesis, wave 38 of chanpak_waves_20261113.
  (i) FINITE-EXACT: the census is complete below 4,190,211, since blocks of order >= 11 have P2 >= 2^11 - 2. There,
  26 of the 73,940 primes = 3 (mod 8) are diagonal joins: 3, 11, 227, 35963, 43427, 56891, ... The realised fraction is
  0.0096, 0.0024, 0.0002, 0.0003 in the decades from 10^4 to 4.2*10^6. HEURISTIC: the angle argument below.
  (ii) FINITE-EXACT through order 12: all 261,888 two-apex sources, total variation 0.005, 0.007, 0.011, 0.013 from
  uniform for l = 5, 7, 11, 13. TV roughly halves per spine vertex from order 8, and at order 12 it is near the
  sampling floor.
  Concurrent and credited: klein-2026-11-15-synthesis (chanpak_mod8_synthesis_span3_20261115/SYNTHESIS.md s.5) names
  the mod-8 reciprocity obstruction of the diagonal sub-family and the local freedom of the full family (via the
  residue densities of chanpak_generic_lines_llt_20261115/LLT.md s.1). This hypothesis adds the archimedean
  obstruction, the exact census and the explicit two-apex family.
session: mac-mini-2026-11-15-synthesis
date: 2026-11-15
source: 05-knowledge/results/chanpak_waves_20261113/README.md (wave 38), GRAND_SYNTHESIS.md
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/w38_norm_bridges.py (census (6))
  - 05-knowledge/results/chanpak_waves_20261113/scripts/w38_two_apex_equidistribution.py
related:
  - THM-5242, THM-5222, THM-5172 (C), chanpak_pair_line_census_20261115/TWO_APEX_COARSE_PLACEMENT.md, HYP-9772
---

# HYP-9792 - norms are thin, the bilinear family is locally free

**Setting.**
- For a tournament A put x = H(A) + P₂(A) = Q(A) − H(A), y = H(A) and β_A = x + y√−2.
- A diagonal join cyc(A, A′, ·) counts N(β_A) = x² + 2y² ≡ 3 (mod 8) (THM-5172 (C), THM-5242 (B)).
- A point-cyclic join cyc(A, B, ·) counts Re(β_A β̄_B) = x_A x_B + 2y_A y_B.

**(i) Diagonal joins are a density-zero source of primes.**

A prime p ≡ 3 (mod 8) has a unique representation p = x² + 2y² with x, y > 0. Cornacchia finds it. So p is a
diagonal join iff the single pair (H, P₂) = (y, x − y) is realised.

The angle obstruction (HEURISTIC, with exact inputs):
1. *x/y grows linearly with the order.* The ratio x/y = Q/H − 1 is at least about 2n − 5 at order n. The exact
   minima for n = 4, …, 10 are 5.40, 6.20, 8.51, 9.67 (attained by the Paley tournament P₇), 11.98, 13.51 and
   15.75. The sampled minima at n = 11–16 are 2n − 4.45 to 2n − 3.93. The bulk sits near 2n − 0.8: the mean of P₂/H
   is 2(n − 1) + 0.2, and the conditional sd of P₂/H falls from 0.435 at n = 10 to 0.229 at n = 16.
2. *Size forces a large order.* p = x² + 2y² ≤ (x/y)²·H_max(n)²·(1 + o(1)) forces n ≥ n(p), with
   n(p) ≍ log p / log log p.
3. *Hence a shrinking angle.* A reachable prime has representation angle arctan(√2·y/x) ≤ arctan(√2/(2n(p) − 5)).
   That is an angular share O(log log p / log p).
4. *Conclusion.* By Hecke's equidistribution of prime angles in ℤ[√−2], this is a vanishing proportion of the primes
   ≡ 3 (mod 8). Within the reachable sector the pair spectrum is itself sparse (fibre fill 3–4% in two dimensions at
   order 10, wave 36).

**(ii) The explicit two-apex family is locally free.**

The family: a spine 0 < … < m − 1, apexes A → B, A beating spine vertex i iff a_i = 1, B iff b_i = 1, with a₀ = 1,
a_(m−1) = 0, b ≠ 0. The conjecture: for every prime ℓ ≥ 5,

    TV( (H, Q) mod ℓ over the family of order m + 2 ,  uniform on F_ℓ² )  ≤  C_ℓ c^m

for some c < 1.

| m (order) | ℓ = 5 | ℓ = 7 | ℓ = 11 | ℓ = 13 |
|---|---:|---:|---:|---:|
| 6 (8) | 0.111 | 0.139 | 0.249 | 0.243 (4 empty cells) |
| 7 (9) | 0.044 | 0.062 | 0.114 | 0.112 |
| 8 (10) | 0.030 | 0.029 | 0.041 | 0.060 |
| 9 (11) | 0.018 | 0.014 | 0.022 | 0.025 |
| 10 (12) | 0.005 | 0.007 | 0.011 | 0.013 |

At ℓ = 2 and 3 the family is constrained by the proved laws: Q ≡ (−1)^n H (mod 3), Q ≡ 0 (mod 4), and Q ≡ 0 (mod 8)
at even order.

**Why it matters.**
- (i) and (ii) together say the norm (diagonal) route to prime counts is dead. The bilinear route survives.
- The point-cyclic join is a bilinear form over ℤ[√−2], and (ii) is the residue-class input a bilinear sieve needs, as
  in Friedlander–Iwaniec (X² + Y⁴) and Heath-Brown (x³ + 2y³).
- The two-apex family alone is dimension-critical: about 4^m/4 sources against a value range of about 6·4^m. So it
  cannot cover every target by itself. Coverage must come from combining families (GRAND_SYNTHESIS.md §6).

**Predictions.**
- *P103.* For every n, (H + P₂)/H ≥ 2n − 5 over all n-tournaments. It is exact for 4 ≤ n ≤ 10. Local search at
  n = 11–13 finds no tournament below 2n − 5.
- *P104.* For the two-apex family of order 14 (m = 12), TV ≤ 0.004 for ℓ = 5, 7, 11, 13, and no cell is empty
  for ℓ ≤ 31.
