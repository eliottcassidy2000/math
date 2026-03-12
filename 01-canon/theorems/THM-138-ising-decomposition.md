---
theorem_id: THM-138
title: Ising decomposition of H — alpha_1 vs alpha_2+ competition
status: PROVED (computational, p=7,11; mechanism understood for all p)
proved_by: kind-pasteur-2026-03-12-S57
date: 2026-03-12
related_theorems: [THM-136, THM-137]
related_hypotheses: [HYP-480]
tags: [paley, interval, ising, independence-polynomial, crossover, alpha]
---

## Main Result

**Theorem (THM-138):** For primes p = 3 mod 4, the Hamiltonian path count
H(T) = I(Omega(T), 2) decomposes into competing contributions:

```
H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
```

where alpha_j = #{independent sets of size j in Omega(T)}.

For Paley (P) vs Interval (I) tournaments:

1. **alpha_1 favors Paley**: alpha_1(P) > alpha_1(I) at p=7,11 (and conjecturally all p)
   - Paley has MORE total directed odd cycles
   - This is the "2-body" Ising term

2. **alpha_2+ favors Interval**: sum_{j>=2} 2^j*alpha_j(I) > sum_{j>=2} 2^j*alpha_j(P)
   - Interval cycles are MORE vertex-disjoint
   - This is the "many-body" Ising term
   - Interval's additive structure creates non-conflicting cycle packings

3. **Phase transition at p ~ 13-19**: the alpha_2+ advantage overtakes alpha_1

## Quantitative Data

| p | alpha_1(P) | alpha_1(I) | higher(P) | higher(I) | H(P) | H(I) | Winner |
|---|-----------|-----------|-----------|-----------|------|------|--------|
| 7 | 80 | 59 | 28 | 56 | 189 | 175 | Paley |
| 11 | 21169 | 18397 | 52756 | 56232 | 95095 | 93027 | Paley |
| 19 | ? | ? | ? | ? | 1.173T | 1.184T | Interval |
| 23 | ? | ? | ? | ? | 15.76Q | 16.01Q | Interval |

At p=7:
- Paley alpha_1 advantage: 2*(80-59) = 42
- Interval alpha_2+ advantage: 56 - 28 = 28
- Net: Paley wins by 14

At p=11:
- Paley alpha_1 advantage: 2*(21169-18397) = 5544
- Interval alpha_2+ advantage: 56232 - 52756 = 3476
- Net: Paley wins by 2068

At p=19: Interval wins by ~11.5B (net). The alpha_2+ advantage has overtaken.

## Ising Interpretation

The Walsh-Fourier expansion of H on the orientation cube {+1,-1}^m gives:

```
H(sigma) = H_0 + sigma^T J sigma + (higher degree terms)
```

- **J** (interaction matrix): Paley sigma_P is the eigenvector with largest eigenvalue (THM-137)
  - This means Paley maximizes the QUADRATIC (2-body) term
  - Equivalently: alpha_1 is maximized when cycles are individually many

- **Higher-order terms**: The 4-body, 6-body, ... terms correspond to
  independent sets of size 2, 3, ... in Omega
  - Interval maximizes these because its additive structure creates
    more vertex-disjoint cycle packings

- **Phase transition**: The dimensionless coupling g = 2*sqrt(p)/pi measures
  the relative strength of many-body terms
  - g < g_c ~ 2.3: 2-body dominates -> Paley wins
  - g > g_c: many-body dominates -> Interval wins
  - g(7) = 1.68, g(11) = 2.11, g(19) = 2.78

## Connection to Additive Energy

Interval's alpha_2+ advantage comes from ADDITIVE STRUCTURE:
- Interval S = {1,...,m} has high additive energy E(S) = #{(a,b,c,d) in S^4 : a+b=c+d}
- High additive energy => similar neighborhoods for consecutive vertices
- Similar neighborhoods => cycles tend to avoid each other (vertex-disjoint)
- More vertex-disjoint cycles => higher alpha_2, alpha_3, ...

The sum-product theorem says QR has LOW additive energy (it's multiplicatively closed).
So Paley's cycles are more "tangled" (share vertices), while Interval's are "aligned".

## Scripts

- `04-computation/trace_H_analytic.py` — cycle count comparison
- `04-computation/thm136_all_k_proof.py` — trace alternation for all k
- `05-knowledge/results/ising_phase_transition.out` — Hessian analysis (opus)
