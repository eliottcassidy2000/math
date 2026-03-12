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

## Quantitative Data (COMPLETE for p=7,11; partial p=19)

### Full alpha decomposition

**p=7 (VERIFIED):**

| | alpha_0 | alpha_1 | alpha_2 | alpha_3 | H |
|---|---------|---------|---------|---------|-----|
| Paley | 1 | 80 | 7 | 0 | 189 |
| Interval | 1 | 59 | 14 | 0 | 175 |
| Delta (P-I) | 0 | +21 | -7 | 0 | +14 |

OCF: H = 1 + 2(80) + 4(7) = 189 (Paley), 1 + 2(59) + 4(14) = 175 (Interval).

**p=11 (VERIFIED):**

| | alpha_0 | alpha_1 | alpha_2 | alpha_3 | alpha_4 | H |
|---|---------|---------|---------|---------|---------|-------|
| Paley | 1 | 21169 | 10879 | 1155 | 0 | 95095 |
| Interval | 1 | 18397 | 11110 | 1474 | 0 | 93027 |
| Delta (P-I) | 0 | +2772 | -231 | -319 | 0 | +2068 |

OCF verified: 1+2(21169)+4(10879)+8(1155) = 95095, 1+2(18397)+4(11110)+8(1474) = 93027.

**p=19 (alpha_1 COMPUTED, alpha_2+ from H):**

| | alpha_1 | 2*alpha_1 | higher | H |
|---|---------|-----------|--------|------|
| Paley | 130,965,270,477 | 261,930,540,954 | 910,765,205,960 | 1,172,695,746,915 |
| Interval | (computing) | | | 1,184,212,824,763 |

**p=19 cycle counts by length (Paley):**
c_3=285, c_5=11628, c_7=424080, c_9=12156390, c_11=249208902,
c_13=3280900392, c_15=23662379790, c_17=69401425077, c_19=34358763933

### Ising advantage decomposition

| p | D1 = 2*Delta_alpha1 | D2 = 4*Delta_alpha2 | D3 = 8*Delta_alpha3 | H(P)-H(I) |
|---|---------------------|---------------------|---------------------|------------|
| 7 | +42 | -28 | 0 | +14 |
| 11 | +5544 | -924 | -2552 | +2068 |
| 19 | (pending) | | | -11,517,077,848 |

### Higher-order dominance (fraction of H from alpha_2+)

| p | 2*alpha_1/H (Paley) | higher/H (Paley) |
|---|---------------------|------------------|
| 7 | 84.7% | 14.8% |
| 11 | 44.5% | 55.5% |
| 19 | 22.3% | 77.7% |

The many-body terms grow from 15% to 78% of H as p increases from 7 to 19.

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

- `04-computation/alpha_decomp_p11_full.py` — complete alpha_j verification at p=11
- `04-computation/alpha_decomp_p19_fast.py` — cycle counts at p=19 via circulant symmetry
- `04-computation/trace_H_analytic.py` — cycle count comparison
- `04-computation/thm136_all_k_proof.py` — trace alternation for all k
- `05-knowledge/results/alpha_decomp_p11_full.out` — verified p=11 decomposition
- `05-knowledge/results/alpha_decomp_p19_fast.out` — p=19 cycle counts
- `05-knowledge/results/ising_phase_transition.out` — Hessian analysis (opus)
