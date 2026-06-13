---
theorem_id: THM-132
title: Z_11 OCF structure — alpha_1 ordering does NOT match H ordering
status: PROVED (computational, exhaustive over all 32 circulants)
proved_by: opus-2026-03-12-S58
date: 2026-03-12
related_theorems: [THM-126, THM-129, THM-130]
related_hypotheses: [HYP-466]
tags: [paley, ocf, circulant, Z11, independence-polynomial]
---

## Main Result

Among all 32 circulant tournaments on Z_11, Paley T_11 (S = QR_11 = {1,3,4,5,9})
uniquely maximizes H = 95095. There are 4 distinct H values and the full OCF
decomposition H = 1 + 2α₁ + 4α₂ + 8α₃ is:

| S representative | H | α₁ | α₂ | α₃ | Paley? |
|-----------------|--------|--------|--------|------|--------|
| {1,3,4,5,9} | 95095 | 21169 | 10879 | 1155 | YES |
| {1,2,4,5,8} | 93467 | 19541 | 11220 | 1188 | |
| {1,2,3,4,5} | 93027 | 18397 | 11110 | 1474 | |
| {2,3,4,5,10} | 92411 | 19629 | 10912 | 1188 | |

Note: α₄ = 0 for all Z_11 circulants (4 disjoint odd cycles require ≥ 12 vertices).

## Key Finding: α₁ Ordering ≠ H Ordering

The α₁ (total odd cycle count) ordering is:
```
21169 (Paley) > 19629 > 19541 > 18397
```

But the H ordering is:
```
95095 (Paley) > 93467 > 93027 > 92411
```

The tournament with α₁ = 19629 (second-highest total cycles) has the LOWEST H!
This is because its α₂ = 10912 is low — having many cycles but fewer disjoint
pairs. Conversely, the cyclic interval (α₁ = 18397, lowest) compensates with
α₃ = 1474 (highest), and its 8α₃ = 11792 contribution rescues its H.

## Contribution Breakdown

| Tournament | 2α₁ (%) | 4α₂ (%) | 8α₃ (%) |
|-----------|---------|---------|---------|
| Paley | 42338 (44.5%) | 43516 (45.8%) | 9240 (9.7%) |
| H=93467 | 39082 (41.8%) | 44880 (48.0%) | 9504 (10.2%) |
| H=93027 | 36794 (39.6%) | 44440 (47.8%) | 11792 (12.7%) |
| H=92411 | 39258 (42.5%) | 43648 (47.2%) | 9504 (10.3%) |

Paley wins the 2α₁ contribution overwhelmingly (+3256 over second place), and
this single-term advantage is enough to overcome any α₂ or α₃ deficit.

## Comparison with p=7

At p=7 (THM-126 + dihedral_hp_analysis):
- H = 1 + 2α₁ + 4α₂ (α₃ = 0 forced by vertex count)
- Paley: α₁ = 80, α₂ = 7 → H = 1 + 160 + 28 = 189
- Non-Paley: α₁ = 59, α₂ = 14 → H = 1 + 118 + 56 = 175

At p=7, Paley wins via 2α₁ dominance (+42) despite losing on 4α₂ (-28).

**Pattern (HYP-466):** For Paley primes p ≡ 3 mod 4, the Paley tournament wins
among circulants primarily via the 2α₁ (total odd cycle count) term of OCF.
The advantage comes from having MORE odd cycles (especially c_5 — see THM-130),
not more disjoint pairs.

## Cycle Count Comparison

| Cycle length | Paley | H=93467 | H=93027 | H=92411 |
|-------------|-------|---------|---------|---------|
| 3-cycles | 55 | 55 | 55 | 55 |
| 5-cycles | 594 | 550 | 484 | 572 |
| 7-cycles | 3960 | 3586 | 3399 | 3729 |
| 9-cycles | 11055 | 10197 | 9350 | 10274 |
| 11-cycles | 5505 | 5153 | 5109 | 4999 |

3-cycles universal (= p(p²-1)/24 = 55 for all regular tournaments on 11 vertices).
Paley maximizes c_5, c_7, c_9, c_11 individually. The phase alignment principle
(THM-130) explains the c_5 advantage via Gauss sum structure.

## Scripts

Script: `04-computation/z11_ocf_structure.py`
Output: `05-knowledge/results/z11_ocf_structure.out`
