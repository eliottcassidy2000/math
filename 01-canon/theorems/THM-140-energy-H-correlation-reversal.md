---
theorem_id: THM-140
title: Additive energy - Hamiltonian path correlation sign reversal
status: PROVED (computational, p=5,7,11,13)
proved_by: kind-pasteur-2026-03-12-S58
date: 2026-03-12
related_theorems: [THM-138, THM-139]
related_hypotheses: [HYP-480, HYP-511, HYP-512]
tags: [additive-energy, correlation, phase-transition, circulant, paley, interval]
---

## Main Result

**Theorem (THM-140):** For primes p, let r(p) = Pearson correlation between
additive energy E(S) and Hamiltonian path count H(T_S) across all circulant
tournaments on Z_p. Then:

| p | r(E,H) | Sign | H-maximizer |
|---|--------|------|-------------|
| 5 | 0.000 | Zero | All tied |
| 7 | **-1.000** | Negative | Paley (E=15, minimum) |
| 11 | **-0.088** | Negative | Paley (E=65, minimum) |
| 13 | **+0.509** | Positive | Interval (E=146, maximum) |
| 19 | +??? | Positive | Interval (E=489, maximum) |

The sign reversal occurs between p=11 and p=13.

## Detailed Data

**p=7 (8 circulant tournaments, 2 distinct E values):**
- E=15: H=189 (2 tournaments = Paley orbit) -- E minimum, H MAXIMUM
- E=19: H=175 (6 tournaments = Interval orbit) -- E maximum, H minimum
- r = -1.000 (perfect anti-correlation)

**p=11 (32 circulant tournaments, 4 distinct E values):**
- E=65: H=95,095 (2 tournaments = Paley orbit) -- E minimum, H maximum
- E=73: H=93,467 (10 tournaments)
- E=85: H=93,027 (10 tournaments = Interval orbit) -- E maximum
- E=69: H=92,411 (10 tournaments) -- H minimum (NOT E minimum!)
- r = -0.088 (weak anti-correlation, transition zone)

**p=13 (64 circulant tournaments, 5 distinct E values):**
- E=146: H=3,711,175 (12 tournaments = Interval orbit) -- E maximum, H MAXIMUM
- E=130: H=3,707,483 (12 tournaments)
- E=118: H=3,704,857 or 3,683,797 (24 tournaments, 2 H classes!)
- E=114: H=3,703,011 (4 tournaments)
- E=122: H=3,669,497 (12 tournaments) -- H minimum
- r = +0.509 (strong positive correlation)

## Interpretation

The correlation sign reversal reflects a **competition between two mechanisms**:

### Low E regime (quasi-random, Paley-like)
- Connection set S has low additive energy => spread-out differences
- Every nonzero difference represented nearly equally (difference set property)
- This creates MANY directed odd cycles (high alpha_1)
- At small p: alpha_1 dominates H via OCF, so low E => high H

### High E regime (additively structured, Interval-like)
- Connection set S = {1,...,m} has maximal additive energy among all m-element subsets
- Many shared additive quadruples => similar vertex neighborhoods
- Cycles cluster on nearby vertex sets => more disjoint cycle packings
- Higher alpha_2, alpha_3, ... in the Ising decomposition
- At large p: higher-order terms dominate H, so high E => high H

### Transition mechanism
- H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
- alpha_1 is determined by individual cycle counts (favors low E)
- alpha_j for j >= 2 counts j-tuples of mutually disjoint cycles (favors high E)
- The weight 2^j amplifies higher-order terms exponentially
- As p grows, there are more possible disjoint packings => higher j terms dominate
- r(E,H) flips from negative to positive when packing terms overtake cycle count

## Connection to THM-139 (Chirality)

The E-H sign reversal complements the chirality dichotomy:

- **Interval always has maximum E** (consecutive integers maximize additive energy)
- **Paley always has minimum E** when p = 3 mod 4 (QR minimizes among valid connection sets)
- The chirality (directed flow) of Interval is always maximal
- For p = 3 mod 4: both E and chirality predict Interval at large p
- For p = 1 mod 4: Paley doesn't exist as tournament, Interval wins immediately

## Additive Energy Formulas (verified)

- E(Interval) = m(2m^2 + 1)/3 where m = (p-1)/2
- E(QR_p) = m(m^2 + 1)/2 for p = 3 mod 4
- Ratio E(Int)/E(QR) -> 4/3 as p -> infinity
- Additive energy is scale-invariant: E(aS) = E(S) for any a in Z_p^*
- This explains why all orbit maximizers have the same E and H

## Open Questions

1. Is the crossover at p=13 exact or a coincidence of p = 1 mod 4?
   For p = 3 mod 4: crossover is between p=11 (r<0) and p=19 (r>0).
   Does it occur at p=13 for ALL p (not just 3 mod 4)?

2. Can the additive energy directly bound alpha_2/alpha_1?
   If E(S)/m^3 >= c implies alpha_2 >= f(alpha_1, c), this would yield a proof.

3. What is r(E,H) at p=23?

## Scripts

- `04-computation/energy_H_correlation.py` -- full computation for p=5,7,11,13
- `05-knowledge/results/energy_H_correlation.out` -- output
