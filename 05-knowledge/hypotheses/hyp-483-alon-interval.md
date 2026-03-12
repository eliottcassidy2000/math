---
hypothesis_id: HYP-483
title: "Alon's interval tournament becomes H-maximizer at large n"
status: PARTIALLY_CONFIRMED
tested_by: opus-2026-03-12-S60
date: 2026-03-12
related: [THM-135, HYP-464, HYP-466]
---

## Statement

The interval tournament C_n (connection set S={1,...,⌊(n-1)/2⌋}) becomes the global
H-maximizer among all tournaments on n vertices for sufficiently large n.

## Evidence

### Small n (interval does NOT maximize):
| n | H(interval) | P(n) = max | Winner |
|---|-------------|------------|--------|
| 3 | 3 | 3 | TIE |
| 5 | 15 | 15 | TIE |
| 7 | 175 | 189 | Paley |
| 9 | 3267 | 3357 | Unknown |
| 11 | 93027 | 95095 | Paley |

### Large n (interval wins among circulants):
| p | H(Paley) | H(Interval) | Winner (circulants) |
|---|----------|-------------|---------------------|
| 19 | 1,172,695,746,915 | 1,184,212,824,763 | Interval |

### Alon (1990) proved:
- P(T_n) ≥ n!/(2+o(1))^n (interval is asymptotically optimal)
- lim P(T_n)^{1/n} = lim P(n)^{1/n} = n/(2e)
- Conjectured C(T_n) = C(n) (maximizes Hamiltonian cycles too)

### OCF mechanism for crossover:
At p=7: Paley has more cycles (α_1=80 vs 59) but fewer disjoint pairs (α_2=7 vs 14).
The α_1 term dominates at small n. At large n, higher α_k (disjoint cycle packings)
dominate due to 2^k weighting, and interval's local structure creates better packings.

## Status

The interval tournament is CONFIRMED as asymptotically optimal (Alon 1990).
The crossover to finite-n dominance is CONFIRMED at p=19 among circulants.
Global maximality at p=19 (among ALL tournaments) is UNVERIFIED.

## Key reference

Noga Alon, "The Maximum Number of Hamiltonian Paths in Tournaments",
Combinatorica 10(4), 319-324 (1990). Uses Brégman's theorem (Minc's conjecture).
