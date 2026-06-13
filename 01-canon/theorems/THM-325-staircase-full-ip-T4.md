---
id: THM-325
title: Full Independence Polynomial of T_4
status: VERIFIED
session: opus-2026-05-23-S5
tags: [staircase, independence-polynomial, real-rooted, TRRT]
related: [THM-318, THM-319, THM-322, THM-323, HYP-1732]
---

## Statement

The conflict graph Ω(T_4) of the all-0 staircase T_4 (n=8) has independence polynomial:

$$I(\Omega(T_4), x) = 1 + 68x + 24x^2$$

## Structure

T_4 has n=8 vertices. Odd cycle breakdown:
- 12 triangles (3-cycles)
- 28 pentagons (5-cycles)  
- 28 heptagons (7-cycles)
- Total: m = 68 odd cycles

The independence number α(Ω(T_4)) = 2, with α_2 = 24 disjoint pairs.

Pair breakdown by cycle-length combination:
- (3,3) pairs: 16 (two disjoint triangles)
- (3,5) pairs: 8 (one triangle + one pentagon)
- (3,7) pairs: 0 (no disjoint triangle-heptagon pair)
- (5,5) pairs: 0 
- (5,7) pairs: 0
- (7,7) pairs: 0

## Key Properties

- **Degree:** d(Ω(T_4)) = 2 (consistent with THM-318: floor(2·4/3) = 2)
- **α(Ω(T_4)) = 2:** maximum IS has size 2
- **Real-rooted:** discriminant = 68² - 4·24 = 4624 - 96 = 4528 > 0; roots ≈ -0.3548 and -67.645
- **Roots are both negative** ✓

## HYP-1732 Check

For T_4, α(Ω) = 2, so HYP-1732 applies. Among 68 cycles with α_2=24:

Full pair-partner testing (552 tests) showed 0 violations. All valid C* selections satisfy:
$$\alpha_2(\Omega) \leq p(m - p)$$

where p = #{B-cycles: disjoint from C*}. See 05-knowledge/results/hyp1732_full_investigation.out.

## Comparison with 3-Cycle Truncation

I_3(T_4, x) = 1 + 12x + 16x² (3-cycles only, from THM-319)

The full IP I(Ω(T_4), x) = 1 + 68x + 24x² includes contributions from 5- and 7-cycles:
- m grows from 12 to 68 (factor ~5.7)
- α_2 shrinks from 16 to 24... wait: I_3 has α_{3,2}=16, full IP has α_2=24. The 5- and 7-cycles ADD disjoint pairs (net +8).

Note: Even though 5- and 7-cycles add new IS members, their conflict structure is complex — the net effect on α_2 is modest (+8 pairs).

## Files

- Output: 05-knowledge/results/hyp1732_full_investigation.out
- Script: 04-computation/hyp1732_full_investigation.py and cycle_count_formula.py
