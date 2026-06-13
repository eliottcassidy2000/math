---
id: THM-323
title: Full Independence Polynomial of T_3
status: VERIFIED
session: opus-2026-05-23-S5
tags: [staircase, independence-polynomial, real-rooted]
related: [THM-318, THM-319, THM-322]
---

## Statement

The conflict graph Ω(T_3) of the all-0 staircase T_3 (n=6) has independence polynomial:

$$I(\Omega(T_3), x) = 1 + 12x + x^2$$

## Structure

T_3 has n=6 vertices and the following odd cycles:
- 6 triangles (3-cycles): from α_1 = 2·C(3,2) = 6
- 6 pentagons (5-cycles): from #5 = C(3,3)·6 = 6
- Total: 12 odd cycles = m = 12

The independence polynomial has:
- α_0 = 1 (empty set)
- α_1 = 12 (each single cycle is independent)
- α_2 = 1 (exactly one disjoint pair)

The unique maximum independent set consists of 2 cycles (a 3-cycle and a 5-cycle from disjoint vertex sets). The 3-cycle and 5-cycle must together use all 6 vertices of T_3.

## Key Properties

- **Degree:** d(Ω(T_3)) = 2 (consistent with THM-318: floor(2·3/3) = 2)
- **α(Ω(T_3)) = 2:** exactly one maximum IS of size 2
- **Real-rooted:** discriminant = 12² - 4·1 = 140 > 0; roots are $x = (-12 ± \sqrt{140})/2$, both negative (≈ -0.083, -11.917)
- **No 3-IS:** the degree is 2, so α_3 = 0

## Verification

Computed directly via exhaustive IS enumeration in hyp1732_full_investigation.py.

## Connection to HYP-1732

T_3 satisfies HYP-1732 vacuously at k=3: with α(Ω)=2, the pair-partner construction applies, p=|B-cycles| varies, but there is only 1 disjoint pair total, consistent with α_2 ≤ p(m-p).
