---
id: THM-879
title: THE EXACT COVERING LOCUS OF THE TIGHT AP — the 46/58-chamber case analysis run: between Farey-13 breakpoints every gap of {i·x} is linear in x, so {x : maxgap ≤ 1/7} is an exact union of TWELVE closed rational intervals of total measure 601/1078 ≈ 0.5575 (witness/lonely measure 477/1078 ≈ 0.4425); the locus starts EXACTLY at the clock 1/14 (the whole far tail x < 1/14 is lonely); every one of the 29 order-level chamber words carries covering mass
status: PROVED-EXACT (linear programming per cell, exact rationals; refereed against a 3000-point maxgap scan, 0 mismatches)
source: klein-2026-07-16-S314; runs the chamber case analysis promised in cont.5 (T1536)
depends_on: [LEM-020 (frame; mid-band law), the Farey-14/13 chamber structure (cont.5)]
verification: 04-computation/ap13_chamber_covering_locus_klein_S314.py -> 05-knowledge/results/ap13_chamber_covering_locus_klein_S314.out (4/4)
---

# THM-879 — the exact covering locus of {1,…,13}

Between consecutive Farey-13 fractions the circular order of {frac(i·x)} is constant, each
gap is a linear function of x, and covering (maxgap ≤ 1/7 ⟺ μ₀ = 0) is a one-variable
linear system per cell. Solving all 58 cells exactly:

> **Cov = [1/14, 8/49] ∪ [6/35, 4/21] ∪ [3/14, 5/21] ∪ [13/49, 11/35] ∪ [5/14, 11/28] ∪
> [20/49, 36/77] ∪ (mirror images under x ↦ 1−x),  measure 601/1078.**

Consequences and structure:
- **The far tail is all-lonely:** x ∈ (0, 1/14) is entirely witness — the locus begins at
  the clock 1/14 exactly (the mid-band law's boundary made sharp).
- **Witness measure = 477/1078 ≈ 0.4425**, symmetric under negation; widths include 2/105
  (numerically = THM-863's c₃ — flagged coincidence, unproved), 32/539, 9/98, 1/28, 1/42,
  12/245.
- **Chamber census:** 58 order-cells carrying 29 words (each on exactly 2 mirror cells);
  ALL 29 words contain covering mass — no chamber type is fully lonely; the finer 46-word
  count (cont.5) is the gap-VALUE-refined stratification (cells subdivide where two of the
  three gap values coincide).
- The tight AP is "lonely at 44.2% of slow times and covering at 55.8%" — exact; every
  future argument about the extremal instance can now cite intervals instead of sampling.
