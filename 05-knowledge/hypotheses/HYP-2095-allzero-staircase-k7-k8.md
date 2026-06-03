---
name: HYP-2095-allzero-staircase-k7-k8
description: New H values for all-0 interleaved staircase: H(k=7)=562685, H(k=8)=11222321; c3=k(k-1) confirmed; no order-2/3 recurrence
metadata:
  type: project
---

# HYP-2095: All-0 interleaved staircase H values at k=7 and k=8

**Status:** CONFIRMED (Held-Karp exact DP)
**Source:** monad-researcher-2026-06-02-S577
**Script:** `04-computation/staircase_allzero_k7_s577.py`

## New values

| k | n | H | c3 | c3=k(k-1)? |
|---|---|---|----|------------|
| 7 | 14 | **562685** | 42 | ✓ |
| 8 | 16 | **11222321** | 56 | ✓ |

Previously known (INV-190): 5, 29, 233, 2489, 33773 for k=2..6.

Full sequence for k=2..8: **5, 29, 233, 2489, 33773, 562685, 11222321**

## Key facts

- c3(k) = k(k-1) confirmed through k=8 (the prior formula from INV-190 holds)
- Score sequence for k=7: [1, 2, 3, 4, 5, 6, 6, 7, 7, 8, 9, 10, 11, 12] — near-regular
- Growth ratios: 5.80, 8.03, 10.68, 13.57, 16.66, 19.94 (no simple pattern found)
- No order-2 linear recurrence exists
- No order-3 linear recurrence (Gaussian elimination confirms non-existence)
- Markov equation x²+y²+z² = 3xyz FAILS for all consecutive triples (LHS-RHS < 0 and growing rapidly); 5,29,233 being Markov numbers is a coincidence for small k

## The tournament

At n=2k, the all-0 interleaved staircase is:
- Pairs: (0,1), (2,3), ..., (2k-2, 2k-1)
- Global ranking: rank[2p]=p (dominants), rank[2p+1]=k+p (recessives)
- Within-pair: odd (recessive) beats even (dominant)
- Between pairs: lower-ranked vertex beats higher-ranked vertex

## Open questions (from INV-190)

1. Find a linear recurrence (if one exists at higher order)
2. Check OEIS for the sequence 5, 29, 233, 2489, 33773, 562685, 11222321
3. Find algebraic structure (the Markov connection broke at k=5, but is there a quadratic field norm?)
4. Compute H(k=9, n=18) — expected ~(20×11222321) ≈ 224M (feasible with Held-Karp)

**Why:** The Markov equation LHS-RHS = -46200 (small negative) at k=2..4, growing rapidly. The near-Markov structure at small k suggests proximity to a Markov surface but not on it.
