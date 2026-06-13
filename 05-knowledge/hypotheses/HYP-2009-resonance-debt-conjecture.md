---
id: HYP-2009
name: resonance-debt-conjecture
status: OPEN
date: 2026-06-01
session: opus-2026-06-01-S531
---

# HYP-2009: Resonance Debt Conjecture — initial segment maximizes debt/credit

## Statement

The lonely measure decomposes as:

    μ(LONELY) = CREDIT + DEBT

where CREDIT = ((n-2)/n)^{n-1} ≈ e^{-2} and DEBT = Σ_{r≥2} resonance terms.

**Conjecture:** Among all primitive speed sets at a given n, the initial
segment {1,...,n-1} has the maximum |DEBT/CREDIT| ratio, and this maximum
is exactly 1. Equivalently:

    |DEBT(V)| ≤ CREDIT    for all primitive V

with equality iff V = {1,...,n-1} (up to scaling).

## Consequence

If HYP-2009 is true: μ(LONELY) = CREDIT - |DEBT| ≥ 0 for all primitive V.
- For V = {1,...,n-1}: μ = 0 exactly (wall-only, proved by regular polygon)
- For V ≠ {1,...,n-1}: μ > 0 (open lonely times, positive measure)

This would prove LRC for all n.

## Evidence

Verified at n=3,...,8,14 (initial segment: ratio 1.000000 exactly).
Verified for non-initial sets: primes at n=14 has ratio 0.338.
All 3631+ tested speed sets have ratio ≤ 1.

## Structure of the debt

The debt decomposes by resonance order:
- Order 2 (pairwise): closed form via Bernoulli polynomials B₂
- Order ≥ 3 (higher): the "inside debt" — controlled by the polygon's diagonals

The order-2 (pairwise) fraction of total debt varies:
  n=3: 99.8% (pairwise dominates)
  n=5: 7.3% (higher dominates!)
  n=7: 27% (mixed)
  n=14: 83% (pairwise large but higher matters)

## Why the initial segment is extremal

The initial segment speeds {1,...,n-1} form an ARITHMETIC PROGRESSION.
At the lonely time t=k/n, the positions are {k/n, 2k/n, ..., (n-1)k/n}
= the regular n-gon. This is the UNIQUE configuration where:
1. All gaps exactly = 1/n (boundary lonely)
2. The tournament is the regular circular tournament
3. The resonance terms sum to exactly -CREDIT (perfect cancellation)
4. The Gauss sum √n controls the character-sum balance (oracle-S529)

Non-AP speed sets break this perfect balance, reducing the debt.

## Proof strategy

1. Show the pairwise debt is maximized by AP (Bernoulli + arithmetic)
2. Show the higher-order debt is maximized by AP (character sums)
3. Use the Gauss/Kloosterman bound √n per shell to control each order
4. Sum the bounds and show total debt ≤ CREDIT
