# THM-161: Effective Marginal H Coefficients for Circulant Tournaments

**Status:** VERIFIED (p=7, p=11)
**Session:** kind-pasteur-2026-03-13-S60

## Context

For circulant tournaments on Z_p with orientation set S:

    H(T_S) = 1 + 2*N + 4*alpha_2 + 8*alpha_3    (OCF decomposition)

where N = sum c_k, alpha_j counts j-tuples of pairwise-disjoint directed cycles.

## Statement

Define the **effective marginal H coefficient** b_k as the NET effect on H per unit change in c_k, across the family of circulant orientations:

    H = const + sum_k b_k * c_k

At p=11 (with b_11 = 2 fixed by construction):
- b_5 = -36.9 (NEGATIVE — 5-cycles REDUCE H on net)
- b_7 = -1.5 (NEGATIVE — 7-cycles slightly reduce H)
- b_9 = +3.6 (POSITIVE — 9-cycles are super-valuable, more than Ham cycles)
- b_11 = +2.0 (NEUTRAL — Hamiltonian cycles are "free" benefit)

## Mechanism

Each k-cycle C on vertex set V contributes:
- +1 to N (adding 2 to H directly)
- +C_odd(T[comp(V)]) to alpha_2 (adding 4*C_odd(comp) to H)
- Higher-order alpha_3 terms

The "per-cycle disjoint contribution" at Paley p=11:
- k=3: C_odd(comp) = 172 (constant!) — each 3-cycle is disjoint from 172 others
- k=5: C_odd(comp) = 7-16 (avg 12.56) — moderate disjoint contribution
- k=7: C_odd(comp) = 0-2 (avg 1.22) — minimal disjoint contribution
- k=9: C_odd(comp) = 0 — complement too small for any cycle
- k=11: C_odd(comp) = 0 — complement is empty

## Why Negative b_5?

The effective b_k is NOT the per-cycle marginal H. It captures the cross-orientation effect: when one orientation has more 5-cycles than another, it typically has MORE disjoint pairs (higher alpha_2), but the alpha_2 increase is correlated with LOSS of 9-cycles and 11-cycles. The effective coefficient captures this total effect including cross-length interactions.

## Why Paley Maximizes H

Paley maximizes EVERY c_k simultaneously (Savchenko's theorem for DRTs). The large-cycle advantage (c_9, c_11) is the dominant mechanism:

- Paley's c_9 advantage = +858 over Class B, contributing 3.6 * 858 = +3089 to dH
- Paley's c_11 advantage = +352, contributing 2.0 * 352 = +704
- Paley's c_5 "excess" = +44, "costing" 37 * 44 = -1628 (but offset by c_9, c_11)
- Paley's c_7 "excess" = +374, "costing" 1.5 * 374 = -561

Net: +3089 + 704 - 1628 - 561 = +1604, close to actual dH = 1628.

## Key Structural Facts

1. **c3 is constant for ALL regular tournaments**: c3 = n(n-1)(n+1)/24 at odd n (depends only on degree sequence, not orientation)
2. **4 distinct H classes** at p=11 with sizes 2:10:10:10
3. **alpha_2 is a function of N alone** within the circulant family at both p=7 and p=11
4. **Overlap coefficient** |d(alpha_2)/dN| < 0.5 at both p=7 (0.33) and p=11 (0.07), meaning more cycles always increases H
5. **H-class ordering follows c_11** (Hamiltonian cycles), not c_5/c_7/c_9

## Verification

- alpha_decomposition_all_orientations.out: all 32 orientations, exact H decomposition
- N_maximization_paley.out: Paley maximizes every c_k at p=7 and p=11
- ham_cycle_dominance.out: effective coefficient computation with exact verification
- H_maximization_mechanism.out: overlap coefficient analysis
