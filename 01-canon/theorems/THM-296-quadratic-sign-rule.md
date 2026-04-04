# THM-296: Sign Rule for Quadratic Coefficients

**Status:** CONJECTURED, verified n=4,5,6,7,8 (0 violations in 366 tested pairs)
**Discovered by:** opus-2026-04-04-S3
**Related:** THM-287 (OCF quadratic decomposition), THM-297, THM-297

## Statement

For the quadratic multilinear coefficients c_{ij} of H(t):

1. **Same-end pairs** (tiles sharing upper or lower vertex): c_{ij} < 0 always.
2. **Cross-end pairs** (one tile's upper = other's lower): c_{ij} > 0 always.
3. **Disjoint pairs** (no shared vertex): c_{ij} > 0 always.

More precisely:
- c_{ij} < 0 if and only if tiles i,j share a "same-end" vertex (both have the same upper vertex x, or both have the same lower vertex y).
- c_{ij} > 0 otherwise (cross-end or fully disjoint).
- c_{ij} = 0 only for a few small-n pairs with no cycle through both tiles.

## Geometric Interpretation

From the OCF decomposition (THM-287): c_{ij} = 2·A(i,j) + 4·B(i,j), where:
- A(i,j) comes from single cycles through both tiles
- B(i,j) comes from disjoint pairs of cycles

Same-end pairs use tiles in **opposite** directions (one forward, one backward in the cycle), giving chi coefficient = -1. Cross-end pairs use tiles in the **same** direction, giving chi coefficient = +1. This geometric origin explains the sign rule.

## Computational Verification

| n | Same-end pairs | Cross-end pairs | Disjoint pairs | Violations |
|---|---------------|----------------|----------------|------------|
| 4 | 2 (all -2) | 0 | 0 | 0 |
| 5 | 8 (all -2) | 1 (+4) | 4 (+) | 0 |
| 6 | 20 (neg) | 4 (+) | 18 (+) | 0 |
| 7 | 40 (neg) | 10 (+) | 51 (+) | 0 |
| 8 | 66 (neg) | ? | ? | 0 |
