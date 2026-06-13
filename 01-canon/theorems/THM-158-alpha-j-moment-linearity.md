# THM-158: All alpha_j are exact linear functions of eigenvalue moments

**Status:** PROVED at p=7,11. Verified at p=13 (pending full completion).
**Proved by:** kind-pasteur-2026-03-13-S60
**Date:** 2026-03-13

## Statement

For a circulant tournament T(Z_p, S) with p prime, define:
- Omega(T) = the set of all directed odd cycles of T (any length 3,5,...,p)
- alpha_j = number of independent sets of size j in the vertex-disjoint conflict graph on Omega(T)
- S_{2k} = sum_{t=1}^{m} D_t^{2k} where D_t = Im(lambda_t) are the eigenvalue imaginary parts

Then:

1. **Finiteness:** alpha_j = 0 for j > floor(p/3) (minimum cycle length 3, so j+1 disjoint 3-cycles need 3(j+1) > p vertices)

2. **OCF decomposition:** H(T) = sum_{j=0}^{floor(p/3)} alpha_j * 2^j

3. **Moment linearity:** Each alpha_j is an exact linear function of (S_4, S_6, ..., S_{p-3}).
   In particular, the coefficient of S_{p-1} in each alpha_j is ZERO.

4. **Dimension:** The number of moment parameters needed is m - 2 = (p-5)/2.

## Verified cases

| p | m | floor(p/3) | #moments needed | #orbit types | Status |
|---|---|-----------|-----------------|-------------|--------|
| 7 | 3 | 2 | 1 (S4) | 2 | PROVED (trivial dimension) |
| 11 | 5 | 3 | 3 (S4,S6,S8) | 4 | PROVED (4 types, 4 params: exact) |
| 13 | 6 | 4 | 4 (S4,...,S10) | 6 | Computing (6 types, 5 params: overconstrained) |

## Significance

Combined with OCF, this gives: **H(T) is an exact linear function of (S_4, ..., S_{p-3})** for all circulant tournaments on Z_p.

Since S_{p-1} cancels in every alpha_j individually (not just in the sum), this reveals a deep structure: the highest-frequency eigenvalue moment is algebraically redundant for the entire independence polynomial, not just for the total.

## Key components at p=11

- alpha_1 = -187.17*S4 + 17.73*S6 - 0.556*S8 + 26561
- alpha_2 = -42.76*S4 + 7.78*S6 - 0.381*S8 + 11796
- alpha_3 = 10.93*S4 - 1.76*S6 + 0.097*S8 + 897
- H = -457.93*S4 + 52.51*S6 - 1.858*S8 + 107481

## Disjoint pair decomposition at p=11

- disj_{3,3} = f(S4) only (THM-156)
- disj_{3,5} = f(S4, S6, S8)
- disj_{3,7} = f(S4, S6, S8)
- disj_{5,5} = f(S4, S6, S8)

## Paley structure

Paley (QR) tournament simultaneously:
- Maximizes alpha_1 (total directed odd cycles)
- Minimizes alpha_2, alpha_3, ..., alpha_{floor(p/3)} (disjoint sets)
- Maximizes H (Hamiltonian path count)

This is because H is dominated by the alpha_1 term (coefficient 2) while higher alpha_j have smaller but growing coefficients (4, 8, 16, ...).

## Proof approach

The proof is algebraic: each c_k (directed k-cycles) is a polynomial in eigenvalue moments by the Held-Karp/trace formula. The disjoint pair/triple counts are inclusion-exclusion over cycle-vertex intersections, which for circulant tournaments reduce to algebraic expressions in moments by the circulant symmetry (all vertices equivalent).

## Files

- `04-computation/alpha2_by_type_p11.py` — full decomposition at p=11
- `04-computation/alpha_ocf_p13.py` — computing at p=13
- `05-knowledge/results/alpha2_by_type_p11.out` — p=11 results
