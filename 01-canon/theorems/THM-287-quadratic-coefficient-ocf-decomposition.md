# THM-287: Quadratic Multilinear Coefficient via OCF Decomposition

**Status:** PROVED (algebraic for all n, verified n=4,5,6)
**Filed by:** opus-2026-04-04-S2
**Depends on:** OCF (H=I(Ω,2)), THM-286

## Statement

The quadratic multilinear coefficient c_{ij} in H(t₁,...,tₘ) = Σ_S c_S Π_{k∈S} t_k
decomposes as:

  **c_{ij} = 2·A(i,j) + 4·B(i,j)**

where:
- A(i,j) = Σ_C [coeff of t_i·t_j in χ_C(t)]  (sum over all directed odd cycles C)
- B(i,j) = Σ_{(C₁,C₂) vertex-disjoint} [coeff of t_i in χ_{C₁}] · [coeff of t_j in χ_{C₂}]

No higher independent set sizes (|I| ≥ 3) contribute to the quadratic coefficient.

## Proof

By OCF: H(t) = Σ_{I independent in Ω(t)} 2^|I| = Σ_{I vertex-disjoint odd cycles} 2^|I| · Π_{C∈I} χ_C(t)

where χ_C(t) is the indicator that cycle C exists in tournament T(t).

**Key lemma:** χ_C(0,...,0) = 0 for every directed odd cycle C.

*Proof of lemma:* T(0,...,0) is the transitive tournament (all tiles forward, all arcs go from higher to lower). The transitive tournament has no directed cycles. ∎

**Consequence:** The constant term of every χ_C is 0. Therefore every cycle indicator polynomial starts at degree ≥ 1.

**For independent sets of size |I| = k:** The product Π_{C∈I} χ_C has minimum degree ≥ k (each factor contributes degree ≥ 1).

For the quadratic coefficient c_{ij} (degree 2), we need contributions from terms of degree exactly 2. This requires k ≤ 2:
- k = 0: contributes only to degree 0 (constant = 1). No.
- k = 1: single cycle C can contribute if χ_C has a t_i·t_j term. Weight: 2¹ = 2.
- k = 2: pair (C₁, C₂) contributes if χ_{C₁} has t_i and χ_{C₂} has t_j (or vice versa). Weight: 2² = 4. (For χ_{C₁} to have t_i·t_j, we'd need C₂ to contribute constant 0. So only the split works.)
- k ≥ 3: each of ≥ 3 factors has degree ≥ 1, so product has degree ≥ 3 > 2. **No contribution.** ∎

## Corollaries

1. **Linear coefficients:** c_i = 2·Σ_C [coeff of t_i in χ_C(t)]. Only single cycles contribute.

2. **Shared-vertex pairs:** If tiles i and j share a vertex, then B(i,j) = 0 (no vertex-disjoint cycle pair can include both tiles). So c_{ij} = 2·A(i,j).

3. **n ≤ 5:** B(i,j) = 0 for ALL pairs (not enough vertices for disjoint cycle pairs). All quadratic coefficients come from single cycles only.

4. **n ≥ 6:** B(i,j) can be nonzero for vertex-disjoint tile pairs, adding 4× corrections that create non-power-of-2 values (e.g., c₂ = 10 = 2·3 + 4·1).

5. **Degree-k generalization:** The degree-k coefficient c_S (with |S|=k) gets contributions only from independent sets of size |I| ≤ k. Weight 2^|I|. This gives a finite hierarchical decomposition of ALL multilinear coefficients.

## Computational Verification

PERFECT MATCH at n=4, n=5, n=6 between OCF decomposition and direct Möbius inversion.

n=6 examples of the decomposition:
| Tile pair | Type | Skips | c₂ | TypeA (×2) | TypeB (×4) |
|-----------|------|-------|-----|-----------|-----------|
| (6,1)-(5,2) | disjoint | (5,3) | 16 | 8 | 8 |
| (5,1)-(6,3) | disjoint | (4,3) | 10 | 6 | 4 |
| (3,1)-(6,4) | disjoint | (2,2) | 4 | 0 | 4 |
| (6,1)-(3,1) | shared | (5,2) | -4 | -4 | 0 |
| (3,1)-(5,3) | cross-end | (2,2) | 4 | 4 | 0 |

## See Also
- THM-284 (linear coefficient = 2^(skip-1))
- THM-286 (all coefficients even)
- THM-260 (Walsh degree bound)
- Scripts: h_ocf_multilinear_v2.py, h_ocf_quadratic_anatomy.py
