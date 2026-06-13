# THM-142: 3-Cycle Disjointness Excess Formula

**Status:** PROVED
**Author:** kind-pasteur-2026-03-12-S59
**Date:** 2026-03-12

## Statement

For primes p ≡ 3 (mod 4), the number of vertex-disjoint 3-3 cycle pairs satisfies:

    disjoint_3(Interval) - disjoint_3(Paley) = p(p-1)(p+1)(p-3) / 192

This quantity is always positive for p ≥ 5 and grows as Θ(p⁴).

## Proof

Both Interval and Paley are regular tournaments on Z_p with score m = (p-1)/2, so both have exactly c₃ = p(p²-1)/24 three-cycle vertex sets (standard formula for regular tournaments).

By circulant symmetry, every vertex v belongs to n_v = 3c₃/p = (p²-1)/8 three-cycle vertex sets.

**Key identity (inclusion-exclusion):**

    disjoint_pairs = C(c₃, 2) - Σ_v C(n_v, 2) + #{overlap ≥ 2 pairs}

The first two terms are IDENTICAL for Paley and Interval (same c₃, same n_v by regularity). Therefore:

    Δ_disjoint = #{ov≥2}_Interval - #{ov≥2}_Paley

**Computing #{ov≥2} (pairs sharing ≥ 2 vertices):**

For a circulant tournament, #{ov=2} = Σ_{d=1}^{m} p · C(co_occ_3(d), 2), where the sum is over distinct gap values.

For **Paley** (co_occ = c = (p+1)/4 constant):
    #{ov=2}_P = m · p · C(c, 2) = m · p · c(c-1)/2

For **Interval** (co_occ(d) = d by THM-141):
    #{ov=2}_I = p · Σ_{d=1}^{m} C(d, 2) = p · Σ_{d=1}^{m} d(d-1)/2 = p · m(m+1)(m-1)/6

**Computing the excess:**

    Δ = #{ov=2}_I - #{ov=2}_P
      = p · [(m²-1)/6 - (p+1)(p-3)/32] · m

Substituting m = (p-1)/2:
    m²-1 = (p²-2p-3)/4, so (m²-1)/6 = (p²-2p-3)/24
    (p+1)(p-3)/32 = (p²-2p-3)/32

    Bracket = (p²-2p-3) · [1/24 - 1/32] = (p²-2p-3) · 8/768 = (p²-2p-3)/96

    Δ = p · (p-1)/2 · (p²-2p-3)/96 = p(p-1)(p+1)(p-3)/192  ∎

## Verification

| p  | Formula  | Computed (Paley → Int) | ✓ |
|----|----------|------------------------|---|
| 7  | 7·6·8·4/192 = 7    | 7 → 14 (Δ=7)       | ✓ |
| 11 | 11·10·12·8/192 = 55 | 495 → 550 (Δ=55)   | ✓ |
| 19 | 19·18·20·16/192 = 570 | 23370 → 23940 (Δ=570) | ✓ |

## Exact formulas for absolute counts

**Paley disjoint 3-3 pairs:**
    D_P = C(c₃,2) - p(p²-1)(p²-9)/128 + m·p·(p+1)(p-3)/32

**Interval disjoint 3-3 pairs:**
    D_I = D_P + p(p-1)(p+1)(p-3)/192

## Interpretation: The Phase Transition Mechanism

The excess disjoint pairs arise from Interval's **bimodal overlap distribution**:
- Interval has MORE overlap=2 pairs (heavily-overlapping cycles using nearby vertices)
- Interval has MORE overlap=0 pairs (disjoint cycles in different neighborhoods)
- Interval has FEWER overlap=1 pairs (moderate overlap)

By inclusion-exclusion, more ov=2 pairs create more disjoint pairs (the double-counting correction adds back disjoint pairs).

This mechanism underpins the Paley → Interval H-maximizer phase transition at p ≈ 13:
- Small p: Paley's cycle count advantage (alpha_1) dominates H = Σ α_j · 2^j
- Large p: Interval's disjointness advantage (alpha_2, alpha_3, ...) dominates
- The crossover occurs when the O(p⁴) disjointness excess overcomes the O(p²) cycle count deficit

## Related

- THM-141 (Interval co-occurrence formula — key ingredient)
- THM-027 (BIBD cycle analysis)
- HYP-480 (Interval is global H-maximizer for p ≥ 13)
- THM-139 (Chirality dichotomy for Paley)
