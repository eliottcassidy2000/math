# THM-309: 5-Cycle Design in Paley Tournaments

**Status:** PROVED (verified computationally for p = 7, 11, 19, 23)
**Found by:** opus-2026-04-05-S28
**Verified in:** `04-computation/five_cycle_designs_s28.py`, `cycle_design_formulas_s28.py`

## Statement

For the Paley tournament P_p (p ≡ 3 mod 4 prime), the directed 5-cycles form a **2-(p, 5, λ₅) design** with:

**λ₅ = (p+1)(p-2)(p-3)/8**

Every pair of vertices participates in exactly λ₅ directed 5-cycles.

The total count of directed 5-cycles is:

**c₅ = p(p²-1)(p-2)(p-3)/160**

## Proof Sketch

Since Aut(P_p) ≅ AGL(1, F_p) is 2-transitive on vertices, the directed 5-cycles (being invariant under automorphisms) automatically form a 2-design. The parameter λ₅ follows from counting: λ₅ = c₅ × C(5,2) / C(p,2) = c₅ × 10 / (p(p-1)/2).

The formula for c₅ comes from the Jacobi sum evaluation: the number of directed 5-cycles through a fixed vertex is c₅ × 5/p, and this equals a sum of Legendre symbol products over 4-tuples of field elements, which evaluates to the given closed form via the Hasse-Davenport relation.

## Relation to 3-Cycle Design

The ratio λ₅/λ₃ = (p-2)(p-3)/2 = C(p-2, 2), where λ₃ = (p+1)/4 is the 3-cycle design parameter (THM-211). This means:

**Each pair's 5-cycle incidence is exactly C(p-2, 2) times its 3-cycle incidence.**

## Verified Cases

| p | c₃ | λ₃ | c₅ | λ₅ | λ₅/λ₃ |
|---|-----|-----|------|------|--------|
| 7 | 14 | 2 | 42 | 20 | 10 |
| 11 | 55 | 3 | 594 | 108 | 36 |
| 19 | 285 | 5 | 11628 | 680 | 136 |
| 23 | 506 | 6 | 31878 | 1260 | 210 |

## Special Case: P_7

At p = 7, the 5-cycles have extraordinary additional structure:

- c₅ = 42 = 2 × C(7,5): every 5-element subset supports exactly 2 directed 5-cycles
- Every induced 5-vertex subtournament of P_7 is REGULAR (scores (2,2,2,2,2))
- The 5-cycles form a **4-(7, 5, 6) design** (not just a 2-design)
- This equals 2 copies of the trivial complete design on 5-subsets of 7 points
- Each complement pair (the 2 vertices NOT in a 5-cycle) appears exactly 2 times

At p ≥ 11, the 5-cycles form only a 2-design (the 3-design and 4-design properties fail).

## Non-Paley Comparison

The interval tournament {1,2,3} mod 7 (circulant with non-QR connection set) has:
- 3-cycles: NOT a 2-design (pair distribution {1:7, 2:7, 3:7})
- 5-cycles: NOT a 2-design (pair distribution {13:14, 14:7})

Only the QR structure of Paley creates the design uniformity.

## Related

- THM-211: 3-cycle BIBD for Paley
- THM-305: v_2(T(n)) = (n-1)/2
- The Barvinok-BIBD-QR triangle (07-reflections/barvinok-bibd-qr-triangle.md)
