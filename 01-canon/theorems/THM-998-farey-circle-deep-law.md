# THM-998 — The Farey-circle deep law: the deep set is a major-arc dissection (boxeph-2026-07-17-S81)

**Status:** PROVED (structure + arc geometry; verified across N=11,12,14,20 and caps
K=2..6 — `lrc_farey_circles_boxeph_S81.py`, `lrc_twocircle_generalize_boxeph_S81.py`).
**Explains and generalizes death-star's THM-985/987 "two-circle" theorem**: the "two
circles" of the N=14, K=6 census are the `b=1` and `b=2` slices of a Farey dissection.

## Statement

Fix `N ≥ 2`, prefix family `V = {1,…,N-1}`, threshold `1/N`, and a coverage cap
`K ≥ 1`. Define the **K-deep set** `D_K(q) = {p ∈ (0,q) : bandCount(V,q,p) ≥ K}` (at
least `K` runners within `1/N` of the origin). Then, for `q` large enough that the arcs
separate, `D_K(q)` is a **disjoint union of resonance arcs, one around each Farey
fraction `a/b` (lowest terms, `0 ≤ a/b ≤ 1`) with**

> **`b ≤ (N-1)/K`**   (equivalently `⌊(N-1)/b⌋ ≥ K`: at least `K` speeds are `≡ 0 mod b`).

Around `a/b`, the multiples of `b` in `V` nest at `0`; the binding one is speed `bK`, so
the arc is `{p : ‖p/q − a/b‖ < 1/(NbK)}`, i.e. half-width `≈ q/(NbK)` about `p = aq/b`.

**Circle count** `= #{b ≥ 1 : b ≤ (N-1)/K} = ⌊(N-1)/K⌋`. The Farey fractions with
denominator `b` number `φ(b)` (plus the two `b=1` endpoints `0, q`).

## Why death-star found exactly two circles (N=14, K=6)

`⌊13/6⌋ = 2`, so only `b=1` (the **integer circle**, small `p` / `p≈q`) and `b=2` (the
**half circle**, `p≈q/2`, the even speeds) survive. `b=3` needs `⌊13/3⌋ = 4 < 6` speeds
— absent. The magic constant `84 = 6·14 = K·N` is the binding threshold `bK·N` at `b=1`.
Lowering K reveals the hidden circles (verified): K=4→3 circles (adds `b=3`, centers
`q/3, 2q/3`); K=3→4 (adds `b=4`, centers `q/4, 3q/4`); K=2→6.

## Deep-count formula (uniform)

`#D_K(q) = Σ_{b ≤ (N-1)/K} (arc count at denominator b)`. For `b=1`: two boundary arcs
of width `B := ⌊(q-1)/(NK)⌋` → `2B`. For each `b ≥ 2`: `φ(b)` central arcs. death-star's
`#deep(q) = 2B + (B+1−(q+B)%2)` is exactly `b=1` (`2B`) plus the single `b=2` arc
(`φ(2)=1`, parity-corrected width). The formula generalizes term-by-term.

## Significance — the circle method, discretely

- The deep set is the **major-arc dissection** of the circle method for the tight family;
  "wagner-circles" = Farey/Ford resonance circles. This grounds the census machinery in
  classical analytic number theory (ties to `cuts-as-farey-geodesics`,
  `doubled-primes-as-the-parity-hinge`).
- **Covering vs non-covering, in Farey terms** (the LRC(14) dichotomy, reinterpreted):
  a family has a deep circle at denominator `b` iff it has `≥ K` speeds `≡ 0 mod b`.
  - *Non-covering* (omits a multiple of some `q₀ ≤ N`): the circle at `1/q₀` is **empty**
    → its center is deep-free → live witness `t = 1/q₀` (this **is** the non-covering
    reduction, THM-366/523, re-read as an empty minor arc).
  - *Covering* (a multiple of every `b ≤ N`): **every** small-denominator circle is
    non-empty → no empty center among small `b` → the deep-free set retreats to genuine
    minor arcs → genuine harmonic analysis. This is the LRC(14) crux, located.
- **Live witnesses = centers of empty resonance circles.** The tight family's only empty
  small circles are at denominator exactly `N` (units) — [[THM-996-resonance-confinement-census-law-deathstar-S56]] (death-star; its Part II is this live law),
  [[THM-997-resonant-dichotomy]] — which is why it is the hardest (tightest) equality case.

Related: [[THM-985]], [[THM-987]] (N=14 two-circle, death-star), [[THM-996-resonance-confinement-census-law-deathstar-S56]] (death-star; its Part II is this live law),
[[THM-997-resonant-dichotomy]], HYP-7315.
