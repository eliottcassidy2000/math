---
source: opus-2026-06-06-S699 (χ(ℝ²) narrowing attempts)
status: NO ELIMINATION of {5,6,7} (a major open problem) — but a useful rigorous META-statement + reframings. (1) χ_f(ℝ²)=1/m₁ ≤ 4.36 < 5, so the fractional/spectral/measure bound CANNOT reach 5: narrowing {5,6,7} is IRREDUCIBLY COMBINATORIAL = the LRC Vitali-wall integrality gap (rigorous, given known m₁). (2) LOESCHIAN reframing: 5,6 are not norms a²+ab+b², so no Eisenstein-hexagonal periodic coloring uses 5 or 6; 7=N(2+ζ6) is the smallest ⟹ the periodic upper bound is exactly 7, and a 6-coloring must be non-Eisenstein. (3) Field-tower: χ climbs ℚ²=2, ℤ[ζ6]=3, +√−11≥4, →[5,7]; both ends of [5,7] escape the cyclotomic floor 3. (4) Heegner roadmap: χ = 2 + (#independent Heegner rotations).
tags: [hadwiger-nelson, chromatic-plane, fractional-chromatic, m1, integrality-gap, vitali-wall, loeschian, eisenstein, 7-coloring, field-tower, heegner, honest-negative]
---

# χ(ℝ²) ∈ {5,6,7}: narrowing attempts, and why it's a combinatorial wall

**Prompt (user):** spend a long session on χ(ℝ²); attempt any novel statement, even small, like
eliminating one of {5,6,7}.

**Honest headline: I did not eliminate any of 5, 6, 7.** Doing so means exhibiting a 6-coloring
(kills 7), a 6-chromatic unit-distance graph (kills 5), or a 5-coloring (kills 7+6) — each a major
open result. What the session *did* produce: one rigorous meta-statement on *why* it's hard and
*which techniques can't work*, plus reframings and a roadmap.

## 1. The useful rigorous statement: narrowing is irreducibly combinatorial

> **`χ_f(ℝ²) = 1/m₁ ≤ 4.36 < 5 ≤ χ(ℝ²)`**, where `m₁ ∈ [0.2293, 0.2598]` is the plane's
> 1-avoiding density (known bounds). So the **fractional chromatic number is strictly below 5**,
> and the integrality gap `χ − χ_f ≥ 0.64`.

Consequences (rigorous, given the `m₁` bounds):
- **No fractional / spectral / measure / density argument can even reach `χ ≥ 5`** — they are all
  bounded by `χ_f ≤ 4.36`. The `χ ≥ 5` lower bound (de Grey) is *necessarily combinatorial*
  (a finite gadget), and **distinguishing 5 vs 6 vs 7 lives entirely in the integrality gap.**
- **This integrality gap IS the LRC "Vitali wall"** (S699g, THM-406 M2): the gap between the
  *measurable/fractional* density bound and the *true combinatorial* value. So HN's `{5,6,7}`
  uncertainty and LRC's worry-set both sit in the same place — the measure-blind combinatorial
  residual. *Any narrowing of `{5,6,7}` must be combinatorial, not analytic* — a clean honest
  constraint on the problem.

## 2. The Loeschian reframing of the upper bound (why 7, and what a 6 would require)

The Isbell 7-coloring is the norm-`7` Eisenstein coloring `ℤ[ζ₆]/(2+ζ₆) ≅ 𝔽₇` (S699m). The
periodic-coloring index must be a **Loeschian number** `a²+ab+b²` (the Eisenstein norms):
`1,3,4,7,9,12,13,…`. **Verified: `5` and `6` are NOT Loeschian; `7 = N(2+ζ₆)` is the smallest
> 4.** Therefore:

> **No Eisenstein-hexagonal periodic coloring uses 5 or 6 colors; the smallest is 7.** Equivalently,
> **a 6-coloring of the plane (if one exists) must be non-Eisenstein** — aperiodic, or on a
> different lattice. This is a structural constraint on any upper-bound improvement: *to beat 7 you
> must leave the hexagonal/cyclotomic colorings.* (It does not bound the true `χ` — non-lattice
> colorings are not excluded — so it is a statement about the *natural* coloring family, not an
> elimination.)

## 3. The field-tower: both ends of [5,7] escape the cyclotomic floor

Chromatic number climbs along the number-field tower (S687 extended):
```
   χ(ℚ²) = 2   (rational plane bipartite)
   χ(ℤ[ζ₆]) = 3   (Eisenstein lattice; W₆ forces 3, 3-colors via mod √−3 = 𝔽₃)   ← the cyclotomic FLOOR
   χ(ℚ(√−3,√−11)) ≥ 4   (contains the Moser spindle; √−11 = the non-cyclotomic rotation)
   χ(ℝ²) ∈ {5,6,7}
```
> **The cyclotomic floor is `3`; `[5,7]` is the post-cyclotomic regime, and BOTH ends escape the
> lattice.** The *lower* bound (`χ≥5`) escapes by forcing Heegner rotations (S687); the *upper*
> bound (`≤6`) would require a non-Eisenstein coloring (§2). The lattice pins only `χ=3`; the true
> value is an *escape* statement on both sides — which is exactly why it sits in the
> integrality-gap/Vitali wall (§1).

## 4. The Heegner roadmap (the only concrete path to narrowing)

S699m's conjecture: `χ = 3 ↦ √−3`, `χ = 4 ↦ √−11`, `χ = 5 ↦ √−19?` (Heegner, class number 1).
> **If a finite unit-distance graph forces at most `k` independent Heegner imaginary-quadratic
> rotations, then `χ(ℝ²) = 2 + k`.** Narrowing `{5,6,7}` ⟺ **bounding the maximal Heegner-rotation
> rank of a finite unit-distance graph.** A 6-chromatic graph (kills 5) would need rank-4 Heegner
> rotations; ruling that out (kills 6,7, leaving 5) would need a rank cap of 3. This is a concrete,
> arithmetic reformulation of the lower-bound side — the one place a *finite* computation could
> bite. (Conjectural; the verified content is the `√−3, √−11` realizations.)

## 5. Honest status

- **Not done:** eliminating any of `5,6,7` (the open problem).
- **Rigorous (given known `m₁`):** `χ_f ≤ 4.36 < 5` ⟹ narrowing is irreducibly combinatorial = the
  Vitali-wall integrality gap; `5,6` not Loeschian ⟹ no Eisenstein-hexagonal 5/6-coloring.
- **Structural framings:** the field-tower chromatic growth (`2,3,4,…,[5,7]`); both ends of `[5,7]`
  escape the cyclotomic floor `3`; a 6-coloring must be non-Eisenstein.
- **Roadmap (conjectural):** `χ = 2 + #Heegner rotations`; narrowing ⟺ a Heegner-rank cap on finite
  unit-distance graphs.
- **The takeaway:** the value of `χ(ℝ²)` is locked in the *combinatorial integrality gap* — no
  analytic/spectral method can touch it (§1) — so progress must come from finite gadgets, and the
  cleanest finite handle is the Heegner-rotation rank (§4). This both explains the difficulty and
  points at the only door.

**Artifacts:** `04-computation/hn_narrowing_attempts_s699n.py` (+`.out`). Builds on S699g (spectral
unification / Vitali wall = integrality gap), S699m (`z⁷−z`/`𝔽₇`/Heegner), S687 (field tower),
THM-406 (Vitali wall), Falconer/de Grey/Isbell. New: **HYP-2278**.
