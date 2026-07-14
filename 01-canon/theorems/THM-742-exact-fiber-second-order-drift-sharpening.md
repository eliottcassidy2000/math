---
id: THM-742
title: The exact-fiber / second-order drift sharpening of THM-739/740 — L(V) is EXACTLY the slope-W geodesic measure of the polygonal torus region R (strand identity, verified as Fraction equality); wedge/sawtooth analysis replaces the adversarial constants — first-order constant 146 → 7.25 (shape 1), W₀ 1948 → 452 and 2676 → 584; THM-740's inner-ride C_in 44.7 → 21.8
status: PROVED (proof below; stated constants carry a deliberate 2× margin over the derived ones) + VERIFIED-EXACT (strand identity EQUAL as Fractions; zero bound violations over the W-spread at both derived and stated constants, both shapes)
source: opus-2026-07-13-S275 (owner prompt: prove the second-order drift sharpening for THM-739 and 740)
depends_on:
  - THM-739 (the object being sharpened; its Area and A_J machinery)
related:
  - THM-740 (inner ride sharpened here), THM-737, kps THM-735
  - the density route's Q_s (THM-729): the REMAINING gap (observed |L−Area| is ~50–500× below even the sharpened bound) is cross-line joint equidistribution — the same cancellation class
---

# THM-742 — the exact-fiber / second-order drift sharpening

## 1. The reformulation (an identity, not an estimate)

For `V = B ∪ (W+J)`, write `s = {Wt}` and define the **fixed** polygonal torus region

> `R = {(u, σ) ∈ T² : u ∈ G_B, ‖σ + ju‖ ≥ 1/14 ∀j ∈ J}`,  `Area(R) = ∫_{G_B} A_J = Area(B,J)`.

Then **exactly** (change of variables, no approximation):

> `L(V) = (1/W) Σ_{m=0}^{W−1} f(m)`,  `f(m) = |{σ : ((m+σ)/W, σ) ∈ R}|`

— the measure of the slope-`W` geodesic inside `R`. (Verified as an exact Fraction equality at
`W = 30` for both shapes: engine `L` = strand sum, EQUAL.) THM-739's error is therefore precisely
a strand-vs-area discrepancy for a region bounded by segments of the `2|J|` lines
`σ = −ju ∓ 1/14` (slope `−j`) plus `2#comp(G_B)` vertical segments.

## 2. The sharpened bound

Let the **arrangement data** of `R` be: the maximal *exposed segments* per sloped line (a
segment = a maximal `u`-run where the arc-`j` endpoint lies on `∂(∪ arcs)` inside `G_B`), with
`Ξ = Σ_segments min(j·Δu_seg/2, 1/4)`; and `V_R` = the number of arrangement breakpoints inside
`G_B`. Then with

- `C₁ = 2·(#comp(G_B) + Ξ)`   (stated; derived `#comp + Ξ`),
- `C₂ = 6·Σ_{j∈J} j² + 4·max(J)·V_R`   (stated; derived `3Σj² + 2·max(J)·V_R`),

> **`|L(V) − Area(B,J)| ≤ C₁/W + C₂/W²`.**

## 3. Proof

Compare `f(m)` with the strip average `g(m) = W·∫_{strip m} |R_u| du` (`Σ_m g = W·Area`).

1. **Strips crossed by no boundary:** `R ∩ strip` is a union of full-width bands ⟹ `f = g`
   exactly. The `V(A_J)` Riemann term of THM-739 is *gone* (we never sample `A_J`).
2. **Vertical edges** (at the `2#comp(G_B)` points of `∂G_B`): the fiber measure jumps by
   `≤ 1`; the strand-vs-average error in that strip is `≤ ½·jump`. Total `≤ #comp(G_B)`.
   The base-fattening term `2Σ_B b` of THM-739 is *replaced* by this (the `Σ_B b` many danger
   arcs interior to the complement never mattered — only the surviving boundary does).
3. **Sloped edges — the signed wedge.** A slope-`−j` edge crossing strip `m` at entry height
   `s₀` contributes `f − g = (j/W)·ψ(s₀) + O((j/W)²)`, `ψ(x) = 1/2 − x` (exact computation of
   the diagonal crossing vs the strip average of a linear boundary). The quadratic remainders
   total `≤ Σ_j 2j²/W` over each line's `≤ W+j` crossings (`→ C₂`).
4. **Sawtooth cancellation along a line.** Over the crossings of one maximal segment, the entry
   heights march by `j/W` steps: `Σψ` over a *full wrap* collapses by the distribution relation
   `Σ_{k mod q}ψ({x+k/q}) = ψ({qx})` (cost `≤ gcd(j,W)/2` per wrap-set, `→ C₂` as `≤ j²`);
   a *partial wrap* (segment end) costs at most the maximal partial sum of a monotone sweep:
   `(j/W)·max_r [r(1/2−x) − r²j/2W] ≤ (j/W)·W/(8j) = 1/8`. Hence per segment
   `≤ min(j·Δu/2, 1/4)` (short segments: crossing-count × `|ψ| ≤ 1/2`; long: two `1/8`-caps).
   Summing: `Ξ`.
5. **Vertices** (breakpoints of the arrangement): interval-endpoint annihilation makes the
   per-edge linear model overcount by at most the local wedge sizes, `≤ 2max(J)/W` per vertex;
   total `≤ 2max(J)·V_R/W` (`→ C₂`). ∎

## 4. Verified numbers

| shape | `Ξ` (#segs) | derived `C₁` | stated `C₁` | crude 739 `C` | `C₂` (stated) | **new `W₀`** | old `W₀` |
|---|---|---|---|---|---|---|---|
| 1: `{1}∪{W..W+11}` | 6.246 (172) | **7.25** | 14.49 | 146 | 8712 | **452** | 1948 |
| 2: `{1,2,3}∪{W..W+9}` | 5.569 (154) | **9.57** | 19.14 | 120 | 4086 | **584** | 2676 |

Zero violations over `W ∈ {10,…,800}` at both constant sets. The first-order constant improved
**20×**; the overall `W₀` improved **4.3–4.6×** and is now *bound by the second-order vertex
term* `4max(J)V_R` — the next constant to attack, along with cross-line cancellation.

**THM-740 inner ride sharpened:** the same strand-vs-strip analysis with the continuous weight
`A_{J₂}` (adds only its variation, no jumps): `C_in^sharp = 2(#comp(G_B)+Ξ(J₁)) + V(A_{J₂}) =
21.81` vs crude `44.7` (shape A; verified against the exact S274 values, all valid). The outer
ride's `O(W₁)` terms (cluster-1 base-fattening and `#comp(G₁)`) are *genuine* (many components
⟹ many real jumps), so THM-740's linear cone ratio stands.

## 5. Honest gap

Observed `|L − Area|` is still ~50–500× below even the *derived* bound: the remaining
cancellation is **across lines** (the `|J|` sawtooth families at different `j` are jointly
equidistributed, not adversarially aligned) — a joint/multi-frequency estimate of exactly the
`Q_s`/Weyl class the density route owns (THM-729). The per-line analysis here is the limit of
what single-frequency arguments give.

## Files

`04-computation/lrc14_second_order_drift_thm742_opus_S275.py` (+ `.out`): strand-identity
Fraction check, exact arrangement/Ξ/V_R computation, bound verification over the W-spread
(zero violations), new-W₀ solve, THM-740 inner-ride check.
