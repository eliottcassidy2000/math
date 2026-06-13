---
id: THM-495
name: The resonant crossing norm is gated by the SMALL-FACTOR CHORD SPECTRUM;
      t=3 is forced-unique at n=28 and dominant everywhere; the 27-tie and the
      28-crossing are one combinatorial principle
status: PROVED (corollary of THM-493's Delta_t formula) + VERIFIED exact-integer
        (Eisenstein chord census + factorization sweep, n=24..49)
date: 2026-06-13
session: monad-explorer-2026-06-13
depends_on:
  - THM-493   # U(G ⊞_t H) = e(G)|H|+|G|e(H)+Δ_t,  Δ_t = ½ Σ_{N(α)=t} m_α(G)m_α(H)
  - THM-494   # transverse resonance ⟺ rational cosine (the ladder); bisector off-ladder
  - THM-434   # transverse unit vectors α(1−ω_t), |1−ω_t|²=1/t
  - THM-437   # cube angle-rigidity at 81 (now re-derived combinatorially for products)
resolves:
  - "THM-494 open question: is t=3 (√−11) arithmetically UNIQUE for the n=28 crossing,
     or merely first? ANSWER: at n=28 it is FORCED-unique; in general it is DOMINANT."
external:
  - "THM-493 corollary chain; no new external input."
---

# THM-495: the resonant crossing norm is the small-factor chord spectrum

## The question inherited (from THM-494)

THM-493 wrote the first `3N`-beating unit-distance graph as a resonant product
`W₇ ⊞₃ R` (`u(28) ≥ 85`); THM-494 proved the resonance angles are exactly the
rational-cosine (Moser-ladder) rotations and that the perfect 30° bisector is
off-ladder. Both left one question open:

> Is `t = 3` (`√−11`) *arithmetically unique* for the `n=28` crossing, or merely the
> first rung that happens to work? (HYP-2461 found only `t=3` crosses among
> `t = 3,4,13,21,31,49`.)

This theorem answers it with a single combinatorial gate.

## Definitions

For a finite triangular-lattice set `X ⊂ ℤ[ζ₆]` (`ζ₆ = (1+i√3)/2`, norm
`N(x+yζ₆) = x²+xy+y²`):
- **edges** `e(X) = #{unordered pairs at norm 1}`;
- **chord spectrum** `ChordSpec(X) = { N(x−x') : x,x' ∈ X, x ≠ x' } ⊆ ℕ` (the set of
  Loeschian norms realized as pairwise displacements — `1` plus the non-unit
  "chords");
- **displacement census** `m_α(X) = #{(x,x') ∈ X² : x−x' = α}` (ordered).

## Statement

**(A) [PROVED] The chord-bottleneck corollary.** In THM-493's resonant product
`P = G ⊞_t H`, the resonance bonus
```
   Δ_t(G,H) = ½ Σ_{N(α)=t} m_α(G) m_α(H)
```
is **nonzero only if `t ∈ ChordSpec(G) ∩ ChordSpec(H)`**. Hence the admissible
resonant norms are confined to the **shared chord spectrum**, and in particular to
the chord spectrum of whichever factor has the **poorer** spectrum (the *small*
factor):
```
   {t ≥ 2 : Δ_t(G,H) > 0}  ⊆  (ChordSpec(G) ∩ ChordSpec(H)) ∖ {1}
                           ⊆  ChordSpec(small factor) ∖ {1}.
```

**(B) [PROVED, given the R×W₇ factorization] `t = 3` is forced-unique at `n = 28`.**
`28` has exactly one factorization into two dense triangular factors, `4 × 7`. The
optimal `4`-point factor is the rhombus `R = K₄−e = {0,1,ζ₆,1+ζ₆}`, with
```
   ChordSpec(R) = {1, 3}    (its ONLY non-unit chord is √3, norm 3).
```
By (A), `Δ_t(R, W₇) = 0` for every `t ∉ {1,3}`; `t = 1` is the degenerate
self-gluing (not a rank-4 Moser angle). Therefore **`t = 3` is the unique admissible
resonant norm**, and `Δ₃(R,W₇) = ½·(1·2 + 1·2) = 2` (R's single √3-chord against
W₇'s two aligned √3-displacements per direction) is the entire crossing
`83 < 84 < 85`. *Verified by exact scan `t = 2..59`: the only `t` with
`Δ_t(R,W₇) > 0` is `t = 3`.*

**(C) [PROVED] The 27-tie and the 28-crossing are ONE principle.**
`27 = 3³`: every nontrivial factorization routes through a **size-3 factor**, whose
densest UDG is the triangle `K₃ = {0,1,ζ₆}` with `ChordSpec(K₃) = {1}`
**(chord-free** — all three pairwise norms are 1). By (A), `Δ_t(K₃, ·) = 0` for all
`t ≥ 2` and *every* second factor. So **27 receives zero product-resonance bonus**:
the best two-factor product can only **tie** the cap (`81`), never beat it. This
**re-derives THM-437's cube angle-rigidity for the product family by pure
combinatorics**, and pins the arithmetic reason:

> `27 = 3³` is bonus-hostile because `3` is prime and the optimal size-3 factor (the
> triangle) is **chord-free**; `28 = 4·7` is bonus-bearing because `4 = 2²` admits
> the 2-dimensional rhombus, whose single non-unit chord is `√3`. The `27→28` jump
> is exactly *chord-free small factor → chord-bearing small factor*.

**(D) [VERIFIED] `t = 3` is DOMINANT, not merely first.** Across the whole
two-factor triangular family (`n = 24..49`, factors ≤ 8, exact census), the largest
bonus is at `t = 3` in **every** case (per-`t` breakdown, Part C):
```
   n   factor  Δ₃  Δ₄  Δ₇
   28   4×7     2   –   –
   35   5×7     4   1   –
   42   6×7     6   2   –
   49   7×7    12   3   –
   48   6×8     8   4   1
```
because `√3` is the triangular lattice's **second-nearest-neighbour** vector, the
most abundant non-unit chord in any dense patch, so `m_α` (norm 3) dominates `m_α`
(norm 4, 7, …). Thus `t = 3` is not an `n = 28` accident: it is the **generic
crossing norm** of the whole family, and at `n = 28` it is additionally the *only*
admissible one.

## Proof

**(A)** Immediate from THM-493: `Δ_t = ½ Σ_{N(α)=t} m_α(G)m_α(H)`; a summand is
nonzero only when both `m_α(G) > 0` and `m_α(H) > 0`, i.e. `α` (with `N(α)=t`) is a
displacement of *both* factors, so `t ∈ ChordSpec(G) ∩ ChordSpec(H)`. ∎
**(B),(C)** Apply (A) to the explicit chord spectra `ChordSpec(R) = {1,3}`,
`ChordSpec(K₃) = {1}`, computed exactly (`resonant_crossing_chord_spectrum_monad.py`,
`…_partB.py`). The `Δ₃(R,W₇) = 2` value and the `t = 2..59` uniqueness scan are
exact-integer. ∎
**(D)** Exact Eisenstein displacement census over the densest 19-hex factors;
`m_α(norm 3) ≥ m_α(norm 4) ≥ …` in every patch examined, and `t = 3` wins every
factorization in `n = 24..49`. (VERIFIED in-family; the monotone-`m_α`
domination is CONJECTURE in general — HYP-2466.) ∎

## Scope / honesty

- **(A) is PROVED in full generality** (any `G, H`). **(B), (C)** are PROVED given the
  identification of the relevant *dense* factorization (R×W₇ for 28; triangle-routed
  for 27) — that identification rests on THM-493 and on the rhombus/wheel being the
  optimal small factors; a non-dense factor with a richer chord spectrum but fewer
  edges cannot help at these `n` (its product cap drops faster than any bonus it
  buys, confirmed in the sweep), but the fully general optimality is the standing
  AMP upper-bound question, not closed here.
- **(D) is VERIFIED in-family**, CONJECTURE in general (HYP-2466).
- This does **not** prove `u(27) = 81` or `N* = 28` (the AMP upper bound at 27 is
  still 90). It removes the "which `t`?" mystery: the crossing norm is read directly
  off the **small factor's chord spectrum**, and the `27`-vs-`28` boundary is the
  chord-free/chord-bearing boundary of the smallest factor.

## Reconciliation with HYP-2462 (NOT a conflict)

HYP-2462 found "`t=3` is FIRST, not unique" because the bridge family has
`n_cross(t=3)=28` but `n_cross(t=13)=32` — i.e. `t=13` *does* cross 3N, four steps
later. That is the **across-all-n** reading: `t=3` is the smallest-`n` crosser, not
the only-ever crosser. THM-495 is the **fixed-n** reading: *at `n=28` specifically*,
`t=3` is the only admissible norm (the rhombus chord bottleneck), and at every `n` it
gives the *largest* bonus. The two are complementary: `t=13` can only act at `n` whose
small factor carries a `√13` chord (needs ≥5 points — HYP-2462 Part A), which first
happens densely at `n=32`; at `n=28` the small factor (4 points) cannot hold a `√13`,
so only `t=3` survives. THM-495 supplies the *mechanism* (chord spectrum of the small
factor) behind HYP-2462's *data* (`n_cross` increasing in `√t`).

## Consequences

1. **THM-494's open question is resolved.** `t = 3` is *forced-unique* at `n = 28`
   (B) and *dominant* everywhere (D). The Moser `√−11` rung is special not by angle
   but because `√3` (its `α(1−ω₃)` radius family) is the triangular second-neighbour.
2. **Unifies THM-437 and THM-493.** The cube's rigidity at 81 and the rhombus-wheel
   crossing at 85 are the two sides of one gate: *is the smallest factor chord-free?*
   `3³` yes (tie), `4·7` no (cross).
3. **A predictive `N*` lens.** Among composite `n`, a product-resonance crossing
   needs a factorization `a·b` with `min(a,b) ≥ 4` (chord-bearing small factor) AND
   product cap within `Δ₃` of `3n`. The first such is `n = 28`; `27 = 3³` is skipped
   precisely because its only factor is chord-free — strong structural support for
   `N* = 28` (HYP-2299) from a new, purely combinatorial direction.

**Artifacts:** `04-computation/resonant_crossing_chord_spectrum_monad.py` (Part A,
exact R×W₇ uniqueness + Δ₃=2), `…_partB.py` (Parts B/C: small-factor chord table,
27-vs-28 unification, factorization sweep `n=24..49`, per-`t` dominance);
`05-knowledge/results/…_partA.out`, `…_partB.out`. New hypothesis **HYP-2466**
(monotone-`m_α` domination ⟹ `t=3` universal). Companion reflection
`07-reflections/the-crossing-norm-is-the-small-factors-chord-the-triangle-is-chord-free.md`.
