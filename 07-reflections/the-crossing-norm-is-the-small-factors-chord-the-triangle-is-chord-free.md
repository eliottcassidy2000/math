---
source: monad-explorer-2026-06-13 (deep-research; OPEN-Q-057 frontier; THM-494 handoff)
status: REFLECTION grounded in PROVED arithmetic (THM-495(A), a corollary of the
  THM-493 Δ_t formula) + exact-integer Eisenstein chord census and factorization
  sweep (n=24..49, every count exact).
tags: [unit-distance, N-star, 3N-crossover, resonant-product, THM-493, THM-494,
  THM-495, THM-437, chord-spectrum, eisenstein, moser-lattice, sqrt3, second-neighbour,
  OPEN-Q-057, HYP-2461, HYP-2466]
---

# The crossing norm is the small factor's chord — and the triangle is chord-free

## The pebble in the shoe

For three same-day theorems the `u(28) ≥ 85` unit-distance crossing kept being
pinned to one arithmetic accident: it crosses **only** in the Moser lattice `L₃`
(`√−11`). THM-493 wrote `85` as the product `W₇ ⊞₃ R` plus a resonance bonus
`Δ₃ = 2`. THM-494 proved the resonance angles are exactly the rational-cosine
rotations and that even the geometrically perfect 30° bisector is off the ladder and
barren. Every result said *t = 3 is special* — but none said **why 3, and whether 3
is unique or merely first**. That was the pebble. It turns out to dissolve in one
line.

## The one line

THM-493's bonus is `Δ_t(G,H) = ½ Σ_{N(α)=t} m_α(G) m_α(H)`. A summand survives only
when the displacement `α` (Eisenstein norm `t`) occurs in **both** factors. So:

> **The resonance bonus at norm `t` is nonzero only if `t` is a *shared chord* of
> both factors.** The admissible crossing norms are confined to the chord spectrum of
> the *smaller* factor.

That is the whole content. No angle calculus, no class numbers — just "does the
small factor have a `√t` chord?" Everything else is reading off a tiny table of small
triangular-lattice point sets:

```
   factor size   densest UDG       ChordSpec (Loeschian norms of displacements)
       2          edge K₂          {1}                 ← chord-free
       3          triangle K₃      {1}                 ← chord-free
       4          rhombus K₄−e     {1, 3}              ← first non-unit chord: √3
       5          (rhombus+1)      {1, 3, 4}
       6          T₃ / 2×3         {1, 3, 4, 7}
       7          wheel W₇         {1, 3, 4}
```

## What it explains, all at once

**Why `t = 3` at `n = 28`, uniquely.** `28 = 4·7` is the only dense factorization.
The 4-factor is the rhombus, whose *only* non-unit chord is `√3` (norm 3). So
`Δ_t(R, W₇) = 0` for every `t ≠ 3` — not "t=3 happens to win," but *t=3 is the only
admissible norm at all*. Exact scan `t = 2..59`: a single survivor, `t = 3`,
`Δ₃ = 2`. Forced-unique.

**Why `n = 27` cannot cross, ever (via products).** `27 = 3³`. Every nontrivial
factorization routes through a size-3 factor, and the densest size-3 UDG is the
**triangle**, whose three chords are *all unit edges*: `ChordSpec(K₃) = {1}`,
chord-free. So the bonus is identically zero through any triangle factor — `27` can
only tie the product cap (`81`), never beat it. This is **THM-437's cube
angle-rigidity, re-derived by counting chords** instead of solving
`cos u + cos w + cos(u−w) = −1`. And it names the reason: `3` is prime, and the
prime-3 optimal factor is chord-free.

So the entire `27 → 28` boundary — the exact location of `N*` — is one combinatorial
dichotomy:

> **chord-free smallest factor (tie) vs chord-bearing smallest factor (cross).**
> `3³`: the smallest factor is the chord-free triangle. `4·7`: the smallest factor is
> the chord-bearing rhombus. The crossing happens at the first `n` whose smallest
> dense factor carries a `√3`.

## The deeper shape: `t = 3` is the lattice's own second neighbour

The follow-up — *is 3 merely first?* — has a sharper answer than "first." Across the
whole two-factor family (`n = 24..49`), the bonus is **largest at `t = 3` in every
case** (`Δ₃ ≥ Δ₄ ≥ Δ₇`). The reason is geometric: `√3` is the **second-nearest-
neighbour distance** of the triangular lattice (after the unit edge). Any dense patch
is thick with `√3` pairs, so `m_α(norm 3)` dominates `m_α(norm 4), m_α(norm 7), …`.
The Moser `√−11` rung is privileged because its transverse radius family `α(1−ω₃)`
has Eisenstein norm `N(α) = 3` — it resonates with the lattice's *own* second
neighbour. `t = 3` is therefore not the first rung of a ladder we climb; it is the
**ground tone** of the triangular lattice, and the higher rungs are overtones that
never sound as loud.

## The transcending pattern

This is the same lesson the project keeps meeting, sharpened once more. The
threshold-breaking structure looked *angular* (a Moser rotation, a magic `33.6°`),
then *arithmetic* (a rational cosine, a class-number texture). It is really
**neighbourly**: the crossing rides the lattice's shortest non-trivial chord, and the
obstruction at `27` is simply that its arithmetic (`3³`) forces a factor too small to
*have* a chord. Geometry → arithmetic → **combinatorics of the second neighbour**.
The hardest-looking gate in the `N*` problem is, at bottom, the question *"is your
smallest piece big enough to hold a √3?"*

Two pieces of unfinished business this opens (HYP-2466):
1. Prove the `m_α`-domination (`norm 3` displacement count ≥ all higher norms) for
   *every* dense triangular patch — that would upgrade "`t=3` dominant in-family
   (VERIFIED)" to a theorem.
2. The chord-bottleneck is a *product* statement. The free (non-product) crossings
   AMP/Engel actually use are denser; do they still ride `√3`? HYP-2461's free-patch
   data (only `t=3` crosses) says yes — formalize the bridge from product chords to
   free-patch chords.

Cross-links: [[THM-493]] (the Δ_t formula this corollarizes), [[THM-494]] (the ladder
is rational-cosine; bisector off-ladder), [[THM-437]] (cube rigidity, now a chord
count), `the-perfect-bisector-is-off-the-ladder-rational-cosine-not-rational-angle.md`,
`the-unit-distance-tie-is-carrier-robust-the-crossing-is-resonant.md`.
