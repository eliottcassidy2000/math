# HYP-7011 — The 6617 identity: the clock flat set is the DOUBLING IMAGE of the corridor

**Status:** PROVED-EXACT (death-star-2026-07-16-S17; exact rational interval identity,
`lrc14_flat_vs_corridor_6617_deathstar_S17.py`, referee output in results/).
**Resolves:** mac-mini-S113 backlog item 2 / kps cont.28 handoff ("flat 6617/97020 =
2 × corridor 6617/194040 — candidate mechanism: the ±/2-fold S2-kink symmetry").

## Statement

Let `S2(x) = Σ_{d=1}^{12} (13−d)·max(0, 1/7 − ‖dx‖)` (pair-overlap energy of the tight AP
{1..13}; THM-878), `F = {x : S2(x) = 6/7}` its flat set, and
`G = {y : ‖vy‖ ≥ 1/14 ∀v ∈ {1..12}}` the deep-well corridor good set (THM-853, value
`m(G) = 6617/194040`). Then, at the level of positive-length components:

> **`F = 2·G (mod 1)`** — the doubling map `y ↦ 2y mod 1` sends the 12 components of `G`
> bijectively onto the 12 components of `F` (verified exact, no overlaps), and `G` contains
> no pair `{y, y+1/2}`, so the map is injective on `G` and
> **`λ(F) = 2·m(G) = 6617/97020`** — the factor 2 is the JACOBIAN of the doubling map.

Even the four isolated (measure-zero) points of `G` — `3/14, 5/14, 9/14, 11/14` — double to
isolated flat points of `S2` (`3/7, 5/7, 2/7, 4/7`, all with `S2 = 6/7` exactly).

## Component table (exact)

| G component | len | ×2 → F component | len |
|---|---|---|---|
| [1/14, 13/168] | 1/168 | [1/7, 13/84] | 1/84 |
| [15/98, 13/84] | 1/588 | [15/49, 13/42] | 1/294 |
| [29/126, 13/56] | 1/504 | [29/63, 13/28] | 1/252 |
| [43/140, 13/42] | 1/420 | [43/70, 13/21] | 1/210 |
| [43/112, 27/70] | 1/560 | [43/56, 27/35] | 1/280 |
| [71/154, 13/28] | 1/308 | [71/77, 13/14] | 1/154 |
| [15/28, 83/154] | 1/308 | [1/14, 6/77] | 1/154 |
| [43/70, 69/112] | 1/560 | [8/35, 13/56] | 1/280 |
| [29/42, 97/140] | 1/420 | [8/21, 27/70] | 1/210 |
| [43/56, 97/126] | 1/504 | [15/28, 34/63] | 1/252 |
| [71/84, 83/98] | 1/588 | [29/42, 34/49] | 1/294 |
| [155/168, 13/14] | 1/168 | [71/84, 6/7] | 1/84 |

(Right endpoints of G-components are `13/(14k)`-type; F right-endpoints `13/(7k)`-type —
the doubling in plain sight.)

## Mechanism reading

`‖d·(2y)‖ = ‖(2d)·y‖`: the flat-set condition at `x = 2y` reads the EVEN sublattice of
speeds at `y`, and the corridor condition supplies exactly the clearances needed — the
identity is the **2-adic descent deck map** (THM-580's `t → 2t`) connecting the Fejes-Tóth
world (THM-878) to the corridor world (THM-853). NOTE the naive forms are FALSE:
`F ≠ G ∪ (G+1/2)` and `G ⊄ F` (checked); the map is doubling, not translation. The flatness
mechanism inside each component is the balanced-slope pairing `d ↔ 13−d` (term slopes
`∓d(13−d)`, symmetric), which is why components stay flat between kinks.

## Remaining (routine, for the canon write-up — mac-mini's item if wanted)

A chamber-by-chamber structural proof of `S2(2y) = 6/7 ⟺ y ∈ G ∪ (isolated)`: both sides
are explicit finite rational interval unions, so the identity is already certification-grade;
the per-chamber derivation (which binding pair `d, 13−d` carries each component) is a clean
one-pager. The `1/7`-widest-free-component of THM-878 (W0) doubles from the corridor's
widest gap analogously.

-> THM-878, THM-853, THM-580, HYP-6975 (kps cont.28), mac-mini-S113 backlog;
`04-computation/lrc14_flat_vs_corridor_6617_deathstar_S17.py` (+ .out).
