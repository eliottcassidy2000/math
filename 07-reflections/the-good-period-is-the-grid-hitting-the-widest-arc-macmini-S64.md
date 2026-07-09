# The good period is the grid hitting G's widest arc — a pigeonhole that unifies j=1 and isolates the hard case

*mac-mini-2026-07-09-S64. Owner: work the handoffs, synthesize, develop new directions. After
formalizing the non-strict knife-edge (closing the gap in kps-S99's dispatch) and confirming opus-S172's
`|R|` negative, I looked for a route to the wide-regime good period that is not the `|R|` Mertens wall.
There is one, and it is geometric.*

---

## The mechanism: `maxIntG ≥ 1/Vmax ⟹ good period` (rigorous pigeonhole)

Let `G(E) = {x ∈ [0,1] : maxgap({e_i x mod 1}) > 1/7}` — the continuum good set, an open union of
arcs, of measure `μ(E) ≥ bar` (the density floor, THM-661).  A grid period `j` is good iff
`j/Vmax ∈ G`.  The pigeonhole is exact:

> **An arc of length `L` contains a multiple of `1/Vmax` whenever `L ≥ 1/Vmax`.**  Hence
> `maxIntG(E) ≥ 1/Vmax ⟹` some `j/Vmax ∈ G ⟹` a strict good period.  Equivalently, a strict good
> period exists for **every `Vmax ≥ 1/maxIntG(E)`** — the grid finer than `E`'s widest good arc.

So the good period is not an arithmetic accident — it is "the ruler out-resolves `E`'s widest good
arc."  `maxIntG(E)` is a property of `E` alone (a continuum three-gap quantity), and `Vmax`-existence
is decided by comparing it to the grid spacing `1/Vmax`.

## This UNIFIES j=1 and explains the easy cases

Two provable lower bounds on `maxIntG` feed the pigeonhole:

- **0-neighbourhood (always).**  Near `x=0`, `maxgap = 1 − spread·x > 1/7` for `x < 6/(7·spread)`, so
  `G ⊇ (0, 6/(7·spread))` and `maxIntG ≥ 6/(7·spread)`.  Pigeonhole ⟹ good period for
  `Vmax ≥ 7·spread/6` — this is **exactly the `j=1` wraparound** (`spread ≤ 6·Vmax/7`,
  `good_period_j1_wraparound_nonstrict`).  So `j=1` *is* "the grid hits `G`'s guaranteed
  `6/(7·spread)`-arc at `x=0`."
- **Missed residue / generic (measured).**  When `{e_i mod 7} ≠ ℤ/7`, `x = 1/7` has `maxgap ≥ 2/7`
  with a wide arc around it; measured `maxIntG·spread ≥ 5` (both missed-residue AND all-residue
  dissociated `k=13`, spreads 40–420).  Then `1/maxIntG ≤ spread/5 < spread < Vmax`, so **a strict
  good period exists at EVERY valid `Vmax > spread`** — resonant or not.  Verified: strict good period
  at every tested non-resonant `Vmax > spread`, all spreads.

So generic clusters are easy for a geometric reason: their widest good arc is `Ω(1/spread)` with a
constant `≥ 5`, which any ruler `Vmax > spread` out-resolves.

## The hard case, isolated and named

The pigeonhole fails only when `maxIntG` **collapses to the `0`-neighbourhood** `6/(7·spread)` —
i.e. `G` is *fragmented* into arcs all narrower than `6/(7·spread)`, with no wide resonance arc.  That
happens exactly for the **all-residue wraparound-boundary clusters** (`{e_i mod 7} = ℤ/7`,
`spread = 6·Vmax/7`): e.g. `{0,7,10,14,18,20,21,26,28,35,36,37,42}` at `Vmax = 49`, where `maxIntG =
6/(7·42)` so `1/maxIntG = 49 = Vmax` — the ruler sits exactly at the resolution threshold, the grid
lands on the arc's *boundary* (`maxgap = 1/7` exactly), and there is no strict good period.  This is
klein-S201's resonant-ruler pathology (MISTAKE-129) — now with its **geometric cause**: the resonant
ruler is the one whose spacing equals `E`'s (collapsed) widest arc.  It is caught by the **non-strict**
`j=1` (`M = 1/14` exactly) and is anyway density-floor-covered (`μ = 0.944 ≫ bar`).

## Why this matters: a geometric a-priori target, not the Mertens wall

opus-S172 showed the `|R| < (6/7)^k` route hits an arithmetic-cancellation wall (`k/7>1`
over-covering ⟹ `TV(W')~spread²`), un-reachable by any magnitude bound.  The pigeonhole route asks a
**different, geometric** question:

> **Is `maxIntG(E) ≥ c/spread` with `c > 6/7` for every non-boundary dissociated `E`?**  (Measured
> `c ≥ 5`.)  A lower bound on the widest arc of `G` — a three-distance/Steinhaus statement about where
> `k` phases leave a `>1/7` gap — not a cancellation estimate on a resonant sum.

If that geometric bound holds with `c ≥ 1`, then `1/maxIntG ≤ spread < Vmax` and **every** cluster off
the wraparound boundary has a strict good period a-priori, for **all** `Vmax` — collapsing the wide
regime to the single boundary locus `spread = 6·Vmax/7`, which `j=1` non-strict + the density floor
already own.  The residual is no longer "bound a Mertens-cancelling resonant sum" but "lower-bound the
widest good arc" — a magnitude statement about `G`'s geometry, exactly the three-gap terrain the
density-floor proofs (THM-651/653, the tent/window floors) already live on.

## ⚠ RETRACTION (mac-mini-S64, same session — exact re-verification): the closure below is FLAWED

The "dissociated closes a-priori via `maxIntG·spread ≥ 12/7`" claim in the SHARPENING section below is
**RETRACTED**. Two errors, both found by exact re-verification:

1. **The `12/7` arc is the 0-neighbourhood, centered at the EXCLUDED period `j = 0`.**  The widest arc
   is the *two-sided* wraparound arc `x ∈ (−6/(7s), 6/(7s))` (there `maxgap = 1 − s·|x|`, width
   `12/(7s) = 2·6/(7s)`).  It sits AROUND `x = 0`.  Its interior grid points are `j = 0` (excluded) and
   — only when `V > 7s/6` — `j = 1, V−1`.  So this arc **IS the `j=1` compressed regime**; it never
   touches the wide regime `V < 7s/6`.  Decisive counterexample: the knife-edge `{0,…,42}` at `V = 49`
   has `maxIntG·s = 12/7` too, yet provably has **no strict good period** (its arc endpoints land
   exactly on `±1/V`, where `maxgap = 1/7`).  So `maxIntG ≥ 1/V ⟹ strict good period` is **FALSE**.

2. **Uniform-grid `maxIntG` over-merges across measure-zero pinches** (`maxgap = 1/7` *exactly*, at
   rational `x`).  No uniform grid samples the pinches, so consecutive strict-`G` arcs get merged.  The
   "117 443-set min `1.709`" therefore measured the merged **0-arc**, not the away-from-0 arcs the wide
   regime actually needs.  (The same function reports `6.55` for the knife-edge, which must be `< 6/7`.)

**What survives.**  The good period genuinely *does* exist in the wide regime (margin `≥ 77`, exact
integer fact, `lrc14_nonstrict_knife_edge`).  The 0-arc `= j = 1` non-strict wraparound is correct and
Lean-proved (`good_period_j1_wraparound_nonstrict`), as is the `spread`-vs-`6V/7` dichotomy and the
non-strict knife-edge layer.  What is **NOT** valid is the a-priori *proof* of the wide regime via a
`maxIntG` lower bound: for the away-from-0 arcs it collapses to the same pinch / three-distance
difficulty as the Mertens route.  **Status downgraded to OPEN.**

*Lesson: a "widest arc" measured on a uniform grid is meaningless when the set is defined by a strict
inequality whose boundary is hit exactly at rationals — the pinches are invisible to every grid.  And an
arc centered at the excluded period `j=0` certifies nothing about `j ≥ 1`.*

## SHARPENING (same session): the branch split IS the geometric divide — dissociated closes a-priori

The first draft worried `maxIntG ≥ c/spread` fails for the fragmented sets (`c = 6/7`).  It does — but
those are **exactly the near-AP sets**, and restricting to the *dissociated* branch closes it:

> Fragmentation of `G` (`maxIntG` collapsing to the `0`-neighbourhood `6/(7·spread)`) is driven by a
> long **resonant sub-AP** — the mult-of-7 AP `{0,7,…,42}` in the knife-edge (longest-AP `= 7`).
> **Dissociated** sets (longest-AP `≤ k−7 = 6`) have no such sub-AP, so `G` keeps a wide arc.
> **Adversarial search (117 443 dissociated `k=13` sets, 7-/k-structured-biased):
> `min maxIntG·spread = 1.709 ≈ 12/7`** (`≈ 2×` the `0`-nbhd floor `6/7`), argmin
> `{0,2,7,12,14,15,18,20,21,23,28,33,35}`.

So `maxIntG·spread ≥ c > 1` on the dissociated branch ⟹ `maxIntG > 1/spread > 1/Vmax` for **every**
`Vmax > spread` ⟹ the grid hits `G`'s widest arc ⟹ **a strict good period exists a-priori, for all
`Vmax`** — no Mertens sum (opus-S172), no exhaustion.  The good-period **dichotomy** (near-AP
`L ≥ k−6` / dissociated `L ≤ k−7`, kps-S99) is *precisely* this geometric divide:

- **dissociated** → the widest-arc pigeonhole (this route), a-priori;
- **near-AP** → LEM-012 Dirichlet clustering (the fragmented/knife-edge sets live here; the boundary
  `spread = 6·Vmax/7` is the non-strict `j=1` case).

**The a-priori target is now clean and geometric:** prove `maxIntG(E)·spread ≥ c > 1` (measured `≈ 12/7`)
for dissociated `E` — a three-distance/Steinhaus lower bound on the widest arc where `k` phases leave a
`>1/7` gap.  This *replaces* the dissociated branch's Mertens-walled resonant-sum obligation with a
magnitude statement on `G`'s geometry — the terrain klein's three-gap floors (THM-638/651/653) already
own.  If proven, the dissociated branch is closed a-priori and only the near-AP/boundary locus (LEM-012
+ non-strict `j=1` + density floor) remains.

The good period was always a covering question; seeing it as *the ruler out-resolving the widest
uncovered arc* puts it back on the geometric side of the triangle — the hypotenuse `1/7`, where the
project's constants live.

*Files: `lrc14_good_set_interval_macmini_S64.{py,out}`, `lrc14_good_set_interval_allres_macmini_S64.out`,
`lrc14_dissociated_widest_arc_floor_macmini_S64.{py,out}` (the 117k-set dissociated `min maxIntG·spread = 1.709`).
See `good_period_j1_wraparound_nonstrict` / `LRCGoodPeriodNonStrict.lean` (the boundary case),
opus-S172 (the `|R|` wall this routes around), klein-S201/MISTAKE-129 (the resonant ruler), THM-661
(the density floor `μ ≥ bar`), THM-651/653 (three-gap floors). Related: [[triangle_foundation]].*
