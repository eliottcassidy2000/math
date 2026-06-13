---
source: oracle-2026-06-01-S544o
status: reframe + computation (how global spread can guarantee local emptiness: the space/frequency duality)
tags:
  - lonely-runner
  - reframing
  - spread
  - decorrelation
  - space-frequency
  - resonance
  - equidistribution
---

# When Global Spread Guarantees Local Emptiness: the Space/Frequency Reframe

LRC is "local emptiness": a hole of clearance `1/n` at the observer. The instinct
"global spread should force a local hole" runs straight into a paradox (S543): the
**regular polygon — maximal spatial spread — is the tight, hardest case**. So
naively, spread is the *enemy* of emptiness. The resolution is to **choose the right
notion of spread**, and it turns on a space/frequency duality.

## The paradox, made precise (instantaneous spatial spread is the enemy)

At a fixed time, "global spread" = an even, high-entropy configuration = small gaps =
**no hole**. Verified (`lrc_global_spread_local_emptiness_s544.py`): the instantaneous
`corr(S_H(t), observer-hole(t))` is **negative for generic speeds** (`-0.22, -0.28`),
i.e. more spatial spread ⇒ smaller observer hole. **So "global spread" must NOT mean
instantaneous spatial evenness** — that choice makes the implication false.

## The right choice #1: spread = DECORRELATION (the independence main term)

Reframe the danger arcs `B_i = {t : ‖v_i t‖ < 1/n}` (each of measure `2/n`). LRC fails
iff they cover `[0,1)`. If the arcs are **statistically independent** (the genuine
"global spread" = decorrelation), inclusion–exclusion gives a safe (lonely) set of
measure

> **`(1 − 2/n)^{n-1} > 0`  — decorrelation GUARANTEES local emptiness.**

Verified: generic primitive speeds realize lonely *measure* `≈ (1-2/n)^{n-1}` (`0.14,
0.11, 0.14` vs main term `0.13` at `n=5,6,7`) — **spread wins**. The regular polygon
(AP) gives lonely measure **exactly `0`** — the arithmetic resonances cancel the
spread term precisely (tight). So:

> **Decorrelated (spread) danger arcs ⇒ a guaranteed hole; the ONLY obstruction is
> the arithmetic resonances (the inside debt, S529).**

## The resolution: spatial spread ⟺ FREQUENCY concentration (the duality)

Why is the regular polygon both maximally *spatially* spread and the obstruction? By
the space/frequency duality (S536 sectors ↔ S537 characters): the regular polygon has
**commensurate (AP) speeds = maximal FREQUENCY concentration** (all resonances align).
It is spatially spread precisely *because* it is frequency-concentrated. So:

> **The correct "global spread" is FREQUENCY spread = incommensurability = no
> resonances = decorrelation.** Frequency-spread ⇒ the arcs decorrelate ⇒ the positive
> main term ⇒ local emptiness. The regular polygon is frequency-*concentrated*
> (commensurate), which is exactly why it is the tight obstruction.

This dissolves the paradox: spatial spread was the wrong axis; **frequency spread is
the friend of local emptiness.** Equidistribution (Q-linearly independent speeds =
maximal frequency spread) is the extreme case where the orbit is dense and the hole is
hit with room to spare (S521o).

## The right choice #2: spread in TIME + ROTATION (the mechanism)

For the hard (commensurate) case, the mechanism is temporal: the runner orbit is a
**moving geodesic** — it cannot sit at the regular polygon. Verified: the max-gap
`max_t g(t)` reaches `≥ 2/n` for every set, and the **hole's position rotates to cover
all `n` sectors**. So:

> **Temporal spread (a moving orbit must DIP to the bunched/holed extreme) + ROTATION
> (the hole's position equidistributes) ⇒ the hole visits the observer.**

Local emptiness is global *anti*-spread (bunching), and what guarantees it is that the
orbit, being globally spread *in time*, must pass through the bunched extreme, and the
bunch's complementary hole rotates onto the observer. For generic speeds this gives
positive-measure loneliness; for the AP it occurs only at the measure-zero boundary
`t = k/n` (tight).

## The reframed choices (the answer)

| our original choice | reframed so spread ⇒ emptiness |
|---|---|
| "spread" = spatial evenness (snapshot) | **wrong axis** — anti-correlated with the hole |
| "spread" = DECORRELATION of danger arcs | ⇒ safe measure `(1-2/n)^{n-1} > 0` (main term) |
| "spread" = FREQUENCY spread (incommensurability) | the *cause* of decorrelation; resonance-free ⇒ hole |
| "spread" = TEMPORAL (the orbit moves) | ⇒ must dip to the bunched/holed extreme |
| observer FIXED at 0 | ride the ROTATION: the hole's position equidistributes onto 0 |

> **One-line reframe:** *local emptiness is guaranteed by frequency spread
> (decorrelation), not spatial spread; the positive independence term `(1-2/n)^{n-1}`
> is the engine, temporal-spread-plus-rotation is the mechanism, and the sole
> obstruction is arithmetic frequency-concentration — maximal, and exactly
> canceling, at the regular polygon.*

So LRC is precisely the claim that **frequency concentration can never cancel the
spatial-decorrelation hole below set-nonemptiness** — the resonances push the lonely
*measure* to `0` (the regular polygon) but, conjecturally, never empty the lonely
*set*.

## Verdict / next
- The paradox (spatial spread = enemy) is resolved by the space/frequency duality:
  frequency spread (decorrelation) is the friend; the regular polygon is
  frequency-concentrated, hence the obstruction.
- Engine = the positive independence main term `(1-2/n)^{n-1}` (verified realized for
  generic speeds, cancelled to 0 at the AP). Mechanism = temporal spread + hole
  rotation (verified: hole covers all sectors).
- Concrete next: (1) quantify "frequency spread" as the smallest resonance / the
  speed set's additive energy, and bound the lonely measure below by
  `(1-2/n)^{n-1} − (energy term)`; (2) make the hole-rotation equidistribution a
  three-gap/CF statement at the observer; (3) the set-vs-measure gap (measure 0 but
  set nonempty) as the true residual at the regular polygon.

## Artifacts
```
04-computation/lrc_global_spread_local_emptiness_s544.py
05-knowledge/results/lrc_global_spread_local_emptiness_s544.out
```
Related: S543 (regular polygon = max H-entropy = tight), S529 (main term + inside
debt), S536/S537 (space/frequency duality), S521o (equidistribution), S530 (apex =
largest gap), S527 (cascade/independence).
