# The LRC(14) crux is a lower bound, and its sharp core is a mod-25 covering fact

*kps-2026-07-06-S41 — working toward closing LRC(14), math first. Pinning the
Freiman residual (opus S120) as a lower bound, reducing its sharp core to a finite
mod-25 covering statement, and formalizing the lower-bound certificate.*

## Where the proof stands

opus S120's architecture audit is the frame. LRC(14) needs two open things:

1. **The math crux** — the *Freiman-stability step*: a first-gap member is an
   `(N−2)`-AP + exactly 2 defects; the ≤2-defect world is swept empty, so the
   residual is **rule out ≥3 defects**.
2. **The Lean bridge** — `(A) → Mreach` is unwired (the gap thread doesn't yet
   reach the Fin-13 top-level).

This session works (1), math first.

## The crux is a lower bound

Since **LRC(13) is proved** (owner directive), every 12-speed family has
`M ≥ 1/13`. A first-gap member has `M ∈ (1/13, 2/25)`. So the Freiman step is
*equivalent* to a clean lower bound:

> **longest-AP-subset(v) ≤ 9 (≥3 defects) ⟹ M(v) ≥ 2/25.**

This is the **hard direction** — a lower bound, won by exhibiting a witness, not by
magnitude (the "lower-bound is the hard direction" theme, opus S48). Tested on
7,785 structured `≥3`-defect families at N=12: **0 land in the gap** — the bound
holds robustly.

## An honest redirect: the harmonic/S30 tool is spacing-1 only

I formalized (S30) that a sorted sequence is an AP iff all its second differences
vanish. Tempting to phrase "2 defects" as "≤2 nonzero second differences." But the
dilated member `{1,5,6,11,16,17}` (sub-AP `{1,6,11,16}`, spacing 5) has **3**
nonzero *sorted* second differences — because its sub-AP is interleaved with the
defects when sorted. So the sorted-second-difference measure detects only
*spacing-1* APs; the Freiman step is genuinely about **subset APs** (any spacing),
which is Freiman's inverse-sumset territory, not the discrete Laplacian. Recorded
so the S30 harmonic tool isn't over-applied here.

## The sharp core is a mod-25 covering fact

`M ≥ 2/25` is *sufficiently* witnessed at `t = c/25`: if some `c ∈ (ℤ/25)*` rotates
every speed off the forbidden band `{0, ±1}` mod 25 (i.e. `2 ≤ vᵢc mod 25 ≤ 23`),
then `M ≥ 2/25`. Testing the `≥3`-defect families:

- **~all clear this way** (a rotation exists), **except** those containing a
  **multiple of 25** — which sits at residue `0` for *every* `c`, so no rotation
  helps.
- The multiple-of-25 families clear instead at **small denominators** with `M`
  *far* above `2/25` (e.g. `2/11, 2/17, 3/19, 9/46` — all ≫ `2/25`): being
  divisible by 25 makes them "loose" the easy way.

So the residual splits cleanly:

> **near-tight `≥3`-defect families (no multiple of 25) are mod-25-clearable ⟹
> `M ≥ 2/25`** (a *finite* covering-system statement — klein's covering leg),
> **and** the multiple-of-25 families are loose by a small-denominator witness.

This connects the Freiman step to the project's `2/25`-floor / covering machinery
(klein's `≥ 1/12` covering band, `LRCGapLadder`, `gap_tower_step`) instead of a
bespoke inverse-sumset theorem.

## Lean: the lower-bound certificate (GREEN)

`LRCMod25Floor.lean` (kernel-pure, `[propext, Classical.choice, Quot.sound]`):

- `mod25_covering_floor` — `∀i, 2 ≤ vᵢc mod 25 ≤ 23  ⟹  ∀ i m, 2/25 ≤ |vᵢ·(c/25) − m|`.
- `loose_of_mod25_covering` — the existential loose form `∃ t, ∀ i m, 2/25 ≤ |vᵢ t − m|`.

Both are one-line instances of the `rational_point_margin` atom at `s = 25`,
`μ = 2`. So *once* the mod-25 covering fact is established (the residue-pattern
count), the loose conclusion is machine-checkable from a decidable integer
hypothesis — the lower-bound half of (G) is formally in hand.

## What remains (sharpened)

- **Math:** prove the mod-25 covering fact — every `≥3`-defect 12-speed family
  without a multiple of 25 admits a rotation off `{0, ±1}` mod 25. This is finite
  (residue patterns mod 25) and is klein's covering leg; plus the easy
  multiple-of-25 case (small-denominator witness). Together with the swept
  ≤2-defect world, this closes (G).
- **Assembly:** opus's `(A) → Mreach` bridge (the rank/density case-split) — still
  the gate between a closed (C) and a closed LRC(14) in Lean.

## Pointers

- `lrc_freiman_lowerbound_mod25_kps_S41.out` (lower-bound test + mod-25 reduction),
  `LRCMod25Floor.lean` (the certificate, GREEN).
- opus S120 (Freiman step, architecture audit), HYP-4516 (mod-30 mediant gate);
  mac-mini THM-632 / HYP-4602 (2-defect + order-3 empty); klein `LRCGapLadder`
  (`2/25` window machinery); kps S39 (dilated APs), S40 (per-order gauntlet), S30
  (harmonic — spacing-1 only, redirected here).
