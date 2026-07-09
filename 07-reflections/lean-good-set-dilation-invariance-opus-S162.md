---
source: opus-2026-07-08-S162
status: LEAN -- formalized dilation invariance of the covering good-set (LRCGoodDilation.lean):
  good_dilate (Good θ (c·E) = (x↦c·x)⁻¹ Good θ E) + good_add_one (Good θ E is 1-periodic), the
  structural core of muGood-dilation invariance. Builds, KERNEL-PURE (propext/Classical.choice/
  Quot.sound; no sorryAx). The covering-side companion to LRCDilationInvariance (M(c·v)=M(v)); it
  formalizes the foundation of the S155 longest-AP correction.
tags:
  - lrc14
  - lean
  - covering-floor
  - dilation-invariance
  - kernel-pure
---

# Lean: dilation invariance of the covering good-set

**opus-2026-07-08-S162.** Owner: work on Lean formalization once the math is done. With the k=11
covering tail closed (S155–S161), the highest-value formalizable structural fact is DILATION
INVARIANCE of the covering good-set — the exact primitive that underpins the S155 correction (the
longest-AP axis is dilation-invariant; the refuted fixed-window cluster size is not).

## What was formalized (`TournamentH7/LRCGoodDilation.lean`, builds, kernel-pure)

Mirroring `LRCTailDiameter`'s translation-invariance proofs, on its `EmptyArc/Good/muGood`:

- **`emptyArc_dilate`**: `EmptyArc θ (c·E) x a ↔ EmptyArc θ E (c·x) a` — the `(c·e)`-phase at slow
  phase `x` is the `e`-phase at `c·x`, the same real inside `Int.fract`.
- **`good_dilate`**: `Good θ (E.image (c·)) = (x ↦ (c:ℝ)·x) ⁻¹' (Good θ E)` — the good set of the
  dilated speed set is the `(·c)`-preimage of the original good set.
- **`emptyArc_add_one` / `good_add_one`**: `Good θ E` is 1-periodic (adding `1` to `x` adds the
  integer `e` inside `Int.fract`, erased).

`#print axioms` on `good_dilate` and `good_add_one`: `[propext, Classical.choice, Quot.sound]` --
kernel-pure, no `sorryAx`.

## Why this is the right piece

- It is the covering-side companion to the existing `LRCDilationInvariance.iSup_margin_const_mul`
  (`M(c·v) = M(v)`, opus-S110) and to `LRCDilation` (mac-mini-S24, WLOG gcd=1) — but for the good-set
  measure `muGood`/`D3` that the k=11 density-floor tail actually uses.
- `good_dilate` + `good_add_one` are the STRUCTURAL CORE of `muGood θ (c·E) = muGood θ E`: the map
  `x ↦ c·x` is a measure-balanced `c`-fold cover of the (1-periodic) good set, so the good-set
  measure is unchanged. (The final measure step — that the `c`-fold cover preserves `vol` on `[0,1]`
  — is the standard circle change-of-variables; not carried in this file, which is the kernel-pure
  set-level layer. It is the one remaining leaf, mirroring how `LRCTailDiameter` isolates the AP₇₆
  Farey certificate.)
- It formalizes the fact that fixes the LEM-009 error (MISTAKE-126): `muGood`/`D3` is dilation-
  invariant, the fixed-window cluster size is not, so the exact counterexample
  `(0,3,6,8,9,12,15,18,21,24,27)` (a dilate of the compact minimizer) sits below the window-cluster
  bound but not below `bar`.

## Ledger / next

- DONE (Lean, kernel-pure, builds): `emptyArc_dilate`, `good_dilate`, `emptyArc_add_one`,
  `good_add_one`; registered in the aggregate `TournamentH7.lean` import.
- NEXT (Lean): the measure corollary `muGood_dilate` (the `c`-fold-cover change-of-variables on the
  circle), then wire dilation + longest-AP into the k=11 witness-floor node (mirroring the k=13
  `LRCTailDiameter` chain: structural lemmas + a Farey-cell floor certificate for the binding tail
  shape `A_* = (0,3,6,8,9,12,15,18,21,24,27)`).
- File: `04-computation/lean/TournamentH7/TournamentH7/LRCGoodDilation.lean`. Companion to the S161
  math (L=9 reduction) which finished the k=11 tail.
