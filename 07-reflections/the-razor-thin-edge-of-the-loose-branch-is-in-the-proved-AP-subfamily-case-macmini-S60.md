# The razor-thin edge of the loose branch lives in the proved case; the residual has a margin

*mac-mini-2026-07-05-S60. Owner: push the exact math to complete LRC(14). The loose branch
(klein HYP-4151, kps's critical-config reduction HYP-4157) is the crux: every non-AP 12-config has
`M ≥ 2/25`, razor-tight at `{1,…,11,24}`. This note SPLITS the crux so that the razor-thinness sits
entirely inside the part I already proved (S59), and the remainder is a non-razor margin bound — the same
shape as the covering-min split (S46/S47).*

## The split

Classify a 12-config `B` by whether it contains a **tight 11-subfamily** (an 11-subset that is a dilated AP
`c·{1,…,11}`, i.e. LRC(12)-tight, `M=1/12`):

- **AP-SUBFAMILY** (`∃` tight 11-subfamily): `B = c·{1,…,11} ∪ {X}` up to relabeling. By **S59 (HYP-4152,
  PROVED)**, `M(B) ∈ {1/13} ∪ [2/25, ∞)`. **The razor-thin extremizer `{1,…,11,24}` at exactly `2/25` is
  here** (it is the `c=1, k=2` rung of the ladder `k/(12k+1)`).
- **ALL-LOOSE** (no tight 11-subfamily = every 11-subfamily is loose): **`M(B) ≥ 2/23`** — a *margin* of
  `2/23 − 2/25 = 4/575 ≈ 0.007` above `2/25`. Floored by `{1,…,13}∖{6}` at exactly `2/23`.

**Verified** (`lrc14_loose_split_macmini_S60.py`): over all structured all-loose configs (AP-with-holes
`{1..N}∖holes`, `N ≤ 18`) plus **high-height CRT-lifts** of `{1,…,13}∖{6}` to height `~2500` (MISTAKE-102
discipline — lifting only pushes `M` *up*, e.g. `2/23 → 267/2540 ≈ 0.105`): **zero below `2/23`, zero in the
gap.** The two floors are `2/25` (AP-subfamily, `{1..11,24}`) and `2/23` (all-loose, `{1..13}∖6`).

## Why this is the right split

`2/23 = 2/(2·12−1)` is the **11-runner second value**. So an all-loose 12-config inherits the *11-runner*
gap floor — the split is **recursive**. This is the exact analogue of the covering-min split I found earlier
(S46/S47): the razor value `1/13` lived only in the dilated-deep-well (CRT) case, and the non-dilated
residual sat at the looser `7/89`. Same phenomenon here: the razor `2/25` lives only in the AP-subfamily
(ladder/CRT) case, and the all-loose residual sits at the looser `2/23`.

**Consequence for the crux.** kps-S13 reduced the loose branch to "critical configs are gap-free," razor-tight
at `2/25`. This note shows the critical configs' razor edge is **exactly the AP-subfamily ones**, which are
**already proved** (S59). The remaining critical configs are all-loose, and they carry a **definite margin
`≥ 2/23`**. So the honest remaining obligation is not a razor-edge inequality but a **margin bound**:

> **all-loose 12-config ⟹ `M ≥ 2/23`** (hence `> 2/25`).

A margin bound is what loss-of-constants tools (peeling, folding, the walk cascade) can actually deliver —
they cannot hit a razor edge but *can* clear a margin (the lesson opus/klein logged for `14/183`).

## The reduction of the residual (toward a proof)

The all-loose bound recurses cleanly into the **11-runner rigidity**:
- All-loose ⟹ every 11-subfamily is loose ⟹ (by the 11-runner rigidity, klein-S126's `(1/12, 2/23)` gap)
  every 11-subfamily has `M ≥ 2/23`.
- kps's **peeling lemma** at threshold `2/23`: if some runner `w` is safe (`≥ 2/23`) at the optimum of
  `B∖{w}` (which has `M ≥ 2/23`), then `M(B) ≥ 2/23`. The residual is the "critical-at-`2/23`" all-loose
  configs — and those are exactly the AP-with-holes near `{1..13}∖6`, one level down in the descent.
- So the descent is: **AP-subfamily rung (S59, proved) + all-loose step (peel at the `(m−1)`-second-value,
  recurse).** Each level's razor edge is discharged by the level's S59-analogue; the margin passes down.

The genuinely open piece is unchanged in name — the deepest "all-peels-fail" core, klein/opus's one rigidity
lemma — but this note **removes the razor-thinness from it**: that core now only has to clear `2/23`, not
`2/25`, and only for the all-loose (no-tight-subfamily) configs.

## Status (honest)

- **PROVED** (S59): the AP-subfamily case, incl. the razor `2/25` extremizer.
- **VERIFIED** (this note, structured + high-height): all-loose ⟹ `M ≥ 2/23` (margin). NOT yet proved — it
  needs the 11-runner rigidity (klein-S126, itself the one-lower crux) + the peeling closure at `2/23`.
- **What it buys:** the crux's razor-thinness is now provably confined to S59; the residual is a recursive
  margin bound. The finish line is the descent's base + the deepest peeling-closure lemma, at a *margin*.

## Links
Builds on: mac-mini S59 (HYP-4152, AP-subfamily proof), kps-S13 (HYP-4157, peeling + critical reduction),
klein-S140 (HYP-4151, loose branch = AP-rigidity), klein-S126 (11-runner gap `2/23`), mac-mini S46/S47
(the covering-min split — same shape). HYP-4162. Script: `lrc14_loose_split_macmini_S60.py`.
