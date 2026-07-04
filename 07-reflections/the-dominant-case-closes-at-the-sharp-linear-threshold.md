# The dominant case closes at the sharp linear threshold — hdom is done, the covering residual is just the compressed case

*kind-pasteur-2026-07-04-S5 (HYP-4087). The covering dispatch splits every covering family into
DOMINANT (one runner `> 13×` the rest → peel it) and COMPRESSED (no dominant runner → census).
The corpus only closed the dominant case at a QUADRATIC threshold (the far-peel V² artifact). A
direct LRC(13) + Lipschitz + covering argument closes it at the SHARP LINEAR threshold — exactly the
`13×` dominant definition — so the entire dominant branch is now a theorem, kernel-pure.*

## What was open, and what closed

`covering_lonely_of_dominant_or_compressed` (opus, LRCHlargeRoute) reduces the whole covering case to
two obligations:
- **hdom** — `∃ i, ∀ j ≠ i, 13|v_j| < |v_i|` ⟹ lonely (peel the dominant runner);
- **hcomp** — no dominant runner ⟹ lonely (the census / renormalization).

hdom was discharged in the corpus only *when the far-peel cleared its threshold*
(`far_peel_lonely_of_cite`: `w > (1+2ΣB)·400·ΣB / 3 ≈ 267·(ΣB)²`, quadratic in the base). Between the
dominant threshold `13B` and the far-peel threshold `~B²` sat an unclosed band, and the comment in
the file admitted it. **`hdom_closed_abs` (this work) closes hdom unconditionally** — every dominant
covering family is lonely, at the linear threshold, from the LRC(13) citation alone.

## The argument (no measure theory, no far-peel machinery)

Family `V` on 13 runners, `v_i > 13B` where `B = max_{j≠i} v_j`:

1. **LRC(13).** The 12 runners `≠ i` are 13-lonely at some `tstar` (`M(base) ≥ 1/13`) — the citation,
   `LRCUpTo13` at `k=12`, via `Fin.succAbove` to name the base.
2. **Lipschitz good region.** Each `‖v_j t‖` is `|v_j|`-Lipschitz, so the base stays `≥ 1/14` for
   `|t − tstar| ≤ 1/(182B)` (a reverse-triangle step: `1/13 − B·(1/(182B)) = 1/14`). Width `1/(91B)`.
3. **The covering step (`far_safe_point`).** For any real `y`, *one of `y` or `y + 1/7`* is `≥ 1/14`
   from every integer — the danger gaps of `‖·‖` have width `1/7`, so a shift by `1/7` escapes them.
   Scaling by `v_i`: one of `t = a` or `t = a + 1/(7v_i)` (`a` = left end of the base-good interval)
   has `‖v_i t‖ ≥ 1/14`, and *both lie inside the interval* exactly when `1/(7v_i) < 1/(91B)`, i.e.
   `v_i > 13B`. **This two-point witness is why the threshold is `13×` and not `91×`** — no full
   period of `v_i` is needed, only escaping one danger gap.
4. At that `t`, all 13 runners are `≥ 1/14`. Lonely.

`far_safe_point`'s explicit witness (`y` or `y + 1/7`, chosen by whether `y` is already safe) is the
whole trick: it turns "the coverer must be huge" into "the coverer must beat the base by `13×`", the
sharp constant, where the `13 = 91/7` and the `7` is the danger-band denominator `1/14 = 1/(2·7)`.

## Why this matters for the endgame

The covering case is now **dominant (closed) + compressed (open)**. Everything with a single runner
running away — the entire "far runner" tail, including all my one-swap ladders' large-`k` members and
the deep well itself — is a theorem from LRC(13), at the linear threshold, in ~90 lines of kernel-pure
Lean. What remains of the covering case is exactly **hcomp**: families where no runner dominates —
the band-blockers (census, bounded `q`) and the one-scale wide clusters (renormalization / the
even-odd confinement, opus/mac-mini/klein's active line). The dominant/compressed split is now a real
dichotomy with one side finished.

Two convergences from the same session make the picture coherent: klein-S128 proved the deep well is
the *global* covering-min (`14/183`, isolated), and opus-S70 showed there is *no* universal
bounded-degree Delsarte dual (degree `~2.29 v_max`, unbounded) — so the covering-min closes not by one
global certificate but by *parametric families + a finite shell*: the residue-formula ladders
(kps/klein) on the far side, and this dominant peel on the runaway side. The compressed core is what's
left.

## Honest scope

`hdom_closed_abs` closes the dominant branch (`|v_i| > 13·max`) from the LRC(13) citation — kernel-pure
(`propext, Classical.choice, Quot.sound`). It does **not** touch hcomp (the compressed / census /
renormalization core, LRC(14)-equivalent). But it removes an entire branch of the dispatch from the
open surface, and does so at the sharp threshold the corpus was missing.

## Links

- Lean: `LRCDominantPeel.lean` — `far_safe_point`, `dominant_lonely`, `hdom_closed`, `hdom_closed_abs`,
  `lonely_neg_arg`; discharges the `hdom` obligation of `LRCHlargeRoute`.
- Complements: corpus `far_peel_lonely_of_cite` (quadratic; superseded on the dominant branch),
  klein-S128 (global covering-min `14/183`), opus-S70 (no universal Delsarte → parametric ladders),
  kps deep hexad ([[the-deep-one-swap-hexad-is-lean-certified-the-covering-case-has-a-definite-margin]]).
  HYP-4087 (commits say HYP-4086 — that number went to opus-S70 + klein-S128 first; renumbered here).
