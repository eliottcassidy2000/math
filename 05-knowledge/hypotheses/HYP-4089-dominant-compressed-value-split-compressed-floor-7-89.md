# HYP-4089 — The dominant/compressed split is a VALUE split: the razor-thin covering-min lives in the DISCHARGED dominant branch; the open compressed leaf floors at 7/89 (the lcm(12,14) shadow of the deep well)

**Status:** VERIFIED (exact-M, two independent descents + structured scan). Empirical, not yet proved.
**Source:** mac-mini-2026-07-04-S45
**Type:** CONFIRMATION + Lean-dispatch bridge (not a novel result — see credit)

> **CREDIT / OVERLAP (honest).** This substantially re-derives concurrent same-day fleet results:
> **klein-S129 (HYP-4090)** — deep well `{1..12,182}` is the UNIQUE covering family at `14/183`, every
> other is `≥ 7/89 = {1..11,13,84}` (drop-12 residue-liar, PROVED by kps `LRCResidueLiar`), residual
> non-sharp (`35/16287` slack); **kps-S5 (HYP-4087)** — "the far-runner tail ... the deep well" is in
> the DISCHARGED dominant branch; **kps-S6 / opus-S72 (HYP-4091)** — the open leaf is exactly `hcomp`
> (compressed), with "razor-sharp `14/183` only at the already-proved deep well." **My value-add is
> only:** (a) an *independent* confirmation of the `7/89` floor computed **directly on opus's Lean
> `compressed` predicate** (`largest ≤ 13·second-largest`), not on klein's "minimal-tightener" partition
> — so the two agree on the floor for the exact set kps must discharge; (b) the **lcm-shadow** structural
> "why `7/89`" below.

---

## Setup — the covering dispatch (opus, formalized by kps HYP-4091)
`covering_lonely_of_dominant_or_compressed` splits every primitive covering 13-family `V` by whether it
has a **dominant runner**:
- **DOMINANT:** `∃ i, ∀ j≠i: 13·|v_j| < |v_i|` — the largest runner exceeds `13×` the second-largest.
  **DISCHARGED in Lean** (kps-S5 sharp dominant peel `hdom_closed_abs`, HYP-4087): the giant runner's danger
  arcs (`1/(7v_i) < 1/(91B)`) are too thin to spoil the 12-runner base's good interval (LRC(13)).
- **COMPRESSED:** `largest ≤ 13·(second-largest)` — no dominant runner. **The SOLE remaining open leaf.**

## The statement (confirmed) — the dispatch split is also a VALUE split
> **The razor-thin covering-min `14/183` lives ENTIRELY in the DISCHARGED dominant branch.
> The open COMPRESSED leaf (`hcomp`) is bounded away: every compressed covering family has
> `M ≥ 7/89 = 0.078652`, a margin of `35/16287 ≈ 0.00215` above `14/183 = 0.076503`.**
> (klein-S129 established this for "every non-deep-well family"; here it is verified on opus's exact
> `compressed` predicate — the two sets share the floor `7/89`.)

- **Deep well `{1,…,12, 182}`** (the covering-min extremizer, `M = 14/183`) is **DOMINANT**:
  `182 > 13·12 = 156`. So the extremizer is in the **closed** branch. This is not a coincidence —
  it is *forced*: to cover the hard pair `{13,14}` with a **single** killer needs `lcm(13,14)=182 | killer`,
  and with base second-max `≤ 12` that killer `≥ 182 > 156` is always dominant. **A single-killer covering
  of `{13,14}` cannot be compressed.**
- **Compressed extremizer `{1,…,11, 13, 84}`** (`M = 7/89 = [0;12,7]`) is the **lcm(12,14) shadow**:
  unable to use `182`, the compressed family covers `{12,14}` together with `84 = lcm(12,14) = 12·7`
  (`84 ≤ 13·13`, compressed), leaving `13` as a separate small killer over the tight base `{1,…,11}`.
  This is the exact structural parallel of the deep well, one lcm down.

## The parallel (why 7/89 and why 14/183)
| branch | cover `{a,b}` via one killer | killer | base | `M` | CF |
|---|---|---|---|---|---|
| **dominant** (closed) | `{13,14}` | `182 = lcm(13,14)` | `{1..12}` | `14/183` | `[0;13,14]` |
| **compressed** (open) | `{12,14}` | `84 = lcm(12,14)` | `{1..11,13}` | `7/89` | `[0;12,7]` |

`84 = lcm(12,14)`, and covering `14` forces `7 ∣ K` on the `{1..11,13,12K}` ladder `M = K/(12K+1)`, whose
minimum over covering `K` is `K=7 → 7/89`. Both extremizers are single-"big-lcm-killer" + tight base;
the covering constraint on `{13,14}` (needing `182`) is what makes the tight one dominant.

## Verification (exact-M, `Fraction` breakpoints)
- **`0` compressed covering families below `14/183`** over `>15,000` sampled/structured (random broad sample
  `0/7939`; structured drop-1/drop-2 ladders `0/…`).
- **Compressed floor `= 7/89`**, reached at `{1,…,11,13,84}`, robust across: broad random + local descent
  (two independent runs, wide range `[1,260]`, single- and multi-restart) both converged to `7/89`;
  structured exhaustive over `{1..13}\{≤2 elts} + killers`. Values seen climbing the ladder above it:
  `7/89 < 2/23 < 9/109 < 1/12 < …` (Ostrowski rungs `[0;12,k]`, `[0;11,k]`).
- `7/89 > 1/13 = 0.076923` (so the floor even clears the 12-runner LRC bound) `> 14/183`.

## Consequence for the Lean endgame (the point)
kps's remaining obligation `hcomp` (**every compressed covering family is lonely**, i.e. `M ≥ 1/14`) is
**not razor-thin**: the true value in the compressed branch is `≥ 7/89`, a `~2.8%` margin above the
covering-min and a `~10%` margin above the LRC target `1/14`. So the compressed leaf can be discharged with
a **slack** bound (any `c ≤ 7/89`, e.g. a clean `M ≥ 1/13` or even `M ≥ 1/14 + ε`), **not** the sharp
`14/183` equioscillation. The sharp `14/183` argument (THM-618 killer-offset) is only needed in the
dominant branch — which is already closed. **This relocates all remaining razor-thinness out of the open
leaf.**

## Resonance
`89 = 12·7+1` and killer `84 = 12·7` coincide with **kps's `residueLiar84_lonely` (K=7, 37/89)** — the
compressed floor and kps's residue-liar family are the same `12·7` arithmetic. The compressed branch is the
`[0;12,k]`-ladder world (klein HYP-4080 spectral gap, Ostrowski); the dominant branch is the `[0;13,k]`
world (THM-618, deep well). See [[the-covering-min-is-an-ostrowski-ladder-and-the-ap-and-deep-well-are-its-ends]].

## Not claimed / open
The floor `7/89` is **empirical** (exact-M over a large structured + random family, not a proof). A proof
that *compressed covering ⟹ M ≥ 1/13* (or `≥ 7/89`) would fully de-sharpen kps's open leaf — likely a
"the extra `{13,14}`-coverage slot costs a runner, dropping to the 12-runner LRC floor `1/13`" peel,
complementary to opus's dominant peel. Related: THM-617, THM-618, HYP-4082, HYP-4091.

**Scripts:** `04-computation/compressed_covering_min_macmini_20260704.py`,
`05-knowledge/results/compressed_covering_min_macmini_20260704.out`.
