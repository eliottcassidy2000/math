---
source: kind-pasteur-2026-07-06-S53
status: GREEN formalization (the S52 coarse reduction, kernel-pure) + the single-scale
  density-floor structure it leaves for Route 1's analytic node
tags:
  - lonely-runner
  - LRC14
  - route-1
  - coarse-reduction
  - density-floor
  - lean
  - three-gap
---

# The coarse reduction is formalized (Route 1), and the density node is reduced to single-scale

Owner: *work more on Route 1 density and any other remaining tasks.* Formalization is back
in scope, so this session turns the S52 coarse/scale reduction — the tool that survived
opus-S130's Route-2 retirement because it bounds the **loneliness threshold** (the right
object) — into GREEN Lean, and then reads off what it leaves for the density node.

## What I formalized (GREEN, kernel-pure — `LRCCoarseReduction.lean`, HYP-4707)

Five theorems, all `[propext, Classical.choice, Quot.sound]`, no `sorry`, no `native_decide`,
integrated into the manifest (corpus builds, 8727 jobs):

- **`reach_transfer_coarse`** — the core, one 1-Lipschitz triangle inequality: if every
  coarse phase `k i · s₀` is `≥ μ` from every integer, `|a i| ≤ A`, `|s₀| ≤ 1`, `0 < L`, then
  every fine phase `(a i + L·k i)·(s₀/L)` is `≥ μ − A/L` from every integer. This *is*
  `M(v) ≥ M(K) − A/L` at the level of a fixed witness.
- **`lonely14_of_coarse_of_lonely13`** — with `μ = 1/13` and the budget `A/L ≤ 1/13 − 1/14`,
  the transferred margin clears `1/14`.
- **`lonely13_of_card_le12`** — a 13-tuple whose distinct nonzero values number `≤ 12` is
  `1/13`-lonely, by the `LRCUpTo13` citation (the coarse-family input).
- **`lonely14_of_coarse_le12`** (the headline) — a 13-speed family `v` with a scale
  decomposition `v i = a i + L·k i` whose coarse part `k` has `≤ 12` distinct nonzero values
  (`v` clusters into `≤ 12` groups at scale `L`), `|a i| ≤ A`, `A/L ≤ 1/13 − 1/14`, is
  **Lonely at 1/14** — from the settled LRC(≤13), with **no new analysis**.

So the entire **multi-scale branch** of Route 1 is now discharged in Lean against the settled
cases. This is the honest content of "peeling/compression," re-grounded on the supremum
(where it survives opus-S130) and made a quantitative certificate rather than a hand-wave.

## What it leaves the density node: single-scale families (the residue)

Route 1 wants the good-period density `ρ(v) = Leb{t ∈ [0,1) : min_i ‖v_i t‖ ≥ 1/14}` to
witness `M(v) ≥ 1/14`. The coarse reduction removes the multi-scale families, so the density
node only needs **single-scale (bounded-ratio) families**. I mapped that residue directly on
the LRC(14) object (`lrc_singlescale_density_floor_kps_S53`):

- **The AP `{1,…,13}` is the unique single-scale tight family:** `M = 1/14` exactly, optimal
  witness `t* = 1/14`, and `ρ = 0` — the good set is a single **isolated** point, yet the
  family is **lonely** (equality). So there is **no uniform positive density floor** (any
  "2/7 uniform ρ" claim is refuted by the AP's zero) — the floor is a **dichotomy**, not a
  constant.
- **Every perturbation jumps off the threshold:** bumping one AP speed by `+1, +2, …` sends
  `M` from `1/14` to `1/13`, with `ρ` jumping to `0.024–0.034 > 0`. So there is a **rung
  gap** `(1/14, 1/13)` around the AP — the direct-LRC(14) analogue of the `(G)` gap — and
  the tight locus is the **isolated** AP.
- **Spread single-scale families have positive density and are easily lonely:** random
  bounded-ratio 13-families give `ρ ≈ 0.08–0.13` and `M ≫ 1/14`.

**The clean readout for the density node.** Over single-scale families the floor splits:
- **near-AP / arithmetic corner** — `ρ → 0`, witness is an isolated arithmetic orbit `{k/14}`,
  `M` pinned to a continued-fraction rung. Closed by **rigidity / three-gap** (the metric
  face of AP-uniqueness; OPEN-Q-110's "the problem is arithmetic", mac-mini HYP-4412), **not**
  by measure.
- **spread bulk** — `ρ > 0`. Closed by **decorrelation** (the Riesz/Fourier witness-density
  estimate, Bedert-style).

(Honest caveat: the raw "distinct gap-count `g`" I logged is coarse — at a small-denominator
witness the phases land on a coarse grid, so `g` is small for *loose* families too. `g` is a
clean tight/loose separator only in mac-mini's *near-tight* regime, which I did not reproduce
here; the `ρ`-dichotomy above is the robust finding.)

## The net

Route 1 now reads: **multi-scale ⟹ LRC(≤13) [FORMALIZED GREEN this session] ; single-scale
⟹ {AP: isolated-witness rigidity} + {spread: positive-density decorrelation}**. The
formalizable structural half is done; the remaining analytic core is exactly the single-scale
density floor, and it is genuinely two problems (rigidity + decorrelation), not one uniform
estimate. That is the sharpest honest statement of where LRC(14) stands on the correctly-aimed
route.

## Ledger

- **GREEN this session:** `LRCCoarseReduction.lean` (5 theorems, kernel-pure, in manifest).
- **Verified:** the single-scale density dichotomy (AP unique zero-density tight; perturbations
  jump to `1/13`; spread `ρ > 0`).
- **Does NOT claim:** a proof of LRC(14). The single-scale density floor (rigidity +
  decorrelation) is open.
- **Files:** `LRCCoarseReduction.lean`; `lrc_singlescale_density_floor_kps_S53.py` (+ `.out`).
  HYP-4707.
- **Pointers:** kps-S52 (the coarse bound); opus-S130 (Route 2 retired, Route 1 the target);
  the Route-1 skeleton (`LRCFourteenSkeleton`, `LRCWitnessPartA`); OPEN-Q-110 (the arithmetic
  near-AP corner); mac-mini HYP-4412 (three-gap); Bedert arXiv:2511.16636 (Riesz products).
