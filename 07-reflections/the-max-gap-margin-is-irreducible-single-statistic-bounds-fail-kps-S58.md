---
source: kind-pasteur-2026-07-07-S58
status: honest delineation — attempting inf_E E[maxgap] > 1/7 (the S57 reduced target); the
  clean single-statistic bounds fall short, so the margin is irreducibly in the max (a
  discrepancy statement with the same AP-minimality core as opus's μ_{1/7} route)
tags:
  - lonely-runner
  - LRC14
  - density-floor
  - max-gap
  - three-distance
  - discrepancy
---

# The max-gap margin is irreducible: the single-statistic lower bounds fall short

Owner: *work the former* — attempt `inf_E E[maxgap] > 1/7` (the reduced density-floor target
from S57, `μ_{1/7} ≥ (7/6)(E[maxgap] − 1/7)`) via discrepancy theory. I chased the two clean
rigorous lower bounds; both fall below `1/7`, which pins *why* the target resists an easy proof
and where its content really is.

## A clean identity (reflection symmetry)

The map `x → 1−x` sends `frac(e_i x) → 1 − frac(e_i x)`, so the phase configuration is symmetric
under `θ → 1−θ`. Hence `min_i frac(e_i x)` and `1 − max_i frac(e_i x)` are equidistributed, and
the **origin gap** `gap_0 = min_i frac(e_i x) + (1 − max_i frac(e_i x))` satisfies

> **`E[gap_0] = 2·E[min_i frac(e_i x)]`** (exact; verified to 4 digits).

So `E[maxgap] ≥ E[gap_0] = 2·E[min]`, reducing the origin-gap piece to the **expected minimum
fractional part**.

## Both single-statistic bounds fall below 1/7 — the max is irreducible

Two natural rigorous lower bounds on `E[maxgap]`, and their adversarial infima (over
maxgap-minimizing 13-families):

| bound | value at inf | vs `1/7 = 0.1429` |
|---|---|---|
| **`E[maxgap]`** (target) | **0.211** | **> 1/7, margin +0.069** ✓ |
| `E[gap_0] = 2·E[min]` (origin gap) | 0.137 | **< 1/7** ✗ |
| `E[Σgᵢ²]` (length-biased, `≥ 1/k`) | 0.135 | **< 1/7** ✗ |

`inf_E E[min] = 0.065 < 1/(k+1) = 0.0714`, so `E[gap_0]` genuinely drops below `1/7`; and the
length-biased mean does too. **Yet `E[maxgap]` stays at 0.21** — because when the origin gap (or
the length-biased mean) is small, *other* gaps are large. So the comfortable margin (~0.07) lives
**irreducibly in the max over all gaps**; it cannot be captured by any single gap or by the
length-biased mean. This rules out the two tempting one-statistic reductions and says the target
is a genuine **order-statistic / discrepancy** statement.

## The crux structure: same AP-minimality core as opus's μ_{1/7}, cleaner functional

The `E[maxgap]`-minimizers are **AP-like**: `inf E[maxgap] ≈ 0.211`, exactly `E[maxgap]({1..13})`.
So, empirically, the **AP minimizes `E[maxgap]`** — the same shape as opus-S131's `μ_{1/7}`
AP-minimality (A′). Thus `inf E[maxgap] > 1/7` factors as

1. **AP-minimality of `E[maxgap]`** — `E[maxgap](E) ≥ E[maxgap](AP)` for all `E` (the hard
   extremal core, shared with the `μ_{1/7}` route);
2. **`E[maxgap](AP) > 1/7`** — and for the **AP orbit** `{frac(i·α)}` this is a *classical
   three-distance* quantity: `E[maxgap](AP_13) = E_α[max gap of {iα : i≤13}] = 0.211`, computable
   in closed form from the continued-fraction gap structure (the max gap is `δ + δ'`, consecutive
   convergent contributions). `0.211 > 1/7` with comfortable margin.

So `E[maxgap]` does **not** dodge the AP-minimality difficulty (part 1) — but it is a **cleaner
functional** than `μ_{1/7}`: a single order statistic, a comfortable `+0.07` margin, and an
AP value (part 2) that is a *classical three-distance average*, not the bespoke `477/1078`
piecewise-linear computation. If AP-minimality is going to be proved for *some* gap functional,
`E[maxgap]` — monotone under equidistribution, three-distance-computable at the AP — looks like
the friendlier one to attack. This dovetails with the governing three-gap frame (opus-S132).

## Status

- **Rigorous & clean:** `μ_{1/7} ≥ (7/6)(E[maxgap] − 1/7)` (S57); `E[gap_0] = 2·E[min]` (S58).
- **Ruled out (dead ends):** origin-gap and length-biased single-statistic bounds — both dip
  below `1/7`; the margin is not in any single gap.
- **The target, pinned:** `inf_E E[maxgap] > 1/7 = ` AP-minimality of `E[maxgap]` (hard, shared
  with opus's route) `+ E[maxgap](AP) = ` three-distance average `≈ 0.211 > 1/7` (computable).
- **Does NOT prove LRC(14).** A delineation of the reduced crux, not its resolution.

## Ledger

- **Files:** verification inline (see SESSION-LOG); builds on `lrc_maxgap_reduction_kps_S57.py`.
  HYP-4757.
- **Pointers:** kps-S57 (max-gap reduction, HYP-4747); opus-S131 (`μ_{1/7}` AP-minimality / PZ);
  opus-S132 (three-gap governing frame); three-distance theorem (max gap of `{iα}`).
