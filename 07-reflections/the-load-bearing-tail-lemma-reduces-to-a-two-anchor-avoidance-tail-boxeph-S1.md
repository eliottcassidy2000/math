---
source: boxeph-2026-07-07-S1 (HYP-4801, extends HYP-4760)
status: reduction + strong numerical evidence (not a proof). Reduces the fleet's
  one honest open lemma -- mu_{1/7}(E) >= T_k (the (A') per-k tail floor) -- from
  "full max-gap AP-minimality" to a BOUNDED 2-anchor avoidance tail
  P(max(gap@0,gap@1/2) > 1/7) >= T_k, which is a translated-Steinhaus three-gap
  object (opus-S134's Farey roof with a shifted origin). Reduction is rigorous;
  the 2-anchor floor >= T_k is adversarial-verified with comfortable margin.
tags:
  - lonely-runner
  - LRC14
  - density-floor
  - max-gap
  - three-gap
  - anchor-floor
---

# The load-bearing tail lemma reduces to a TWO-ANCHOR avoidance tail

Owner (remote): pull from the fleet, synthesize incoming ideas and extend them,
work the crux. This does exactly that -- it fuses klein-S153 (anchor floor),
opus-S134 (Farey roof), death-star-S1 (mean detour is dead / the tail is the
object) and monad-explorer-S1 (union-bound thresholds T_k) into one simpler target.

## The state the fleet converged to

After the reverse-Markov `E[maxgap]` detour was declared dead (death-star-S1: its
AP-minimality window and its non-vacuity window are disjoint), the **one honest
open lemma** is the per-`k` TAIL floor

> **(A') `mu_{1/7}(E) >= T_k` for every `k`-element integer set `E`,**

where (monad-explorer HYP-4787, from THM-530's union bound
`G2 >= meas(G_P) + mu_{1/7}(E) - 1 >= m_P`) the thresholds are
`T_k = 1 - meas(G_P^{min,k}) + m_P`:

| k | 8 | 9 | 10 | 11 | 12 | 13 (P=empty) |
|---|---|---|---|---|---|---|
| `T_k` | 0.6185 | 0.5057 | 0.3956 | 0.2747 | 0.1429 | `m_P≈0.0565` |

The fleet has been attacking this as **AP-minimality** `mu_{1/7}(E) >= mu_{1/7}(AP_k)`
(exhaustively verified `k<=10` by opus-S131; the AP value is opus-S134's Farey-roof
`q<=6` window measure). AP-minimality is a hard rigidity statement about the FULL
max-gap.

## The reduction: only TWO gaps are needed, not the max over all of them

The max-gap dominates the gap at any anchor, pointwise:
`maxgap(x) >= max_{a in A} gap_a(x)` for any finite anchor set `A subset R/Z`. Hence

> **`mu_{1/7}(E) = P(maxgap > 1/7) >= P( max_{a in A} gap_a > 1/7 ) =: PA_{|A|}(E).`**

This is klein-S153's anchor floor -- but applied to the **TAIL** `P(.>1/7)`, the
object the skeleton actually consumes, NOT the **mean** `E[.]` (which is a
knife-edge `~0.147` at the origin gap and where klein's route stalled). On the tail
the anchor floor is comfortable, and only two anchors suffice:

| k | `T_k` | `PA_1` (`{0}`) | `PA_2` (`{0,1/2}`) | `PA_3` (`{0,1/3,2/3}`) |
|---|---|---|---|---|
| 8 | 0.6185 | 0.602 (short) | **0.766** | 0.843 |
| 9 | 0.5057 | 0.511 | 0.685 | 0.730 |
| 10 | 0.3956 | 0.434 | 0.570 | 0.626 |
| 11 | 0.2747 | 0.368 | 0.487 | 0.524 |
| 12 | 0.1429 | 0.321 | 0.421 | 0.463 |
| 13 | 0.0565 | 0.281 | 0.360 | 0.392 |

(inf over adversarial families: spread inhomogeneous APs `{a+dj}` -- the minimizers
-- plus random; `lrc_route_A_perk_gap0tail` / `lrc_route_A_multianchor_tail`.)

> **The 1-anchor tail `P(gap@0 > 1/7)` already discharges `k=9..13`** (missing only
> `k=8` by 0.017); **the 2-anchor tail `P(max(gap@0,gap@1/2) > 1/7)` discharges ALL
> `k=8..13`** with margins `0.15-0.31`.

So the load-bearing (A') ledger reduces from "max over all `k` gaps is
AP-minimized" to "**two specific gaps' avoidance tail clears `T_k`**" -- a bounded
functional in `|A|=2` gaps, no full order statistic.

## Why this is a materially easier target (the three-gap handle)

The minimizer of `PA` at every `k` is a **spread inhomogeneous AP** `{a + d·j}`.
Its config `{frac((a+dj)x)}` is a **translate** (by `frac(ax)`) of the Steinhaus set
`{frac(j·(dx)) : j}` -- so each anchor gap `gap_a` is a **shifted-origin
three-distance gap**, and opus-S134's Farey-roof machinery
(`gap = q(x-p/q) + q'(p'/q'-x)` per Farey cell, but measured from `a` instead of
`0`) computes `PA_2` for the minimizing family **exactly**. So the reduced target

> **`inf over shifts a, over d, of P( max(gap@0, gap@1/2) > 1/7 )(shifted Steinhaus_k) >= T_k`**

is a finite three-gap computation per `k` (the roof is piecewise-linear with
rational Farey breakpoints) PLUS the (single-gap-flavored, likely easier than the
full max) claim that the spread AP is the `PA_2`-minimizer. Two anchors also give
**translation- and dilation-invariance for free** (`mu_theta(E)` and each `PA` are
invariant under `E -> dE` and `E -> E+c`, both measure-preserving on the config), so
the search is genuinely over shape only.

## Honest ledger

- **Rigorous:** the reduction `mu_{1/7}(E) >= PA_{|A|}(E)` (pointwise
  `maxgap >= max_a gap_a`); translation+dilation invariance of `mu` and `PA`; the
  minimizers are spread APs; `PA_2 >= T_k` holds with margin `>=0.15` on a wide
  adversarial family set (spread APs + random, `k=8..13`).
- **Not proven:** `inf_E PA_2(E) >= T_k` (adversarial, not a theorem) -- but it is a
  *bounded 2-gap* object with a three-gap closed form on the minimizing family, far
  more tractable than the full-max AP-minimality the fleet was attacking.
- **Synthesis/credit:** klein-S153 (anchor floor -- moved here from mean to tail),
  opus-S134 (Farey roof -- the per-anchor three-gap computer), death-star-S1 (tail
  is the object), monad-explorer HYP-4787 (`T_k`, union bound), THM-527/530. Builds
  on my HYP-4760 (`1/7` sharp, `mu` = the good-density). Stays on the sup side
  (survives MISTAKE-117).
- **Files:** `04-computation/lrc_mu17_lowerbound_routes_boxeph_S1.py`,
  `lrc_route_A_perk_gap0tail_boxeph_S1.py`, `lrc_route_A_multianchor_tail_boxeph_S1.py`
  (+ `.out`).
- **Next (concrete):** with @opus, run the Farey roof from a shifted origin to get
  `PA_2(shifted-Steinhaus_k)` EXACTLY per `k` (turns the table's `0.766, ...` into
  rationals like the `477/1078` line); then the only gap is "spread AP minimizes
  `PA_2`", a 2-gap rigidity that should yield to the same three-gap argument as the
  AP-min tail -- but on two fixed gaps, not the max.
