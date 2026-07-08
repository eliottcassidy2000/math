---
source: opus-2026-07-08-S148
status: exact closed-form degree-3 covering floor (a sharpening within mac-mini's THM-661);
  de-risks the binding k=11 leg 4.6x; block-minimized exhaustively. Concurrent with
  mac-mini's THM-661 (degree-<=4 LP) and klein's THM-660 compact checks — reconciled, cited.
tags:
  - lrc14
  - covering-floor
  - paley-zygmund
  - moment-method
  - closed-form
  - k11
---

# The degree-3 covering floor has a closed form, and it de-risks k=11

**opus-2026-07-08-S148 (HYP-5327).** Owner: keep working the LRC(14) bleeding edge. The
bleeding edge tonight is the covering reformulation (mac-mini THM-657: `μ_{1/7}(E) = P(W>0)`,
`W = Σ(g_i−1/7)_+` the uncovered measure) and the moment floors on it — klein's
Paley–Zygmund `B_2 = E[W]²/E[W²]` (THM-660) and mac-mini's unified degree-`d` bound `B_d`
(THM-661). The six density-floor legs `k=8..13` are all `μ ≥ B_d(E)`; the single **binding**
one is `k=11`, where `B_2` clears the honest bar by a razor-thin `+0.0159`.

## The point

`B_d` is defined as an LP (best degree-`d` polynomial `p ≤ 1_{w>0}` on `[0,6/7]`). For
`d = 3` the LP is unnecessary — the optimum is a **closed form**. The optimal minorant is

> `p(t) = 1 − (1 − t/r)²(1 − t/M)`,  `M = 6/7`.

It is feasible for **any** `r`: `p(0) = 0`, and `p(t) ≤ 1 ⟺ (1−t/r)²(1−t/M) ≥ 0 ⟺ t ≤ M`
(the square is automatic; the second factor is `≥ 0` on `[0,M]`). Its expectation is a
rational function of `r` maximized at `r* = (m2 − m3/M)/(m1 − m2/M)`, giving

> **`B_3(E) = E[W]/M + (E[W] − E[W²]/M)² / (E[W²] − E[W³]/M)`**  (`m_i = E[W^i]`; valid when
> `m2 − m3/M > 0`).

The LP confirmed it to the digit. It's the exact three-moment Markov–Krein bound, written
as one line — the degree-3 sibling of `B_2 = m1²/m2` (whose optimal minorant is the same
form with the `(1−t/M)` factor dropped: `1 − (1−t/r)²`).

## Why it matters: the binding leg

THM-660/THM-661 discharge `k ≥ 11` with `B_2`. At `k=11` that is `+0.0159` — the tightest
point in the entire density-floor program. The exact `B_3` lifts it to **`+0.0735` — 4.6×
thicker** (block value `54912120381817/135668932727076 = 0.404751 ≥ bar 83549/252252`), and
exhaustively (diam ≤ 14) the **block is the exact `B_3`-minimizer**, so `min B_3 = 0.404751`
over the whole compact regime — the same `+0.0735` uniformly. The razor-thin margin that made
`k=11` uncomfortable is gone at essentially no cost: same moments `E[W], E[W²], E[W³]` (mac-mini
already computes them via Farey-cell integration), one more closed-form line.

The closed form also clears all six legs for the block on its own
(`+0.00058/+0.0055/+0.033/+0.0735/+0.159/+0.257` at `k=8..13`), but `k=8,9` are thin at
degree 3 — mac-mini's `B_4` is the right tool there. The clean division: **`B_4` for the
small-`k` legs (thick), `B_3` for the binding `k=11` (the exact sharpening `B_2` needed),
`B_2` fine for `k=12,13`** (already comfortable).

## Honest positioning

This does not close a new leg — `k=11` was already discharged (thinly) by the fleet's `B_2`
compact check (`+0.0156`, klein-S178/S179; me extending the exhaustive to diam ≤ 16) plus the
decorrelation tail (`+0.055`, mac-mini). It is a **robustness upgrade**: the one uncomfortable
margin in the proof becomes comfortable, with an exact closed form rather than an LP. Two small
infrastructure gains came with it: mac-mini's block-only exact PZ integrator is now generalized
to **arbitrary families** (`pz_exact`, validated against the block values to the digit), which
is what let me run the exhaustive `B_3`-minimizer check; and klein's exhaustive PZ compact check
is extended one diameter (≤ 16, min `0.3468`, all clear).

## The lesson worth keeping

When a moment bound is defined by an LP over feasible polynomials, check whether the extremal
polynomial has a **fixed factored shape** first. Here every degree had one:
`1 − (1−t/r)²·Π(boundary factors)` — the double root at the free point `r`, times linear
factors pinning the endpoint behavior (`1−t/M` for degree 3, appearing because the LP's slope
constraint `s = −1/M` is active at `t = M`). The LP is then a one-parameter optimization with a
rational critical point, not a numeric solve — and the bound becomes an exact certificate. The
Markov–Krein principal representation is exactly this factored shape; recognizing it turns
"verified by LP" into "here is the rational."

## Ledger

- EXACT/NEW: the closed-form `B_3 = m1/M + (m1−m2/M)²/(m2−m3/M)`; optimal minorant
  `1−(1−t/r*)²(1−t/M)`; k=11 block/compact-min `B_3 = 0.404751`, margin `+0.0735` (4.6× PZ).
- INFRA: `pz_exact` (arbitrary-family exact covering-moment integrator); exhaustive PZ to
  diam ≤ 16; exhaustive `B_3`-min (block-minimized).
- HONEST: a robustness sharpening of the binding leg, not a newly-closed leg; k=8,9 want `B_4`.
- Files: `lrc14_pz_general_integrator`, `lrc14_pz_degree3_floor`,
  `lrc14_degree3_closed_form_floor`, `lrc14_pz_k11_exhaustive` `_opus_S148` (+outs);
  THM-661 addendum; HYP-5327.
- Builds on / cites: mac-mini THM-657 (covering reformulation) + THM-661 (the `B_d` umbrella),
  klein THM-660 (PZ = `B_2`) + S179 (compact checks, `Var(W)~R2`); Markov–Krein moment problem.
