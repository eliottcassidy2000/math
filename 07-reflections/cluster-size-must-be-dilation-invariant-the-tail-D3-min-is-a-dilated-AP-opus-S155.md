---
source: opus-2026-07-08-S155
status: working the cluster-monotonicity step, found an EXACT counterexample to klein-S186/S187 /
  LEM-009 / kps-S86 -- the k=11 tail D3-minimum is NOT the block+outlier (0.4587) nor >= D3_10
  (0.4646); the dilated AP + interior point (0,3,6,8,9,12,15,18,21,24,27) has D3 = 0.452986. Root:
  D3 is dilation-invariant, the fixed-window cluster size is not. Corrected axis = longest AP. The
  k=11 CLOSURE survives (tail min ~0.4530 >= bar, margin +0.12); the extremal argument must be
  re-derived. Court case filed.
tags:
  - lrc14
  - covering-floor
  - D3
  - dilation-invariance
  - cluster-monotonicity
  - counterexample
  - extremal
---

# Cluster size must be dilation-invariant: the k=11 tail D3-minimum is a dilated AP, not the block+outlier

**opus-2026-07-08-S155 (HYP-5467).** Owner: work the cluster-monotonicity step. The step (klein-S186/
S187, LEM-009, kps-S86) reduces the k=11 prim-diam ≥ 25 tail to "`D3(E) ≥ D3_{c(E)} ≥ D3_10 = 0.4646`,
`c(E)` = max points in a length-9 window, `D3_c` decreasing; global tail min = block+outlier `{0..9,25}`
= 0.4587." Working it turned up an exact counterexample, and the counterexample points at a clean
general principle.

## The counterexample (exact, by klein's own D3)

> **`A = (0, 3, 6, 8, 9, 12, 15, 18, 21, 24, 27)`** = the AP `3·{0..9}` (common difference 3) plus the
> interior point `8`. Primitive (`gcd = 1`), **prim-diam 27**, and
> **`D3(A) = 0.452986 < D3_10 = 0.4646`** — below the claimed *global* tail lower bound, and below the
> claimed minimizer 0.4587.

Both klein-S184's exact Farey `D3` and opus-S148's independent `moments_exact` return the identical
rational `88747403972619401646021583/195916463945506515076905312`. A 56 840-shape search finds the true
tail min ≈ 0.452986 (at `A` and its reflection `(0,3,6,9,12,15,18,19,21,24,27)`).

## The principle: a dilation-invariant functional needs a dilation-invariant axis

`μ` and its degree-3 floor `D3` are **dilation-invariant**: `W_{cE}(x) = W_E(cx)`, so all moments of `W`
(hence `D3`) are unchanged by `E → cE`; and `prim-diam` is dilation-invariant too. Any extremal or
monotonicity statement for `D3` must therefore be phrased on a **dilation-invariant** axis. klein's
"cluster size" = max points in a fixed **length-9 window** is *not* dilation-invariant, and that is
exactly the crack:

- `A` has window-cluster `c = 5` (predicting `D3 ≥ D3_5 ≈ 0.6`), yet contains a **length-10 AP**
  (`0,3,…,27`). Its dilation-invariant "cluster" is 10.
- `A` is the tail image of the **exhaustive minimizer** `2·{0..9}∪{9} = (0,2,4,6,8,9,10,12,14,16,18)`
  (`D3 = 0.4356`, prim-diam 18). Both are **"AP₁₀ (energy 570) + one point (+20) = R2 = 590"** — the AP
  sitting at scale 2 (compact) or 3 (tail). Same `R2 = 590` as the block+outlier, **different D3** — a
  clean exact witness that `D3` is not a function of `R2`, and the max-`R2` shape is not the min-`D3`
  shape (the fleet knew `R2` "scatters"; here is the sharp pair).

## The fix: longest AP is the clean, monotone, dilation-invariant axis

Replacing the window count by the **longest AP** contained in `E` restores monotonicity: stratified
min `D3` runs `0.76 / 0.67 / 0.61 / 0.62 / 0.62 / 0.54 / 0.52 / 0.467 / 0.453` at longest-AP `= 2..10`,
descending to the extremal **"AP₁₀ + 1 point"** family. The minimum over that family (and the whole
tail search) is ≈ **0.4530**, attained at scale `d = 3` with an interior point — not the `d = 1` block +
far point. klein's block-decorrelation *limits* (`D3_10 = 0.4646`, the `D3_c` table) are correct for
their families; what fails is the identification of those families as the tail *minimizers*.

## What this does to the k=11 leg (honest)

- **The closure is NOT threatened.** Every tail shape searched has `D3 ≥ 0.4530 ≥ bar = 0.3312`
  (margin **+0.12**); the exhaustive prim-diam ≤ 24 (klein-S184) and the `R2 ≤ 590` bound (THM-662)
  both stand.
- **The extremal argument regresses.** "k=11 closes modulo cluster-monotonicity" is not valid as
  stated — the cluster-monotonicity is false. Honest status: **k=11 closes IF the tail infimum ≥ bar
  — strongly evidenced (≈ 0.4530, margin +0.12) — via the corrected longest-AP / AP-extremal picture,
  not window-cluster monotonicity.** A rigorous closure now wants either (a) a dilation-invariant
  extremal lemma ("AP₁₀ + 1 point minimizes `D3` over the tail, and its min ≥ bar"), or (b) the L²
  variance route (Var is near-dominated = local AP overlaps; the longest AP controls the binding near
  part — the dilation-invariant restatement of "cluster controls the floor").

## Ledger

- REFUTED (exact): the block+outlier as tail D3-minimizer; `D3 ≥ D3_10 = 0.4646`; fixed-window
  cluster-monotonicity `D3 ≥ D3_{c(E)}` (LEM-009 / klein-S186/S187 / kps-S86). Court case filed.
- CORRECTED: dilation-invariant axis = longest AP; extremal family = "AP₁₀ + 1 point"; true tail min
  ≈ 0.4530; closure survives (+0.12).
- GENERAL LESSON (MISTAKE-126): state extremal/monotonicity claims for dilation-invariant functionals
  on dilation-invariant axes; test candidate extremizers against their dilates.
- Files: `lrc14_cluster_monotonicity_opus_S155.py`, `lrc14_tail_true_min_opus_S155.py` (+outs);
  re-verified via `lrc14_d3_exact_verify_klein_S184.D3`. Connects to the S154 finding that `Var(W)` is
  the dilation-invariant L² object (`the-discrepancy-route-is-L2-not-L1-...-opus-S154`).
