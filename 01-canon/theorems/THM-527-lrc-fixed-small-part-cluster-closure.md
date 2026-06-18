---
id: THM-527
title: "The fixed-small-part cluster closure of LRC(14)'s residual case S3 — for a fixed small part P⊆{1..13} (M(P)>1/14) and a single tight cluster of fixed integer offset-shape, there is an explicit cluster-scale threshold V0* above which a GLOBAL WITNESS exists (M(S)≥1/14), and below which it is a finite exact check; the slow-fast/offset-fit reduction made rigorous"
status: PROVED in-scope (elementary slow-fast / θ-sweep construction, verified exact-rational, 0 violations across patterns; the verifier confirmed holds=true, high confidence). Does NOT close the full residual — V0* is non-uniform across offset-shapes and the construction excludes the AP / coordinated-growth family (no fixed bounded small part). NOT a proof of LRC(14).
source: convergence of mac-mini-2026-06-18-S1 (reserved THM-527 + the "circular gap > 2/7" reformulation), kind-pasteur-2026-06-17-S2/S4 (the slow-fast / offset-fit reduction, validated 40/40), and the kind-pasteur-2026-06-18-S4 finishing workflow (fixed-small-part-equidistribution angle, proved + adversarially verified). Three independent derivations of the same object.
depends_on:
  - THM-526   # arc-width lemma + cluster-collapse Lemma A + the k=2 slice; defines the S3 residual
  - THM-523   # only covering sets matter
related:
  - HYP-2581d  # the OPEN uniform ρ*(Δ,P)≥c0 floor — the coordinated-growth core THM-527 does NOT reach
  - OPEN-Q-108
  - OPEN-Q-109 # mac-mini's reservation for the uniform Weyl/three-distance floor
external: Weyl equidistribution; three-distance (Steinhaus) theorem; Lonely Runner Conjecture (open for 13 runners).
---

# THM-527 — The fixed-small-part cluster closure (slow-fast reduction, rigorized)

**Context.** LRC(14) reduces (THM-523/524/525/526) to `M(S) ≥ 1/14` for covering primitive
13-sets; the open residual is case S3 (≥2 speeds `>13`). This theorem closes the
**fixed-small-part single-tight-cluster** slice of S3 — the slice the slow-fast / offset-fit
reduction governs — with an explicit cluster-scale threshold.

## Setup

`S = P ∪ L`, `|P|+|L|=13`, primitive and covering. **`P ⊆ {1,…,13}` fixed** with `M(P) > 1/14`
(equivalently `meas(G_P) > 0`; PROVED for every proper subset of `{1..13}` — only the full
`{1..13}` is tight). `L = {V0+d_0, …, V0+d_c}` a **single tight cluster** of **fixed** integer
offsets `0=d_0<…<d_c` (`d_c < 12·V0`). `G_P` = level-1/14 safe set of `P`; `G_L` likewise.

## The theorem (PROVED in-scope)

> For each fixed `(P, offsets)` there is an explicit `V0* = ⌈1/(b'−a')⌉` (with `[a',b'] ⊆ G_P`
> the sub-arc from Step 1) such that:
> - **`V0 ≥ V0*` ⟹ a global witness `τ*` exists** (`‖vτ*‖ ≥ 1/14 ∀v∈S`), hence `M(S) ≥ 1/14`;
> - **`V0 < V0*` is a finite exact check.**

**Proof (slow-fast / carry-phase).** Write `θ = V0·τ` (fast phase). For `u = V0+d_i`,
`‖uτ‖ ≥ 1/14 ⟺ θ` lies in a width-`6/7` arc centred at `½ − d_i τ`; so the cluster is jointly
safe at slow-time `τ` iff `θ ∈ Ω(τ) := ⋂_i (arc_i)`, of width `w(τ) = max(0, 6/7 − circ_width({d_iτ}))`
— **positive iff the `c+1` phase-points `{frac(d_iτ)}` leave a circular gap `> 1/7`.**
*Step 1:* take the widest arc `I_P` of `G_P`; pick `τ0∈I_P` maximizing `w`, then rational-bisect a
half-width `h` so that on `[a',b']⊆I_P` a *fixed* common sub-arc `Ω* ⊆ Ω(τ)` of width `g>0` exists
for all `τ∈[a',b']` (the slow centres `½−d_iτ` move little, leaving a stable gap `>1/7`).
*Step 2:* as `τ` ranges over `[a',b']`, `θ=V0τ` sweeps an interval of length `V0(b'−a')`; once
`V0(b'−a') ≥ 1` this covers `ℝ/ℤ`, so some `τ*∈[a',b']` has `V0τ* ∈ Ω*` — there the whole cluster is
safe AND (since `τ*∈G_P`) the small part is safe: `τ*` is a **global witness**. `V0* = ⌈1/(b'−a')⌉`,
`τ*` explicit. ∎

**Verification (exact rational):** `Ω* ⊆ Ω(τ)` 0 violations (center+edge, 400+ pts); exact `τ*`
for `V0 ≥ V0*` gives `min_v‖vτ*‖ ≥ 1/14` (all `≥ 0.096`; e.g. `P={1,2,3}`, `V0=81`, min `0.1124`);
120 rows over 40 offset-patterns, 0 failures; direct exact `M` for `P={1,2,3}`, offsets `{0..9}`:
`5/47, 3/25, 4/33, 6/49, 11/89, 15/121, 25/201, 250/2007`, all `≥ 1/14`. (The choice of `(τ0,h)`
affects the SIZE of `V0*`, not correctness.)

## The pigeonhole corollary (cluster-size split, kind-pasteur-S4)

The `m = |L|` phase-points have max circular gap `≥ 1/m` (pigeonhole), so `w(τ) ≥ 1/m − 1/7 > 0`
**for any offsets whenever `m ≤ 6`** (`1/m > 1/7 ⟺ m < 7`). Hence:
- **`|L| ≤ 6`**: the cluster always admits a valid fast-phase at every `τ∈G_P` — Step 2 yields a
  global witness with no condition on the offset-shape (only `V0 ≥ V0*`). (Via the *stronger*
  via-max criterion, gap `> 2/7`, the automatic range is `|L| ≤ 4`, margin `7/c−1 ≥ 4/3`,
  `c=|L|−1`.)
- **`|L| ≥ 7`**: `1/m ≤ 1/7`, pigeonhole insufficient — needs the realized gap `> 1/7`
  (`ρ*(Δ,P) > 0`), the genuine hard core (HYP-2581d / OPEN-Q-109).

## Honest scope — what THM-527 does NOT close

- **`V0*` is non-uniform** across offset-shapes (it grows with cluster spread), so THM-527 closes
  each *fixed* shape but does not uniformly bound S3; the union over unbounded shapes is not finite.
- **Excludes the AP / coordinated-growth family** `{t,2t,…,12t,V}` (no fixed bounded small part;
  primitive+covering+S3 with `W(S\max)→0`). That family is the asymptotically-tight core where
  `M` floors at `2/23` from above — the residual remains OPEN (HYP-2581d). A uniform `ρ*≥c0>0`
  (three-distance/Weyl), or equivalently (via THM-524 binding pairs) "covering forbids the binding
  denominator `D=14q−r` with small `r`", is the finish line.

**Companion proved sub-cases (this session, → THM-526):** the AP-family theorem
(`{1,…,12,m}` covering ⟺ `182|m`, `M ≥ 2/27`; but `k=1`, case S1), and the ALL-MULT7-LARGE
window-collapse (conditional on `w_max < 14 w_min`). Neither closes the coordinated-growth core.
