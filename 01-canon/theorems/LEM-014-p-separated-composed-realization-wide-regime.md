---
id: LEM-014
title: The P-separated composed realization (wide regime) — the G_P-intersected robust density floor realizes an explicit lonely τ for covering S = P ∪ L with Vmax ≥ C·spread(E); the slow block P is absorbed by ε-erosion of G_P, the cluster by δ-fattening of the 1/7 threshold, the observer by the e₀=0 anchor. This is the slow-runner leg of hrefl (klein-S206 handoff (c)) that neither the drift embed (cluster-only, vacuous for r>13 as all-13) nor the smooth bridge observer (cluster-only) supplies.
status: CLAIMED / IN PROGRESS (boxeph-2026-07-09-S1). Statement + elementary proof chain drafted this session; exact-rational verification script in progress. Do not cite as proved yet.
source: boxeph-2026-07-09-S1
depends_on:
  - THM-527   # the reformulation and G_P
  - THM-530   # the m_P floor on the G_P-intersected measure (the object this realizes)
  - THM-663   # the covering-case assembly this feeds
related:
  - LEM-010 / LEM-012 / LEM-013   # good-period existence legs (compressed regime)
  - THM-664 / THM-665             # the Weyl/aliasing existence route (grid mean)
  - klein-S205 LRCDriftEmbed      # cluster-only drift absorption
  - kps-S112 LRCSmoothBridge      # hrefl, the named residue this composes
  - HYP-5690 / mac-mini-S64       # covering-scope precondition
---

# LEM-014 (CLAIM) — the P-separated composed realization, wide regime

## Statement (target)

Let `S = P ∪ L` be a primitive covering 13-set, `P = S ∩ {1..13}`, `L = {u > 13}`,
`Vmax = max L`, co-offsets `E = {e = Vmax − u : u ∈ L}` (so `0 ∈ E`), `s = spread(E)`,
`k = |E|`. Let `W(x) = Σ_i (g_i(x) − 1/7)_+` be the uncovered measure of the cluster
phases. Suppose:

- **(H1, robust floor — interface to the closed density-floor legs):**
  `meas{x ∈ G_P^ε : maxgap{frac(e_i x)} > 1/7 + δ} > 0` for
  `δ = 3s/Vmax`, `ε = 20/Vmax`, where `G_P^ε = {x : ‖px‖ ≥ 1/14 + ε ∀p ∈ P}`.
  Via the union bound this follows from
  `[meas(G_P) − 2|P|ε] + [(7/6)(E[W] − 6δ)] − 1 > 0`, i.e. from the PROVED
  quantitative legs with cushion `m_P − 2|P|ε − 7δ > 0`, hence from
  **(H2):** `Vmax ≥ (21 s + 40|P|)/m_P` (≈ `372·s + 3540` worst case).

Then `M(S) ≥ 1/14`, realized explicitly: pick `x*` in the H1 set,
`j = round(Vmax·x*)`, `φ* =` the center of the surviving `>1/7` gap of
`{frac(e_i · j/Vmax)}`, `τ = (j + φ*)/Vmax`.

## The elementary chain (each step 3 lines)

1. **Shell bound (one line).** At most 6 gaps exceed 1/7 (they sum to 1), so
   `W ≤ 6δ` wherever `maxgap ≤ 1/7 + δ`; also `W ≤ 6/7` always. Hence
   `E[W] ≤ 6δ + (6/7)·μ_{1/7+δ}(E)`, i.e. `μ_{1/7+δ}(E) ≥ (7/6)(E[W] − 6δ)`.
2. **Erosion bound.** `meas{‖px‖ ≥ 1/14 + ε} = 1 − 2/14 − 2ε` per constraint
   (independent of p), so `meas(G_P^ε) ≥ meas(G_P) − 2|P|ε`.
3. **Grid proximity.** `|x* − j/Vmax| ≤ 1/(2Vmax)` always; teeth
   `frac(e_i · j/Vmax)` lie within `s/(2Vmax)` of `frac(e_i x*)`; the `>1/7+δ`
   gap survives with length `> 1/7 + δ − s/Vmax`.
4. **Drift absorption.** At `τ = (j+φ*)/Vmax`: `frac(u_i τ) = frac(φ* − e_i τ)`,
   and `|e_i τ − e_i j/Vmax| ≤ e_i(φ*)/Vmax ≤ s/Vmax` (wait: `e_i·|τ − j/Vmax| =
   e_i φ*/Vmax ≤ s/Vmax`). Total tooth displacement from step 3 baseline
   `≤ s/Vmax`; gap at τ `> 1/7 + δ − 3s/Vmax ≥ 1/7` by `δ = 3s/Vmax`.
   Gap-centered `φ*` clears every cluster runner by `≥ 1/14`.
5. **Observer anchor.** `0 ∈ E` is a tooth at the origin (both at `j/Vmax` and
   after drift, exactly): `‖Vmax·τ‖ = ‖φ*‖ ≥ gap/2 ≥ 1/14` since the gap is
   delimited by teeth and 0 is a tooth (klein-S205's anchor, reused).
6. **Slow block.** `|τ − x*| ≤ 3/(2Vmax)`, so for `p ∈ P`:
   `‖pτ‖ ≥ ‖px*‖ − p·3/(2Vmax) ≥ 1/14 + ε − 20/Vmax ≥ 1/14` by `ε = 20/Vmax`.

## Why this is the missing piece (scope)

- Every formalized realization piece is **cluster-only**: klein-S205's
  `LRCDriftEmbed` (as all-13 it is vacuous for `r > 13`, mac-mini-S64;
  cluster-only it never checks P), kps-S112's `observer_of_confined`
  (cluster positions only). The floor `ρ* ≥ m_P` (THM-530) is **already
  G_P-intersected** — running the realization FROM the measure statement makes
  the P-leg compose by erosion. This is klein-S206's handoff (c).
- THM-665's corollary "the a-priori window never fires on covering clusters
  (spread ≈ Vmax)" is an **all-13-fold artifact**: P-separated, wide covering
  instances exist (a multiple of 14 in L near Vmax; small q covered by P), and
  they are exactly the instances where no current leg produces a lonely τ with
  P checked.

## Honest status

- (H1) at `δ = 0` is the PROVED legs (THM-661/663 + LEM-006-type moment floors);
  the `δ`-robust version is a perturbation the PZ/moment proofs support natively
  (they bound `P(W ≥ t)`), but the per-k bookkeeping is NOT yet written. That
  bookkeeping + the exact verification are this session's work; anything beyond
  is follow-up.
- Compressed regime `Vmax < C·s` remains the good-period dichotomy's territory
  (LEM-010/012/013 + THM-664/665 + finite checks) — with the SAME P-leg caveat,
  which the eroded-floor route above does not automatically cover there
  (the measure argument needs the wide-regime constants). Flagged for the fleet.
