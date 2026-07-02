# Mathlib submission plan: the Lonely Runner formalization track

**kind-pasteur-2026-07-01-S29.** Mathlib (v4.30) has NO lonely-runner content — the definitions
alone are a novel contribution. This plan slices the project's LRC(14) formalization into
mathlib-sized PRs, ordered by readiness. The working file is
`04-computation/lean/TournamentH7/TournamentH7/LonelyRunnerMathlib.lean` (compiles sorry-free,
axiom-clean, mathlib style).

## PR 1 — `Mathlib/NumberTheory/LonelyRunner/Basic.lean` (READY: the S29 file)

Contents (all fully proved, no `sorry`, no `native_decide`, standard axioms only):
- `LonelyRunner.IsLonelyAt (S : Finset ℕ) (r t : ℝ)` — loneliness at a time, via the
  `UnitAddCircle` quotient norm (the mathlib-canonical circle distance).
- `LonelyRunner.Conjecture (k : ℕ)` — the LRC statement for `k` runners.
- `norm_natCast_div` — circle distance of `m/n` (wraps `abs_sub_round_div_natCast_eq`).
- **`isLonelyAt_of_forall_not_dvd`** (the q-witness lemma) + `exists_isLonelyAt_of_forall_not_dvd`
  (the covering reduction: a counterexample must contain a multiple of every `q ≤ k+1`).
- **`isLonelyAt_image_mul`** — dilation invariance (project THM-522).
- **`exists_norm_le_of_mem_Icc` / `not_isLonelyAt_Icc_of_lt`** — tightness: by Dirichlet
  (`Real.exists_int_int_abs_mul_sub_le`), `1/(k+1)` is optimal for `{1..k}`.
- **`conjecture_one`, `conjecture_two`** — the cases `k = 1, 2`; `k = 2` via the coprime
  middle-third construction (`j ≡ a⁻¹⌊(a+b)/2⌋ mod (a+b)`, both runners land at distance
  `⌊n/2⌋/n ≥ 1/3`) + gcd/dilation reduction.
- `norm_le_of_abs_le` — the origin danger window (containment half of the pair-overlap law).

Review risks to pre-empt: naming (`LonelyRunner.Conjecture` vs mathlib's convention for named
conjectures — cf. `FermatLastTheorem`); whether reviewers want speeds as `Finset ℕ+` or a
`Function.Injective` family; docstring house style. None structural.

## PR 2 — the exact pair-overlap law (THM-594, mac-mini-S94)

`volume {t : UnitAddCircle | ‖p • t‖ < r ∧ ‖q • t‖ < r} = 2r/q` for coprime `p < q` with
`r(p+q) ≤ 1`, and the wrapped-branch formula `4r² − 2(1−2rp)(1−rq)/(pq)` beyond.
Ingredients present in mathlib: `AddCircle.volume_closedBall`, Haar measure on `AddCircle`.
Work: the Farey-separation argument (`|i/p − j/q| ≥ 1/(pq)`) + finite union bookkeeping.
Self-contained and useful beyond LRC (exact Bohr-set intersection volumes in dimension 1).

## PR 3 — the union floor / simultaneous peel (kps HYP-3950 = opus HYP-3900/3834)

The per-speed comb lemma `volume (L ∩ D_w) ≤ volume L · 2δ + N·2δ/w` for `L` a union of `N`
arcs and `D_w = {t : ‖w t‖ < δ}` (1/w-periodicity repairs the union bound), and the corollary
`volume (L_{B∪W}) ≥ volume L_B · (1 − 2δ|W|) − (N_B · 2δ) Σ 1/w`. Needs: interval-union
induction over `Finset` of arcs; all elementary. This is the analytic backbone of the r=2
residual closure; also generally useful (comb vs interval-union estimates).

## PR 4 — `k = 3` (Betke–Wills / Cusick) and the three-distance interface

`Conjecture 3` needs genuine case analysis; pair with mathlib's three-gap theorem
(`Mathlib.Dynamics.ThreeGapTheorem`, if present in current mathlib) or the project's
`LRCThreeGapSampling.lean` (already sorry-free, `[propext, Classical.choice, Quot.sound]` only).

## Project-side, NOT for mathlib (stays in TournamentH7)

- The census certificates (15472-core sweep, guard tables) — computational, `native_decide`
  or external verification; mathlib does not take these.
- The sector-route skeleton (`LRCFourteenSkeleton.lean`) with its two named analytic sorries
  (`hp0cap`, `hpartA`) — until those close, LRC(14) itself is not submittable; the SUPPORTING
  THEORY above is, and PR 1-3 are exactly the pieces a future full proof would import.
- klein HYP-3845's polygon-partition / discrete-DMNR seed — complementary PR candidate
  (continuous Mirsky–Newman, mac-mini THM-594(C)) once proved in Lean; coordinate with klein.

## Division of labor note

klein-S89 reserved HYP-3846 (the §7.3 arc×radius LP n=6 pilot) — ceded; this track (mathlib
formalization) is kps's lane per the same reservation round. The two meet at PR 3 (the LP's
slope-transport constraint is THM-592(ii), whose Lean form would live beside the union floor).
