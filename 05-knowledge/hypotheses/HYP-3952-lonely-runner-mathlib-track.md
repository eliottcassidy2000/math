# HYP-3952: The Lonely Runner mathlib track — Basic.lean compiles sorry-free and axiom-clean

**Status:** VERIFIED (machine-checked) — kind-pasteur-2026-07-01-S29. kps block (3950+).

## Claim
The project's LRC work now has a MATHLIB-SUBMITTABLE core:
`04-computation/lean/TournamentH7/TournamentH7/LonelyRunnerMathlib.lean` — compiles sorry-free
through `lake build`, axioms exactly `[propext, Classical.choice, Quot.sound]` (the mathlib
standard three), written to mathlib conventions over `UnitAddCircle`. Mathlib v4.30 has NO
lonely-runner content, so the definitions alone are novel.

## Contents (all fully proved)
1. `IsLonelyAt (S : Finset ℕ) (r t : ℝ)` and `Conjecture (k : ℕ)` — the LRC objects on the
   mathlib-canonical circle (`UnitAddCircle` quotient norm).
2. **q-witness lemma** `isLonelyAt_of_forall_not_dvd` (+ covering reduction
   `exists_isLonelyAt_of_forall_not_dvd`): no multiple of `q` in `S` ⟹ time `1/q` is
   `1/q`-lonely — THM-523's easy half; kills every non-covering configuration.
3. **Dilation invariance** `isLonelyAt_image_mul` — THM-522.
4. **Tightness** `exists_norm_le_of_mem_Icc` / `not_isLonelyAt_Icc_of_lt`: via mathlib's
   Dirichlet (`Real.exists_int_int_abs_mul_sub_le`), some speed in `{1..k}` is always within
   `1/(k+1)` — the conjectured constant is optimal.
5. **`conjecture_one` and `conjecture_two`**: k=2 by the coprime middle-third construction
   (`j ≡ a⁻¹·⌊(a+b)/2⌋ mod (a+b)` puts BOTH runners at distance `⌊n/2⌋/n ≥ 1/3`; the two
   circle-distances agree because `a ≡ −b mod (a+b)`), then gcd/dilation reduction. ZMod-based;
   the key mathlib pieces: `ZMod.coe_mul_inv_eq_one`, `abs_sub_round_div_natCast_eq`.
6. `norm_le_of_abs_le` — origin danger-window containment (THM-594's easy half).

## The submission plan
`03-artifacts/drafts/mathlib-submission-plan-lonely-runner-kps.md`: PR 1 = this file (ready);
PR 2 = exact pair-overlap law (THM-594; `AddCircle.volume_closedBall` exists); PR 3 = the union
floor / simultaneous peel (HYP-3950 = opus HYP-3900) as comb-vs-arc-union volume estimates;
PR 4 = k=3 + three-gap interface. Census certificates and the `hp0cap`/`hpartA` skeleton stay
project-side.

## Coordination
klein-S89 reserved HYP-3846 (the §7.3 arc×radius LP n=6 pilot) minutes before I started it —
CEDED to klein; this mathlib track is the differentiated lane. klein HYP-3845's
polygon-partition/DMNR mathlib seed is complementary (targets THM-594(C), a different lemma
family). mac-mini THM-594's correction to my HYP-3950 pair-law gloss is integrated in the
HYP-3950 file.

## Artifacts
- 04-computation/lean/TournamentH7/TournamentH7/LonelyRunnerMathlib.lean (registered in root)
- 03-artifacts/drafts/mathlib-submission-plan-lonely-runner-kps.md

## Depends on / relates to
THM-522, THM-523, THM-594, HYP-3950, klein HYP-3845/3846, LRCThreeGapSampling.lean,
LRCFourteenSkeleton.lean (hp0cap/hpartA), mathlib Mathlib.NumberTheory.DiophantineApproximation.
