/-
  TournamentH7.LRCDeepWellReach  (kind-pasteur-2026-07-11-S127 cont.60)

  THE SINGLE-KILLER BALANCE BOUND IN THE reach/margin API — the covering-min VALUE.

  klein-S119 (LRCDeepWellWitness / LRCDeepWellLonely) formalized the deep-well witness
  as per-runner distance bounds and assembled them into `Lonely 14 deepWell14 (14/183)`
  — i.e. every runner is `≥ 1/14` from the origin at `t* = 14/183`. That is the LRC(14)
  (`≥ 1/n`) conclusion. This file assembles the SAME witness into the project's
  `reach = sSup (margin · '' [0,1])` API at the TIGHT value:

        reach(deepWell) ≥ 14/183,   deepWell = {1,…,12, 182},

  keeping the full margin `14/183` (not just `1/14`). This is the covering-min VALUE in
  the sSup framework the covering-min lower-bound crux is stated in — the reach-API
  companion to klein's `Lonely` statement, and the LOWER half of the equality
  `reach(deepWell) = 14/183` (the single-killer balance `[1/13]·182/183 = 14/183` at its
  tight extremal `v_f = 182 = 13·14 = lcm(13,14)`, kps cont.55/59).

  It goes through the `Fin 13` twin of the covering-reach atom
  `LRCCoveringReach.reach_ge_of_covering` (the atom's `rational_point_margin` core is
  polymorphic in the index type, so it copies verbatim from `Fin 12` to `Fin 13`), plus
  a single decidable clearing certificate at the odd base `q = 183` (kps cont.56: the
  covering doorway is odd) — residues `(vᵢ·14) mod 183 = {14,28,…,168,169} ⊆ [14,169]`.

  Kernel-pure (standard trio: propext, Classical.choice, Quot.sound).
-/
import TournamentH7.LRCCoveringReach
import TournamentH7.LRCDecorrelation13

open TournamentH7.LRCWitness LonelyRunner.HarmonicGate TournamentH7.LRCDecorr13

namespace TournamentH7.LRCDeepWellReach

/-- **The covering-reach atom for 13 speeds.**  If a rotation `c` (`0 ≤ c ≤ q`) takes every
speed off the forbidden band `{0,…,±(μ−1)} mod q` — `μ ≤ (vᵢ·c) % q ≤ q − μ` for all `i` —
then `reach v = sSup (margin v '' [0,1]) ≥ μ/q`.  The `Fin 13` twin of
`LRCCoveringReach.reach_ge_of_covering` (whose `rational_point_margin` core is polymorphic). -/
theorem reach_ge_of_covering13 (v : Fin 13 → ℤ) (q c μ : ℤ)
    (hq : 0 < q) (hc0 : 0 ≤ c) (hcq : c ≤ q)
    (hclear : ∀ i, μ ≤ (v i * c) % q ∧ (v i * c) % q ≤ q - μ) :
    (μ : ℝ) / q ≤ sSup (margin v '' Set.Icc (0 : ℝ) 1) := by
  have hmargin : (μ : ℝ) / q ≤ margin v ((c : ℝ) / q) :=
    (le_margin_iff v ((μ : ℝ) / q) ((c : ℝ) / q)).2
      (rational_point_margin v c q μ hq hclear)
  have hqR : (0 : ℝ) < q := by exact_mod_cast hq
  have ht01 : (c : ℝ) / q ∈ Set.Icc (0 : ℝ) 1 := by
    refine Set.mem_Icc.mpr ⟨?_, ?_⟩
    · exact div_nonneg (by exact_mod_cast hc0) (le_of_lt hqR)
    · rw [div_le_one hqR]; exact_mod_cast hcq
  have hmem : margin v ((c : ℝ) / q) ∈ margin v '' Set.Icc (0 : ℝ) 1 := ⟨_, ht01, rfl⟩
  have hbdd : BddAbove (margin v '' Set.Icc (0 : ℝ) 1) :=
    ⟨1 / 2, by rintro _ ⟨t, _, rfl⟩; exact margin_le_half13 v t⟩
  exact le_trans hmargin (le_csSup hbdd hmem)

/-- The **deep well** `{1,…,12, 182}` (`= klein's deepWell14`), the covering-min minimizer
at `n = 14` (`182 = 13·14 = lcm(13,14)`, the minimal single killer repairing `d = 13, 14`). -/
def deepWell : Fin 13 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 182]

/-- **The deep well reaches the covering-min floor:** `reach(deepWell) ≥ 14/183`.
Witness `t* = 14/183` (`q = 183`, rotation `c = 14`, `μ = 14`): the residues
`(vᵢ·14) mod 183 = {14,28,…,168,169}` all lie in the band `[14,169] = [μ, q−μ]`, with
runners `1` and `182` at the two edges (the tight band-fit). This keeps the full `14/183`
margin — the covering-min VALUE — where klein's `Lonely 14` keeps only `1/14`. -/
theorem deepWell_reach : (14 : ℝ) / 183 ≤ sSup (margin deepWell '' Set.Icc (0 : ℝ) 1) := by
  have hcov : ∀ i, (14 : ℤ) ≤ (deepWell i * 14) % 183 ∧ (deepWell i * 14) % 183 ≤ 183 - 14 := by
    decide
  have h := reach_ge_of_covering13 deepWell 183 14 14 (by norm_num) (by norm_num) (by norm_num) hcov
  norm_num at h
  exact h

/-- The deep well strictly clears the tight LRC(14) bound `1/14`: `reach(deepWell) > 1/14`
(margin surplus `14/183 − 1/14 = 13/2562`). The covering minimizer is *loose*. -/
theorem deepWell_reach_gt_tight : (1 : ℝ) / 14 < sSup (margin deepWell '' Set.Icc (0 : ℝ) 1) := by
  have h := deepWell_reach
  have hlt : (1 : ℝ) / 14 < (14 : ℝ) / 183 := by norm_num
  linarith

#print axioms reach_ge_of_covering13
#print axioms deepWell_reach
#print axioms deepWell_reach_gt_tight

end TournamentH7.LRCDeepWellReach
