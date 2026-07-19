/-
  TournamentH7.LRCMod19LedgerBridge — wiring the mod-19 antipodal-spread lemma into the
  LRC(14) ledger (boxeph-2026-07-19-S129).

  Connects `LonelyRunner.antipodal_cover` (LRCMod19Spread) to the ledger's loneliness-margin
  framework (`TournamentH7.LRCWitness.margin`, the `L(t) = min_i dist(v_i t, ℤ)` used by the
  uniqueness rung `LRCLadderD1`).

  In the ledger's terms: any family in the mod-19 GAP REGIME — margin `< 2/19` at every rational
  time `b/19` (which holds whenever `M(C) = sSup (margin v '' [0,1]) < 2/19`, and `2/19` spans the
  whole uniqueness gap and beyond, `2/19 > 2/25 > 3/38 > 1/13`) — with no speed divisible by 19,
  has residues covering EVERY antipodal unit-pair of `ℤ/19`.  This registers the mod-19
  antipodal-spread constraint as a proved NECESSARY CONDITION on any family populating the n=12
  uniqueness gap (problem (C)); it is translation-sensitive (opus THM-1185/1220 triage) and lives
  on the same (C)-axis as `LRCLadderD1`'s `2/25` reach bound.
-/
import TournamentH7.LRCMod19Spread
import TournamentH7.LRCWitnessAttainment

open TournamentH7.LRCWitness

namespace LRC14

/-- **Ledger bridge: gap-regime ⟹ antipodal covering (mod 19).**  For `k ≥ 1` speeds `v` with no
speed divisible by 19, if the loneliness margin `min_i dist(v_i t, ℤ)` is `< 2/19` at every rational
time `b/19` (the gap regime — implied by `M(C) < 2/19`), then for every nonzero `u : ZMod 19` some
runner has `v_i ≡ u` or `v_i ≡ -u (mod 19)`: the residues cover every antipodal unit-pair of `ℤ/19`. -/
theorem antipodal_cover_of_margin {k : ℕ} [Nonempty (Fin k)] (v : Fin k → ℤ)
    (hunit : ∀ i, ¬ ((19 : ℤ) ∣ v i))
    (hgap : ∀ b : ℤ, margin v ((b : ℝ) / 19) < 2 / 19)
    (u : ZMod 19) (hu : u ≠ 0) :
    ∃ i, (v i : ZMod 19) = u ∨ (v i : ZMod 19) = -u := by
  refine LonelyRunner.antipodal_cover v hunit (fun b => ?_) u hu
  have hnot : ¬ (∀ i, ∀ m : ℤ, (2 : ℝ) / 19 ≤ |(v i : ℝ) * ((b : ℝ) / 19) - m|) := by
    rw [← le_margin_iff v (2 / 19) ((b : ℝ) / 19)]
    exact not_le.mpr (hgap b)
  simp only [not_forall, not_le] at hnot
  exact hnot

end LRC14

#print axioms LRC14.antipodal_cover_of_margin
