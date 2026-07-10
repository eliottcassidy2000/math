/-
# LRCStrictWitnessFloor — the measure floor needs no measure theory

`SafeMeasureFloor` (kps-S127) asks every residual family for a safe set of positive measure.  This file
proves that obligation follows from a **single point with a strict uniform margin** — by nothing but the
reverse triangle inequality:

  `StrictWitness v`  :  `∃ t₀ ∈ (0,1), ∃ ε > 0, ∀ i m, 1/14 + ε ≤ |vᵢ·t₀ − m|`

  `StrictWitness v → ∃ a < b, Icc a b ⊆ safePeriod v → 0 < volume (safePeriod v)`.

Indeed if every phase clears `1/14` at `t₀` by `ε`, then for `|t − t₀| ≤ δ := min(ε/(M+1), t₀/2,
(1−t₀)/2)` (with `M = Σ|vᵢ|`) each phase still clears `1/14`, since
`|vᵢ·t − m| ≥ |vᵢ·t₀ − m| − |vᵢ|·|t − t₀| ≥ 1/14 + ε − M·ε/(M+1) = 1/14 + ε/(M+1)`.
The whole interval `[t₀−δ, t₀+δ]` is therefore safe, and it has positive length.

**Consequence for the endgame.**  Composing with `residualObligation_of_measureFloor`:

  `LRC14Statement  ⟸  LRCUpTo13  +  (every residual family has a strict-margin point)`,  kernel-pure,

so the entire remaining content of LRC(14) is a *pointwise, elementary* statement — no measure theory, no
continuum estimates.  It is NOT proved here: exhibiting a strict-margin point for every residual family
IS the open case of the Lonely Runner Conjecture (equivalently `Mreach v > 1/14` on the residual class;
the covering census measures `min M = 1/12`, cushion `1/84`, but that is evidence, not proof).  What is
proved here is that nothing *beyond* such a point is needed.

kind-pasteur-2026-07-10-S127.
-/
import Mathlib
import TournamentH7.LRCResidualMeasureFloor

namespace LonelyRunner
namespace LRC14Grand

open MeasureTheory

/-- A **strict lonely witness**: a time strictly inside the period at which every runner clears `1/14`
by a uniform positive margin `ε`. -/
def StrictWitness (v : Fin 13 → ℤ) : Prop :=
  ∃ t₀ : ℝ, t₀ ∈ Set.Ioo (0 : ℝ) 1 ∧ ∃ ε : ℝ, 0 < ε ∧
    ∀ i, ∀ m : ℤ, (1 : ℝ) / 14 + ε ≤ |(v i : ℝ) * t₀ - m|

/-- **A strict-margin point yields a safe interval.**  Pure reverse-triangle: shrinking `t` towards `t₀`
by at most `δ = min(ε/(M+1), t₀/2, (1−t₀)/2)` costs each phase at most `M·δ < ε`. -/
theorem safeInterval_of_strictWitness {v : Fin 13 → ℤ} (h : StrictWitness v) :
    ∃ a b : ℝ, a < b ∧ Set.Icc a b ⊆ safePeriod v := by
  obtain ⟨t₀, ⟨ht0, ht1⟩, ε, hε, hmargin⟩ := h
  set M : ℝ := ∑ i, |(v i : ℝ)| with hMdef
  have hM0 : 0 ≤ M := Finset.sum_nonneg (fun i _ => abs_nonneg _)
  have hMle : ∀ j, |(v j : ℝ)| ≤ M :=
    fun j => Finset.single_le_sum (f := fun i => |(v i : ℝ)|)
      (fun i _ => abs_nonneg _) (Finset.mem_univ j)
  have hM1 : (0 : ℝ) < M + 1 := by linarith
  set δ : ℝ := min (ε / (M + 1)) (min (t₀ / 2) ((1 - t₀) / 2)) with hδdef
  have hδ0 : 0 < δ := by
    refine lt_min (by positivity) (lt_min (by linarith) (by linarith))
  have hδ1 : δ ≤ ε / (M + 1) := min_le_left _ _
  have hδ2 : δ ≤ t₀ / 2 := le_trans (min_le_right _ _) (min_le_left _ _)
  have hδ3 : δ ≤ (1 - t₀) / 2 := le_trans (min_le_right _ _) (min_le_right _ _)
  refine ⟨t₀ - δ, t₀ + δ, by linarith, ?_⟩
  intro t ht
  obtain ⟨hta, htb⟩ := Set.mem_Icc.mp ht
  have habs : |t - t₀| ≤ δ := abs_le.mpr ⟨by linarith, by linarith⟩
  refine ⟨Set.mem_Ico.mpr ⟨by linarith, by linarith⟩, ?_⟩
  intro i m
  -- reverse triangle: |v·t − m| ≥ |v·t₀ − m| − |v|·|t − t₀|
  have hrt : |(v i : ℝ) * t₀ - m| - |(v i : ℝ) * t - m| ≤ |((v i : ℝ) * t₀ - m) - ((v i : ℝ) * t - m)| :=
    abs_sub_abs_le_abs_sub _ _
  have hdiff : ((v i : ℝ) * t₀ - m) - ((v i : ℝ) * t - m) = (v i : ℝ) * (t₀ - t) := by ring
  have hprod : |(v i : ℝ) * (t₀ - t)| = |(v i : ℝ)| * |t₀ - t| := abs_mul _ _
  have habs' : |t₀ - t| ≤ δ := by rwa [abs_sub_comm]
  have hbound : |(v i : ℝ)| * |t₀ - t| ≤ M * δ := by
    apply mul_le_mul (hMle i) habs' (abs_nonneg _) hM0
  have hMδ : M * δ ≤ M * (ε / (M + 1)) := by
    exact mul_le_mul_of_nonneg_left hδ1 hM0
  have hlt : M * (ε / (M + 1)) < ε := by
    have h1 : M * (ε / (M + 1)) = M * ε / (M + 1) := by ring
    rw [h1, div_lt_iff₀ hM1]
    nlinarith
  have hmg := hmargin i m
  have hgoal : (1 : ℝ) / 14 ≤ |(v i : ℝ) * t - m| := by
    rw [hdiff, hprod] at hrt
    linarith
  -- `Lonely 14 v t` unfolds to `(1:ℝ)/(14:ℕ) ≤ …`
  simpa using hgoal

/-- **The measure floor from a strict-margin point.** -/
theorem measureFloor_of_strictWitness {v : Fin 13 → ℤ} (h : StrictWitness v) :
    0 < volume (safePeriod v) := by
  obtain ⟨a, b, hab, hsub⟩ := safeInterval_of_strictWitness h
  exact safePeriod_measure_pos_of_Icc_subset hab hsub

/-- **THE REMAINING OBLIGATION, measure-theory-free.**  Every residual family has a strict-margin point.
Equivalently `Mreach v > 1/14` on the residual class — the open case of LRC(14). -/
def StrictWitnessSupply : Prop := ∀ v : Fin 13 → ℤ, IsResidual v → StrictWitness v

/-- Strict-witness supply gives the measure floor. -/
theorem safeMeasureFloor_of_strictWitnessSupply (h : StrictWitnessSupply) : SafeMeasureFloor :=
  fun v hres => measureFloor_of_strictWitness (h v hres)

/-- **LRC(14) from the citation and a strict-margin point per residual family.**  Kernel-pure; the whole
remaining mathematical content is `StrictWitnessSupply`, a pointwise elementary statement. -/
theorem lrc14_of_strictWitnessSupply (cite : LRCUpTo13) (h : StrictWitnessSupply) :
    LRC14.LRC14Statement :=
  lrc14_of_measureFloor cite (safeMeasureFloor_of_strictWitnessSupply h)

end LRC14Grand
end LonelyRunner
