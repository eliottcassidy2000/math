/-
# LRCResidualMeasureFloor — LRC(14) from the LRC(≤13) citation plus ONE measure floor

The sharpest statement of the LRC(14) formalization's state.

`Lonely 14 v t` unfolds to `∀ i m, 1/14 ≤ |v i · t − m|` — i.e. `t` lies in the **safe set** of `v`
(every phase inside the closed band `[1/14, 13/14]` mod 1).  So the grand assembly's `ResidualObligation`
says precisely: *every residual family has a nonempty safe set*.  By klein's THM-685 the modulus side is
closed by elementary counting, and **the remaining analytic content of the covering case is exactly the
measure floors**.  This file makes that literal in Lean:

  `SafeMeasureFloor`   :  every residual family has `volume {t ∈ [0,1) | Lonely 14 v t} > 0`;
  `SafeIntervalFloor`  :  … contains a nondegenerate closed interval  (⟹ the measure floor),

and then

  **`LRC14Statement  ⟸  LRCUpTo13 citation  +  SafeMeasureFloor`,  kernel-pure.**

Every other branch of `lrc14_grand_assembly` is already discharged with foundational axioms only
(non-covering sieve; ratio ≤ 13 spread13; the window ≤ 22 six-witness pigeonhole LEM-024 — the
`native_decide` there was removed in kps-S127; dominant peel; repeated speeds; detuned; common-residue).
Hence the ENTIRE remaining mathematical content of LRC(14) is the single floor above.

This is a *reduction*, not a proof of the floor: the floor is the open analytic core (mac-mini's
witness-floor bricks, klein's THM-685 Corollary 1, boxeph's `μ_L` chain measures).  The interval form is
the shape mac-mini's `LRCIntervalBridge` (`Ico_subset_safeSet_of_bounds`) already produces.  Stating both
as named Lean hypotheses, with `SafeIntervalFloor → SafeMeasureFloor → ResidualObligation → LRC(14)`, is
what puts the formalization in its best available state.

kind-pasteur-2026-07-10-S127.
-/
import Mathlib
import TournamentH7.LRC14GrandAssembly

namespace LonelyRunner
namespace LRC14Grand

open MeasureTheory

/-- The **safe set** of `v` inside one period: the times at which every runner clears `1/14`.
`t ∈ safePeriod v` is definitionally `t ∈ [0,1)` together with `Lonely 14 v t`. -/
def safePeriod (v : Fin 13 → ℤ) : Set ℝ := {t | t ∈ Set.Ico (0 : ℝ) 1 ∧ Lonely 14 v t}

/-- **The residual class**, bundled: exactly the hypotheses of `ResidualObligation` — covering,
scale-gapped, compressed (no dominant runner), distinct absolute speeds, reaching past the
machine-checked window, no detuned harmonic, not a coarse near-AP, difference-primitive. -/
def IsResidual (v : Fin 13 → ℤ) : Prop :=
  (∀ i, v i ≠ 0) ∧
  LRC14.CoveringFamily v ∧
  GapFamily v ∧
  (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) ∧
  (∀ i j, i ≠ j → |v i| ≠ |v j|) ∧
  (∃ i, 23 ≤ |v i|) ∧
  (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) ∧
  (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
    (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
    (Finset.univ.image k).card ≤ 12) ∧
  (∀ d : ℤ, 2 ≤ d → ∀ a : ℤ, (∀ i, d ∣ (v i - a)) → d ∣ a)

/-- **THE SOLE REMAINING ANALYTIC OBLIGATION (measure form).**  Every residual family has a safe set of
positive measure inside one period.  By klein's THM-685 this is exactly the "measure floor" the covering
case reduces to. -/
def SafeMeasureFloor : Prop := ∀ v : Fin 13 → ℤ, IsResidual v → 0 < volume (safePeriod v)

/-- **The same obligation, interval form** — the shape the analytic machinery actually produces
(mac-mini's `LRCIntervalBridge.Ico_subset_safeSet_of_bounds`, boxeph's `μ_L` interval certificates). -/
def SafeIntervalFloor : Prop :=
  ∀ v : Fin 13 → ℤ, IsResidual v → ∃ a b : ℝ, a < b ∧ Set.Icc a b ⊆ safePeriod v

/-- **A positive-measure safe set yields a lonely time.**  Positive measure ⟹ nonempty, and every point
of the safe set *is* a witness (`Lonely` is exactly safe-set membership). -/
theorem lonely_of_safePeriod_measure_pos {v : Fin 13 → ℤ} (h : 0 < volume (safePeriod v)) :
    ∃ t : ℝ, Lonely 14 v t := by
  obtain ⟨t, ht⟩ := nonempty_of_measure_ne_zero (ne_of_gt h)
  exact ⟨t, ht.2⟩

/-- **A nondegenerate safe interval gives the measure floor.** -/
theorem safePeriod_measure_pos_of_Icc_subset {v : Fin 13 → ℤ} {a b : ℝ}
    (hab : a < b) (hsub : Set.Icc a b ⊆ safePeriod v) : 0 < volume (safePeriod v) := by
  have hmono : volume (Set.Icc a b) ≤ volume (safePeriod v) := measure_mono hsub
  have hIcc : volume (Set.Icc a b) = ENNReal.ofReal (b - a) := Real.volume_Icc
  have hpos : (0 : ENNReal) < ENNReal.ofReal (b - a) := by
    rw [ENNReal.ofReal_pos]; linarith
  exact lt_of_lt_of_le (hIcc ▸ hpos) hmono

/-- **Interval floor ⟹ measure floor.** -/
theorem measureFloor_of_intervalFloor (h : SafeIntervalFloor) : SafeMeasureFloor := by
  intro v hres
  obtain ⟨a, b, hab, hsub⟩ := h v hres
  exact safePeriod_measure_pos_of_Icc_subset hab hsub

/-- **The residual obligation reduces to the measure floor.** -/
theorem residualObligation_of_measureFloor (hfloor : SafeMeasureFloor) : ResidualObligation := by
  intro v hv hcov hgap hcomp hdist hlarge hdiv hcoarse hres
  exact lonely_of_safePeriod_measure_pos
    (hfloor v ⟨hv, hcov, hgap, hcomp, hdist, hlarge, hdiv, hcoarse, hres⟩)

/-- **LRC(14) from the citation and ONE measure floor.**  With `lrc14_grand_assembly` now
foundational-axioms-only (the `winData22` `native_decide` was removed by the LEM-024 six-witness
pigeonhole), this is the complete statement of the formalization: every branch is discharged, and the
entire remaining mathematical content of LRC(14) is `SafeMeasureFloor`. -/
theorem lrc14_of_measureFloor (cite : LRCUpTo13) (hfloor : SafeMeasureFloor) :
    LRC14.LRC14Statement :=
  lrc14_grand_assembly cite (residualObligation_of_measureFloor hfloor)

/-- **LRC(14) from the citation and a safe interval per residual family** — the form directly consumable
by the interval-certificate stack. -/
theorem lrc14_of_intervalFloor (cite : LRCUpTo13) (hfloor : SafeIntervalFloor) :
    LRC14.LRC14Statement :=
  lrc14_of_measureFloor cite (measureFloor_of_intervalFloor hfloor)

end LRC14Grand
end LonelyRunner
