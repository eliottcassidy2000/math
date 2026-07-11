/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: kind-pasteur (LRC multi-agent project, 2026-07-10-S127)
-/
import Mathlib
import TournamentH7.LRCDetunedD3
import TournamentH7.LRCTwoDetunedClearing
import TournamentH7.LRCDissociatedAssembly

/-!
# Wiring the detuned clearings into `MultiDetunedDispatch` (kps-S127)

opus's `MultiDetunedDispatch` (S209) is the THM-678 citation consumed by the dissociated assembly:
*every family with some `g ≥ 2` at detuning level `nonMultCard v g ∈ {2,3}` is lonely.*  Both levels are
now PROVED in the **generic** counting regime — d=2 by opus (`twoDetunedClearing`, S211) and d=3 by kps
(`threeDetunedClearing`, this session).  This file connects `nonMultCard v g = 2` / `= 3` to the concrete
detuned-coordinate hypotheses and **shrinks the citation to its exceptional residual**.

* `lonely14_of_nonMultCard_three` — the requested wire: `nonMultCard v g = 3` + generic counting
  (`Σ_{g∤vᵢ} badCount (vᵢ) g < g`) ⟹ lonely, via `DetunedD3.lonely14_of_three_detuned'`.
* `lonely14_of_nonMultCard_two` — the d=2 analogue, bridging the sum-count to opus's `(q₁,q₂) ≠ (2,2)`.
* `ExceptionalDetunedDispatch` — the dispatch obligation restricted to the NON-generic (`Σ badCount ≥ g`)
  detuned families; and `multiDetunedDispatch_of_exceptional` : the full `MultiDetunedDispatch` follows from
  it plus the citation.  So the cited THM-678 shrinks from "all detuned `d ∈ {2,3}`" to "only the
  half-harmonic exceptional counts" — the `(2,2)` / `(2,2,·)` residual and finitely many small-`q` triples.

Kernel-pure: no `sorry`, no `native_decide`.  Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

/-- The generic counting condition at `g`: the sum of per-coordinate bad-branch counts over the detuned
(non-multiple) coordinates is `< g`.  When it holds, a single branch clears all detuned phases. -/
def genericCount (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  (Finset.univ.filter (fun i => ¬ g ∣ v i)).sum (fun i => DetunedD3.badCount (v i) g) < g.toNat

/-- **The `d = 3` wire.**  A family with exactly three coordinates not divisible by `g ≥ 2`, satisfying the
generic counting `Σ badCount < g`, is lonely.  Extracts the three detuned coordinates from
`nonMultCard v g = 3` and feeds them to `DetunedD3.lonely14_of_three_detuned'`. -/
theorem lonely14_of_nonMultCard_three (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (h3 : nonMultCard v g = 3) (hgen : genericCount v g) :
    ∃ t : ℝ, Lonely 14 v t := by
  rw [nonMultCard] at h3
  obtain ⟨i₁, i₂, i₃, h12, h13, h23, hfilt⟩ := Finset.card_eq_three.mp h3
  have hmem : ∀ j, j ∈ (Finset.univ.filter (fun i => ¬ g ∣ v i)) ↔ ¬ g ∣ v j := by
    intro j; rw [Finset.mem_filter]
    exact ⟨fun h => h.2, fun h => ⟨Finset.mem_univ j, h⟩⟩
  have hi₁ : ¬ g ∣ v i₁ := (hmem i₁).mp (by rw [hfilt]; simp)
  have hi₂ : ¬ g ∣ v i₂ := (hmem i₂).mp (by rw [hfilt]; simp)
  have hi₃ : ¬ g ∣ v i₃ := (hmem i₃).mp (by rw [hfilt]; simp)
  have hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j := by
    intro j hj1 hj2 hj3
    by_contra hnd
    have hj : j ∈ (Finset.univ.filter (fun i => ¬ g ∣ v i)) := (hmem j).mpr hnd
    rw [hfilt] at hj
    simp only [Finset.mem_insert, Finset.mem_singleton] at hj
    rcases hj with h | h | h
    · exact hj1 h
    · exact hj2 h
    · exact hj3 h
  have hsum3 : DetunedD3.badCount (v i₁) g + DetunedD3.badCount (v i₂) g
      + DetunedD3.badCount (v i₃) g < g.toNat := by
    rw [genericCount, hfilt] at hgen
    rw [Finset.sum_insert (by simp [h12, h13]), Finset.sum_insert (by simp [h23]),
      Finset.sum_singleton] at hgen
    omega
  exact DetunedD3.lonely14_of_three_detuned' cite v hv g hg i₁ i₂ i₃ h12 h13 h23 hdvd hi₁ hi₂ hi₃ hsum3

/-- When `q = g/gcd(δ,g) = 2` the bad count is exactly `gcd(δ,g)`, and `g = 2·gcd` — so two `q = 2`
coordinates saturate the count (`badCount₁ + badCount₂ = g`).  This is why the `(2,2)` case is exceptional:
the generic strict bound `Σ badCount < g` cannot hold. -/
theorem badCount_of_q_two {δ g : ℤ} (hg : 0 < g) (hq : g / (Int.gcd δ g : ℤ) = 2) :
    DetunedD3.badCount δ g = Int.gcd δ g ∧ 2 * Int.gcd δ g = g.toNat := by
  have hdvd : ((Int.gcd δ g : ℤ)) ∣ g := Int.gcd_dvd_right δ g
  have htoNat : (g / (Int.gcd δ g : ℤ)).toNat = 2 := by rw [hq]; rfl
  refine ⟨?_, ?_⟩
  · rw [DetunedD3.badCount, htoNat]; norm_num
  · have h := Int.mul_ediv_cancel' hdvd
    rw [hq] at h
    omega

/-- **The `d = 2` wire.**  A family with exactly two coordinates not divisible by `g ≥ 2`, satisfying the
generic counting `Σ badCount < g`, is lonely — via opus's `DetunedD2.lonely14_of_two_detuned'`.  The generic
sum bound supplies opus's `(q₁,q₂) ≠ (2,2)` hypothesis (`badCount_of_q_two`: `(2,2)` would force the sum to
equal `g`, contradicting `< g`). -/
theorem lonely14_of_nonMultCard_two (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (h2 : nonMultCard v g = 2) (hgen : genericCount v g) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hg0 : (0 : ℤ) < g := by omega
  rw [nonMultCard] at h2
  obtain ⟨i₁, i₂, h12, hfilt⟩ := Finset.card_eq_two.mp h2
  have hmem : ∀ j, j ∈ (Finset.univ.filter (fun i => ¬ g ∣ v i)) ↔ ¬ g ∣ v j := by
    intro j; rw [Finset.mem_filter]
    exact ⟨fun h => h.2, fun h => ⟨Finset.mem_univ j, h⟩⟩
  have hi₁ : ¬ g ∣ v i₁ := (hmem i₁).mp (by rw [hfilt]; simp)
  have hi₂ : ¬ g ∣ v i₂ := (hmem i₂).mp (by rw [hfilt]; simp)
  have hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → g ∣ v j := by
    intro j hj1 hj2
    by_contra hnd
    have hj : j ∈ (Finset.univ.filter (fun i => ¬ g ∣ v i)) := (hmem j).mpr hnd
    rw [hfilt] at hj
    simp only [Finset.mem_insert, Finset.mem_singleton] at hj
    rcases hj with h | h
    · exact hj1 h
    · exact hj2 h
  have hsum2 : DetunedD3.badCount (v i₁) g + DetunedD3.badCount (v i₂) g < g.toNat := by
    rw [genericCount, hfilt, Finset.sum_insert (by simp [h12]), Finset.sum_singleton] at hgen
    exact hgen
  have hq : ¬ (g / (Int.gcd (v i₁) g : ℤ) = 2 ∧ g / (Int.gcd (v i₂) g : ℤ) = 2) := by
    rintro ⟨hq1, hq2⟩
    obtain ⟨hbc1, hgt1⟩ := badCount_of_q_two hg0 hq1
    obtain ⟨hbc2, hgt2⟩ := badCount_of_q_two hg0 hq2
    rw [hbc1, hbc2] at hsum2
    omega
  exact DetunedD2.lonely14_of_two_detuned' cite v hv g hg i₁ i₂ h12 hdvd hi₁ hi₂ hq

/-! ### Shrinking the citation: `MultiDetunedDispatch` reduces to its exceptional residual. -/

/-- **The dispatch restricted to the exceptional (non-generic) detuned counts.**  Only families whose
detuned counting is *not* generic (`¬ Σ badCount < g`, i.e. `Σ badCount ≥ g`) need this — the half-harmonic
residual: `(2,2)` at `d = 2`, and `(2,2,·)` plus finitely many small-`q` triples at `d = 3`
(`lrc14_three_detuned_exceptional_kps_S127`).  The generic bulk is discharged unconditionally. -/
def ExceptionalDetunedDispatch : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∀ g : ℤ, 2 ≤ g →
    (nonMultCard v g = 2 ∨ nonMultCard v g = 3) → ¬ genericCount v g →
    ∃ t : ℝ, Lonely 14 v t

/-- **`MultiDetunedDispatch` from its exceptional residual.**  The full THM-678 dispatch (`d ∈ {2,3}`)
follows from the citation and the exceptional-only dispatch: the generic `d = 2` (opus) and `d = 3` (kps)
cases are proved, so only the non-generic half-harmonic counts remain cited.  This is the sharpest form of
opus's `MultiDetunedDispatch` — the near-dilate peel now costs only the `(2,2)` / `(2,2,·)` residual. -/
theorem multiDetunedDispatch_of_exceptional (cite : LRCUpTo13)
    (hexc : ExceptionalDetunedDispatch) : MultiDetunedDispatch := by
  intro v hv hex
  obtain ⟨g, hg, hcard⟩ := hex
  by_cases hgen : genericCount v g
  · rcases hcard with h2 | h3
    · exact lonely14_of_nonMultCard_two cite v hv g hg h2 hgen
    · exact lonely14_of_nonMultCard_three cite v hv g hg h3 hgen
  · exact hexc v hv g hg hcard hgen

/-- **LRC(14) through the dissociated assembly, with the detuned citation shrunk to its residual.**
Composes `multiDetunedDispatch_of_exceptional` with opus's `lrc14_grand_assembly_dissoc`: the THM-678
citation is now needed only for the exceptional half-harmonic detuned counts. -/
theorem lrc14_grand_assembly_dissoc_exceptional (cite : LRCUpTo13)
    (hexc : ExceptionalDetunedDispatch) (hdissoc : ResidualObligationDissoc) :
    LRC14.LRC14Statement :=
  lrc14_grand_assembly_dissoc cite (multiDetunedDispatch_of_exceptional cite hexc) hdissoc

end LRC14Grand
end LonelyRunner
