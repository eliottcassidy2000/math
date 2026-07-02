/-
klein-2026-07-02-S110 (HYP-4015) — the LRC(≤13) CITATION NODE + the repeated-entry
reduction (closing kps-S15's census-class flag).

Owner policy (CLAUDE.md): LRC(n ≤ 13) is PROVED (Sungkawichai–Trakulthongchai
arXiv:2604.23906 for n = 11–13; Rosenfeld n=8; Trakulthongchai n=9,10; classical below)
and enters Lean as a NAMED CITATION HYPOTHESIS — not an axiom, not a sorry. This file
creates that node (`LRCUpTo13`) and uses it to close the flagged census-class gap:
MONOTONE tuples with REPEATED entries (which the pack generators never enumerate) have
≤ 12 distinct speeds, hence are lonely at gap 1/13 ≥ 1/14 by the citation. `Lonely`
depends only on the RANGE of the speed vector, so repeats are redundant by construction.
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner

/-- **The LRC(≤13) citation node** (owner policy): for every k ≤ 12 nonzero speeds there
is a time with all runners at circle-distance ≥ 1/(k+1) from the origin. Consumed as a
hypothesis parameter by the endgame; discharged by citation, never by sorry. -/
def LRCUpTo13 : Prop :=
  ∀ (k : ℕ), k ≤ 12 → ∀ (w : Fin k → ℤ), (∀ i, w i ≠ 0) →
    ∃ t : ℝ, Lonely (k + 1) w t

/-- **The repeated-entry reduction**: a 13-tuple with a repeated entry has ≤ 12 distinct
speeds and is lonely at the 1/14 band by the LRC(≤13) citation. Closes the census-class
gap flagged in kps-S15 (pack generators enumerate distinct tuples only). -/
theorem lonely14_of_repeat (cite : LRCUpTo13) (v : Fin 13 → ℤ)
    (hv : ∀ i, v i ≠ 0) {i j : Fin 13} (hij : i ≠ j) (heq : v i = v j) :
    ∃ t : ℝ, Lonely 14 v t := by
  set S : Finset ℤ := Finset.univ.image v with hS
  have hcard : S.card ≤ 12 := by
    have hle : S.card ≤ 13 := by
      simpa [hS] using Finset.card_image_le (s := (Finset.univ : Finset (Fin 13))) (f := v)
    rcases Nat.lt_or_ge S.card 13 with h | h
    · omega
    · exfalso
      have heq13 : S.card = 13 := le_antisymm hle h
      have hinj : Set.InjOn v (Finset.univ : Finset (Fin 13)) := by
        apply Finset.injOn_of_card_image_eq
        simpa [hS] using heq13
      exact hij (hinj (Finset.mem_coe.mpr (Finset.mem_univ i))
        (Finset.mem_coe.mpr (Finset.mem_univ j)) heq)
  -- enumerate the distinct values
  let e : Fin S.card ≃ S := S.equivFin.symm
  set w : Fin S.card → ℤ := fun m => (e m : ℤ) with hw
  have hwne : ∀ m, w m ≠ 0 := by
    intro m
    have hmem : (e m : ℤ) ∈ S := (e m).2
    simp only [hS, Finset.mem_image] at hmem
    obtain ⟨i0, -, hi0⟩ := hmem
    rw [hw]
    simp only
    rw [← hi0]
    exact hv i0
  obtain ⟨t, hL⟩ := cite S.card hcard w hwne
  refine ⟨t, fun i0 m => ?_⟩
  have hmem : v i0 ∈ S := by
    simp [hS, Finset.mem_image]
  set m0 : Fin S.card := e.symm ⟨v i0, hmem⟩ with hm0
  have hwm0 : w m0 = v i0 := by
    rw [hw, hm0]
    simp only
    rw [Equiv.apply_symm_apply]
  have hstep : (1 : ℝ) / ((S.card : ℝ) + 1) ≤ |(v i0 : ℝ) * t - m| := by
    have := hL m0 m
    rw [hwm0] at this
    have hcast : ((S.card + 1 : ℕ) : ℝ) = (S.card : ℝ) + 1 := by push_cast; ring
    rwa [hcast] at this
  have h14 : (1 : ℝ) / 14 ≤ 1 / ((S.card : ℝ) + 1) := by
    have hpos : (0 : ℝ) < (S.card : ℝ) + 1 := by positivity
    have hle14 : ((S.card : ℝ) + 1) ≤ 14 := by
      have h12 : (S.card : ℝ) ≤ 12 := by exact_mod_cast hcard
      linarith
    exact one_div_le_one_div_of_le hpos hle14
  linarith

end LonelyRunner
