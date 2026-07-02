/-
klein-2026-07-02-S102 (HYP-4007 K1) — the seven-translate uniqueness, membership level.

The ℚ-avatar of module 3's core fact (opus INSIGHT 1): for 7 ∤ P, the seven 1/7-translates
of the P-comb partition the circle at the level of residues — for every x, the value
P·x − (P·j)/7 lands within the comb's window for EXACTLY ONE j < 7, because
j ↦ P·j mod 7 is a bijection. This file proves the arithmetic core (K1-arith): the
residue selector exists and is unique. The interval-level K1 and the K2 cursor lemma
(length_inter_partition) build on this; spec in SESSION-LOG klein-S102.
-/
import Mathlib

namespace LRCSevenTranslate

/-- K1-arith: for 7 ∤ P, every residue class mod 7 is hit by P·j for exactly one j < 7. -/
theorem exists_unique_translate (P : ℕ) (hP : ¬ (7 ∣ P)) (t : ZMod 7) :
    ∃! j : ZMod 7, (P : ZMod 7) * j = t := by
  have hne : (P : ZMod 7) ≠ 0 := by
    rw [Ne, CharP.cast_eq_zero_iff (ZMod 7) 7 P]
    exact hP
  haveI : Fact (Nat.Prime 7) := ⟨by norm_num⟩
  refine ⟨(P : ZMod 7)⁻¹ * t, ?_, fun y hy => ?_⟩
  · field_simp
  · have h1 : (P : ZMod 7)⁻¹ * ((P : ZMod 7) * y) = (P : ZMod 7)⁻¹ * t := by rw [hy]
    rwa [← mul_assoc, inv_mul_cancel₀ hne, one_mul] at h1

end LRCSevenTranslate
