/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Coherent blocker two-cycle alignment (THM-1262)

This dependency-free module consumes the named paper-topology inputs:
ascent-target containment in the low-safe component, danger/safe
disjointness, overlap of consecutive selected teeth, and the binary
mismatch-to-adjacency law.  It checks their ordered-word composition,
existence of a corridor neighbor, exclusion of both marked owners at that
neighbor, and the gcd/lcm seam arithmetic.
-/

namespace LRC14
namespace BlockerTwoCycleAlignment

/-- Two positions are consecutive in the chronological tooth word. -/
def Consecutive (a b : ℕ) : Prop := a + 1 = b ∨ b + 1 = a

/-- `k` lies strictly between `a` and `b`, in either orientation. -/
def StrictBetween (a k b : ℕ) : Prop :=
  (a < k ∧ k < b) ∨ (b < k ∧ k < a)

/-- The reverse high-to-low edge of an oriented two-cycle lies in THM-1248's
binary target range for every positive carrier. -/
theorem reverse_edge_binary_gate
    {carrier low high : ℤ} (hcarrier : 0 < carrier) (hlowHigh : low < high) :
    low < high ∧ high < carrier + high ∧ low < carrier + high := by
  constructor
  · exact hlowHigh
  · constructor <;> linarith

theorem consecutive_ne {a b : ℕ} (h : Consecutive a b) : a ≠ b := by
  unfold Consecutive at h
  omega

/-- Ascent containment and reverse-target danger make the two marked target
teeth disjoint. -/
theorem ascent_protection_disjoint
    {point : Type*} {highTooth lowTooth lowSafe : point → Prop}
    (hHighSafe : ∀ x, highTooth x → lowSafe x)
    (hLowDanger : ∀ x, lowTooth x → ¬ lowSafe x) :
    ¬ ∃ x, highTooth x ∧ lowTooth x := by
  rintro ⟨x, hxHigh, hxLow⟩
  exact hLowDanger x hxLow (hHighSafe x hxHigh)

/-- Since consecutive members of the irredundant word overlap, the protected
high target and reverse low target cannot be consecutive. -/
theorem protected_target_teeth_not_consecutive
    {point : Type*} {highTooth lowTooth lowSafe : point → Prop}
    {highPos lowPos : ℕ}
    (hHighSafe : ∀ x, highTooth x → lowSafe x)
    (hLowDanger : ∀ x, lowTooth x → ¬ lowSafe x)
    (hConsecutiveOverlap : Consecutive highPos lowPos →
      ∃ x, highTooth x ∧ lowTooth x) :
    ¬ Consecutive highPos lowPos := by
  intro hconsecutive
  exact ascent_protection_disjoint hHighSafe hLowDanger
    (hConsecutiveOverlap hconsecutive)

/-- THM-1256's binary landing dichotomy: after nonconsecutivity eliminates
the mismatch/adjacent branch, alignment is forced. -/
theorem nonconsecutive_binary_pair_aligned
    {aligned mismatch : Prop} {highPos lowPos : ℕ}
    (hDichotomy : aligned ∨ mismatch)
    (hMismatchAdjacent : mismatch → Consecutive highPos lowPos)
    (hNonconsecutive : ¬ Consecutive highPos lowPos) :
    aligned := by
  rcases hDichotomy with hAligned | hMismatch
  · exact hAligned
  · exact False.elim (hNonconsecutive (hMismatchAdjacent hMismatch))

/-- Direct composition of ascent protection with the binary reverse-edge
landing law. -/
theorem coherent_two_cycle_forces_alignment
    {point : Type*} {highTooth lowTooth lowSafe : point → Prop}
    {highPos lowPos : ℕ} {aligned mismatch : Prop}
    (hHighSafe : ∀ x, highTooth x → lowSafe x)
    (hLowDanger : ∀ x, lowTooth x → ¬ lowSafe x)
    (hConsecutiveOverlap : Consecutive highPos lowPos →
      ∃ x, highTooth x ∧ lowTooth x)
    (hDichotomy : aligned ∨ mismatch)
    (hMismatchAdjacent : mismatch → Consecutive highPos lowPos) :
    aligned := by
  apply nonconsecutive_binary_pair_aligned hDichotomy hMismatchAdjacent
  exact protected_target_teeth_not_consecutive
    hHighSafe hLowDanger hConsecutiveOverlap

/-- Two distinct nonconsecutive marks have a tooth strictly between them;
it can be chosen as the immediate neighbor of the first mark toward the
second. -/
theorem nonconsecutive_has_toward_neighbor
    {a b : ℕ} (hne : a ≠ b) (hNonconsecutive : ¬ Consecutive a b) :
    ∃ k, StrictBetween a k b ∧ Consecutive a k := by
  by_cases hab : a < b
  · have hgap : a + 1 < b := by
      have hneSucc : a + 1 ≠ b := by
        intro heq
        exact hNonconsecutive (Or.inl heq)
      omega
    refine ⟨a + 1, ?_, ?_⟩
    · exact Or.inl ⟨by omega, hgap⟩
    · exact Or.inl rfl
  · have hba : b < a := by omega
    have hgap : b < a - 1 := by
      have hneSucc : b + 1 ≠ a := by
        intro heq
        exact hNonconsecutive (Or.inr heq)
      omega
    refine ⟨a - 1, ?_, ?_⟩
    · exact Or.inr ⟨hgap, by omega⟩
    · exact Or.inr (by omega)

/-- The two-cycle composition simultaneously forces alignment and a
nonempty chronological corridor. -/
theorem coherent_two_cycle_alignment_with_corridor
    {point : Type*} {highTooth lowTooth lowSafe : point → Prop}
    {highPos lowPos : ℕ} {aligned mismatch : Prop}
    (hPositionsDistinct : highPos ≠ lowPos)
    (hHighSafe : ∀ x, highTooth x → lowSafe x)
    (hLowDanger : ∀ x, lowTooth x → ¬ lowSafe x)
    (hConsecutiveOverlap : Consecutive highPos lowPos →
      ∃ x, highTooth x ∧ lowTooth x)
    (hDichotomy : aligned ∨ mismatch)
    (hMismatchAdjacent : mismatch → Consecutive highPos lowPos) :
    aligned ∧ ∃ k,
      StrictBetween highPos k lowPos ∧ Consecutive highPos k := by
  have hNonconsecutive : ¬ Consecutive highPos lowPos :=
    protected_target_teeth_not_consecutive
      hHighSafe hLowDanger hConsecutiveOverlap
  exact ⟨nonconsecutive_binary_pair_aligned
      hDichotomy hMismatchAdjacent hNonconsecutive,
    nonconsecutive_has_toward_neighbor hPositionsDistinct hNonconsecutive⟩

/-- The first tooth on the aligned corridor has neither marked owner and its
overlap with the high target lies in the low-safe component.  The assumptions
`hSameOwnerDisjoint` and `hLowOwnerDanger` are exactly the same-comb and
ascent-protection paper providers. -/
theorem corridor_neighbor_is_third_owner_and_protected
    {owner point : Type*}
    (tooth : ℕ → point → Prop) (toothOwner : ℕ → owner)
    (lowSafe : point → Prop)
    {highPos neighborPos : ℕ} {lowOwner highOwner : owner}
    (hAdjacent : Consecutive highPos neighborPos)
    (hHighOwner : toothOwner highPos = highOwner)
    (hHighSafe : ∀ x, tooth highPos x → lowSafe x)
    (hLowOwnerDanger : ∀ p x,
      toothOwner p = lowOwner → tooth p x → ¬ lowSafe x)
    (hAdjacentOverlap : ∀ p q, Consecutive p q →
      ∃ x, tooth p x ∧ tooth q x)
    (hSameOwnerDisjoint : ∀ p q, p ≠ q →
      toothOwner p = toothOwner q →
        ¬ ∃ x, tooth p x ∧ tooth q x) :
    toothOwner neighborPos ≠ highOwner ∧
      toothOwner neighborPos ≠ lowOwner ∧
      ∃ x, lowSafe x ∧ tooth highPos x ∧ tooth neighborPos x := by
  obtain ⟨x, hxHigh, hxNeighbor⟩ :=
    hAdjacentOverlap highPos neighborPos hAdjacent
  have hposne : highPos ≠ neighborPos := consecutive_ne hAdjacent
  have hnotHigh : toothOwner neighborPos ≠ highOwner := by
    intro hNeighborOwner
    have hsame : toothOwner highPos = toothOwner neighborPos :=
      hHighOwner.trans hNeighborOwner.symm
    exact hSameOwnerDisjoint highPos neighborPos hposne hsame
      ⟨x, hxHigh, hxNeighbor⟩
  have hnotLow : toothOwner neighborPos ≠ lowOwner := by
    intro hNeighborOwner
    exact hLowOwnerDanger neighborPos x hNeighborOwner hxNeighbor
      (hHighSafe x hxHigh)
  exact ⟨hnotHigh, hnotLow,
    ⟨x, hHighSafe x hxHigh, hxHigh, hxNeighbor⟩⟩

/-- Alignment is equality of the two nonzero orientation signs. -/
theorem aligned_orientation_sign
    {wordDiff phaseDiff : ℝ}
    (hOrientation :
      (0 < wordDiff ∧ 0 < phaseDiff) ∨
        (wordDiff < 0 ∧ phaseDiff < 0)) :
    (0 < wordDiff ↔ 0 < phaseDiff) ∧
      (wordDiff < 0 ↔ phaseDiff < 0) := by
  rcases hOrientation with hpos | hneg
  · constructor <;> constructor <;> intro <;> linarith
  · constructor <;> constructor <;> intro <;> linarith

/-- Circle reflection negates both differences and preserves alignment. -/
theorem reflection_preserves_alignment
    (wordDiff phaseDiff : ℝ) :
    0 < wordDiff * phaseDiff ↔
      0 < (-wordDiff) * (-phaseDiff) := by
  constructor <;> intro h <;> nlinarith

/-- Natural-number conversion of the gcd endpoint quantum to the lcm clock. -/
theorem gcd_mul_lcm_identity (u v : ℕ) :
    Nat.gcd u v * Nat.lcm u v = u * v := by
  exact Nat.gcd_mul_lcm u v

/-- The two equivalent rational forms of the protected seam quantum. -/
theorem gcd_quantum_eq_lcm_quantum
    {u v : ℕ} (hu : 0 < u) (hv : 0 < v) :
    (Nat.gcd u v : ℚ) / (14 * u * v) =
      1 / (14 * Nat.lcm u v) := by
  have hu0 : (u : ℚ) ≠ 0 := by positivity
  have hv0 : (v : ℚ) ≠ 0 := by positivity
  have hl : 0 < Nat.lcm u v := Nat.lcm_pos hu hv
  have hl0 : (Nat.lcm u v : ℚ) ≠ 0 := by positivity
  field_simp
  exact_mod_cast Nat.gcd_mul_lcm u v

/-- A positive compatible seam numerator is at least the gcd quantum's
integer numerator. -/
theorem positive_divisible_numerator_ge_gcd
    {u v numerator : ℕ}
    (hpos : 0 < numerator) (hdiv : Nat.gcd u v ∣ numerator) :
    Nat.gcd u v ≤ numerator := by
  exact Nat.le_of_dvd hpos hdiv

/-- Combining positivity, gcd divisibility, and the gcd/lcm identity gives
the protected handoff's literal lcm lower bound. -/
theorem positive_compatible_seam_pays_lcm_quantum
    {u v numerator : ℕ} (hu : 0 < u) (hv : 0 < v)
    (hpos : 0 < numerator) (hdiv : Nat.gcd u v ∣ numerator) :
    (1 : ℚ) / (14 * Nat.lcm u v) ≤
      numerator / (14 * u * v) := by
  rw [← gcd_quantum_eq_lcm_quantum hu hv]
  have hleNat : Nat.gcd u v ≤ numerator :=
    positive_divisible_numerator_ge_gcd hpos hdiv
  have hle : (Nat.gcd u v : ℚ) ≤ numerator := by
    exact_mod_cast hleNat
  have hden : (0 : ℚ) < 14 * u * v := by positivity
  exact (div_le_div_iff_of_pos_right hden).2 hle

/-- The lcm seam quantum is strictly positive for positive owner speeds. -/
theorem lcm_quantum_positive
    {u v : ℕ} (hu : 0 < u) (hv : 0 < v) :
    (0 : ℚ) < 1 / (14 * Nat.lcm u v) := by
  have hl : 0 < Nat.lcm u v := Nat.lcm_pos hu hv
  positivity

#print axioms ascent_protection_disjoint
#print axioms reverse_edge_binary_gate
#print axioms protected_target_teeth_not_consecutive
#print axioms nonconsecutive_binary_pair_aligned
#print axioms coherent_two_cycle_forces_alignment
#print axioms nonconsecutive_has_toward_neighbor
#print axioms coherent_two_cycle_alignment_with_corridor
#print axioms corridor_neighbor_is_third_owner_and_protected
#print axioms aligned_orientation_sign
#print axioms reflection_preserves_alignment
#print axioms gcd_quantum_eq_lcm_quantum
#print axioms positive_divisible_numerator_ge_gcd
#print axioms positive_compatible_seam_pays_lcm_quantum
#print axioms lcm_quantum_positive

end BlockerTwoCycleAlignment
end LRC14
