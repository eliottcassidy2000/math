/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Closest-return leaf rank and the sharp paid-star obstruction (THM-1266)

The paper topology supplies a deletion-minimal chronological tooth word, a
closest repeated owner, endpoint order, and the exact return-polygon identity.
This dependency-free module checks the arithmetic consumers: the strict `6/5`
leaf ascent, the packet-packing count, the 42-step compact rank, the repeated
low-owner separation which forbids a sixth consecutive-address star rung, and
the exact rational data of the primitive sharp cell.
-/

namespace LRC14
namespace ClosestReturnLeafPaidStar

/-- Once topology bounds a closest return by six edges, positivity of its
return polygon forces the repeated outer owner above `6/5` of the smallest
internal owner.  `ratioSum` is the sum of all internal `a / s_q` terms. -/
theorem closest_return_forces_six_fifths
    {a b ratioSum R : ℚ} {m : ℕ}
    (ha : 0 < a) (hb : 0 < b)
    (hmLower : 2 ≤ m) (hmUpper : m ≤ 6) (hR : 1 ≤ R)
    (hReturn : 7 * R - 1 < ratioSum)
    (hMinimum : ratioSum ≤ (m - 1 : ℕ) * (a / b)) :
    6 * b < 5 * a := by
  have hCount : m - 1 ≤ 5 := by omega
  have hCountQ : ((m - 1 : ℕ) : ℚ) ≤ 5 := by exact_mod_cast hCount
  have hRatioNonnegative : 0 ≤ a / b := by positivity
  have hCountBound : ((m - 1 : ℕ) : ℚ) * (a / b) ≤ 5 * (a / b) := by
    exact mul_le_mul_of_nonneg_right hCountQ hRatioNonnegative
  have hSix : (6 : ℚ) ≤ 7 * R - 1 := by linarith
  have hRatio : (6 : ℚ) < 5 * (a / b) :=
    lt_of_le_of_lt hSix (lt_of_lt_of_le hReturn (hMinimum.trans hCountBound))
  have hCancel : (a / b) * b = a := by field_simp
  nlinarith

/-- The recursive left/right extraction of vertex-disjoint closest packets
uses at most seven vertices per packet and leaves at most `P+1` distinct-owner
blocks of at most six vertices. -/
theorem disjoint_packet_packing_bound
    {N P leaves packetVertices terminalVertices : ℕ}
    (hPartition : N = packetVertices + terminalVertices)
    (hPackets : packetVertices ≤ 7 * P)
    (hLeaves : leaves ≤ P + 1)
    (hTerminals : terminalVertices ≤ 6 * leaves) :
    N ≤ 13 * P + 6 := by
  omega

/-- The weakest closest-return factor fits 42 strict projective ascents in
the compact box, but a 43rd ascent is impossible. -/
theorem six_fifths_rank_cutoff :
    ((6 : ℚ) / 5) ^ 42 < 2345 ∧
      (2345 : ℚ) < ((6 : ℚ) / 5) ^ 43 := by
  norm_num

/-- Endpoint separation for two occurrences of one low owner.  The paper
chain gives `1/j ≤ x < (7 delta + 1)/(7h)`; one r=1 toothpick gives `h>6j`.
Together they force the two low slots at least six high-gap addresses apart. -/
theorem repeated_low_owner_separation
    {h j x : ℚ} {delta : ℕ}
    (hh : 0 < h) (hj : 0 < j) (hAscent : 6 * j < h)
    (hLeft : 1 / j ≤ x)
    (hRight : x < (7 * (delta : ℚ) + 1) / (7 * h)) :
    6 ≤ delta := by
  have hDenominator : 0 < 7 * h := by positivity
  have hFraction : 1 / j < (7 * (delta : ℚ) + 1) / (7 * h) :=
    lt_of_le_of_lt hLeft hRight
  have hCross : 7 * h < j * (7 * (delta : ℚ) + 1) := by
    have := (div_lt_div_iff₀ hj hDenominator).mp hFraction
    nlinarith
  by_contra hNot
  have hDelta : delta ≤ 5 := by omega
  have hDeltaQ : (delta : ℚ) ≤ 5 := by exact_mod_cast hDelta
  nlinarith

/-- Six slots colored by only five non-high owners repeat a color. -/
theorem six_low_slots_repeat (owner : Fin 6 → Fin 5) :
    ∃ i j, i ≠ j ∧ owner i = owner j := by
  obtain ⟨i, j, hne, heq⟩ :=
    Fintype.exists_ne_map_eq_of_card_lt owner (by norm_num)
  exact ⟨i, j, hne, heq⟩

/-- A six-rung consecutive-address star would color six slots by five low
owners, while repeated-low separation requires equal colors six slots apart.
That is impossible inside `Fin 6`. -/
theorem no_six_rung_consecutive_address_star
    (owner : Fin 6 → Fin 5)
    (hSeparated : ∀ i j, i ≠ j → owner i = owner j →
      i.val + 6 ≤ j.val ∨ j.val + 6 ≤ i.val) :
    False := by
  obtain ⟨i, j, hne, heq⟩ := six_low_slots_repeat owner
  have hFar := hSeparated i j hne heq
  have hi := i.isLt
  have hj := j.isLt
  omega

noncomputable def toothLeft (speed address : ℕ) : ℚ :=
  (14 * (address : ℚ) - 1) / (14 * (speed : ℚ))

noncomputable def toothRight (speed address : ℕ) : ℚ :=
  (14 * (address : ℚ) + 1) / (14 * (speed : ℚ))

def sharpRows : List (ℕ × ℕ) :=
  [(1805, 1036), (256, 147), (1805, 1037), (254, 146),
   (1805, 1038), (292, 168), (1805, 1039), (257, 148),
   (1805, 1040), (255, 147), (1805, 1041)]

def AdjacentOverlap : List (ℕ × ℕ) → Prop
  | first :: second :: tail =>
      toothLeft second.1 second.2 < toothRight first.1 first.2 ∧
        AdjacentOverlap (second :: tail)
  | _ => True

def TwoApartSeparated : List (ℕ × ℕ) → Prop
  | first :: second :: third :: tail =>
      toothRight first.1 first.2 < toothLeft third.1 third.2 ∧
        TwoApartSeparated (second :: third :: tail)
  | _ => True

/-- Both endpoint sequences of the sharp eleven-tooth row are strictly
increasing. -/
theorem sharp_star_endpoint_order :
    sharpRows.Pairwise (fun first second =>
      toothLeft first.1 first.2 < toothLeft second.1 second.2) ∧
    sharpRows.Pairwise (fun first second =>
      toothRight first.1 first.2 < toothRight second.1 second.2) := by
  norm_num [sharpRows, toothLeft, toothRight]

/-- Adjacent teeth overlap and every pair two positions apart is strictly
separated, so the row is an irredundant interval chain. -/
theorem sharp_star_chain_geometry :
    AdjacentOverlap sharpRows ∧ TwoApartSeparated sharpRows := by
  norm_num [sharpRows, AdjacentOverlap, TwoApartSeparated, toothLeft, toothRight]

/-- The sharp cell lies strictly inside the complete `c=140,k=80` carrier
gap. -/
theorem sharp_star_inside_carrier :
    (14 * (80 : ℚ) + 1) / (14 * 140) < toothLeft 1805 1036 ∧
      toothRight 1805 1041 < (14 * (80 : ℚ) + 13) / (14 * 140) := by
  norm_num [toothLeft, toothRight]

/-- The five immediate returns have the displayed positive detunings. -/
theorem sharp_star_detunings :
    1805 - 6 * 256 = 269 ∧
    1805 - 6 * 254 = 281 ∧
    1805 - 6 * 292 = 53 ∧
    1805 - 6 * 257 = 263 ∧
    1805 - 6 * 255 = 275 := by
  norm_num

/-- The row is primitive, compact, and lies on the cover-compatible side of
the first harmonic scalar gate. -/
theorem sharp_star_projective_checks :
    Nat.gcd 140 (Nat.gcd 254 (Nat.gcd 255 (Nat.gcd 256
      (Nat.gcd 257 (Nat.gcd 292 1805))))) = 1 ∧
    1805 < 2345 * 140 ∧
    (1 : ℚ) / 140 <
      1 / 254 + 1 / 255 + 1 / 256 + 1 / 257 + 1 / 292 + 1 / 1805 := by
  norm_num

/-- Every selected centered blocker in the sharp row strictly contains its
phase.  The target word positions are respectively `8,3,7,1,1,3`. -/
theorem sharp_centered_blocker_membership :
    toothLeft 1805 1040 < (227 : ℚ) / 394 ∧
      (227 : ℚ) / 394 < toothRight 1805 1040 ∧
    toothLeft 254 146 < (227 : ℚ) / 395 ∧
      (227 : ℚ) / 395 < toothRight 254 146 ∧
    toothLeft 257 148 < (19 : ℚ) / 33 ∧
      (19 : ℚ) / 33 < toothRight 257 148 ∧
    toothLeft 256 147 < (228 : ℚ) / 397 ∧
      (228 : ℚ) / 397 < toothRight 256 147 ∧
    toothLeft 256 147 < (31 : ℚ) / 54 ∧
      (31 : ℚ) / 54 < toothRight 256 147 ∧
    toothLeft 254 146 < (1118 : ℚ) / 1945 ∧
      (1118 : ℚ) / 1945 < toothRight 254 146 := by
  norm_num [toothLeft, toothRight]

/-- The two local blocker cycles are phase/word aligned and nonconsecutive. -/
theorem sharp_two_cycle_alignment_signs :
    (1118 : ℚ) / 1945 - 227 / 394 = -1023 / 766330 ∧
    (228 : ℚ) / 397 - 19 / 33 = -19 / 13101 ∧
    (-1023 : ℚ) / 766330 < 0 ∧
    (-19 : ℚ) / 13101 < 0 ∧
    3 + 1 < 8 ∧ 1 + 1 < 7 := by
  norm_num

/-- The four residual components are complete fastest-safe gaps, each of
length `6/(7*1805)`, after fastest addresses 1034, 1035, 1041, and 1042. -/
theorem sharp_four_tail_safe_gaps :
    toothRight 1805 1034 = (14477 : ℚ) / 25270 ∧
    toothLeft 1805 1035 = (14489 : ℚ) / 25270 ∧
    toothRight 1805 1035 = (14491 : ℚ) / 25270 ∧
    toothLeft 1805 1036 = (14503 : ℚ) / 25270 ∧
    toothRight 1805 1041 = (2915 : ℚ) / 5054 ∧
    toothLeft 1805 1042 = (14587 : ℚ) / 25270 ∧
    toothRight 1805 1042 = (14589 : ℚ) / 25270 ∧
    toothLeft 1805 1043 = (14601 : ℚ) / 25270 ∧
    (14489 : ℚ) / 25270 - 14477 / 25270 = 6 / (7 * 1805) ∧
    (14503 : ℚ) / 25270 - 14491 / 25270 = 6 / (7 * 1805) ∧
    (14587 : ℚ) / 25270 - 2915 / 5054 = 6 / (7 * 1805) ∧
    (14601 : ℚ) / 25270 - 14589 / 25270 = 6 / (7 * 1805) := by
  norm_num [toothLeft, toothRight]

#print axioms closest_return_forces_six_fifths
#print axioms disjoint_packet_packing_bound
#print axioms six_fifths_rank_cutoff
#print axioms repeated_low_owner_separation
#print axioms six_low_slots_repeat
#print axioms no_six_rung_consecutive_address_star
#print axioms sharp_star_endpoint_order
#print axioms sharp_star_chain_geometry
#print axioms sharp_star_inside_carrier
#print axioms sharp_star_detunings
#print axioms sharp_star_projective_checks
#print axioms sharp_centered_blocker_membership
#print axioms sharp_two_cycle_alignment_signs
#print axioms sharp_four_tail_safe_gaps

end ClosestReturnLeafPaidStar
end LRC14
