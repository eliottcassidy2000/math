/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import TournamentH7.LRCCenteredBlockerAddressCompression

/-!
# Binary centered phases land on the irredundant tooth word (THM-1256)

The paper layer combines the interval-chain topology of THM-1253 with the
coherent binary edge of THM-1254.  This module checks the exact endpoint
interpretation of the residual invoice, the local alignment/mismatch interval
lemmas, the same-provider toothpick detuning, the ABAB contradiction, and the
resulting turn-count arithmetic.
-/

namespace LRC14
namespace BinaryPhaseWordLanding

/-- Reciprocal-centre drift is exactly the determinant of the two endpoint
tooth centres. -/
theorem reciprocal_endpoint_identity
    (s₀ sᵣ n₀ nᵣ : ℚ) (hn₀ : n₀ ≠ 0) (hnᵣ : nᵣ ≠ 0) :
    n₀ * nᵣ * (s₀ / n₀ - sᵣ / nᵣ) = nᵣ * s₀ - n₀ * sᵣ := by
  field_simp [hn₀, hnᵣ]

/-- The apparent residual surplus is just the marked endpoint-address
difference. -/
theorem residual_minus_endpoint
    (P n₀ nᵣ s₀ sᵣ k delta : ℚ)
    (hP : P = nᵣ + k + delta) :
    (P * s₀ - n₀ * sᵣ) - (nᵣ * s₀ - n₀ * sᵣ) =
      s₀ * (k + delta) := by
  rw [hP]
  ring

/-- Consequently the THM-1254 invoice is address order, not a second metric
overlap credit. -/
theorem residual_invoice_iff_address_order
    {P n₀ nᵣ s₀ sᵣ : ℚ} (hs₀ : 0 < s₀) :
    nᵣ * s₀ - n₀ * sᵣ ≤ P * s₀ - n₀ * sᵣ ↔ nᵣ ≤ P := by
  constructor <;> intro h <;> nlinarith

/-- Binary digit zero puts the source centered phase after the target phase. -/
theorem digit_zero_phase_order
    {Q theta phaseDiff : ℚ}
    (hQ : 0 < Q) (htheta : 0 < theta)
    (htransport : Q * phaseDiff = theta) :
    0 < phaseDiff := by
  nlinarith

/-- Binary digit one puts the source centered phase before the target phase. -/
theorem digit_one_phase_order
    {Q theta phaseDiff : ℚ}
    (hQ : 0 < Q) (htheta : theta < 1)
    (htransport : Q * phaseDiff = theta - 1) :
    phaseDiff < 0 := by
  nlinarith

/-- THM-1248 doubles the old central-band separation on an actual digit-zero
blocker edge. -/
theorem digit_zero_sharp_phase_separation
    {Q theta phaseDiff : ℚ}
    (hQ : 0 < Q) (htheta : 5/14 < theta)
    (htransport : Q * phaseDiff = theta) :
    5/(14*Q) < phaseDiff := by
  have hden : 0 < 14*Q := by positivity
  rw [div_lt_iff₀ hden]
  nlinarith

/-- The reflected digit-one form of the doubled actual-blocker separation. -/
theorem digit_one_sharp_phase_separation
    {Q theta phaseDiff : ℚ}
    (hQ : 0 < Q) (htheta : theta < 9/14)
    (htransport : Q * phaseDiff = theta-1) :
    5/(14*Q) < -phaseDiff := by
  have hden : 0 < 14*Q := by positivity
  rw [div_lt_iff₀ hden]
  nlinarith

/-- The key landing upgrade: if the target's own centered phase lies outside
the tooth which blocked the source, then phase order and both endpoint orders
cannot mismatch.  In a minimal word this is exactly chronological alignment. -/
theorem coherent_blocker_marks_align
    {aL aR bL bR x y : ℝ}
    (hxA : aL < x) (hxA' : x < aR)
    (hyB : bL < y) (hyB' : y < bR)
    (houtside : y < aL ∨ aR < y) :
    (y < x ∧ bL < aL) ∨ (x < y ∧ aR < bR) := by
  rcases houtside with hleft | hright
  · left
    constructor <;> linarith
  · right
    constructor <;> linarith

/-- The doubled phase clearance is wider than the whole target tooth by more
than `1/(14Q)` whenever `Q=c+j` and `j>c`. -/
theorem sharp_corridor_exceeds_target_tooth
    {c j Q : ℝ} (hc : 0 < c) (hcj : c < j) (hQ : Q = c+j) :
    1/(7*j)+1/(14*Q) < 5/(14*Q) := by
  have hj : 0 < j := lt_trans hc hcj
  have hQpos : 0 < Q := by rw [hQ]; nlinarith
  have hjden : 0 < 7*j := by positivity
  have hQden : 0 < 14*Q := by positivity
  rw [div_add_div, div_lt_div_iff₀] <;> try positivity
  rw [hQ]
  nlinarith

/-- If one following tooth spans the sharp target-wall corridor, its owner
must be smaller than `4j/5`.  Hence this adjacency is impossible when `j` is
the minimum speed on the blocker cycle. -/
theorem adjacent_marked_tooth_forces_speed_drop
    {c j h Q leg : ℝ}
    (hc : 0 < c) (hcj : c < j) (hh : 0 < h)
    (hQ : Q = c+j)
    (hlower : 5/(14*Q) < leg)
    (hspan : leg < 1/(7*h)) :
    5*h < 2*Q ∧ 5*h < 4*j := by
  have hj : 0 < j := lt_trans hc hcj
  have hQpos : 0 < Q := by rw [hQ]; nlinarith
  have hleft : 5/(14*Q) < 1/(7*h) := lt_trans hlower hspan
  have hdenQ : 0 < 14*Q := by positivity
  have hdenh : 0 < 7*h := by positivity
  have hcross := (div_lt_div_iff₀ hdenQ hdenh).mp hleft
  constructor
  · nlinarith
  · rw [hQ] at hcross
    nlinarith

/-- At a strict cycle minimum the preceding lemma contradicts adjacency. -/
theorem cycle_minimum_marked_teeth_not_adjacent
    {c j h Q leg : ℝ}
    (hc : 0 < c) (hcj : c < j) (hjh : j < h)
    (hQ : Q = c+j)
    (hlower : 5/(14*Q) < leg)
    (hspan : leg < 1/(7*h)) : False := by
  have hh : 0 < h := lt_trans (lt_trans hc hcj) hjh
  have hdrop := adjacent_marked_tooth_forces_speed_drop
    hc hcj hh hQ hlower hspan
  nlinarith

/-- The target-free corridor is covered by the other five combs.  Combining
the sharp one-interval discrepancy with its `>5/(14Q)` length gives the
explicit harmonic side invoice `H>5/(6Q)`. -/
theorem five_comb_corridor_harmonic_invoice
    {L H Q : ℝ}
    (hL : 5/(14*Q) < L)
    (hcover : L ≤ 5*L/7+6*H/49) :
    5/(6*Q) < H := by
  calc
    5/(6*Q) = (7/3)*(5/(14*Q)) := by ring
    _ < (7/3)*L := by
      exact mul_lt_mul_of_pos_left hL (by norm_num)
    _ ≤ H := by nlinarith

/-- The sharp binary corridor is shorter than one target safe gap, so after
leaving its selected target tooth it cannot reach the next target tooth. -/
theorem sharp_binary_corridor_shorter_than_target_safe_gap
    {c j Q : ℝ} (hc : 0 < c) (hcj : c < j) (hQ : Q = c+j) :
    23/(28*Q) < 6/(7*j) := by
  have hj : 0 < j := lt_trans hc hcj
  have hQpos : 0 < Q := by rw [hQ]; positivity
  have hdenQ : 0 < 28*Q := by positivity
  have hdenj : 0 < 7*j := by positivity
  apply (div_lt_div_iff₀ hdenQ hdenj).2
  rw [hQ]
  nlinarith

/-- If a full internal handoff is counted twice inside the five-comb
corridor, the one-interval discrepancy yields an additive `49/6` seam term. -/
theorem nested_corridor_tariff_rearrangement
    {L W H : ℝ} (hcover : L+W ≤ 5*L/7+6*H/49) :
    (7/3)*L+(49/6)*W ≤ H := by
  nlinarith

/-- In the binary-ascent two-cycle subbranch, four outside combs cover the
facing gap; the same discrepancy algebra turns a seam quantum into the
coefficient `1/4`. -/
theorem four_comb_facing_gap_invoice
    {L H : ℝ} (hcover : L ≤ 4*L/7+6*H/49) :
    (7/2)*L ≤ H := by
  nlinarith

theorem four_comb_seam_coefficient (q : ℝ) :
    (7/2)*(1/(14*q)) = 1/(4*q) := by
  ring

/-- Three located gcd seams in the protected component change THM-1244's
coefficient from `7/6` to `7/4`. -/
theorem three_located_seams_coefficient (g D : ℝ) :
    (49/6)*(3*g/(14*D^2)) = 7*g/(4*D^2) := by
  ring

/-- Normalized reciprocal obstruction for the range `2<x≤3`, proved without
calculus by splitting once at `x=5/2`. -/
theorem binary_six_cycle_ratio_obstruction_low
    {x : ℝ} (hx : 2 < x) (hx3 : x ≤ 3) :
    7+x < 5*(x/(x+1)+x/(x+2)+x/(x+3)+x/(x+4)) := by
  by_cases hmid : x ≤ 5/2
  · have h1 : (2:ℝ)/3 < x/(x+1) := by
      apply (div_lt_div_iff₀ (by norm_num) (by nlinarith)).2
      nlinarith
    have h2 : (1:ℝ)/2 < x/(x+2) := by
      apply (div_lt_div_iff₀ (by norm_num) (by nlinarith)).2
      nlinarith
    have h3 : (2:ℝ)/5 < x/(x+3) := by
      apply (div_lt_div_iff₀ (by norm_num) (by nlinarith)).2
      nlinarith
    have h4 : (1:ℝ)/3 < x/(x+4) := by
      apply (div_lt_div_iff₀ (by norm_num) (by nlinarith)).2
      nlinarith
    nlinarith
  · have hxmid : 5/2 < x := lt_of_not_ge hmid
    have h1 : (5:ℝ)/7 < x/(x+1) := by
      apply (div_lt_div_iff₀ (by norm_num) (by nlinarith)).2
      nlinarith
    have h2 : (5:ℝ)/9 < x/(x+2) := by
      apply (div_lt_div_iff₀ (by norm_num) (by nlinarith)).2
      nlinarith
    have h3 : (5:ℝ)/11 < x/(x+3) := by
      apply (div_lt_div_iff₀ (by norm_num) (by nlinarith)).2
      nlinarith
    have h4 : (5:ℝ)/13 < x/(x+4) := by
      apply (div_lt_div_iff₀ (by norm_num) (by nlinarith)).2
      nlinarith
    have hconst : (10:ℝ) < 5*(5/7+5/9+5/11+5/13) := by norm_num
    nlinarith

/-- The complementary normalized reciprocal obstruction for `x≥3`. -/
theorem binary_six_cycle_ratio_obstruction_high
    {x : ℝ} (hx : 3 ≤ x) :
    7+2*x/(x-1) < 5*(x/(x+1)+x/(x+2)+x/(x+3)+x/(x+4)) := by
  have h1 : (3:ℝ)/4 ≤ x/(x+1) := by
    apply (div_le_div_iff₀ (by norm_num) (by nlinarith)).2
    nlinarith
  have h2 : (3:ℝ)/5 ≤ x/(x+2) := by
    apply (div_le_div_iff₀ (by norm_num) (by nlinarith)).2
    nlinarith
  have h3 : (1:ℝ)/2 ≤ x/(x+3) := by
    apply (div_le_div_iff₀ (by norm_num) (by nlinarith)).2
    nlinarith
  have h4 : (3:ℝ)/7 ≤ x/(x+4) := by
    apply (div_le_div_iff₀ (by norm_num) (by nlinarith)).2
    nlinarith
  have hratio : 2*x/(x-1) ≤ 3 := by
    apply (div_le_iff₀ (by nlinarith)).2
    nlinarith
  have hconst : (10:ℝ) < 5*(3/4+3/5+1/2+3/7) := by norm_num
  nlinarith

/-- Exact clock-ratio strip when the slowest two-cycle's descent digit is
zero and its ascent digit is `m≥1`. -/
theorem two_cycle_descent_zero_ratio_strip
    {R m thetaUp downMag : ℝ}
    (hm : 1 ≤ m)
    (huLo : 5/28 < thetaUp) (huHi : thetaUp < 23/28)
    (hdLo : 5/14 < downMag) (hdHi : downMag < 23/28)
    (hR : R = (m-thetaUp)/downMag) :
    (28*m-23)/23 < R ∧ R < (28*m-5)/10 := by
  have hdpos : 0 < downMag := by nlinarith
  have hloPos : 0 < (28*m-23)/23 := by nlinarith
  have hupPos : 0 < (28*m-5)/10 := by nlinarith
  constructor
  · rw [hR]
    apply (lt_div_iff₀ hdpos).2
    calc
      ((28*m-23)/23)*downMag
          < ((28*m-23)/23)*(23/28) :=
            mul_lt_mul_of_pos_left hdHi hloPos
      _ = m-23/28 := by ring
      _ < m-thetaUp := by linarith
  · rw [hR]
    apply (div_lt_iff₀ hdpos).2
    calc
      m-thetaUp < m-5/28 := by linarith
      _ = ((28*m-5)/10)*(5/14) := by ring
      _ < ((28*m-5)/10)*downMag :=
        mul_lt_mul_of_pos_left hdLo hupPos

/-- Reflected clock-ratio strip when the descent digit is one and the ascent
digit is `-p`, `p≥0`. -/
theorem two_cycle_descent_one_ratio_strip
    {R p thetaUp downMag : ℝ}
    (hp : 0 ≤ p)
    (huLo : 5/28 < thetaUp) (huHi : thetaUp < 23/28)
    (hdLo : 5/14 < downMag) (hdHi : downMag < 23/28)
    (hR : R = (p+thetaUp)/downMag) :
    (28*p+5)/23 < R ∧ R < (28*p+23)/10 := by
  have hdpos : 0 < downMag := by nlinarith
  have hloPos : 0 < (28*p+5)/23 := by nlinarith
  have hupPos : 0 < (28*p+23)/10 := by nlinarith
  constructor
  · rw [hR]
    apply (lt_div_iff₀ hdpos).2
    calc
      ((28*p+5)/23)*downMag
          < ((28*p+5)/23)*(23/28) :=
            mul_lt_mul_of_pos_left hdHi hloPos
      _ = p+5/28 := by ring
      _ < p+thetaUp := by linarith
  · rw [hR]
    apply (div_lt_iff₀ hdpos).2
    calc
      p+thetaUp < p+23/28 := by linarith
      _ = ((28*p+23)/10)*(5/14) := by ring
      _ < ((28*p+23)/10)*downMag :=
        mul_lt_mul_of_pos_left hdLo hupPos

/-- Nonadjacent marked teeth contribute at least two corridor handoffs; the
opposite wall of the first marked tooth supplies a third occurrence. -/
theorem nonadjacent_marks_export_three_handoffs
    {p q : ℕ} (hsep : p+2 ≤ q) : 3 ≤ (q-p)+1 := by
  omega

/-- In the incidence/reuse-free two-cycle residual, distinct owners on the
two facing walls rule out a unique intermediate tooth.  Including both outer
wall handoffs leaves a canonical path with at least five occurrences. -/
theorem distinct_facing_owners_export_five_handoffs
    {p q : ℕ} (hsep : p+2 ≤ q) (hnotOne : q ≠ p+2) :
    5 ≤ (q-p)+2 := by
  omega

/-- In a four-tooth window, the ABAB exclusion says that a same-provider
target fork is followed by a genuinely different projected edge. -/
theorem same_provider_fork_forces_next_turn
    {α : Type*} {h₀ j h₁ h₂ : α}
    (hnoABAB : ¬ (h₀ = h₁ ∧ j = h₂)) :
    h₀ ≠ h₁ ∨ j ≠ h₂ := by
  tauto

/-- If interval A precedes interval B but their marked interior phases occur
in the reverse order, A and B overlap. -/
theorem phase_position_mismatch_forces_overlap
    {aL aR bL bR x y : ℝ}
    (_hleft : aL < bL)
    (_hxA : aL < x) (hxA' : x < aR)
    (hyB : bL < y) (_hyB' : y < bR)
    (hmismatch : y < x) :
    bL < aR := by
  linarith

/-- THM-1253 separation therefore excludes a mismatch between nonconsecutive
members of the minimal interval word. -/
theorem separated_intervals_cannot_mismatch
    {aL aR bL bR x y : ℝ}
    (hleft : aL < bL)
    (hxA : aL < x) (hxA' : x < aR)
    (hyB : bL < y) (hyB' : y < bR)
    (hseparated : aR ≤ bL) :
    ¬ y < x := by
  intro hmismatch
  have := phase_position_mismatch_forces_overlap
    hleft hxA hxA' hyB hyB' hmismatch
  linarith

/-- Two adjacent phase/position mismatches would reverse the two
nonconsecutive outer marks, contradicting THM-1253 separation.  Thus the
inversion graph on marked teeth is a matching. -/
theorem adjacent_mismatches_cannot_share_a_mark
    {aR cL x y z : ℝ}
    (hxA : x < aR) (hzC : cL < z)
    (hseparated : aR ≤ cL)
    (hxy : y < x) (hyz : z < y) :
    False := by
  linarith

/-- Two overlapping ordered intervals cover every point between an interior
mark of the first and a later interior mark of the second.  Iteration gives
the paper subword statement. -/
theorem aligned_overlapping_pair_covers_segment
    {aL aR bL bR x y z : ℝ}
    (_hleft : aL < bL) (_hright : aR < bR)
    (hoverlap : bL < aR)
    (hxA : aL < x) (_hxA' : x < aR)
    (_hyB : bL < y) (hyB' : y < bR)
    (hxz : x ≤ z) (hzy : z ≤ y) :
    (aL < z ∧ z < aR) ∨ (bL < z ∧ z < bR) := by
  by_cases hzA : z < aR
  · left
    exact ⟨lt_of_lt_of_le hxA hxz, hzA⟩
  · right
    have hzA' : aR ≤ z := le_of_not_gt hzA
    exact ⟨lt_of_lt_of_le hoverlap hzA', lt_of_le_of_lt hzy hyB'⟩

/-- Exact same-provider toothpick law for a backtrack with address rung r. -/
theorem backtrack_detuning_identity
    (h j r : ℚ) (hh : h ≠ 0) (hj : j ≠ 0) :
    1 / (7 * j) - (7 * r - 1) / (7 * h) =
      (h - (7 * r - 1) * j) / (7 * j * h) := by
  field_simp [hh, hj]

/-- A positive backtrack detuning forces the repeated outer owner to be more
than six times the middle owner. -/
theorem backtrack_outer_gt_six_middle
    {h j r : ℤ}
    (hj : 0 < j) (hr : 1 ≤ r)
    (hdetune : 0 < h - (7 * r - 1) * j) :
    6 * j < h := by
  have hc : (6 : ℤ) ≤ 7 * r - 1 := by omega
  have hm : 6 * j ≤ (7 * r - 1) * j :=
    mul_le_mul_of_nonneg_right hc (le_of_lt hj)
  have hd : (7 * r - 1) * j < h := by omega
  exact lt_of_le_of_lt hm hd

/-- Oppositely oriented backtracks on one unordered owner pair are
incompatible.  This is the arithmetic core of the ABAB exclusion. -/
theorem abab_speed_contradiction
    {a b : ℤ} (_ha : 0 < a) (hb : 0 < b)
    (hab : 6 * b < a) (hba : 6 * a < b) :
    False := by
  nlinarith

/-- Five successive toothpick returns cannot fit in the compact LRC14 speed
box: each return multiplies speed by more than six, while `6^5>2345`. -/
theorem no_five_backtrack_ascent_chain
    {c x₀ x₁ x₂ x₃ x₄ x₅ : ℝ}
    (hc : 0 < c) (hx₀ : c < x₀)
    (h01 : 6 * x₀ < x₁) (h12 : 6 * x₁ < x₂)
    (h23 : 6 * x₂ < x₃) (h34 : 6 * x₃ < x₄)
    (h45 : 6 * x₄ < x₅) (hcap : x₅ < 2345 * c) :
    False := by
  nlinarith

/-- If B internal positions are backtracks and T are nonbacktracking turns,
the absence of consecutive backtracks gives the sharp half-density turn
floor. -/
theorem turn_floor_from_no_consecutive_backtracks
    {N B T : ℕ}
    (hN : 2 ≤ N)
    (hledger : B + T = N - 2)
    (hseparated : 2 * B ≤ N - 1) :
    (N - 2) / 2 ≤ T := by
  omega

#print axioms reciprocal_endpoint_identity
#print axioms residual_minus_endpoint
#print axioms residual_invoice_iff_address_order
#print axioms digit_zero_phase_order
#print axioms digit_one_phase_order
#print axioms digit_zero_sharp_phase_separation
#print axioms digit_one_sharp_phase_separation
#print axioms coherent_blocker_marks_align
#print axioms sharp_corridor_exceeds_target_tooth
#print axioms adjacent_marked_tooth_forces_speed_drop
#print axioms cycle_minimum_marked_teeth_not_adjacent
#print axioms five_comb_corridor_harmonic_invoice
#print axioms sharp_binary_corridor_shorter_than_target_safe_gap
#print axioms nested_corridor_tariff_rearrangement
#print axioms four_comb_facing_gap_invoice
#print axioms four_comb_seam_coefficient
#print axioms three_located_seams_coefficient
#print axioms binary_six_cycle_ratio_obstruction_low
#print axioms binary_six_cycle_ratio_obstruction_high
#print axioms two_cycle_descent_zero_ratio_strip
#print axioms two_cycle_descent_one_ratio_strip
#print axioms nonadjacent_marks_export_three_handoffs
#print axioms distinct_facing_owners_export_five_handoffs
#print axioms same_provider_fork_forces_next_turn
#print axioms phase_position_mismatch_forces_overlap
#print axioms separated_intervals_cannot_mismatch
#print axioms adjacent_mismatches_cannot_share_a_mark
#print axioms aligned_overlapping_pair_covers_segment
#print axioms backtrack_detuning_identity
#print axioms backtrack_outer_gt_six_middle
#print axioms abab_speed_contradiction
#print axioms no_five_backtrack_ascent_chain
#print axioms turn_floor_from_no_consecutive_backtracks

end BinaryPhaseWordLanding
end LRC14
