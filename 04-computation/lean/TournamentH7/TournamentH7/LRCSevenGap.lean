/-
  TournamentH7.LRCSevenGap — the SEVEN-POINT GAP-DEFICIT LEMMA
  (kind-pasteur-2026-07-02-S24, HYP-3981, stage A).

  THE CLUSTERED-BLOCK HEART.  Seven arcs of width `1/7` (the 1/14-danger arcs of
  seven runners) can cover the circle ONLY by perfect packing: all seven circular
  spacings exactly `1/7`.  With the phase of runner 1 PINNED at `0`, perfect
  packing forces the offset multiset to be exactly `{0, 1/7, …, 6/7}`, whose sum
  is `3 ≡ 0 (mod 1)`.  So a citation margin `θ` on the SUM-COMBO
  `C₀ = Σ (wⱼ − w₁)` — a SINGLE extra cited runner — quantitatively precludes
  packing: some spacing is `≥ 1/7 + θ/22`, and its midpoint clears every point
  by `≥ 1/14 + θ/44` — STRICTLY above the sharp margin, with room for drift.

  The sum identity behind it: for sorted positions `0 = q₀ ≤ … ≤ q₆ < 1`,
      Σⱼ qⱼ − 3 = Σᵢ (6−i)·(sᵢ − 1/7)     (sᵢ = the spacings),
  so all spacings `< 1/7 + θ/22` forces `|Σ qⱼ − 3| < θ` — contradicting the
  transported citation margin `∀ m, θ ≤ |Σ qⱼ − m|` at `m = 3`.

  Kernel-pure: no native_decide, no sorry.
-/
import Mathlib
import TournamentH7.LRCRealRegions

namespace LonelyRunner
namespace SevenGap

/-- Mid-gap clearance: the midpoint of a point-free arc `(a, b) ⊆ [0, 1]` keeps
integer-circle distance `≥ c` from every point of `[0,1)` outside the open arc,
whenever the arc is `2c` wide. -/
theorem mid_gap_clearance (a b c : ℝ) (hc : 0 ≤ c) (h2c : 2 * c ≤ b - a)
    (ha : 0 ≤ a) (hb : b ≤ 1) :
    ∀ x : ℝ, 0 ≤ x → x < 1 → (x ≤ a ∨ b ≤ x) → ∀ m : ℤ, c ≤ |(a + b) / 2 - x - m| := by
  intro x hx0 hx1 hside m
  rcases hside with hxa | hxb
  · -- x left of the gap: (a+b)/2 - x ∈ [c, 1 - c]
    have h1 : c ≤ (a + b) / 2 - x := by linarith
    have h2 : (a + b) / 2 - x ≤ 1 - c := by linarith
    rcases le_or_gt (m : ℝ) 0 with hm | hm
    · have : c ≤ (a + b) / 2 - x - m := by linarith
      calc c ≤ (a + b) / 2 - x - m := this
        _ ≤ |(a + b) / 2 - x - m| := le_abs_self _
    · have hm1 : (1 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm
      have : c ≤ -((a + b) / 2 - x - m) := by linarith
      calc c ≤ -((a + b) / 2 - x - m) := this
        _ ≤ |(a + b) / 2 - x - m| := neg_le_abs _
  · -- x right of the gap: x - (a+b)/2 ∈ [c, 1 - c)
    have h1 : c ≤ x - (a + b) / 2 := by linarith
    have h2 : x - (a + b) / 2 ≤ 1 - c := by linarith
    rcases le_or_gt 0 (m : ℝ) with hm | hm
    · have : c ≤ -((a + b) / 2 - x - m) := by linarith
      calc c ≤ -((a + b) / 2 - x - m) := this
        _ ≤ |(a + b) / 2 - x - m| := neg_le_abs _
    · have hm1 : (m : ℝ) ≤ -1 := by
        have : (m : ℤ) < 0 := by exact_mod_cast hm
        have : (m : ℤ) ≤ -1 := by omega
        exact_mod_cast this
      have : c ≤ (a + b) / 2 - x - m := by linarith
      calc c ≤ (a + b) / 2 - x - m := this
        _ ≤ |(a + b) / 2 - x - m| := le_abs_self _

/-- **THE SEVEN-POINT GAP-DEFICIT LEMMA**: sorted circle points
`0 = q₀ ≤ … ≤ q₆ < 1` whose SUM keeps integer distance `≥ θ` admit a probe `φ`
with circle clearance `≥ 1/14 + θ/44` from every point.  (Perfect packing has
sum `3`; the margin on the sum forces a spacing surplus `θ/22` somewhere.) -/
theorem seven_gap_deficit (q0 q1 q2 q3 q4 q5 q6 θ : ℝ) (hθ : 0 < θ)
    (h0 : q0 = 0)
    (h01 : q0 ≤ q1) (h12 : q1 ≤ q2) (h23 : q2 ≤ q3) (h34 : q3 ≤ q4)
    (h45 : q4 ≤ q5) (h56 : q5 ≤ q6) (h6 : q6 < 1)
    (hsum : ∀ m : ℤ, θ ≤ |q0 + q1 + q2 + q3 + q4 + q5 + q6 - (m : ℝ)|) :
    ∃ φ : ℝ, ∀ x ∈ [q0, q1, q2, q3, q4, q5, q6], ∀ m : ℤ,
      1 / 14 + θ / 44 ≤ |φ - x - m| := by
  have hq00 : (0 : ℝ) ≤ q0 := le_of_eq h0.symm
  -- the spacing surplus: some circular gap is ≥ 1/7 + θ/22
  have hdisj : 1 / 7 + θ / 22 ≤ q1 - q0 ∨ 1 / 7 + θ / 22 ≤ q2 - q1 ∨ 1 / 7 + θ / 22 ≤ q3 - q2 ∨ 1 / 7 + θ / 22 ≤ q4 - q3 ∨ 1 / 7 + θ / 22 ≤ q5 - q4 ∨ 1 / 7 + θ / 22 ≤ q6 - q5 ∨ 1 / 7 + θ / 22 ≤ (1 : ℝ) - q6 := by
    by_contra hcon
    push_neg at hcon
    obtain ⟨hs0, hs1, hs2, hs3, hs4, hs5, hs6⟩ := hcon
    have h3 := hsum 3
    push_cast at h3
    rcases abs_cases (q0 + q1 + q2 + q3 + q4 + q5 + q6 - (3 : ℝ)) with ⟨heq, _⟩ | ⟨heq, _⟩
    · rw [heq] at h3
      linarith
    · rw [heq] at h3
      linarith
  rcases hdisj with hg | hg | hg | hg | hg | hg | hg
  · -- gap case 0
      refine ⟨(q0 + q1) / 2, ?_⟩
      intro x hx m
      simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
      have hgap := mid_gap_clearance q0 q1 (1 / 14 + θ / 44) (by linarith) (by linarith) (by linarith) (by linarith)
      rcases hx with rfl | rfl | rfl | rfl | rfl | rfl | rfl
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
  · -- gap case 1
      refine ⟨(q1 + q2) / 2, ?_⟩
      intro x hx m
      simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
      have hgap := mid_gap_clearance q1 q2 (1 / 14 + θ / 44) (by linarith) (by linarith) (by linarith) (by linarith)
      rcases hx with rfl | rfl | rfl | rfl | rfl | rfl | rfl
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
  · -- gap case 2
      refine ⟨(q2 + q3) / 2, ?_⟩
      intro x hx m
      simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
      have hgap := mid_gap_clearance q2 q3 (1 / 14 + θ / 44) (by linarith) (by linarith) (by linarith) (by linarith)
      rcases hx with rfl | rfl | rfl | rfl | rfl | rfl | rfl
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
  · -- gap case 3
      refine ⟨(q3 + q4) / 2, ?_⟩
      intro x hx m
      simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
      have hgap := mid_gap_clearance q3 q4 (1 / 14 + θ / 44) (by linarith) (by linarith) (by linarith) (by linarith)
      rcases hx with rfl | rfl | rfl | rfl | rfl | rfl | rfl
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
  · -- gap case 4
      refine ⟨(q4 + q5) / 2, ?_⟩
      intro x hx m
      simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
      have hgap := mid_gap_clearance q4 q5 (1 / 14 + θ / 44) (by linarith) (by linarith) (by linarith) (by linarith)
      rcases hx with rfl | rfl | rfl | rfl | rfl | rfl | rfl
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
  · -- gap case 5
      refine ⟨(q5 + q6) / 2, ?_⟩
      intro x hx m
      simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
      have hgap := mid_gap_clearance q5 q6 (1 / 14 + θ / 44) (by linarith) (by linarith) (by linarith) (by linarith)
      rcases hx with rfl | rfl | rfl | rfl | rfl | rfl | rfl
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inr (by linarith)) m
  · -- gap case 6
      refine ⟨(q6 + (1 : ℝ)) / 2, ?_⟩
      intro x hx m
      simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
      have hgap := mid_gap_clearance q6 (1 : ℝ) (1 / 14 + θ / 44) (by linarith) (by linarith) (by linarith) (by linarith)
      rcases hx with rfl | rfl | rfl | rfl | rfl | rfl | rfl
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m
      · exact hgap x (by linarith) (by linarith) (Or.inl (by linarith)) m

end SevenGap
end LonelyRunner
