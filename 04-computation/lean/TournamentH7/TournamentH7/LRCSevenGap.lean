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

set_option maxHeartbeats 1600000

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
    ∃ φ : ℝ, (0 ≤ φ ∧ φ ≤ 1) ∧ ∀ x ∈ [q0, q1, q2, q3, q4, q5, q6], ∀ m : ℤ,
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
      refine ⟨(q0 + q1) / 2, ⟨by linarith, by linarith⟩, ?_⟩
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
      refine ⟨(q1 + q2) / 2, ⟨by linarith, by linarith⟩, ?_⟩
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
      refine ⟨(q2 + q3) / 2, ⟨by linarith, by linarith⟩, ?_⟩
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
      refine ⟨(q3 + q4) / 2, ⟨by linarith, by linarith⟩, ?_⟩
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
      refine ⟨(q4 + q5) / 2, ⟨by linarith, by linarith⟩, ?_⟩
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
      refine ⟨(q5 + q6) / 2, ⟨by linarith, by linarith⟩, ?_⟩
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
      refine ⟨(q6 + (1 : ℝ)) / 2, ⟨by linarith, by linarith⟩, ?_⟩
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

/-! ## Stage B1 — the unsorted wrapper -/

/-- **Gap deficit for raw phase lists**: any 7 phases in `[0,1)` containing `0`
whose sum keeps integer distance `≥ θ` admit a probe clearing every phase by
`1/14 + θ/44`.  (Sorts with `insertionSort`, transfers pin/sum/membership along
the permutation, applies `seven_gap_deficit`.) -/
theorem gap_deficit_of_mem (L : List ℝ) (hlen : L.length = 7) (h0mem : (0 : ℝ) ∈ L)
    (hrange : ∀ x ∈ L, 0 ≤ x ∧ x < 1) (θ : ℝ) (hθ : 0 < θ)
    (hsum : ∀ m : ℤ, θ ≤ |L.sum - (m : ℝ)|) :
    ∃ φ : ℝ, (0 ≤ φ ∧ φ ≤ 1) ∧ ∀ x ∈ L, ∀ m : ℤ, 1 / 14 + θ / 44 ≤ |φ - x - m| := by
  classical
  set S : List ℝ := L.insertionSort (· ≤ ·) with hS
  have hperm : S.Perm L := List.perm_insertionSort _ L
  have hsort : List.Pairwise (· ≤ ·) S := List.pairwise_insertionSort _ L
  have hlenS : S.length = 7 := by rw [hperm.length_eq, hlen]
  -- destructure the sorted list
  rcases S with _ | ⟨q0, S1⟩
  · simp at hlenS
  rcases S1 with _ | ⟨q1, S2⟩
  · simp at hlenS
  rcases S2 with _ | ⟨q2, S3⟩
  · simp at hlenS
  rcases S3 with _ | ⟨q3, S4⟩
  · simp at hlenS
  rcases S4 with _ | ⟨q4, S5⟩
  · simp at hlenS
  rcases S5 with _ | ⟨q5, S6⟩
  · simp at hlenS
  rcases S6 with _ | ⟨q6, S7⟩
  · simp at hlenS
  rcases S7 with _ | ⟨q7, S8⟩
  swap
  · simp at hlenS
  -- sorted chain
  simp only [List.pairwise_cons, List.mem_cons, List.not_mem_nil, or_false,
    forall_eq_or_imp, forall_eq, List.Pairwise.nil, and_true] at hsort
  obtain ⟨⟨h01, h02, h03, h04, h05, h06⟩, ⟨h12, h13, h14, h15, h16⟩,
    ⟨h23, h24, h25, h26⟩, ⟨h34, h35, h36⟩, ⟨h45, h46⟩, h56⟩ := hsort
  -- the pin: 0 is a member and q0 is the sorted minimum; members are ≥ 0
  have hmemS : ∀ x ∈ L, x ∈ ([q0, q1, q2, q3, q4, q5, q6] : List ℝ) := fun x hx =>
    hperm.mem_iff.mpr hx
  have h0S : (0 : ℝ) ∈ ([q0, q1, q2, q3, q4, q5, q6] : List ℝ) := hmemS 0 h0mem
  have hq0L : q0 ∈ L := hperm.mem_iff.mp (by simp)
  have hq6L : q6 ∈ L := hperm.mem_iff.mp (by simp)
  have hq0nn : 0 ≤ q0 := (hrange q0 hq0L).1
  have hq0z : q0 = 0 := by
    simp only [List.mem_cons, List.not_mem_nil, or_false] at h0S
    rcases h0S with h | h | h | h | h | h | h <;> linarith
  have hq6lt : q6 < 1 := (hrange q6 hq6L).2
  -- the sum along the permutation
  have hsumS : q0 + q1 + q2 + q3 + q4 + q5 + q6 = L.sum := by
    rw [← hperm.sum_eq]
    simp only [List.sum_cons, List.sum_nil]
    ring
  obtain ⟨φ, hφb, hφ⟩ := seven_gap_deficit q0 q1 q2 q3 q4 q5 q6 θ hθ hq0z
    h01 h12 h23 h34 h45 h56.1 hq6lt (by rw [hsumS]; exact hsum)
  exact ⟨φ, hφb, fun x hx m => hφ x (hmemS x hx) m⟩

/-! ## Stage B2 — the frozen-drift sweep: a clustered 7-block owns a good point
in any window holding two `w1`-periods, given the transported sum-combo margin -/

/-- **THE CLUSTER SWEEP STEP**: seven ordered speeds within drift budget
(`1232·(w7−w1) ≤ w1`), a window `[τ, τ+L]` holding two base periods, and the
transported citation margin `1/14` on the SUM-COMBO `Σ(wⱼ−w1)` at `τ` — then some
`t'` in the window is `1/14`-good for ALL SEVEN runners.  The probe phase comes
from the gap-deficit lemma at the frozen offsets; sweeping `w1` to the probe
costs at most `2/w1` of time, and the offsets drift at most `1/616` — exactly
the gap-deficit surplus. -/
theorem cluster_sweep_step (w1 w2 w3 w4 w5 w6 w7 : ℤ) (hw1 : 0 < w1)
    (h12 : w1 ≤ w2) (h23 : w2 ≤ w3) (h34 : w3 ≤ w4) (h45 : w4 ≤ w5)
    (h56 : w5 ≤ w6) (h67 : w6 ≤ w7)
    (hclu : 1232 * (w7 - w1) ≤ w1)
    (τ L : ℝ) (hsweep : 2 ≤ (w1 : ℝ) * L)
    (hmargin : ∀ m : ℤ, (1 : ℝ) / 14 ≤ |((w2 + w3 + w4 + w5 + w6 + w7 - 6 * w1 : ℤ) : ℝ) * τ - m|) :
    ∃ t' : ℝ, τ ≤ t' ∧ t' ≤ τ + L ∧
      ∀ u ∈ ([w1, w2, w3, w4, w5, w6, w7] : List ℤ), ∀ m : ℤ,
        (1 : ℝ) / 14 ≤ |(u : ℝ) * t' - m| := by
  have hw1R : (0 : ℝ) < (w1 : ℝ) := by exact_mod_cast hw1
  have hz : Int.fract (-(((w1 - w1 : ℤ) : ℝ)) * τ) = 0 := by
    simp
  set P : List ℝ := [Int.fract (-(((w1 - w1 : ℤ) : ℝ)) * τ), Int.fract (-(((w2 - w1 : ℤ) : ℝ)) * τ), Int.fract (-(((w3 - w1 : ℤ) : ℝ)) * τ), Int.fract (-(((w4 - w1 : ℤ) : ℝ)) * τ), Int.fract (-(((w5 - w1 : ℤ) : ℝ)) * τ), Int.fract (-(((w6 - w1 : ℤ) : ℝ)) * τ), Int.fract (-(((w7 - w1 : ℤ) : ℝ)) * τ)] with hP
  have hlen : P.length = 7 := rfl
  have h0mem : (0 : ℝ) ∈ P := by
    rw [hP]
    simp only [List.mem_cons]
    exact Or.inl hz.symm
  have hrange : ∀ x ∈ P, 0 ≤ x ∧ x < 1 := by
    intro x hx
    rw [hP] at hx
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
    rcases hx with rfl | rfl | rfl | rfl | rfl | rfl | rfl <;>
      exact ⟨Int.fract_nonneg _, Int.fract_lt_one _⟩
  obtain ⟨K, hKdef⟩ : ∃ k : ℤ, k = ⌊-(((w2 - w1 : ℤ) : ℝ)) * τ⌋ + ⌊-(((w3 - w1 : ℤ) : ℝ)) * τ⌋ + ⌊-(((w4 - w1 : ℤ) : ℝ)) * τ⌋ + ⌊-(((w5 - w1 : ℤ) : ℝ)) * τ⌋ + ⌊-(((w6 - w1 : ℤ) : ℝ)) * τ⌋ + ⌊-(((w7 - w1 : ℤ) : ℝ)) * τ⌋ := ⟨_, rfl⟩
  have hPsum : P.sum = -(((w2 + w3 + w4 + w5 + w6 + w7 - 6 * w1 : ℤ) : ℝ) * τ) - (K : ℝ) := by
    rw [hP, hKdef]
    have hf2 := Int.floor_add_fract (-(((w2 - w1 : ℤ) : ℝ)) * τ)
    have hf3 := Int.floor_add_fract (-(((w3 - w1 : ℤ) : ℝ)) * τ)
    have hf4 := Int.floor_add_fract (-(((w4 - w1 : ℤ) : ℝ)) * τ)
    have hf5 := Int.floor_add_fract (-(((w5 - w1 : ℤ) : ℝ)) * τ)
    have hf6 := Int.floor_add_fract (-(((w6 - w1 : ℤ) : ℝ)) * τ)
    have hf7 := Int.floor_add_fract (-(((w7 - w1 : ℤ) : ℝ)) * τ)
    simp only [List.sum_cons, List.sum_nil]
    rw [hz]
    push_cast at *
    linarith
  have hsumMargin : ∀ m : ℤ, (1 : ℝ) / 14 ≤ |P.sum - (m : ℝ)| := by
    intro m
    rw [hPsum]
    have h := hmargin (-(K + m))
    have heq : -((((w2 + w3 + w4 + w5 + w6 + w7 - 6 * w1 : ℤ) : ℝ)) * τ) - (K : ℝ) - (m : ℝ)
        = -(((w2 + w3 + w4 + w5 + w6 + w7 - 6 * w1 : ℤ) : ℝ) * τ - ((-(K + m) : ℤ) : ℝ)) := by
      push_cast
      ring
    rw [heq, abs_neg]
    exact h
  obtain ⟨φ, ⟨hφ0, hφ1⟩, hφ⟩ := gap_deficit_of_mem P hlen h0mem hrange (1 / 14)
    (by norm_num) hsumMargin
  obtain ⟨t', ht'def⟩ : ∃ x : ℝ, x = ((⌊(w1 : ℝ) * τ⌋ : ℝ) + 1 + φ) / (w1 : ℝ) := ⟨_, rfl⟩
  have hw1t' : (w1 : ℝ) * t' = (⌊(w1 : ℝ) * τ⌋ : ℝ) + 1 + φ := by
    rw [ht'def, mul_div_cancel₀ _ (ne_of_gt hw1R)]
  have hfloor_le := Int.floor_le ((w1 : ℝ) * τ)
  have hlt_floor := Int.lt_floor_add_one ((w1 : ℝ) * τ)
  have hτt' : τ ≤ t' := by
    have h : (w1 : ℝ) * τ ≤ (w1 : ℝ) * t' := by
      rw [hw1t']
      linarith
    exact le_of_mul_le_mul_left h hw1R
  have ht'L : t' ≤ τ + L := by
    have h : (w1 : ℝ) * t' ≤ (w1 : ℝ) * (τ + L) := by
      have e : (w1 : ℝ) * (τ + L) = (w1 : ℝ) * τ + (w1 : ℝ) * L := by ring
      rw [e, hw1t']
      linarith
    exact le_of_mul_le_mul_left h hw1R
  have ht'τ0 : 0 ≤ t' - τ := by linarith
  have ht'τ2 : (w1 : ℝ) * (t' - τ) ≤ 2 := by
    have e : (w1 : ℝ) * (t' - τ) = (w1 : ℝ) * t' - (w1 : ℝ) * τ := by ring
    rw [e, hw1t']
    linarith
  have hcluR : (1232 : ℝ) * ((w7 : ℝ) - (w1 : ℝ)) ≤ (w1 : ℝ) := by exact_mod_cast hclu
  refine ⟨t', hτt', ht'L, ?_⟩
  intro u hu m
  simp only [List.mem_cons, List.not_mem_nil, or_false] at hu
  rcases hu with hu | hu | hu | hu | hu | hu | hu
  · rw [hu]
    have hfl := Int.floor_add_fract (-(((w1 - w1 : ℤ) : ℝ)) * τ)
    have hfr : Int.fract (-(((w1 - w1 : ℤ) : ℝ)) * τ) = -(((w1 - w1 : ℤ) : ℝ)) * τ - (⌊-(((w1 - w1 : ℤ) : ℝ)) * τ⌋ : ℝ) := by linarith
    have hpt := hφ (Int.fract (-(((w1 - w1 : ℤ) : ℝ)) * τ)) (by rw [hP]; simp) (m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w1 - w1 : ℤ) : ℝ)) * τ⌋)
    have hdposZ : (0 : ℤ) ≤ w1 - w1 := by linarith
    have hdpos : (0 : ℝ) ≤ ((w1 - w1 : ℤ) : ℝ) := by exact_mod_cast hdposZ
    have hj7Z : (w1 - w1 : ℤ) ≤ w7 - w1 := by linarith
    have hj7 : (((w1 - w1 : ℤ)) : ℝ) ≤ ((w7 - w1 : ℤ) : ℝ) := by exact_mod_cast hj7Z
    have hw7w1 : ((w7 - w1 : ℤ) : ℝ) = (w7 : ℝ) - (w1 : ℝ) := by push_cast; ring
    have hdrift : (((w1 - w1 : ℤ)) : ℝ) * (t' - τ) ≤ 1 / 616 := by
      have hA := mul_le_mul_of_nonneg_right hj7 ht'τ0
      have hB := mul_le_mul_of_nonneg_right hcluR ht'τ0
      rw [hw7w1] at hA
      nlinarith [ht'τ2]
    have habs_drift : |(((w1 - w1 : ℤ)) : ℝ) * (t' - τ)| ≤ 1 / 616 := by
      rw [abs_of_nonneg (mul_nonneg hdpos ht'τ0)]
      exact hdrift
    have heq : φ - Int.fract (-(((w1 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w1 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)
        = ((w1 : ℝ) * t' - m) - (((w1 - w1 : ℤ)) : ℝ) * (t' - τ) := by
      rw [hfr]
      push_cast
      linear_combination -hw1t'
    have htri : |φ - Int.fract (-(((w1 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w1 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)|
        ≤ |(w1 : ℝ) * t' - m| + |(((w1 - w1 : ℤ)) : ℝ) * (t' - τ)| := by
      rw [heq, sub_eq_add_neg]
      calc |((w1 : ℝ) * t' - m) + -((((w1 - w1 : ℤ)) : ℝ) * (t' - τ))|
          ≤ |(w1 : ℝ) * t' - m| + |(-((((w1 - w1 : ℤ)) : ℝ) * (t' - τ)))| := abs_add_le _ _
        _ = |(w1 : ℝ) * t' - m| + |(((w1 - w1 : ℤ)) : ℝ) * (t' - τ)| := by rw [abs_neg]
    linarith [hpt, htri, habs_drift]
  · rw [hu]
    have hfl := Int.floor_add_fract (-(((w2 - w1 : ℤ) : ℝ)) * τ)
    have hfr : Int.fract (-(((w2 - w1 : ℤ) : ℝ)) * τ) = -(((w2 - w1 : ℤ) : ℝ)) * τ - (⌊-(((w2 - w1 : ℤ) : ℝ)) * τ⌋ : ℝ) := by linarith
    have hpt := hφ (Int.fract (-(((w2 - w1 : ℤ) : ℝ)) * τ)) (by rw [hP]; simp) (m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w2 - w1 : ℤ) : ℝ)) * τ⌋)
    have hdposZ : (0 : ℤ) ≤ w2 - w1 := by linarith
    have hdpos : (0 : ℝ) ≤ ((w2 - w1 : ℤ) : ℝ) := by exact_mod_cast hdposZ
    have hj7Z : (w2 - w1 : ℤ) ≤ w7 - w1 := by linarith
    have hj7 : (((w2 - w1 : ℤ)) : ℝ) ≤ ((w7 - w1 : ℤ) : ℝ) := by exact_mod_cast hj7Z
    have hw7w1 : ((w7 - w1 : ℤ) : ℝ) = (w7 : ℝ) - (w1 : ℝ) := by push_cast; ring
    have hdrift : (((w2 - w1 : ℤ)) : ℝ) * (t' - τ) ≤ 1 / 616 := by
      have hA := mul_le_mul_of_nonneg_right hj7 ht'τ0
      have hB := mul_le_mul_of_nonneg_right hcluR ht'τ0
      rw [hw7w1] at hA
      nlinarith [ht'τ2]
    have habs_drift : |(((w2 - w1 : ℤ)) : ℝ) * (t' - τ)| ≤ 1 / 616 := by
      rw [abs_of_nonneg (mul_nonneg hdpos ht'τ0)]
      exact hdrift
    have heq : φ - Int.fract (-(((w2 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w2 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)
        = ((w2 : ℝ) * t' - m) - (((w2 - w1 : ℤ)) : ℝ) * (t' - τ) := by
      rw [hfr]
      push_cast
      linear_combination -hw1t'
    have htri : |φ - Int.fract (-(((w2 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w2 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)|
        ≤ |(w2 : ℝ) * t' - m| + |(((w2 - w1 : ℤ)) : ℝ) * (t' - τ)| := by
      rw [heq, sub_eq_add_neg]
      calc |((w2 : ℝ) * t' - m) + -((((w2 - w1 : ℤ)) : ℝ) * (t' - τ))|
          ≤ |(w2 : ℝ) * t' - m| + |(-((((w2 - w1 : ℤ)) : ℝ) * (t' - τ)))| := abs_add_le _ _
        _ = |(w2 : ℝ) * t' - m| + |(((w2 - w1 : ℤ)) : ℝ) * (t' - τ)| := by rw [abs_neg]
    linarith [hpt, htri, habs_drift]
  · rw [hu]
    have hfl := Int.floor_add_fract (-(((w3 - w1 : ℤ) : ℝ)) * τ)
    have hfr : Int.fract (-(((w3 - w1 : ℤ) : ℝ)) * τ) = -(((w3 - w1 : ℤ) : ℝ)) * τ - (⌊-(((w3 - w1 : ℤ) : ℝ)) * τ⌋ : ℝ) := by linarith
    have hpt := hφ (Int.fract (-(((w3 - w1 : ℤ) : ℝ)) * τ)) (by rw [hP]; simp) (m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w3 - w1 : ℤ) : ℝ)) * τ⌋)
    have hdposZ : (0 : ℤ) ≤ w3 - w1 := by linarith
    have hdpos : (0 : ℝ) ≤ ((w3 - w1 : ℤ) : ℝ) := by exact_mod_cast hdposZ
    have hj7Z : (w3 - w1 : ℤ) ≤ w7 - w1 := by linarith
    have hj7 : (((w3 - w1 : ℤ)) : ℝ) ≤ ((w7 - w1 : ℤ) : ℝ) := by exact_mod_cast hj7Z
    have hw7w1 : ((w7 - w1 : ℤ) : ℝ) = (w7 : ℝ) - (w1 : ℝ) := by push_cast; ring
    have hdrift : (((w3 - w1 : ℤ)) : ℝ) * (t' - τ) ≤ 1 / 616 := by
      have hA := mul_le_mul_of_nonneg_right hj7 ht'τ0
      have hB := mul_le_mul_of_nonneg_right hcluR ht'τ0
      rw [hw7w1] at hA
      nlinarith [ht'τ2]
    have habs_drift : |(((w3 - w1 : ℤ)) : ℝ) * (t' - τ)| ≤ 1 / 616 := by
      rw [abs_of_nonneg (mul_nonneg hdpos ht'τ0)]
      exact hdrift
    have heq : φ - Int.fract (-(((w3 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w3 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)
        = ((w3 : ℝ) * t' - m) - (((w3 - w1 : ℤ)) : ℝ) * (t' - τ) := by
      rw [hfr]
      push_cast
      linear_combination -hw1t'
    have htri : |φ - Int.fract (-(((w3 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w3 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)|
        ≤ |(w3 : ℝ) * t' - m| + |(((w3 - w1 : ℤ)) : ℝ) * (t' - τ)| := by
      rw [heq, sub_eq_add_neg]
      calc |((w3 : ℝ) * t' - m) + -((((w3 - w1 : ℤ)) : ℝ) * (t' - τ))|
          ≤ |(w3 : ℝ) * t' - m| + |(-((((w3 - w1 : ℤ)) : ℝ) * (t' - τ)))| := abs_add_le _ _
        _ = |(w3 : ℝ) * t' - m| + |(((w3 - w1 : ℤ)) : ℝ) * (t' - τ)| := by rw [abs_neg]
    linarith [hpt, htri, habs_drift]
  · rw [hu]
    have hfl := Int.floor_add_fract (-(((w4 - w1 : ℤ) : ℝ)) * τ)
    have hfr : Int.fract (-(((w4 - w1 : ℤ) : ℝ)) * τ) = -(((w4 - w1 : ℤ) : ℝ)) * τ - (⌊-(((w4 - w1 : ℤ) : ℝ)) * τ⌋ : ℝ) := by linarith
    have hpt := hφ (Int.fract (-(((w4 - w1 : ℤ) : ℝ)) * τ)) (by rw [hP]; simp) (m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w4 - w1 : ℤ) : ℝ)) * τ⌋)
    have hdposZ : (0 : ℤ) ≤ w4 - w1 := by linarith
    have hdpos : (0 : ℝ) ≤ ((w4 - w1 : ℤ) : ℝ) := by exact_mod_cast hdposZ
    have hj7Z : (w4 - w1 : ℤ) ≤ w7 - w1 := by linarith
    have hj7 : (((w4 - w1 : ℤ)) : ℝ) ≤ ((w7 - w1 : ℤ) : ℝ) := by exact_mod_cast hj7Z
    have hw7w1 : ((w7 - w1 : ℤ) : ℝ) = (w7 : ℝ) - (w1 : ℝ) := by push_cast; ring
    have hdrift : (((w4 - w1 : ℤ)) : ℝ) * (t' - τ) ≤ 1 / 616 := by
      have hA := mul_le_mul_of_nonneg_right hj7 ht'τ0
      have hB := mul_le_mul_of_nonneg_right hcluR ht'τ0
      rw [hw7w1] at hA
      nlinarith [ht'τ2]
    have habs_drift : |(((w4 - w1 : ℤ)) : ℝ) * (t' - τ)| ≤ 1 / 616 := by
      rw [abs_of_nonneg (mul_nonneg hdpos ht'τ0)]
      exact hdrift
    have heq : φ - Int.fract (-(((w4 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w4 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)
        = ((w4 : ℝ) * t' - m) - (((w4 - w1 : ℤ)) : ℝ) * (t' - τ) := by
      rw [hfr]
      push_cast
      linear_combination -hw1t'
    have htri : |φ - Int.fract (-(((w4 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w4 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)|
        ≤ |(w4 : ℝ) * t' - m| + |(((w4 - w1 : ℤ)) : ℝ) * (t' - τ)| := by
      rw [heq, sub_eq_add_neg]
      calc |((w4 : ℝ) * t' - m) + -((((w4 - w1 : ℤ)) : ℝ) * (t' - τ))|
          ≤ |(w4 : ℝ) * t' - m| + |(-((((w4 - w1 : ℤ)) : ℝ) * (t' - τ)))| := abs_add_le _ _
        _ = |(w4 : ℝ) * t' - m| + |(((w4 - w1 : ℤ)) : ℝ) * (t' - τ)| := by rw [abs_neg]
    linarith [hpt, htri, habs_drift]
  · rw [hu]
    have hfl := Int.floor_add_fract (-(((w5 - w1 : ℤ) : ℝ)) * τ)
    have hfr : Int.fract (-(((w5 - w1 : ℤ) : ℝ)) * τ) = -(((w5 - w1 : ℤ) : ℝ)) * τ - (⌊-(((w5 - w1 : ℤ) : ℝ)) * τ⌋ : ℝ) := by linarith
    have hpt := hφ (Int.fract (-(((w5 - w1 : ℤ) : ℝ)) * τ)) (by rw [hP]; simp) (m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w5 - w1 : ℤ) : ℝ)) * τ⌋)
    have hdposZ : (0 : ℤ) ≤ w5 - w1 := by linarith
    have hdpos : (0 : ℝ) ≤ ((w5 - w1 : ℤ) : ℝ) := by exact_mod_cast hdposZ
    have hj7Z : (w5 - w1 : ℤ) ≤ w7 - w1 := by linarith
    have hj7 : (((w5 - w1 : ℤ)) : ℝ) ≤ ((w7 - w1 : ℤ) : ℝ) := by exact_mod_cast hj7Z
    have hw7w1 : ((w7 - w1 : ℤ) : ℝ) = (w7 : ℝ) - (w1 : ℝ) := by push_cast; ring
    have hdrift : (((w5 - w1 : ℤ)) : ℝ) * (t' - τ) ≤ 1 / 616 := by
      have hA := mul_le_mul_of_nonneg_right hj7 ht'τ0
      have hB := mul_le_mul_of_nonneg_right hcluR ht'τ0
      rw [hw7w1] at hA
      nlinarith [ht'τ2]
    have habs_drift : |(((w5 - w1 : ℤ)) : ℝ) * (t' - τ)| ≤ 1 / 616 := by
      rw [abs_of_nonneg (mul_nonneg hdpos ht'τ0)]
      exact hdrift
    have heq : φ - Int.fract (-(((w5 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w5 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)
        = ((w5 : ℝ) * t' - m) - (((w5 - w1 : ℤ)) : ℝ) * (t' - τ) := by
      rw [hfr]
      push_cast
      linear_combination -hw1t'
    have htri : |φ - Int.fract (-(((w5 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w5 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)|
        ≤ |(w5 : ℝ) * t' - m| + |(((w5 - w1 : ℤ)) : ℝ) * (t' - τ)| := by
      rw [heq, sub_eq_add_neg]
      calc |((w5 : ℝ) * t' - m) + -((((w5 - w1 : ℤ)) : ℝ) * (t' - τ))|
          ≤ |(w5 : ℝ) * t' - m| + |(-((((w5 - w1 : ℤ)) : ℝ) * (t' - τ)))| := abs_add_le _ _
        _ = |(w5 : ℝ) * t' - m| + |(((w5 - w1 : ℤ)) : ℝ) * (t' - τ)| := by rw [abs_neg]
    linarith [hpt, htri, habs_drift]
  · rw [hu]
    have hfl := Int.floor_add_fract (-(((w6 - w1 : ℤ) : ℝ)) * τ)
    have hfr : Int.fract (-(((w6 - w1 : ℤ) : ℝ)) * τ) = -(((w6 - w1 : ℤ) : ℝ)) * τ - (⌊-(((w6 - w1 : ℤ) : ℝ)) * τ⌋ : ℝ) := by linarith
    have hpt := hφ (Int.fract (-(((w6 - w1 : ℤ) : ℝ)) * τ)) (by rw [hP]; simp) (m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w6 - w1 : ℤ) : ℝ)) * τ⌋)
    have hdposZ : (0 : ℤ) ≤ w6 - w1 := by linarith
    have hdpos : (0 : ℝ) ≤ ((w6 - w1 : ℤ) : ℝ) := by exact_mod_cast hdposZ
    have hj7Z : (w6 - w1 : ℤ) ≤ w7 - w1 := by linarith
    have hj7 : (((w6 - w1 : ℤ)) : ℝ) ≤ ((w7 - w1 : ℤ) : ℝ) := by exact_mod_cast hj7Z
    have hw7w1 : ((w7 - w1 : ℤ) : ℝ) = (w7 : ℝ) - (w1 : ℝ) := by push_cast; ring
    have hdrift : (((w6 - w1 : ℤ)) : ℝ) * (t' - τ) ≤ 1 / 616 := by
      have hA := mul_le_mul_of_nonneg_right hj7 ht'τ0
      have hB := mul_le_mul_of_nonneg_right hcluR ht'τ0
      rw [hw7w1] at hA
      nlinarith [ht'τ2]
    have habs_drift : |(((w6 - w1 : ℤ)) : ℝ) * (t' - τ)| ≤ 1 / 616 := by
      rw [abs_of_nonneg (mul_nonneg hdpos ht'τ0)]
      exact hdrift
    have heq : φ - Int.fract (-(((w6 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w6 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)
        = ((w6 : ℝ) * t' - m) - (((w6 - w1 : ℤ)) : ℝ) * (t' - τ) := by
      rw [hfr]
      push_cast
      linear_combination -hw1t'
    have htri : |φ - Int.fract (-(((w6 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w6 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)|
        ≤ |(w6 : ℝ) * t' - m| + |(((w6 - w1 : ℤ)) : ℝ) * (t' - τ)| := by
      rw [heq, sub_eq_add_neg]
      calc |((w6 : ℝ) * t' - m) + -((((w6 - w1 : ℤ)) : ℝ) * (t' - τ))|
          ≤ |(w6 : ℝ) * t' - m| + |(-((((w6 - w1 : ℤ)) : ℝ) * (t' - τ)))| := abs_add_le _ _
        _ = |(w6 : ℝ) * t' - m| + |(((w6 - w1 : ℤ)) : ℝ) * (t' - τ)| := by rw [abs_neg]
    linarith [hpt, htri, habs_drift]
  · rw [hu]
    have hfl := Int.floor_add_fract (-(((w7 - w1 : ℤ) : ℝ)) * τ)
    have hfr : Int.fract (-(((w7 - w1 : ℤ) : ℝ)) * τ) = -(((w7 - w1 : ℤ) : ℝ)) * τ - (⌊-(((w7 - w1 : ℤ) : ℝ)) * τ⌋ : ℝ) := by linarith
    have hpt := hφ (Int.fract (-(((w7 - w1 : ℤ) : ℝ)) * τ)) (by rw [hP]; simp) (m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w7 - w1 : ℤ) : ℝ)) * τ⌋)
    have hdposZ : (0 : ℤ) ≤ w7 - w1 := by linarith
    have hdpos : (0 : ℝ) ≤ ((w7 - w1 : ℤ) : ℝ) := by exact_mod_cast hdposZ
    have hj7Z : (w7 - w1 : ℤ) ≤ w7 - w1 := by linarith
    have hj7 : (((w7 - w1 : ℤ)) : ℝ) ≤ ((w7 - w1 : ℤ) : ℝ) := by exact_mod_cast hj7Z
    have hw7w1 : ((w7 - w1 : ℤ) : ℝ) = (w7 : ℝ) - (w1 : ℝ) := by push_cast; ring
    have hdrift : (((w7 - w1 : ℤ)) : ℝ) * (t' - τ) ≤ 1 / 616 := by
      have hA := mul_le_mul_of_nonneg_right hj7 ht'τ0
      have hB := mul_le_mul_of_nonneg_right hcluR ht'τ0
      rw [hw7w1] at hA
      nlinarith [ht'τ2]
    have habs_drift : |(((w7 - w1 : ℤ)) : ℝ) * (t' - τ)| ≤ 1 / 616 := by
      rw [abs_of_nonneg (mul_nonneg hdpos ht'τ0)]
      exact hdrift
    have heq : φ - Int.fract (-(((w7 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w7 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)
        = ((w7 : ℝ) * t' - m) - (((w7 - w1 : ℤ)) : ℝ) * (t' - τ) := by
      rw [hfr]
      push_cast
      linear_combination -hw1t'
    have htri : |φ - Int.fract (-(((w7 - w1 : ℤ) : ℝ)) * τ) - ((m - (⌊(w1 : ℝ) * τ⌋ + 1) + ⌊-(((w7 - w1 : ℤ) : ℝ)) * τ⌋ : ℤ) : ℝ)|
        ≤ |(w7 : ℝ) * t' - m| + |(((w7 - w1 : ℤ)) : ℝ) * (t' - τ)| := by
      rw [heq, sub_eq_add_neg]
      calc |((w7 : ℝ) * t' - m) + -((((w7 - w1 : ℤ)) : ℝ) * (t' - τ))|
          ≤ |(w7 : ℝ) * t' - m| + |(-((((w7 - w1 : ℤ)) : ℝ) * (t' - τ)))| := abs_add_le _ _
        _ = |(w7 : ℝ) * t' - m| + |(((w7 - w1 : ℤ)) : ℝ) * (t' - τ)| := by rw [abs_neg]
    linarith [hpt, htri, habs_drift]

/-! ## Stage C — the citation assembly: a clustered 7-far block is lonely -/

open LRC14

/-- The sum-combo of a 7-tuple: `Σⱼ(wⱼ − w₁) = w₂+…+w₇ − 6w₁`. -/
def sumCombo (w1 w2 w3 w4 w5 w6 w7 : ℤ) : ℤ := w2 + w3 + w4 + w5 + w6 + w7 - 6 * w1

/-- **THE CLUSTERED-BLOCK CLOSER**: a covering family of 13 nonzero speeds whose far
block is exactly seven *distinct* clustered integers `w₁ < … < w₇`
(`1232·(w₇−w₁) ≤ w₁`, so `w₁ ≥ 7392`) and whose remaining `k ≤ 11` runners are bounded
by `B` with `364·B ≤ w₁` — is lonely.

Cite the `k` bounded runners together with the sum-combo `C₀ = Σ(wⱼ−w₁)` (one extra
slot, `k+1 ≤ 12`); this pins `C₀` off the origin at `τ = t₀`, which by
`cluster_sweep_step` forces a good sweep point `t'` for all seven far runners in a
window of width `2/w₁`.  The bounded runners drift only `B·2/w₁ ≤ 1/182`, so their
citation margin `1/(k+2) ≥ 1/13` survives to exactly `1/14`. -/
theorem cite_cluster7_lonely (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (k : ℕ) (hk : k ≤ 11) (B : ℤ) (hB : 0 < B)
    (hcited : ∀ i : Fin 13, (i : ℕ) < k → |v i| ≤ B)
    (w1 w2 w3 w4 w5 w6 w7 : ℤ) (hw1 : 0 < w1)
    (h12 : w1 ≤ w2) (h23 : w2 ≤ w3) (h34 : w3 ≤ w4) (h45 : w4 ≤ w5)
    (h56 : w5 ≤ w6) (h67 : w6 ≤ w7) (hlt : w1 < w7)
    (hclu : 1232 * (w7 - w1) ≤ w1) (hbig : 364 * B ≤ w1)
    (hsplit : ∀ i : Fin 13, (i : ℕ) < k ∨
      |v i| ∈ ([w1, w2, w3, w4, w5, w6, w7] : List ℤ)) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hw1R : (0 : ℝ) < (w1 : ℝ) := by exact_mod_cast hw1
  have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB
  have hbigR : (364 : ℝ) * (B : ℝ) ≤ (w1 : ℝ) := by exact_mod_cast hbig
  set C0 : ℤ := sumCombo w1 w2 w3 w4 w5 w6 w7 with hC0
  have hC0pos : 0 < C0 := by
    rw [hC0, sumCombo]; omega
  -- the citation tuple: k bounded runners, then C0
  have hk13 : k ≤ 13 := by omega
  set wc : Fin (k + 1) → ℤ :=
    Fin.snoc (fun j : Fin k => v (Fin.castLE hk13 j)) C0 with hwc
  have hwcne : ∀ i, wc i ≠ 0 := by
    intro i
    refine Fin.lastCases ?_ ?_ i
    · rw [hwc, Fin.snoc_last]; omega
    · intro j
      rw [hwc, Fin.snoc_castSucc]; exact hv _
  obtain ⟨t₀, hcite⟩ := cite (k + 1) (by omega) wc hwcne
  -- the sweep window
  set L : ℝ := 2 / (w1 : ℝ) with hL
  have hsweep : 2 ≤ (w1 : ℝ) * L := by
    rw [hL, mul_div_cancel₀ _ (ne_of_gt hw1R)]
  -- margin on C0 at τ = t₀ (last citation entry), weakened to 1/14
  have hcombo : ∀ m : ℤ, (1 : ℝ) / 14 ≤ |((sumCombo w1 w2 w3 w4 w5 w6 w7 : ℤ) : ℝ) * t₀ - m| := by
    intro m
    have hlast := hcite (Fin.last k) m
    rw [hwc, Fin.snoc_last] at hlast
    have hkk : (1 : ℝ) / 14 ≤ 1 / ((k + 1 : ℕ) + 1 : ℕ) := by
      apply one_div_le_one_div_of_le (by positivity)
      have : ((k + 1 : ℕ) + 1 : ℕ) ≤ 14 := by omega
      exact_mod_cast this
    rw [← hC0]
    calc (1 : ℝ) / 14 ≤ 1 / ((k + 1 : ℕ) + 1 : ℕ) := hkk
      _ ≤ |(C0 : ℝ) * t₀ - m| := hlast
  obtain ⟨t', ht'1, ht'2, hgood⟩ := cluster_sweep_step w1 w2 w3 w4 w5 w6 w7 hw1
    h12 h23 h34 h45 h56 h67 hclu t₀ L hsweep hcombo
  refine ⟨t', ?_⟩
  intro i m
  rcases hsplit i with hnear | hfar
  · -- bounded runner: citation margin 1/(k+2) minus drift B·L ≥ 1/14
    set j : Fin k := ⟨(i : ℕ), hnear⟩ with hj
    have hidx : Fin.castLE hk13 j = i := by
      apply Fin.ext; rfl
    have hcj := hcite (Fin.castSucc j) m
    rw [hwc, Fin.snoc_castSucc, hidx] at hcj
    -- hcj : 1/(k+2) ≤ |v i · t₀ - m|
    have hvB : |(v i : ℝ)| ≤ (B : ℝ) := by
      rw [← Int.cast_abs]; exact_mod_cast hcited i hnear
    have htri : |(v i : ℝ) * t₀ - m| ≤ |(v i : ℝ) * t' - m| + |(v i : ℝ)| * |t₀ - t'| := by
      calc |(v i : ℝ) * t₀ - m|
          = |((v i : ℝ) * t' - m) + (v i : ℝ) * (t₀ - t')| := by congr 1; ring
        _ ≤ |(v i : ℝ) * t' - m| + |(v i : ℝ) * (t₀ - t')| := abs_add_le _ _
        _ = |(v i : ℝ) * t' - m| + |(v i : ℝ)| * |t₀ - t'| := by rw [abs_mul]
    have hdt : |t₀ - t'| ≤ L := by
      rw [abs_le]
      constructor <;> [linarith; linarith]
    have hLpos : 0 ≤ L := by positivity
    have hlip : |(v i : ℝ)| * |t₀ - t'| ≤ (B : ℝ) * L := by
      apply mul_le_mul hvB hdt (abs_nonneg _) hBR.le
    have hBL : (B : ℝ) * L ≤ 1 / 182 := by
      rw [hL]
      have he : (B : ℝ) * (2 / w1) = 2 * B / w1 := by ring
      rw [he, div_le_iff₀ hw1R]
      nlinarith [hbigR]
    have hmarg : (1 : ℝ) / ((k + 1 : ℕ) + 1 : ℕ) ≤ |(v i : ℝ) * t₀ - m| := hcj
    have h13 : (1 : ℝ) / 13 ≤ 1 / ((k + 1 : ℕ) + 1 : ℕ) := by
      apply one_div_le_one_div_of_le (by positivity)
      have : ((k + 1 : ℕ) + 1 : ℕ) ≤ 13 := by omega
      exact_mod_cast this
    -- 1/13 - 1/182 = 1/14
    have hfinal : (1 : ℝ) / 14 ≤ |(v i : ℝ) * t' - m| := by
      have hchain : (1 : ℝ) / 13 ≤ |(v i : ℝ) * t' - m| + 1 / 182 := by
        calc (1 : ℝ) / 13 ≤ 1 / ((k + 1 : ℕ) + 1 : ℕ) := h13
          _ ≤ |(v i : ℝ) * t₀ - m| := hmarg
          _ ≤ |(v i : ℝ) * t' - m| + |(v i : ℝ)| * |t₀ - t'| := htri
          _ ≤ |(v i : ℝ) * t' - m| + (B : ℝ) * L := by linarith [hlip]
          _ ≤ |(v i : ℝ) * t' - m| + 1 / 182 := by linarith [hBL]
      linarith
    exact hfinal
  · -- far runner: sign bridge onto the seven-block good point
    have hgd := hgood |v i| hfar
    rcases abs_cases ((v i : ℝ)) with ⟨habs, _⟩ | ⟨habs, _⟩
    · have hcast : ((|v i| : ℤ) : ℝ) = (v i : ℝ) := by rw [Int.cast_abs, habs]
      have := hgd m
      rw [hcast] at this
      exact this
    · have hcast : ((|v i| : ℤ) : ℝ) = -(v i : ℝ) := by rw [Int.cast_abs, habs]
      have := hgd (-m)
      rw [hcast] at this
      have heq : |-(v i : ℝ) * t' - ((-m : ℤ) : ℝ)| = |(v i : ℝ) * t' - m| := by
        rw [show -(v i : ℝ) * t' - ((-m : ℤ) : ℝ) = -((v i : ℝ) * t' - m) by push_cast; ring,
          abs_neg]
      rw [heq] at this
      exact this

end SevenGap
end LonelyRunner
