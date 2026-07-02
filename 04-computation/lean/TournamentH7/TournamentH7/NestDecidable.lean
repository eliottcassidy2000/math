/-
klein-2026-07-01-S93 (HYP-4000) — the Lean-ready restatement of THM-601's verification.

THE RESTATEMENT. A non-origin point common to the three dangers D_u, D_v, D_w at
r = 1/14 is a point x with |x - a/u| < r/u, |x - b/v| < r/v, |x - c/w| < r/w for integer
indices (a,b,c) outside the two origin classes. By 1D HELLY (three open intervals share
a point iff they pairwise intersect — `helly_three_Ioo` below), existence is equivalent
to the three PAIRWISE integer inequalities

    14·|a·v − b·u| ≤ u+v−1,   14·|b·w − c·v| ≤ v+w−1,   14·|a·w − c·u| ≤ u+w−1,

a BOUNDED INTEGER statement with no measure theory. THM-601's 8100-case verification
reduces (by monotonicity of intersections) to: no gcd-reduced triple of {1..13} admits a
non-origin solution — `cap_universe_no_coincidence` below, discharged by `decide`
(cross-validated against exact interval arithmetic on all 1140 triples ≤ 20 in
helly_integer_criterion_klein.py: zero mismatches).
-/
import Mathlib

namespace NestDecidable

/-- The pairwise closeness condition at r = 1/14: the index-(a,b) intervals of speeds
u, v overlap iff `14|av − bu| ≤ u + v − 1`. -/
def pairOK (u v a b : ℕ) : Prop :=
  14 * ((a * v : ℤ) - (b * u : ℤ)).natAbs ≤ u + v - 1

instance (u v a b : ℕ) : Decidable (pairOK u v a b) := by
  unfold pairOK; infer_instance

/-- A non-origin triple coincidence for speeds u < v < w: bounded indices, both origin
classes excluded, all three pairwise conditions. -/
def hasCoincidence (u v w : ℕ) : Prop :=
  ∃ a ≤ u, ∃ b ≤ v, ∃ c ≤ w,
    (¬(a = 0 ∧ b = 0 ∧ c = 0)) ∧ (¬(a = u ∧ b = v ∧ c = w)) ∧
    pairOK u v a b ∧ pairOK v w b c ∧ pairOK u w a c

instance (u v w : ℕ) : Decidable (hasCoincidence u v w) := by
  unfold hasCoincidence; infer_instance

/-- **The cap-universe check (THM-601's enumeration, restated):** no gcd-reduced triple
of `{1..13}` admits a non-origin coincidence. 286 bounded-integer instances (compiled
evaluation; cross-validated against exact interval arithmetic, zero mismatches). -/
theorem cap_universe_no_coincidence :
    ∀ w ≤ 13, ∀ v < w, ∀ u < v, 0 < u →
      ¬ hasCoincidence (u / ((u.gcd v).gcd w)) (v / ((u.gcd v).gcd w))
        (w / ((u.gcd v).gcd w)) := by
  native_decide

/-- **1D Helly for three open intervals** (densely ordered): pairwise intersecting open
intervals have a common point. This is the bridge that turns the triple-danger
non-existence into the three pairwise integer conditions. -/
theorem helly_three_Ioo {α : Type*} [LinearOrder α] [DenselyOrdered α]
    {a₁ b₁ a₂ b₂ a₃ b₃ : α}
    (h₁₂ : (Set.Ioo a₁ b₁ ∩ Set.Ioo a₂ b₂).Nonempty)
    (h₂₃ : (Set.Ioo a₂ b₂ ∩ Set.Ioo a₃ b₃).Nonempty)
    (h₁₃ : (Set.Ioo a₁ b₁ ∩ Set.Ioo a₃ b₃).Nonempty) :
    (Set.Ioo a₁ b₁ ∩ Set.Ioo a₂ b₂ ∩ Set.Ioo a₃ b₃).Nonempty := by
  obtain ⟨x, ⟨hx1, hx1'⟩, hx2, hx2'⟩ := h₁₂
  obtain ⟨y, ⟨hy2, hy2'⟩, hy3, hy3'⟩ := h₂₃
  obtain ⟨z, ⟨hz1, hz1'⟩, hz3, hz3'⟩ := h₁₃
  have h11 : a₁ < b₁ := hx1.trans hx1'
  have h12 : a₁ < b₂ := hx1.trans hx2'
  have h13 : a₁ < b₃ := hz1.trans hz3'
  have h21 : a₂ < b₁ := hx2.trans hx1'
  have h22 : a₂ < b₂ := hx2.trans hx2'
  have h23 : a₂ < b₃ := hy2.trans hy3'
  have h31 : a₃ < b₁ := hz3.trans hz1'
  have h32 : a₃ < b₂ := hy3.trans hy2'
  have h33 : a₃ < b₃ := hy3.trans hy3'
  have hLR : max a₁ (max a₂ a₃) < min b₁ (min b₂ b₃) := by
    rw [max_lt_iff, max_lt_iff]
    refine ⟨?_, ?_, ?_⟩ <;> rw [lt_min_iff, lt_min_iff]
    · exact ⟨h11, h12, h13⟩
    · exact ⟨h21, h22, h23⟩
    · exact ⟨h31, h32, h33⟩
  obtain ⟨t, htL, htR⟩ := exists_between hLR
  refine ⟨t, ⟨⟨?_, ?_⟩, ?_, ?_⟩, ?_, ?_⟩
  · exact lt_of_le_of_lt (le_max_left _ _) htL
  · exact lt_of_lt_of_le htR (min_le_left _ _)
  · exact lt_of_le_of_lt ((le_max_left a₂ a₃).trans (le_max_right a₁ _)) htL
  · exact lt_of_lt_of_le htR ((min_le_right b₁ _).trans (min_le_left b₂ b₃))
  · exact lt_of_le_of_lt ((le_max_right a₂ a₃).trans (le_max_right a₁ _)) htL
  · exact lt_of_lt_of_le htR ((min_le_right b₁ _).trans (min_le_right b₂ b₃))

end NestDecidable
