/-
  TournamentH7.LRCPinch — the grid-invisible pinch abscissae (kind-pasteur-2026-07-09-S107).

  A "pinch" of the LRC good set is a tooth collision `frac(a·x) = frac(b·x)` — two teeth coincide, the
  local normal form of the lemniscate node (two crossing branches).  This file pins the exact fact that
  the pinches sit at RATIONAL abscissae `x = m/(a−b)` (`d = a−b` a cluster difference), hence at a
  FINITE, a-priori set whose denominators are the cluster differences — NOT the ruler `Vmax` — which is
  why they are invisible to every uniform ruler grid `x_j = (j+½)/Vmax` (mac-mini-S64 MISTAKE-130).
  Self-contained (imports Mathlib).
-/
import Mathlib

namespace LRC14Pinch

/-- **The pinch characterization.**  Two teeth `frac(a·x)`, `frac(b·x)` collide exactly when
`(a−b)·x` is an integer — i.e. `x = m/(a−b)` for some `m ∈ ℤ`.  The collision abscissae are the
cluster-difference rationals; for `a ≠ b` and `x ∈ [0,1)` there are exactly `|a−b|` of them. -/
theorem collision_iff (a b : ℤ) (x : ℝ) :
    Int.fract ((a : ℝ) * x) = Int.fract ((b : ℝ) * x) ↔ ∃ m : ℤ, ((a - b : ℤ) : ℝ) * x = (m : ℝ) := by
  rw [Int.fract_eq_fract]
  constructor
  · rintro ⟨z, hz⟩
    exact ⟨z, by push_cast; linear_combination hz⟩
  · rintro ⟨m, hm⟩
    exact ⟨m, by push_cast at hm ⊢; linarith⟩

/-- **A pinch is a rational with denominator a cluster difference.**  If teeth `a, b` (`a ≠ b`) collide
at `x`, then `x = m/(a−b)` for some `m ∈ ℤ` — a grid-invisible rational (denominator `= a−b`, not
`Vmax`). -/
theorem pinch_rational (a b : ℤ) (hab : a ≠ b) (x : ℝ)
    (h : Int.fract ((a : ℝ) * x) = Int.fract ((b : ℝ) * x)) :
    ∃ m : ℤ, x = (m : ℝ) / ((a - b : ℤ) : ℝ) := by
  obtain ⟨m, hm⟩ := (collision_iff a b x).mp h
  have hd : ((a - b : ℤ) : ℝ) ≠ 0 := by
    have : a - b ≠ 0 := sub_ne_zero.mpr hab
    exact_mod_cast this
  exact ⟨m, by field_simp; linarith [hm]⟩

end LRC14Pinch
