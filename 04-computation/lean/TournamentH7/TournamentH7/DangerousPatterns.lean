/-
  TournamentH7.DangerousPatterns  (mac-mini-2026-07-01-S101)

  THM-601(i), the constructive half, Lean-checked: the dangerous-pattern
  characterization of the free-phase pair layer (THM-598/599/601).

  For a coprime pattern (P, Q) and radius r, if the phase θ satisfies
  ‖P·θ‖ ≥ r(P+Q) (distance to ℤ), then NO time x is double-covered:
  ¬(‖Px‖ < r ∧ ‖Qx − θ‖ < r).  Reason: with a, b the near integers,
  Q(Px − a) − P(Qx − θ − b) = Pθ − (Qa − Pb), so dist(Pθ, ℤ) < r(P+Q).
  Consequently, whenever r(P+Q) ≤ 1/2 the witness θ = 1/(2P) gives a phase
  with EMPTY pattern overlap — the "dangerous" (zero-min) regime of THM-601:
  at r = 1/14 exactly the nine coprime patterns with P + Q ≤ 7.

  This eliminates the Fourier layer of THM-598(A) from the formal chain:
  the avoidance direction is pure arithmetic.  Sorry-free.
-/
import Mathlib.Tactic

namespace TournamentH7.DangerousPatterns

/-- **Avoidance from a protected phase.**  If every integer is at distance
`≥ r*(P+Q)` from `P*θ`, then no `x` has both `|P*x − a| < r` and
`|Q*x − θ − b| < r` for integers `a, b`. -/
theorem no_double_cover (P Q : ℕ) (r θ : ℝ) (hP : 0 < P) (hQ : 0 < Q)
    (hθ : ∀ m : ℤ, r * (P + Q) ≤ |P * θ - m|) :
    ∀ x : ℝ, ¬ ((∃ a : ℤ, |P * x - a| < r) ∧ (∃ b : ℤ, |Q * x - θ - b| < r)) := by
  rintro x ⟨⟨a, ha⟩, ⟨b, hb⟩⟩
  have hPr : (0 : ℝ) < P := by exact_mod_cast hP
  have hQr : (0 : ℝ) < Q := by exact_mod_cast hQ
  -- the key algebraic identity
  have key : (P : ℝ) * θ - (Q * a - P * b : ℤ)
      = (P : ℝ) * (Q * x - θ - b) * (-1) + (Q : ℝ) * (P * x - a) := by
    push_cast
    ring
  have hcombo : |(P : ℝ) * θ - (Q * a - P * b : ℤ)| < r * (P + Q) := by
    rw [key]
    calc |(P : ℝ) * (Q * x - θ - b) * (-1) + (Q : ℝ) * (P * x - a)|
        ≤ |(P : ℝ) * (Q * x - θ - b) * (-1)| + |(Q : ℝ) * (P * x - a)| := abs_add_le _ _
      _ = (P : ℝ) * |Q * x - θ - b| + (Q : ℝ) * |P * x - a| := by
          rw [abs_mul, abs_mul, abs_mul, abs_neg, abs_one, mul_one,
            abs_of_pos hPr, abs_of_pos hQr]
      _ < (P : ℝ) * r + (Q : ℝ) * r := by
          apply add_lt_add
          · exact mul_lt_mul_of_pos_left hb hPr
          · exact mul_lt_mul_of_pos_left ha hQr
      _ = r * (P + Q) := by ring
  exact absurd hcombo (not_lt.mpr (hθ (Q * a - P * b)))

/-- **The dangerous-phase witness.**  If `r*(P+Q) ≤ 1/2`, the phase
`θ = 1/(2P)` satisfies the protection hypothesis: `P*θ = 1/2`, whose distance
to every integer is `1/2 ≥ r(P+Q)`. -/
theorem protected_phase_exists (P Q : ℕ) (r : ℝ) (hP : 0 < P)
    (hr : r * (P + Q) ≤ 1 / 2) :
    ∀ m : ℤ, r * (P + Q) ≤ |(P : ℝ) * (1 / (2 * P)) - m| := by
  intro m
  have hPr : (0 : ℝ) < P := by exact_mod_cast hP
  have hval : (P : ℝ) * (1 / (2 * P)) = 1 / 2 := by
    field_simp
  rw [hval]
  -- distance from 1/2 to any integer is at least 1/2
  have : (1 : ℝ) / 2 ≤ |1 / 2 - (m : ℝ)| := by
    by_cases h : (m : ℝ) ≤ 0
    · rw [abs_of_nonneg (by linarith)]
      linarith
    · push_neg at h
      have hm1 : (1 : ℝ) ≤ (m : ℝ) := by
        have h0 : (0 : ℤ) < m := by exact_mod_cast h
        exact_mod_cast h0
      rw [abs_of_nonpos (by linarith)]
      linarith
  linarith

/-- **THM-601(i), avoidance direction, packaged.**  For any coprime-or-not
pattern `(P, Q)` with `r*(P+Q) ≤ 1/2` there is a phase whose pattern overlap
is EMPTY: the pattern is dangerous.  (At `r = 1/14`: `P + Q ≤ 7`.) -/
theorem dangerous_of_small_sum (P Q : ℕ) (r : ℝ) (hP : 0 < P) (hQ : 0 < Q)
    (hr : r * (P + Q) ≤ 1 / 2) :
    ∃ θ : ℝ, ∀ x : ℝ,
      ¬ ((∃ a : ℤ, |P * x - a| < r) ∧ (∃ b : ℤ, |Q * x - θ - b| < r)) :=
  ⟨1 / (2 * P), no_double_cover P Q r _ hP hQ (protected_phase_exists P Q r hP hr)⟩

end TournamentH7.DangerousPatterns
