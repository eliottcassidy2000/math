/-
  TournamentH7.IntervalEscape  (mac-mini-2026-07-04-S50)

  The loose branch's geometric kernel (THM-619(i)) + a consolidation.

  PRIMITIVE `interval_escape`: any interval of length > 2r (0 < r ≤ 1/2)
  contains a point at distance ≥ r from every integer.  Proof geometry: with
  m = ⌊A + 1/2⌋, either the half-integer m + 1/2 lies in [A,B] (distance to
  ℤ is 1/2 ≥ r), or [A,B] sits inside the unit cell (m − 1/2, m + 1/2), where
  length > 2r forces a poke-out past m ± r; the poked endpoint clears m by the
  poke and every other integer by the triangle inequality (≥ 1 − 1/2 = 1/2).

  COROLLARY `no_cover_of_long` (THM-619(i) in Lean): a single radius-r comb of
  speed w cannot cover any interval J with w·|J| > 2r — one-tooth containment
  is the only way.  (The S16 antipode dodge is the |J| = 1/w instance, with
  constant r in place of its sharper 1/2 — consolidation note in the log.)

  Sorry-free.  All-ℚ.
-/
import Mathlib.Tactic

namespace TournamentH7.IntervalEscape

/-- distance to every integer other than `m` is large when `x` is within the
unit cell of `m`. -/
private lemma far_from_others {x : ℚ} {m k : ℤ} (hk : k ≠ m)
    (hcell : |x - m| ≤ 1 / 2) : (1 : ℚ) / 2 ≤ |x - k| := by
  have h0 : (m - k : ℤ) ≠ 0 := sub_ne_zero.mpr (Ne.symm hk)
  have h1 : (1 : ℚ) ≤ |(m : ℚ) - k| := by
    have hz : (1 : ℤ) ≤ |m - k| := Int.one_le_abs h0
    calc (1 : ℚ) = ((1 : ℤ) : ℚ) := by norm_num
      _ ≤ ((|m - k| : ℤ) : ℚ) := by exact_mod_cast hz
      _ = |(m : ℚ) - k| := by rw [Int.cast_abs]; push_cast; ring_nf
  have htri : |(m : ℚ) - k| ≤ |x - m| + |x - k| := by
    calc |(m : ℚ) - k| = |((m : ℚ) - x) + (x - k)| := by ring_nf
      _ ≤ |(m : ℚ) - x| + |x - k| := abs_add_le _ _
      _ = |x - m| + |x - k| := by rw [abs_sub_comm]
  linarith

/-- **The escape primitive.** -/
theorem interval_escape (r A B : ℚ) (hr : 0 < r) (hr2 : r ≤ 1 / 2)
    (hlen : 2 * r < B - A) :
    ∃ x : ℚ, A ≤ x ∧ x ≤ B ∧ ∀ k : ℤ, r ≤ |x - k| := by
  set m : ℤ := ⌊A + 1 / 2⌋ with hm
  have hm1 : (m : ℚ) ≤ A + 1 / 2 := Int.floor_le _
  have hm2 : A + 1 / 2 < (m : ℚ) + 1 := Int.lt_floor_add_one _
  by_cases hin : (m : ℚ) + 1 / 2 ≤ B
  · -- the half-integer m + 1/2 ∈ [A, B]
    refine ⟨(m : ℚ) + 1 / 2, by linarith, hin, ?_⟩
    intro k
    have h1 : ((m : ℚ) + 1 / 2) - k = ((m - k : ℤ) : ℚ) + 1 / 2 := by push_cast; ring
    rw [h1]
    by_cases h : (0 : ℚ) ≤ ((m - k : ℤ) : ℚ)
    · rw [abs_of_nonneg (by linarith)]; linarith
    · push_neg at h
      have hle : ((m - k : ℤ) : ℚ) ≤ -1 := by
        have h0 : (m - k : ℤ) < 0 := by exact_mod_cast h
        have h1' : (m - k : ℤ) ≤ -1 := by omega
        exact_mod_cast h1'
      rw [abs_of_nonpos (by linarith)]; linarith
  · -- whole interval inside the unit cell (m − 1/2, m + 1/2)
    push_neg at hin
    have hAcell : (m : ℚ) - 1 / 2 ≤ A := by linarith
    have hBcell : B < (m : ℚ) + 1 / 2 := hin
    -- length > 2r forces a poke-out past m − r or m + r
    have hpoke : A < (m : ℚ) - r ∨ (m : ℚ) + r < B := by
      by_contra hcon
      push_neg at hcon
      obtain ⟨h1, h2⟩ := hcon
      linarith
    rcases hpoke with hA | hB
    · -- x = A: clears m by the poke; others by the cell
      refine ⟨A, le_refl _, by linarith, ?_⟩
      intro k
      have hcell : |A - (m : ℚ)| ≤ 1 / 2 := by
        rw [abs_le]; constructor <;> linarith
      rcases eq_or_ne k m with rfl | hk
      · have hAm : A - (m : ℚ) ≤ 0 := by linarith
        rw [abs_of_nonpos hAm]
        linarith
      · have hf := far_from_others hk hcell
        linarith
    · -- x = B symmetric
      refine ⟨B, by linarith, le_refl _, ?_⟩
      intro k
      have hcell : |B - (m : ℚ)| ≤ 1 / 2 := by
        rw [abs_le]; constructor <;> linarith
      rcases eq_or_ne k m with rfl | hk
      · have hBm : (0 : ℚ) ≤ B - (m : ℚ) := by linarith
        rw [abs_of_nonneg hBm]
        linarith
      · have hf := far_from_others hk hcell
        linarith

/-- **THM-619(i) in Lean: one-tooth containment.**  A single radius-`r` comb of
speed `w` cannot cover an interval `J = [a,b]` with `w·(b − a) > 2r`: some
`t ∈ J` has `|w·t − k| ≥ r` for every integer `k`. -/
theorem no_cover_of_long (w : ℕ) (hw : 0 < w) (r a b : ℚ) (hr : 0 < r)
    (hr2 : r ≤ 1 / 2) (hlen : 2 * r < (w : ℚ) * (b - a)) :
    ∃ t : ℚ, a ≤ t ∧ t ≤ b ∧ ∀ k : ℤ, r ≤ |(w : ℚ) * t - k| := by
  have hw' : (0 : ℚ) < (w : ℚ) := by exact_mod_cast hw
  obtain ⟨x, hx1, hx2, hx3⟩ := interval_escape r ((w : ℚ) * a) ((w : ℚ) * b) hr hr2
    (by nlinarith)
  refine ⟨x / w, ?_, ?_, ?_⟩
  · rw [le_div_iff₀ hw']; linarith
  · rw [div_le_iff₀ hw']; linarith
  · intro k
    have hval : (w : ℚ) * (x / w) = x := by field_simp
    rw [hval]
    exact hx3 k

end TournamentH7.IntervalEscape
