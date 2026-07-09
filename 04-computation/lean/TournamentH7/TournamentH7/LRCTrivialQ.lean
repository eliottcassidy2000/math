/-
  TournamentH7.LRCTrivialQ — the trivial-`q` loneliness lemma (mac-mini-2026-07-09-S64).

  For any `q` with `0 < q ≤ 14`, the rational time `t = 1/q` is lonely for `v` as soon as `q` divides
  NO speed.  Reason: `|v_i/q − m| = |v_i − m·q| / q`, and `v_i − m·q` is a NONZERO integer (that is
  exactly `q ∤ v_i`), so its absolute value is `≥ 1`; hence `|v_i/q − m| ≥ 1/q ≥ 1/14`.

  COROLLARY (`lonely_of_not_dvd_14`): if `14` divides no speed, LRC(14) holds outright at `t = 1/14`.
  So a HARD instance must have, for EVERY `q ∈ {2,…,14}`, some speed divisible by `q` — a covering
  condition on the speed set.  (This is why the equality extremal `{1,…,13}` is settled by `t = 1/14`:
  no speed there is a multiple of `14`, and the margin is exactly `1/14`.)

  Companion to kps-S28's `spread13_lonely` (which closes every ratio `Vmax ≤ 13·Vmin`): together they
  strip away everything except speed sets that are simultaneously wide (`Vmax > 13·Vmin`) and
  small-modulus-covering.  Self-contained on `LonelyRunner`.
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LRC14

open LonelyRunner

/-- **The trivial-`q` loneliness lemma.**  If `0 < q ≤ 14` and `q` divides none of the speeds, then
every runner is at least `1/14` from the origin at time `t = 1/q`.  (Margin is in fact `≥ 1/q`.) -/
theorem lonely_of_not_dvd {ι : Type*} (v : ι → ℤ) (q : ℤ)
    (hq : 0 < q) (hq14 : q ≤ 14) (hnd : ∀ i, ¬ (q ∣ v i)) :
    Lonely 14 v (1 / (q : ℝ)) := by
  intro i m
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hq0 : (q : ℝ) ≠ 0 := ne_of_gt hqR
  -- `v i * (1/q) − m = (v i − m·q)/q`
  have key : (v i : ℝ) * (1 / (q : ℝ)) - (m : ℝ) = (((v i - m * q : ℤ)) : ℝ) / (q : ℝ) := by
    push_cast
    field_simp
  -- `v i − m·q ≠ 0`, i.e. it is a nonzero integer, so `|v i − m·q| ≥ 1`
  have hne : (v i - m * q) ≠ 0 := by
    intro h
    exact hnd i ⟨m, by rw [mul_comm]; linarith⟩
  have h0 : 0 < |v i - m * q| := abs_pos.mpr hne
  have h1 : (1 : ℤ) ≤ |v i - m * q| := by omega
  have h1R : (1 : ℝ) ≤ |(((v i - m * q : ℤ)) : ℝ)| := by
    rw [← Int.cast_abs]; exact_mod_cast h1
  have hq14R : (q : ℝ) ≤ 14 := by exact_mod_cast hq14
  -- `1/14 ≤ 1/q ≤ |v i − m·q|/q`
  have step1 : (1 : ℝ) / 14 ≤ 1 / (q : ℝ) := one_div_le_one_div_of_le hqR hq14R
  have step2 : (1 : ℝ) / (q : ℝ) ≤ |(((v i - m * q : ℤ)) : ℝ)| / (q : ℝ) := by gcongr
  have goal : (1 : ℝ) / 14 ≤ |(((v i - m * q : ℤ)) : ℝ)| / (q : ℝ) := le_trans step1 step2
  rw [key, abs_div, abs_of_pos hqR]
  simpa using goal

/-- **Corollary: LRC(14) is trivial unless some speed is a multiple of 14.**  If `14 ∤ v i` for every
`i`, the time `t = 1/14` is lonely (margin exactly `1/14` in the worst case).  In particular the
LRC(14) equality extremal `{1,…,13}` is settled here. -/
theorem lonely_of_not_dvd_14 {ι : Type*} (v : ι → ℤ) (hnd : ∀ i, ¬ ((14 : ℤ) ∣ v i)) :
    Lonely 14 v (1 / (14 : ℝ)) := by
  have h := lonely_of_not_dvd v 14 (by norm_num) (by norm_num) hnd
  simpa using h

/-! ## Axiom audit -/
#print axioms lonely_of_not_dvd
#print axioms lonely_of_not_dvd_14

end LRC14
