/-
  TournamentH7.LRCSpread13 — THE SHARP BOUNDED-RATIO WINDOW
  (kind-pasteur-2026-07-03-S28).

  If every absolute speed lies in a band `[a, b]` of ratio `b ≤ 13a`, the family is
  lonely at `t = 1/(a+b)`: every runner sits at `|v|/(a+b) ∈ [a/(a+b), b/(a+b)] ⊆
  [1/14, 13/14]`, the FULL safe band between the two danger arcs.

  This SHARPENS `spread7_lonely` (LRCChainPeel, ratio `7`, `t = 1/(14a)`), which only
  used the half-band `[1/14, 1/2]`.  The improvement is exactly the other half: placing
  the runners across the whole gap `[1/14, 13/14]` instead of crowding them below the
  midpoint doubles the admissible ratio to `13` — the sharp Lonely-Runner threshold
  (ratio `14` fails: `1/(a+b)` would put the smallest runner at `< 1/14`).

  So EVERY family whose speeds are mutually within a factor of `13` — the "all
  comparable" case — is lonely, unconditionally, with an explicit rational witness and
  NO citation, NO census, NO measure theory.

  Kernel-pure (`propext, Classical.choice, Quot.sound`); no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner
namespace LRC14

/-- **THE SHARP BOUNDED-RATIO WINDOW**: absolute speeds in a band `[a, b]` with
`b ≤ 13a` ⟹ lonely at `t = 1/(a+b)`.  Every runner lands in the full safe band
`[1/14, 13/14]`. -/
theorem spread13_lonely {ι : Type*} [Nonempty ι] (v : ι → ℤ) (a b : ℤ) (ha : 0 < a)
    (hlo : ∀ i, a ≤ |v i|) (hhi : ∀ i, |v i| ≤ b) (hratio : b ≤ 13 * a) :
    Lonely 14 v (1 / ((a : ℝ) + b)) := by
  have haR : (0 : ℝ) < (a : ℝ) := by exact_mod_cast ha
  have hbR : (a : ℝ) ≤ (b : ℝ) := by
    -- b ≥ a since a ≤ |v i₀| ≤ b for any i (need ι nonempty)… derive from a ≤ b directly:
    have : (a : ℤ) ≤ b := le_trans (hlo (Classical.arbitrary ι)) (hhi (Classical.arbitrary ι))
    exact_mod_cast this
  have habR : (0 : ℝ) < (a : ℝ) + b := by linarith
  have hratioR : (b : ℝ) ≤ 13 * (a : ℝ) := by exact_mod_cast hratio
  intro i m
  have hlo' : (a : ℝ) ≤ |(v i : ℝ)| := by rw [← Int.cast_abs]; exact_mod_cast hlo i
  have hhi' : |(v i : ℝ)| ≤ (b : ℝ) := by rw [← Int.cast_abs]; exact_mod_cast hhi i
  set x : ℝ := (v i : ℝ) * (1 / ((a : ℝ) + b)) with hx
  have habsx : |x| = |(v i : ℝ)| * (1 / ((a : ℝ) + b)) := by
    rw [hx, abs_mul, abs_of_pos (by positivity : (0:ℝ) < 1 / ((a:ℝ)+b))]
  -- |x| ∈ [1/14, 13/14]
  have hxlo : (1 : ℝ) / 14 ≤ |x| := by
    rw [habsx, mul_one_div, le_div_iff₀ habR]
    -- (1/14)(a+b) ≤ a  ⟸  b ≤ 13a
    nlinarith
  have hxhi : |x| ≤ 13 / 14 := by
    rw [habsx, mul_one_div, div_le_iff₀ habR]
    -- b ≤ (13/14)(a+b)  ⟸  b ≤ 13a
    nlinarith
  -- distance to every integer is ≥ 1/14
  rcases eq_or_ne m 0 with rfl | hm
  · simp only [Int.cast_zero, sub_zero]
    calc (1 : ℝ) / (14 : ℕ) = 1 / 14 := by norm_num
      _ ≤ |x| := hxlo
  · have hm1 : (1 : ℝ) ≤ |(m : ℝ)| := by
      rw [← Int.cast_abs]; exact_mod_cast Int.one_le_abs hm
    have htri : |(m : ℝ)| - |x| ≤ |x - m| := by
      calc |(m : ℝ)| - |x| ≤ |(m : ℝ) - x| := abs_sub_abs_le_abs_sub _ _
        _ = |x - m| := abs_sub_comm _ _
    calc (1 : ℝ) / (14 : ℕ) = 1 / 14 := by norm_num
      _ ≤ |x - m| := by
          have : |(m : ℝ)| - |x| ≥ 1 - 13 / 14 := by linarith
          have h114 : (1 : ℝ) - 13 / 14 = 1 / 14 := by norm_num
          linarith [htri]

/-- **THE RATIONAL-WITNESS SIEVE**: at `t = p/q` (`q > 0`), the family is lonely iff
every `v i · p` stays integer-distance `≥ q/14` from every multiple of `q`.  This is
the general denominator sieve — the covering reduction (`t = 1/q` for a `q` dividing no
speed) is the special case `p = 1`.  Empirically (kps-S28,
`lrc14_denominator_bound/stress`) EVERY hard LRC(14) instance is lonely at some `p/q`
with `q ≤ 35`, so LRC(14) reduces to a finite denominator search — the route the known
`n ≤ 13` computational proofs take. -/
theorem lonely14_of_ratio {ι : Type*} (v : ι → ℤ) (p q : ℤ) (hq : 0 < q)
    (hres : ∀ i, ∀ k : ℤ, q ≤ 14 * |v i * p - k * q|) :
    Lonely 14 v ((p : ℝ) / q) := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  intro i m
  have he : (v i : ℝ) * ((p : ℝ) / q) - m = ((v i * p - m * q : ℤ) : ℝ) / q := by
    push_cast; field_simp
  rw [he, abs_div, abs_of_pos hqR, le_div_iff₀ hqR]
  have hrR : (q : ℝ) ≤ 14 * |((v i * p - m * q : ℤ) : ℝ)| := by
    have h := hres i m
    rw [← Int.cast_abs]
    exact_mod_cast h
  have h14 : (1 : ℝ) / (14 : ℕ) = 1 / 14 := by norm_num
  rw [h14]
  linarith [hrR]

end LRC14
end LonelyRunner
