/-
  TournamentH7.LRCDilatedSieve — THM-1013, the dilated-AP sieve (boxeph-2026-07-18-S86).

  The dilated generalization of the empty-circle sieve (`sieve_frac`).  If every speed
  lies at distance `≥ d` from the lattice `n·d·ℤ`, then `t = 1/(n·d)` is `1/n`-lonely:
      ‖v·t‖ = dist(v, n·d·ℤ)/(n·d) ≥ d/(n·d) = 1/n.
  For `n = 13` this is the WITNESS the compact LRC(14) minimizers were missing: the
  families `d·{1..12} ∪ {killer}` (which defeat fill-1, far-element descent, the sharp
  recursion, the ordinary sieve, and the measure bound) all satisfy the band condition
  at their dilation `d`, so `M ≥ 1/13 > 1/14` by this one explicit rational witness.

  Elementary: one division inequality on ℝ/ℤ.  Companion to `sieve_frac`,
  `fill1_perturbation`, `descent_general`.
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner

variable {ι : Type*}

/-- **THM-1013 (dilated-AP sieve).**  If every speed is at distance `≥ d` from the
lattice `n·d·ℤ` (i.e. `|v i − n·d·m| ≥ d` for all integers `m`), then `t = 1/(n·d)` is
`1/n`-lonely. -/
theorem dilated_sieve (n d : ℕ) (hn : 1 ≤ n) (hd : 1 ≤ d) (v : ι → ℤ)
    (hband : ∀ i, ∀ m : ℤ, (d : ℤ) ≤ |v i - (n : ℤ) * d * m|) :
    Lonely n v (1 / ((n : ℝ) * d)) := by
  have hn0 : (0 : ℝ) < n := by exact_mod_cast (by omega : 0 < n)
  have hd0 : (0 : ℝ) < d := by exact_mod_cast (by omega : 0 < d)
  have hnd : (0 : ℝ) < (n : ℝ) * d := by positivity
  intro i m
  -- put v_i·t − m over the common denominator n·d
  have key : (v i : ℝ) * (1 / ((n : ℝ) * d)) - m
      = ((v i - (n : ℤ) * d * m : ℤ) : ℝ) / ((n : ℝ) * d) := by
    push_cast; field_simp
  rw [key, abs_div, abs_of_pos hnd]
  -- numerator ≥ d, denominator n·d, so ratio ≥ d/(n·d) = 1/n
  have hb : (d : ℝ) ≤ |((v i - (n : ℤ) * d * m : ℤ) : ℝ)| := by
    rw [← Int.cast_abs]; exact_mod_cast hband i m
  rw [le_div_iff₀ hnd]
  -- goal: 1/n * (n*d) ≤ |num| ;  1/n*(n*d) = d ≤ |num|
  have hsimp : (1 : ℝ) / n * ((n : ℝ) * d) = d := by field_simp
  rw [hsimp]; exact hb

/-- The LRC(14) reading: the same witness clears the `1/14` band with margin
(`1/13 > 1/14`). -/
theorem dilated_sieve_lonely14 (d : ℕ) (hd : 1 ≤ d) (v : ι → ℤ)
    (hband : ∀ i, ∀ m : ℤ, (d : ℤ) ≤ |v i - 13 * d * m|) :
    Lonely 14 v (1 / (13 * d : ℝ)) := by
  have h13 := dilated_sieve 13 d (by norm_num) hd v (by
    intro i m; simpa using hband i m)
  intro i m
  have := h13 i m
  have h14le13 : (1 : ℝ) / 14 ≤ 1 / 13 := by norm_num
  have hcast : ((13 : ℕ) : ℝ) * d = (13 * d : ℝ) := by push_cast; ring
  have := h13 i m
  -- h13 gives 1/13 ≤ |...|; 1/14 ≤ 1/13
  calc (1 : ℝ) / 14 ≤ 1 / 13 := by norm_num
    _ ≤ |(v i : ℝ) * (1 / ((13 : ℕ) * (d : ℝ))) - m| := this
    _ = |(v i : ℝ) * (1 / (13 * (d : ℝ))) - m| := by norm_num

end LonelyRunner

#print axioms LonelyRunner.dilated_sieve
#print axioms LonelyRunner.dilated_sieve_lonely14
