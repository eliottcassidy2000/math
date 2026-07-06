/-
  TournamentH7.LRCMidrangeWitness — THE MIDRANGE WITNESS (general lower bound for the covering-min)
  (klein-2026-07-05-S141, HYP-4161).

  For ANY family of speeds all lying in `[m, Mx]` (`0 < m ≤ Mx`), the single rational time
  `t = 1/(m+Mx)` clears every speed by `m/(m+Mx)`:

        ∀ i, ∀ k : ℤ,   m/(m+Mx) ≤ ‖v_i · t − k‖.

  So `M(S) := max_t min_i ‖v_i t‖ ≥ v_min/(v_min+v_max)` for every family — the elementary
  midrange lower bound (THM-526; opus-S85's "phases cluster symmetrically around ½, the extreme
  speeds bind").  Full consecutive combs SATURATE it (`M({a,…,a+r−1}) = a/(2a+r−1)`, opus-S85), so
  it is exact on the compressed extremal.  It is the `≥`-half every gap/rigidity argument needs, and
  e.g. gives `v_max ≤ 13 v_min ⟹ M ≥ 1/14` immediately (`midrange_margin_compressed`).

  Reduction: `m/(m+Mx) ≤ |v_i/(m+Mx) − k|  ⟺  m ≤ |v_i − k·(m+Mx)|`, the latter a pure integer
  inequality (cases on `k ≤ 0` / `k ≥ 1`).

  Kernel-pure (propext, Classical.choice, Quot.sound); no sorry, no native_decide.
-/
import Mathlib

namespace LonelyRunner

/-- **The integer core.**  For `m ≤ v ≤ Mx` with `0 < m`, every integer `k` has
`m ≤ |v − k·(m+Mx)|`. -/
lemma midrange_int_core (v m Mx : ℤ) (hm : 0 < m) (hlo : m ≤ v) (hhi : v ≤ Mx) (k : ℤ) :
    m ≤ |v - k * (m + Mx)| := by
  have hw : 0 ≤ m + Mx := by omega
  by_cases hk : k ≤ 0
  · -- k ≤ 0 : v − k·(m+Mx) ≥ v ≥ m
    have hkw : k * (m + Mx) ≤ 0 := mul_nonpos_of_nonpos_of_nonneg hk hw
    have hge : m ≤ v - k * (m + Mx) := by omega
    exact le_trans hge (le_abs_self _)
  · -- k ≥ 1 : v − k·(m+Mx) ≤ Mx − (m+Mx) = −m
    have hk1 : 1 ≤ k := by omega
    have hkw : m + Mx ≤ k * (m + Mx) := le_mul_of_one_le_left hw hk1
    have hle : m ≤ -(v - k * (m + Mx)) := by omega
    exact le_trans hle (neg_le_abs _)

/-- **The midrange witness.**  Every speed in `[m, Mx]` is integer-distance `≥ m/(m+Mx)` from every
integer at the time `t = 1/(m+Mx)`.  General `ι`; the extreme speeds `m`, `Mx` bind. -/
theorem midrange_margin {ι : Type*} (v : ι → ℤ) (m Mx : ℤ) (hm : 0 < m) (hmMx : m ≤ Mx)
    (hb : ∀ i, m ≤ v i ∧ v i ≤ Mx) (i : ι) (k : ℤ) :
    (m : ℝ) / ((m : ℝ) + Mx) ≤ |(v i : ℝ) * (1 / ((m : ℝ) + Mx)) - k| := by
  have hmR : (0 : ℝ) < (m : ℝ) := by exact_mod_cast hm
  have hMxR : (m : ℝ) ≤ (Mx : ℝ) := by exact_mod_cast hmMx
  have hs : (0 : ℝ) < (m : ℝ) + Mx := by linarith
  have hs0 : ((m : ℝ) + Mx) ≠ 0 := ne_of_gt hs
  obtain ⟨hlo, hhi⟩ := hb i
  have hrw : (v i : ℝ) * (1 / ((m : ℝ) + Mx)) - k
      = (((v i - k * (m + Mx) : ℤ)) : ℝ) / ((m : ℝ) + Mx) := by
    push_cast
    field_simp
  rw [hrw, abs_div, abs_of_pos hs]
  gcongr
  rw [← Int.cast_abs]
  exact_mod_cast midrange_int_core (v i) m Mx hm hlo hhi k

/-- **Compressed corollary**: if `v_max ≤ 13 v_min` then the midrange time clears every runner by
`≥ 1/14`.  (The `≥`-half of the compressed leaf; equality iff the AP, opus-S85.) -/
theorem midrange_margin_compressed {ι : Type*} (v : ι → ℤ) (m Mx : ℤ) (hm : 0 < m)
    (hmMx : m ≤ Mx) (hcomp : Mx ≤ 13 * m) (hb : ∀ i, m ≤ v i ∧ v i ≤ Mx) (i : ι) (k : ℤ) :
    (1 : ℝ) / 14 ≤ |(v i : ℝ) * (1 / ((m : ℝ) + Mx)) - k| := by
  refine le_trans ?_ (midrange_margin v m Mx hm hmMx hb i k)
  have hmR : (0 : ℝ) < (m : ℝ) := by exact_mod_cast hm
  have hs : (0 : ℝ) < (m : ℝ) + Mx := by
    have : (m : ℝ) ≤ Mx := by exact_mod_cast hmMx
    linarith
  have hcompR : (Mx : ℝ) ≤ 13 * m := by exact_mod_cast hcomp
  rw [div_le_div_iff₀ (by norm_num) hs]
  linarith

#print axioms midrange_margin
#print axioms midrange_margin_compressed

end LonelyRunner
