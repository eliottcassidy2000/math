/-
  TournamentH7.LRCBandMargin — THE RECIPROCAL-WINDOW MARGIN + the LRC(13) analog of
  spread13, and the FIRST NECESSARY CONDITION of the tight-locus rigidity
  (kind-pasteur-2026-07-05-S1, HYP-4096).

  `band_margin_reciprocal`: if every speed lies in a band `[a,b]` (`0<a`), then at
  `t = 1/(a+b)` EVERY runner sits at circle-distance `≥ a/(a+b)` from `ℤ`.  This is the
  explicit-margin refinement of klein's `spread13_lonely` (which only extracts the
  `1/14` threshold at ratio `13`).  The margin `a/(a+b)` is exact (attained by both the
  min and the max speed) and feeds `lonely_of_window_margin` with `β = a/(a+b)` for ANY
  band ratio — not just `≤ 13`.

  `spread12_lonely13`: the LRC(13) twin of `spread13_lonely` — band ratio `b ≤ 12a`
  ⟹ `Lonely 13` at `t=1/(a+b)` (margin `≥ 1/13`).

  RIGIDITY CONNECTION (HYP-4096).  The tight-locus rigidity — `W` a primitive 12-set,
  `M(W)=max_t min_w ||wt|| = 1/13`  ⟺  `W = {1,…,12}` — is the uncited lemma klein-S132
  `tight_ap_free_rider` needs to cover the WHOLE tight case of hcomp.  It is VERIFIED
  DECISIVELY (kps-S1: exhaustive over all 1820 primitive 12-subsets of [1,16] AP-unique;
  531376 residue-system adversaries AP-unique; 400k random to height 200; second value
  2/25) but NOT proved.  `band_margin_reciprocal` supplies its FIRST NECESSARY CONDITION:
  since `M(W) ≥ a/(a+b) = w_min/(w_min+w_max)`, tightness forces `w_min/(w_min+w_max) ≤
  1/13`, i.e. **`w_max ≥ 12·w_min`** (an EQUALITY for the AP: `12 = 12·1`).  Equivalently,
  every NARROW band (`w_max < 12·w_min`) is LOOSE — a provable slice of the rigidity's
  contrapositive, with the explicit loose witness `t=1/(w_min+w_max)`.
  (A second necessary condition, `w_max ≤ 12·w_2nd` via LRC(12) on the peeled 11-subset,
  is proved in the reflection; formalizing it needs the interval-covering machinery.)

  Kernel-pure (`propext, Classical.choice, Quot.sound`); no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner
namespace BandMargin

/-- **THE RECIPROCAL-WINDOW MARGIN.**  Speeds in a band `[a,b]` with `0<a` ⟹ at
`t = 1/(a+b)` every runner is at circle-distance `≥ a/(a+b)` from every integer.
The bound is exact: the min speed lands at `a/(a+b)` and the max at `b/(a+b) =
1 - a/(a+b)`, both at distance exactly `a/(a+b)`. -/
theorem band_margin_reciprocal {ι : Type*} [Nonempty ι] (v : ι → ℤ) (a b : ℤ) (ha : 0 < a)
    (hlo : ∀ i, a ≤ |v i|) (hhi : ∀ i, |v i| ≤ b) :
    ∀ i, ∀ m : ℤ, (a : ℝ) / ((a : ℝ) + b) ≤ |(v i : ℝ) * (1 / ((a : ℝ) + b)) - m| := by
  have haR : (0 : ℝ) < (a : ℝ) := by exact_mod_cast ha
  have hab : (a : ℤ) ≤ b := le_trans (hlo (Classical.arbitrary ι)) (hhi (Classical.arbitrary ι))
  have hbR : (a : ℝ) ≤ (b : ℝ) := by exact_mod_cast hab
  have habR : (0 : ℝ) < (a : ℝ) + b := by linarith
  have hne : ((a : ℝ) + b) ≠ 0 := ne_of_gt habR
  intro i m
  have hlo' : (a : ℝ) ≤ |(v i : ℝ)| := by rw [← Int.cast_abs]; exact_mod_cast hlo i
  have hhi' : |(v i : ℝ)| ≤ (b : ℝ) := by rw [← Int.cast_abs]; exact_mod_cast hhi i
  set x : ℝ := (v i : ℝ) * (1 / ((a : ℝ) + b)) with hx
  have habsx : |x| = |(v i : ℝ)| * (1 / ((a : ℝ) + b)) := by
    rw [hx, abs_mul, abs_of_pos (one_div_pos.mpr habR)]
  -- |x| · (a+b) = |v i|  (clear the denominator once, robustly)
  have hxs : |x| * ((a : ℝ) + b) = |(v i : ℝ)| := by
    rw [habsx, mul_one_div, div_mul_cancel₀ _ hne]
  -- a/(a+b) ≤ |x| ≤ b/(a+b)
  have hxlo : (a : ℝ) / ((a : ℝ) + b) ≤ |x| := by
    rw [div_le_iff₀ habR, hxs]; exact hlo'
  have hxhi : |x| ≤ (b : ℝ) / ((a : ℝ) + b) := by
    rw [le_div_iff₀ habR, hxs]; exact hhi'
  -- a/(a+b) and b/(a+b) are complementary in [0,1]
  have hcompl : (a : ℝ) / ((a : ℝ) + b) + (b : ℝ) / ((a : ℝ) + b) = 1 := by
    rw [← add_div, div_self hne]
  rcases eq_or_ne m 0 with rfl | hm
  · simpa using hxlo
  · have hm1 : (1 : ℝ) ≤ |(m : ℝ)| := by
      rw [← Int.cast_abs]; exact_mod_cast Int.one_le_abs hm
    have htri : |(m : ℝ)| - |x| ≤ |x - m| := by
      calc |(m : ℝ)| - |x| ≤ |(m : ℝ) - x| := abs_sub_abs_le_abs_sub _ _
        _ = |x - m| := abs_sub_comm _ _
    -- |x-m| ≥ |m| - |x| ≥ 1 - b/(a+b) = a/(a+b)
    linarith [htri, hm1, hxhi, hcompl]

/-- **THE LRC(13) BOUNDED-RATIO WINDOW** — the `spread13_lonely` twin at threshold `12`.
Speeds in `[a,b]` with `b ≤ 12a` ⟹ `Lonely 13` at `t = 1/(a+b)` (margin `≥ 1/13`).
CONTRAPOSITIVE = the tight-locus rigidity's first necessary condition: a primitive
12-set with `M(W)=1/13` must have `w_max ≥ 12·w_min`. -/
theorem spread12_lonely13 {ι : Type*} [Nonempty ι] (v : ι → ℤ) (a b : ℤ) (ha : 0 < a)
    (hlo : ∀ i, a ≤ |v i|) (hhi : ∀ i, |v i| ≤ b) (hratio : b ≤ 12 * a) :
    Lonely 13 v (1 / ((a : ℝ) + b)) := by
  have haR : (0 : ℝ) < (a : ℝ) := by exact_mod_cast ha
  have hbR : (b : ℝ) ≤ 12 * (a : ℝ) := by exact_mod_cast hratio
  have habR : (0 : ℝ) < (a : ℝ) + b := by
    have h1 : (a : ℤ) ≤ b := le_trans (hlo (Classical.arbitrary ι)) (hhi (Classical.arbitrary ι))
    have h2 : (a : ℝ) ≤ (b : ℝ) := by exact_mod_cast h1
    linarith
  -- 1/13 ≤ a/(a+b) since (a+b) ≤ 13a
  have hkey : (1 : ℝ) / 13 ≤ (a : ℝ) / ((a : ℝ) + b) := by
    rw [le_div_iff₀ habR]; linarith
  intro i m
  calc (1 : ℝ) / (13 : ℕ) = 1 / 13 := by norm_num
    _ ≤ (a : ℝ) / ((a : ℝ) + b) := hkey
    _ ≤ |(v i : ℝ) * (1 / ((a : ℝ) + b)) - m| :=
        band_margin_reciprocal v a b ha hlo hhi i m

#print axioms band_margin_reciprocal
#print axioms spread12_lonely13

end BandMargin
end LonelyRunner
