/-
  TournamentH7.LRCAlignedResonance

  The arithmetic residue of the zero-color branch.  If a stalk has common
  speed factor `h`, then its potentially aligned multipliers are precisely
  the nonzero residues `p < q` for which `q ∣ h * p`.  Cancelling the common
  gcd identifies one reduced congruence class and gives the exact count
  `gcd h q - 1`.

  Tournament/Kakeya carrier audit: zero determinant collapses the runner
  slopes to one projective direction, but it does not erase the multiplier
  (the needle intercept/phase).  Thus the faithful vertices are still
  multiplier events in a zero-color fiber, and the tie path is their natural
  order in `p`.  The divisibility predicate is the binary switch; its score
  histogram consists of one reduced residue class.  A static tournament on
  slopes destroys the intercept and therefore cannot recover the exact
  activity count.

  What is proved here is only the arithmetic resonance filter and the
  integrality atom turning sufficiently strong closeness into exact
  proportionality.  It does not prove that zero overlap color supplies the
  required common witness parameter, nor that local proportionalities glue
  across a connected seven-stalk with bounded aggregation multiplicity.

  No `sorry`; no `native_decide`.
-/

import Mathlib.Algebra.Order.Ring.Int
import Mathlib.Data.Nat.Factorization.Basic

namespace LonelyRunner
namespace LRC14Grand

open Finset

/-- The nontrivial multipliers whose common speed factor can resonate with
the modulus. -/
def alignedResonantMultipliers (h q : ℕ) : Finset ℕ :=
  (Finset.Ioo 0 q).filter fun p => q ∣ h * p

/-- Cancel the common speed/modulus gcd.  The zero-color resonance condition
is one reduced divisibility condition on the multiplier. -/
theorem dvd_mul_iff_reduced_modulus_dvd (h q p : ℕ) (hq : 0 < q) :
    q ∣ h * p ↔ q / Nat.gcd h q ∣ p := by
  let g := Nat.gcd h q
  have hg : 0 < g := by
    simpa [g] using Nat.gcd_pos_of_pos_right h hq
  have hg_left : g ∣ h := by
    simpa [g] using Nat.gcd_dvd_left h q
  have hg_right : g ∣ q := by
    simpa [g] using Nat.gcd_dvd_right h q
  have hrhs : (h / g * p) * g = h * p := by
    calc
      (h / g * p) * g = (h / g * g) * p := by ac_rfl
      _ = h * p := by rw [Nat.div_mul_cancel hg_left]
  have hcoprime : Nat.Coprime (q / g) (h / g) :=
    (Nat.coprime_div_gcd_div_gcd hg).symm
  change q ∣ h * p ↔ q / g ∣ p
  calc
    q ∣ h * p ↔ (q / g) * g ∣ (h / g * p) * g := by
      rw [Nat.div_mul_cancel hg_right, hrhs]
    _ ↔ q / g ∣ (h / g) * p := mul_dvd_mul_iff_right hg.ne'
    _ ↔ q / g ∣ p := hcoprime.dvd_mul_left

/-- There are exactly `gcd h q - 1` nonzero resonant multipliers below `q`.
The missing one from the closed interval `(0,q]` is its right endpoint `q`.
-/
theorem alignedResonantMultipliers_card (h q : ℕ) (hq : 0 < q) :
    (alignedResonantMultipliers h q).card = Nat.gcd h q - 1 := by
  let g := Nat.gcd h q
  let d := q / g
  have hg : 0 < g := by
    simpa [g] using Nat.gcd_pos_of_pos_right h hq
  have hgq : g ∣ q := by
    simpa [g] using Nat.gcd_dvd_right h q
  have hdpos : 0 < d := by
    apply Nat.div_pos
    · exact Nat.le_of_dvd hq hgq
    · exact hg
  have hdvdq : d ∣ q := by
    refine ⟨g, ?_⟩
    simpa [d] using (Nat.div_mul_cancel hgq).symm
  have hquot : q / d = g := by
    have hfactor : q = d * g := by
      simpa [d] using (Nat.div_mul_cancel hgq).symm
    rw [hfactor, mul_comm, Nat.mul_div_left _ hdpos]
  have hset :
      alignedResonantMultipliers h q =
        (Finset.Ioo 0 q).filter (fun p => d ∣ p) := by
    ext p
    simp only [alignedResonantMultipliers, Finset.mem_filter,
      Finset.mem_Ioo]
    rw [dvd_mul_iff_reduced_modulus_dvd h q p hq]
  rw [hset]
  change ((Finset.Ioo 0 q).filter (fun p => d ∣ p)).card = g - 1
  have hqnot : q ∉ (Finset.Ioo 0 q).filter (fun p => d ∣ p) := by
    simp
  have hinsert :
      insert q ((Finset.Ioo 0 q).filter (fun p => d ∣ p)) =
        (Finset.Ioc 0 q).filter (fun p => d ∣ p) := by
    rw [← Finset.Ioo_insert_right hq]
    rw [Finset.filter_insert, if_pos hdvdq]
  have hcard :
      ((Finset.Ioo 0 q).filter (fun p => d ∣ p)).card + 1 = g := by
    calc
      ((Finset.Ioo 0 q).filter (fun p => d ∣ p)).card + 1 =
          (insert q ((Finset.Ioo 0 q).filter (fun p => d ∣ p))).card :=
        (Finset.card_insert_of_notMem hqnot).symm
      _ = ((Finset.Ioc 0 q).filter (fun p => d ∣ p)).card := by
        rw [hinsert]
      _ = q / d := Nat.Ioc_filter_dvd_card_eq_div q d
      _ = g := hquot
  omega

/-- **Integrality atom for an aligned needle.**  Once the scaled error is
strictly smaller than the common-factor window and that window is at most the
scale, the integer cross-product error has absolute value below one and hence
vanishes. -/
theorem eq_of_aligned_resonance_closeness
    (M h q p r : ℤ)
    (hM : 0 < M) (_hh : 0 ≤ h)
    (hwindow : h * q ≤ 14 * M)
    (hclose : 14 * M * |h * p - r * q| < h * q) :
    h * p = r * q := by
  have hscale : (0 : ℤ) < 14 * M := mul_pos (by omega) hM
  have hscaled :
      14 * M * |h * p - r * q| < 14 * M * 1 :=
    lt_of_lt_of_le hclose (by simpa using hwindow)
  have habs : |h * p - r * q| < 1 :=
    (mul_lt_mul_iff_right₀ hscale).mp hscaled
  exact sub_eq_zero.mp (Int.abs_lt_one_iff.mp habs)

end LRC14Grand
end LonelyRunner

#print axioms LonelyRunner.LRC14Grand.dvd_mul_iff_reduced_modulus_dvd
#print axioms LonelyRunner.LRC14Grand.alignedResonantMultipliers_card
#print axioms LonelyRunner.LRC14Grand.eq_of_aligned_resonance_closeness
