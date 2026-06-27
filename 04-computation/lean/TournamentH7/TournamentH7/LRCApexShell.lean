/-
  TournamentH7.LRCApexShell -- exact arithmetic core of the LRC(14)
  apex-binding step.

  This file formalizes the part of the THM-568 / HYP-+2909 apex-denominator
  discussion that is genuinely elementary and sorry-free.  At a tight
  `1/14` crossing, an active runner at `+1/14` gives

      14 * (u*a - m*D) = D,

  and an active runner at `-1/14` gives

      14 * (v*a - n*D) = -D.

  The safe conclusion is:

    * `14 ∣ D`;
    * `D ∣ (u+v)` when `a` and `D` are coprime;
    * hence `14 ∣ (u+v)`.

  This does NOT prove `D = 14` for primitive rows.  It leaves the shell-collapse
  theorem `D/14 = 1` as a separate mathematical obligation.
-/

import Mathlib.Tactic

namespace LonelyRunner
namespace ApexShell

/-- A `+1/14` binding equation forces the denominator to be a multiple of `14`. -/
theorem fourteen_dvd_den_of_pos_binding (u a m D : ℤ)
    (h : 14 * (u * a - m * D) = D) : (14 : ℤ) ∣ D :=
  ⟨u * a - m * D, by simpa [mul_comm] using h.symm⟩

/-- A `-1/14` binding equation also forces the denominator to be a multiple of `14`. -/
theorem fourteen_dvd_den_of_neg_binding (u a m D : ℤ)
    (h : 14 * (u * a - m * D) = -D) : (14 : ℤ) ∣ D := by
  refine ⟨-(u * a - m * D), ?_⟩
  linarith

/-- Opposite-side binders at `±1/14` have pair sum whose product with the
numerator is divisible by the denominator. -/
theorem den_dvd_pair_sum_mul_num_of_opposite_bindings (u v a m n D : ℤ)
    (hpos : 14 * (u * a - m * D) = D)
    (hneg : 14 * (v * a - n * D) = -D) :
    D ∣ (u + v) * a := by
  refine ⟨m + n, ?_⟩
  have hsum : u * a - m * D + (v * a - n * D) = 0 := by
    have h : 14 * (u * a - m * D + (v * a - n * D)) = 0 := by
      linarith
    nlinarith
  nlinarith

/-- If the rational time `a/D` is in lowest terms, opposite-side binders have
pair sum divisible by `D`. -/
theorem den_dvd_pair_sum_of_opposite_bindings (u v a m n D : ℤ)
    (hcop : IsCoprime D a)
    (hpos : 14 * (u * a - m * D) = D)
    (hneg : 14 * (v * a - n * D) = -D) :
    D ∣ u + v := by
  exact hcop.dvd_of_dvd_mul_right
    (den_dvd_pair_sum_mul_num_of_opposite_bindings u v a m n D hpos hneg)

/-- The HYP-+2909 forward readout: opposite-side tight binders are antipodal
modulo `14`. -/
theorem fourteen_dvd_pair_sum_of_opposite_bindings (u v a m n D : ℤ)
    (hcop : IsCoprime D a)
    (hpos : 14 * (u * a - m * D) = D)
    (hneg : 14 * (v * a - n * D) = -D) :
    (14 : ℤ) ∣ u + v := by
  exact dvd_trans (fourteen_dvd_den_of_pos_binding u a m D hpos)
    (den_dvd_pair_sum_of_opposite_bindings u v a m n D hcop hpos hneg)

/-- Same pair-sum divisibility with the two signs swapped. -/
theorem fourteen_dvd_pair_sum_of_opposite_bindings_swap (u v a m n D : ℤ)
    (hcop : IsCoprime D a)
    (hneg : 14 * (u * a - m * D) = -D)
    (hpos : 14 * (v * a - n * D) = D) :
    (14 : ℤ) ∣ u + v := by
  have h : (14 : ℤ) ∣ v + u :=
    fourteen_dvd_pair_sum_of_opposite_bindings v u a n m D hcop hpos hneg
  simpa [add_comm] using h

/-! ## Axiom audit -/

#print axioms fourteen_dvd_den_of_pos_binding
#print axioms fourteen_dvd_den_of_neg_binding
#print axioms den_dvd_pair_sum_of_opposite_bindings
#print axioms fourteen_dvd_pair_sum_of_opposite_bindings

end ApexShell
end LonelyRunner
