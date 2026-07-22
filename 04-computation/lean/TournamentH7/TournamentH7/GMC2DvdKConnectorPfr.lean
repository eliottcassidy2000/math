import TournamentH7.GMC2DvdKConnector
import TournamentH7.GMC2DvdKUnitOrigin
import TournamentH7.GMC2DvdKMultiplicativeClosing
import TournamentH7.GMC2DvdKFrameHSide
import TournamentH7.GMC2DvdKFrameDegree
import TournamentH7.GMC2DvdKHderivAssembly

/-!
# The concrete φ-connections and the closed `hderiv` (GMC(2) formalization capstone)

This module glues death-star's transpose `φ = map(ofPowerSeries) ∘ tauHom` to the *concrete* Weierstrass
factors of `Φ = Xᴹ − t·R` (mac-mini's `smallRootFactor`, `weierstrassUnit`), producing the six per-object
inputs that kps's `GMC2DvdKHderivAssembly.hderiv_of_transpose_glue` consumes:

* `Pfr := φ(↑smallRootFactor)` — the transported distinguished factor;
* `coeff_zero_Pfr` — its `t⁰` term is `xᴹ` (`single M 1`), so `isUnit_Pfr` (unit, `xᴹ ≠ 0` in `F⸨x⸩`);
* `hlt_Pfr` — every positive-`t` term of `Pfr` has `x`-degree `< M` (monic, degree `M`), feeding the (c)
  degree lemma `GMC2DvdKFrameDegree.xCoeff0_logDeriv_eq_zero_of_monic`;
* `phiFrame_factor` — `PhiFrame = Pfr · φ(weierstrassUnit)` (from the connector `phi_Phi` + `Φ = P·h`);
* `xCoeff0_phi_weierstrassUnit` — the **bridge** `xCoeff0(φ h) = unitCoeff0` (`h(0,t)`), so `hg` is a unit
  (`h(0,0)=1`) and the assembly's conclusion is the literal DvdK derivative.

Composing these with the assembly and mac-mini's monic (c)-lemma gives **`hderiv_concrete`**: the
`t`-derivative of `h(0,t)` vanishes, with the *only* remaining hypothesis being `hvanish`, the vanishing of
every positive frame moment `[x^{M·m}] R(x)^m = 0` (the frame image of `D_m = 0`).  Feeding that into
mac-mini's `smallRootFactor_coeff0_eq_of_derivative_vanishes'` closes the multiplicative DvdK crux
(`P.coeff 0 = −t·r₀`) in `CharZero`, exp/log/Puiseux-free.
-/

open PowerSeries GMC2DvdKTranspose GMC2DvdKFrame GMC2DvdKWeierstrass

namespace GMC2DvdKConnectorPfr

variable {F : Type*} [Field F]

/-- `coeff k (φ H) = ofPowerSeries (mk (fun n => coeff k (coeff n H)))`. -/
theorem coeff_phi (H : PowerSeries (PowerSeries F)) (k : ℕ) :
    PowerSeries.coeff (R := LaurentSeries F) k (phi H)
      = HahnSeries.ofPowerSeries ℤ F
          (PowerSeries.mk fun n => PowerSeries.coeff (R := F) k
            (PowerSeries.coeff (R := PowerSeries F) n H)) := by
  rw [phi, RingHom.comp_apply, coeff_map]
  congr 1
  refine PowerSeries.ext fun n => ?_
  rw [PowerSeries.coeff_mk]
  exact GMC2DvdKTranspose.coeff_coeff_tau H k n

/-- The transported distinguished factor `Pfr := φ(smallRootFactor)`. -/
noncomputable def Pfr (R : Polynomial F) (M : ℕ) : PowerSeries (LaurentSeries F) :=
  phi ((smallRootFactor R M : PowerSeries (PowerSeries F)))

/-- `coeff 0 (Pfr) = single M 1` (the `t⁰` term is `xᴹ`), from `smallRootFactor ≡ Xᴹ mod t`. -/
theorem coeff_zero_Pfr (R : Polynomial F) (M : ℕ) :
    PowerSeries.coeff (R := LaurentSeries F) 0 (Pfr R M) = HahnSeries.single (M : ℤ) 1 := by
  rw [Pfr, coeff_phi]
  have hmodt : PowerSeries.map (PowerSeries.constantCoeff (R := F))
      ((smallRootFactor R M : (PowerSeries F)⟦X⟧)) = PowerSeries.X ^ M :=
    GMC2DvdKUnitOrigin.map_constantCoeff_smallRootFactor R M
  have hmk : (PowerSeries.mk fun n => PowerSeries.coeff (R := F) 0
        (PowerSeries.coeff (R := PowerSeries F) n (smallRootFactor R M : (PowerSeries F)⟦X⟧)))
      = (PowerSeries.X ^ M : PowerSeries F) := by
    refine PowerSeries.ext fun n => ?_
    rw [PowerSeries.coeff_mk, PowerSeries.coeff_zero_eq_constantCoeff, ← coeff_map, hmodt]
  rw [hmk, HahnSeries.ofPowerSeries_X_pow]

/-- `Pfr` is a unit (its `t⁰` term `xᴹ ≠ 0` in the field `F⸨x⸩`). -/
theorem isUnit_Pfr (R : Polynomial F) (M : ℕ) : IsUnit (Pfr R M) := by
  rw [isUnit_iff_constantCoeff_ne_zero, ← PowerSeries.coeff_zero_eq_constantCoeff, coeff_zero_Pfr]
  simp

/-- **The frame factorization** `PhiFrame = Pfr · φ(h)`, from `φ(Φ)=PhiFrame` + `Φ = P·h`. -/
theorem phiFrame_factor (R : Polynomial F) (M : ℕ) :
    PhiFrame (Polynomial.aeval (HahnSeries.single (1 : ℤ) (1 : F)) R) M
      = Pfr R M * phi (PowerSeries.weierstrassUnit (Phi R M) (phi_residue_ne_zero R M)) := by
  have heq := PowerSeries.isWeierstrassFactorization_weierstrassDistinguished_weierstrassUnit
    (phi_residue_ne_zero R M)
  rw [← GMC2DvdKConnector.phi_Phi, Pfr, ← map_mul]
  congr 1
  exact heq.eq_mul

/-- **The bridge** `xCoeff0(φ h) = unitCoeff0` (the frame `x⁰` of the transported unit is `h(0,t)`).
Via `xCoeff0 ∘ (map ofPowerSeries) = map constantCoeff` (kps) + the `τ` coefficient swap. -/
theorem xCoeff0_phi_weierstrassUnit (R : Polynomial F) (M : ℕ) :
    xCoeff0 (phi (PowerSeries.weierstrassUnit (Phi R M) (phi_residue_ne_zero R M)))
      = GMC2DvdKMultiplicativeClosing.unitCoeff0 R M := by
  rw [phi, RingHom.comp_apply, GMC2DvdKFrameHSide.xCoeff0_map_ofPowerSeries]
  refine PowerSeries.ext fun k => ?_
  rw [coeff_map, GMC2DvdKMultiplicativeClosing.unitCoeff0,
    ← PowerSeries.coeff_zero_eq_constantCoeff, ← PowerSeries.coeff_zero_eq_constantCoeff]
  exact GMC2DvdKTranspose.coeff_coeff_tau
    (PowerSeries.weierstrassUnit (Phi R M) (phi_residue_ne_zero R M)) k 0

/-- For `t^n` with `n ≥ 1`, `Pfr` has `x`-degree `< M` (the `x^M` term is only at `t^0`, monic). -/
theorem hlt_Pfr (R : Polynomial F) (M : ℕ) :
    ∀ n, 1 ≤ n → (PowerSeries.coeff (R := LaurentSeries F) n (Pfr R M)).support ⊆ Set.Iio (M : ℤ) := by
  intro n hn j hj
  rw [Set.mem_Iio]
  by_contra hle
  push_neg at hle
  rw [HahnSeries.mem_support, Pfr, coeff_phi] at hj
  apply hj
  obtain ⟨m, rfl⟩ : ∃ m : ℕ, (m : ℤ) = j := ⟨j.toNat, Int.toNat_of_nonneg (by omega)⟩
  rw [HahnSeries.ofPowerSeries_apply_coeff, PowerSeries.coeff_mk, Polynomial.coeff_coe]
  have hmM : M ≤ m := by exact_mod_cast hle
  rcases eq_or_lt_of_le hmM with hm | hm
  · have hlc : (smallRootFactor R M).coeff M = 1 := by
      have hmon := smallRootFactor_monic R M
      rwa [Polynomial.Monic, Polynomial.leadingCoeff, smallRootFactor_natDegree] at hmon
    rw [← hm, hlc, PowerSeries.coeff_one, if_neg (by omega)]
  · rw [Polynomial.coeff_eq_zero_of_natDegree_lt (by rw [smallRootFactor_natDegree]; omega), map_zero]

/-- **The closed `hderiv`.**  Composing the six φ-connections with kps's assembly and mac-mini's monic
(c)-lemma: given only the vanishing of every positive frame moment `[x^{M·m}] R(x)^m = 0` (`hvanish`, the
frame image of `D_m = 0`), the `t`-derivative of the Weierstrass unit's `x⁰`-coefficient `h(0,t)` vanishes.
This is exactly the `hderiv` that mac-mini's multiplicative closing consumes. -/
theorem hderiv_concrete (R : Polynomial F) (M : ℕ)
    (hvanish : ∀ m : ℕ, 1 ≤ m →
      ((Polynomial.aeval (HahnSeries.single (1 : ℤ) (1 : F)) R) ^ m).coeff ((M : ℤ) * m) = 0) :
    PowerSeries.derivativeFun (GMC2DvdKMultiplicativeClosing.unitCoeff0 R M) = 0 := by
  have hg : IsUnit (xCoeff0 (phi (PowerSeries.weierstrassUnit (Phi R M) (phi_residue_ne_zero R M)))) := by
    rw [xCoeff0_phi_weierstrassUnit R M, PowerSeries.isUnit_iff_constantCoeff,
      GMC2DvdKUnitOrigin.unitCoeff0_constantCoeff_eq_one R M]
    exact isUnit_one
  exact GMC2DvdKHderivAssembly.hderiv_of_transpose_glue
    (Polynomial.aeval (HahnSeries.single (1 : ℤ) (1 : F)) R) M
    (PowerSeries.weierstrassUnit (Phi R M) (phi_residue_ne_zero R M))
    (PowerSeries.isWeierstrassFactorization_weierstrassDistinguished_weierstrassUnit
      (phi_residue_ne_zero R M)).isUnit
    (Pfr R M) (phiFrame_factor R M) (isUnit_Pfr R M)
    (GMC2DvdKFrameDegree.xCoeff0_logDeriv_eq_zero_of_monic (Pfr R M) (isUnit_Pfr R M) M
      (coeff_zero_Pfr R M) (hlt_Pfr R M))
    hg hvanish (GMC2DvdKMultiplicativeClosing.unitCoeff0 R M) (xCoeff0_phi_weierstrassUnit R M)

/-- **The multiplicative DvdK crux, closed modulo `hvanish` (`D_m = 0`).**  In `CharZero`, feeding the
closed `hderiv` into mac-mini's `smallRootFactor_coeff0_eq_of_derivative_vanishes'`: the small-root
factor's constant coefficient is `P.coeff 0 = −t·r₀` (`r₀ = R.coeff 0`).  No exp, log, Puiseux, or
Fredholm determinant — the entire multiplicative route now rests on exactly the frame moment vanishing. -/
theorem smallRootFactor_coeff0_closed [CharZero F] (R : Polynomial F) (M : ℕ) (hM : 1 ≤ M)
    (hvanish : ∀ m : ℕ, 1 ≤ m →
      ((Polynomial.aeval (HahnSeries.single (1 : ℤ) (1 : F)) R) ^ m).coeff ((M : ℤ) * m) = 0) :
    (smallRootFactor R M).coeff 0
      = - PowerSeries.X * (algebraMap F (PowerSeries F)) (R.coeff 0) :=
  GMC2DvdKUnitOrigin.smallRootFactor_coeff0_eq_of_derivative_vanishes' R M hM
    (hderiv_concrete R M hvanish)

end GMC2DvdKConnectorPfr

#print axioms GMC2DvdKConnectorPfr.hderiv_concrete
#print axioms GMC2DvdKConnectorPfr.smallRootFactor_coeff0_closed
