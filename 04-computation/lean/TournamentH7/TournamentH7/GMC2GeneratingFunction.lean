import TournamentH7.GMC2ConstantTermRelations

/-!
# The DvdK generating function is trivial when all constant terms vanish

The DvdK constant terms are `D_m = aeval c (constantTermRelation q m)` (codex's Check A,
`GMC2LaurentShiftCheckA.constantCoeff_pow_eq_aeval_constantTermRelation`), and `D_0 = 1`.  The DvdK1
hypothesis being refuted is `D_m = 0` for all `m ≥ 1`, which makes the generating function
`F(t) = ∑_m D_m tᵐ` the constant power series `1`.  This is the elementary "the generating function is
trivial" step that both the additive (`b = F(t) = 1`) and multiplicative (`t·Π'/Π = F(t) = 1 ⟹ Π = c·t`)
endgames consume after the deep residue/Weierstrass identity `∑_{small} = F(t)` / `P.coeff 0 ~ F`.
-/

open MvPolynomial

namespace GMC2GeneratingFunction

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- `D_0 = aeval c (constantTermRelation q 0) = 1`: the only size-`0` composition is the zero
composition, which is balanced with multinomial coefficient `1` and empty coefficient product `1`. -/
theorem aeval_constantTermRelation_zero (q : ι → ℤ) (c : ι → ℂ) :
    MvPolynomial.aeval c (GMC2ConstantTermRelations.constantTermRelation q 0) = 1 := by
  rw [GMC2ConstantTermRelations.aeval_constantTermRelation]
  rw [Finset.piAntidiag_zero]
  rw [Finset.sum_singleton]
  have hbal : GMC2ConstantTermRelations.totalCharge q (0 : ι → ℕ) = 0 := by
    simp [GMC2ConstantTermRelations.totalCharge]
  rw [if_pos hbal]
  simp [Nat.multinomial]

/-- **The generating function is trivial.**  If every positive-power constant term
`aeval c (constantTermRelation q m)` (`m ≥ 1`) vanishes, then the DvdK generating function
`F(t) = ∑_m D_m tᵐ` equals the constant power series `1`. -/
theorem generatingFunction_eq_one (q : ι → ℤ) (c : ι → ℂ)
    (hvanish : ∀ m : ℕ, 1 ≤ m →
      MvPolynomial.aeval c (GMC2ConstantTermRelations.constantTermRelation q m) = 0) :
    (PowerSeries.mk (fun m : ℕ =>
        MvPolynomial.aeval c (GMC2ConstantTermRelations.constantTermRelation q m)) : PowerSeries ℂ)
      = 1 := by
  ext m
  rw [PowerSeries.coeff_mk]
  rcases Nat.eq_zero_or_pos m with hm | hm
  · subst hm
    rw [aeval_constantTermRelation_zero, PowerSeries.coeff_zero_eq_constantCoeff, map_one]
  · rw [hvanish m hm, PowerSeries.coeff_one, if_neg (by omega)]

end GMC2GeneratingFunction

#print axioms GMC2GeneratingFunction.aeval_constantTermRelation_zero
#print axioms GMC2GeneratingFunction.generatingFunction_eq_one
