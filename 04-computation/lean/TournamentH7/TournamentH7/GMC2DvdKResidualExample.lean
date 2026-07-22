import TournamentH7.GMC2DvdKMonomialCertificate

/-!
# A residual coincident-channel support with DvdK discharged by a two-mass certificate

`{-2,-1,1,2}` is the paradigm support the unique-channel bypass cannot reach: the involution
`u ↦ -1/u` (`f(-1/u) = -f(u)`) pairs balanced compositions, so every mass carries ≥2 of them and no
single power certifies non-vanishing coefficient-independently.  Yet DvdK is discharged elementarily
by the **two-mass monomial certificate**

  `12 · (X₁X₂)² = (3 X₀X₃ + 9 X₁X₂) · CT(f²) − CT(f⁴)`,

so on the coefficient torus `CT(f²)` and `CT(f⁴)` cannot both vanish (`CT(f²)=0 ⟹ CT(f⁴) = -12 c₁²c₂²
≠ 0`).  This is `dvdk1_neg2_neg1_1_2` below — DvdK1 for `{-2,-1,1,2}`, kernel-pure, no Galois input.
The 14 residual straddling supports of size 3–4 in `[-4,4]` each carry such a certificate (degree ≤ 6,
verified in `04-computation/dvdk_monomial_certificate_residual_boxeph_S231.py`); this file mechanizes
the representative.
-/

open MvPolynomial Finset GMC2ConstantTermRelations

namespace GMC2DvdKResidualExample

/-- The charge vector of the paradigm residual support. -/
def q4 : Fin 4 → ℤ := ![(-2 : ℤ), -1, 1, 2]

/-- `CT(f²) = 2(c₀c₃ + c₁c₂)` — two coincident channels `(±2)` and `(±1)`. -/
theorem ct2_eval (c : Fin 4 → ℂ) :
    aeval c (constantTermRelation q4 2) = 2 * (c 0 * c 3) + 2 * (c 1 * c 2) := by
  rw [aeval_constantTermRelation, ← Finset.sum_filter]
  rw [show (Finset.piAntidiag (Finset.univ : Finset (Fin 4)) 2).filter
        (fun r => totalCharge q4 r = 0) = {![1,0,0,1], ![0,1,1,0]} from by decide]
  rw [Finset.sum_pair (by decide)]
  simp only [Fin.prod_univ_four, Fin.sum_univ_four, Nat.multinomial,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val, Fin.isValue]
  norm_num

/-- `CT(f⁴) = 6c₀²c₃² + 24c₀c₁c₂c₃ + 6c₁²c₂²`. -/
theorem ct4_eval (c : Fin 4 → ℂ) :
    aeval c (constantTermRelation q4 4)
      = 6 * (c 0 ^ 2 * c 3 ^ 2) + 24 * (c 0 * c 1 * c 2 * c 3) + 6 * (c 1 ^ 2 * c 2 ^ 2) := by
  rw [aeval_constantTermRelation, ← Finset.sum_filter]
  rw [show (Finset.piAntidiag (Finset.univ : Finset (Fin 4)) 4).filter
        (fun r => totalCharge q4 r = 0) = {![2,0,0,2], ![1,1,1,1], ![0,2,2,0]} from by decide]
  rw [Finset.sum_insert (by decide), Finset.sum_pair (by decide)]
  simp only [Fin.prod_univ_four, Fin.sum_univ_four, Nat.multinomial,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val, Fin.isValue]
  norm_num
  ring

/-- **DvdK1 for the residual support `{-2,-1,1,2}`** — no DvdK premise, no Galois.  The two-mass
certificate `{2,4}`: either `CT(f²) ≠ 0`, or `CT(f²)=0` forces `CT(f⁴) = -12 c₁²c₂² ≠ 0`. -/
theorem dvdk1_neg2_neg1_1_2 (c : Fin 4 → ℂ) (hc : ∀ i, c i ≠ 0) :
    ∃ m : ℕ, 1 ≤ m ∧ aeval c (constantTermRelation q4 m) ≠ 0 := by
  by_cases h2 : aeval c (constantTermRelation q4 2) = 0
  · -- CT(f²) = 0 ⟹ CT(f⁴) = -12 c₁²c₂² ≠ 0, so m = 4 works
    rw [ct2_eval] at h2
    refine ⟨4, by norm_num, ?_⟩
    rw [ct4_eval]
    have hE : c 0 * c 3 + c 1 * c 2 = 0 := by linear_combination h2 / 2
    have hval : 6 * (c 0 ^ 2 * c 3 ^ 2) + 24 * (c 0 * c 1 * c 2 * c 3) + 6 * (c 1 ^ 2 * c 2 ^ 2)
        = -12 * (c 1 ^ 2 * c 2 ^ 2) := by
      linear_combination (6 * (c 0 * c 3) + 18 * (c 1 * c 2)) * hE
    rw [hval]
    exact mul_ne_zero (by norm_num)
      (mul_ne_zero (pow_ne_zero 2 (hc 1)) (pow_ne_zero 2 (hc 2)))
  · -- CT(f²) ≠ 0, so m = 2 works
    exact ⟨2, by norm_num, h2⟩

end GMC2DvdKResidualExample

#print axioms GMC2DvdKResidualExample.dvdk1_neg2_neg1_1_2
