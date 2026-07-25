import TournamentH7.GMC2ParityGuardrail

/-!
# Cubic polarization of the parity-fiber obstruction

`GMC2ParityGuardrail` proves that the even and odd fibers of the Boolean
three-cube agree on every Walsh observable of degree at most two.  This file
identifies the lost coordinate exactly: for an arbitrary integer-valued
statistic, the difference of its two parity-fiber integrals is simultaneously
its top Walsh pairing and its third mixed finite difference.

Thus a statistic separates the two fibers if and only if it sees the cubic
polarization.  This is a structural guardrail for moment compression, not a
missing premise of the already unconditional `GMC2Main.gmc2` proof.  In
particular, it is the three-cube version of THM-2237's missing top atom.

All identities are exact on eight points.  The file uses ordinary kernel
reduction only: no `sorry`, no `native_decide`, and no extra axiom.
-/

namespace GMC2.ParityPolarization

open ParityGuardrail

/-- The cube vertex with the displayed three Boolean coordinates. -/
def cubePoint (a b c : Bool) : Cube3 :=
  Fin.cases a (Fin.cases b (fun _ => c))

/-- Coordinate presentation of the three-cube.  This lets proofs expand the
eight vertices without asking the kernel to enumerate an entire function
type. -/
def cubeCoordinatesEquiv : (Bool × Bool × Bool) ≃ Cube3 where
  toFun abc := cubePoint abc.1 abc.2.1 abc.2.2
  invFun x := (x 0, x 1, x 2)
  left_inv abc := by
    rcases abc with ⟨a, b, c⟩
    rfl
  right_inv x := by
    funext i
    fin_cases i <;> rfl

/-- Pairing with the top Walsh character. -/
def topWalshPairing (f : Cube3 → ℤ) : ℤ :=
  ∑ x : Cube3, walshCharacter Finset.univ x * f x

/-- The third mixed finite difference, with the convention
`Δᵢ f = f|_(xᵢ=0) - f|_(xᵢ=1)`. -/
def cubicFiniteDifference (f : Cube3 → ℤ) : ℤ :=
    f (cubePoint false false false)
  - f (cubePoint true  false false)
  - f (cubePoint false true  false)
  - f (cubePoint false false true)
  + f (cubePoint true  true  false)
  + f (cubePoint true  false true)
  + f (cubePoint false true  true)
  - f (cubePoint true  true  true)

/-- Pointwise, even-fiber weight minus odd-fiber weight is exactly the top
Walsh character. -/
theorem parity_weight_difference_eq_top_character :
    ∀ x : Cube3,
      parityWeight false x - parityWeight true x =
        walshCharacter Finset.univ x := by
  intro x
  obtain ⟨⟨a, b, c⟩, rfl⟩ := cubeCoordinatesEquiv.surjective x
  cases a <;> cases b <;> cases c <;>
    decide

/-- The difference between the even- and odd-fiber integrals of any statistic
is its top Walsh coefficient. -/
theorem parity_expectation_difference_eq_top_walsh (f : Cube3 → ℤ) :
    parityExpectation false f - parityExpectation true f =
      topWalshPairing f := by
  unfold parityExpectation topWalshPairing
  rw [← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro x _
  rw [← sub_mul, parity_weight_difference_eq_top_character]

/-- On the three-cube, the top Walsh pairing is the third mixed finite
difference. -/
theorem top_walsh_eq_cubic_finite_difference (f : Cube3 → ℤ) :
    topWalshPairing f = cubicFiniteDifference f := by
  unfold topWalshPairing
  rw [← cubeCoordinatesEquiv.sum_comp
    (fun x => walshCharacter Finset.univ x * f x)]
  simp [Fintype.sum_prod_type, cubicFiniteDifference,
    cubeCoordinatesEquiv, cubePoint, walshCharacter, bitSign,
    Fin.prod_univ_succ]
  ring

/-- Exact polarization identity: parity-fiber separation is entirely cubic. -/
theorem parity_expectation_difference_eq_cubic (f : Cube3 → ℤ) :
    parityExpectation false f - parityExpectation true f =
      cubicFiniteDifference f := by
  rw [parity_expectation_difference_eq_top_walsh,
    top_walsh_eq_cubic_finite_difference]

/-- A statistic has equal integrals on the two parity fibers if and only if
its third mixed finite difference vanishes. -/
theorem parity_expectations_agree_iff_cubic_vanishes (f : Cube3 → ℤ) :
    parityExpectation false f = parityExpectation true f ↔
      cubicFiniteDifference f = 0 := by
  rw [← sub_eq_zero, parity_expectation_difference_eq_cubic]

/-- Equivalently, a statistic separates the two parity fibers exactly when
its cubic polarization is nonzero. -/
theorem parity_expectations_differ_iff_cubic_nonzero (f : Cube3 → ℤ) :
    parityExpectation false f ≠ parityExpectation true f ↔
      cubicFiniteDifference f ≠ 0 := by
  rw [ne_eq, ne_eq, parity_expectations_agree_iff_cubic_vanishes]

end GMC2.ParityPolarization

#print axioms GMC2.ParityPolarization.parity_weight_difference_eq_top_character
#print axioms GMC2.ParityPolarization.parity_expectation_difference_eq_top_walsh
#print axioms GMC2.ParityPolarization.top_walsh_eq_cubic_finite_difference
#print axioms GMC2.ParityPolarization.parity_expectation_difference_eq_cubic
#print axioms GMC2.ParityPolarization.parity_expectations_agree_iff_cubic_vanishes
#print axioms GMC2.ParityPolarization.parity_expectations_differ_iff_cubic_nonzero
