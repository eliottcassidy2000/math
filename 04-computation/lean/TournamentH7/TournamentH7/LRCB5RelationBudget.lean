import Mathlib

/-!
# The exact depth-five relation budget

THM-935 groups the level-five singular series by exact relation support.  Once
that analytic identity is supplied, its numerical consumer is completely
elementary: the equilibrium mass is `2052 / 16807`, while exact supports
`2,3,4,5` enter with signed coefficients `24/343`, `-24/49`, `-2/7`, and
`-1`.

This module checks both the resulting absolute-debt implication and the sharp
one-sided budget split at pair horizon `H = 30`.  It intentionally does not assert the
singular-series identity for the concrete discrete `B5`; that analytic bridge
remains a separate proof obligation.
-/

namespace LonelyRunner
namespace LRCB5RelationBudget

noncomputable section

/-- The independent/equilibrium depth-five mass for thirteen runners. -/
def equilibrium : ℝ := 2052 / 16807

/-- Absolute coefficient of the exact-support-two relation mass. -/
def pairWeight : ℝ := 24 / 343

/-- Absolute coefficient of the exact-support-three relation mass. -/
def tripleWeight : ℝ := 24 / 49

/-- Absolute coefficient of the exact-support-four relation mass. -/
def quadWeight : ℝ := 2 / 7

/-- The exact-support decomposition's signed algebraic model. -/
def relationModel (mass2 mass3 mass4 mass5 : ℝ) : ℝ :=
  equilibrium + pairWeight * mass2 - tripleWeight * mass3 -
    quadWeight * mass4 - mass5

/-- The sign-blind debt used to certify positivity uniformly. -/
def relationDebt (mass2 mass3 mass4 mass5 : ℝ) : ℝ :=
  pairWeight * |mass2| + tripleWeight * |mass3| +
    quadWeight * |mass4| + |mass5|

/-- The support-three/four/five portion of the absolute debt. -/
def higherRelationDebt (mass3 mass4 mass5 : ℝ) : ℝ :=
  tripleWeight * |mass3| + quadWeight * |mass4| + |mass5|

/-- The signed support-three/four/five contribution which actually harms the
lower bound.  Negative values help the depth-five model and consume no budget. -/
def harmfulHigherContribution (mass3 mass4 mass5 : ℝ) : ℝ :=
  tripleWeight * mass3 + quadWeight * mass4 + mass5

theorem equilibrium_pos : 0 < equilibrium := by
  norm_num [equilibrium]

theorem weights_nonnegative :
    0 ≤ pairWeight ∧ 0 ≤ tripleWeight ∧ 0 ≤ quadWeight := by
  norm_num [pairWeight, tripleWeight, quadWeight]

/-- Triangle inequality in the exact coefficient frame. -/
theorem equilibrium_sub_debt_le_model (mass2 mass3 mass4 mass5 : ℝ) :
    equilibrium - relationDebt mass2 mass3 mass4 mass5 ≤
      relationModel mass2 mass3 mass4 mass5 := by
  dsimp [equilibrium, relationDebt, relationModel, pairWeight, tripleWeight, quadWeight]
  nlinarith [neg_abs_le mass2, le_abs_self mass3, le_abs_self mass4,
    le_abs_self mass5]

/-- Any absolute relation debt below the equilibrium budget forces a positive
depth-five model. -/
theorem relationModel_pos_of_debt_lt (mass2 mass3 mass4 mass5 : ℝ)
    (hdebt : relationDebt mass2 mass3 mass4 mass5 < equilibrium) :
    0 < relationModel mass2 mass3 mass4 mass5 := by
  have hlower := equilibrium_sub_debt_le_model mass2 mass3 mass4 mass5
  linarith

/-- Contrapositive budget form: a failed depth-five certificate must spend at
least the full equilibrium mass in exact-support relation debt. -/
theorem equilibrium_le_debt_of_model_nonpos (mass2 mass3 mass4 mass5 : ℝ)
    (hmodel : relationModel mass2 mass3 mass4 mass5 ≤ 0) :
    equilibrium ≤ relationDebt mass2 mass3 mass4 mass5 := by
  have hlower := equilibrium_sub_debt_le_model mass2 mass3 mass4 mass5
  linarith

theorem relationDebt_eq_pair_add_higher (mass2 mass3 mass4 mass5 : ℝ) :
    relationDebt mass2 mass3 mass4 mass5 =
      pairWeight * |mass2| + higherRelationDebt mass3 mass4 mass5 := by
  simp [relationDebt, higherRelationDebt, add_assoc]

/-- THM-935's proved support-two tail at horizon `30` consumes strictly less
than one quarter of the exact equilibrium budget. -/
theorem pair_tail_horizon_thirty_lt_quarter :
    pairWeight * (13 / 30 : ℝ) < equilibrium / 4 := by
  norm_num [pairWeight, equilibrium]

/-- Exact pair-tail expenditure at horizon thirty. -/
theorem pair_tail_horizon_thirty_eq :
    pairWeight * (13 / 30 : ℝ) = 52 / 1715 := by
  norm_num [pairWeight]

/-- Exact higher-support budget remaining after the full proved pair tail. -/
theorem remaining_after_pair_tail_eq :
    equilibrium - pairWeight * (13 / 30 : ℝ) = 7712 / 84035 := by
  norm_num [equilibrium, pairWeight]

/-- The exact remainder is slightly stronger than the convenient
three-quarter socket. -/
theorem remaining_after_pair_tail_eq_three_quarters_add_margin :
    equilibrium - pairWeight * (13 / 30 : ℝ) =
      3 * equilibrium / 4 + 17 / 84035 := by
  norm_num [equilibrium, pairWeight]

/-- Sharp one-sided horizon-thirty consumer.  The pair coefficient is positive,
so only a lower bound on `mass2` is needed; the higher supports enter through
their signed harmful combination rather than through absolute values. -/
theorem relationModel_pos_of_signed_horizon_thirty_split
    (mass2 mass3 mass4 mass5 : ℝ)
    (hpair : -(13 / 30 : ℝ) ≤ mass2)
    (hhigher : harmfulHigherContribution mass3 mass4 mass5 < 7712 / 84035) :
    0 < relationModel mass2 mass3 mass4 mass5 := by
  have hweight : (0 : ℝ) ≤ pairWeight := weights_nonnegative.1
  have hpairWeighted : -(52 / 1715 : ℝ) ≤ pairWeight * mass2 := by
    have := mul_le_mul_of_nonneg_left hpair hweight
    norm_num [pairWeight] at this ⊢
    exact this
  dsimp [relationModel, harmfulHigherContribution, equilibrium, pairWeight,
    tripleWeight, quadWeight] at hpairWeighted hhigher ⊢
  linarith

/-- A direct socket for the remaining `T_s(H)`, `s=3,4,5`, estimates: the
proved pair-tail quarter plus a strict three-quarter higher-support bound
forces positivity. -/
theorem relationModel_pos_of_quarter_threeQuarter_split
    (mass2 mass3 mass4 mass5 : ℝ)
    (hpair : pairWeight * |mass2| ≤ equilibrium / 4)
    (hhigher : higherRelationDebt mass3 mass4 mass5 < 3 * equilibrium / 4) :
    0 < relationModel mass2 mass3 mass4 mass5 := by
  apply relationModel_pos_of_debt_lt
  rw [relationDebt_eq_pair_add_higher]
  linarith

/-- Strongest direct horizon-thirty consumer: the actual pair estimate leaves
the exact higher-support allowance `7712/84035`. -/
theorem relationModel_pos_of_exact_horizon_thirty_split
    (mass2 mass3 mass4 mass5 : ℝ)
    (hpair : |mass2| ≤ 13 / 30)
    (hhigher : higherRelationDebt mass3 mass4 mass5 < 7712 / 84035) :
    0 < relationModel mass2 mass3 mass4 mass5 := by
  apply relationModel_pos_of_debt_lt
  rw [relationDebt_eq_pair_add_higher]
  have hpairWeighted : pairWeight * |mass2| ≤ 52 / 1715 := by
    have hweight : (0 : ℝ) ≤ pairWeight := weights_nonnegative.1
    calc
      pairWeight * |mass2| ≤ pairWeight * (13 / 30) :=
        mul_le_mul_of_nonneg_left hpair hweight
      _ = 52 / 1715 := pair_tail_horizon_thirty_eq
  norm_num [equilibrium]
  linarith

/-! ## Axiom audit -/

#print axioms equilibrium_sub_debt_le_model
#print axioms relationModel_pos_of_debt_lt
#print axioms equilibrium_le_debt_of_model_nonpos
#print axioms pair_tail_horizon_thirty_lt_quarter
#print axioms relationModel_pos_of_signed_horizon_thirty_split
#print axioms relationModel_pos_of_quarter_threeQuarter_split
#print axioms relationModel_pos_of_exact_horizon_thirty_split

end
end LRCB5RelationBudget
end LonelyRunner
