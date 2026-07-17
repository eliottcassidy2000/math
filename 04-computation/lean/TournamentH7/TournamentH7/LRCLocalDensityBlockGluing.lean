/-
  TournamentH7.LRCLocalDensityBlockGluing

  Algebraic formalization core of THM-933, the sharp local-density
  block-gluing theorem (codex-2026-07-16-S21, HYP-7152).

  The geometric input to THM-933 says that every component J of an earlier
  survivor retains at least

      density * length(J) - discrepancy

  inside the next certified block.  This file proves the reusable pieces
  around that input:

  * a bounded centered primitive gives the one-interval discrepancy inequality,
    with equality when an interval joins an upper extremum to a lower extremum;
  * summing the local inequality over components pays exactly card * q, and a
    component cap turns this into M * q;
  * if each newly deleted tooth raises component count by at most one, N teeth
    leave at most N components;
  * any sequence of nonnegative-density recurrence steps is bounded below by
    the recursive product-minus-weighted-debt ledger.

  It also kernel-checks the R = 7 LRC(14) arithmetic and the exact three-block
  certificate from the THM-933 referee.  No native_decide and no sorry.
-/

import Mathlib.Tactic

open scoped BigOperators

namespace LRC14.LocalDensityBlockGluing

section PrimitiveBridge

variable {Point : Type*}

/-- The analytic core of THM-933 Lemma 1.  Once retained interval measure is
represented as a difference of a bounded centered primitive, its oscillation
is a valid additive discrepancy loss. -/
theorem primitive_interval_lower
    (primitive : Point → ℝ) (startPoint endPoint : Point)
    (density intervalLength retained lowerValue upperValue : ℝ)
    (hlower : ∀ point, lowerValue ≤ primitive point)
    (hupper : ∀ point, primitive point ≤ upperValue)
    (hidentity :
      retained - density * intervalLength =
        primitive endPoint - primitive startPoint) :
    density * intervalLength - (upperValue - lowerValue) ≤ retained := by
  have hstart := hupper startPoint
  have hend := hlower endPoint
  linarith

/-- The primitive loss is sharp when the oriented interval starts at an upper
extremum and ends at a lower extremum. -/
theorem primitive_interval_sharp
    (primitive : Point → ℝ) (startPoint endPoint : Point)
    (density intervalLength retained lowerValue upperValue : ℝ)
    (hstart : primitive startPoint = upperValue)
    (hend : primitive endPoint = lowerValue)
    (hidentity :
      retained - density * intervalLength =
        primitive endPoint - primitive startPoint) :
    retained = density * intervalLength - (upperValue - lowerValue) := by
  rw [hstart, hend] at hidentity
  linarith

/-- An attained fixed-scale deficit lies below the primitive oscillation.
Taking the supremum over scales and using `primitive_interval_sharp` is the
paper proof of `q = sup_ell ell * (delta - eta ell)`. -/
theorem fixedScale_deficit_le_discrepancy
    (primitive : Point → ℝ) (startPoint endPoint : Point)
    (density eta intervalLength retained lowerValue upperValue : ℝ)
    (hlower : ∀ point, lowerValue ≤ primitive point)
    (hupper : ∀ point, primitive point ≤ upperValue)
    (heta : retained = eta * intervalLength)
    (hidentity :
      retained - density * intervalLength =
        primitive endPoint - primitive startPoint) :
    intervalLength * (density - eta) ≤ upperValue - lowerValue := by
  have hstart := hupper startPoint
  have hend := hlower endPoint
  have hprimitive :
      primitive startPoint - primitive endPoint ≤ upperValue - lowerValue :=
    sub_le_sub hstart hend
  calc
    intervalLength * (density - eta)
        = density * intervalLength - eta * intervalLength := by ring
    _ = density * intervalLength - retained := by rw [heta]
    _ = primitive startPoint - primitive endPoint := by linarith
    _ ≤ upperValue - lowerValue := hprimitive

/-- At a fixed-scale minimizer whose arc joins primitive extrema, the
fixed-scale deficit equals the primitive discrepancy exactly. -/
theorem fixedScale_extremizer_eq_discrepancy
    (primitive : Point → ℝ) (startPoint endPoint : Point)
    (density eta intervalLength retained lowerValue upperValue : ℝ)
    (heta : retained = eta * intervalLength)
    (hstart : primitive startPoint = upperValue)
    (hend : primitive endPoint = lowerValue)
    (hidentity :
      retained - density * intervalLength =
        primitive endPoint - primitive startPoint) :
    intervalLength * (density - eta) = upperValue - lowerValue := by
  have hsharp := primitive_interval_sharp primitive startPoint endPoint
    density intervalLength retained lowerValue upperValue hstart hend hidentity
  calc
    intervalLength * (density - eta)
        = density * intervalLength - eta * intervalLength := by ring
    _ = density * intervalLength - retained := by rw [heta]
    _ = upperValue - lowerValue := by linarith

end PrimitiveBridge

section ComponentSum

variable {ι : Type*}

/-- Summing a one-component local-density certificate pays one discrepancy
unit per component.  This is the algebraic content of THM-933 Lemma 2. -/
theorem local_to_component_sum
    (components : Finset ι) (length kept : ι → ℝ) (density discrepancy : ℝ)
    (hlocal : ∀ component ∈ components,
      density * length component - discrepancy ≤ kept component) :
    density * (∑ component ∈ components, length component)
        - (components.card : ℝ) * discrepancy
      ≤ ∑ component ∈ components, kept component := by
  calc
    density * (∑ component ∈ components, length component)
          - (components.card : ℝ) * discrepancy
        = ∑ component ∈ components,
            (density * length component - discrepancy) := by
              simp only [Finset.sum_sub_distrib, Finset.sum_const,
                nsmul_eq_mul, Finset.mul_sum]
    _ ≤ ∑ component ∈ components, kept component :=
      Finset.sum_le_sum fun component hcomponent => hlocal component hcomponent

/-- If the survivor has at most `complexity` components and discrepancy is
nonnegative, replace the exact cardinality debt by the advertised block
complexity debt. -/
theorem local_to_complexity_sum
    (components : Finset ι) (length kept : ι → ℝ)
    (density discrepancy complexity : ℝ)
    (hlocal : ∀ component ∈ components,
      density * length component - discrepancy ≤ kept component)
    (hcard : (components.card : ℝ) ≤ complexity)
    (hdiscrepancy : 0 ≤ discrepancy) :
    density * (∑ component ∈ components, length component)
        - complexity * discrepancy
      ≤ ∑ component ∈ components, kept component := by
  calc
    density * (∑ component ∈ components, length component)
          - complexity * discrepancy
        ≤ density * (∑ component ∈ components, length component)
          - (components.card : ℝ) * discrepancy := by
            exact sub_le_sub_left
              (mul_le_mul_of_nonneg_right hcard hdiscrepancy) _
    _ ≤ ∑ component ∈ components, kept component :=
      local_to_component_sum components length kept density discrepancy hlocal

/-- Fixed-scale sampling form (Opus S333 G1, after each component is tiled by
length-`ell` test intervals): summing the component losses pays
`eta * card * ell`.  THM-933's primitive form and this fixed-scale form are
dual certificate interfaces. -/
theorem fixedScale_sampling_sum
    (components : Finset ι) (length kept : ι → ℝ) (eta ell : ℝ)
    (hlocal : ∀ component ∈ components,
      eta * (length component - ell) ≤ kept component) :
      eta * ((∑ component ∈ components, length component)
        - (components.card : ℝ) * ell)
      ≤ ∑ component ∈ components, kept component := by
  have hsumConst :
      (∑ _component ∈ components, ell) = (components.card : ℝ) * ell := by
    rw [Finset.sum_const, nsmul_eq_mul]
  calc
    eta * ((∑ component ∈ components, length component)
          - (components.card : ℝ) * ell)
        = eta * ((∑ component ∈ components, length component)
            - ∑ _component ∈ components, ell) := by rw [hsumConst]
    _ = eta * (∑ component ∈ components,
          (length component - ell)) := by rw [Finset.sum_sub_distrib]
    _ = ∑ component ∈ components,
          eta * (length component - ell) := by rw [Finset.mul_sum]
    _ ≤ ∑ component ∈ components, kept component :=
      Finset.sum_le_sum fun component hcomponent => hlocal component hcomponent

/-- Opus's `eta * mass - loss` ledger is a valid weakening of the sharper
`eta * (mass - loss)` bound whenever `eta ≤ 1` and the loss is nonnegative. -/
theorem fixedScale_weaker_loss
    (mass boundaryLoss eta : ℝ) (heta : eta ≤ 1) (hloss : 0 ≤ boundaryLoss) :
    eta * mass - boundaryLoss ≤ eta * (mass - boundaryLoss) := by
  have hnonnegative : 0 ≤ (1 - eta) * boundaryLoss :=
    mul_nonneg (sub_nonneg.mpr heta) hloss
  calc
    eta * mass - boundaryLoss
        = eta * (mass - boundaryLoss) - (1 - eta) * boundaryLoss := by ring
    _ ≤ eta * (mass - boundaryLoss) := sub_le_self _ hnonnegative

end ComponentSum

section ComponentCap

/-- Abstract geometric ledger behind THM-933 Lemma 3.  The complement of one
deleted tooth has at most one component; if every additional tooth raises the
count by at most one, the complement of `toothCount` teeth has at most
`toothCount` components. -/
theorem component_count_le_tooth_count
    (componentCount : ℕ → ℕ)
    (hfirst : componentCount 1 ≤ 1)
    (hstep : ∀ toothCount, 1 ≤ toothCount →
      componentCount (toothCount + 1) ≤ componentCount toothCount + 1)
    (toothCount : ℕ) :
    1 ≤ toothCount → componentCount toothCount ≤ toothCount := by
  induction toothCount with
  | zero => omega
  | succ previousCount inductionHypothesis =>
      intro _
      by_cases hzero : previousCount = 0
      · subst previousCount
        simpa using hfirst
      · have hpositive : 1 ≤ previousCount := Nat.one_le_iff_ne_zero.mpr hzero
        calc
          componentCount (previousCount + 1)
              ≤ componentCount previousCount + 1 :=
                hstep previousCount hpositive
          _ ≤ previousCount + 1 :=
            Nat.add_le_add_right (inductionHypothesis hpositive) 1

end ComponentCap

section Recurrence

/-- Recursive lower ledger: multiply the previous lower bound by the next
local density and subtract the next boundary debt. -/
def lowerBound (initial : ℝ) (density debt : ℕ → ℝ) : ℕ → ℝ
  | 0 => initial
  | n + 1 => density n * lowerBound initial density debt n - debt n

/-- Recursive accumulated debt.  Earlier debts acquire every later density
factor, exactly as in the closed THM-933 formula. -/
def weightedDebt (density debt : ℕ → ℝ) : ℕ → ℝ
  | 0 => 0
  | n + 1 => density n * weightedDebt density debt n + debt n

/-- Product of all density factors used through stage `n`. -/
def densityProduct (density : ℕ → ℝ) (n : ℕ) : ℝ :=
  ∏ k ∈ Finset.range n, density k

/-- The recursive ledger is exactly product times the initial mass minus the
later-density-weighted debt. -/
theorem lowerBound_eq_product_sub_weightedDebt
    (initial : ℝ) (density debt : ℕ → ℝ) (n : ℕ) :
    lowerBound initial density debt n =
      densityProduct density n * initial - weightedDebt density debt n := by
  induction n with
  | zero => simp [lowerBound, densityProduct, weightedDebt]
  | succ n ih =>
      rw [lowerBound, weightedDebt, ih]
      simp only [densityProduct, Finset.prod_range_succ]
      ring

/-- Closed suffix-product form of the debt ledger.  The debt charged at stage
`k` is multiplied by every density factor from `k+1` through `n-1`. -/
theorem weightedDebt_eq_suffix_sum
    (density debt : ℕ → ℝ) (n : ℕ) :
    weightedDebt density debt n =
      ∑ k ∈ Finset.Ico 0 n,
        debt k * ∏ j ∈ Finset.Ico (k + 1) n, density j := by
  induction n with
  | zero => simp [weightedDebt]
  | succ n ih =>
      rw [weightedDebt, ih, Finset.sum_Ico_succ_top (Nat.zero_le n),
        Finset.Ico_self, Finset.prod_empty, mul_one]
      refine congr_arg (· + debt n) ?_
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro k hk
      rw [Finset.prod_Ico_succ_top (by
        have := Finset.mem_Ico.mp hk
        omega)]
      ring

/-- Fully unrolled product-minus-suffix-debt formula used in THM-933 (BG). -/
theorem lowerBound_eq_closed
    (initial : ℝ) (density debt : ℕ → ℝ) (n : ℕ) :
    lowerBound initial density debt n =
      densityProduct density n * initial
        - ∑ k ∈ Finset.Ico 0 n,
            debt k * ∏ j ∈ Finset.Ico (k + 1) n, density j := by
  rw [lowerBound_eq_product_sub_weightedDebt,
    weightedDebt_eq_suffix_sum]

/-- Soundness of the entire gluing recurrence.  Nonnegative density factors
allow every previously proved lower bound to be multiplied forward. -/
theorem lowerBound_le_actual
    (initial : ℝ) (density debt actual : ℕ → ℝ)
    (hinitial : initial ≤ actual 0)
    (hdensity : ∀ n, 0 ≤ density n)
    (hstep : ∀ n, density n * actual n - debt n ≤ actual (n + 1)) :
    ∀ n, lowerBound initial density debt n ≤ actual n := by
  intro n
  induction n with
  | zero => simpa [lowerBound] using hinitial
  | succ n ih =>
      calc
        lowerBound initial density debt (n + 1)
            = density n * lowerBound initial density debt n - debt n := rfl
        _ ≤ density n * actual n - debt n := by
          exact sub_le_sub_right (mul_le_mul_of_nonneg_left ih (hdensity n)) _
        _ ≤ actual (n + 1) := hstep n

/-- The explicit two-transition / three-block form consumed by THM-933:
`d₁d₂d₃ - (e₂d₃+e₃)`. -/
theorem three_block_gluing
    (d₁ d₂ d₃ e₂ e₃ w₁ w₂ w₃ : ℝ)
    (hd₂ : 0 ≤ d₂) (hd₃ : 0 ≤ d₃)
    (h₁ : d₁ ≤ w₁)
    (h₂ : d₂ * w₁ - e₂ ≤ w₂)
    (h₃ : d₃ * w₂ - e₃ ≤ w₃) :
    d₁ * d₂ * d₃ - (e₂ * d₃ + e₃) ≤ w₃ := by
  have h₂' : d₂ * d₁ - e₂ ≤ w₂ := by
    calc
      d₂ * d₁ - e₂ ≤ d₂ * w₁ - e₂ := by
        exact sub_le_sub_right (mul_le_mul_of_nonneg_left h₁ hd₂) _
      _ ≤ w₂ := h₂
  calc
    d₁ * d₂ * d₃ - (e₂ * d₃ + e₃) = d₃ * (d₂ * d₁ - e₂) - e₃ := by ring
    _ ≤ d₃ * w₂ - e₃ := by
      exact sub_le_sub_right (mul_le_mul_of_nonneg_left h₂' hd₃) _
    _ ≤ w₃ := h₃

end Recurrence

section ExactArithmetic

/-- Exact arithmetic behind the improved `R = 7` LRC(14) corollary. -/
theorem six_pow_twelve_gt_seven_pow_eleven :
    7 ^ 11 < 6 ^ 12 := by norm_num

/-- The explicit positive measure floor in THM-933 equation (10). -/
theorem ratio_seven_margin_pos :
    (0 : ℚ) < 199455593 / 13841287201 := by norm_num

/-- Kernel check of the three-block product-minus-debt identity. -/
theorem exact_three_block_ledger :
    (53 / 105 : ℚ) * (3 / 5) * (386 / 735)
        - (193 / 8575 + 193 / 6174)
      = 81253 / 771750 := by norm_num

/-- Kernel check of the stronger ledger using exact prefix component counts
`10` and `132` instead of tooth caps `15` and `375`. -/
theorem exact_three_block_component_ledger :
    (53 / 105 : ℚ) * (3 / 5) * (386 / 735)
        - (386 / 25725 + 4246 / 385875)
      = 7334 / 55125 := by norm_num

/-- The exact THM-933 block-gluing floor is positive. -/
theorem exact_three_block_margin_pos :
    (0 : ℚ) < 81253 / 771750 := by norm_num

/-- Exact arithmetic for the sharp fixed-scale Opus S333 7+6 witness. -/
theorem opus_fixedScale_two_block_ledger :
    (274025490881738650 / 1001359472502594621 : ℚ)
        * (12208485893419843 / 38882590600758450 - 1724 / 45000)
      = 60354211840721383388269695262412
        / 800043501647462161192289496375975 := by norm_num

/-- The fixed-scale 7+6 composition margin is strictly positive. -/
theorem opus_fixedScale_two_block_margin_pos :
    (0 : ℚ) < 60354211840721383388269695262412
      / 800043501647462161192289496375975 := by norm_num

/-- The direct exact safe measure dominates the proved block-gluing floor. -/
theorem exact_measure_dominates_ledger :
    (81253 / 771750 : ℚ) ≤ 55063 / 330750 := by norm_num

end ExactArithmetic

/-! ## Axiom audit -/

#print axioms local_to_component_sum
#print axioms local_to_complexity_sum
#print axioms fixedScale_sampling_sum
#print axioms primitive_interval_lower
#print axioms primitive_interval_sharp
#print axioms fixedScale_deficit_le_discrepancy
#print axioms fixedScale_extremizer_eq_discrepancy
#print axioms fixedScale_weaker_loss
#print axioms component_count_le_tooth_count
#print axioms lowerBound_eq_product_sub_weightedDebt
#print axioms weightedDebt_eq_suffix_sum
#print axioms lowerBound_eq_closed
#print axioms lowerBound_le_actual
#print axioms three_block_gluing
#print axioms six_pow_twelve_gt_seven_pow_eleven
#print axioms ratio_seven_margin_pos
#print axioms exact_three_block_ledger
#print axioms exact_three_block_component_ledger
#print axioms exact_three_block_margin_pos
#print axioms opus_fixedScale_two_block_ledger
#print axioms opus_fixedScale_two_block_margin_pos
#print axioms exact_measure_dominates_ledger

end LRC14.LocalDensityBlockGluing
