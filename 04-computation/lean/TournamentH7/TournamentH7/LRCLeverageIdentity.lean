import Mathlib

/-!
# The leverage identity (THM-930/THM-935 algebraic core)

kind-pasteur-S128c37.  The Bonferroni truncation error in closed form, at the
level of an arbitrary finite cell decomposition — no measure theory:

Given cells `c` with nonnegative weights `w c` and depths `D c`, and the
binomial moments `S k = ∑_c C(D c, k) * w c`, the alternating truncation obeys

  `∑_{k ≤ m} (-1)^k S k  =  μ₀ + (-1)^m ∑_{D c ≥ 1} C(D c - 1, m) * w c`,

where `μ₀ = ∑_{D c = 0} w c` is the good mass.  The engine is the partial
alternating row sum `∑_{k ≤ m} (-1)^k C(d,k) = (-1)^m C(d-1,m)` (Pascal
telescoping).  Corollaries: the two-sided Bonferroni inequalities with exact
error, and THE CERTIFICATE THEOREM — `0 < B_m (m odd) → 0 < μ₀` — the
statement the level-5 program rests on.  Constants anchors tie the
exact-support weights `E_s` (24/343, 24/49, -2/7), the equilibrium
`2052/16807`, the depth-13 leverage `C(12,5) = 792`, and the kill threshold
to `LRCB5RelationBudget`'s interface.

The remaining analytic step for a concrete packet is only the sweep encoding:
that the circle decomposes into finitely many cells of these depths — the
identity itself is kernel-pure algebra.
-/

namespace LonelyRunner
namespace LRCLeverageIdentity

open Finset

/-- Partial alternating row sums of the zeroth binomial row are all `1`. -/
theorem alternating_partial_choose_zero (m : ℕ) :
    ∑ k ∈ range (m + 1), ((-1 : ℚ)) ^ k * ((0 : ℕ).choose k : ℚ) = 1 := by
  induction m with
  | zero => simp
  | succ m ih =>
      rw [sum_range_succ, ih]
      have hz : (0 : ℕ).choose (m + 1) = 0 :=
        Nat.choose_eq_zero_of_lt (Nat.succ_pos m)
      rw [hz]
      push_cast
      ring

/-- THE PARTIAL ALTERNATING ROW SUM (Pascal telescoping): for `1 ≤ d`,
    `∑_{k=0}^{m} (-1)^k C(d,k) = (-1)^m C(d-1,m)`. -/
theorem alternating_partial_choose (d : ℕ) (hd : 1 ≤ d) (m : ℕ) :
    ∑ k ∈ range (m + 1), ((-1 : ℚ)) ^ k * (d.choose k : ℚ)
      = ((-1 : ℚ)) ^ m * ((d - 1).choose m : ℚ) := by
  induction m with
  | zero => simp
  | succ m ih =>
      have pas : d.choose (m + 1) = (d - 1).choose m + (d - 1).choose (m + 1) := by
        conv_lhs => rw [← Nat.sub_add_cancel hd]
        exact Nat.choose_succ_succ (d - 1) m
      rw [sum_range_succ, ih]
      push_cast [pas]
      ring

variable {N : ℕ}

/-- The `k`-th binomial moment of a weighted depth profile
    (`S_k` of the depth spectrum). -/
def binomMoment (w : Fin N → ℚ) (D : Fin N → ℕ) (k : ℕ) : ℚ :=
  ∑ c, ((D c).choose k : ℚ) * w c

/-- The good mass: total weight at depth `0`. -/
def goodMass (w : Fin N → ℚ) (D : Fin N → ℕ) : ℚ :=
  ∑ c ∈ univ.filter (fun c => D c = 0), w c

/-- The `m`-leveraged tail: `∑ C(D c - 1, m) * w c` over cells of depth ≥ 1.
    (Cells with `1 ≤ D c ≤ m` contribute `0` automatically.) -/
def leveragedTail (w : Fin N → ℚ) (D : Fin N → ℕ) (m : ℕ) : ℚ :=
  ∑ c ∈ univ.filter (fun c => 1 ≤ D c), ((D c - 1).choose m : ℚ) * w c

/-- **THE LEVERAGE IDENTITY**: the Bonferroni truncation error in closed form.
    `∑_{k ≤ m} (-1)^k S_k = μ₀ + (-1)^m · leveragedTail`. -/
theorem leverage_identity (w : Fin N → ℚ) (D : Fin N → ℕ) (m : ℕ) :
    ∑ k ∈ range (m + 1), ((-1 : ℚ)) ^ k * binomMoment w D k
      = goodMass w D + ((-1 : ℚ)) ^ m * leveragedTail w D m := by
  have swap : ∑ k ∈ range (m + 1), ((-1 : ℚ)) ^ k * binomMoment w D k
      = ∑ c, (∑ k ∈ range (m + 1), ((-1 : ℚ)) ^ k * ((D c).choose k : ℚ)) * w c := by
    unfold binomMoment
    calc
      ∑ k ∈ range (m + 1), ((-1 : ℚ)) ^ k * ∑ c, ((D c).choose k : ℚ) * w c
          = ∑ k ∈ range (m + 1), ∑ c,
              ((-1 : ℚ)) ^ k * (((D c).choose k : ℚ) * w c) := by
            refine Finset.sum_congr rfl fun k _ => ?_
            rw [Finset.mul_sum]
      _ = ∑ c, ∑ k ∈ range (m + 1),
              ((-1 : ℚ)) ^ k * (((D c).choose k : ℚ) * w c) := Finset.sum_comm
      _ = ∑ c, (∑ k ∈ range (m + 1), ((-1 : ℚ)) ^ k * ((D c).choose k : ℚ)) * w c := by
            refine Finset.sum_congr rfl fun c _ => ?_
            rw [Finset.sum_mul]
            exact Finset.sum_congr rfl fun k _ => (mul_assoc _ _ _).symm
  rw [swap]
  rw [← Finset.sum_filter_add_sum_filter_not
    (Finset.univ : Finset (Fin N)) (fun c => D c = 0)]
  unfold goodMass leveragedTail
  congr 1
  · refine Finset.sum_congr rfl fun c hc => ?_
    have h0 : D c = 0 := (Finset.mem_filter.mp hc).2
    rw [h0, alternating_partial_choose_zero, one_mul]
  · have hset : (Finset.univ.filter (fun c : Fin N => ¬ D c = 0))
        = Finset.univ.filter (fun c => 1 ≤ D c) := by
      apply Finset.filter_congr
      intro c _
      simp [Nat.one_le_iff_ne_zero]
    rw [hset, Finset.mul_sum]
    refine Finset.sum_congr rfl fun c hc => ?_
    have h1 : 1 ≤ D c := (Finset.mem_filter.mp hc).2
    rw [alternating_partial_choose (D c) h1 m, mul_assoc]

/-- Bonferroni, odd truncation: a lower bound with the error exact. -/
theorem bonferroni_odd_le (w : Fin N → ℚ) (D : Fin N → ℕ) (m : ℕ)
    (hm : Odd m) (hw : ∀ c, 0 ≤ w c) :
    ∑ k ∈ range (m + 1), ((-1 : ℚ)) ^ k * binomMoment w D k ≤ goodMass w D := by
  rw [leverage_identity]
  have htail : 0 ≤ leveragedTail w D m := by
    unfold leveragedTail
    exact Finset.sum_nonneg fun c _ =>
      mul_nonneg (Nat.cast_nonneg _) (hw c)
  have hsign : ((-1 : ℚ)) ^ m = -1 := Odd.neg_one_pow hm
  rw [hsign]
  linarith

/-- Bonferroni, even truncation: an upper bound with the error exact. -/
theorem bonferroni_even_ge (w : Fin N → ℚ) (D : Fin N → ℕ) (m : ℕ)
    (hm : Even m) (hw : ∀ c, 0 ≤ w c) :
    goodMass w D ≤ ∑ k ∈ range (m + 1), ((-1 : ℚ)) ^ k * binomMoment w D k := by
  rw [leverage_identity]
  have htail : 0 ≤ leveragedTail w D m := by
    unfold leveragedTail
    exact Finset.sum_nonneg fun c _ =>
      mul_nonneg (Nat.cast_nonneg _) (hw c)
  have hsign : ((-1 : ℚ)) ^ m = 1 := Even.neg_one_pow hm
  rw [hsign]
  linarith

/-- **THE CERTIFICATE THEOREM**: a positive odd-level Bonferroni truncation
    certifies positive good mass.  (`m = 5` is the level-5 wall's use.) -/
theorem goodMass_pos_of_bonferroni_pos (w : Fin N → ℚ) (D : Fin N → ℕ)
    (m : ℕ) (hm : Odd m) (hw : ∀ c, 0 ≤ w c)
    (hpos : 0 < ∑ k ∈ range (m + 1), ((-1 : ℚ)) ^ k * binomMoment w D k) :
    0 < goodMass w D :=
  lt_of_lt_of_le hpos (bonferroni_odd_le w D m hm hw)

/-! ## Constants anchors (THM-935's exact-support weights, THM-930's budget) -/

/-- `E₂ = ∑_{j≤3} (-1)^j C(11,j) (1/7)^j = 24/343`. -/
theorem E2_eq :
    ∑ j ∈ range 4, ((-1 : ℚ)) ^ j * ((11 : ℕ).choose j : ℚ) * (1 / 7) ^ j
      = 24 / 343 := by
  have c1 : (11 : ℕ).choose 1 = 11 := by decide
  have c2 : (11 : ℕ).choose 2 = 55 := by decide
  have c3 : (11 : ℕ).choose 3 = 165 := by decide
  simp only [Finset.sum_range_succ, Finset.sum_range_zero,
    Nat.choose_zero_right, c1, c2, c3]
  norm_num

/-- `E₃ = ∑_{j≤2} (-1)^j C(10,j) (1/7)^j = 24/49`. -/
theorem E3_eq :
    ∑ j ∈ range 3, ((-1 : ℚ)) ^ j * ((10 : ℕ).choose j : ℚ) * (1 / 7) ^ j
      = 24 / 49 := by
  have c1 : (10 : ℕ).choose 1 = 10 := by decide
  have c2 : (10 : ℕ).choose 2 = 45 := by decide
  simp only [Finset.sum_range_succ, Finset.sum_range_zero,
    Nat.choose_zero_right, c1, c2]
  norm_num

/-- `E₄ = ∑_{j≤1} (-1)^j C(9,j) (1/7)^j = -2/7`. -/
theorem E4_eq :
    ∑ j ∈ range 2, ((-1 : ℚ)) ^ j * ((9 : ℕ).choose j : ℚ) * (1 / 7) ^ j
      = -2 / 7 := by
  have c1 : (9 : ℕ).choose 1 = 9 := by decide
  simp only [Finset.sum_range_succ, Finset.sum_range_zero,
    Nat.choose_zero_right, c1]
  norm_num

/-- The level-5 equilibrium mass:
    `∑_{k≤5} (-1)^k C(13,k) (1/7)^k = 2052/16807`. -/
theorem B5_equilibrium_eq :
    ∑ k ∈ range 6, ((-1 : ℚ)) ^ k * ((13 : ℕ).choose k : ℚ) * (1 / 7) ^ k
      = 2052 / 16807 := by
  have c1 : (13 : ℕ).choose 1 = 13 := by decide
  have c2 : (13 : ℕ).choose 2 = 78 := by decide
  have c3 : (13 : ℕ).choose 3 = 286 := by decide
  have c4 : (13 : ℕ).choose 4 = 715 := by decide
  have c5 : (13 : ℕ).choose 5 = 1287 := by decide
  simp only [Finset.sum_range_succ, Finset.sum_range_zero,
    Nat.choose_zero_right, c1, c2, c3, c4, c5]
  norm_num

/-- The depth-13 leverage: `C(12,5) = 792`. -/
theorem depth13_leverage : (12 : ℕ).choose 5 = 792 := by decide

/-- The kill threshold `2052/16807/792 = 57/369754` in lowest terms: a
    full-depth atom of larger mass alone breaks the level-5 certificate. -/
theorem kill_threshold_eq : (2052 : ℚ) / 16807 / 792 = 57 / 369754 := by
  norm_num

/-! ## Axiom audit -/

#print axioms alternating_partial_choose
#print axioms leverage_identity
#print axioms goodMass_pos_of_bonferroni_pos
#print axioms B5_equilibrium_eq

end LRCLeverageIdentity
end LonelyRunner
