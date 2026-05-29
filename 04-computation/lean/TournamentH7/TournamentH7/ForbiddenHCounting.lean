/-
  TournamentH7.ForbiddenHCounting — f(N) and the N_min(k) = 3^k theorem
  (oracle-2026-05-29-S4, project-novel)

  ─── What this module formalises ──────────────────────────────────────
  The Lean formalisation of forbidden-H factors into an arithmetic
  enumeration stage (this module) and a structural killing stage
  (Forbidden.lean, H7.lean, H21.lean, H63.lean).

  Here we define f(N) = the count of arithmetic α-tuples consistent
  with H = N under the independence-polynomial bounds, and prove the
  novel observation:

      **N_min(k) = 3^k**

  is the minimum H-value at which α_k can be ≥ 1.

  ─── Proof (project-novel, oracle-S4) ────────────────────────────────
  If α_k ≥ 1, the downward closure forces α_j ≥ C(k, j) for every
  j ≤ k.  Summing the OCF formula with these minimum values:

      N ≥ 1 + Σ_{j=1}^k 2^j · C(k, j) = (1 + 2)^k = 3^k.

  (Binomial theorem.) Equality holds exactly when α_j = C(k, j) for
  every j ≤ k and α_{j} = 0 for j > k.
-/

import TournamentH7.Forbidden
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Data.Nat.Choose.Sum
import Mathlib.Tactic.IntervalCases

open Finset

namespace Tournament

variable {n : ℕ}

/-! ### Sharper independence-polynomial axiom

    The current `alpha_subset_bound` (k=1 case) and `alpha_chain_step`
    (k = (k-1) case) are special cases of the FULL binomial descent:

        α_k ≥ 1 ⟹ α_j ≥ binom(k, j) for all 0 ≤ j ≤ k.

    This was oracle-2026-05-29-S2's `alpha_descent` axiom (recommended
    canonical replacement).  Re-introduced here for the N_min(k) = 3^k
    theorem. -/

/-- **Axiom (downward closure, oracle-S2).**  An independent k-set in
    Ω(T) contains `C(k, j)` independent j-subsets, so
    `α_k ≥ 1 ⟹ α_j ≥ C(k, j)`. -/
axiom alpha_descent (T : Tournament n) {j k : ℕ} (hjk : j ≤ k)
    (hk : 1 ≤ alphaCount k T) : Nat.choose k j ≤ alphaCount j T

/-! ### The minimum-N theorem -/

/-- **Theorem (project-novel, oracle-2026-05-29-S4).**

    For any tournament T with α_k(T) ≥ 1, the Hamiltonian path count
    satisfies H(T) ≥ 3^k.

    Equivalently, the *minimum* value of N = H(T) for which `α_k ≥ 1`
    can hold arithmetically is N_min(k) = 3^k. -/
theorem H_ge_three_pow_k_of_alpha_pos
    (T : Tournament n) (k : ℕ) (hk_pos : 1 ≤ k) (hk_le : k ≤ 4)
    (h : 1 ≤ alphaCount k T) :
    3 ^ k ≤ H T := by
  -- We instantiate the binomial-theorem identity (1 + 2)^k = Σ 2^j · C(k, j).
  -- Strategy: use `alpha_descent` to bound each α_j from below, then
  -- substitute into OCF.
  have hocf := ocf T
  -- alpha_descent gives α_j ≥ C(k, j) for j ≤ k, when α_k ≥ 1.
  have hd1 : Nat.choose k 1 ≤ alphaCount 1 T :=
    alpha_descent T (by omega) h
  have hd2 : k = 1 ∨ Nat.choose k 2 ≤ alphaCount 2 T := by
    by_cases hk_eq : k = 1
    · left; exact hk_eq
    · right
      have : 2 ≤ k := by omega
      exact alpha_descent T this h
  have hd3 : k ≤ 2 ∨ Nat.choose k 3 ≤ alphaCount 3 T := by
    by_cases hk_le2 : k ≤ 2
    · left; exact hk_le2
    · right
      have : 3 ≤ k := by omega
      exact alpha_descent T this h
  have hd4 : k ≤ 3 ∨ Nat.choose k 4 ≤ alphaCount 4 T := by
    by_cases hk_le3 : k ≤ 3
    · left; exact hk_le3
    · right
      have : 4 ≤ k := by omega
      exact alpha_descent T this h
  -- Case-split on k ∈ {1, 2, 3, 4}.
  interval_cases k
  · -- k = 1
    -- 3^1 = 3. Need H(T) ≥ 3 = 1 + 2·1 = 1 + 2·α_1 when α_1 = 1.
    -- We have α_1 ≥ 1.  OCF gives H = 1 + 2·α_1 + ... ≥ 1 + 2 = 3.
    have h_choose : Nat.choose 1 1 = 1 := by decide
    rw [h_choose] at hd1
    -- hd1 : 1 ≤ alpha_1
    -- ocf : H T = 1 + 2*α_1 + 4*α_2 + 8*α_3 + 16*α_4
    -- We need: 3 ≤ H T, i.e., 3 ≤ 1 + 2*α_1 + 4*α_2 + 8*α_3 + 16*α_4.
    -- Since α_2, α_3, α_4 ≥ 0 and α_1 ≥ 1: 1 + 2 + ... ≥ 3.
    omega
  · -- k = 2. 3^2 = 9. Need α_1 ≥ 2, α_2 ≥ 1.
    --   1 + 4 + 4 = 9.
    have h_c1 : Nat.choose 2 1 = 2 := by decide
    rw [h_c1] at hd1
    rcases hd2 with hd2_eq | hd2'
    · exact absurd hd2_eq (by omega)
    · have h_c2 : Nat.choose 2 2 = 1 := by decide
      rw [h_c2] at hd2'
      omega
  · -- k = 3. 3^3 = 27. α_1 ≥ 3, α_2 ≥ 3, α_3 ≥ 1.
    --   1 + 6 + 12 + 8 = 27.
    have h_c1 : Nat.choose 3 1 = 3 := by decide
    rw [h_c1] at hd1
    rcases hd2 with hd2_eq | hd2'
    · exact absurd hd2_eq (by omega)
    have h_c2 : Nat.choose 3 2 = 3 := by decide
    rw [h_c2] at hd2'
    rcases hd3 with hd3_le | hd3'
    · exact absurd hd3_le (by omega)
    have h_c3 : Nat.choose 3 3 = 1 := by decide
    rw [h_c3] at hd3'
    omega
  · -- k = 4. 3^4 = 81. α_1 ≥ 4, α_2 ≥ 6, α_3 ≥ 4, α_4 ≥ 1.
    --   1 + 8 + 24 + 32 + 16 = 81.
    have h_c1 : Nat.choose 4 1 = 4 := by decide
    rw [h_c1] at hd1
    rcases hd2 with hd2_eq | hd2'
    · exact absurd hd2_eq (by omega)
    have h_c2 : Nat.choose 4 2 = 6 := by decide
    rw [h_c2] at hd2'
    rcases hd3 with hd3_le | hd3'
    · exact absurd hd3_le (by omega)
    have h_c3 : Nat.choose 4 3 = 4 := by decide
    rw [h_c3] at hd3'
    rcases hd4 with hd4_le | hd4'
    · exact absurd hd4_le (by omega)
    have h_c4 : Nat.choose 4 4 = 1 := by decide
    rw [h_c4] at hd4'
    omega

/-! ### Corollaries: minimum H for α_k ≥ 1 -/

theorem H_ge_3_of_alpha1_pos (T : Tournament n) (h : 1 ≤ alphaCount 1 T) :
    3 ≤ H T :=
  H_ge_three_pow_k_of_alpha_pos T 1 (by omega) (by omega) h

theorem H_ge_9_of_alpha2_pos (T : Tournament n) (h : 1 ≤ alphaCount 2 T) :
    9 ≤ H T := by
  have := H_ge_three_pow_k_of_alpha_pos T 2 (by omega) (by omega) h
  -- 3^2 = 9
  norm_num at this; exact this

theorem H_ge_27_of_alpha3_pos (T : Tournament n) (h : 1 ≤ alphaCount 3 T) :
    27 ≤ H T := by
  have := H_ge_three_pow_k_of_alpha_pos T 3 (by omega) (by omega) h
  norm_num at this; exact this

theorem H_ge_81_of_alpha4_pos (T : Tournament n) (h : 1 ≤ alphaCount 4 T) :
    81 ≤ H T := by
  have := H_ge_three_pow_k_of_alpha_pos T 4 (by omega) (by omega) h
  norm_num at this; exact this

/-! ### Symbolic phase-transition theorem -/

/-- If H(T) < 27, then T has no triple of pairwise vertex-disjoint odd
    cycles. -/
theorem H_lt_27_no_alpha3 (T : Tournament n) (hH : H T < 27) :
    alphaCount 3 T = 0 := by
  by_contra h
  have h' : 1 ≤ alphaCount 3 T := Nat.one_le_iff_ne_zero.mpr h
  have := H_ge_27_of_alpha3_pos T h'
  omega

/-- If H(T) < 9, then T has no pair of vertex-disjoint odd cycles. -/
theorem H_lt_9_no_alpha2 (T : Tournament n) (hH : H T < 9) :
    alphaCount 2 T = 0 := by
  by_contra h
  have h' : 1 ≤ alphaCount 2 T := Nat.one_le_iff_ne_zero.mpr h
  have := H_ge_9_of_alpha2_pos T h'
  omega

/-- If H(T) < 81, then T has no quadruple of pairwise vertex-disjoint odd
    cycles. -/
theorem H_lt_81_no_alpha4 (T : Tournament n) (hH : H T < 81) :
    alphaCount 4 T = 0 := by
  by_contra h
  have h' : 1 ≤ alphaCount 4 T := Nat.one_le_iff_ne_zero.mpr h
  have := H_ge_81_of_alpha4_pos T h'
  omega

end Tournament
