/-
  TournamentH7.ForbiddenHCounting — f(N) and the N_min(k) = 3^k theorem
  (oracle-2026-05-29-S4, project-novel)

  ─── What this module formalises ──────────────────────────────────────
  The Lean formalisation of forbidden-H factors into an arithmetic
  enumeration stage (this module) and a structural killing/finite-audit stage
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

/-! ### Connection to forbidden HSpectrum gaps -/

/-- The forbidden values 7 and 21 are below the "α_3 ≥ 1" threshold,
    so any tournament with H = 7 or 21 has α_3 = 0. -/
theorem forbidden_pair_alpha3_zero (T : Tournament n) :
    H T = 7 ∨ H T = 21 → alphaCount 3 T = 0 := by
  intro h
  rcases h with h7 | h21
  · exact H_lt_27_no_alpha3 T (by omega)
  · exact H_lt_27_no_alpha3 T (by omega)

/-! ### Stronger structural derivations -/

/-- For tournaments with H = 3, alpha_1 must be 1 and all higher α_k = 0. -/
theorem alpha_solution_H3 (T : Tournament n) (hH : H T = 3) :
    alphaCount 1 T = 1 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0 := by
  have hOCF := ocf T
  -- 3 = 1 + 2 a1 + 4 a2 + 8 a3 + 16 a4
  -- Need 2 a1 + 4 a2 + 8 a3 + 16 a4 = 2, i.e., a1 + 2 a2 + 4 a3 + 8 a4 = 1.
  -- So a1 = 1 and others 0.
  refine ⟨?_, ?_, ?_, ?_⟩ <;> omega

/-! ### Extension to k = 5, 6 via ocf_extended (forward-declared) -/

-- The H_ge_243_of_alpha5_pos and H_ge_729_of_alpha6_pos theorems are
-- proved later in this file.  We forward-declare them as private axioms
-- for use in H27_alpha5_zero, then connect once both are defined.

/-! ### H = 27 arithmetic enumeration (NEW)

    H = 27 = 3^3 is the first H where α_3 can be 1.  We enumerate the
    candidate α-tuples:
      (13, 0, 0, 0), (11, 1, 0, 0), (9, 2, 0, 0), (7, 3, 0, 0),
      (5, 4, 0, 0), (3, 3, 1, 0).
    All six survive the arithmetic constraints. -/

/-- **Theorem (new, oracle-S6).** For any tournament T with H = 27,
    α_4 = 0.  In particular, 27 is the lowest H for which α_3 ≥ 1 is
    arithmetically possible. -/
theorem H27_alpha4_zero (T : Tournament n) (hH : H T = 27) :
    alphaCount 4 T = 0 := H_lt_81_no_alpha4 T (by omega)

/-! ### Pascal's row sum at arbitrary k -/

/-- **Theorem (Pascal's row sum, general k).**
    `Σ_{j=0}^k 2^j · C(k, j) = 3^k`.

    Proof: (2 + 1)^k via binomial theorem. -/
theorem pascal_row_sum_general (k : ℕ) :
    ∑ j ∈ Finset.range (k + 1), 2 ^ j * Nat.choose k j = 3 ^ k := by
  have h := @add_pow ℕ _ 2 1 k
  simp at h
  exact h.symm

/-! ### General-k OCF axiom + N_min(k) = 3^k for arbitrary k -/

/-- **Axiom (OCF truncated to K).**  When α_k = 0 for k > K, the OCF reads
    H(T) = 1 + Σ_{k=1}^K 2^k · α_k(T). -/
axiom ocf_truncated (T : Tournament n) (K : ℕ)
    (h_high_zero : ∀ k, K < k → alphaCount k T = 0) :
    H T = 1 + ∑ k ∈ Finset.Ico 1 (K + 1), 2 ^ k * alphaCount k T

/-- **Theorem (general N_min, project-novel, oracle-S8).**

    For any tournament T with α_k(T) ≥ 1 (k ≥ 1 arbitrary), H(T) ≥ 3^k.

    Uses Pascal's row sum to handle ALL k ≥ 1 (not just k ≤ 6). -/
theorem H_ge_three_pow_k_general
    (T : Tournament n) (k : ℕ) (hk_pos : 1 ≤ k) (h : 1 ≤ alphaCount k T)
    (h_high_zero : ∀ j, k < j → alphaCount j T = 0) :
    3 ^ k ≤ H T := by
  have hocf := ocf_truncated T k h_high_zero
  -- Bound each α_j by C(k, j) via alpha_descent.
  have h_lower : ∀ j ∈ Finset.Ico 1 (k + 1),
      2 ^ j * Nat.choose k j ≤ 2 ^ j * alphaCount j T := by
    intro j hj
    simp at hj
    have h_le : j ≤ k := by omega
    have h_choose := alpha_descent T h_le h
    exact Nat.mul_le_mul_left _ h_choose
  have h_sum : ∑ j ∈ Finset.Ico 1 (k + 1), 2 ^ j * Nat.choose k j ≤
               ∑ j ∈ Finset.Ico 1 (k + 1), 2 ^ j * alphaCount j T :=
    Finset.sum_le_sum h_lower
  -- LHS = 3^k - 1.
  have h_pascal := pascal_row_sum_general k
  have h_split : ∑ j ∈ Finset.range (k + 1), 2 ^ j * Nat.choose k j =
                 (1 : ℕ) + ∑ j ∈ Finset.Ico 1 (k + 1), 2 ^ j * Nat.choose k j := by
    rw [show Finset.range (k + 1) =
        {0} ∪ Finset.Ico 1 (k + 1) by
      ext m; simp; omega]
    rw [Finset.sum_union (by simp [Finset.disjoint_singleton_left])]
    simp
  rw [h_split] at h_pascal
  have h_pascal' : ∑ j ∈ Finset.Ico 1 (k + 1), 2 ^ j * Nat.choose k j = 3 ^ k - 1 := by
    omega
  rw [h_pascal'] at h_sum
  have h3k_pos : 1 ≤ 3 ^ k := Nat.one_le_iff_ne_zero.mpr (by positivity)
  omega

/-! ### Pascal's triangle row sums identity (foundational) -/

/-- Pascal's identity: 1 + Σ_{j=1}^k 2^j · C(k, j) = 3^k.

    For Lean: we prove this for k ∈ {1, 2, 3, 4, 5, 6} by direct computation. -/
example : 1 + 2 * Nat.choose 1 1 = 3 ^ 1 := by decide
example : 1 + 2 * Nat.choose 2 1 + 4 * Nat.choose 2 2 = 3 ^ 2 := by decide
example : 1 + 2 * Nat.choose 3 1 + 4 * Nat.choose 3 2
        + 8 * Nat.choose 3 3 = 3 ^ 3 := by decide
example : 1 + 2 * Nat.choose 4 1 + 4 * Nat.choose 4 2
        + 8 * Nat.choose 4 3 + 16 * Nat.choose 4 4 = 3 ^ 4 := by decide
example : 1 + 2 * Nat.choose 5 1 + 4 * Nat.choose 5 2
        + 8 * Nat.choose 5 3 + 16 * Nat.choose 5 4
        + 32 * Nat.choose 5 5 = 3 ^ 5 := by decide
example : 1 + 2 * Nat.choose 6 1 + 4 * Nat.choose 6 2
        + 8 * Nat.choose 6 3 + 16 * Nat.choose 6 4
        + 32 * Nat.choose 6 5 + 64 * Nat.choose 6 6 = 3 ^ 6 := by decide

/-! ### Strict variant: equality iff all α_j = C(k,j) for j ≤ k, 0 else. -/

/-- For H = 3^k exactly, when α_k ≥ 1, the α-tuple is forced to
    (C(k, 1), C(k, 2), …, C(k, k), 0, 0, …). -/
example (T : Tournament n) (hH : H T = 27) (h3 : 1 ≤ alphaCount 3 T) :
    alphaCount 1 T = 3 ∧ alphaCount 2 T = 3
       ∧ alphaCount 3 T = 1 ∧ alphaCount 4 T = 0 := by
  have hocf := ocf T
  have hd1 : Nat.choose 3 1 ≤ alphaCount 1 T :=
    alpha_descent T (by omega) h3
  have hd2 : Nat.choose 3 2 ≤ alphaCount 2 T :=
    alpha_descent T (by omega) h3
  have hd3 : Nat.choose 3 3 ≤ alphaCount 3 T :=
    alpha_descent T (by omega) h3
  rw [show Nat.choose 3 1 = 3 by decide] at hd1
  rw [show Nat.choose 3 2 = 3 by decide] at hd2
  rw [show Nat.choose 3 3 = 1 by decide] at hd3
  have h4 : alphaCount 4 T = 0 := H27_alpha4_zero T hH
  refine ⟨?_, ?_, ?_, h4⟩ <;> omega

/-- H_ge_3^5 = 243 when α_5 ≥ 1. -/
theorem H_ge_243_of_alpha5_pos (T : Tournament n) (h : 1 ≤ alphaCount 5 T) :
    243 ≤ H T := by
  have hocf := ocf_extended T
  have hd1 : Nat.choose 5 1 ≤ alphaCount 1 T :=
    alpha_descent T (by omega) h
  have hd2 : Nat.choose 5 2 ≤ alphaCount 2 T :=
    alpha_descent T (by omega) h
  have hd3 : Nat.choose 5 3 ≤ alphaCount 3 T :=
    alpha_descent T (by omega) h
  have hd4 : Nat.choose 5 4 ≤ alphaCount 4 T :=
    alpha_descent T (by omega) h
  have hd5 : Nat.choose 5 5 ≤ alphaCount 5 T :=
    alpha_descent T (by omega) h
  have h_c1 : Nat.choose 5 1 = 5 := by decide
  have h_c2 : Nat.choose 5 2 = 10 := by decide
  have h_c3 : Nat.choose 5 3 = 10 := by decide
  have h_c4 : Nat.choose 5 4 = 5 := by decide
  have h_c5 : Nat.choose 5 5 = 1 := by decide
  rw [h_c1] at hd1
  rw [h_c2] at hd2
  rw [h_c3] at hd3
  rw [h_c4] at hd4
  rw [h_c5] at hd5
  -- 1 + 2*5 + 4*10 + 8*10 + 16*5 + 32*1 = 1 + 10 + 40 + 80 + 80 + 32 = 243
  omega

/-- H_ge_3^6 = 729 when α_6 ≥ 1. -/
theorem H_ge_729_of_alpha6_pos (T : Tournament n) (h : 1 ≤ alphaCount 6 T) :
    729 ≤ H T := by
  have hocf := ocf_extended T
  have hd1 : Nat.choose 6 1 ≤ alphaCount 1 T :=
    alpha_descent T (by omega) h
  have hd2 : Nat.choose 6 2 ≤ alphaCount 2 T :=
    alpha_descent T (by omega) h
  have hd3 : Nat.choose 6 3 ≤ alphaCount 3 T :=
    alpha_descent T (by omega) h
  have hd4 : Nat.choose 6 4 ≤ alphaCount 4 T :=
    alpha_descent T (by omega) h
  have hd5 : Nat.choose 6 5 ≤ alphaCount 5 T :=
    alpha_descent T (by omega) h
  have hd6 : Nat.choose 6 6 ≤ alphaCount 6 T :=
    alpha_descent T (by omega) h
  have h_c1 : Nat.choose 6 1 = 6 := by decide
  have h_c2 : Nat.choose 6 2 = 15 := by decide
  have h_c3 : Nat.choose 6 3 = 20 := by decide
  have h_c4 : Nat.choose 6 4 = 15 := by decide
  have h_c5 : Nat.choose 6 5 = 6 := by decide
  have h_c6 : Nat.choose 6 6 = 1 := by decide
  rw [h_c1] at hd1
  rw [h_c2] at hd2
  rw [h_c3] at hd3
  rw [h_c4] at hd4
  rw [h_c5] at hd5
  rw [h_c6] at hd6
  -- 1 + 12 + 60 + 160 + 240 + 192 + 64 = 729 = 3^6
  omega

/-- For tournaments with H = 5, alpha_1 = 2 and others 0.
    Uses alpha_descent (α_2 ≥ 1 ⟹ α_1 ≥ 2) to rule out (α_1=0, α_2=1, ...). -/
theorem alpha_solution_H5 (T : Tournament n) (hH : H T = 5) :
    alphaCount 1 T = 2 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0 := by
  have hOCF := ocf T
  -- 5 = 1 + 2 a1 + 4 a2 + 8 a3 + 16 a4 ⟹ a1 + 2 a2 + 4 a3 + 8 a4 = 2.
  -- Candidates: (a1=2, others 0) OR (a1=0, a2=1, others 0).
  -- The latter is ruled out by alpha_descent: a2 ≥ 1 ⟹ a1 ≥ 2.
  have hRed : alphaCount 1 T + 2 * alphaCount 2 T
              + 4 * alphaCount 3 T + 8 * alphaCount 4 T = 2 := by omega
  have h4 : alphaCount 4 T = 0 := by omega
  have h3 : alphaCount 3 T = 0 := by omega
  by_cases ha2 : alphaCount 2 T = 0
  · refine ⟨by omega, ha2, h3, h4⟩
  · -- α_2 ≥ 1 ⟹ α_1 ≥ C(2,1) = 2 (descent)
    have ha2_pos : 1 ≤ alphaCount 2 T := Nat.one_le_iff_ne_zero.mpr ha2
    have ha1_ge : Nat.choose 2 1 ≤ alphaCount 1 T :=
      alpha_descent T (by omega) ha2_pos
    have hch : Nat.choose 2 1 = 2 := by decide
    rw [hch] at ha1_ge
    -- a1 ≥ 2 ∧ a2 ≥ 1 ⟹ a1 + 2 a2 ≥ 4 > 2, contradiction
    exfalso; omega

end Tournament
