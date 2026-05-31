/-
  TournamentH7.SmallHEnumerations — α-tuple enumerations for small H

  For each small odd H value, we enumerate the (α_1, α_2, α_3, α_4)
  candidates consistent with OCF + binomial bounds.  This complements
  `Forbidden.lean` (forbidden values) by showing the structure of
  REALIZABLE small H.

  ─── Verified via Python brute-force (oracle-S15) ─────────────────────
    H = 9:  candidates = {(4,0,0,0), (2,1,0,0)}                 [2 cands]
    H = 11: {(5,0,0,0), (3,1,0,0)}                              [2 cands]
    H = 13: {(6,0,0,0), (4,1,0,0)}                              [2 cands]
    H = 15: {(7,0,0,0), (5,1,0,0), (3,2,0,0)}                   [3 cands]
    H = 17: {(8,0,0,0), (6,1,0,0), (4,2,0,0)}                   [3 cands]
    H = 19: {(9,0,0,0), (7,1,0,0), (5,2,0,0), (3,3,0,0)}        [4 cands]
    H = 23: {(11,0,0,0), (9,1,0,0), (7,2,0,0), (5,3,0,0)}       [4 cands]
    H = 25: {(12,0,0,0), (10,1,0,0), (8,2,0,0), (6,3,0,0), (4,4,0,0)}  [5 cands]

  For H ≤ 25, ALL candidates have α_3 = α_4 = 0 (the phase transition
  to α_3 ≥ 1 starts at H = 27 = 3³ — proven generally by N_min(k)=3^k).
-/

import TournamentH7.ForbiddenHCounting

namespace Tournament

variable {n : ℕ}

/-! ### H = 9 -/

/-- For H = 9, the α-tuple is (4, 0, 0, 0) or (2, 1, 0, 0). -/
theorem alpha_solution_H9 (T : Tournament n) (hH : H T = 9) :
    (alphaCount 1 T = 4 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 2 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) := by
  have hocf := ocf T
  have h3_zero : alphaCount 3 T = 0 := H_lt_27_no_alpha3 T (by omega)
  have h4_zero : alphaCount 4 T = 0 := H_lt_81_no_alpha4 T (by omega)
  -- After substituting: 1 + 2α₁ + 4α₂ = 9 ⟹ α₁ + 2α₂ = 4.
  have hRed : alphaCount 1 T + 2 * alphaCount 2 T = 4 := by omega
  by_cases h2_zero : alphaCount 2 T = 0
  · left; exact ⟨by omega, h2_zero, h3_zero, h4_zero⟩
  · -- α₂ ≥ 1; descent: α₁ ≥ 2; binomial: α₂ ≤ C(α₁, 2).
    have h2_pos : 1 ≤ alphaCount 2 T := Nat.one_le_iff_ne_zero.mpr h2_zero
    have h1_ge : Nat.choose 2 1 ≤ alphaCount 1 T :=
      alpha_descent T (by omega) h2_pos
    rw [show Nat.choose 2 1 = 2 by decide] at h1_ge
    -- α₂ ≤ C(α₁, 2). For α₁ = 2: C(2,2) = 1. For α₁ ≥ 3: more allowed.
    -- α₁ + 2α₂ = 4 with α₁ ≥ 2: (2, 1) or (4, 0). But we're in α₂ ≥ 1 branch.
    -- So α₁ = 2, α₂ = 1.
    right
    refine ⟨?_, ?_, h3_zero, h4_zero⟩ <;> omega

/-! ### H = 11 -/

/-- For H = 11, the α-tuple is (5, 0, 0, 0) or (3, 1, 0, 0). -/
theorem alpha_solution_H11 (T : Tournament n) (hH : H T = 11) :
    (alphaCount 1 T = 5 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 3 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) := by
  have hocf := ocf T
  have h3_zero : alphaCount 3 T = 0 := H_lt_27_no_alpha3 T (by omega)
  have h4_zero : alphaCount 4 T = 0 := H_lt_81_no_alpha4 T (by omega)
  have hRed : alphaCount 1 T + 2 * alphaCount 2 T = 5 := by omega
  by_cases h2_zero : alphaCount 2 T = 0
  · left; exact ⟨by omega, h2_zero, h3_zero, h4_zero⟩
  · have h2_pos : 1 ≤ alphaCount 2 T := Nat.one_le_iff_ne_zero.mpr h2_zero
    have h1_ge : Nat.choose 2 1 ≤ alphaCount 1 T :=
      alpha_descent T (by omega) h2_pos
    rw [show Nat.choose 2 1 = 2 by decide] at h1_ge
    right
    refine ⟨?_, ?_, h3_zero, h4_zero⟩ <;> omega

/-! ### H = 13 -/

/-- For H = 13, the α-tuple is (6, 0, 0, 0) or (4, 1, 0, 0). -/
theorem alpha_solution_H13 (T : Tournament n) (hH : H T = 13) :
    (alphaCount 1 T = 6 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 4 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) := by
  have hocf := ocf T
  have h3_zero : alphaCount 3 T = 0 := H_lt_27_no_alpha3 T (by omega)
  have h4_zero : alphaCount 4 T = 0 := H_lt_81_no_alpha4 T (by omega)
  have hRed : alphaCount 1 T + 2 * alphaCount 2 T = 6 := by omega
  by_cases h2_zero : alphaCount 2 T = 0
  · left; exact ⟨by omega, h2_zero, h3_zero, h4_zero⟩
  · have h2_pos : 1 ≤ alphaCount 2 T := Nat.one_le_iff_ne_zero.mpr h2_zero
    have h1_ge : Nat.choose 2 1 ≤ alphaCount 1 T :=
      alpha_descent T (by omega) h2_pos
    rw [show Nat.choose 2 1 = 2 by decide] at h1_ge
    -- α₁ + 2α₂ = 6 with α₁ ≥ 2, α₂ ≥ 1: (2,2), (4,1), (6,0 excluded).
    -- Need also α₂ ≤ C(α₁, 2). For α₁ = 2: C=1, but α₂=2 would violate.
    have h_binom := alpha_binomial_bound T 2
    rcases Nat.lt_or_ge (alphaCount 1 T) 3 with h_lt3 | h_ge3
    · -- α₁ = 2 (since α₁ ≥ 2 from descent)
      have h_eq2 : alphaCount 1 T = 2 := by omega
      rw [h_eq2] at h_binom
      rw [show Nat.choose 2 2 = 1 by decide] at h_binom
      -- α₂ ≤ 1, and α₂ = (6 - 2)/2 = 2. Contradiction.
      exfalso; omega
    · -- α₁ ≥ 3
      right
      refine ⟨?_, ?_, h3_zero, h4_zero⟩ <;> omega

/-! ### H = 15 -/

/-- For H = 15, the α-tuple is (7, 0, 0, 0) or (5, 1, 0, 0) or (3, 2, 0, 0). -/
theorem alpha_solution_H15 (T : Tournament n) (hH : H T = 15) :
    (alphaCount 1 T = 7 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 5 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 3 ∧ alphaCount 2 T = 2
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) := by
  have hocf := ocf T
  have h3_zero : alphaCount 3 T = 0 := H_lt_27_no_alpha3 T (by omega)
  have h4_zero : alphaCount 4 T = 0 := H_lt_81_no_alpha4 T (by omega)
  have hRed : alphaCount 1 T + 2 * alphaCount 2 T = 7 := by omega
  by_cases h2_zero : alphaCount 2 T = 0
  · left; exact ⟨by omega, h2_zero, h3_zero, h4_zero⟩
  have h2_pos : 1 ≤ alphaCount 2 T := Nat.one_le_iff_ne_zero.mpr h2_zero
  have h1_ge : Nat.choose 2 1 ≤ alphaCount 1 T :=
    alpha_descent T (by omega) h2_pos
  rw [show Nat.choose 2 1 = 2 by decide] at h1_ge
  have h_binom := alpha_binomial_bound T 2
  -- α₁ + 2α₂ = 7. Options: (5,1), (3,2), (1,3 ruled out by α₁ ≥ 2).
  -- (3, 2) requires α₂ ≤ C(3, 2) = 3, OK.
  -- (1, 3) ruled out by α₁ ≥ 2.
  rcases Nat.lt_or_ge (alphaCount 1 T) 4 with h_lt4 | h_ge4
  · -- α₁ ∈ {2, 3}.
    rcases Nat.lt_or_ge (alphaCount 1 T) 3 with h_lt3 | h_ge3
    · -- α₁ = 2.
      have h_eq2 : alphaCount 1 T = 2 := by omega
      -- 2 + 2α₂ = 7 ⟹ α₂ = 5/2, not integer. Contradiction.
      exfalso; omega
    · -- α₁ = 3.
      right; right
      refine ⟨by omega, ?_, h3_zero, h4_zero⟩; omega
  · -- α₁ ≥ 4.
    rcases Nat.lt_or_ge (alphaCount 1 T) 6 with h_lt6 | h_ge6
    · -- α₁ ∈ {4, 5}.
      rcases Nat.lt_or_ge (alphaCount 1 T) 5 with h_lt5 | h_ge5
      · -- α₁ = 4. 4 + 2α₂ = 7 ⟹ α₂ = 3/2, not integer.
        exfalso; omega
      · -- α₁ = 5.
        right; left
        refine ⟨by omega, ?_, h3_zero, h4_zero⟩; omega
    · -- α₁ ≥ 6. 6 + 2α₂ = 7 ⟹ α₂ negative or not integer. Contradiction.
      exfalso; omega

/-! ### H = 17 -/

/-- For H = 17, the α-tuple is (8, 0, 0, 0) or (6, 1, 0, 0) or (4, 2, 0, 0). -/
theorem alpha_solution_H17 (T : Tournament n) (hH : H T = 17) :
    (alphaCount 1 T = 8 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 6 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 4 ∧ alphaCount 2 T = 2
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) := by
  have hocf := ocf T
  have h3_zero : alphaCount 3 T = 0 := H_lt_27_no_alpha3 T (by omega)
  have h4_zero : alphaCount 4 T = 0 := H_lt_81_no_alpha4 T (by omega)
  have hRed : alphaCount 1 T + 2 * alphaCount 2 T = 8 := by omega
  by_cases h2_zero : alphaCount 2 T = 0
  · left; exact ⟨by omega, h2_zero, h3_zero, h4_zero⟩
  have h2_pos : 1 ≤ alphaCount 2 T := Nat.one_le_iff_ne_zero.mpr h2_zero
  have h1_ge : Nat.choose 2 1 ≤ alphaCount 1 T :=
    alpha_descent T (by omega) h2_pos
  rw [show Nat.choose 2 1 = 2 by decide] at h1_ge
  have h_binom := alpha_binomial_bound T 2
  -- α₁ + 2α₂ = 8 with α₁ ≥ 2, α₂ ≥ 1.
  -- Options: (2, 3), (4, 2), (6, 1).
  -- (2, 3): α₂ ≤ C(2, 2) = 1. ✗
  rcases Nat.lt_or_ge (alphaCount 1 T) 4 with h_lt4 | h_ge4
  · -- α₁ = 2 or 3.
    rcases Nat.lt_or_ge (alphaCount 1 T) 3 with h_lt3 | h_ge3
    · -- α₁ = 2. α₂ = 3. But α₂ ≤ C(2,2) = 1. Contradiction.
      have h_eq2 : alphaCount 1 T = 2 := by omega
      rw [h_eq2] at h_binom
      rw [show Nat.choose 2 2 = 1 by decide] at h_binom
      exfalso; omega
    · -- α₁ = 3. 3 + 2α₂ = 8 ⟹ α₂ = 5/2, not integer.
      exfalso; omega
  · -- α₁ ≥ 4.
    rcases Nat.lt_or_ge (alphaCount 1 T) 6 with h_lt6 | h_ge6
    · -- α₁ = 4 or 5.
      rcases Nat.lt_or_ge (alphaCount 1 T) 5 with h_lt5 | h_ge5
      · -- α₁ = 4. α₂ = 2.
        right; right
        refine ⟨by omega, ?_, h3_zero, h4_zero⟩; omega
      · -- α₁ = 5. 5 + 2α₂ = 8 ⟹ α₂ = 3/2, not integer.
        exfalso; omega
    · -- α₁ ≥ 6.
      rcases Nat.lt_or_ge (alphaCount 1 T) 7 with h_lt7 | h_ge7
      · -- α₁ = 6. α₂ = 1.
        right; left
        refine ⟨by omega, ?_, h3_zero, h4_zero⟩; omega
      · -- α₁ ≥ 7. 7 + 2α₂ ≥ 7 + 2 = 9 > 8.
        exfalso; omega

/-! ### H = 19 -/

/-- For H = 19, the α-tuple is one of (9,0), (7,1), (5,2), (3,3). -/
theorem alpha_solution_H19 (T : Tournament n) (hH : H T = 19) :
    (alphaCount 1 T = 9 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 7 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 5 ∧ alphaCount 2 T = 2
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 3 ∧ alphaCount 2 T = 3
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) := by
  have hocf := ocf T
  have h3_zero : alphaCount 3 T = 0 := H_lt_27_no_alpha3 T (by omega)
  have h4_zero : alphaCount 4 T = 0 := H_lt_81_no_alpha4 T (by omega)
  have hRed : alphaCount 1 T + 2 * alphaCount 2 T = 9 := by omega
  by_cases h2_zero : alphaCount 2 T = 0
  · left; exact ⟨by omega, h2_zero, h3_zero, h4_zero⟩
  have h2_pos : 1 ≤ alphaCount 2 T := Nat.one_le_iff_ne_zero.mpr h2_zero
  have h1_ge : Nat.choose 2 1 ≤ alphaCount 1 T :=
    alpha_descent T (by omega) h2_pos
  rw [show Nat.choose 2 1 = 2 by decide] at h1_ge
  have h_binom := alpha_binomial_bound T 2
  -- α₁ ≥ 2, α₂ ≥ 1, α₁ + 2α₂ = 9 (odd).
  -- So α₁ is odd. Options: α₁ ∈ {3, 5, 7, 9}.
  -- α₁ = 3, α₂ = 3: C(3, 2) = 3, OK.
  -- α₁ = 9, α₂ = 0: but we're in α₂ ≥ 1, excluded.
  -- So (3, 3), (5, 2), (7, 1).
  rcases Nat.lt_or_ge (alphaCount 1 T) 5 with h_lt5 | h_ge5
  · -- α₁ < 5 with α₁ ≥ 2 and odd (so α₁ = 3).
    rcases Nat.lt_or_ge (alphaCount 1 T) 3 with h_lt3 | h_ge3
    · -- α₁ = 2. 2 + 2α₂ = 9 ⟹ α₂ = 7/2.
      exfalso; omega
    · rcases Nat.lt_or_ge (alphaCount 1 T) 4 with h_lt4 | h_ge4
      · -- α₁ = 3.
        right; right; right
        refine ⟨by omega, ?_, h3_zero, h4_zero⟩; omega
      · -- α₁ = 4. 4 + 2α₂ = 9 ⟹ α₂ = 5/2.
        exfalso; omega
  · -- α₁ ≥ 5.
    rcases Nat.lt_or_ge (alphaCount 1 T) 7 with h_lt7 | h_ge7
    · -- α₁ ∈ {5, 6}.
      rcases Nat.lt_or_ge (alphaCount 1 T) 6 with h_lt6 | h_ge6
      · -- α₁ = 5.
        right; right; left
        refine ⟨by omega, ?_, h3_zero, h4_zero⟩; omega
      · -- α₁ = 6. 6 + 2α₂ = 9 ⟹ α₂ = 3/2.
        exfalso; omega
    · rcases Nat.lt_or_ge (alphaCount 1 T) 8 with h_lt8 | h_ge8
      · -- α₁ = 7.
        right; left
        refine ⟨by omega, ?_, h3_zero, h4_zero⟩; omega
      · -- α₁ ≥ 8. 8 + 2α₂ ≥ 8 + 2 = 10 > 9.
        exfalso; omega

/-! ### H = 23 -/

/-- For H = 23, the α-tuple is one of (11,0), (9,1), (7,2), (5,3). -/
theorem alpha_solution_H23 (T : Tournament n) (hH : H T = 23) :
    (alphaCount 1 T = 11 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 9 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 7 ∧ alphaCount 2 T = 2
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 5 ∧ alphaCount 2 T = 3
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) := by
  have hocf := ocf T
  have h3_zero : alphaCount 3 T = 0 := H_lt_27_no_alpha3 T (by omega)
  have h4_zero : alphaCount 4 T = 0 := H_lt_81_no_alpha4 T (by omega)
  have hRed : alphaCount 1 T + 2 * alphaCount 2 T = 11 := by omega
  by_cases h2_zero : alphaCount 2 T = 0
  · left; exact ⟨by omega, h2_zero, h3_zero, h4_zero⟩
  have h2_pos : 1 ≤ alphaCount 2 T := Nat.one_le_iff_ne_zero.mpr h2_zero
  have h1_ge : Nat.choose 2 1 ≤ alphaCount 1 T :=
    alpha_descent T (by omega) h2_pos
  rw [show Nat.choose 2 1 = 2 by decide] at h1_ge
  have h_binom := alpha_binomial_bound T 2
  -- α₁ + 2α₂ = 11. α₁ is odd.
  -- Options: (3, 4), (5, 3), (7, 2), (9, 1), (11, 0 excluded).
  -- (3, 4): C(3,2) = 3, so α₂ ≤ 3. But α₂ = 4 > 3. ✗
  rcases Nat.lt_or_ge (alphaCount 1 T) 5 with h_lt5 | h_ge5
  · rcases Nat.lt_or_ge (alphaCount 1 T) 4 with h_lt4 | h_ge4
    · rcases Nat.lt_or_ge (alphaCount 1 T) 3 with h_lt3 | h_ge3
      · -- α₁ = 2. 2 + 2α₂ = 11 ⟹ α₂ = 9/2.
        exfalso; omega
      · -- α₁ = 3. α₂ = 4. C(3,2) = 3, violation.
        have h_eq3 : alphaCount 1 T = 3 := by omega
        rw [h_eq3] at h_binom
        rw [show Nat.choose 3 2 = 3 by decide] at h_binom
        exfalso; omega
    · -- α₁ = 4. 4 + 2α₂ = 11 ⟹ α₂ = 7/2.
      exfalso; omega
  · rcases Nat.lt_or_ge (alphaCount 1 T) 7 with h_lt7 | h_ge7
    · rcases Nat.lt_or_ge (alphaCount 1 T) 6 with h_lt6 | h_ge6
      · -- α₁ = 5, α₂ = 3.
        right; right; right
        refine ⟨by omega, ?_, h3_zero, h4_zero⟩; omega
      · -- α₁ = 6. 6 + 2α₂ = 11 ⟹ α₂ = 5/2.
        exfalso; omega
    · rcases Nat.lt_or_ge (alphaCount 1 T) 9 with h_lt9 | h_ge9
      · rcases Nat.lt_or_ge (alphaCount 1 T) 8 with h_lt8 | h_ge8
        · -- α₁ = 7, α₂ = 2.
          right; right; left
          refine ⟨by omega, ?_, h3_zero, h4_zero⟩; omega
        · -- α₁ = 8. 8 + 2α₂ = 11 ⟹ α₂ = 3/2.
          exfalso; omega
      · rcases Nat.lt_or_ge (alphaCount 1 T) 10 with h_lt10 | h_ge10
        · -- α₁ = 9, α₂ = 1.
          right; left
          refine ⟨by omega, ?_, h3_zero, h4_zero⟩; omega
        · -- α₁ ≥ 10. α₁ ≥ 10 + 2*1 = 12 > 11.
          exfalso; omega

/-! ### Summary: realizable small H values -/

/-- For H ∈ {1, 3, 5, 9, 11, 13}, ALL α-tuples have α_3 = α_4 = 0
    (the "α_3-zero phase"). -/
theorem small_H_alpha34_zero (T : Tournament n) (hH : H T ≤ 26) :
    alphaCount 3 T = 0 ∧ alphaCount 4 T = 0 := by
  refine ⟨H_lt_27_no_alpha3 T (by omega), H_lt_81_no_alpha4 T (by omega)⟩

end Tournament
