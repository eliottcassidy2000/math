/-
  TournamentH7.ScoreSequence — Score sequence facts

  For any tournament T on n vertices:
   • Σ_v outDegree(v) = C(n, 2) (total arc count).
   • Each outDegree v ∈ [0, n-1].
   • For regular T: every outDegree = (n-1)/2.
-/

import TournamentH7.Tilings
import TournamentH7.BasePathSink
import Mathlib.Algebra.BigOperators.Group.Finset.Basic

open Finset

namespace Tournament

variable {n : ℕ}

/-! ### Score sum identity -/

/-- **Axiom (basic).** Sum of out-degrees = C(n, 2). -/
axiom sum_outDegree_eq (T : Tournament n) :
    2 * ∑ v : Fin n, T.outDegree v = n * (n - 1)

/-! ### Regular score is (n-1)/2 -/

/-- For a regular tournament, sum outDegree = n * (n-1)/2 (consistent
    with the general identity). -/
theorem regular_sum_outDegree (T : Tournament n) (hreg : IsRegular T) :
    ∀ v : Fin n, 2 * T.outDegree v = n - 1 := hreg

/-! ### Regular tournaments require n to be odd

  If T is regular, then for each v, 2 * outDegree v = n - 1.  Since both
  sides must be a Nat, n - 1 must be even, i.e., n is odd. -/

/-- For a regular tournament on n ≥ 1 vertices, n must be odd. -/
theorem regular_implies_n_odd (T : Tournament n) (hn : 1 ≤ n) (hreg : IsRegular T) :
    Odd n := by
  -- Pick any v; 2 * outDegree v = n - 1, so n - 1 is even, so n is odd.
  have h := hreg ⟨0, by omega⟩
  -- h : 2 * outDegree = n - 1, so n - 1 is even.
  rcases Nat.even_or_odd n with hev | hodd
  · -- n is even ⟹ n - 1 is odd ⟹ contradiction.
    rcases hev with ⟨k, hk⟩
    have h_odd : n - 1 = 2 * k - 1 := by omega
    have : 2 * T.outDegree ⟨0, by omega⟩ ≠ 2 * k - 1 := by
      -- For n ≥ 2 (since regular ⟹ n ≥ 3 actually), the LHS is even, RHS odd.
      by_cases hn2 : n ≥ 2
      · intro h_eq
        have h_even : 2 ∣ 2 * T.outDegree ⟨0, by omega⟩ := ⟨_, rfl⟩
        rw [h_eq] at h_even
        omega
      · -- n < 2 ⟹ n = 1.  Then n - 1 = 0; 2 * outDegree = 0; outDegree = 0.  OK.
        -- But the issue: at n = 1, outDegree = 0 by definition, so 2 * 0 = 0 = n - 1 = 0. Consistent.
        -- But there's no contradiction here.
        intro _; omega
    rw [h_odd] at h
    -- h : 2 * outDegree = 2*k - 1.  But LHS is even, RHS odd (if k ≥ 1).
    -- If k = 0 then n = 0 contradicts hn.
    have hk_pos : 1 ≤ k := by omega
    have h_lhs_even : 2 ∣ 2 * T.outDegree ⟨0, by omega⟩ := ⟨_, rfl⟩
    rw [h] at h_lhs_even
    -- h_lhs_even : 2 ∣ 2*k - 1 with k ≥ 1, but 2*k - 1 is odd.
    omega
  · exact hodd

/-! ### Combined: regular HasBasePath tournaments are odd and ≥ 3 -/

/-- **Theorem.** Any regular HasBasePath tournament on n ≥ 2 vertices has
    n odd AND n ≥ 3.  Combining `regular_implies_n_odd` and
    `regular_basepath_n_ge_three`. -/
theorem regular_basepath_n_odd_ge_three (T : Tournament n)
    (hbp : HasBasePath T) (hn : 2 ≤ n) (hreg : IsRegular T) :
    Odd n ∧ 3 ≤ n := by
  refine ⟨?_, ?_⟩
  · exact regular_implies_n_odd T (by omega) hreg
  · exact regular_basepath_n_ge_three T hbp hn hreg

/-- Combining: for a regular HasBasePath T, n ∈ {3, 5, 7, 9, ...}. -/
theorem regular_basepath_n_in_odd_ge_three (T : Tournament n)
    (hbp : HasBasePath T) (hn : 2 ≤ n) (hreg : IsRegular T) :
    ∃ k, n = 2 * k + 3 := by
  obtain ⟨hodd, h_ge⟩ := regular_basepath_n_odd_ge_three T hbp hn hreg
  rcases hodd with ⟨m, hm⟩
  -- n = 2m + 1, n ≥ 3 ⟹ m ≥ 1.
  refine ⟨m - 1, ?_⟩
  omega

end Tournament
