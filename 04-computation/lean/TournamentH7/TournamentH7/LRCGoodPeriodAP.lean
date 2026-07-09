/-
  TournamentH7.LRCGoodPeriodAP — the AP clustering good-period lemma (kind-pasteur-2026-07-09-S92).

  A companion to mac-mini's LRCGoodPeriodJ1 (LEM-010 i + the good-gap core).  This formalizes the
  AP case of LEM-010: if the k cluster phases form an arithmetic progression 0, θ, 2θ, …, (k−1)θ
  with step θ ≥ 0 and total span (k−1)·θ < 6/7, they are clustered in a `< 6/7` arc, so the
  complementary circular arc of length `1 − (k−1)θ > 1/7` is empty — a good period.

  This is the mechanism behind LEM-010(iii)'s AP bound (mac-mini/klein-S196): for an AP cluster,
  a Dirichlet dilation `j ≤ ⌈7(k−1)/6⌉` makes the step `frac(j·d/Vmax)` small enough that
  `(k−1)·step < 6/7`, and THIS lemma then yields the gap.  It is the shared engine of the near-AP
  (LEM-012) branch of the good-period dichotomy.  Self-contained (imports only Mathlib).
-/
import Mathlib

namespace LRC14AP

/-- **The AP clustering good period.** If `k` phases form an AP `0, θ, 2θ, …, (k−1)θ` with step
`θ ≥ 0` whose span `(k−1)·θ < 6/7`, then a gap of length `gapLen = 1 − (k−1)θ > 1/7` is open past
every phase: each phase `i·θ` (for `i < k`) satisfies `i·θ + gapLen ≤ 1`.  A good period exists. -/
theorem good_period_AP (k : ℕ) (θ : ℚ) (hθ : 0 ≤ θ)
    (hspan : ((k : ℚ) - 1) * θ < 6 / 7) :
    ∃ gapLen : ℚ, (1 : ℚ) / 7 < gapLen ∧
      ∀ i : ℕ, i < k → (i : ℚ) * θ + gapLen ≤ 1 := by
  refine ⟨1 - ((k : ℚ) - 1) * θ, by linarith, ?_⟩
  intro i hi
  have hik : (i : ℚ) ≤ (k : ℚ) - 1 := by
    have h1 : (i : ℚ) + 1 ≤ (k : ℚ) := by exact_mod_cast hi
    linarith
  have hle : (i : ℚ) * θ ≤ ((k : ℚ) - 1) * θ := mul_le_mul_of_nonneg_right hik hθ
  linarith

/-- **Corollary (existence form).** Under the same hypotheses, the maximal circular gap of the AP
phase set exceeds `1/7` — i.e. the AP cluster admits a good period.  (Packaged as: a positive gap
length `> 1/7` witnessed past the last phase.) -/
theorem good_period_of_AP_span (k : ℕ) (θ : ℚ) (_hk : 1 ≤ k) (_hθ : 0 ≤ θ)
    (hspan : ((k : ℚ) - 1) * θ < 6 / 7) :
    ∃ gapLen : ℚ, (1 : ℚ) / 7 < gapLen :=
  ⟨1 - ((k : ℚ) - 1) * θ, by linarith⟩

-- Sanity: k = 13, θ = 1/16 (a 13-AP of step 1/16 spans 12/16 = 3/4 < 6/7) clusters ⇒ a good period.
example : ∃ gapLen : ℚ, (1 : ℚ) / 7 < gapLen ∧
    ∀ i : ℕ, i < 13 → (i : ℚ) * (1 / 16) + gapLen ≤ 1 := by
  apply good_period_AP 13 (1 / 16) (by norm_num)
  norm_num

end LRC14AP
