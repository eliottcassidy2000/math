/-
  TournamentH7.LRCINVcovCounterexample — an exact counterexample to the literal
  covering inverse target `INVcov`.

  The doubled tight arithmetic progression

      v_i = 2(i+1),  i : Fin 13,

  is positive and covers every modulus 2..14.  Dilation invariance reduces any
  hypothetical 1/13-lonely time for this family to one for the undilated AP,
  contradicting the sharp Dirichlet ceiling `minReach_AP_le <= 1/14`.  Yet no
  entry is 13-fold dominant: speed 4 defeats every proposed dominant entry
  except itself, and speed 6 defeats that remaining proposal.

  Thus the literal `INVcov` predicate in `LRCMSplit` is false.  Any viable
  noncircular inverse interface must exclude tight AP dilations (for example by
  a primitive/gcd normalization plus a separate equality branch), or weaken
  its conclusion.
-/
import Mathlib
import TournamentH7.LRCMSplit
import TournamentH7.LRCAPTight
import TournamentH7.LRCDilation

namespace LonelyRunner
namespace INVcovCounterexample

/-- The doubled 13-term arithmetic progression `2,4,...,26`. -/
def doubledAP (i : Fin 13) : ℤ := 2 * ((i.val : ℤ) + 1)

/-- Every speed in the doubled AP is positive. -/
theorem doubledAP_pos (i : Fin 13) : 0 < doubledAP i := by
  unfold doubledAP
  omega

/-- The doubled AP covers every modulus in the divisor-sieve window `2..14`.
For `n <= 13`, the entry with index `n-1` equals `2n`; for `n=14`, use the
entry `14` at index `6`. -/
theorem doubledAP_covering : Covering doubledAP := by
  intro n hn2 hn14
  by_cases hn13 : n ≤ 13
  · let i : Fin 13 := ⟨n - 1, by omega⟩
    refine ⟨i, 2, ?_⟩
    simp only [doubledAP, i]
    omega
  · have hn : n = 14 := by omega
    subst n
    refine ⟨⟨6, by norm_num⟩, 1, ?_⟩
    norm_num [doubledAP]

/-- A `1/13`-lonely time forces the corresponding min-reach to be at least
`1/13`, for the AP family used below. -/
private theorem one_thirteenth_le_minReach_of_lonely
    (t : ℝ)
    (h : Lonely 13 (fun i : Fin 13 => (i.val : ℤ) + 1) t) :
    (1 : ℝ) / 13 ≤ LRC14Concrete.minReach (fun i : Fin 13 => (i.val : ℤ) + 1) t := by
  unfold LRC14Concrete.minReach
  refine le_ciInf (fun i => ?_)
  exact LRC14Concrete.le_nearInt_of_forall_int _ _ (h i)

/-- The doubled AP has no `1/13`-lonely time.  By `lonely_smul`, such a time
would give a `1/13`-lonely time for the undilated AP at twice that time; its
min-reach would then be at least `1/13`, contradicting `minReach_AP_le`. -/
theorem doubledAP_no_lonely13 : ¬ ∃ t : ℝ, Lonely 13 doubledAP t := by
  rintro ⟨t, ht⟩
  have hscaled : Lonely 13 (fun i : Fin 13 => (i.val : ℤ) + 1) ((2 : ℝ) * t) := by
    apply (lonely_smul 13 (fun i : Fin 13 => (i.val : ℤ) + 1) 2 t).mp
    simpa [doubledAP] using ht
  have hlo := one_thirteenth_le_minReach_of_lonely ((2 : ℝ) * t) hscaled
  have hhi := LRC14Concrete.minReach_AP_le ((2 : ℝ) * t)
  norm_num at hlo hhi
  linarith

/-- No entry of the doubled AP dominates every other entry by a factor of 13.
Speed `4` is already too large for any proposed dominant entry (all are at
most `26`), unless that proposed entry is speed `4` itself; speed `6` handles
that one exceptional index. -/
theorem doubledAP_no_thirteen_fold_dominance :
    ¬ ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * doubledAP i ≤ doubledAP vstar := by
  rintro ⟨vstar, hdom⟩
  let i4 : Fin 13 := ⟨1, by norm_num⟩
  let i6 : Fin 13 := ⟨2, by norm_num⟩
  by_cases hs : i4 = vstar
  · have hne : i6 ≠ vstar := by
      intro h
      have : i6 = i4 := h.trans hs.symm
      simp [i4, i6] at this
    have hbad := hdom i6 hne
    have hvstar : (vstar.val : ℤ) ≤ 12 := by omega
    simp only [doubledAP, i6] at hbad
    norm_num at hbad
    omega
  · have hbad := hdom i4 hs
    have hvstar : (vstar.val : ℤ) ≤ 12 := by omega
    simp only [doubledAP, i4] at hbad
    norm_num at hbad
    omega

/-- **Exact refutation of the literal covering inverse target.** -/
theorem not_INVcov : ¬ INVcov := by
  intro hinv
  obtain ⟨vstar, hdom⟩ :=
    hinv doubledAP doubledAP_pos doubledAP_covering doubledAP_no_lonely13
  exact doubledAP_no_thirteen_fold_dominance ⟨vstar, hdom⟩

end INVcovCounterexample
end LonelyRunner

#print axioms LonelyRunner.INVcovCounterexample.doubledAP_pos
#print axioms LonelyRunner.INVcovCounterexample.doubledAP_covering
#print axioms LonelyRunner.INVcovCounterexample.doubledAP_no_lonely13
#print axioms LonelyRunner.INVcovCounterexample.doubledAP_no_thirteen_fold_dominance
#print axioms LonelyRunner.INVcovCounterexample.not_INVcov
