/-
  TournamentH7.LRCHighTailIdentity -- THM-676 in Lean (klein-2026-07-09-S219).

  THE HIGH-TAIL IDENTITY: for odd depth D, the Bonferroni functional equals the
  exact live count minus a penalty supported ENTIRELY on high coverage:

      bonf D v q = liveCount v q - penalty D v q,
      penalty D v q = Σ_{p ∈ (0,q)} (bandCount v q p - 1).choose D.

  Nat.choose kills every coverage c ≤ D automatically ((c-1).choose D = 0 for
  c - 1 < D), so the penalty is exactly the mass at coverage > D -- THM-676's
  "the whole Bonferroni error is the high-coverage (discrete apex-7) mass".

  Engine: the partial alternating binomial identity
      Σ_{d ≤ D} (-1)^d C(c,d) = (-1)^D C(c-1, D)   (c ≥ 1),
  which at odd D and with ℕ-truncated subtraction becomes the UNIFORM per-p
  closed form  Σ = 1_{c=0} - C(c-1, D)  (the c = 0 case rides on 0-1 = 0 in ℕ).

  Corollaries: bonf_le_liveCount recovered with the error explicit, and THE
  DEPTH LADDER bonf D ≤ bonf (D+2) for odd 5 ≤ D ≤ 11 -- pointwise via
  C(x, D+2) ≤ C(x, D) for x ≤ 12 (kernel decide), the 13-runner domain fact
  (THM-675's escalation B5 ≤ B7 ≤ ... ≤ B13 = LM, now formal end to end
  together with death-star's bonf13_eq_liveCount).

  Kernel-pure: no native_decide, no sorry.
-/

import Mathlib
import TournamentH7.LRCDiscreteBonferroni

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- Partial alternating binomial sum, closed form (`c ≥ 1`):
`Σ_{d ≤ D} (-1)^d C(c,d) = (-1)^D C(c-1, D)`. Pascal telescopes. -/
theorem partial_alternating_choose (c D : ℕ) (hc : 1 ≤ c) :
    (∑ d ∈ range (D + 1), (-1 : ℤ) ^ d * (c.choose d : ℤ))
      = (-1) ^ D * ((c - 1).choose D : ℤ) := by
  induction D with
  | zero => simp
  | succ D ih =>
      rw [Finset.sum_range_succ, ih]
      obtain ⟨c', rfl⟩ : ∃ c', c = c' + 1 := ⟨c - 1, (Nat.succ_pred_eq_of_pos hc).symm⟩
      have hpascal : (c' + 1).choose (D + 1) = c'.choose D + c'.choose (D + 1) :=
        Nat.choose_succ_succ c' D
      have hsub : c' + 1 - 1 = c' := rfl
      rw [hsub, hpascal]
      push_cast
      ring

/-- The uniform per-multiplier closed form at odd depth: with ℕ-subtraction,
`Σ_{d ≤ D} (-1)^d C(c,d) = 1_{c=0} - C(c-1, D)` for every `c` (the `c = 0`
case rides on `0 - 1 = 0` and `C(0, D) = 0` for `D ≥ 1`). -/
theorem odd_truncation_closed_form (c D : ℕ) (hD : D % 2 = 1) :
    (∑ d ∈ range (D + 1), (-1 : ℤ) ^ d * (c.choose d : ℤ))
      = (if c = 0 then (1 : ℤ) else 0) - ((c - 1).choose D : ℤ) := by
  have hD1 : 1 ≤ D := Nat.one_le_iff_ne_zero.mpr (by rintro rfl; simp at hD)
  rcases Nat.eq_zero_or_pos c with rfl | hc
  · have hchoose : (0 - 1 : ℕ).choose D = 0 := by
      simpa using Nat.choose_eq_zero_of_lt hD1
    rw [hchoose, Finset.sum_range_succ']
    simp [Nat.choose_zero_succ]
  · have hneg : ((-1 : ℤ)) ^ D = -1 := Odd.neg_one_pow ⟨D / 2, by omega⟩
    rw [partial_alternating_choose c D hc, hneg, if_neg (Nat.pos_iff_ne_zero.mp hc)]
    ring

/-- The high-coverage penalty: mass at coverage strictly above `D`
(`Nat.choose` vanishes below, so no explicit filter is needed). -/
def penalty (D : ℕ) (v : Fin 13 → ℤ) (q : ℕ) : ℕ :=
  ∑ p ∈ Finset.Ioo 0 q, (bandCount v q p - 1).choose D

/-- **THE HIGH-TAIL IDENTITY (THM-676)**: at odd depth,
`bonf D = liveCount - penalty D` EXACTLY -- the entire Bonferroni error is
the high-coverage mass. -/
theorem bonf_eq_liveCount_sub_penalty (D : ℕ) (hD : D % 2 = 1)
    (v : Fin 13 → ℤ) (q : ℕ) :
    bonf D v q = (liveCount v q : ℤ) - (penalty D v q : ℤ) := by
  have hswap : bonf D v q
      = ∑ p ∈ Finset.Ioo 0 q,
          ∑ d ∈ range (D + 1), (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) := by
    unfold bonf momentS
    calc ∑ d ∈ range (D + 1), (-1 : ℤ) ^ d
            * ((∑ p ∈ Finset.Ioo 0 q, (bandCount v q p).choose d : ℕ) : ℤ)
        = ∑ d ∈ range (D + 1), ∑ p ∈ Finset.Ioo 0 q,
            (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) := by
          refine Finset.sum_congr rfl fun d _ => ?_
          rw [Nat.cast_sum, Finset.mul_sum]
      _ = ∑ p ∈ Finset.Ioo 0 q,
            ∑ d ∈ range (D + 1), (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) :=
          Finset.sum_comm
  rw [hswap]
  have hpt : ∀ p ∈ Finset.Ioo 0 q,
      (∑ d ∈ range (D + 1), (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ))
        = (if bandCount v q p = 0 then (1 : ℤ) else 0)
          - ((bandCount v q p - 1).choose D : ℤ) :=
    fun p _ => odd_truncation_closed_form (bandCount v q p) D hD
  rw [Finset.sum_congr rfl hpt, Finset.sum_sub_distrib]
  congr 1
  · unfold liveCount
    rw [Finset.sum_boole]
  · unfold penalty
    rw [Nat.cast_sum]

/-- The 13-runner domain fact behind the escalation ladder:
`C(x, D+2) ≤ C(x, D)` for all `x ≤ 12` and `5 ≤ D ≤ 13` (kernel decide). -/
theorem choose_ladder_dom :
    ∀ x ∈ Finset.range 13, ∀ D ∈ Finset.range 14,
      5 ≤ D → x.choose (D + 2) ≤ x.choose D := by decide

/-- **THE DEPTH LADDER (THM-675/676)**: for odd `5 ≤ D ≤ 11`,
`bonf D ≤ bonf (D+2)` -- the escalation `B5 ≤ B7 ≤ B9 ≤ B11 (≤ B13 = LM)` is
pointwise on the 13-runner domain. -/
theorem bonf_le_bonf_next (D : ℕ) (hD : D % 2 = 1) (hD5 : 5 ≤ D) (hD11 : D ≤ 11)
    (v : Fin 13 → ℤ) (q : ℕ) :
    bonf D v q ≤ bonf (D + 2) v q := by
  have hD2 : (D + 2) % 2 = 1 := by omega
  rw [bonf_eq_liveCount_sub_penalty D hD v q,
      bonf_eq_liveCount_sub_penalty (D + 2) hD2 v q]
  have hpen : penalty (D + 2) v q ≤ penalty D v q := by
    unfold penalty
    refine Finset.sum_le_sum fun p _ => ?_
    have hx : bandCount v q p - 1 ∈ Finset.range 13 := by
      have := bandCount_le_thirteen v q p
      simp only [Finset.mem_range]
      omega
    have hDmem : D ∈ Finset.range 14 := by simp; omega
    exact choose_ladder_dom _ hx D hDmem hD5
  have : (penalty (D + 2) v q : ℤ) ≤ (penalty D v q : ℤ) := by exact_mod_cast hpen
  omega

/-- The certificate reading of the identity: positivity of any odd-depth
Bonferroni functional says the live count strictly exceeds the high-coverage
penalty (and in particular is positive -- the THM-671 pipe). -/
theorem penalty_lt_liveCount_of_bonf_pos (D : ℕ) (hD : D % 2 = 1)
    (v : Fin 13 → ℤ) (q : ℕ) (h : 0 < bonf D v q) :
    (penalty D v q : ℤ) < (liveCount v q : ℤ) := by
  have := bonf_eq_liveCount_sub_penalty D hD v q
  omega

end LRC14Concrete
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.LRC14Concrete.partial_alternating_choose
#print axioms LonelyRunner.LRC14Concrete.odd_truncation_closed_form
#print axioms LonelyRunner.LRC14Concrete.bonf_eq_liveCount_sub_penalty
#print axioms LonelyRunner.LRC14Concrete.choose_ladder_dom
#print axioms LonelyRunner.LRC14Concrete.bonf_le_bonf_next
#print axioms LonelyRunner.LRC14Concrete.penalty_lt_liveCount_of_bonf_pos
