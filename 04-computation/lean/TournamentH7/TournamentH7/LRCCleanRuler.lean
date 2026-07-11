/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: kind-pasteur (LRC multi-agent project, 2026-07-11-S127 cont.28)
-/
import TournamentH7.LRCDiscreteBonferroni

/-!
# THM-707 — the clean-ruler lemma: shallow coverage collapses `B5` to `liveCount`

The single remaining LRC(14) obligation `hB5` asks each residual covering family for a modulus `q` with
`0 < B5 v q`, where `B5 = S₀ − S₁ + S₂ − S₃ + S₄ − S₅` is the depth-5 alternating Bonferroni functional of
the coverage histogram (`LRCDiscreteBonferroni`).  THM-671 gives `B5 ≤ liveCount`.  This file proves the
**exact** companion on *clean* rulers, via the elementary identity
`∑_{d<6} (−1)^d C(n,d) = 0` for `1 ≤ n ≤ 5` (and `= 1` at `n = 0`):

> **`b5_pos_of_clean`**: if every multiplier covers at most `5` runners (`bandCount ≤ 5`) and some multiplier
> is live, then `B5 v q = liveCount v q > 0`.

So on a clean ruler the opaque signed Bonferroni is exactly the count of loneliness witnesses — no
deep-coverage penalty.  This reduces `hB5` to the transparent seven-sector statement "every residual family
has a ruler with a live multiplier and no multiplier covering ≥ 6 runners".
-/

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- **The depth-5 truncation on shallow overlaps.**  For `n ≤ 5`, the alternating binomial partial sum
`∑_{d<6} (−1)^d C(n,d)` is the live indicator: `1` if `n = 0`, else `0`.  (The general identity is
`∑_{d≤k}(−1)^d C(n,d) = (−1)^k C(n−1,k)`, which vanishes for `1 ≤ n ≤ k`.) -/
lemma trunc5_of_le_five (n : ℕ) (h : n ≤ 5) :
    ∑ d ∈ range 6, (-1 : ℤ) ^ d * ((n.choose d : ℤ)) = (if n = 0 then 1 else 0) := by
  interval_cases n <;> decide

/-- **The clean-ruler lemma (THM-707).**  If the ruler `q` has shallow coverage (`bandCount ≤ 5` at every
multiplier) and at least one live multiplier, then `B5 v q = liveCount v q`, hence `0 < B5 v q`. -/
theorem b5_pos_of_clean (v : Fin 13 → ℤ) (q : ℕ)
    (hclean : ∀ p ∈ Finset.Ioo 0 q, bandCount v q p ≤ 5)
    (hlive : 0 < liveCount v q) :
    0 < B5 v q := by
  -- swap the depth/multiplier sums (verbatim from `B5_le_liveCount`)
  have hswap : B5 v q
      = ∑ p ∈ Finset.Ioo 0 q,
          ∑ d ∈ range 6, (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) := by
    unfold B5 momentS
    calc ∑ d ∈ range 6, (-1 : ℤ) ^ d
            * ((∑ p ∈ Finset.Ioo 0 q, (bandCount v q p).choose d : ℕ) : ℤ)
        = ∑ d ∈ range 6, ∑ p ∈ Finset.Ioo 0 q,
            (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) := by
          refine Finset.sum_congr rfl fun d _ => ?_
          rw [Nat.cast_sum, Finset.mul_sum]
      _ = ∑ p ∈ Finset.Ioo 0 q,
            ∑ d ∈ range 6, (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ) :=
          Finset.sum_comm
  -- pointwise: shallow coverage ⟹ each multiplier contributes its live indicator
  have hpt : ∀ p ∈ Finset.Ioo 0 q,
      (∑ d ∈ range 6, (-1 : ℤ) ^ d * ((bandCount v q p).choose d : ℤ))
        = (if bandCount v q p = 0 then (1 : ℤ) else 0) := by
    intro p hp
    have h5 : bandCount v q p ≤ 5 := hclean p hp
    revert h5
    generalize bandCount v q p = n
    intro h5
    exact trunc5_of_le_five n h5
  -- so `B5 = ∑ live-indicator = liveCount`
  have heq : B5 v q = (liveCount v q : ℤ) := by
    rw [hswap, Finset.sum_congr rfl hpt]
    unfold liveCount
    rw [Finset.sum_boole]
  rw [heq]
  exact_mod_cast hlive

/-- **A clean ruler for the family `v`**: a positive modulus with shallow coverage everywhere and at least
one live multiplier.  This is the transparent seven-sector condition the architecture reduces `hB5` to. -/
def CleanRuler (v : Fin 13 → ℤ) : Prop :=
  ∃ q : ℕ, 0 < q ∧ (∀ p ∈ Finset.Ioo 0 q, bandCount v q p ≤ 5) ∧ 0 < liveCount v q

/-- **The reduction (THM-707).**  A clean ruler supplies exactly the per-family `hB5` witness
`∃ q, 0 < q ∧ 0 < B5 v q`.  Hence a clean-ruler supply over all residual families discharges `hB5`, and with
it the last analytic obligation of the LRC(14) grand assembly. -/
theorem exists_B5_pos_of_cleanRuler (v : Fin 13 → ℤ) (h : CleanRuler v) :
    ∃ q : ℕ, 0 < q ∧ 0 < B5 v q := by
  obtain ⟨q, hq, hclean, hlive⟩ := h
  exact ⟨q, hq, b5_pos_of_clean v q hclean hlive⟩

end LRC14Concrete
end LonelyRunner

#print axioms LonelyRunner.LRC14Concrete.b5_pos_of_clean
#print axioms LonelyRunner.LRC14Concrete.trunc5_of_le_five
#print axioms LonelyRunner.LRC14Concrete.exists_B5_pos_of_cleanRuler
