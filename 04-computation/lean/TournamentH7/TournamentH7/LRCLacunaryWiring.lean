/-
  TournamentH7.LRCLacunaryWiring — remove the 7/3-lacunary families from the
  chain-dichotomy dense core.

  `LRCLacunaryNest.lonely_of_pos_lacunary` closes every positive ordered family
  whose consecutive ratios are at least 7/3.  This module transports that theorem
  through absolute values and a sorting permutation, then exposes the genuinely
  smaller grand-assembly residual: a sorted `ChainDenseCore` which is *not*
  7/3-lacunary.

  Assumption challenge / Tournament Analysis: the natural vertices in this branch
  are the twelve adjacent scale gaps, not the thirteen runners or their arcs.  A
  binary tournament relation would retain only which speed is larger; it destroys
  the quantitative 7/3 threshold and the nested safe-interval widths used by the
  proof.  Accordingly this is recorded as an ordered ratio-chain certificate rather
  than forcing a lossy tournament quotient.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import TournamentH7.LRCLacunaryNest
import TournamentH7.LRCChainDichotomy

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner

/-- The twelve adjacent ratios of an ordered thirteen-speed family are all at
least `7/3`.  The cross-multiplied integer form avoids division and sign issues. -/
def SevenThreeLacunary (w : Fin 13 → ℤ) : Prop :=
  ∀ i : Fin 12, 7 * w i.castSucc ≤ 3 * w i.succ

/-- The universal lacunary theorem transported through a relabeling and arbitrary
speed signs.  This is the interface consumed by residual reductions. -/
theorem lonely_of_abs_perm_sevenThreeLacunary (v : Fin 13 → ℤ)
    (hv : ∀ i, v i ≠ 0) (σ : Equiv.Perm (Fin 13))
    (hlac : SevenThreeLacunary (fun i ↦ |v (σ i)|)) :
    ∃ t : ℝ, Lonely 14 v t := by
  let w : Fin 13 → ℤ := fun i ↦ |v (σ i)|
  have hwpos : ∀ i, 0 < w i := fun i ↦ abs_pos.mpr (hv (σ i))
  have hwchain : ∀ i : Fin 12, 7 * w i.castSucc ≤ 3 * w i.succ := by
    simpa [SevenThreeLacunary, w] using hlac
  obtain ⟨t, ht⟩ := LRC14.lonely_of_pos_lacunary w hwpos hwchain
  refine ⟨t, (LRC14.lonely_abs_iff 14 v t).mp ?_⟩
  have ht' : Lonely 14 ((fun i ↦ |v i|) ∘ σ) t := by
    simpa [w, Function.comp_apply] using ht
  exact (LRC14.lonely_comp_equiv 14 (fun i ↦ |v i|) t σ).mp ht'

/-- **The sharpened dense-core obligation.**  Besides all grand-assembly residual
clauses and the chain-dichotomy normal form, the sorted absolute speeds must have an
adjacent ratio below `7/3`.  Thus the former overlap
`7/3 ≤ ratio < 3` at the last dense pair is no longer left to the residual. -/
def NonLacunaryDenseCoreObligation : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
    LRC14.CoveringFamily v →
    GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    (∃ i, 23 ≤ |v i|) →
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
      (Finset.univ.image k).card ≤ 12) →
    (∀ d : ℤ, 2 ≤ d → ∀ a : ℤ, (∀ i, d ∣ (v i - a)) → d ∣ a) →
    (∃ σ : Equiv.Perm (Fin 13), Monotone (fun i ↦ |v (σ i)|) ∧
      ChainDenseCore (fun i ↦ |v (σ i)|) ∧
      ¬ SevenThreeLacunary (fun i ↦ |v (σ i)|)) →
    ∃ t : ℝ, Lonely 14 v t

/-- A solver for only the non-7/3-lacunary dense core suffices for the old dense-core
obligation: the complementary branch is closed by nested safe intervals. -/
theorem denseCoreObligation_of_nonLacunary
    (hnonlac : NonLacunaryDenseCoreObligation) : DenseCoreObligation := by
  intro v hv hcov hgap hcomp hdist hlarge hdiv hcoarse hcres hcert
  obtain ⟨σ, hmono, hcore⟩ := hcert
  by_cases hlac : SevenThreeLacunary (fun i ↦ |v (σ i)|)
  · exact lonely_of_abs_perm_sevenThreeLacunary v hv σ hlac
  · exact hnonlac v hv hcov hgap hcomp hdist hlarge hdiv hcoarse hcres
      ⟨σ, hmono, hcore, hlac⟩

/-- **Grand-assembly residual after both chain reductions.**  The citation node
closes citable ratio-3 tails; nested gaps close every remaining 7/3-lacunary dense
core.  Only a `ChainDenseCore` with some adjacent ratio below 7/3 remains. -/
theorem residualObligation_of_nonLacunaryDenseCore (cite : LRCUpTo13)
    (hnonlac : NonLacunaryDenseCoreObligation) : ResidualObligation :=
  residualObligation_of_denseCore cite
    (denseCoreObligation_of_nonLacunary hnonlac)

/-- **LRC(14) from the sharpened non-lacunary dense-core obligation.** -/
theorem lrc14_of_nonLacunaryDenseCore (cite : LRCUpTo13)
    (hnonlac : NonLacunaryDenseCoreObligation) : LRC14.LRC14Statement :=
  lrc14_of_denseCore cite (denseCoreObligation_of_nonLacunary hnonlac)

/-! ## Axiom audit -/
#print axioms lonely_of_abs_perm_sevenThreeLacunary
#print axioms denseCoreObligation_of_nonLacunary
#print axioms residualObligation_of_nonLacunaryDenseCore
#print axioms lrc14_of_nonLacunaryDenseCore

end LRC14Grand
end LonelyRunner
