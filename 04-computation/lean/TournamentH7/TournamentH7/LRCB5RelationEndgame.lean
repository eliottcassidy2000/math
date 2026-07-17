import TournamentH7.LRCEndgameParameterDischargeTwoThree
import TournamentH7.LRCB5NormalizedBridge
import TournamentH7.LRCMomentCertificates
import TournamentH7.LRCArcWire
import TournamentH7.LRCDeepCount
import TournamentH7.LRCB5DeviationBudget

/-!
# Relation-budget certificates feed the chain-dense LRC(14) endgame

This module is the formal consumer for THM-935.  A certificate supplies a
modulus with at least one sampled phase, exact-support relation masses, the
exact `(q-1)`-scaled identity between their signed model and the concrete
integer `B5`, and the sharp one-sided horizon-thirty tail split.  The checked
budget then gives `B5 > 0`, which feeds the chain-dense endgame without any
additional analytic assumption.
-/

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open LRCB5RelationBudget
open LRCB5NormalizedBridge
open LRCB5CleanModulus
open scoped Classical

/-- A proof-producing THM-935 certificate for one concrete speed tuple. -/
structure B5RelationBudgetCertificate (v : Fin 13 → ℤ) where
  q : ℕ
  one_lt_q : 1 < q
  mass2 : ℝ
  mass3 : ℝ
  mass4 : ℝ
  mass5 : ℝ
  b5_eq_scaled_model :
    (LRC14Concrete.B5 v q : ℝ) =
      ((q : ℝ) - 1) * relationModel mass2 mass3 mass4 mass5
  pair_mass_lower_bound : -(13 / 30 : ℝ) ≤ mass2
  harmful_higher_budget :
    harmfulHigherContribution mass3 mass4 mass5 < 7712 / 84035

/-- Every relation-budget certificate produces the integer positivity witness
consumed by the discrete Bonferroni endgame. -/
theorem B5RelationBudgetCertificate.b5_pos {v : Fin 13 → ℤ}
    (certificate : B5RelationBudgetCertificate v) :
    0 < LRC14Concrete.B5 v certificate.q := by
  have hmodel : 0 < relationModel certificate.mass2 certificate.mass3
      certificate.mass4 certificate.mass5 :=
    relationModel_pos_of_signed_horizon_thirty_split _ _ _ _
      certificate.pair_mass_lower_bound certificate.harmful_higher_budget
  have hqR : (1 : ℝ) < certificate.q := by
    exact_mod_cast certificate.one_lt_q
  have hscale : (0 : ℝ) < (certificate.q : ℝ) - 1 := by
    linarith
  have hreal : (0 : ℝ) < (LRC14Concrete.B5 v certificate.q : ℝ) := by
    rw [certificate.b5_eq_scaled_model]
    positivity
  exact_mod_cast hreal

/-- Concrete THM-940 certificate at the canonical clean modulus.  The modulus,
coprimality, singleton sign, and finite-height modular-to-exact relation bridge
are automatic; only the second factorial-moment floor and harmful depth budget
remain. -/
structure NormalizedB5RelationBudgetCertificate (v : Fin 13 → ℤ) where
  height : ℕ
  coverage_excess_budget : 443 / 1470 ≤
    ((LRC14Concrete.liveCount v (cleanModulus v height) : ℚ) +
      (shiftedPairDepthMoment v (cleanModulus v height) : ℚ)) /
        ((cleanModulus v height : ℚ) - 1)
  harmful_depth_budget :
    harmfulDepthMoment v (cleanModulus v height) /
      ((cleanModulus v height : ℚ) - 1) < -(65218 / 84035)

/-- The normalized certificate feeds the concrete integer B5 directly; no
separate relation-mass identification premise remains. -/
theorem NormalizedB5RelationBudgetCertificate.b5_pos {v : Fin 13 → ℤ}
    (certificate : NormalizedB5RelationBudgetCertificate v)
    (hv : ∀ i, v i ≠ 0) :
    0 < LRC14Concrete.B5 v (cleanModulus v certificate.height) :=
  B5_pos_at_cleanModulus_of_coverage_excess_and_depth_budget
    v hv certificate.height certificate.coverage_excess_budget
      certificate.harmful_depth_budget

/-- **Exact usable-ruler window for a depth-six cap.**  At any capped modulus
`q > 1`, at most six speeds can lie strictly below `q/14`.  Otherwise the first
multiplier alone is simultaneously bad for seven runners.  Equivalently, a
cap-six search must keep `q` at most fourteen times at least seven of the
absolute speeds; cofinal rulers are excluded before relation analysis begins. -/
theorem small_speed_card_le_six_of_coverageCapped
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q)
    (hcap : LRC14Concrete.CoverageCapped v q 6) :
    ((Finset.univ : Finset (Fin 13)).filter fun i =>
      14 * (v i).natAbs < q).card ≤ 6 := by
  have hp : 1 ∈ Finset.Ioo 0 q := by
    simp only [Finset.mem_Ioo]
    exact ⟨by omega, hq⟩
  have hsubset :
      ((Finset.univ : Finset (Fin 13)).filter fun i => 14 * (v i).natAbs < q) ⊆
        (Finset.univ.filter fun i : Fin 13 => ¬ LRC14Concrete.inBand v q 1 i) := by
    intro i hi
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hi ⊢
    exact LRC14Concrete.not_inBand_one_of_fourteen_natAbs_lt v q i hi
  have hcard := Finset.card_le_card hsubset
  have hdepth := hcap 1 hp
  unfold LRC14Concrete.bandCount at hdepth
  omega

/-- **Clean-ruler near-zero catastrophe.**  At the cofinal clean modulus the
first multiplier is bad for all thirteen runners.  Consequently the clean
ruler can never be coverage-capped at depth six.  The clean relation bridge
and the cap-six bridge therefore require genuinely different modulus choices. -/
theorem not_coverageCapped_six_at_cleanModulus (v : Fin 13 → ℤ)
    (hv : ∀ i, v i ≠ 0) (height : ℕ) :
    ¬ LRC14Concrete.CoverageCapped v (cleanModulus v height) 6 := by
  intro hcap
  have hp : 1 ∈ Finset.Ioo 0 (cleanModulus v height) := by
    simp only [Finset.mem_Ioo]
    exact ⟨by omega, one_lt_cleanModulus v height hv⟩
  have hdepth := hcap 1 hp
  have hthirteen : LRC14Concrete.bandCount v (cleanModulus v height) 1 = 13 :=
    LRC14Concrete.bandCount_one_eq_thirteen_of_fourteen_natAbs_lt
      v (cleanModulus v height) (fun i => fourteen_natAbs_lt_cleanModulus v height hv i)
  omega

/-- THM-945 certificate at an arbitrary **usable** modulus.  On a depth-six
capped stratum, concrete B5 is exactly live multipliers minus depth-six
multipliers, so this strict census is both transparent and sharp.  The modulus
is deliberately not forced to be the cofinal clean ruler: the theorem above
shows that choice would make this structure empty. -/
structure CoverageCappedB5Certificate (v : Fin 13 → ℤ) where
  q : ℕ
  q_pos : 0 < q
  capped : LRC14Concrete.CoverageCapped v q 6
  live_gt_deep_six :
    ((Finset.Ioo 0 q).filter fun p =>
      LRC14Concrete.bandCount v q p = 6).card <
      LRC14Concrete.liveCount v q

/-- THM-947 constructor in the killer-arc language: exclude one rational
multiplier meeting seven bad arcs, then supply the exact THM-945 census. -/
def CoverageCappedB5Certificate.of_noSeven
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    (hnoSeven : ¬ ∃ p ∈ Finset.Ioo 0 q,
      7 ≤ LRC14Concrete.bandCount v q p)
    (hlive :
      ((Finset.Ioo 0 q).filter fun p =>
        LRC14Concrete.bandCount v q p = 6).card <
        LRC14Concrete.liveCount v q) :
    CoverageCappedB5Certificate v := by
  exact ⟨q, hq, (LRC14Concrete.coverageCapped_iff_no_seven v q).2 hnoSeven,
    hlive⟩

theorem CoverageCappedB5Certificate.b5_pos {v : Fin 13 → ℤ}
    (certificate : CoverageCappedB5Certificate v) :
    0 < LRC14Concrete.B5 v certificate.q := by
  have hcast :
      (((Finset.Ioo 0 certificate.q).filter fun p =>
        LRC14Concrete.bandCount v certificate.q p = 6).card : ℤ) <
        (LRC14Concrete.liveCount v certificate.q : ℤ) :=
    Int.ofNat_lt.mpr certificate.live_gt_deep_six
  rw [LRC14Concrete.B5_eq_live_sub_deepSix v
    certificate.q certificate.capped]
  omega
/-- A nonvacuous certificate stated directly in THM-940's concrete subset-
deviation ledger.  Unlike the relation-mass interface above, every term here
is defined from `jointFail`; favorable signed deviation is retained instead of
being charged through an absolute value. -/
structure B5DeviationBudgetCertificate (v : Fin 13 → ℤ) where
  q : ℕ
  one_lt_q : 1 < q
  deviation_surplus :
    -(((q : ℚ) - 1) * (2052 / 16807)) <
      LRC14Concrete.signedDeviationLedger v q

/-- THM-940 turns the concrete deviation inequality into integer `B5`
positivity without a relation-mass normalization or singular-series bridge. -/
theorem B5DeviationBudgetCertificate.b5_pos {v : Fin 13 → ℤ}
    (certificate : B5DeviationBudgetCertificate v) :
    0 < LRC14Concrete.B5 v certificate.q :=
  LRC14Concrete.B5_pos_of_signedDeviationLedger_lower v certificate.q
    certificate.deviation_surplus

/-- Unconditional THM-950 census certificate.  No coverage cap is assumed:
all multipliers of depth at least six are charged at the sharp universal
pointwise cost `792`. -/
structure CensusB5Certificate (v : Fin 13 → ℤ) where
  q : ℕ
  q_pos : 0 < q
  live_beats_deep :
    792 * (((Finset.Ioo 0 q).filter fun p =>
      6 ≤ LRC14Concrete.bandCount v q p).card : ℤ) <
      (LRC14Concrete.liveCount v q : ℤ)

/-- The unconditional census certificate produces the positive B5 modulus
consumed by the dense-core endgame. -/
theorem CensusB5Certificate.b5_pos {v : Fin 13 → ℤ}
    (certificate : CensusB5Certificate v) :
    0 < LRC14Concrete.B5 v certificate.q :=
  LRC14Concrete.B5_pos_of_live_beats_deep
    v certificate.q certificate.live_beats_deep

/-- THM-935 certificate supply only on the primitive, dissociated,
chain-dense core isolated by the current endgame. -/
def DenseCoreRelationBudgetSupply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.tupleGcd v = 1 →
    LRC14.CoveringFamily v → GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    (∃ i, 23 ≤ |v i|) →
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13,
      (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ),
      (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1 / 13 - 1 / 14 ∧
      (∀ i, k i ≠ 0) ∧ (Finset.univ.image k).card ≤ 12) →
    (∀ g : ℤ, 2 ≤ g →
      nonMultCard v g ≠ 2 ∧ nonMultCard v g ≠ 3) →
    (∃ σ : Equiv.Perm (Fin 13), Monotone (fun i => |v (σ i)|) ∧
      ChainDenseCore (fun i => |v (σ i)|)) →
    Nonempty (B5RelationBudgetCertificate v)

/-- Stronger concrete supplier on the same primitive, dissociated,
chain-dense core.  Its open mathematics is exactly a coprime modulus together
with the normalized pair and signed higher-support estimates. -/
def DenseCoreNormalizedRelationBudgetSupply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.tupleGcd v = 1 →
    LRC14.CoveringFamily v → GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    (∃ i, 23 ≤ |v i|) →
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13,
      (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ),
      (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1 / 13 - 1 / 14 ∧
      (∀ i, k i ≠ 0) ∧ (Finset.univ.image k).card ≤ 12) →
    (∀ g : ℤ, 2 ≤ g →
      nonMultCard v g ≠ 2 ∧ nonMultCard v g ≠ 3) →
    (∃ σ : Equiv.Perm (Fin 13), Monotone (fun i => |v (σ i)|) ∧
      ChainDenseCore (fun i => |v (σ i)|)) →
    Nonempty (NormalizedB5RelationBudgetCertificate v)

/-- THM-945 supply on the same dense core.  Its remaining mathematics is the
seven-wall cap together with the exact live-versus-depth-six census. -/
def DenseCoreCoverageCappedB5Supply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.tupleGcd v = 1 →
    LRC14.CoveringFamily v → GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    (∃ i, 23 ≤ |v i|) →
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13,
      (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ),
      (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1 / 13 - 1 / 14 ∧
      (∀ i, k i ≠ 0) ∧ (Finset.univ.image k).card ≤ 12) →
    (∀ g : ℤ, 2 ≤ g →
      nonMultCard v g ≠ 2 ∧ nonMultCard v g ≠ 3) →
    (∃ σ : Equiv.Perm (Fin 13), Monotone (fun i => |v (σ i)|) ∧
      ChainDenseCore (fun i => |v (σ i)|)) →
    Nonempty (CoverageCappedB5Certificate v)

/-- THM-950 supply on the same dense core.  This version removes the
seven-wall cap entirely; its sole combinatorial obligation is the explicit
live-versus-all-deep census. -/
def DenseCoreCensusB5Supply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.tupleGcd v = 1 →
    LRC14.CoveringFamily v → GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    (∃ i, 23 ≤ |v i|) →
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13,
      (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ),
      (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1 / 13 - 1 / 14 ∧
      (∀ i, k i ≠ 0) ∧ (Finset.univ.image k).card ≤ 12) →
    (∀ g : ℤ, 2 ≤ g →
      nonMultCard v g ≠ 2 ∧ nonMultCard v g ≠ 3) →
    (∃ σ : Equiv.Perm (Fin 13), Monotone (fun i => |v (σ i)|) ∧
      ChainDenseCore (fun i => |v (σ i)|)) →
    Nonempty (CensusB5Certificate v)

/- The direct THM-940 alternative on the same primitive, dissociated,
chain-dense core.  It asks only for the literal signed subset-deviation
surplus, with no relation-mass identification. -/
def DenseCoreDeviationBudgetSupply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.tupleGcd v = 1 →
    LRC14.CoveringFamily v → GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    (∃ i, 23 ≤ |v i|) →
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13,
      (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ),
      (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1 / 13 - 1 / 14 ∧
      (∀ i, k i ≠ 0) ∧ (Finset.univ.image k).card ≤ 12) →
    (∀ g : ℤ, 2 ≤ g →
      nonMultCard v g ≠ 2 ∧ nonMultCard v g ≠ 3) →
    (∃ σ : Equiv.Perm (Fin 13), Monotone (fun i => |v (σ i)|) ∧
      ChainDenseCore (fun i => |v (σ i)|)) →
    Nonempty (B5DeviationBudgetCertificate v)

/-- The structured relation-budget supplier discharges the raw positive-B5
supplier used by `LRCDenseCoreEndgame`. -/
theorem denseCoreDissociatedB5Supply_of_relationBudget
    (hsupply : DenseCoreRelationBudgetSupply) :
    DenseCoreDissociatedB5Supply := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  obtain ⟨certificate⟩ :=
    hsupply v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  exact ⟨certificate.q, lt_trans Nat.zero_lt_one certificate.one_lt_q,
    certificate.b5_pos⟩

/-- The normalized THM-940 supplier discharges the raw positive-B5 supplier
without routing through the older abstract equality certificate. -/
theorem denseCoreDissociatedB5Supply_of_normalizedRelationBudget
    (hsupply : DenseCoreNormalizedRelationBudgetSupply) :
    DenseCoreDissociatedB5Supply := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  obtain ⟨certificate⟩ :=
    hsupply v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  exact ⟨cleanModulus v certificate.height,
    lt_of_lt_of_le (by omega) (fourteen_le_cleanModulus v certificate.height hv),
    certificate.b5_pos hv⟩

/-- The exact cap-six census also discharges the raw dense-core B5 supplier. -/
theorem denseCoreDissociatedB5Supply_of_coverageCapped
    (hsupply : DenseCoreCoverageCappedB5Supply) :
    DenseCoreDissociatedB5Supply := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  obtain ⟨certificate⟩ :=
    hsupply v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  exact ⟨certificate.q, certificate.q_pos, certificate.b5_pos⟩

/-- The unconditional THM-950 census supplier also discharges the raw
dense-core B5 interface. -/
theorem denseCoreDissociatedB5Supply_of_census
    (hsupply : DenseCoreCensusB5Supply) :
    DenseCoreDissociatedB5Supply := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  obtain ⟨certificate⟩ :=
    hsupply v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  exact ⟨certificate.q, certificate.q_pos, certificate.b5_pos⟩

/-- The literal signed-deviation supply feeds the same positive-B5 consumer
without passing through normalized relation coordinates. -/
theorem denseCoreDissociatedB5Supply_of_deviationBudget
    (hsupply : DenseCoreDeviationBudgetSupply) :
    DenseCoreDissociatedB5Supply := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  obtain ⟨certificate⟩ :=
    hsupply v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  exact ⟨certificate.q, lt_trans Nat.zero_lt_one certificate.one_lt_q,
    certificate.b5_pos⟩

/-- **THM-935-shaped current capstone.**  The sanctioned LRC(≤13) citation,
the sharply reduced `q=2,3`/two-adic exceptional dispatch, and relation-budget
certificates only on the primitive dissociated chain-dense core imply LRC(14). -/
theorem lrc14_from_four_detuned_and_relationBudget
    (cite : LRCUpTo13) (hdeep : DeepExceptionalDetunedDispatchFour)
    (hsupply : DenseCoreRelationBudgetSupply) : LRC14.LRC14Statement :=
  lrc14_from_four_detuned_and_denseCore_dissociated_B5 cite hdeep
    (denseCoreDissociatedB5Supply_of_relationBudget hsupply)

/-- **Sharpest THM-935-shaped capstone.**  Degree arithmetic reduces the
three-detuned residual to a `q = 2` row with a distinct `q ≤ 8` companion or
the uniform `(3,3,3)` pattern; relation-budget certificates are required only
on the primitive dissociated chain-dense core. -/
theorem lrc14_from_twoThree_detuned_and_relationBudget
    (cite : LRCUpTo13) (hdeep : DeepExceptionalDetunedDispatchTwoThree)
    (hsupply : DenseCoreRelationBudgetSupply) : LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5 cite hdeep
    (denseCoreDissociatedB5Supply_of_relationBudget hsupply)

/-- Concrete normalized-THM-940 capstone.  The B5 side now exposes only the
sharp coverage-excess floor and signed higher-depth budget at a clean
modulus. -/
theorem lrc14_from_twoThree_detuned_and_normalizedRelationBudget
    (cite : LRCUpTo13) (hdeep : DeepExceptionalDetunedDispatchTwoThree)
    (hsupply : DenseCoreNormalizedRelationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5 cite hdeep
    (denseCoreDissociatedB5Supply_of_normalizedRelationBudget hsupply)

/-- THM-945 capstone.  On the dense core, a depth-six cap and the strict
live-versus-depth-six census suffice for the entire B5 branch. -/
theorem lrc14_from_twoThree_detuned_and_coverageCappedB5
    (cite : LRCUpTo13) (hdeep : DeepExceptionalDetunedDispatchTwoThree)
    (hsupply : DenseCoreCoverageCappedB5Supply) :
    LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5 cite hdeep
    (denseCoreDissociatedB5Supply_of_coverageCapped hsupply)

/-- THM-950 capstone.  The B5 branch is reduced to one explicit modulus whose
live multipliers outnumber `792` times every multiplier of depth at least six;
no separate coverage-cap hypothesis remains. -/
theorem lrc14_from_twoThree_detuned_and_censusB5
    (cite : LRCUpTo13) (hdeep : DeepExceptionalDetunedDispatchTwoThree)
    (hsupply : DenseCoreCensusB5Supply) :
    LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5 cite hdeep
    (denseCoreDissociatedB5Supply_of_census hsupply)

/-- Direct THM-940 capstone using the signed subset-deviation surplus. -/
theorem lrc14_from_twoThree_detuned_and_deviationBudget
    (cite : LRCUpTo13) (hdeep : DeepExceptionalDetunedDispatchTwoThree)
    (hsupply : DenseCoreDeviationBudgetSupply) : LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5 cite hdeep
    (denseCoreDissociatedB5Supply_of_deviationBudget hsupply)

/-! ## Axiom audit -/

#print axioms B5RelationBudgetCertificate.b5_pos
#print axioms NormalizedB5RelationBudgetCertificate.b5_pos
#print axioms small_speed_card_le_six_of_coverageCapped
#print axioms not_coverageCapped_six_at_cleanModulus
#print axioms CoverageCappedB5Certificate.b5_pos
#print axioms CoverageCappedB5Certificate.of_noSeven
#print axioms CensusB5Certificate.b5_pos
#print axioms denseCoreDissociatedB5Supply_of_relationBudget
#print axioms denseCoreDissociatedB5Supply_of_normalizedRelationBudget
#print axioms denseCoreDissociatedB5Supply_of_coverageCapped
#print axioms denseCoreDissociatedB5Supply_of_census
#print axioms lrc14_from_four_detuned_and_relationBudget
#print axioms lrc14_from_twoThree_detuned_and_relationBudget
#print axioms lrc14_from_twoThree_detuned_and_normalizedRelationBudget
#print axioms lrc14_from_twoThree_detuned_and_coverageCappedB5
#print axioms lrc14_from_twoThree_detuned_and_censusB5
#print axioms B5DeviationBudgetCertificate.b5_pos
#print axioms denseCoreDissociatedB5Supply_of_deviationBudget
#print axioms lrc14_from_twoThree_detuned_and_deviationBudget

end LRC14Grand
end LonelyRunner
