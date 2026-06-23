/-
  TournamentH7.Verify — Axiom dependency audit

  Inspect after `lake build` — the `#print axioms` output should list
  exactly the axioms documented in OCF.lean, Forbidden.lean, H21.lean,
  H63.lean, and Redei.lean (plus Lean's foundational `propext`,
  `Classical.choice`, `Quot.sound`).
-/

import TournamentH7.HSpectrum
import TournamentH7.RootSigns
import TournamentH7.RootPackets
import TournamentH7.NaturalOperationDigraphs
import TournamentH7.Tilings
import TournamentH7.GridReflection
import TournamentH7.StaircaseModel
import TournamentH7.SelfComplementary
import TournamentH7.AntiAutomorphism
import TournamentH7.HPIPIdentity
import TournamentH7.Iso
import TournamentH7.IsoProperties
import TournamentH7.SCCounts
import TournamentH7.SmallTournaments
import TournamentH7.ForbiddenHCounting
import TournamentH7.GoodCuts
import TournamentH7.StaircaseConnectivity
import TournamentH7.BucketBalance
import TournamentH7.PaleyAxiomatic
import TournamentH7.StaircaseTileModel
import TournamentH7.HSpectrumExtended
import TournamentH7.IsomorphismClasses
import TournamentH7.ApexBridge
import TournamentH7.StaircaseBucketTransport
import TournamentH7.RedeiFromOCF
import TournamentH7.HSpectrumClean
import TournamentH7.HSpectrumSmallN
import TournamentH7.SmallHEnumerations
import TournamentH7.BasePathSink
import TournamentH7.IsoCharacterizations
import TournamentH7.ScoreSequence
import TournamentH7.Paley3
import TournamentH7.TransitiveH
import TournamentH7.ProductSum
import TournamentH7.LRCDeathChain
import TournamentH7.LRCFactorialAtom
import TournamentH7.LRCBooleanTypeCut
import TournamentH7.LRCPeriodmaxCertificate
import TournamentH7.LRCGenuineWideCorrection
import TournamentH7.LRCQ6Contraction
import TournamentH7.LRCGk8SingleFar
import TournamentH7.LRCDoubletWitnessFloor
import TournamentH7.LRCLowerThresholdNeighborhood
import TournamentH7.LRCMreachConcrete
import TournamentH7.LRCApex7Floor
import TournamentH7.LRCUnitGrid14
import TournamentH7.LRCBindingPair
import TournamentH7.LRCGapReach
import TournamentH7.LRCWitnessAttainment
import TournamentH7.LRCWitnessAttainmentBridge
import TournamentH7.LRCMaxGapPigeonhole
import TournamentH7.LRCDenseCovers
import TournamentH7.LRCCoverBound
import TournamentH7.LRCMarginalUniform
import TournamentH7.LRCCoverBoxes
import TournamentH7.LRCArcComplexity
import TournamentH7.LRCGoodSet
import TournamentH7.LRCBonferroniMeasure
import TournamentH7.LRCEventMeasureBridge
import TournamentH7.LRCWitnessFloorConcrete
import TournamentH7.LRCWitnessBonferroni
import TournamentH7.LRCWitnessPartA
import TournamentH7.LRCL7Discrepancy
import TournamentH7.LRCApexShell
import TournamentH7.LRCTournamentStateLift

open Tournament

/-! ### LRC14 denominator-14 binding-pair arithmetic -/

theorem lrc_binding_pair_dvd_audit (a si sj : ℤ)
    (ha : IsCoprime (14 : ℤ) a)
    (hi : si * a ≡ 1 [ZMOD 14])
    (hj : sj * a ≡ -1 [ZMOD 14]) :
    (14 : ℤ) ∣ si + sj :=
  LonelyRunner.BindingPair.binding_pair_dvd a si sj ha hi hj
#print axioms lrc_binding_pair_dvd_audit

/-! ### LRC14 apex-shell arithmetic kernel -/

theorem lrc_apex_shell_fourteen_dvd_pair_sum_audit (u v a m n D : ℤ)
    (hcop : IsCoprime D a)
    (hpos : 14 * (u * a - m * D) = D)
    (hneg : 14 * (v * a - n * D) = -D) :
    (14 : ℤ) ∣ u + v :=
  LonelyRunner.ApexShell.fourteen_dvd_pair_sum_of_opposite_bindings
    u v a m n D hcop hpos hneg
#print axioms lrc_apex_shell_fourteen_dvd_pair_sum_audit

/-! ### LRC14 tournament-state lift closure endpoint -/

theorem lrc14_no_tournament_state_lift_audit
    (L : LonelyRunner.TournamentStateLift) : False :=
  LonelyRunner.no_tournament_state_lift L
#print axioms lrc14_no_tournament_state_lift_audit

theorem lrc14_not_bad_of_tournament_state_lift_audit {Bad : Prop}
    (hLift : Bad → LonelyRunner.TournamentStateLift) : ¬ Bad :=
  LonelyRunner.not_bad_of_tournament_state_lift hLift
#print axioms lrc14_not_bad_of_tournament_state_lift_audit

theorem lrc14_not_bad_of_H_eq_seven_lift_audit {Bad : Prop}
    (hLift : Bad → ∃ n : ℕ, ∃ T : Tournament n, H T = 7) : ¬ Bad :=
  LonelyRunner.not_bad_of_H_eq_seven_lift hLift
#print axioms lrc14_not_bad_of_H_eq_seven_lift_audit

/-! ### LRC14 HYP-2805 genuine-wide correction kernel -/

theorem lrc_genuine_wide_correction_all_reported_rows_below_cap_audit :
    ∀ r : Fin 5,
      0 < LonelyRunner.GenuineWideCorrection.marginNum
        (LonelyRunner.GenuineWideCorrection.maxRow r) :=
  LonelyRunner.GenuineWideCorrection.all_reported_rows_below_cap
#print axioms lrc_genuine_wide_correction_all_reported_rows_below_cap_audit

theorem lrc_genuine_wide_correction_k10_smallest_margin_audit :
    LonelyRunner.GenuineWideCorrection.marginLE LonelyRunner.GenuineWideCorrection.k10
        LonelyRunner.GenuineWideCorrection.k8 = true ∧
      LonelyRunner.GenuineWideCorrection.marginLE LonelyRunner.GenuineWideCorrection.k10
        LonelyRunner.GenuineWideCorrection.k9 = true ∧
      LonelyRunner.GenuineWideCorrection.marginLE LonelyRunner.GenuineWideCorrection.k10
        LonelyRunner.GenuineWideCorrection.k10 = true ∧
      LonelyRunner.GenuineWideCorrection.marginLE LonelyRunner.GenuineWideCorrection.k10
        LonelyRunner.GenuineWideCorrection.k11 = true ∧
      LonelyRunner.GenuineWideCorrection.marginLE LonelyRunner.GenuineWideCorrection.k10
        LonelyRunner.GenuineWideCorrection.k12 = true :=
  LonelyRunner.GenuineWideCorrection.k10_is_smallest_reported_margin
#print axioms lrc_genuine_wide_correction_k10_smallest_margin_audit

theorem lrc_genuine_wide_correction_robust_margin_flags_audit :
    0 <= LonelyRunner.GenuineWideCorrection.robustNum
        LonelyRunner.GenuineWideCorrection.k8 ∧
      0 <= LonelyRunner.GenuineWideCorrection.robustNum
        LonelyRunner.GenuineWideCorrection.k9 ∧
      LonelyRunner.GenuineWideCorrection.robustNum
        LonelyRunner.GenuineWideCorrection.k10 < 0 ∧
      0 <= LonelyRunner.GenuineWideCorrection.robustNum
        LonelyRunner.GenuineWideCorrection.k11 ∧
      0 <= LonelyRunner.GenuineWideCorrection.robustNum
        LonelyRunner.GenuineWideCorrection.k12 :=
  LonelyRunner.GenuineWideCorrection.robust_margin_flags
#print axioms lrc_genuine_wide_correction_robust_margin_flags_audit

theorem lrc_genuine_wide_correction_nonprimitive_base_guardrail_audit :
    LonelyRunner.GenuineWideCorrection.k9.basePrimitive = false ∧
      LonelyRunner.GenuineWideCorrection.k10.basePrimitive = false :=
  LonelyRunner.GenuineWideCorrection.nonprimitive_base_guardrail
#print axioms lrc_genuine_wide_correction_nonprimitive_base_guardrail_audit

/-! ### LRC14 doublet rho*/witness floor arithmetic checksum -/

theorem lrc_doublet_witness_floor_positive_and_separated_audit :
    0 < LonelyRunner.DoubletWitnessFloor.rhoFloorNum ∧
      0 < LonelyRunner.DoubletWitnessFloor.witnessFloorNum ∧
      LonelyRunner.DoubletWitnessFloor.consecRhoFloorNum *
          LonelyRunner.DoubletWitnessFloor.rhoFloorDen <
        LonelyRunner.DoubletWitnessFloor.rhoFloorNum *
          LonelyRunner.DoubletWitnessFloor.consecRhoFloorDen ∧
      LonelyRunner.DoubletWitnessFloor.thm530WitnessFloorNum *
          LonelyRunner.DoubletWitnessFloor.witnessFloorDen <
        LonelyRunner.DoubletWitnessFloor.witnessFloorNum *
          LonelyRunner.DoubletWitnessFloor.thm530WitnessFloorDen ∧
      LonelyRunner.DoubletWitnessFloor.rhoFloorNum *
          LonelyRunner.DoubletWitnessFloor.witnessFloorDen <
        LonelyRunner.DoubletWitnessFloor.witnessFloorNum *
          LonelyRunner.DoubletWitnessFloor.rhoFloorDen :=
  LonelyRunner.DoubletWitnessFloor.audited_doublet_floors_positive_and_separated
#print axioms lrc_doublet_witness_floor_positive_and_separated_audit

/-! ### LRC lower-threshold Farey-neighborhood arithmetic checksum -/

theorem lrc_lower_threshold_neighborhood_readouts_audit :
    LonelyRunner.LowerThresholdNeighborhood.properBoundedN6MinNum = 0 ∧
      LonelyRunner.LowerThresholdNeighborhood.properBoundedN8MinNum = 0 ∧
      LonelyRunner.LowerThresholdNeighborhood.properBoundedN10MinNum = 0 ∧
      LonelyRunner.LowerThresholdNeighborhood.properBoundedN12MinNum = 0 ∧
      LonelyRunner.LowerThresholdNeighborhood.properBoundedN14MinNum = 0 ∧
      0 < LonelyRunner.LowerThresholdNeighborhood.originBoundedN6MinNum ∧
      0 < LonelyRunner.LowerThresholdNeighborhood.originBoundedN8MinNum ∧
      0 < LonelyRunner.LowerThresholdNeighborhood.originBoundedN10MinNum ∧
      0 < LonelyRunner.LowerThresholdNeighborhood.originBoundedN12MinNum ∧
      0 < LonelyRunner.LowerThresholdNeighborhood.originBoundedN14MinNum ∧
      LonelyRunner.LowerThresholdNeighborhood.originBoundedN14MinNum *
          LonelyRunner.LowerThresholdNeighborhood.originBoundedN6MinDen ≤
        LonelyRunner.LowerThresholdNeighborhood.originBoundedN6MinNum *
          LonelyRunner.LowerThresholdNeighborhood.originBoundedN14MinDen ∧
      LonelyRunner.LowerThresholdNeighborhood.originBoundedN14MinNum *
          LonelyRunner.LowerThresholdNeighborhood.originBoundedN8MinDen ≤
        LonelyRunner.LowerThresholdNeighborhood.originBoundedN8MinNum *
          LonelyRunner.LowerThresholdNeighborhood.originBoundedN14MinDen ∧
      LonelyRunner.LowerThresholdNeighborhood.originBoundedN14MinNum *
          LonelyRunner.LowerThresholdNeighborhood.originBoundedN10MinDen ≤
        LonelyRunner.LowerThresholdNeighborhood.originBoundedN10MinNum *
          LonelyRunner.LowerThresholdNeighborhood.originBoundedN14MinDen ∧
      LonelyRunner.LowerThresholdNeighborhood.originBoundedN14MinNum *
          LonelyRunner.LowerThresholdNeighborhood.originBoundedN12MinDen ≤
        LonelyRunner.LowerThresholdNeighborhood.originBoundedN12MinNum *
          LonelyRunner.LowerThresholdNeighborhood.originBoundedN14MinDen :=
  LonelyRunner.LowerThresholdNeighborhood.audited_lower_threshold_neighborhood_readouts
#print axioms lrc_lower_threshold_neighborhood_readouts_audit

/-! ### LRC14 HYP-2829 gK8 single-far arithmetic checksum -/

theorem lrc_gk8_singlefar_binding_checks_audit :
    (2633 : Nat) * 588 < 2243 * 735 ∧
      (3259 : Nat) * 2002 < 9895 * 735 ∧
      (37 : Nat) * 7 < 40 * 7 ∧
      (2323 : Nat) * 588 < 2243 * 980 ∧
      (2876 : Nat) * 2002 < 9895 * 735 ∧
      (62267 : Nat) * 7 < 40 * 12936 ∧
      (2323 : Nat) * 735 < 2633 * 980 ∧
      (2876 : Nat) < 3259 ∧
      (62267 : Nat) * 7 < 37 * 12936 :=
  LonelyRunner.Gk8SingleFar.all_binding_checks
#print axioms lrc_gk8_singlefar_binding_checks_audit

/-! ### LRC14 concrete Mreach compactness bridge skeleton -/

theorem lrc_concrete_lonely_of_mreach_ge_audit
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hM : (1 : ℝ) / 14 ≤ LonelyRunner.LRC14Concrete.Mreach v) :
    ∃ t : ℝ, LonelyRunner.Lonely 14 v t :=
  LonelyRunner.LRC14Concrete.lonely_of_Mreach_ge v hv hM
#print axioms lrc_concrete_lonely_of_mreach_ge_audit

/-! ### LRC14 denominator-14 unit-grid sieve and apex floor ### -/

theorem lrc14_denom14_unit_lonely_iff_no_multiple_audit
    {ι : Type*} (a : ℤ) (hcop : IsCoprime (14 : ℤ) a) (v : ι → ℤ) :
    LonelyRunner.Lonely 14 v ((a : ℝ) / 14) ↔
      ∀ i, ¬ ((14 : ℤ) ∣ v i) :=
  LonelyRunner.UnitGrid14.denom14_unit_lonely_iff_no_multiple a hcop v
#print axioms lrc14_denom14_unit_lonely_iff_no_multiple_audit

theorem lrc14_counterexample_needs_multiple14_audit
    {ι : Type*} (v : ι → ℤ) (hcex : ∀ t : ℝ, ¬ LonelyRunner.Lonely 14 v t) :
    ∃ i, (14 : ℤ) ∣ v i :=
  LonelyRunner.UnitGrid14.counterexample_needs_multiple14 v hcex
#print axioms lrc14_counterexample_needs_multiple14_audit

theorem lrc14_apex7_D14_never_certifies_audit
    (S : Finset ℤ) (s : ℤ) (hsS : s ∈ S) (hs : (14 : ℤ) ∣ s) (a : ℤ) :
    ∃ t ∈ S,
      LonelyRunner.LRC14Concrete.nearInt ((t : ℝ) * ((a : ℝ) / 14)) < 1 / 14 :=
  LonelyRunner.Apex7.D14_never_certifies S s hsS hs a
#print axioms lrc14_apex7_D14_never_certifies_audit

/-! ### LRC14 geometric gap-to-reach core -/

theorem lrc_gapreach_margin_ge_of_free_interval_audit
    (C : Finset ℝ) (a g : ℝ)
    (hfree : ∀ c ∈ C, ∀ n : ℤ,
      (c + (n : ℝ)) ∉ Set.Ioo a (a + g)) :
    ∀ c ∈ C, ∀ n : ℤ,
      g / 2 ≤ |(a + g / 2) - (c + (n : ℝ))| :=
  LonelyRunner.GapReach.margin_ge_of_free_interval C a g hfree
#print axioms lrc_gapreach_margin_ge_of_free_interval_audit

theorem lrc_gapreach_margin_gt_one_div_14_of_gap_audit
    (C : Finset ℝ) (a g : ℝ) (hg : 1 / 7 < g)
    (hfree : ∀ c ∈ C, ∀ n : ℤ,
      (c + (n : ℝ)) ∉ Set.Ioo a (a + g)) :
    ∀ c ∈ C, ∀ n : ℤ,
      1 / 14 < |(a + g / 2) - (c + (n : ℝ))| :=
  LonelyRunner.GapReach.margin_gt_one_div_14_of_gap C a g hg hfree
#print axioms lrc_gapreach_margin_gt_one_div_14_of_gap_audit

theorem lrc_gapreach_exists_lonely_phase_of_gap_audit
    (C : Finset ℝ) (a g : ℝ) (hg : 1 / 7 < g)
    (hfree : ∀ c ∈ C, ∀ n : ℤ,
      (c + (n : ℝ)) ∉ Set.Ioo a (a + g)) :
    ∃ φ : ℝ, ∀ c ∈ C, ∀ n : ℤ,
      1 / 14 < |φ - (c + (n : ℝ))| :=
  LonelyRunner.GapReach.exists_lonely_phase_of_gap C a g hg hfree
#print axioms lrc_gapreach_exists_lonely_phase_of_gap_audit

theorem lrc_gapreach_nearInt_gt_of_forall_int_audit {y r : ℝ}
    (h : ∀ n : ℤ, r < |y - (n : ℝ)|) :
    r < LonelyRunner.LRC14Concrete.nearInt y :=
  LonelyRunner.GapReach.nearInt_gt_of_forall_int h
#print axioms lrc_gapreach_nearInt_gt_of_forall_int_audit

theorem lrc_gapreach_exists_nearInt_margin_of_gap_audit
    (C : Finset ℝ) (a g : ℝ) (hg : 1 / 7 < g)
    (hfree : ∀ c ∈ C, ∀ n : ℤ,
      (c + (n : ℝ)) ∉ Set.Ioo a (a + g)) :
    ∃ φ : ℝ, ∀ c ∈ C,
      1 / 14 < LonelyRunner.LRC14Concrete.nearInt (φ - c) :=
  LonelyRunner.GapReach.exists_nearInt_margin_of_gap C a g hg hfree
#print axioms lrc_gapreach_exists_nearInt_margin_of_gap_audit

/-! ### LRC14 witness-attainment margin/Mreach bridge -/

theorem lrc_witness_attainment_distZ_eq_nearInt_audit (x : ℝ) :
    TournamentH7.LRCWitness.distZ x = LonelyRunner.LRC14Concrete.nearInt x :=
  TournamentH7.LRCWitness.distZ_eq_nearInt x
#print axioms lrc_witness_attainment_distZ_eq_nearInt_audit

theorem lrc_witness_attainment_margin_sup_audit
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hmargin :
      (1 : ℝ) / 14 ≤
        sSup (TournamentH7.LRCWitness.margin v '' Set.Icc (0 : ℝ) 1)) :
    ∃ t : ℝ, LonelyRunner.Lonely 14 v t :=
  TournamentH7.LRCWitness.exists_lonely_of_margin_sSup_ge v hv hmargin
#print axioms lrc_witness_attainment_margin_sup_audit

/-! ### LRC14 max-gap pigeonhole hnu1 kernel -/

theorem lrc_maxgap_exists_one_div_card_le_audit
    {k : ℕ} (hk : 0 < k) (g : Fin k → ℝ)
    (hsum : ∑ i, g i = 1) :
    ∃ i, (1 : ℝ) / k ≤ g i :=
  TournamentH7.MaxGapPigeonhole.exists_one_div_card_le hk g hsum
#print axioms lrc_maxgap_exists_one_div_card_le_audit

theorem lrc_maxgap_exists_gap_gt_one_seventh_audit
    {k : ℕ} (hk : 0 < k) (hk6 : k ≤ 6)
    (g : Fin k → ℝ) (hsum : ∑ i, g i = 1) :
    ∃ i, (1 : ℝ) / 7 < g i :=
  TournamentH7.MaxGapPigeonhole.exists_gap_gt_one_seventh hk hk6 g hsum
#print axioms lrc_maxgap_exists_gap_gt_one_seventh_audit

theorem lrc_maxgap_all_eq_one_seventh_of_le_audit
    (g : Fin 7 → ℝ) (hsum : ∑ i, g i = 1)
    (hle : ∀ i, g i ≤ (1 : ℝ) / 7) :
    ∀ i, g i = (1 : ℝ) / 7 :=
  TournamentH7.MaxGapPigeonhole.all_eq_one_seventh_of_le g hsum hle
#print axioms lrc_maxgap_all_eq_one_seventh_of_le_audit

theorem lrc_maxgap_exists_gap_gt_or_all_eq_one_seventh_audit
    (g : Fin 7 → ℝ) (hsum : ∑ i, g i = 1) :
    (∃ i, (1 : ℝ) / 7 < g i) ∨ ∀ i, g i = (1 : ℝ) / 7 :=
  TournamentH7.MaxGapPigeonhole.exists_gap_gt_or_all_eq_one_seventh g hsum
#print axioms lrc_maxgap_exists_gap_gt_or_all_eq_one_seventh_audit

theorem lrc_goodSet_measurableSet_arc_audit (c : ℤ) :
    MeasurableSet {x : ℝ | Int.fract ((c : ℝ) * x) ∉ Set.Ioc (0 : ℝ) (1 / 7)} :=
  TournamentH7.GoodSet.measurableSet_arc c
#print axioms lrc_goodSet_measurableSet_arc_audit

theorem lrc_goodSet_measurableSet_goodSet_audit (E : List ℤ) :
    MeasurableSet (TournamentH7.GoodSet.goodSet E) :=
  TournamentH7.GoodSet.measurableSet_goodSet E
#print axioms lrc_goodSet_measurableSet_goodSet_audit

theorem lrc_goodSet_fract_sub_eq_fract_fract_sub_fract_audit (u v : ℝ) :
    Int.fract (u - v) = Int.fract (Int.fract u - Int.fract v) :=
  TournamentH7.GoodSet.fract_sub_eq_fract_fract_sub_fract u v
#print axioms lrc_goodSet_fract_sub_eq_fract_fract_sub_fract_audit

theorem lrc_goodSet_phaseGapSet_subset_goodSet_audit (E : List ℤ) :
    LonelyRunner.DenseCovers.phaseGapSet E ⊆ TournamentH7.GoodSet.goodSet E :=
  TournamentH7.GoodSet.phaseGapSet_subset_goodSet E
#print axioms lrc_goodSet_phaseGapSet_subset_goodSet_audit

theorem lrc_goodSet_denseSet_compl_subset_goodSet_audit (E : List ℤ) :
    (LonelyRunner.DenseCovers.denseSet E)ᶜ ⊆ TournamentH7.GoodSet.goodSet E :=
  TournamentH7.GoodSet.denseSet_compl_subset_goodSet E
#print axioms lrc_goodSet_denseSet_compl_subset_goodSet_audit

/-! ### LRC14 witness/p0 event-level elementary inclusions -/

theorem lrc_dense_covers_all_inner_audit
    (S : Finset ℝ) (hsub : ∀ s ∈ S, 0 ≤ s ∧ s < 1) (h0 : (0 : ℝ) ∈ S)
    (hdense : LonelyRunner.DenseCovers.Dense17 S) :
    ∀ j : ℕ, 1 ≤ j → j ≤ 6 →
      ∃ s ∈ S, (j : ℝ) / 7 ≤ s ∧ s < ((j : ℝ) + 1) / 7 :=
  LonelyRunner.DenseCovers.dense_covers_all_inner S hsub h0 hdense
#print axioms lrc_dense_covers_all_inner_audit

theorem lrc_exists_phase_arc_empty_of_not_dense_audit
    (S : Finset ℝ) (hsub : ∀ s ∈ S, 0 ≤ s ∧ s < 1)
    (hnd : ¬ LonelyRunner.DenseCovers.Dense17 S) :
    ∃ a ∈ S, ∀ c ∈ S,
      Int.fract (c - a) ∉ Set.Ioc (0 : ℝ) (1 / 7) :=
  LonelyRunner.DenseCovers.exists_phase_arc_empty_of_not_dense S hsub hnd
#print axioms lrc_exists_phase_arc_empty_of_not_dense_audit

theorem lrc_denseSet_subset_coverSet_audit
    (E : List ℤ) (hE : (0 : ℤ) ∈ E) :
    LonelyRunner.DenseCovers.denseSet E ⊆ LonelyRunner.DenseCovers.coverSet E :=
  LonelyRunner.DenseCovers.denseSet_subset_coverSet E hE
#print axioms lrc_denseSet_subset_coverSet_audit

theorem lrc_coverSet_compl_subset_denseSet_compl_audit
    (E : List ℤ) (hE : (0 : ℤ) ∈ E) :
    (LonelyRunner.DenseCovers.coverSet E)ᶜ ⊆
      (LonelyRunner.DenseCovers.denseSet E)ᶜ :=
  LonelyRunner.DenseCovers.coverSet_compl_subset_denseSet_compl E hE
#print axioms lrc_coverSet_compl_subset_denseSet_compl_audit

theorem lrc_denseSet_compl_subset_phaseGapSet_audit (E : List ℤ) :
    (LonelyRunner.DenseCovers.denseSet E)ᶜ ⊆
      LonelyRunner.DenseCovers.phaseGapSet E :=
  LonelyRunner.DenseCovers.denseSet_compl_subset_phaseGapSet E
#print axioms lrc_denseSet_compl_subset_phaseGapSet_audit

theorem lrc_volume_denseSet_le_coverSet_audit
    (E : List ℤ) (hE : (0 : ℤ) ∈ E) :
    MeasureTheory.volume
        (LonelyRunner.DenseCovers.denseSet E ∩ Set.Ico (0 : ℝ) 1) ≤
      MeasureTheory.volume
        (LonelyRunner.DenseCovers.coverSet E ∩ Set.Ico (0 : ℝ) 1) :=
  LonelyRunner.DenseCovers.volume_denseSet_le_coverSet E hE
#print axioms lrc_volume_denseSet_le_coverSet_audit

theorem lrc_volume_coverSet_inter_le_one_audit (E : List ℤ) :
    MeasureTheory.volume
        (LonelyRunner.DenseCovers.coverSet E ∩ Set.Ico (0 : ℝ) 1) ≤ 1 :=
  LonelyRunner.DenseCovers.volume_coverSet_inter_le_one E
#print axioms lrc_volume_coverSet_inter_le_one_audit

theorem lrc_measurableSet_coverSet_audit (E : List ℤ) :
    MeasurableSet (LonelyRunner.DenseCovers.coverSet E) :=
  LonelyRunner.DenseCovers.measurableSet_coverSet E
#print axioms lrc_measurableSet_coverSet_audit

theorem lrc_measurableSet_safeSet_audit (P : List ℤ) :
    MeasurableSet (LonelyRunner.DenseCovers.safeSet P) :=
  LonelyRunner.DenseCovers.measurableSet_safeSet P
#print axioms lrc_measurableSet_safeSet_audit

theorem lrc_slowμ_denseSet_le_coverSet_audit
    (E : List ℤ) (hE : (0 : ℤ) ∈ E) :
    LonelyRunner.DenseCovers.slowμ (LonelyRunner.DenseCovers.denseSet E) ≤
      LonelyRunner.DenseCovers.slowμ (LonelyRunner.DenseCovers.coverSet E) :=
  LonelyRunner.DenseCovers.slowμ_denseSet_le_coverSet E hE
#print axioms lrc_slowμ_denseSet_le_coverSet_audit

theorem lrc_bonferroni_toReal_audit
    {α : Type*} [MeasurableSpace α]
    (μ : MeasureTheory.Measure α) [MeasureTheory.IsProbabilityMeasure μ]
    (A B : Set α) (hB : MeasurableSet B) :
    (μ A).toReal + (μ B).toReal - 1 ≤ (μ (A ∩ B)).toReal :=
  LonelyRunner.BonferroniMeasure.toReal_bonferroni μ A B hB
#print axioms lrc_bonferroni_toReal_audit

theorem lrc_event_measure_bridge_bonferroni_audit
    {α : Type*} [MeasurableSpace α]
    (μ : MeasureTheory.Measure α) [MeasureTheory.IsProbabilityMeasure μ]
    (GOOD GP : LonelyRunner.LRC14.Shape → Set α)
    (hGP : ∀ s, MeasurableSet (GP s))
    (nuShape measGP : LonelyRunner.LRC14.Shape → ℝ)
    (hnu : ∀ s, nuShape s = (μ (GOOD s)).toReal)
    (hgp : ∀ s, measGP s = (μ (GP s)).toReal)
    (hwitness :
      ∀ s, LonelyRunner.LRC14.witnessG2 s =
        (μ ((GOOD s) ∩ (GP s))).toReal) :
    ∀ s, nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s :=
  LonelyRunner.LRC14.EventMeasureBridge.shape_bonferroni_handoff
    μ GOOD GP hGP nuShape measGP hnu hgp hwitness
#print axioms lrc_event_measure_bridge_bonferroni_audit

theorem lrc_event_measure_bridge_D_le_p0_audit
    {α : Type*} [MeasurableSpace α]
    (μ : MeasureTheory.Measure α) [MeasureTheory.IsProbabilityMeasure μ]
    (Dset P0set : LonelyRunner.LRC14.Shape → Set α)
    (nuShape DShape p0Shape : LonelyRunner.LRC14.Shape → ℝ)
    (hsub : ∀ s, Dset s ⊆ P0set s)
    (hDmeasure : ∀ s, DShape s = (μ (Dset s)).toReal)
    (hp0measure : ∀ s, p0Shape s = (μ (P0set s)).toReal)
    (hDdef : ∀ s, DShape s = 1 - nuShape s) :
    ∀ s, (1 - nuShape s) ≤ p0Shape s :=
  LonelyRunner.LRC14.EventMeasureBridge.shape_D_le_p0_handoff
    μ Dset P0set nuShape DShape p0Shape hsub hDmeasure hp0measure hDdef
#print axioms lrc_event_measure_bridge_D_le_p0_audit

theorem lrc_event_measure_bridge_safeSet_bonferroni_audit
    (GOOD : LonelyRunner.LRC14.Shape → Set ℝ)
    (Pof : LonelyRunner.LRC14.Shape → List ℤ)
    (nuShape measGP : LonelyRunner.LRC14.Shape → ℝ)
    (hnu : ∀ s, nuShape s =
      (LonelyRunner.DenseCovers.slowμ (GOOD s)).toReal)
    (hgp : ∀ s, measGP s =
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (hwitness :
      ∀ s, LonelyRunner.LRC14.witnessG2 s =
        (LonelyRunner.DenseCovers.slowμ
          ((GOOD s) ∩ LonelyRunner.DenseCovers.safeSet (Pof s))).toReal) :
    ∀ s, nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s :=
  LonelyRunner.LRC14.EventMeasureBridge.shape_bonferroni_safeSet_handoff
    GOOD Pof nuShape measGP hnu hgp hwitness
#print axioms lrc_event_measure_bridge_safeSet_bonferroni_audit

theorem lrc_event_measure_bridge_denseCovers_D_le_p0_audit
    (Eof : LonelyRunner.LRC14.Shape → List ℤ)
    (hanchor : ∀ s, (0 : ℤ) ∈ Eof s)
    (nuShape DShape p0Shape : LonelyRunner.LRC14.Shape → ℝ)
    (hDmeasure : ∀ s, DShape s =
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.denseSet (Eof s))).toReal)
    (hp0measure : ∀ s, p0Shape s =
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet (Eof s))).toReal)
    (hDdef : ∀ s, DShape s = 1 - nuShape s) :
    ∀ s, (1 - nuShape s) ≤ p0Shape s :=
  LonelyRunner.LRC14.EventMeasureBridge.shape_D_le_p0_denseCovers_handoff
    Eof hanchor nuShape DShape p0Shape hDmeasure hp0measure hDdef
#print axioms lrc_event_measure_bridge_denseCovers_D_le_p0_audit

theorem lrc_event_measure_bridge_goodSet_safeSet_bonferroni_audit
    (Eof Pof : LonelyRunner.LRC14.Shape → List ℤ)
    (nuShape measGP : LonelyRunner.LRC14.Shape → ℝ)
    (hnu : ∀ s, nuShape s =
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet (Eof s))).toReal)
    (hgp : ∀ s, measGP s =
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (hwitness : ∀ s, LonelyRunner.LRC14.witnessG2 s =
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet (Eof s) ∩
          LonelyRunner.DenseCovers.safeSet (Pof s))).toReal) :
    ∀ s, nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s :=
  LonelyRunner.LRC14.EventMeasureBridge.shape_bonferroni_goodSet_safeSet_handoff
    Eof Pof nuShape measGP hnu hgp hwitness
#print axioms lrc_event_measure_bridge_goodSet_safeSet_bonferroni_audit

theorem lrc_event_measure_bridge_goodSet_witness_pos_from_strict_cover_bound_audit
    (Eof Pof : LonelyRunner.LRC14.Shape → List ℤ)
    (cap : LonelyRunner.LRC14.Shape → ℝ)
    (hanchor : ∀ s, (0 : ℤ) ∈ Eof s)
    (hwitness : ∀ s, LonelyRunner.LRC14.witnessG2 s =
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet (Eof s) ∩
          LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (hwide : ∀ s,
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet (Eof s))).toReal < cap s)
    (hdual : ∀ s, cap s ≤
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.safeSet (Pof s))).toReal) :
    ∀ s, 0 < LonelyRunner.LRC14.witnessG2 s :=
  LonelyRunner.LRC14.EventMeasureBridge.shape_goodSet_witness_pos_from_strict_cover_bound
    Eof Pof cap hanchor hwitness hwide hdual
#print axioms lrc_event_measure_bridge_goodSet_witness_pos_from_strict_cover_bound_audit

theorem lrc_event_measure_bridge_goodSet_witness_margin_from_wide_bound_audit
    (Eof Pof : LonelyRunner.LRC14.Shape → List ℤ)
    (cap delta : LonelyRunner.LRC14.Shape → ℝ)
    (hanchor : ∀ s, (0 : ℤ) ∈ Eof s)
    (hwitness : ∀ s, LonelyRunner.LRC14.witnessG2 s =
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet (Eof s) ∩
          LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (hwide : ∀ s,
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet (Eof s))).toReal ≤
          cap s - delta s)
    (hdual : ∀ s, cap s ≤
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.safeSet (Pof s))).toReal) :
    ∀ s, delta s ≤ LonelyRunner.LRC14.witnessG2 s :=
  LonelyRunner.LRC14.EventMeasureBridge.shape_goodSet_witness_margin_from_wide_bound
    Eof Pof cap delta hanchor hwitness hwide hdual
#print axioms lrc_event_measure_bridge_goodSet_witness_margin_from_wide_bound_audit

/-! ### LRC14 p0 cover-bound elementary cores -/

theorem lrc_cover_bound_coverSet_mono_audit {E E' : List ℤ}
    (h : ∀ e ∈ E, e ∈ E') :
    LonelyRunner.DenseCovers.coverSet E ⊆
      LonelyRunner.DenseCovers.coverSet E' :=
  LonelyRunner.DenseCovers.coverSet_mono h
#print axioms lrc_cover_bound_coverSet_mono_audit

theorem lrc_cover_bound_slowμ_coverSet_mono_audit {E E' : List ℤ}
    (h : ∀ e ∈ E, e ∈ E') :
    LonelyRunner.DenseCovers.slowμ (LonelyRunner.DenseCovers.coverSet E) ≤
      LonelyRunner.DenseCovers.slowμ (LonelyRunner.DenseCovers.coverSet E') :=
  LonelyRunner.DenseCovers.slowμ_coverSet_mono h
#print axioms lrc_cover_bound_slowμ_coverSet_mono_audit

theorem lrc_cover_bound_six_le_card_of_coverSet_mem_audit
    {E : List ℤ} {x : ℝ}
    (hx : x ∈ LonelyRunner.DenseCovers.coverSet E) :
    6 ≤ E.toFinset.card :=
  LonelyRunner.DenseCovers.six_le_card_of_coverSet_mem hx
#print axioms lrc_cover_bound_six_le_card_of_coverSet_mem_audit

theorem lrc_cover_bound_coverSet_eq_empty_of_card_lt_six_audit
    {E : List ℤ} (h : E.toFinset.card < 6) :
    LonelyRunner.DenseCovers.coverSet E = ∅ :=
  LonelyRunner.DenseCovers.coverSet_eq_empty_of_card_lt_six h
#print axioms lrc_cover_bound_coverSet_eq_empty_of_card_lt_six_audit

theorem lrc_cover_bound_slowμ_coverSet_eq_zero_of_card_lt_six_audit
    {E : List ℤ} (h : E.toFinset.card < 6) :
    LonelyRunner.DenseCovers.slowμ
      (LonelyRunner.DenseCovers.coverSet E) = 0 :=
  LonelyRunner.DenseCovers.slowμ_coverSet_eq_zero_of_card_lt_six h
#print axioms lrc_cover_bound_slowμ_coverSet_eq_zero_of_card_lt_six_audit

theorem lrc_cover_bound_slowμ_coverSet_lt_cap_of_decorrelation_audit
    (E : List ℤ) (p0decorr Q capv : ℝ)
    (hresonance :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal ≤ p0decorr)
    (hfinite : p0decorr ≤ Q)
    (hmargin : Q < capv) :
    (LonelyRunner.DenseCovers.slowμ
      (LonelyRunner.DenseCovers.coverSet E)).toReal < capv :=
  LonelyRunner.DenseCovers.slowμ_coverSet_lt_cap_of_decorrelation
    E p0decorr Q capv hresonance hfinite hmargin
#print axioms lrc_cover_bound_slowμ_coverSet_lt_cap_of_decorrelation_audit

theorem lrc_cover_bound_slowμ_coverSet_lt_cap_audit
    (E : List ℤ) (p0decorr Q capv : ℝ) (hcap : 0 < capv)
    (hbig : 6 ≤ E.toFinset.card →
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal ≤ p0decorr ∧
        p0decorr ≤ Q ∧ Q < capv) :
    (LonelyRunner.DenseCovers.slowμ
      (LonelyRunner.DenseCovers.coverSet E)).toReal < capv :=
  LonelyRunner.DenseCovers.slowμ_coverSet_lt_cap
    E p0decorr Q capv hcap hbig
#print axioms lrc_cover_bound_slowμ_coverSet_lt_cap_audit

theorem lrc_marginal_uniform_slowμ_fract_Ico_le_audit
    (w : ℕ) (hw : 1 ≤ w) {a b : ℝ}
    (ha : 0 ≤ a) (hab : a ≤ b) (hb : b ≤ 1) :
    LonelyRunner.DenseCovers.slowμ
        {x | Int.fract ((w : ℝ) * x) ∈ Set.Ico a b} ≤
      ENNReal.ofReal (b - a) :=
  LonelyRunner.MarginalUniform.slowμ_fract_Ico_le w hw ha hab hb
#print axioms lrc_marginal_uniform_slowμ_fract_Ico_le_audit

theorem lrc_marginal_uniform_slowμ_fract_sector_le_audit
    (w : ℕ) (hw : 1 ≤ w) (j : ℕ) (hj : j ≤ 6) :
    LonelyRunner.DenseCovers.slowμ
        {x | Int.fract ((w : ℝ) * x) ∈
          Set.Ico ((j : ℝ) / 7) (((j : ℝ) + 1) / 7)} ≤
      ENNReal.ofReal (1 / 7) :=
  LonelyRunner.MarginalUniform.slowμ_fract_sector_le w hw j hj
#print axioms lrc_marginal_uniform_slowμ_fract_sector_le_audit

theorem lrc_cover_boxes_coverSet_subset_iUnion_sectorBox_audit
    (E : List ℤ) :
    LonelyRunner.DenseCovers.coverSet E ⊆
      ⋃ σ : LonelyRunner.CoverBoxes.Assignment E,
        LonelyRunner.CoverBoxes.sectorBox E σ :=
  LonelyRunner.CoverBoxes.coverSet_subset_iUnion_sectorBox E
#print axioms lrc_cover_boxes_coverSet_subset_iUnion_sectorBox_audit

theorem lrc_cover_boxes_slowμ_coverSet_le_tsum_sectorBox_audit
    (E : List ℤ) :
    LonelyRunner.DenseCovers.slowμ (LonelyRunner.DenseCovers.coverSet E) ≤
      ∑' σ : LonelyRunner.CoverBoxes.Assignment E,
        LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.CoverBoxes.sectorBox E σ) :=
  LonelyRunner.CoverBoxes.slowμ_coverSet_le_tsum_sectorBox E
#print axioms lrc_cover_boxes_slowμ_coverSet_le_tsum_sectorBox_audit

theorem lrc_arc_complexity_occupiedCount_le_cells_audit
    {ι α : Type*} [DecidableEq α]
    (I : Finset ι) (cells : Finset α) (B : ι → Finset α)
    (hsub : ∀ i ∈ I, B i ⊆ cells)
    (hdisj : (I : Set ι).PairwiseDisjoint B) :
    LonelyRunner.ArcComplexity.FiniteCell.occupiedCount I B ≤ cells.card :=
  LonelyRunner.ArcComplexity.FiniteCell.occupiedCount_le_cells
    I cells B hsub hdisj
#print axioms lrc_arc_complexity_occupiedCount_le_cells_audit

theorem lrc_arc_complexity_occupiedCount_le_seven_mul_sum_audit
    {ι α : Type*} [DecidableEq α]
    (I : Finset ι) (cells : Finset α) (B : ι → Finset α) (sumE : ℕ)
    (hsub : ∀ i ∈ I, B i ⊆ cells)
    (hdisj : (I : Set ι).PairwiseDisjoint B)
    (hcells : cells.card ≤ 7 * sumE) :
    LonelyRunner.ArcComplexity.FiniteCell.occupiedCount I B ≤ 7 * sumE :=
  LonelyRunner.ArcComplexity.FiniteCell.occupiedCount_le_seven_mul_sum
    I cells B sumE hsub hdisj hcells
#print axioms lrc_arc_complexity_occupiedCount_le_seven_mul_sum_audit

/-! ### LRC14 concrete slow-time witness floor -/

theorem lrc_slowμ_compl_toReal_audit {A : Set ℝ} (hA : MeasurableSet A) :
    (LonelyRunner.DenseCovers.slowμ Aᶜ).toReal =
      1 - (LonelyRunner.DenseCovers.slowμ A).toReal :=
  LonelyRunner.DenseCovers.slowμ_compl_toReal hA
#print axioms lrc_slowμ_compl_toReal_audit

theorem lrc_witness_floor_concrete_audit (E P : List ℤ) :
    (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.safeSet P)).toReal -
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal
      ≤ (LonelyRunner.DenseCovers.slowμ
        ((LonelyRunner.DenseCovers.coverSet E)ᶜ ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.witness_floor_concrete E P
#print axioms lrc_witness_floor_concrete_audit

theorem lrc_witness_pos_of_p0_lt_measGP_audit (E P : List ℤ)
    (h :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal <
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    0 <
      (LonelyRunner.DenseCovers.slowμ
        ((LonelyRunner.DenseCovers.coverSet E)ᶜ ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.witness_pos_of_p0_lt_measGP E P h
#print axioms lrc_witness_pos_of_p0_lt_measGP_audit

theorem lrc_witness_carrier_subset_safe_audit (E P : List ℤ) :
    (LonelyRunner.DenseCovers.coverSet E)ᶜ ∩
        LonelyRunner.DenseCovers.safeSet P ⊆
      LonelyRunner.DenseCovers.safeSet P :=
  LonelyRunner.DenseCovers.witness_carrier_subset_safe E P
#print axioms lrc_witness_carrier_subset_safe_audit

theorem lrc_witness_carrier_subset_dense_compl_audit
    (E P : List ℤ) (hE : (0 : ℤ) ∈ E) :
    (LonelyRunner.DenseCovers.coverSet E)ᶜ ∩
        LonelyRunner.DenseCovers.safeSet P ⊆
      (LonelyRunner.DenseCovers.denseSet E)ᶜ ∩
        LonelyRunner.DenseCovers.safeSet P :=
  LonelyRunner.DenseCovers.witness_carrier_subset_dense_compl E P hE
#print axioms lrc_witness_carrier_subset_dense_compl_audit

theorem lrc_witness_carrier_le_dense_compl_measure_audit
    (E P : List ℤ) (hE : (0 : ℤ) ∈ E) :
    (LonelyRunner.DenseCovers.slowμ
        ((LonelyRunner.DenseCovers.coverSet E)ᶜ ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal ≤
      (LonelyRunner.DenseCovers.slowμ
        ((LonelyRunner.DenseCovers.denseSet E)ᶜ ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.witness_carrier_le_dense_compl_measure E P hE
#print axioms lrc_witness_carrier_le_dense_compl_measure_audit

theorem lrc_dense_compl_witness_subset_phaseGap_audit (E P : List ℤ) :
    (LonelyRunner.DenseCovers.denseSet E)ᶜ ∩
        LonelyRunner.DenseCovers.safeSet P ⊆
      LonelyRunner.DenseCovers.phaseGapSet E ∩
        LonelyRunner.DenseCovers.safeSet P :=
  LonelyRunner.DenseCovers.dense_compl_witness_subset_phaseGap E P
#print axioms lrc_dense_compl_witness_subset_phaseGap_audit

theorem lrc_dense_compl_witness_le_phaseGap_measure_audit (E P : List ℤ) :
    (LonelyRunner.DenseCovers.slowμ
        ((LonelyRunner.DenseCovers.denseSet E)ᶜ ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal ≤
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.phaseGapSet E ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.dense_compl_witness_le_phaseGap_measure E P
#print axioms lrc_dense_compl_witness_le_phaseGap_measure_audit

theorem lrc_phaseGap_witness_subset_goodSet_audit (E P : List ℤ) :
    LonelyRunner.DenseCovers.phaseGapSet E ∩
        LonelyRunner.DenseCovers.safeSet P ⊆
      TournamentH7.GoodSet.goodSet E ∩
        LonelyRunner.DenseCovers.safeSet P :=
  LonelyRunner.DenseCovers.phaseGap_witness_subset_goodSet E P
#print axioms lrc_phaseGap_witness_subset_goodSet_audit

theorem lrc_phaseGap_witness_le_goodSet_measure_audit (E P : List ℤ) :
    (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.phaseGapSet E ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal ≤
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet E ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.phaseGap_witness_le_goodSet_measure E P
#print axioms lrc_phaseGap_witness_le_goodSet_measure_audit

theorem lrc_witness_pos_from_wide_bound_audit (E P : List ℤ) (capk : ℝ)
    (hwide :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal ≤ capk)
    (hdual :
      capk <
        (LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    0 <
      (LonelyRunner.DenseCovers.slowμ
        ((LonelyRunner.DenseCovers.coverSet E)ᶜ ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.witness_pos_from_wide_bound E P capk hwide hdual
#print axioms lrc_witness_pos_from_wide_bound_audit

theorem lrc_witness_pos_from_strict_cover_bound_audit (E P : List ℤ) (capk : ℝ)
    (hwide :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal < capk)
    (hdual :
      capk ≤
        (LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    0 <
      (LonelyRunner.DenseCovers.slowμ
        ((LonelyRunner.DenseCovers.coverSet E)ᶜ ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.witness_pos_from_strict_cover_bound E P capk hwide hdual
#print axioms lrc_witness_pos_from_strict_cover_bound_audit

theorem lrc_dense_compl_witness_pos_from_strict_cover_bound_audit
    (E P : List ℤ) (capk : ℝ) (hE : (0 : ℤ) ∈ E)
    (hwide :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal < capk)
    (hdual :
      capk ≤
        (LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    0 <
      (LonelyRunner.DenseCovers.slowμ
        ((LonelyRunner.DenseCovers.denseSet E)ᶜ ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.dense_compl_witness_pos_from_strict_cover_bound
    E P capk hE hwide hdual
#print axioms lrc_dense_compl_witness_pos_from_strict_cover_bound_audit

theorem lrc_phaseGap_witness_pos_from_strict_cover_bound_audit
    (E P : List ℤ) (capk : ℝ) (hE : (0 : ℤ) ∈ E)
    (hwide :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal < capk)
    (hdual :
      capk ≤
        (LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    0 <
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.phaseGapSet E ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.phaseGap_witness_pos_from_strict_cover_bound
    E P capk hE hwide hdual
#print axioms lrc_phaseGap_witness_pos_from_strict_cover_bound_audit

theorem lrc_goodSet_witness_pos_from_strict_cover_bound_audit
    (E P : List ℤ) (capk : ℝ) (hE : (0 : ℤ) ∈ E)
    (hwide :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal < capk)
    (hdual :
      capk ≤
        (LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    0 <
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet E ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.goodSet_witness_pos_from_strict_cover_bound
    E P capk hE hwide hdual
#print axioms lrc_goodSet_witness_pos_from_strict_cover_bound_audit

theorem lrc_witness_margin_from_wide_bound_audit
    (E P : List ℤ) (capk delta : ℝ)
    (hwide :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal ≤ capk - delta)
    (hdual :
      capk ≤
        (LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    delta ≤
      (LonelyRunner.DenseCovers.slowμ
        ((LonelyRunner.DenseCovers.coverSet E)ᶜ ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.witness_margin_from_wide_bound
    E P capk delta hwide hdual
#print axioms lrc_witness_margin_from_wide_bound_audit

theorem lrc_dense_compl_witness_margin_from_wide_bound_audit
    (E P : List ℤ) (capk delta : ℝ) (hE : (0 : ℤ) ∈ E)
    (hwide :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal ≤ capk - delta)
    (hdual :
      capk ≤
        (LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    delta ≤
      (LonelyRunner.DenseCovers.slowμ
        ((LonelyRunner.DenseCovers.denseSet E)ᶜ ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.dense_compl_witness_margin_from_wide_bound
    E P capk delta hE hwide hdual
#print axioms lrc_dense_compl_witness_margin_from_wide_bound_audit

theorem lrc_phaseGap_witness_margin_from_wide_bound_audit
    (E P : List ℤ) (capk delta : ℝ) (hE : (0 : ℤ) ∈ E)
    (hwide :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal ≤ capk - delta)
    (hdual :
      capk ≤
        (LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    delta ≤
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.phaseGapSet E ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.phaseGap_witness_margin_from_wide_bound
    E P capk delta hE hwide hdual
#print axioms lrc_phaseGap_witness_margin_from_wide_bound_audit

theorem lrc_goodSet_witness_margin_from_wide_bound_audit
    (E P : List ℤ) (capk delta : ℝ) (hE : (0 : ℤ) ∈ E)
    (hwide :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal ≤ capk - delta)
    (hdual :
      capk ≤
        (LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    delta ≤
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet E ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.goodSet_witness_margin_from_wide_bound
    E P capk delta hE hwide hdual
#print axioms lrc_goodSet_witness_margin_from_wide_bound_audit

theorem lrc_witness_pos_from_wide_bound_margin_audit
    (E P : List ℤ) (capk delta : ℝ)
    (hdelta : 0 < delta)
    (hwide :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal ≤ capk - delta)
    (hdual :
      capk ≤
        (LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    0 <
      (LonelyRunner.DenseCovers.slowμ
        ((LonelyRunner.DenseCovers.coverSet E)ᶜ ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.witness_pos_from_wide_bound_margin
    E P capk delta hdelta hwide hdual
#print axioms lrc_witness_pos_from_wide_bound_margin_audit

theorem lrc_dense_compl_witness_pos_from_wide_bound_margin_audit
    (E P : List ℤ) (capk delta : ℝ) (hE : (0 : ℤ) ∈ E)
    (hdelta : 0 < delta)
    (hwide :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal ≤ capk - delta)
    (hdual :
      capk ≤
        (LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    0 <
      (LonelyRunner.DenseCovers.slowμ
        ((LonelyRunner.DenseCovers.denseSet E)ᶜ ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.dense_compl_witness_pos_from_wide_bound_margin
    E P capk delta hE hdelta hwide hdual
#print axioms lrc_dense_compl_witness_pos_from_wide_bound_margin_audit

theorem lrc_phaseGap_witness_pos_from_wide_bound_margin_audit
    (E P : List ℤ) (capk delta : ℝ) (hE : (0 : ℤ) ∈ E)
    (hdelta : 0 < delta)
    (hwide :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal ≤ capk - delta)
    (hdual :
      capk ≤
        (LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    0 <
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.phaseGapSet E ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.phaseGap_witness_pos_from_wide_bound_margin
    E P capk delta hE hdelta hwide hdual
#print axioms lrc_phaseGap_witness_pos_from_wide_bound_margin_audit

theorem lrc_goodSet_witness_pos_from_wide_bound_margin_audit
    (E P : List ℤ) (capk delta : ℝ) (hE : (0 : ℤ) ∈ E)
    (hdelta : 0 < delta)
    (hwide :
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet E)).toReal ≤ capk - delta)
    (hdual :
      capk ≤
        (LonelyRunner.DenseCovers.slowμ
          (LonelyRunner.DenseCovers.safeSet P)).toReal) :
    0 <
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet E ∩
          LonelyRunner.DenseCovers.safeSet P)).toReal :=
  LonelyRunner.DenseCovers.goodSet_witness_pos_from_wide_bound_margin
    E P capk delta hE hdelta hwide hdual
#print axioms lrc_goodSet_witness_pos_from_wide_bound_margin_audit

/-! ### LRC14 Bonferroni positive-p0 route -/

theorem lrc_bonferroni_witness_margin_from_p0_wide_bound_shapes_audit
    (nuShape measGP p0Shape cap delta : LonelyRunner.LRC14.Shape → ℝ)
    (hbonf :
      ∀ s, nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, cap s ≤ measGP s)
    (s : LonelyRunner.LRC14.Shape) :
    delta s ≤ LonelyRunner.LRC14.witnessG2 s :=
  LonelyRunner.LRC14.Bonferroni.witness_margin_from_p0_wide_bound_shapes
    nuShape measGP p0Shape cap delta hbonf hDp0 hp0cap hmeasGP s
#print axioms lrc_bonferroni_witness_margin_from_p0_wide_bound_shapes_audit

theorem lrc_bonferroni_large_witness_pos_from_p0_wide_bound_shapes_audit
    (nuShape measGP p0Shape cap delta : LonelyRunner.LRC14.Shape → ℝ)
    (hδpos : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 → 0 < delta s)
    (hbonf : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s)
    (hDp0  : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      cap s ≤ measGP s) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      8 ≤ LonelyRunner.LRC14.clusterSize (LonelyRunner.LRC14.shapeOf v) →
      LonelyRunner.LRC14.clusterSize (LonelyRunner.LRC14.shapeOf v) ≤ 13 →
      0 < LonelyRunner.LRC14.witnessG2 (LonelyRunner.LRC14.shapeOf v) :=
  LonelyRunner.LRC14.Bonferroni.large_witness_pos_from_p0_wide_bound_shapes
    nuShape measGP p0Shape cap delta hδpos hbonf hDp0 hp0cap hmeasGP
#print axioms lrc_bonferroni_large_witness_pos_from_p0_wide_bound_shapes_audit

theorem lrc_bonferroni_lrc14_from_p0_positive_wide_bound_split_nodes_audit
    (nuShape measGP p0Shape cap delta : LonelyRunner.LRC14.Shape → ℝ)
    (hsmall : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      LonelyRunner.LRC14.clusterSize (LonelyRunner.LRC14.shapeOf v) ≤ 7 →
      (LonelyRunner.LRC14.witnessMP : ℝ) ≤
        LonelyRunner.LRC14.witnessG2 (LonelyRunner.LRC14.shapeOf v))
    (hδpos : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 → 0 < delta s)
    (hbonf : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s)
    (hDp0  : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      cap s ≤ measGP s)
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      LonelyRunner.LRC14.clusterSize (LonelyRunner.LRC14.shapeOf v) ≤ 13)
    (hpartA : ∀ v : Fin 13 → ℤ,
      0 < LonelyRunner.LRC14.witnessG2 (LonelyRunner.LRC14.shapeOf v) →
      (1 : ℝ) / 14 ≤ LonelyRunner.LRC14.Mreach v) :
    LonelyRunner.LRC14.LRC14Statement :=
  LonelyRunner.LRC14.Bonferroni.lrc14_from_p0_positive_wide_bound_split_nodes
    nuShape measGP p0Shape cap delta hsmall hδpos hbonf hDp0 hp0cap hmeasGP
    hsize hpartA
#print axioms lrc_bonferroni_lrc14_from_p0_positive_wide_bound_split_nodes_audit

theorem lrc_witness_attainment_exists_lonely_audit
    (v : Fin 13 → ℤ)
    (h : ∃ t : ℝ,
      (1 : ℝ) / 14 ≤ TournamentH7.LRCWitness.margin v t) :
    ∃ t : ℝ, LonelyRunner.Lonely 14 v t :=
  TournamentH7.LRCWitness.exists_lonely_of_margin_ge v 14 h
#print axioms lrc_witness_attainment_exists_lonely_audit

/-! ### LRC14 Bonferroni witness-route bridge -/

theorem lrc_bonferroni_floor_ge_mP_audit :
    ∀ k ∈ [8, 9, 10, 11, 12, 13],
      LonelyRunner.LRC14.witnessMP ≤
        LonelyRunner.LRC14.Bonferroni.nuConsec k +
          LonelyRunner.LRC14.capRat k - 1 :=
  LonelyRunner.LRC14.Bonferroni.bonferroni_floor_ge_mP
#print axioms lrc_bonferroni_floor_ge_mP_audit

theorem lrc_bonferroni_witnessG2_pos_from_p0_wide_bound_shapes_audit
    (nuShape measGP p0Shape capShape deltaShape : LonelyRunner.LRC14.Shape → ℝ)
    (hδ : ∀ s, 0 < deltaShape s)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s)
    (hDp0 : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ capShape s - deltaShape s)
    (hmeasGP : ∀ s, capShape s ≤ measGP s)
    (s : LonelyRunner.LRC14.Shape) :
    0 < LonelyRunner.LRC14.witnessG2 s :=
  LonelyRunner.LRC14.Bonferroni.witnessG2_pos_from_p0_wide_bound_shapes
    nuShape measGP p0Shape capShape deltaShape hδ hbonf hDp0 hp0cap hmeasGP s
#print axioms lrc_bonferroni_witnessG2_pos_from_p0_wide_bound_shapes_audit

theorem lrc14_from_p0_wide_bound_given_partA_audit
    (nuShape measGP p0Shape capShape deltaShape : LonelyRunner.LRC14.Shape → ℝ)
    (hδ : ∀ s, 0 < deltaShape s)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s)
    (hDp0 : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ capShape s - deltaShape s)
    (hmeasGP : ∀ s, capShape s ≤ measGP s)
    (hpartA : ∀ v : Fin 13 → ℤ,
      0 < LonelyRunner.LRC14.witnessG2 (LonelyRunner.LRC14.shapeOf v) →
      (1 : ℝ) / 14 ≤ LonelyRunner.LRC14.Mreach v) :
    LonelyRunner.LRC14.LRC14Statement :=
  LonelyRunner.LRC14.Bonferroni.lrc14_from_p0_wide_bound_given_partA
    nuShape measGP p0Shape capShape deltaShape hδ hbonf hDp0 hp0cap hmeasGP hpartA
#print axioms lrc14_from_p0_wide_bound_given_partA_audit

theorem lrc_bonferroni_witness_floor_all_clusters_audit
    (nuShape measGP : LonelyRunner.LRC14.Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s)
    (hnu1 : ∀ s, LonelyRunner.LRC14.clusterSize s ≤ 7 → nuShape s = 1)
    (hA : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      (LonelyRunner.LRC14.Bonferroni.nuConsec
        (LonelyRunner.LRC14.clusterSize s) : ℝ) ≤ nuShape s)
    (hBsmall : ∀ s, LonelyRunner.LRC14.clusterSize s ≤ 7 →
      (LonelyRunner.LRC14.capRat (LonelyRunner.LRC14.clusterSize s) : ℝ) ≤ measGP s)
    (hBlarge : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      (LonelyRunner.LRC14.capRat (LonelyRunner.LRC14.clusterSize s) : ℝ) ≤ measGP s)
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      LonelyRunner.LRC14.clusterSize (LonelyRunner.LRC14.shapeOf v) ≤ 13) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      (LonelyRunner.LRC14.witnessMP : ℝ) ≤
        LonelyRunner.LRC14.witnessG2 (LonelyRunner.LRC14.shapeOf v) :=
  LonelyRunner.LRC14.Bonferroni.witness_floor_from_bonferroni_nodes_all_clusters
    nuShape measGP hbonf hnu1 hA hBsmall hBlarge hsize
#print axioms lrc_bonferroni_witness_floor_all_clusters_audit

theorem lrc14_from_bonferroni_nodes_given_partA_audit
    (nuShape measGP : LonelyRunner.LRC14.Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s)
    (hnu1 : ∀ s, LonelyRunner.LRC14.clusterSize s ≤ 7 → nuShape s = 1)
    (hA : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      (LonelyRunner.LRC14.Bonferroni.nuConsec
        (LonelyRunner.LRC14.clusterSize s) : ℝ) ≤ nuShape s)
    (hBsmall : ∀ s, LonelyRunner.LRC14.clusterSize s ≤ 7 →
      (LonelyRunner.LRC14.capRat (LonelyRunner.LRC14.clusterSize s) : ℝ) ≤ measGP s)
    (hBlarge : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      (LonelyRunner.LRC14.capRat (LonelyRunner.LRC14.clusterSize s) : ℝ) ≤ measGP s)
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      LonelyRunner.LRC14.clusterSize (LonelyRunner.LRC14.shapeOf v) ≤ 13)
    (hpartA : ∀ v : Fin 13 → ℤ,
      0 < LonelyRunner.LRC14.witnessG2 (LonelyRunner.LRC14.shapeOf v) →
      (1 : ℝ) / 14 ≤ LonelyRunner.LRC14.Mreach v) :
    LonelyRunner.LRC14.LRC14Statement :=
  LonelyRunner.LRC14.Bonferroni.lrc14_from_bonferroni_nodes_given_partA
    nuShape measGP hbonf hnu1 hA hBsmall hBlarge hsize hpartA
#print axioms lrc14_from_bonferroni_nodes_given_partA_audit

/-! ### LRC14 Part-A finite-ruler arc-bound budgets -/

theorem lrc_partA_arc_div_lt_delta_of_le_bound_audit
    (arcCount arcBound vmax : ℕ) (deltaFloor delta : ℝ)
    (hvmax : 0 < vmax)
    (harc : arcCount ≤ arcBound)
    (hdelta : deltaFloor ≤ delta)
    (hbudget : (arcBound : ℝ) < deltaFloor * (vmax : ℝ)) :
    (arcCount : ℝ) / (vmax : ℝ) < delta :=
  LonelyRunner.LRC14.PartA.arc_div_lt_delta_of_le_bound
    arcCount arcBound vmax deltaFloor delta hvmax harc hdelta hbudget
#print axioms lrc_partA_arc_div_lt_delta_of_le_bound_audit

theorem lrc_partA_finite_witness_pos_from_uniform_arc_bound_shapes_audit
    (finiteRho delta : LonelyRunner.LRC14.Shape → ℝ)
    (arcCount vmax : LonelyRunner.LRC14.Shape → ℕ)
    (arcBound : ℕ) (deltaFloor : ℝ)
    (hfloor : ∀ s, delta s ≤ LonelyRunner.LRC14.witnessG2 s)
    (hdelta : ∀ s, deltaFloor ≤ delta s)
    (hvmax : ∀ s, 0 < vmax s)
    (harc : ∀ s, arcCount s ≤ arcBound)
    (herror : ∀ s, |finiteRho s - LonelyRunner.LRC14.witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcBound : ℝ) < deltaFloor * (vmax s : ℝ))
    (s : LonelyRunner.LRC14.Shape) :
    0 < finiteRho s :=
  LonelyRunner.LRC14.PartA.finite_witness_pos_from_uniform_arc_bound_shapes
    finiteRho delta arcCount vmax arcBound deltaFloor
    hfloor hdelta hvmax harc herror hbudget s
#print axioms lrc_partA_finite_witness_pos_from_uniform_arc_bound_shapes_audit

theorem lrc_partA_finite_witness_pos_from_p0_margin_uniform_arc_bound_shapes_audit
    (nuShape measGP p0Shape cap delta finiteRho : LonelyRunner.LRC14.Shape → ℝ)
    (arcCount vmax : LonelyRunner.LRC14.Shape → ℕ)
    (arcBound : ℕ) (deltaFloor : ℝ)
    (hbonf :
      ∀ s, nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, cap s ≤ measGP s)
    (hdelta : ∀ s, deltaFloor ≤ delta s)
    (hvmax : ∀ s, 0 < vmax s)
    (harc : ∀ s, arcCount s ≤ arcBound)
    (herror : ∀ s, |finiteRho s - LonelyRunner.LRC14.witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcBound : ℝ) < deltaFloor * (vmax s : ℝ))
    (s : LonelyRunner.LRC14.Shape) :
    0 < finiteRho s :=
  LonelyRunner.LRC14.PartA.finite_witness_pos_from_p0_margin_uniform_arc_bound_shapes
    nuShape measGP p0Shape cap delta finiteRho arcCount vmax arcBound deltaFloor
    hbonf hDp0 hp0cap hmeasGP hdelta hvmax harc herror hbudget s
#print axioms lrc_partA_finite_witness_pos_from_p0_margin_uniform_arc_bound_shapes_audit

theorem lrc_partA_finite_witness_pos_from_goodSet_margin_shapes_audit
    (Eof Pof : LonelyRunner.LRC14.Shape → List ℤ)
    (cap delta finiteRho : LonelyRunner.LRC14.Shape → ℝ)
    (arcCount vmax : LonelyRunner.LRC14.Shape → ℕ)
    (hanchor : ∀ s, (0 : ℤ) ∈ Eof s)
    (hwitness : ∀ s, LonelyRunner.LRC14.witnessG2 s =
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet (Eof s) ∩
          LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (hwide : ∀ s,
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet (Eof s))).toReal ≤
          cap s - delta s)
    (hdual : ∀ s, cap s ≤
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (herror : ∀ s, |finiteRho s - LonelyRunner.LRC14.witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcCount s : ℝ) / (vmax s : ℝ) < delta s)
    (s : LonelyRunner.LRC14.Shape) :
    0 < finiteRho s :=
  LonelyRunner.LRC14.PartA.finite_witness_pos_from_goodSet_margin_shapes
    Eof Pof cap delta finiteRho arcCount vmax hanchor hwitness
    hwide hdual herror hbudget s
#print axioms lrc_partA_finite_witness_pos_from_goodSet_margin_shapes_audit

theorem lrc_partA_finite_witness_pos_from_goodSet_margin_uniform_arc_bound_shapes_audit
    (Eof Pof : LonelyRunner.LRC14.Shape → List ℤ)
    (cap delta finiteRho : LonelyRunner.LRC14.Shape → ℝ)
    (arcCount vmax : LonelyRunner.LRC14.Shape → ℕ)
    (arcBound : ℕ) (deltaFloor : ℝ)
    (hanchor : ∀ s, (0 : ℤ) ∈ Eof s)
    (hwitness : ∀ s, LonelyRunner.LRC14.witnessG2 s =
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet (Eof s) ∩
          LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (hwide : ∀ s,
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet (Eof s))).toReal ≤
          cap s - delta s)
    (hdual : ∀ s, cap s ≤
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (hdelta : ∀ s, deltaFloor ≤ delta s)
    (hvmax : ∀ s, 0 < vmax s)
    (harc : ∀ s, arcCount s ≤ arcBound)
    (herror : ∀ s, |finiteRho s - LonelyRunner.LRC14.witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcBound : ℝ) < deltaFloor * (vmax s : ℝ))
    (s : LonelyRunner.LRC14.Shape) :
    0 < finiteRho s :=
  LonelyRunner.LRC14.PartA.finite_witness_pos_from_goodSet_margin_uniform_arc_bound_shapes
    Eof Pof cap delta finiteRho arcCount vmax arcBound deltaFloor
    hanchor hwitness hwide hdual hdelta hvmax harc herror hbudget s
#print axioms lrc_partA_finite_witness_pos_from_goodSet_margin_uniform_arc_bound_shapes_audit

theorem lrc_partA_finite_witness_pos_from_goodSet_p0_margin_shapes_audit
    (Eof Pof : LonelyRunner.LRC14.Shape → List ℤ)
    (cap delta finiteRho : LonelyRunner.LRC14.Shape → ℝ)
    (arcCount vmax : LonelyRunner.LRC14.Shape → ℕ)
    (hanchor : ∀ s, (0 : ℤ) ∈ Eof s)
    (hwitness : ∀ s, LonelyRunner.LRC14.witnessG2 s =
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet (Eof s) ∩
          LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (hp0cap : ∀ s, LonelyRunner.DenseCovers.p0 (Eof s) ≤ cap s - delta s)
    (hdual : ∀ s, cap s ≤
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (herror : ∀ s, |finiteRho s - LonelyRunner.LRC14.witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcCount s : ℝ) / (vmax s : ℝ) < delta s)
    (s : LonelyRunner.LRC14.Shape) :
    0 < finiteRho s :=
  LonelyRunner.LRC14.PartA.finite_witness_pos_from_goodSet_p0_margin_shapes
    Eof Pof cap delta finiteRho arcCount vmax hanchor hwitness
    hp0cap hdual herror hbudget s
#print axioms lrc_partA_finite_witness_pos_from_goodSet_p0_margin_shapes_audit

theorem lrc_partA_finite_witness_pos_from_goodSet_p0_margin_uniform_arc_bound_shapes_audit
    (Eof Pof : LonelyRunner.LRC14.Shape → List ℤ)
    (cap delta finiteRho : LonelyRunner.LRC14.Shape → ℝ)
    (arcCount vmax : LonelyRunner.LRC14.Shape → ℕ)
    (arcBound : ℕ) (deltaFloor : ℝ)
    (hanchor : ∀ s, (0 : ℤ) ∈ Eof s)
    (hwitness : ∀ s, LonelyRunner.LRC14.witnessG2 s =
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet (Eof s) ∩
          LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (hp0cap : ∀ s, LonelyRunner.DenseCovers.p0 (Eof s) ≤ cap s - delta s)
    (hdual : ∀ s, cap s ≤
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (hdelta : ∀ s, deltaFloor ≤ delta s)
    (hvmax : ∀ s, 0 < vmax s)
    (harc : ∀ s, arcCount s ≤ arcBound)
    (herror : ∀ s, |finiteRho s - LonelyRunner.LRC14.witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcBound : ℝ) < deltaFloor * (vmax s : ℝ))
    (s : LonelyRunner.LRC14.Shape) :
    0 < finiteRho s :=
  LonelyRunner.LRC14.PartA.finite_witness_pos_from_goodSet_p0_margin_uniform_arc_bound_shapes
    Eof Pof cap delta finiteRho arcCount vmax arcBound deltaFloor
    hanchor hwitness hp0cap hdual hdelta hvmax harc herror hbudget s
#print axioms lrc_partA_finite_witness_pos_from_goodSet_p0_margin_uniform_arc_bound_shapes_audit

theorem lrc_partA_lrc14_from_finite_partA_goodSet_margin_shapes_audit
    (Eof Pof : LonelyRunner.LRC14.Shape → List ℤ)
    (cap delta finiteRho : LonelyRunner.LRC14.Shape → ℝ)
    (arcCount vmax : LonelyRunner.LRC14.Shape → ℕ)
    (hanchor : ∀ s, (0 : ℤ) ∈ Eof s)
    (hwitness : ∀ s, LonelyRunner.LRC14.witnessG2 s =
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet (Eof s) ∩
          LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (hwide : ∀ s,
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.coverSet (Eof s))).toReal ≤
          cap s - delta s)
    (hdual : ∀ s, cap s ≤
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (herror : ∀ s, |finiteRho s - LonelyRunner.LRC14.witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcCount s : ℝ) / (vmax s : ℝ) < delta s)
    (hfinitePartA : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      0 < finiteRho (LonelyRunner.LRC14.shapeOf v) →
      (1 : ℝ) / 14 ≤ LonelyRunner.LRC14.Mreach v) :
    LonelyRunner.LRC14.LRC14Statement :=
  LonelyRunner.LRC14.PartA.lrc14_from_finite_partA_goodSet_margin_shapes
    Eof Pof cap delta finiteRho arcCount vmax hanchor hwitness
    hwide hdual herror hbudget hfinitePartA
#print axioms lrc_partA_lrc14_from_finite_partA_goodSet_margin_shapes_audit

theorem lrc_partA_lrc14_from_finite_partA_goodSet_p0_margin_shapes_audit
    (Eof Pof : LonelyRunner.LRC14.Shape → List ℤ)
    (cap delta finiteRho : LonelyRunner.LRC14.Shape → ℝ)
    (arcCount vmax : LonelyRunner.LRC14.Shape → ℕ)
    (hanchor : ∀ s, (0 : ℤ) ∈ Eof s)
    (hwitness : ∀ s, LonelyRunner.LRC14.witnessG2 s =
      (LonelyRunner.DenseCovers.slowμ
        (TournamentH7.GoodSet.goodSet (Eof s) ∩
          LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (hp0cap : ∀ s, LonelyRunner.DenseCovers.p0 (Eof s) ≤ cap s - delta s)
    (hdual : ∀ s, cap s ≤
      (LonelyRunner.DenseCovers.slowμ
        (LonelyRunner.DenseCovers.safeSet (Pof s))).toReal)
    (herror : ∀ s, |finiteRho s - LonelyRunner.LRC14.witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcCount s : ℝ) / (vmax s : ℝ) < delta s)
    (hfinitePartA : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      0 < finiteRho (LonelyRunner.LRC14.shapeOf v) →
      (1 : ℝ) / 14 ≤ LonelyRunner.LRC14.Mreach v) :
    LonelyRunner.LRC14.LRC14Statement :=
  LonelyRunner.LRC14.PartA.lrc14_from_finite_partA_goodSet_p0_margin_shapes
    Eof Pof cap delta finiteRho arcCount vmax hanchor hwitness
    hp0cap hdual herror hbudget hfinitePartA
#print axioms lrc_partA_lrc14_from_finite_partA_goodSet_p0_margin_shapes_audit

theorem lrc_partA_lrc14_from_finite_partA_p0_margin_split_shapes_audit
    (nuShape measGP p0Shape cap delta finiteRho : LonelyRunner.LRC14.Shape → ℝ)
    (arcCount vmax : LonelyRunner.LRC14.Shape → ℕ)
    (hsmall : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      LonelyRunner.LRC14.clusterSize (LonelyRunner.LRC14.shapeOf v) ≤ 7 →
      (LonelyRunner.LRC14.witnessMP : ℝ) ≤
        LonelyRunner.LRC14.witnessG2 (LonelyRunner.LRC14.shapeOf v))
    (hbonf : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s)
    (hDp0  : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      cap s ≤ measGP s)
    (herror : ∀ s, |finiteRho s - LonelyRunner.LRC14.witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hsmallBudget : ∀ s, LonelyRunner.LRC14.clusterSize s ≤ 7 →
      (arcCount s : ℝ) / (vmax s : ℝ) <
        (LonelyRunner.LRC14.witnessMP : ℝ))
    (hlargeBudget : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      (arcCount s : ℝ) / (vmax s : ℝ) < delta s)
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      LonelyRunner.LRC14.clusterSize (LonelyRunner.LRC14.shapeOf v) ≤ 13)
    (hfinitePartA : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      0 < finiteRho (LonelyRunner.LRC14.shapeOf v) →
      (1 : ℝ) / 14 ≤ LonelyRunner.LRC14.Mreach v) :
    LonelyRunner.LRC14.LRC14Statement :=
  LonelyRunner.LRC14.PartA.lrc14_from_finite_partA_p0_margin_split_shapes
    nuShape measGP p0Shape cap delta finiteRho arcCount vmax hsmall hbonf
    hDp0 hp0cap hmeasGP herror hsmallBudget hlargeBudget hsize hfinitePartA
#print axioms lrc_partA_lrc14_from_finite_partA_p0_margin_split_shapes_audit

theorem lrc14_partA_p0_margin_shapes_mul_audit
    (nuShape measGP p0Shape cap delta finiteRho : LonelyRunner.LRC14.Shape → ℝ)
    (arcCount vmax : LonelyRunner.LRC14.Shape → ℕ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s)
    (hDp0 : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, cap s ≤ measGP s)
    (hvmax : ∀ s, 0 < vmax s)
    (herror : ∀ s, |finiteRho s - LonelyRunner.LRC14.witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcCount s : ℝ) < delta s * (vmax s : ℝ))
    (hfinitePartA : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      0 < finiteRho (LonelyRunner.LRC14.shapeOf v) →
      (1 : ℝ) / 14 ≤ LonelyRunner.LRC14.Mreach v) :
    LonelyRunner.LRC14.LRC14Statement :=
  LonelyRunner.LRC14.PartA.lrc14_from_finite_partA_p0_margin_shapes_mul
    nuShape measGP p0Shape cap delta finiteRho arcCount vmax
    hbonf hDp0 hp0cap hmeasGP hvmax herror hbudget hfinitePartA
#print axioms lrc14_partA_p0_margin_shapes_mul_audit

theorem lrc14_partA_finite_bonferroni_nodes_mul_audit
    (nuShape measGP finiteRho : LonelyRunner.LRC14.Shape → ℝ)
    (arcCount vmax : LonelyRunner.LRC14.Shape → ℕ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ LonelyRunner.LRC14.witnessG2 s)
    (hnu1 : ∀ s, LonelyRunner.LRC14.clusterSize s ≤ 7 → nuShape s = 1)
    (hA : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      (LonelyRunner.LRC14.Bonferroni.nuConsec
        (LonelyRunner.LRC14.clusterSize s) : ℝ) ≤ nuShape s)
    (hBsmall : ∀ s, LonelyRunner.LRC14.clusterSize s ≤ 7 →
      (LonelyRunner.LRC14.capRat (LonelyRunner.LRC14.clusterSize s) : ℝ) ≤ measGP s)
    (hBlarge : ∀ s, 8 ≤ LonelyRunner.LRC14.clusterSize s →
      LonelyRunner.LRC14.clusterSize s ≤ 13 →
      (LonelyRunner.LRC14.capRat (LonelyRunner.LRC14.clusterSize s) : ℝ) ≤ measGP s)
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      LonelyRunner.LRC14.clusterSize (LonelyRunner.LRC14.shapeOf v) ≤ 13)
    (hvmax : ∀ s, 0 < vmax s)
    (herror : ∀ s, |finiteRho s - LonelyRunner.LRC14.witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcCount s : ℝ) <
      (LonelyRunner.LRC14.witnessMP : ℝ) * (vmax s : ℝ))
    (hfinitePartA : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      0 < finiteRho (LonelyRunner.LRC14.shapeOf v) →
      (1 : ℝ) / 14 ≤ LonelyRunner.LRC14.Mreach v) :
    LonelyRunner.LRC14.LRC14Statement :=
  LonelyRunner.LRC14.PartA.lrc14_from_finite_partA_bonferroni_nodes_mul
    nuShape measGP finiteRho arcCount vmax hbonf hnu1 hA hBsmall hBlarge
    hsize hvmax herror hbudget hfinitePartA
#print axioms lrc14_partA_finite_bonferroni_nodes_mul_audit

/-! ### LRC14 THM-563 period-max certificate kernel -/

theorem lrc_periodmax_all_worst_rows_headroom_positive_audit :
    ∀ r : Fin 6,
      0 < LonelyRunner.PeriodmaxCertificate.headroomNum
        (LonelyRunner.PeriodmaxCertificate.worstRow r) :=
  LonelyRunner.PeriodmaxCertificate.all_worst_rows_headroom_positive
#print axioms lrc_periodmax_all_worst_rows_headroom_positive_audit

theorem lrc_periodmax_k9_global_worst_among_worst_rows_audit :
    LonelyRunner.PeriodmaxCertificate.ratioLE LonelyRunner.PeriodmaxCertificate.k8
        LonelyRunner.PeriodmaxCertificate.k9 = true ∧
      LonelyRunner.PeriodmaxCertificate.ratioLE LonelyRunner.PeriodmaxCertificate.k9
        LonelyRunner.PeriodmaxCertificate.k9 = true ∧
      LonelyRunner.PeriodmaxCertificate.ratioLE LonelyRunner.PeriodmaxCertificate.k10
        LonelyRunner.PeriodmaxCertificate.k9 = true ∧
      LonelyRunner.PeriodmaxCertificate.ratioLE LonelyRunner.PeriodmaxCertificate.k11
        LonelyRunner.PeriodmaxCertificate.k9 = true ∧
      LonelyRunner.PeriodmaxCertificate.ratioLE LonelyRunner.PeriodmaxCertificate.k12
        LonelyRunner.PeriodmaxCertificate.k9 = true ∧
      LonelyRunner.PeriodmaxCertificate.ratioLE LonelyRunner.PeriodmaxCertificate.k13
        LonelyRunner.PeriodmaxCertificate.k9 = true :=
  LonelyRunner.PeriodmaxCertificate.k9_is_global_worst_among_worst_rows
#print axioms lrc_periodmax_k9_global_worst_among_worst_rows_audit

theorem lrc_periodmax_count_totals_audit :
    LonelyRunner.PeriodmaxCertificate.totalBases = 12805 ∧
      LonelyRunner.PeriodmaxCertificate.totalTrivial = 3995 ∧
      LonelyRunner.PeriodmaxCertificate.totalScanned = 8810 ∧
      LonelyRunner.PeriodmaxCertificate.totalSkipped = 0 ∧
      LonelyRunner.PeriodmaxCertificate.totalPassed = 12805 ∧
      LonelyRunner.PeriodmaxCertificate.totalFailed = 0 :=
  LonelyRunner.PeriodmaxCertificate.count_totals
#print axioms lrc_periodmax_count_totals_audit

theorem lrc_periodmax_every_count_row_passes_audit :
    ∀ r : Fin 6,
      (LonelyRunner.PeriodmaxCertificate.countRow r).passed =
          (LonelyRunner.PeriodmaxCertificate.countRow r).bases ∧
        (LonelyRunner.PeriodmaxCertificate.countRow r).failed = 0 ∧
        (LonelyRunner.PeriodmaxCertificate.countRow r).skipped = 0 :=
  LonelyRunner.PeriodmaxCertificate.every_count_row_passes
#print axioms lrc_periodmax_every_count_row_passes_audit

theorem lrc_periodmax_k8_periodmax_is_two_audit :
    LonelyRunner.PeriodmaxCertificate.k8.pmNum = 2 ∧
      LonelyRunner.PeriodmaxCertificate.k8.pmDen = 1 :=
  LonelyRunner.PeriodmaxCertificate.k8_periodmax_is_two
#print axioms lrc_periodmax_k8_periodmax_is_two_audit

/-! ### LRC14 HYP-2829 gK8 single-far arithmetic kernel -/

theorem lrc_gk8_singlefar_arithmetic_audit :
    ((2633 : Nat) * 588 < 2243 * 735) ∧
    ((3259 : Nat) * 2002 < 9895 * 735) ∧
    ((37 : Nat) * 7 < 40 * 7) ∧
    ((2323 : Nat) * 588 < 2243 * 980) ∧
    ((2876 : Nat) * 2002 < 9895 * 735) ∧
    ((62267 : Nat) * 7 < 40 * 12936) ∧
    ((2323 : Nat) * 735 < 2633 * 980) ∧
    ((2876 : Nat) < 3259) ∧
    ((62267 : Nat) * 7 < 37 * 12936) :=
  ⟨LonelyRunner.Gk8SingleFar.bounded_k8_below_cap,
   LonelyRunner.Gk8SingleFar.bounded_k9_below_cap,
   LonelyRunner.Gk8SingleFar.bounded_k10_below_cap,
   LonelyRunner.Gk8SingleFar.singlefar_k8_below_cap,
   LonelyRunner.Gk8SingleFar.singlefar_k9_below_cap,
   LonelyRunner.Gk8SingleFar.singlefar_k10_below_cap,
   LonelyRunner.Gk8SingleFar.singlefar_k8_below_bounded,
   LonelyRunner.Gk8SingleFar.singlefar_k9_below_bounded,
   LonelyRunner.Gk8SingleFar.singlefar_k10_below_bounded⟩
#print axioms lrc_gk8_singlefar_arithmetic_audit

/-! ### LRC14 Boolean/type signed aggregate cut -/

theorem lrc_boolean_type_cut_optimalCoeff_sum_audit :
    LonelyRunner.BooleanTypeCut.coeffSum LonelyRunner.BooleanTypeCut.optimalCoeff = 231049 :=
  LonelyRunner.BooleanTypeCut.optimalCoeff_sum
#print axioms lrc_boolean_type_cut_optimalCoeff_sum_audit

theorem lrc_boolean_type_cut_compactCoeff_sum_audit :
    LonelyRunner.BooleanTypeCut.coeffSum LonelyRunner.BooleanTypeCut.compactCoeff = 80 :=
  LonelyRunner.BooleanTypeCut.compactCoeff_sum
#print axioms lrc_boolean_type_cut_compactCoeff_sum_audit

theorem lrc_boolean_type_cut_optimal_active_rows_equal_audit :
    ∀ r : Fin 3,
      LonelyRunner.BooleanTypeCut.optimalMarginDen *
          LonelyRunner.BooleanTypeCut.evalNum LonelyRunner.BooleanTypeCut.optimalCoeff
            (LonelyRunner.BooleanTypeCut.active r) =
        LonelyRunner.BooleanTypeCut.optimalMarginNum *
          (LonelyRunner.BooleanTypeCut.active r).den :=
  LonelyRunner.BooleanTypeCut.optimal_active_rows_equal
#print axioms lrc_boolean_type_cut_optimal_active_rows_equal_audit

theorem lrc_boolean_type_cut_compact_active_rows_above_audit :
    ∀ r : Fin 3,
      LonelyRunner.BooleanTypeCut.compactMarginDen *
          LonelyRunner.BooleanTypeCut.evalNum LonelyRunner.BooleanTypeCut.compactCoeff
            (LonelyRunner.BooleanTypeCut.active r) >=
        LonelyRunner.BooleanTypeCut.compactMarginNum *
          (LonelyRunner.BooleanTypeCut.active r).den :=
  LonelyRunner.BooleanTypeCut.compact_active_rows_above
#print axioms lrc_boolean_type_cut_compact_active_rows_above_audit

theorem lrc_boolean_type_cut_compact_active2_equal_audit :
    LonelyRunner.BooleanTypeCut.compactMarginDen *
        LonelyRunner.BooleanTypeCut.evalNum LonelyRunner.BooleanTypeCut.compactCoeff
          LonelyRunner.BooleanTypeCut.active2 =
      LonelyRunner.BooleanTypeCut.compactMarginNum *
        LonelyRunner.BooleanTypeCut.active2.den :=
  LonelyRunner.BooleanTypeCut.compact_active2_equal
#print axioms lrc_boolean_type_cut_compact_active2_equal_audit

theorem lrc_boolean_type_cut_margins_positive_audit :
    0 < LonelyRunner.BooleanTypeCut.optimalMarginNum ∧
      0 < LonelyRunner.BooleanTypeCut.optimalMarginDen ∧
      0 < LonelyRunner.BooleanTypeCut.compactMarginNum ∧
      0 < LonelyRunner.BooleanTypeCut.compactMarginDen :=
  LonelyRunner.BooleanTypeCut.margins_positive
#print axioms lrc_boolean_type_cut_margins_positive_audit

/-! ### LRC14 L7 discrepancy integer core -/

theorem lrc_l7_cell_apex_necessity_audit (c S : Int)
    (h : LonelyRunner.L7Discrepancy.iabs (7 * c - S) = 0) :
    (7 : Int) ∣ S :=
  LonelyRunner.L7Discrepancy.cell_apex_necessity c S h
#print axioms lrc_l7_cell_apex_necessity_audit

theorem lrc_l7_matrix_apex_necessity_audit
    (c : LonelyRunner.L7Discrepancy.Mat) (S : Int)
    (h : LonelyRunner.L7Discrepancy.Ddef c S = 0) :
    (7 : Int) ∣ S :=
  LonelyRunner.L7Discrepancy.matrix_apex_necessity c S h
#print axioms lrc_l7_matrix_apex_necessity_audit

theorem lrc_l7_c32_disc_audit :
    LonelyRunner.L7Discrepancy.Ddef LonelyRunner.L7Discrepancy.c32 6 = 252 :=
  LonelyRunner.L7Discrepancy.c32_disc
#print axioms lrc_l7_c32_disc_audit

theorem lrc_l7_c87_disc_zero_audit :
    LonelyRunner.L7Discrepancy.Ddef LonelyRunner.L7Discrepancy.c87 56 = 0 :=
  LonelyRunner.L7Discrepancy.c87_disc_zero
#print axioms lrc_l7_c87_disc_zero_audit

theorem lrc_l7_c21_disc_audit :
    LonelyRunner.L7Discrepancy.Ddef LonelyRunner.L7Discrepancy.c21 2 = 140 :=
  LonelyRunner.L7Discrepancy.c21_disc
#print axioms lrc_l7_c21_disc_audit

theorem lrc_l7_c32_rowSums_audit :
    ∀ i, LonelyRunner.L7Discrepancy.rowSum LonelyRunner.L7Discrepancy.c32 i = 6 :=
  LonelyRunner.L7Discrepancy.c32_rowSums
#print axioms lrc_l7_c32_rowSums_audit

theorem lrc_l7_c32_colSums_audit :
    ∀ j, LonelyRunner.L7Discrepancy.colSum LonelyRunner.L7Discrepancy.c32 j = 6 :=
  LonelyRunner.L7Discrepancy.c32_colSums
#print axioms lrc_l7_c32_colSums_audit

/-! ### LRC14 factorial atom / Q0 boundary quotient -/

theorem lrc_factorial_basis_moment_delta_audit :
    ∀ i j : Fin LonelyRunner.FactorialAtom.atomCount,
      LonelyRunner.FactorialAtom.moment (LonelyRunner.FactorialAtom.basis j) i =
        if i = j then 1 else 0 :=
  LonelyRunner.FactorialAtom.basis_moment_delta
#print axioms lrc_factorial_basis_moment_delta_audit

theorem lrc_factorial_originCoeff_delta_audit :
    ∀ t : Fin LonelyRunner.FactorialAtom.atomCount,
      LonelyRunner.FactorialAtom.originCoeff t = if t.val = 0 then 1 else 0 :=
  LonelyRunner.FactorialAtom.originCoeff_delta
#print axioms lrc_factorial_originCoeff_delta_audit

theorem lrc_factorial_basis_q0_sign_audit :
    ∀ j : Fin LonelyRunner.FactorialAtom.atomCount,
      LonelyRunner.FactorialAtom.q0 (LonelyRunner.FactorialAtom.basis j) =
        LonelyRunner.FactorialAtom.altSign j.val :=
  LonelyRunner.FactorialAtom.basis_q0_sign
#print axioms lrc_factorial_basis_q0_sign_audit

theorem lrc_factorial_U4_basis_audit :
    ∀ j : Fin LonelyRunner.FactorialAtom.atomCount,
      LonelyRunner.FactorialAtom.U4 (LonelyRunner.FactorialAtom.basis j) =
        if j.val <= 4 then LonelyRunner.FactorialAtom.altSign j.val else 0 :=
  LonelyRunner.FactorialAtom.U4_basis
#print axioms lrc_factorial_U4_basis_audit

theorem lrc_factorial_low12_basis_audit :
    ∀ j : Fin LonelyRunner.FactorialAtom.atomCount,
      LonelyRunner.FactorialAtom.low12 (LonelyRunner.FactorialAtom.basis j) =
        if j.val = 1 ∨ j.val = 2 then 1 else 0 :=
  LonelyRunner.FactorialAtom.low12_basis
#print axioms lrc_factorial_low12_basis_audit

theorem lrc_factorial_cheapScaled_q0_audit :
    ∀ r : Fin 5,
      LonelyRunner.FactorialAtom.q0 (LonelyRunner.FactorialAtom.cheapScaled r) =
        LonelyRunner.FactorialAtom.cheapScale r :=
  LonelyRunner.FactorialAtom.cheapScaled_q0
#print axioms lrc_factorial_cheapScaled_q0_audit

theorem lrc_factorial_cheapScaled_tail45_audit :
    ∀ r : Fin 5,
      LonelyRunner.FactorialAtom.tail45 (LonelyRunner.FactorialAtom.cheapScaled r) =
        LonelyRunner.FactorialAtom.cheapTailTarget r :=
  LonelyRunner.FactorialAtom.cheapScaled_tail45
#print axioms lrc_factorial_cheapScaled_tail45_audit

theorem lrc_factorial_cheapScaled_outside_tailStrip_bool_audit :
    ∀ r : Fin 5,
      LonelyRunner.FactorialAtom.outsideTailStripBool
          (LonelyRunner.FactorialAtom.cheapScaled r)
          (LonelyRunner.FactorialAtom.cheapScale r) = true :=
  LonelyRunner.FactorialAtom.cheapScaled_outside_tailStrip_bool
#print axioms lrc_factorial_cheapScaled_outside_tailStrip_bool_audit

theorem lrc_factorial_tailStrip_constants_order_audit :
    LonelyRunner.FactorialAtom.tailFloorNum * LonelyRunner.FactorialAtom.tailCeilDen <
      LonelyRunner.FactorialAtom.tailCeilNum * LonelyRunner.FactorialAtom.tailFloorDen :=
  LonelyRunner.FactorialAtom.tailStrip_constants_order
#print axioms lrc_factorial_tailStrip_constants_order_audit

theorem lrc_factorial_cheapScaled_tailStripSide_audit :
    ∀ r : Fin 5,
      LonelyRunner.FactorialAtom.tailStripSide
          (LonelyRunner.FactorialAtom.cheapScaled r)
          (LonelyRunner.FactorialAtom.cheapScale r) =
        LonelyRunner.FactorialAtom.cheapTailSideTarget r :=
  LonelyRunner.FactorialAtom.cheapScaled_tailStripSide
#print axioms lrc_factorial_cheapScaled_tailStripSide_audit

theorem lrc_factorial_gK8_values_audit :
    ∀ t : Fin LonelyRunner.FactorialAtom.atomCount,
      LonelyRunner.FactorialAtom.gK8 t =
        (if t.val = 0 then 10 else if t.val = 3 then 1
          else if t.val = 6 then 10 else 0) :=
  LonelyRunner.FactorialAtom.gK8_values
#print axioms lrc_factorial_gK8_values_audit

theorem lrc_factorial_gK8_dominates_audit :
    ∀ t : Fin LonelyRunner.FactorialAtom.atomCount,
      (if t.val = 0 then (10 : Int) else 0) <= LonelyRunner.FactorialAtom.gK8 t :=
  LonelyRunner.FactorialAtom.gK8_dominates
#print axioms lrc_factorial_gK8_dominates_audit

theorem lrc_factorial_LyK8_readout_audit (q : LonelyRunner.FactorialAtom.Atom) :
    LonelyRunner.FactorialAtom.LyK8 q =
      10 * q ⟨0, by decide⟩ + q ⟨3, by decide⟩ + 10 * q ⟨6, by decide⟩ :=
  LonelyRunner.FactorialAtom.LyK8_readout q
#print axioms lrc_factorial_LyK8_readout_audit

theorem lrc_factorial_LyK8_extremeMass_readout_audit
    (q : LonelyRunner.FactorialAtom.Atom) :
    LonelyRunner.FactorialAtom.LyK8 q =
      10 * LonelyRunner.FactorialAtom.extremeMass q + q ⟨3, by decide⟩ :=
  LonelyRunner.FactorialAtom.LyK8_extremeMass_readout q
#print axioms lrc_factorial_LyK8_extremeMass_readout_audit

theorem lrc_factorial_LyK8_moment_form_audit
    (q : LonelyRunner.FactorialAtom.Atom) :
    LonelyRunner.FactorialAtom.LyK8 q =
      10 * LonelyRunner.FactorialAtom.moment q ⟨0, by decide⟩
        - 10 * LonelyRunner.FactorialAtom.moment q ⟨1, by decide⟩
        + 10 * LonelyRunner.FactorialAtom.moment q ⟨2, by decide⟩
        - 9 * LonelyRunner.FactorialAtom.moment q ⟨3, by decide⟩
        + 6 * LonelyRunner.FactorialAtom.moment q ⟨4, by decide⟩ :=
  LonelyRunner.FactorialAtom.LyK8_moment_form q
#print axioms lrc_factorial_LyK8_moment_form_audit

theorem lrc_factorial_LyK8_probability_moment_form_audit
    (q : LonelyRunner.FactorialAtom.Atom)
    (h0 : LonelyRunner.FactorialAtom.moment q ⟨0, by decide⟩ = 1) :
    LonelyRunner.FactorialAtom.LyK8 q =
      10
        - 10 * LonelyRunner.FactorialAtom.moment q ⟨1, by decide⟩
        + 10 * LonelyRunner.FactorialAtom.moment q ⟨2, by decide⟩
        - 9 * LonelyRunner.FactorialAtom.moment q ⟨3, by decide⟩
        + 6 * LonelyRunner.FactorialAtom.moment q ⟨4, by decide⟩ :=
  LonelyRunner.FactorialAtom.LyK8_probability_moment_form q h0
#print axioms lrc_factorial_LyK8_probability_moment_form_audit

theorem lrc_factorial_LyK8_moment_extremeMass_identity_audit
    (q : LonelyRunner.FactorialAtom.Atom) :
    10 * LonelyRunner.FactorialAtom.moment q ⟨0, by decide⟩
        - 10 * LonelyRunner.FactorialAtom.moment q ⟨1, by decide⟩
        + 10 * LonelyRunner.FactorialAtom.moment q ⟨2, by decide⟩
        - 9 * LonelyRunner.FactorialAtom.moment q ⟨3, by decide⟩
        + 6 * LonelyRunner.FactorialAtom.moment q ⟨4, by decide⟩ =
      10 * LonelyRunner.FactorialAtom.extremeMass q + q ⟨3, by decide⟩ :=
  LonelyRunner.FactorialAtom.LyK8_moment_extremeMass_identity q
#print axioms lrc_factorial_LyK8_moment_extremeMass_identity_audit

theorem lrc_factorial_delsarte_bound_k8_audit
    (q : LonelyRunner.FactorialAtom.Atom)
    (hq : ∀ t : Fin LonelyRunner.FactorialAtom.atomCount, 0 <= q t) :
    10 * LonelyRunner.FactorialAtom.q0 q <= LonelyRunner.FactorialAtom.LyK8 q :=
  LonelyRunner.FactorialAtom.delsarte_bound_k8 q hq
#print axioms lrc_factorial_delsarte_bound_k8_audit

theorem lrc_factorial_LyK9_readout_audit (q : LonelyRunner.FactorialAtom.Atom) :
    LonelyRunner.FactorialAtom.LyK9 q =
      18 * q ⟨0, by decide⟩ + 5 * q ⟨1, by decide⟩
        + 2 * q ⟨4, by decide⟩ + 3 * q ⟨5, by decide⟩ :=
  LonelyRunner.FactorialAtom.LyK9_readout q
#print axioms lrc_factorial_LyK9_readout_audit

theorem lrc_factorial_gK9_values_audit :
    ∀ t : Fin LonelyRunner.FactorialAtom.atomCount,
      LonelyRunner.FactorialAtom.gK9 t =
        (if t.val = 0 then 18 else if t.val = 1 then 5
          else if t.val = 4 then 2 else if t.val = 5 then 3 else 0) :=
  LonelyRunner.FactorialAtom.gK9_values
#print axioms lrc_factorial_gK9_values_audit

theorem lrc_factorial_gK9_dominates_audit :
    ∀ t : Fin LonelyRunner.FactorialAtom.atomCount,
      (if t.val = 0 then (18 : Int) else 0) <= LonelyRunner.FactorialAtom.gK9 t :=
  LonelyRunner.FactorialAtom.gK9_dominates
#print axioms lrc_factorial_gK9_dominates_audit

theorem lrc_factorial_delsarte_bound_k9_audit
    (q : LonelyRunner.FactorialAtom.Atom)
    (hq : ∀ t : Fin LonelyRunner.FactorialAtom.atomCount, 0 <= q t) :
    18 * LonelyRunner.FactorialAtom.q0 q <= LonelyRunner.FactorialAtom.LyK9 q :=
  LonelyRunner.FactorialAtom.delsarte_bound_k9 q hq
#print axioms lrc_factorial_delsarte_bound_k9_audit

theorem lrc_factorial_LyK11_readout_audit (q : LonelyRunner.FactorialAtom.Atom) :
    LonelyRunner.FactorialAtom.LyK11 q =
      6 * q ⟨0, by decide⟩ + 3 * q ⟨1, by decide⟩
        + q ⟨2, by decide⟩ + q ⟨5, by decide⟩ + 3 * q ⟨6, by decide⟩ :=
  LonelyRunner.FactorialAtom.LyK11_readout q
#print axioms lrc_factorial_LyK11_readout_audit

theorem lrc_factorial_gK11_values_audit :
    ∀ t : Fin LonelyRunner.FactorialAtom.atomCount,
      LonelyRunner.FactorialAtom.gK11 t =
        (if t.val = 0 then 6 else if t.val = 1 then 3
          else if t.val = 2 then 1 else if t.val = 5 then 1
          else if t.val = 6 then 3 else 0) :=
  LonelyRunner.FactorialAtom.gK11_values
#print axioms lrc_factorial_gK11_values_audit

theorem lrc_factorial_gK11_dominates_audit :
    ∀ t : Fin LonelyRunner.FactorialAtom.atomCount,
      (if t.val = 0 then (6 : Int) else 0) <= LonelyRunner.FactorialAtom.gK11 t :=
  LonelyRunner.FactorialAtom.gK11_dominates
#print axioms lrc_factorial_gK11_dominates_audit

theorem lrc_factorial_delsarte_bound_k11_audit
    (q : LonelyRunner.FactorialAtom.Atom)
    (hq : ∀ t : Fin LonelyRunner.FactorialAtom.atomCount, 0 <= q t) :
    6 * LonelyRunner.FactorialAtom.q0 q <= LonelyRunner.FactorialAtom.LyK11 q :=
  LonelyRunner.FactorialAtom.delsarte_bound_k11 q hq
#print axioms lrc_factorial_delsarte_bound_k11_audit

theorem lrc_factorial_gK8_all_binding_rows_audit :
    ((2633 : Int) * 588 <= 2243 * 735) ∧
      ((3259 : Int) * 2002 <= 9895 * 735) ∧
      ((37 : Int) * 91 <= 550 * 7) ∧
      ((26603 : Int) * 91 <= 660 * 4410) ∧
      ((29287 : Int) * 7 <= 60 * 4410) ∧
      ((61529 : Int) <= 10 * 8820) :=
  LonelyRunner.FactorialAtom.capClear_gK8_all_binding_rows
#print axioms lrc_factorial_gK8_all_binding_rows_audit

/-! ### LRC14 q6 contraction arithmetic kernel -/

theorem lrc_q6_consecK9_bound_exact_audit :
    LonelyRunner.Q6Contraction.consecK9BoundNum * 5 =
      3 * LonelyRunner.Q6Contraction.commonDen :=
  LonelyRunner.Q6Contraction.consecK9_bound_exact
#print axioms lrc_q6_consecK9_bound_exact_audit

theorem lrc_q6_consecK10_bound_exact_audit :
    LonelyRunner.Q6Contraction.consecK10BoundNum * 35 =
      23 * LonelyRunner.Q6Contraction.commonDen :=
  LonelyRunner.Q6Contraction.consecK10_bound_exact
#print axioms lrc_q6_consecK10_bound_exact_audit

theorem lrc_q6_consecK9_strict_contraction_audit : (3 : Nat) < 5 :=
  LonelyRunner.Q6Contraction.consecK9_strict_contraction
#print axioms lrc_q6_consecK9_strict_contraction_audit

theorem lrc_q6_consecK10_strict_contraction_audit : (23 : Nat) < 35 :=
  LonelyRunner.Q6Contraction.consecK10_strict_contraction
#print axioms lrc_q6_consecK10_strict_contraction_audit

theorem lrc_q6_generalScout_strict_contraction_audit : (33 : Nat) < 35 :=
  LonelyRunner.Q6Contraction.generalScout_strict_contraction
#print axioms lrc_q6_generalScout_strict_contraction_audit

/-! ### LRC14 death-chain/live-depth quotient -/

theorem lrc_cover_oneFar_live_iff_audit (t : Nat) (ht : t <= 6) :
    LonelyRunner.LiveDepth LonelyRunner.coverCurrency 1 t <-> t = 1 :=
  LonelyRunner.cover_oneFar_live_iff t ht
#print axioms lrc_cover_oneFar_live_iff_audit

theorem lrc_survival_oneFar_live_iff_audit (t : Nat) (ht : t <= 6) :
    LonelyRunner.LiveDepth LonelyRunner.survivalCurrency 1 t <->
      t = 1 \/ t = 5 \/ t = 6 :=
  LonelyRunner.survival_oneFar_live_iff t ht
#print axioms lrc_survival_oneFar_live_iff_audit

theorem lrc_survival_twoFar_live_iff_audit (t : Nat) (ht : t <= 6) :
    LonelyRunner.LiveDepth LonelyRunner.survivalCurrency 2 t <->
      t = 1 \/ t = 2 \/ t = 5 \/ t = 6 :=
  LonelyRunner.survival_twoFar_live_iff t ht
#print axioms lrc_survival_twoFar_live_iff_audit

theorem lrc_survival_threeFar_live_iff_audit (t : Nat) (ht : t <= 6) :
    LonelyRunner.LiveDepth LonelyRunner.survivalCurrency 3 t <->
      t = 1 \/ t = 2 \/ t = 3 \/ t = 5 \/ t = 6 :=
  LonelyRunner.survival_threeFar_live_iff t ht
#print axioms lrc_survival_threeFar_live_iff_audit

theorem lrc_survival_fourFar_live_iff_audit (t : Nat) (ht : t <= 6) :
    LonelyRunner.LiveDepth LonelyRunner.survivalCurrency 4 t <-> 1 <= t :=
  LonelyRunner.survival_fourFar_live_iff t ht
#print axioms lrc_survival_fourFar_live_iff_audit

theorem lrc_survival_twoFar_middle_silent_audit :
    LonelyRunner.SilentDepth LonelyRunner.survivalCurrency 2 3 /\
      LonelyRunner.SilentDepth LonelyRunner.survivalCurrency 2 4 :=
  LonelyRunner.survival_twoFar_middle_silent
#print axioms lrc_survival_twoFar_middle_silent_audit

/-! ### Product-sum defect normal form (THM-361 list core) -/

theorem product_sum_iff_core_audit (xs : List Nat) :
    ProductSum.IsProductSum xs ↔
      (ProductSum.core xs).prod = (ProductSum.core xs).sum + ProductSum.ones xs :=
  ProductSum.product_sum_iff_core xs
#print axioms product_sum_iff_core_audit

theorem pad_core_product_sum_audit {d : Nat} {c : List Nat}
    (h : c.prod = c.sum + d) :
    ProductSum.IsProductSum ((List.replicate d 1) ++ c) :=
  ProductSum.pad_core_product_sum h
#print axioms pad_core_product_sum_audit

theorem core_defect_eq_ones_of_product_sum_audit {xs : List Nat}
    (h : ProductSum.IsProductSum xs) :
    (ProductSum.core xs).prod - (ProductSum.core xs).sum = ProductSum.ones xs :=
  ProductSum.core_defect_eq_ones_of_product_sum h
#print axioms core_defect_eq_ones_of_product_sum_audit

theorem two_entry_product_sum_audit {a b : Nat} (ha : 0 < a) (hb : 0 < b)
    (h : ProductSum.IsProductSum [a, b]) :
    a = 2 ∧ b = 2 :=
  ProductSum.two_entry_product_sum ha hb h
#print axioms two_entry_product_sum_audit

theorem no_three_ge_two_product_sum_audit {a b c : Nat}
    (ha : 2 ≤ a) (hb : 2 ≤ b) (hc : 2 ≤ c)
    (h : ProductSum.IsProductSum [a, b, c]) : False :=
  ProductSum.no_three_ge_two_product_sum ha hb hc h
#print axioms no_three_ge_two_product_sum_audit

theorem one_cons_two_entry_product_sum_audit {a b : Nat}
    (ha : 0 < a) (hb : 0 < b) (ha1 : a ≠ 1) (hb1 : b ≠ 1)
    (h : ProductSum.IsProductSum [1, a, b]) :
    (a = 2 ∧ b = 3) ∨ (a = 3 ∧ b = 2) :=
  ProductSum.one_cons_two_entry_product_sum ha hb ha1 hb1 h
#print axioms one_cons_two_entry_product_sum_audit

theorem three_entry_distinct_product_sum_audit {a b c : Nat}
    (ha : 0 < a) (hb : 0 < b) (hc : 0 < c)
    (hab : a ≠ b) (hac : a ≠ c) (hbc : b ≠ c)
    (h : ProductSum.IsProductSum [a, b, c]) :
    List.Perm [a, b, c] [1, 2, 3] :=
  ProductSum.three_entry_distinct_product_sum ha hb hc hab hac hbc h
#print axioms three_entry_distinct_product_sum_audit

/-! ### Type-A root-sign atoms -/

theorem typeA_root_self_audit {n : ℕ} (i : Fin n) :
    TypeA.root i i = 0 :=
  TypeA.root_self i
#print axioms typeA_root_self_audit

theorem typeA_root_swap_audit {n : ℕ} (i j : Fin n) :
    TypeA.root j i = -TypeA.root i j :=
  TypeA.root_swap i j
#print axioms typeA_root_swap_audit

theorem typeA_root_add_root_audit {n : ℕ} (i j k : Fin n) :
    TypeA.root i j + TypeA.root j k = TypeA.root i k :=
  TypeA.root_add_root i j k
#print axioms typeA_root_add_root_audit

theorem typeA_root_eq_zero_iff_audit {n : ℕ} (i j : Fin n) :
    TypeA.root i j = 0 ↔ i = j :=
  TypeA.root_eq_zero_iff i j
#print axioms typeA_root_eq_zero_iff_audit

theorem typeA_root_cycle_sum_audit {n : ℕ} (i j k : Fin n) :
    TypeA.root i j + TypeA.root j k + TypeA.root k i = 0 :=
  TypeA.root_cycle_sum i j k
#print axioms typeA_root_cycle_sum_audit

theorem typeA_walkRootSum_append_single_audit {n : ℕ}
    (a b : Fin n) (middle : List (Fin n)) :
    TypeA.walkRootSum (a :: middle ++ [b]) = TypeA.root a b :=
  TypeA.walkRootSum_append_single a b middle
#print axioms typeA_walkRootSum_append_single_audit

theorem typeA_walkRootSum_closed_audit {n : ℕ}
    (a : Fin n) (middle : List (Fin n)) :
    TypeA.walkRootSum (a :: middle ++ [a]) = 0 :=
  TypeA.walkRootSum_closed a middle
#print axioms typeA_walkRootSum_closed_audit

theorem typeA_rootWalk_rootTotal_audit {n : ℕ} (W : TypeA.RootWalk n) :
    W.rootTotal = TypeA.root W.source W.target :=
  TypeA.RootWalk.rootTotal_eq_boundary W
#print axioms typeA_rootWalk_rootTotal_audit

theorem typeA_rootPacket_rootTotal_audit {n : ℕ} (P : TypeA.RootPacket n) :
    P.rootTotal = 0 :=
  TypeA.RootPacket.rootTotal_eq_zero P
#print axioms typeA_rootPacket_rootTotal_audit

theorem directedCycle_toRootPacket_rootTotal_audit {n k : ℕ} {T : Tournament n}
    (C : Tournament.DirectedCycle T k) :
    C.toRootPacket.rootTotal = 0 :=
  Tournament.DirectedCycle.toRootPacket_rootTotal C
#print axioms directedCycle_toRootPacket_rootTotal_audit

/-! ### Natural-number operation shadows -/

theorem natOperation_addShadow_iff_lt_audit (x z : ℕ) :
    NatOperation.AddShadow x z ↔ x < z :=
  NatOperation.addShadow_iff_lt x z
#print axioms natOperation_addShadow_iff_lt_audit

theorem natOperation_mulUnitShadow_iff_dvd_audit {x z : ℕ} (hz : 1 ≤ z) :
    NatOperation.MulUnitShadow x z ↔ x ∣ z :=
  NatOperation.mulUnitShadow_iff_dvd hz
#print axioms natOperation_mulUnitShadow_iff_dvd_audit

theorem natOperation_mulShadow_iff_dvd_and_lt_audit {x z : ℕ} (hx : 1 ≤ x) :
    NatOperation.MulShadow x z ↔ x ∣ z ∧ x < z :=
  NatOperation.mulShadow_iff_dvd_and_lt hx
#print axioms natOperation_mulShadow_iff_dvd_and_lt_audit

theorem natOperation_shifted_binary_collision_iff_audit (a b : ℕ) :
    (a + 1) + (b + 1) = (a + 1) * (b + 1) ↔ a * b = 1 :=
  NatOperation.shifted_binary_collision_iff a b
#print axioms natOperation_shifted_binary_collision_iff_audit

theorem natOperation_twoFactor_productSum_iff_audit (r a b : ℕ) :
    r + (a + 1) + (b + 1) = (a + 1) * (b + 1) ↔ a * b = r + 1 :=
  NatOperation.twoFactor_productSum_iff r a b
#print axioms natOperation_twoFactor_productSum_iff_audit

theorem natOperation_trivial_twoFactor_productSum_audit (r : ℕ) :
    r + (1 + 1) + ((r + 1) + 1) = (1 + 1) * ((r + 1) + 1) :=
  NatOperation.trivial_twoFactor_productSum r
#print axioms natOperation_trivial_twoFactor_productSum_audit

/-! ### Audit each individual theorem -/

theorem H_ne_seven_audit {n : ℕ} (T : Tournament n) : H T ≠ 7 := H_ne_seven T
#print axioms H_ne_seven_audit

theorem H_ne_twentyone_audit {n : ℕ} (T : Tournament n) : H T ≠ 21 := H_ne_twentyone T
#print axioms H_ne_twentyone_audit

theorem H_ne_sixtythree_le_seven_audit {n : ℕ} (hn : n ≤ 7) (T : Tournament n) :
    H T ≠ 63 := H_ne_sixtythree_le_seven hn T
#print axioms H_ne_sixtythree_le_seven_audit

theorem H_pos_audit {n : ℕ} (hn : 1 ≤ n) (T : Tournament n) : H T ≠ 0 := H_pos hn T
#print axioms H_pos_audit

theorem forbidden_pair_audit {n : ℕ} (T : Tournament n) :
    H T ≠ 7 ∧ H T ≠ 21 := H_not_in_forbidden_pair T
#print axioms forbidden_pair_audit

theorem forbidden_trio_le_seven_audit {n : ℕ} (hn : n ≤ 7) (T : Tournament n) :
    H T ≠ 7 ∧ H T ≠ 21 ∧ H T ≠ 63 := H_not_in_forbidden_trio_le_seven hn T
#print axioms forbidden_trio_le_seven_audit

/-! ### Project-novel results — audit -/

-- (REMOVED) gridSym_iff_audit was based on the wrong tilde_eq_reversed_op
-- axiom; tile-complement and grid-reflection are different tiling
-- involutions; correct THM-280 formalisation requires a concrete
-- tile-coordinate model and is deferred.

/-- Score-formula corollary: regular tournaments are not self-flip
    (project-novel, oracle-2026-05-11-S1).  Used to prove
    Paley(p) ∉ SF for p ≡ 3 (mod 4). -/
theorem regular_not_SF_audit {n : ℕ} (T : Tournament n)
    (hbp : HasBasePath T) (hn : 3 ≤ n) (hreg : IsRegular T)
    (v0 : Fin n) (hv0 : v0.val = 0) :
    (tilde T).outDegree v0 ≠ T.outDegree v0 :=
  regular_not_SF T hbp hn hreg v0 hv0
#print axioms regular_not_SF_audit

/-! ### THM-330 (SC Cut Theorem) - project-novel, opus-2026-05-27-S1 -/

/-- THM-330: SC iff every cut k ∈ {1, …, n-1} has a crossing-upward arc. -/
theorem thm330_audit {n : ℕ} (T : Tournament n) (hbp : HasBasePath T) :
    IsStronglyConnected T ↔ ∀ k, 1 ≤ k → k < n → CrossesUpward T k :=
  thm330_SC_iff_all_cuts_crossing T hbp
#print axioms thm330_audit

/-- THM-330 EASY direction (now FULLY PROVED in Lean — no axioms beyond foundations). -/
theorem thm330_easy_audit {n : ℕ} (T : Tournament n) (hbp : HasBasePath T)
    (h : ∀ k, 1 ≤ k → k < n → CrossesUpward T k) : IsStronglyConnected T :=
  crossesUpward_all_implies_SC T hbp h
#print axioms thm330_easy_audit

/-- Base-path descent: any vertex u reaches any v with v.val ≤ u.val. PROVED. -/
theorem reaches_descent_audit {n : ℕ}
    (T : Tournament n) (hbp : HasBasePath T) (u v : Fin n) (h : v.val ≤ u.val) :
    Reaches T u v :=
  reaches_descent T hbp u v h
#print axioms reaches_descent_audit

/-- Every vertex reaches 0 (via base path). PROVED. -/
theorem reaches_zero_audit {n : ℕ}
    (T : Tournament n) (hbp : HasBasePath T) (hn : 0 < n) (u : Fin n) :
    Reaches T u ⟨0, hn⟩ :=
  reaches_zero T hbp hn u
#print axioms reaches_zero_audit

/-- THM-333 (apex tile is SC): if the apex arc 0 → (n-1) is present, T is SC. -/
theorem apex_implies_SC_audit {n : ℕ} (T : Tournament n)
    (hbp : HasBasePath T) (hn : 3 ≤ n) (hv0 : 0 < n) (hvn : n - 1 < n)
    (h_apex : T.arc ⟨0, hv0⟩ ⟨n - 1, hvn⟩ = true) :
    IsStronglyConnected T :=
  apex_implies_SC T hbp hn hv0 hvn h_apex
#print axioms apex_implies_SC_audit

/-! ### Self-complementary (clean isomorphism formulation) -/

/-- IsSelfComplementary ↔ T ≅ op T (via the new TournamentIso structure). -/
theorem isSelfComplementary_iff_iso_op_audit {n : ℕ} (T : Tournament n) :
    IsSelfComplementary T ↔ T ≅ op T :=
  isSelfComplementary_iff_iso_op T
#print axioms isSelfComplementary_iff_iso_op_audit

/-! ### Regular ⟹ ¬ SF chain (THM-345 candidate) -/

/-- Any regular base-path tournament is not self-flip via identity. -/
theorem regular_not_SF_id_audit {n : ℕ} (T : Tournament n)
    (hbp : HasBasePath T) (hn : 3 ≤ n) (hreg : IsRegular T) :
    ¬ IsSelfFlip_id T :=
  regular_not_SF_id T hbp hn hreg
#print axioms regular_not_SF_id_audit

/-- Any Paley-like tournament is not self-flip via identity. -/
theorem paleyLike_not_SF_audit {n : ℕ} (P : PaleyLike n) :
    ¬ IsSelfFlip_id P.T :=
  paleyLike_not_SF_id P
#print axioms paleyLike_not_SF_audit

/-! ### THM-326 (HP = IP universal identity) -/

theorem hp_ip_truncated_audit {n : ℕ} (T : Tournament n) :
    H T = 1 + 2 * alphaCount 1 T + 4 * alphaCount 2 T
            + 8 * alphaCount 3 T + 16 * alphaCount 4 T :=
  H_eq_independence_poly_at_two_truncated T
#print axioms hp_ip_truncated_audit

/-! ### THM-316 (Abstract anti-palindrome) -/

theorem abstract_anti_palindrome_audit {n : ℕ}
    (T : Tournament n) (hn : 0 < n) (φ : Equiv.Perm (Fin n))
    (hφ : IsAntiAutomorphism T φ) (v : Fin n) :
    epStart T hn v = epEnd T hn (φ v) :=
  abstract_anti_palindrome T hn φ hφ v
#print axioms abstract_anti_palindrome_audit

/-- Endpoint-start fibers partition the Hamiltonian paths. -/
theorem epStart_sum_eq_H_audit {n : ℕ} (T : Tournament n) (hn : 0 < n) :
    ∑ v : Fin n, epStart T hn v = H T :=
  epStart_sum_eq_H T hn
#print axioms epStart_sum_eq_H_audit

/-- Endpoint-end fibers partition the Hamiltonian paths. -/
theorem epEnd_sum_eq_H_audit {n : ℕ} (T : Tournament n) (hn : 0 < n) :
    ∑ v : Fin n, epEnd T hn v = H T :=
  epEnd_sum_eq_H T hn
#print axioms epEnd_sum_eq_H_audit

/-! ### Isomorphism invariants -/

/-- Tournament isomorphism preserves vertex out-degrees (up to relabelling).
    PROVED IN LEAN (no axiom). -/
theorem outDegree_iso_audit {n : ℕ}
    (T₁ T₂ : Tournament n) (f : TournamentIso T₁ T₂) (v : Fin n) :
    T₁.outDegree v = T₂.outDegree (f.perm v) :=
  outDegree_iso T₁ T₂ f v
#print axioms outDegree_iso_audit

/-- Tournament isomorphism preserves the regularity property. PROVED. -/
theorem isRegular_iso_audit {n : ℕ}
    (T₁ T₂ : Tournament n) (h : T₁ ≅ T₂) (hreg : IsRegular T₁) :
    IsRegular T₂ :=
  isRegular_iso T₁ T₂ h hreg
#print axioms isRegular_iso_audit

/-- H is an isomorphism invariant — PROVED IN LEAN. -/
theorem H_iso_invariant_audit {n : ℕ}
    (T₁ T₂ : Tournament n) (h : T₁ ≅ T₂) :
    H T₁ = H T₂ :=
  H_iso_invariant T₁ T₂ h
#print axioms H_iso_invariant_audit

/-- alphaCount is an isomorphism invariant. -/
theorem alphaCount_iso_invariant_audit {n : ℕ}
    (k : ℕ) (T₁ T₂ : Tournament n) (h : T₁ ≅ T₂) :
    alphaCount k T₁ = alphaCount k T₂ :=
  alphaCount_iso_invariant k T₁ T₂ h
#print axioms alphaCount_iso_invariant_audit

/-! ### Bucket-balance half-line conservation (THM-346 core) -/

/-- Ordered half-line bucket balance for any finite quotient and finite move set.
    PROVED. -/
theorem bucketBalance_halfLine_balance_audit
    {X M B : Type*} [Fintype X] [DecidableEq B]
    (q : X → B) (step : M → X → X) (moves : Finset M) (b : B) :
    (BucketBalance.selfHalf q step moves b).card +
      (BucketBalance.crossHalf q step moves b).card =
        (BucketBalance.fiber q b).card * moves.card :=
  BucketBalance.halfLine_balance q step moves b
#print axioms bucketBalance_halfLine_balance_audit

/-! ### Good-cut buckets for staircase tilings -/

/-- Good-cut bucket 0 is exactly the all-down tiling. PROVED. -/
theorem goodCuts_empty_iff_all_down_audit {n : ℕ} (b : StTiling n) :
    b.goodCuts = ∅ ↔ ∀ t : StTile n, b t = false :=
  StTiling.goodCuts_empty_iff_all_down b
#print axioms goodCuts_empty_iff_all_down_audit

/-- Good-cut support is nonempty iff some tile is upward. PROVED. -/
theorem goodCuts_nonempty_iff_exists_upward_tile_audit {n : ℕ} (b : StTiling n) :
    b.goodCuts.Nonempty ↔ ∃ t : StTile n, b t = true :=
  StTiling.goodCuts_nonempty_iff_exists_upward_tile b
#print axioms goodCuts_nonempty_iff_exists_upward_tile_audit

/-- Good-cut bucket 0, cardinality form. PROVED. -/
theorem goodCutCount_eq_zero_iff_all_down_audit {n : ℕ} (b : StTiling n) :
    b.goodCutCount = 0 ↔ ∀ t : StTile n, b t = false :=
  StTiling.goodCutCount_eq_zero_iff_all_down b
#print axioms goodCutCount_eq_zero_iff_all_down_audit

/-- Positive good-cut count iff some tile is upward. PROVED. -/
theorem goodCutCount_pos_iff_exists_upward_tile_audit {n : ℕ} (b : StTiling n) :
    0 < b.goodCutCount ↔ ∃ t : StTile n, b t = true :=
  StTiling.goodCutCount_pos_iff_exists_upward_tile b
#print axioms goodCutCount_pos_iff_exists_upward_tile_audit

/-- Positive good-cut count iff the tiling is not all-down. PROVED. -/
theorem goodCutCount_pos_iff_not_all_down_audit {n : ℕ} (b : StTiling n) :
    0 < b.goodCutCount ↔ ¬ ∀ t : StTile n, b t = false :=
  StTiling.goodCutCount_pos_iff_not_all_down b
#print axioms goodCutCount_pos_iff_not_all_down_audit

/-- One upward tile forces at least two good cuts. PROVED. -/
theorem two_le_goodCutCount_of_upward_tile_audit {n : ℕ}
    {b : StTiling n} {t : StTile n} (ht : b t = true) :
    2 ≤ b.goodCutCount :=
  StTiling.two_le_goodCutCount_of_upward_tile ht
#print axioms two_le_goodCutCount_of_upward_tile_audit

/-- THM-336 Lean core: no tiling has exactly one good cut. PROVED. -/
theorem goodCutCount_ne_one_audit {n : ℕ} (b : StTiling n) :
    b.goodCutCount ≠ 1 :=
  StTiling.goodCutCount_ne_one b
#print axioms goodCutCount_ne_one_audit

/-- THM-336 strengthened: bucket count is 0 or at least 2. PROVED. -/
theorem goodCutCount_eq_zero_or_two_le_audit {n : ℕ} (b : StTiling n) :
    b.goodCutCount = 0 ∨ 2 ≤ b.goodCutCount :=
  StTiling.goodCutCount_eq_zero_or_two_le b
#print axioms goodCutCount_eq_zero_or_two_le_audit

/-- Good-cut set form: empty or cardinality at least two. PROVED. -/
theorem goodCuts_empty_or_two_le_card_audit {n : ℕ} (b : StTiling n) :
    b.goodCuts = ∅ ∨ 2 ≤ b.goodCuts.card :=
  StTiling.goodCuts_empty_or_two_le_card b
#print axioms goodCuts_empty_or_two_le_card_audit

/-- Grid reflection preserves the good-cut bucket. PROVED. -/
theorem goodCutCount_reflect_audit {n : ℕ} (b : StTiling n) :
    b.reflect.goodCutCount = b.goodCutCount :=
  StTiling.goodCutCount_reflect b
#print axioms goodCutCount_reflect_audit

/-- Good cuts are exactly membership in an upward tile interval. PROVED. -/
theorem isGoodCut_interval_union_audit {n : ℕ} {b : StTiling n} {k : ℕ} :
    StTiling.IsGoodCut b k ↔
      ∃ t : StTile n, b t = true ∧ k ∈ t.cutInterval :=
  StTiling.isGoodCut_iff_exists_upward_tile_interval
#print axioms isGoodCut_interval_union_audit

theorem mem_goodCuts_interval_union_audit {n : ℕ} {b : StTiling n} {k : ℕ} :
    k ∈ b.goodCuts ↔
      ∃ t : StTile n, b t = true ∧ k ∈ t.cutInterval :=
  StTiling.mem_goodCuts_iff_exists_upward_tile_interval
#print axioms mem_goodCuts_interval_union_audit

theorem cutInterval_subset_goodCuts_audit {n : ℕ} {b : StTiling n} {t : StTile n}
    (ht : b t = true) :
    t.cutInterval ⊆ b.goodCuts :=
  StTiling.cutInterval_subset_goodCuts_of_upward_tile ht
#print axioms cutInterval_subset_goodCuts_audit

/-- Good-cut count is monotone under turning more tiles upward. PROVED. -/
theorem goodCutCount_mono_audit {n : ℕ} {b c : StTiling n}
    (h : ∀ t : StTile n, b t = true → c t = true) :
    b.goodCutCount ≤ c.goodCutCount :=
  StTiling.goodCutCount_mono h
#print axioms goodCutCount_mono_audit

/-- The only possible buckets are 0 or at least 2, bounded above by n-1. PROVED. -/
theorem goodCutCount_bucket_bounds_audit {n : ℕ} (b : StTiling n) :
    b.goodCutCount = 0 ∨
      (2 ≤ b.goodCutCount ∧ b.goodCutCount ≤ n - 1) :=
  StTiling.goodCutCount_bucket_bounds b
#print axioms goodCutCount_bucket_bounds_audit

/-- The top bucket is equivalent to every legal cut being good. PROVED. -/
theorem goodCutCount_eq_top_iff_all_cuts_good_audit {n : ℕ} (b : StTiling n) :
    b.goodCutCount = n - 1 ↔
      ∀ k, k ∈ cutSet n → StTiling.IsGoodCut b k :=
  StTiling.goodCutCount_eq_top_iff_all_cuts_good b
#print axioms goodCutCount_eq_top_iff_all_cuts_good_audit

/-- For n >= 3, every legal cut is crossed by some staircase tile. PROVED. -/
theorem exists_crossesCut_of_mem_cutSet_audit {n : ℕ} (hn : 3 ≤ n) {k : ℕ}
    (hk : k ∈ cutSet n) :
    ∃ t : StTile n, t.crossesCut k :=
  StTile.exists_crossesCut_of_mem_cutSet hn hk
#print axioms exists_crossesCut_of_mem_cutSet_audit

/-- A single-up tiling has good cuts exactly the interval crossed by that tile. PROVED. -/
theorem goodCuts_singleUp_eq_cutInterval_audit {n : ℕ} (t : StTile n) :
    (StTiling.singleUp t).goodCuts = t.cutInterval :=
  StTiling.goodCuts_singleUp_eq_cutInterval t
#print axioms goodCuts_singleUp_eq_cutInterval_audit

/-- Every allowed nonzero good-cut bucket size is realized. PROVED. -/
theorem exists_goodCutCount_eq_of_allowed_audit {n r : ℕ}
    (hn : 3 ≤ n) (hr2 : 2 ≤ r) (hrn : r ≤ n - 1) :
    ∃ b : StTiling n, b.goodCutCount = r :=
  StTiling.exists_goodCutCount_eq_of_allowed hn hr2 hrn
#print axioms exists_goodCutCount_eq_of_allowed_audit

/-- Exact spectrum of the good-cut bucket abstraction. PROVED. -/
theorem goodCutCount_spectrum_audit {n r : ℕ} (hn : 3 ≤ n) :
    (∃ b : StTiling n, b.goodCutCount = r) ↔
      r = 0 ∨ (2 ≤ r ∧ r ≤ n - 1) :=
  StTiling.goodCutCount_spectrum hn
#print axioms goodCutCount_spectrum_audit

/-- For n >= 3, the all-up tiling is in the top good-cut bucket. PROVED. -/
theorem goodCutCount_allUp_audit {n : ℕ} (hn : 3 ≤ n) :
    (StTiling.allUp n).goodCutCount = n - 1 :=
  StTiling.goodCutCount_allUp hn
#print axioms goodCutCount_allUp_audit

/-- Complementing an all-down tiling puts it in the top bucket. PROVED. -/
theorem goodCutCount_complement_of_all_down_audit {n : ℕ}
    {b : StTiling n} (hn : 3 ≤ n) (h : ∀ t : StTile n, b t = false) :
    b.complement.goodCutCount = n - 1 :=
  StTiling.goodCutCount_complement_of_all_down hn h
#print axioms goodCutCount_complement_of_all_down_audit

/-! ### Concrete staircase tilings as tournaments -/

/-- A concrete staircase tiling induces a valid tournament. PROVED. -/
theorem staircase_toTournament_hasBasePath_audit {n : ℕ} (b : StTiling n) :
    HasBasePath b.toTournament :=
  StTiling.toTournament_hasBasePath b
#print axioms staircase_toTournament_hasBasePath_audit

/-- Good cuts are exactly crossing-upward cuts of the induced tournament. PROVED. -/
theorem isGoodCut_iff_crossesUpward_toTournament_audit {n : ℕ}
    (b : StTiling n) {k : ℕ} :
    StTiling.IsGoodCut b k ↔ CrossesUpward b.toTournament k :=
  StTiling.isGoodCut_iff_crossesUpward_toTournament b
#print axioms isGoodCut_iff_crossesUpward_toTournament_audit

/-- Top good-cut bucket iff strong connectivity of the induced tournament. PROVED. -/
theorem goodCutCount_eq_top_iff_toTournament_SC_audit {n : ℕ}
    (b : StTiling n) :
    b.goodCutCount = n - 1 ↔ IsStronglyConnected b.toTournament :=
  StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected b
#print axioms goodCutCount_eq_top_iff_toTournament_SC_audit

/-- The all-up tiling gives an explicit strongly connected staircase
    tournament for n≥3. PROVED. -/
theorem allUp_toTournament_SC_audit {n : ℕ} (hn : 3 ≤ n) :
    IsStronglyConnected (StTiling.allUp n).toTournament :=
  StTiling.allUp_toTournament_stronglyConnected hn
#print axioms allUp_toTournament_SC_audit

/-- The all-down tiling gives an explicit non-strongly-connected staircase
    tournament for n≥2. PROVED. -/
theorem allDown_toTournament_not_SC_audit {n : ℕ} (hn : 2 ≤ n) :
    ¬ IsStronglyConnected (StTiling.allDown n).toTournament :=
  StTiling.allDown_toTournament_not_stronglyConnected hn
#print axioms allDown_toTournament_not_SC_audit

/-! ### Apex tile bridge -/

/-- The apex tile crosses every legal cut. PROVED. -/
theorem apexTile_cutInterval_eq_cutSet_audit {n : ℕ} (hn : 3 ≤ n) :
    (apexTile n hn).cutInterval = cutSet n :=
  apexTile_cutInterval_eq_cutSet n hn
#print axioms apexTile_cutInterval_eq_cutSet_audit

/-- A single upward apex tile makes every legal cut good. PROVED. -/
theorem singleUp_apex_goodCuts_eq_cutSet_audit {n : ℕ} (hn : 3 ≤ n) :
    (StTiling.singleUp (apexTile n hn)).goodCuts = cutSet n :=
  singleUp_apex_goodCuts_eq_cutSet n hn
#print axioms singleUp_apex_goodCuts_eq_cutSet_audit

/-- A single upward apex tile lies in the top good-cut bucket. PROVED. -/
theorem singleUp_apex_goodCutCount_top_audit {n : ℕ} (hn : 3 ≤ n) :
    (StTiling.singleUp (apexTile n hn)).goodCutCount = n - 1 :=
  singleUp_apex_goodCutCount_top n hn
#print axioms singleUp_apex_goodCutCount_top_audit

/-! ### Abstract bucket balances -/

/-- THM-346 Lean core: oriented half-lines split into internal and escaping
    half-lines. PROVED. -/
theorem bucket_halfLine_balance_audit {alpha beta move : Type}
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    (BucketBalance.selfHalf q step moves b).card +
      (BucketBalance.crossHalf q step moves b).card =
        (BucketBalance.fiber q b).card * moves.card :=
  BucketBalance.halfLine_balance q step moves b
#print axioms bucket_halfLine_balance_audit

/-- Bucket closure criterion: zero escaping half-lines iff every chosen move
    from the fiber stays in the same bucket. PROVED. -/
theorem bucket_crossHalf_card_eq_zero_iff_audit {alpha beta move : Type}
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    (BucketBalance.crossHalf q step moves b).card = 0 <->
      forall x, x ∈ BucketBalance.fiber q b ->
        forall u, u ∈ moves -> q (step u x) = b :=
  BucketBalance.crossHalf_card_eq_zero_iff q step moves b
#print axioms bucket_crossHalf_card_eq_zero_iff_audit

theorem bucket_crossHalf_card_le_total_audit {alpha beta move : Type}
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    (BucketBalance.crossHalf q step moves b).card <=
      (BucketBalance.fiber q b).card * moves.card :=
  BucketBalance.crossHalf_card_le_total q step moves b
#print axioms bucket_crossHalf_card_le_total_audit

/-- Partnering an internal half-line by an involutive move keeps it internal.
    PROVED. -/
theorem bucket_pairHalf_mem_selfHalf_audit {alpha beta move : Type}
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hstep : forall u, u ∈ moves -> Function.Involutive (step u))
    {xu : alpha × move} (hxu : xu ∈ BucketBalance.selfHalf q step moves b) :
    BucketBalance.pairHalf step xu ∈ BucketBalance.selfHalf q step moves b :=
  BucketBalance.pairHalf_mem_selfHalf q step moves b hstep hxu
#print axioms bucket_pairHalf_mem_selfHalf_audit

/-- A fixed-point-free selected move never pairs an internal half-line with
    itself. PROVED. -/
theorem bucket_pairHalf_ne_of_fixedPointFree_audit {alpha beta move : Type}
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hfixed : forall u, u ∈ moves -> forall x, step u x ≠ x)
    {xu : alpha × move} (hxu : xu ∈ BucketBalance.selfHalf q step moves b) :
    BucketBalance.pairHalf step xu ≠ xu :=
  BucketBalance.pairHalf_ne_of_fixedPointFree q step moves b hfixed hxu
#print axioms bucket_pairHalf_ne_of_fixedPointFree_audit

/-- A finite fixed-point-free involution has even cardinality. PROVED. -/
theorem bucket_even_card_of_fixedPointFree_involutiveOn_audit {alpha : Type}
    [DecidableEq alpha] (s : Finset alpha) (f : alpha -> alpha)
    (hmem : forall x, x ∈ s -> f x ∈ s)
    (hinv : forall x, x ∈ s -> f (f x) = x)
    (hfixed : forall x, x ∈ s -> f x ≠ x) :
    Even s.card :=
  BucketBalance.even_card_of_fixedPointFree_involutiveOn s f hmem hinv hfixed
#print axioms bucket_even_card_of_fixedPointFree_involutiveOn_audit

/-- Internal half-lines have even cardinality for fixed-point-free involutive
    move systems. PROVED. -/
theorem bucket_selfHalf_card_even_of_involutive_fixedPointFree_audit
    {alpha beta move : Type}
    [Fintype alpha] [DecidableEq alpha] [DecidableEq beta] [DecidableEq move]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hstep : forall u, u ∈ moves -> Function.Involutive (step u))
    (hfixed : forall u, u ∈ moves -> forall x, step u x ≠ x) :
    Even (BucketBalance.selfHalf q step moves b).card :=
  BucketBalance.selfHalf_card_even_of_involutive_fixedPointFree
    q step moves b hstep hfixed
#print axioms bucket_selfHalf_card_even_of_involutive_fixedPointFree_audit

/-- The unordered balance follows once internal oriented half-lines have even
    cardinality. PROVED. -/
theorem bucket_unordered_balance_of_even_selfHalf_audit {alpha beta move : Type}
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hself : Even (BucketBalance.selfHalf q step moves b).card) :
    2 * BucketBalance.internalLineCount q step moves b +
        (BucketBalance.crossHalf q step moves b).card =
      (BucketBalance.fiber q b).card * moves.card :=
  BucketBalance.unordered_balance_of_even_selfHalf q step moves b hself
#print axioms bucket_unordered_balance_of_even_selfHalf_audit

/-- THM-350 strengthened: fixed-point-free involutive move systems satisfy the
    unordered bucket balance without a separate evenness assumption. PROVED. -/
theorem bucket_unordered_balance_of_involutive_fixedPointFree_audit
    {alpha beta move : Type}
    [Fintype alpha] [DecidableEq alpha] [DecidableEq beta] [DecidableEq move]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hstep : forall u, u ∈ moves -> Function.Involutive (step u))
    (hfixed : forall u, u ∈ moves -> forall x, step u x ≠ x) :
    2 * BucketBalance.internalLineCount q step moves b +
        (BucketBalance.crossHalf q step moves b).card =
      (BucketBalance.fiber q b).card * moves.card :=
  BucketBalance.unordered_balance_of_involutive_fixedPointFree
    q step moves b hstep hfixed
#print axioms bucket_unordered_balance_of_involutive_fixedPointFree_audit

/-- Xor by a Boolean mask is involutive. PROVED. -/
theorem boolCube_xorMask_involutive_audit {index : Type}
    (u : BucketBalance.BoolCube index) :
    Function.Involutive (BucketBalance.xorMask u) :=
  BucketBalance.xorMask_involutive u
#print axioms boolCube_xorMask_involutive_audit

/-- Xor by a nonzero Boolean mask is fixed-point-free. PROVED. -/
theorem boolCube_xorMask_fixedPointFree_audit {index : Type}
    {u : BucketBalance.BoolCube index}
    (hu : BucketBalance.IsNonzeroMask u) :
    ∀ x, BucketBalance.xorMask u x ≠ x :=
  BucketBalance.xorMask_fixedPointFree_of_nonzero hu
#print axioms boolCube_xorMask_fixedPointFree_audit

/-- THM-346 Lean specialization: bucket balance for finite Boolean cubes and
    nonzero xor-mask families. PROVED. -/
theorem bucket_unordered_balance_boolCube_masks_audit
    {index beta : Type} [Fintype index] [DecidableEq index] [DecidableEq beta]
    (q : BucketBalance.BoolCube index -> beta)
    (moves : Finset (BucketBalance.BoolCube index)) (b : beta)
    (hmoves : forall u, u ∈ moves -> BucketBalance.IsNonzeroMask u) :
    2 * BucketBalance.internalLineCount q BucketBalance.xorMask moves b +
        (BucketBalance.crossHalf q BucketBalance.xorMask moves b).card =
      (BucketBalance.fiber q b).card * moves.card :=
  BucketBalance.unordered_balance_boolCube_masks q moves b hmoves
#print axioms bucket_unordered_balance_boolCube_masks_audit

/-- Off-diagonal target-bucket transport rows are exactly escaping half-lines.
    PROVED. -/
theorem bucket_transport_offdiag_eq_crossHalf_audit {alpha beta move : Type}
    [Fintype alpha] [Fintype beta] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    (∑ c ∈ (Finset.univ.erase b),
        (BucketBalance.transportHalf q step moves b c).card) =
      (BucketBalance.crossHalf q step moves b).card :=
  BucketBalance.sum_transportHalf_card_offdiag_eq_crossHalf_card q step moves b
#print axioms bucket_transport_offdiag_eq_crossHalf_audit

/-- Target-bucket transport row checksum for fixed-point-free involutive move
    systems. PROVED. -/
theorem bucket_transport_row_balance_fixedPointFree_audit
    {alpha beta move : Type}
    [Fintype alpha] [Fintype beta] [DecidableEq alpha] [DecidableEq beta]
    [DecidableEq move]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hstep : forall u, u ∈ moves -> Function.Involutive (step u))
    (hfixed : forall u, u ∈ moves -> forall x, step u x ≠ x) :
    2 * BucketBalance.internalLineCount q step moves b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf q step moves b c).card) =
      (BucketBalance.fiber q b).card * moves.card :=
  BucketBalance.transport_row_balance_of_involutive_fixedPointFree
    q step moves b hstep hfixed
#print axioms bucket_transport_row_balance_fixedPointFree_audit

/-- Boolean-cube target-bucket transport row checksum for nonzero xor-mask
    families. PROVED. -/
theorem bucket_transport_row_balance_boolCube_masks_audit
    {index beta : Type} [Fintype index] [DecidableEq index]
    [Fintype beta] [DecidableEq beta]
    (q : BucketBalance.BoolCube index -> beta)
    (moves : Finset (BucketBalance.BoolCube index)) (b : beta)
    (hmoves : forall u, u ∈ moves -> BucketBalance.IsNonzeroMask u) :
    2 * BucketBalance.internalLineCount q BucketBalance.xorMask moves b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf q BucketBalance.xorMask moves b c).card) =
      (BucketBalance.fiber q b).card * moves.card :=
  BucketBalance.transport_row_balance_boolCube_masks q moves b hmoves
#print axioms bucket_transport_row_balance_boolCube_masks_audit

/-! ### Bucket gaps and empty transport rows/columns -/

theorem bucket_fiber_eq_empty_iff_audit {alpha beta : Type}
    [Fintype alpha] [DecidableEq beta] (q : alpha -> beta) (b : beta) :
    BucketBalance.fiber q b = ∅ ↔ ∀ x, q x ≠ b :=
  BucketBalance.fiber_eq_empty_iff q b
#print axioms bucket_fiber_eq_empty_iff_audit

theorem bucket_incidentHalf_empty_of_gap_audit {alpha beta move : Type}
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (moves : Finset move) (b : beta)
    (h : BucketBalance.fiber q b = ∅) :
    BucketBalance.incidentHalf q moves b = ∅ :=
  BucketBalance.incidentHalf_eq_empty_of_fiber_eq_empty q moves b h
#print axioms bucket_incidentHalf_empty_of_gap_audit

theorem bucket_transportHalf_empty_of_source_gap_audit
    {alpha beta move : Type} [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b c : beta)
    (h : BucketBalance.fiber q b = ∅) :
    BucketBalance.transportHalf q step moves b c = ∅ :=
  BucketBalance.transportHalf_eq_empty_of_source_fiber_eq_empty q step moves b c h
#print axioms bucket_transportHalf_empty_of_source_gap_audit

theorem bucket_transportHalf_empty_of_target_gap_audit
    {alpha beta move : Type} [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b c : beta)
    (h : BucketBalance.fiber q c = ∅) :
    BucketBalance.transportHalf q step moves b c = ∅ :=
  BucketBalance.transportHalf_eq_empty_of_target_fiber_eq_empty q step moves b c h
#print axioms bucket_transportHalf_empty_of_target_gap_audit

theorem bucket_transportHalf_card_zero_of_source_gap_audit
    {alpha beta move : Type} [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b c : beta)
    (h : BucketBalance.fiber q b = ∅) :
    (BucketBalance.transportHalf q step moves b c).card = 0 :=
  BucketBalance.transportHalf_card_eq_zero_of_source_fiber_eq_empty
    q step moves b c h
#print axioms bucket_transportHalf_card_zero_of_source_gap_audit

theorem bucket_transportHalf_card_zero_of_target_gap_audit
    {alpha beta move : Type} [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b c : beta)
    (h : BucketBalance.fiber q c = ∅) :
    (BucketBalance.transportHalf q step moves b c).card = 0 :=
  BucketBalance.transportHalf_card_eq_zero_of_target_fiber_eq_empty
    q step moves b c h
#print axioms bucket_transportHalf_card_zero_of_target_gap_audit

theorem bucket_transportHalf_row_sum_zero_of_source_gap_audit
    {alpha beta move : Type} [Fintype alpha] [Fintype beta] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (h : BucketBalance.fiber q b = ∅) :
    (∑ c : beta, (BucketBalance.transportHalf q step moves b c).card) = 0 :=
  BucketBalance.sum_transportHalf_card_eq_zero_of_source_fiber_eq_empty
    q step moves b h
#print axioms bucket_transportHalf_row_sum_zero_of_source_gap_audit

theorem bucket_transportHalf_column_sum_zero_of_target_gap_audit
    {alpha beta move : Type} [Fintype alpha] [Fintype beta] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (c : beta)
    (h : BucketBalance.fiber q c = ∅) :
    (∑ b : beta, (BucketBalance.transportHalf q step moves b c).card) = 0 :=
  BucketBalance.sum_transportHalf_card_eq_zero_of_target_fiber_eq_empty
    q step moves c h
#print axioms bucket_transportHalf_column_sum_zero_of_target_gap_audit

/-! ### Concrete staircase bucket transport -/

theorem stTile_gapPair_roundtrip_audit {n : ℕ} (t : StTile n) :
    StTile.ofGapPair (StTile.toGapPair t) = t :=
  StTile.ofGapPair_toGapPair t
#print axioms stTile_gapPair_roundtrip_audit

theorem stTiling_singleUp_isNonzeroMask_audit {n : ℕ} (t : StTile n) :
    BucketBalance.IsNonzeroMask (StTiling.singleUp t) :=
  StTiling.singleUp_isNonzeroMask t
#print axioms stTiling_singleUp_isNonzeroMask_audit

theorem stTiling_allUp_isNonzeroMask_audit {n : ℕ} (hn : 3 ≤ n) :
    BucketBalance.IsNonzeroMask (StTiling.allUp n) :=
  StTiling.allUp_isNonzeroMask_of_three_le hn
#print axioms stTiling_allUp_isNonzeroMask_audit

theorem stTiling_transport_row_balance_allNonzeroMasks_audit
    {n : ℕ} {beta : Type} [Fintype beta] [DecidableEq beta]
    (q : StTiling n -> beta) (b : beta) :
    2 * BucketBalance.internalLineCount q BucketBalance.xorMask
        (StTiling.nonzeroMasks n) b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf q BucketBalance.xorMask
            (StTiling.nonzeroMasks n) b c).card) =
      (BucketBalance.fiber q b).card * (StTiling.nonzeroMasks n).card :=
  StTiling.transport_row_balance_allNonzeroMasks q b
#print axioms stTiling_transport_row_balance_allNonzeroMasks_audit

theorem stTiling_transport_row_balance_singleTileMasks_audit
    {n : ℕ} {beta : Type} [Fintype beta] [DecidableEq beta]
    (q : StTiling n -> beta) (b : beta) :
    2 * BucketBalance.internalLineCount q BucketBalance.xorMask
        (StTiling.singleTileMasks n) b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf q BucketBalance.xorMask
            (StTiling.singleTileMasks n) b c).card) =
      (BucketBalance.fiber q b).card * (StTiling.singleTileMasks n).card :=
  StTiling.transport_row_balance_singleTileMasks q b
#print axioms stTiling_transport_row_balance_singleTileMasks_audit

theorem stTiling_transport_row_balance_complementMask_audit
    {n : ℕ} {beta : Type} [Fintype beta] [DecidableEq beta]
    (q : StTiling n -> beta) (b : beta) (hn : 3 ≤ n) :
    2 * BucketBalance.internalLineCount q BucketBalance.xorMask
        (StTiling.complementMask n) b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf q BucketBalance.xorMask
            (StTiling.complementMask n) b c).card) =
      (BucketBalance.fiber q b).card * (StTiling.complementMask n).card :=
  StTiling.transport_row_balance_complementMask q b hn
#print axioms stTiling_transport_row_balance_complementMask_audit

theorem stTiling_goodCutBucket_zero_iff_all_down_audit {n : ℕ}
    (u : StTiling n) :
    StTiling.goodCutBucket u = (0 : Fin (n + 1)) ↔
      ∀ t : StTile n, u t = false :=
  StTiling.goodCutBucket_eq_zero_iff_all_down u
#print axioms stTiling_goodCutBucket_zero_iff_all_down_audit

theorem stTiling_goodCutBucket_top_iff_SC_audit {n : ℕ}
    (u : StTiling n) :
    StTiling.goodCutBucket u = StTiling.topGoodCutBucket n ↔
      IsStronglyConnected u.toTournament :=
  StTiling.goodCutBucket_eq_top_iff_toTournament_SC u
#print axioms stTiling_goodCutBucket_top_iff_SC_audit

theorem stTiling_goodCutBucket_image_iff_audit {n : ℕ}
    (hn : 3 ≤ n) (b : Fin (n + 1)) :
    (∃ u : StTiling n, StTiling.goodCutBucket u = b) ↔
      b.val = 0 ∨ (2 ≤ b.val ∧ b.val ≤ n - 1) :=
  StTiling.goodCutBucket_image_iff hn b
#print axioms stTiling_goodCutBucket_image_iff_audit

theorem stTiling_goodCutBucket_ne_one_audit {n : ℕ}
    (hn : 1 ≤ n) (u : StTiling n) :
    StTiling.goodCutBucket u ≠ (⟨1, by omega⟩ : Fin (n + 1)) :=
  StTiling.goodCutBucket_ne_one hn u
#print axioms stTiling_goodCutBucket_ne_one_audit

theorem stTiling_goodCutBucket_ne_overTop_audit {n : ℕ}
    (hn : 1 ≤ n) (u : StTiling n) :
    StTiling.goodCutBucket u ≠ (⟨n, by omega⟩ : Fin (n + 1)) :=
  StTiling.goodCutBucket_ne_overTop hn u
#print axioms stTiling_goodCutBucket_ne_overTop_audit

theorem stTiling_goodCutBucket_fiber_one_empty_audit {n : ℕ}
    (hn : 1 ≤ n) :
    BucketBalance.fiber StTiling.goodCutBucket
      (⟨1, by omega⟩ : Fin (n + 1)) = ∅ :=
  StTiling.goodCutBucket_fiber_one_eq_empty hn
#print axioms stTiling_goodCutBucket_fiber_one_empty_audit

theorem stTiling_goodCutBucket_fiber_overTop_empty_audit {n : ℕ}
    (hn : 1 ≤ n) :
    BucketBalance.fiber StTiling.goodCutBucket
      (⟨n, by omega⟩ : Fin (n + 1)) = ∅ :=
  StTiling.goodCutBucket_fiber_overTop_eq_empty hn
#print axioms stTiling_goodCutBucket_fiber_overTop_empty_audit

theorem stTiling_goodCutBucket_transport_row_singleTileMasks_audit
    {n : ℕ} (b : Fin (n + 1)) :
    2 * BucketBalance.internalLineCount StTiling.goodCutBucket BucketBalance.xorMask
        (StTiling.singleTileMasks n) b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf StTiling.goodCutBucket BucketBalance.xorMask
            (StTiling.singleTileMasks n) b c).card) =
      (BucketBalance.fiber StTiling.goodCutBucket b).card *
        (StTiling.singleTileMasks n).card :=
  StTiling.transport_row_balance_goodCutBucket_singleTileMasks b
#print axioms stTiling_goodCutBucket_transport_row_singleTileMasks_audit

theorem stTiling_goodCutBucket_transport_row_allNonzeroMasks_audit
    {n : ℕ} (b : Fin (n + 1)) :
    2 * BucketBalance.internalLineCount StTiling.goodCutBucket BucketBalance.xorMask
        (StTiling.nonzeroMasks n) b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf StTiling.goodCutBucket BucketBalance.xorMask
            (StTiling.nonzeroMasks n) b c).card) =
      (BucketBalance.fiber StTiling.goodCutBucket b).card *
        (StTiling.nonzeroMasks n).card :=
  StTiling.transport_row_balance_goodCutBucket_allNonzeroMasks b
#print axioms stTiling_goodCutBucket_transport_row_allNonzeroMasks_audit

theorem stTiling_goodCutBucket_transport_row_complementMask_audit
    {n : ℕ} (b : Fin (n + 1)) (hn : 3 ≤ n) :
    2 * BucketBalance.internalLineCount StTiling.goodCutBucket BucketBalance.xorMask
        (StTiling.complementMask n) b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf StTiling.goodCutBucket BucketBalance.xorMask
            (StTiling.complementMask n) b c).card) =
      (BucketBalance.fiber StTiling.goodCutBucket b).card *
        (StTiling.complementMask n).card :=
  StTiling.transport_row_balance_goodCutBucket_complementMask b hn
#print axioms stTiling_goodCutBucket_transport_row_complementMask_audit

/-! ### THM-342 (small diagonal value) -/

example : Qcount 2 1 = 1 := by
  have := thm342_diag0 1 (by omega); simp at this; exact this

/-! ### Concrete tournament examples (no axioms needed) -/

/-- The transitive tournament on 4 vertices has the base path. -/
theorem transitive_4_hasBasePath : HasBasePath (transitiveTournament 4) :=
  transitive_hasBasePath 4
#print axioms transitive_4_hasBasePath

/-- The 3-cycle tournament is regular. -/
theorem threeCycle_regular_audit : IsRegular threeCycle :=
  threeCycle_isRegular
#print axioms threeCycle_regular_audit

/-- The transitive tournament on n ≥ 2 vertices is NOT regular. -/
theorem transitive_not_regular_audit (n : ℕ) (hn : 2 ≤ n) :
    ¬ IsRegular (transitiveTournament n) :=
  transitive_not_regular n hn
#print axioms transitive_not_regular_audit

/-! ### N_min(k) = 3^k theorem -/

/-- For any tournament T with α_k ≥ 1 (k ∈ {1, 2, 3, 4}),
    H(T) ≥ 3^k.  Project-novel (oracle-S4). -/
theorem H_ge_three_pow_k_audit {n : ℕ}
    (T : Tournament n) (k : ℕ) (hk_pos : 1 ≤ k) (hk_le : k ≤ 4)
    (h : 1 ≤ alphaCount k T) :
    3 ^ k ≤ H T :=
  H_ge_three_pow_k_of_alpha_pos T k hk_pos hk_le h
#print axioms H_ge_three_pow_k_audit

/-- H(T) < 27 ⟹ no independent triple of vertex-disjoint odd cycles. -/
theorem H_lt_27_no_alpha3_audit {n : ℕ} (T : Tournament n) (hH : H T < 27) :
    alphaCount 3 T = 0 :=
  H_lt_27_no_alpha3 T hH
#print axioms H_lt_27_no_alpha3_audit

/-- H(T) < 81 ⟹ no independent quadruple of vertex-disjoint odd cycles. -/
theorem H_lt_81_no_alpha4_audit {n : ℕ} (T : Tournament n) (hH : H T < 81) :
    alphaCount 4 T = 0 :=
  H_lt_81_no_alpha4 T hH
#print axioms H_lt_81_no_alpha4_audit

/-! ### PaleyAxiomatic audits -/

theorem paley_7_not_SF_audit : ¬ IsSelfFlip_id paley_7.T := paley_7_not_SF
#print axioms paley_7_not_SF_audit

theorem paley_7_max_audit (T : Tournament 7) : H T ≤ 189 := paley_7_maximises_H T
#print axioms paley_7_max_audit

/-! ### StTile reflection audits (PROVED) -/

theorem stTile_reflect_reflect_audit (t : StTile n) : t.reflect.reflect = t :=
  StTile.reflect_reflect t
#print axioms stTile_reflect_reflect_audit

theorem stTiling_reflect_reflect_audit (b : StTiling n) :
    b.reflect.reflect = b := StTiling.reflect_reflect b
#print axioms stTiling_reflect_reflect_audit

theorem stTiling_complement_complement_audit (b : StTiling n) :
    b.complement.complement = b := StTiling.complement_complement b
#print axioms stTiling_complement_complement_audit

theorem stTiling_reflect_complement_commute_audit (b : StTiling n) :
    b.reflect.complement = b.complement.reflect :=
  StTiling.reflect_complement b
#print axioms stTiling_reflect_complement_commute_audit

/-! ### H = 3, H = 5 arithmetic enumeration (NEW) -/

/-- H = 3 forces α-tuple (1, 0, 0, 0). PROVED. -/
theorem alpha_solution_H3_audit {n : ℕ} (T : Tournament n) (h : H T = 3) :
    alphaCount 1 T = 1 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0 :=
  alpha_solution_H3 T h
#print axioms alpha_solution_H3_audit

/-- H = 5 forces α-tuple (2, 0, 0, 0). PROVED. -/
theorem alpha_solution_H5_audit {n : ℕ} (T : Tournament n) (h : H T = 5) :
    alphaCount 1 T = 2 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0 :=
  alpha_solution_H5 T h
#print axioms alpha_solution_H5_audit

/-! ### N_min(k) for k = 5, 6 (extended) -/

theorem H_ge_243_of_alpha5_pos_audit {n : ℕ} (T : Tournament n)
    (h : 1 ≤ alphaCount 5 T) : 243 ≤ H T :=
  H_ge_243_of_alpha5_pos T h
#print axioms H_ge_243_of_alpha5_pos_audit

theorem H_ge_729_of_alpha6_pos_audit {n : ℕ} (T : Tournament n)
    (h : 1 ≤ alphaCount 6 T) : 729 ≤ H T :=
  H_ge_729_of_alpha6_pos T h
#print axioms H_ge_729_of_alpha6_pos_audit

/-! ### Extended H-spectrum theorems -/

theorem H_spectrum_n7_audit (T : Tournament 7) :
    1 ≤ H T ∧ H T ≤ 189 ∧ Odd (H T) ∧ H T ≠ 7 ∧ H T ≠ 21 ∧ H T ≠ 63 :=
  H_spectrum_n7 T
#print axioms H_spectrum_n7_audit

theorem H_spectrum_n3_audit (T : Tournament 3) :
    1 ≤ H T ∧ Odd (H T) ∧ H T ≠ 7 ∧ H T ≠ 21 :=
  H_spectrum_n3 T
#print axioms H_spectrum_n3_audit

/-! ### Iso class refinement (numNS) -/

theorem numNS_eq_audit (n : ℕ) : numNS n = numIsoClasses n - numSC n :=
  numNS_eq n
#print axioms numNS_eq_audit

theorem numNS_plus_numSC_3_audit : numNS 3 + numSC 3 = 2 :=
  numNS_plus_numSC_3
#print axioms numNS_plus_numSC_3_audit

theorem numNS_plus_numSC_7_audit : numNS 7 + numSC 7 = 456 :=
  numNS_plus_numSC_7
#print axioms numNS_plus_numSC_7_audit

/-! ### Pascal's row sum (generalised) -/

theorem pascal_row_sum_general_audit (k : ℕ) :
    ∑ j ∈ Finset.range (k + 1), 2 ^ j * Nat.choose k j = 3 ^ k :=
  pascal_row_sum_general k
#print axioms pascal_row_sum_general_audit

/-! ### N_min(k) = 3^k for ARBITRARY k -/

/-- General N_min(k) = 3^k for any k. PROVED IN LEAN. -/
theorem H_ge_three_pow_k_general_audit {n : ℕ} (T : Tournament n) (k : ℕ)
    (hk_pos : 1 ≤ k) (h : 1 ≤ alphaCount k T)
    (h_high_zero : ∀ j, k < j → alphaCount j T = 0) :
    3 ^ k ≤ H T :=
  H_ge_three_pow_k_general T k hk_pos h h_high_zero
#print axioms H_ge_three_pow_k_general_audit

/-! ### Apex bridge -/

/-- The single-up-apex tiling at n ≥ 3 is strongly connected. PROVED. -/
theorem singleUp_apex_toTournament_SC_audit (n : ℕ) (hn : 3 ≤ n) :
    IsStronglyConnected (StTiling.singleUp (apexTile n hn)).toTournament :=
  singleUp_apex_toTournament_SC n hn
#print axioms singleUp_apex_toTournament_SC_audit

/-! ### Rédei's existence derived from OCF (no axiom!) -/

theorem H_ge_one_from_ocf_audit {n : ℕ} (T : Tournament n) : 1 ≤ H T :=
  H_ge_one T
#print axioms H_ge_one_from_ocf_audit

theorem H_ne_zero_from_ocf_audit {n : ℕ} (T : Tournament n) : H T ≠ 0 :=
  H_ne_zero T
#print axioms H_ne_zero_from_ocf_audit

theorem H_odd_from_ocf_audit {n : ℕ} (T : Tournament n) : Odd (H T) :=
  H_odd_from_ocf T
#print axioms H_odd_from_ocf_audit

theorem H_mod_two_eq_one_from_ocf_audit {n : ℕ} (T : Tournament n) :
    H T % 2 = 1 :=
  H_mod_two_eq_one_from_ocf T
#print axioms H_mod_two_eq_one_from_ocf_audit

theorem H_not_even_audit {n : ℕ} (T : Tournament n) : ¬ Even (H T) :=
  H_not_even T
#print axioms H_not_even_audit

/-! ### Clean H-spectrum theorems (no redei_existence/parity axiom) -/

theorem H_spectrum_universal_audit {n : ℕ} (T : Tournament n) :
    1 ≤ H T ∧ Odd (H T) ∧ H T ≠ 7 ∧ H T ≠ 21 :=
  H_spectrum_universal T
#print axioms H_spectrum_universal_audit

theorem alpha_solution_H1_audit {n : ℕ} (T : Tournament n) (h : H T = 1) :
    alphaCount 1 T = 0 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0 :=
  alpha_solution_H1 T h
#print axioms alpha_solution_H1_audit

/-! ### Forbidden-set theorem -/

theorem H_in_forbidden_set_audit {n : ℕ} (T : Tournament n) :
    H T ≠ 0 ∧ H T ≠ 2 ∧ H T ≠ 4 ∧ H T ≠ 6 ∧ H T ≠ 7 ∧ H T ≠ 8 ∧
    H T ≠ 10 ∧ H T ≠ 21 :=
  H_in_forbidden_set T
#print axioms H_in_forbidden_set_audit

/-! ### H-spectrum at small n -/

theorem H_n3_eq_one_or_three_audit (T : Tournament 3) : H T = 1 ∨ H T = 3 :=
  H_n3_eq_one_or_three T
#print axioms H_n3_eq_one_or_three_audit

theorem H_n4_eq_135_audit (T : Tournament 4) : H T = 1 ∨ H T = 3 ∨ H T = 5 :=
  H_n4_eq_135 T
#print axioms H_n4_eq_135_audit

theorem H_n5_in_spectrum_audit (T : Tournament 5) :
    H T = 1 ∨ H T = 3 ∨ H T = 5 ∨ H T = 9 ∨ H T = 11 ∨ H T = 13 ∨ H T = 15 :=
  H_n5_in_spectrum T
#print axioms H_n5_in_spectrum_audit

/-! ### Small H arithmetic enumerations -/

theorem alpha_solution_H9_audit {n : ℕ} (T : Tournament n) (h : H T = 9) :
    (alphaCount 1 T = 4 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 2 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) :=
  alpha_solution_H9 T h
#print axioms alpha_solution_H9_audit

theorem alpha_solution_H11_audit {n : ℕ} (T : Tournament n) (h : H T = 11) :
    (alphaCount 1 T = 5 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 3 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) :=
  alpha_solution_H11 T h
#print axioms alpha_solution_H11_audit

theorem alpha_solution_H13_audit {n : ℕ} (T : Tournament n) (h : H T = 13) :
    (alphaCount 1 T = 6 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 4 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) :=
  alpha_solution_H13 T h
#print axioms alpha_solution_H13_audit

theorem alpha_solution_H15_audit {n : ℕ} (T : Tournament n) (h : H T = 15) :
    (alphaCount 1 T = 7 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 5 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 3 ∧ alphaCount 2 T = 2
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) :=
  alpha_solution_H15 T h
#print axioms alpha_solution_H15_audit

theorem small_H_alpha34_zero_audit {n : ℕ} (T : Tournament n) (hH : H T ≤ 26) :
    alphaCount 3 T = 0 ∧ alphaCount 4 T = 0 :=
  small_H_alpha34_zero T hH
#print axioms small_H_alpha34_zero_audit

theorem alpha_solution_H17_audit {n : ℕ} (T : Tournament n) (h : H T = 17) :
    (alphaCount 1 T = 8 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 6 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 4 ∧ alphaCount 2 T = 2
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) :=
  alpha_solution_H17 T h
#print axioms alpha_solution_H17_audit

theorem alpha_solution_H19_audit {n : ℕ} (T : Tournament n) (h : H T = 19) :
    (alphaCount 1 T = 9 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 7 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 5 ∧ alphaCount 2 T = 2
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 3 ∧ alphaCount 2 T = 3
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) :=
  alpha_solution_H19 T h
#print axioms alpha_solution_H19_audit

theorem alpha_solution_H23_audit {n : ℕ} (T : Tournament n) (h : H T = 23) :
    (alphaCount 1 T = 11 ∧ alphaCount 2 T = 0
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 9 ∧ alphaCount 2 T = 1
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 7 ∧ alphaCount 2 T = 2
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) ∨
    (alphaCount 1 T = 5 ∧ alphaCount 2 T = 3
       ∧ alphaCount 3 T = 0 ∧ alphaCount 4 T = 0) :=
  alpha_solution_H23 T h
#print axioms alpha_solution_H23_audit

/-! ### BasePathSink: vertex 0 score bound -/

/-- For any HasBasePath T on n ≥ 2 vertices, vertex 0 has out-degree ≤ n - 2. -/
theorem base_path_sink_outDegree_le_audit {n : ℕ} (T : Tournament n)
    (hbp : HasBasePath T) (hn : 2 ≤ n) :
    T.outDegree ⟨0, by omega⟩ ≤ n - 2 :=
  base_path_sink_outDegree_le T hbp hn
#print axioms base_path_sink_outDegree_le_audit

theorem outDegree_le_n_minus_one_audit {n : ℕ} (T : Tournament n) (hn : 1 ≤ n)
    (v : Fin n) : T.outDegree v ≤ n - 1 :=
  outDegree_le_n_minus_one T hn v
#print axioms outDegree_le_n_minus_one_audit

/-- Vertex n-1 in HasBasePath has out-degree ≥ 1 (base path arc to n-2). -/
theorem base_path_source_outDegree_ge_audit {n : ℕ} (T : Tournament n)
    (hbp : HasBasePath T) (hn : 2 ≤ n) :
    1 ≤ T.outDegree ⟨n - 1, by omega⟩ :=
  base_path_source_outDegree_ge T hbp hn
#print axioms base_path_source_outDegree_ge_audit

/-- Regular HasBasePath tournament: n ≥ 3. -/
theorem regular_basepath_n_ge_three_audit {n : ℕ} (T : Tournament n)
    (hbp : HasBasePath T) (hn : 2 ≤ n) (hreg : IsRegular T) :
    3 ≤ n :=
  regular_basepath_n_ge_three T hbp hn hreg
#print axioms regular_basepath_n_ge_three_audit

/-! ### Iso-class characterizations -/

/-- H(T) = 1 ⟺ T ≅ transitive. -/
theorem H_eq_one_iff_transitive_audit {n : ℕ} (T : Tournament n) (hn : 1 ≤ n) :
    H T = 1 ↔ T ≅ transitiveTournament n :=
  H_eq_one_iff_transitive T hn
#print axioms H_eq_one_iff_transitive_audit

theorem H_transitive_3_eq_one_audit : H (transitiveTournament 3) = 1 :=
  H_transitive_3_eq_one
#print axioms H_transitive_3_eq_one_audit

theorem H_transitive_4_eq_one_audit : H (transitiveTournament 4) = 1 :=
  H_transitive_4_eq_one
#print axioms H_transitive_4_eq_one_audit

theorem H_transitive_5_eq_one_audit : H (transitiveTournament 5) = 1 :=
  H_transitive_5_eq_one
#print axioms H_transitive_5_eq_one_audit

theorem H_threeCycle_eq_three_audit : H threeCycle = 3 :=
  H_threeCycle_eq_three
#print axioms H_threeCycle_eq_three_audit

theorem threeCycle_alpha1_eq_one_audit : alphaCount 1 threeCycle = 1 :=
  threeCycle_alpha1_eq_one
#print axioms threeCycle_alpha1_eq_one_audit

/-! ### Regular tournaments require n to be odd -/

/-- For a regular tournament on n ≥ 1 vertices, n must be odd. PROVED. -/
theorem regular_implies_n_odd_audit {n : ℕ} (T : Tournament n) (hn : 1 ≤ n)
    (hreg : IsRegular T) : Odd n :=
  regular_implies_n_odd T hn hreg
#print axioms regular_implies_n_odd_audit

/-- Regular HasBasePath ⟹ n odd and n ≥ 3. -/
theorem regular_basepath_n_odd_ge_three_audit {n : ℕ} (T : Tournament n)
    (hbp : HasBasePath T) (hn : 2 ≤ n) (hreg : IsRegular T) :
    Odd n ∧ 3 ≤ n :=
  regular_basepath_n_odd_ge_three T hbp hn hreg
#print axioms regular_basepath_n_odd_ge_three_audit

/-- Regular HasBasePath: n = 2k + 3 for some k ≥ 0. -/
theorem regular_basepath_n_in_odd_ge_three_audit {n : ℕ} (T : Tournament n)
    (hbp : HasBasePath T) (hn : 2 ≤ n) (hreg : IsRegular T) :
    ∃ k, n = 2 * k + 3 :=
  regular_basepath_n_in_odd_ge_three T hbp hn hreg
#print axioms regular_basepath_n_in_odd_ge_three_audit

/-! ### Paley(3) ≅ threeCycle -/

theorem paley_3_H_eq_three_audit (P : PaleyType 3) (h : P.T ≅ threeCycle) :
    H P.T = 3 :=
  paley_3_H_eq_three P h
#print axioms paley_3_H_eq_three_audit

/-! ### H(transitive) ≥ 1 for any n -/

theorem H_transitive_ge_one_audit (n : ℕ) (hn : 1 ≤ n) :
    1 ≤ H (transitiveTournament n) :=
  H_transitive_ge_one n hn
#print axioms H_transitive_ge_one_audit

/-- H(transitive_n) = 1 for all n ≥ 1 — PROVED via OCF + transitive_alphaCount_zero. -/
theorem H_transitive_eq_one_general_audit (n : ℕ) (hn : 1 ≤ n) :
    H (transitiveTournament n) = 1 :=
  H_transitive_eq_one_from_ocf n hn
#print axioms H_transitive_eq_one_general_audit
