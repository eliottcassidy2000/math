import TournamentH7.LRCB5PairGridBridge
import TournamentH7.LRCRationalOpenComb

/-!
Exact construction of the strict-open pair ledgers used by the LRC(14)
pair-grid discrepancy bridge.  The normalized rational pair comb supplies one
common component list for exact discrete counts and exact circle volume, with
the always-dangerous point `p = 0` recorded as a separate atom.
-/

namespace LonelyRunner
namespace LRCOpenPairLedger

open Finset RatIntervals LRCStrictInterMerge LRCRationalOpenComb
open LRCOpenDangerComb LRCB5PairGridBridge

noncomputable section

noncomputable def rationalRegionGridCount (q : ℕ) (region : Region) : ℕ := by
  classical
  exact ((Finset.Ico (0 : ℤ) (q : ℤ)).filter fun point : ℤ =>
    strictMem ((point : ℚ) / q) region).card

noncomputable def rationalIntervalGridCount (q : ℕ) (interval : ℚ × ℚ) : ℕ := by
  classical
  exact ((Finset.Ico (0 : ℤ) (q : ℤ)).filter fun point : ℤ =>
    inside ((point : ℚ) / q) interval).card

noncomputable def naturalRegionGridCount (q : ℕ) (region : Region) : ℕ := by
  classical
  exact ((Finset.range q).filter fun point : ℕ =>
    strictMem ((point : ℚ) / q) region).card

theorem naturalRegionGridCount_eq_rationalRegionGridCount
    (q : ℕ) (region : Region) :
    naturalRegionGridCount q region = rationalRegionGridCount q region := by
  classical
  unfold naturalRegionGridCount rationalRegionGridCount
  apply Finset.card_bij (fun (point : ℕ) _ => (point : ℤ))
  · intro point hpoint
    simp only [Finset.mem_filter, Finset.mem_range] at hpoint
    rw [Finset.mem_filter, Finset.mem_Ico]
    refine ⟨⟨by exact_mod_cast (Nat.zero_le point), by exact_mod_cast hpoint.1⟩, ?_⟩
    simpa using hpoint.2
  · intro first hfirst second hsecond heq
    exact_mod_cast heq
  · intro point hpoint
    simp only [Finset.mem_filter, Finset.mem_Ico] at hpoint
    have hpointNonneg : 0 ≤ point := hpoint.1.1
    refine ⟨point.toNat, ?_, ?_⟩
    · rw [Finset.mem_filter, Finset.mem_range]
      constructor
      · have hcast := Int.toNat_of_nonneg hpointNonneg
        omega
      · have hcast : (point.toNat : ℚ) = (point : ℚ) := by
          exact_mod_cast Int.toNat_of_nonneg hpointNonneg
        rw [hcast]
        exact hpoint.2
    · exact Int.toNat_of_nonneg hpointNonneg

theorem singleIntervalGridCount_eq
    (q : ℕ) (interval : ℚ × ℚ) :
    rationalIntervalGridCount q interval =
      openIntervalGridCount q (interval.1 : ℝ) (interval.2 : ℝ) := by
  classical
  unfold rationalIntervalGridCount openIntervalGridCount
  apply congrArg Finset.card
  ext point
  simp only [Finset.mem_filter, Finset.mem_Ico, inside]
  constructor
  · rintro ⟨hpoint, hleft, hright⟩
    refine ⟨hpoint, ?_, ?_⟩
    · have hcast : ((interval.1 : ℚ) : ℝ) <
          ((((point : ℚ) / (q : ℚ)) : ℚ) : ℝ) := by exact_mod_cast hleft
      simpa using hcast
    · have hcast : ((((point : ℚ) / (q : ℚ)) : ℚ) : ℝ) <
          ((interval.2 : ℚ) : ℝ) := by exact_mod_cast hright
      simpa using hcast
  · rintro ⟨hpoint, hleft, hright⟩
    refine ⟨hpoint, ?_, ?_⟩
    · have hcast : ((interval.1 : ℚ) : ℝ) <
          ((((point : ℚ) / (q : ℚ)) : ℚ) : ℝ) := by simpa using hleft
      exact_mod_cast hcast
    · have hcast : ((((point : ℚ) / (q : ℚ)) : ℚ) : ℝ) <
          ((interval.2 : ℚ) : ℝ) := by simpa using hright
      exact_mod_cast hcast

theorem rationalRegionGridCount_cons
    (q : ℕ) (head : ℚ × ℚ) (tail : Region)
    (hnorm : Norm (head :: tail)) :
    rationalRegionGridCount q (head :: tail) =
      openIntervalGridCount q (head.1 : ℝ) (head.2 : ℝ) +
        rationalRegionGridCount q tail := by
  classical
  let grid := Finset.Ico (0 : ℤ) (q : ℤ)
  let headGrid := grid.filter fun point : ℤ => inside ((point : ℚ) / q) head
  let tailGrid := grid.filter fun point : ℤ =>
    strictMem ((point : ℚ) / q) tail
  have hunion :
      grid.filter (fun point : ℤ =>
        strictMem ((point : ℚ) / q) (head :: tail)) =
        headGrid ∪ tailGrid := by
    ext point
    simp only [Finset.mem_filter, strictMem_cons_iff, Finset.mem_union,
      headGrid, tailGrid]
    tauto
  have hdisjoint : Disjoint headGrid tailGrid := by
    rw [Finset.disjoint_left]
    intro point hhead htail
    simp only [headGrid, tailGrid, Finset.mem_filter] at hhead htail
    obtain ⟨other, hotherMem, hother⟩ := htail.2
    have horder := norm_head_le hnorm other hotherMem
    exact (not_lt_of_ge horder) (hother.1.trans hhead.2.2)
  unfold rationalRegionGridCount
  rw [hunion, Finset.card_union_of_disjoint hdisjoint]
  change rationalIntervalGridCount q head + rationalRegionGridCount q tail = _
  rw [singleIntervalGridCount_eq]
  rfl

theorem rationalRegionGridCount_eq_sum
    (q : ℕ) (region : Region) (hnorm : Norm region) :
    rationalRegionGridCount q region =
      (region.map fun interval =>
        openIntervalGridCount q (interval.1 : ℝ) (interval.2 : ℝ)).sum := by
  classical
  induction region with
  | nil => simp [rationalRegionGridCount, strictMem]
  | cons head tail ih =>
      rw [rationalRegionGridCount_cons q head tail hnorm]
      simp only [List.map_cons, List.sum_cons]
      rw [ih (norm_tail hnorm)]

theorem naturalRegionGridCount_eq_sum
    (q : ℕ) (region : Region) (hnorm : Norm region) :
    naturalRegionGridCount q region =
      (region.map fun interval =>
        openIntervalGridCount q (interval.1 : ℝ) (interval.2 : ℝ)).sum := by
  rw [naturalRegionGridCount_eq_rationalRegionGridCount]
  exact rationalRegionGridCount_eq_sum q region hnorm

theorem ratOpenCombTail_mem_bounds
    (w start count : ℕ) (interval : ℚ × ℚ)
    (hinterval : interval ∈ ratOpenCombTail w start count) :
    0 ≤ interval.1 ∧ interval.2 ≤ 1 := by
  induction count generalizing start with
  | zero => simp [ratOpenCombTail] at hinterval
  | succ count ih =>
      simp only [ratOpenCombTail, List.mem_cons] at hinterval
      rcases hinterval with rfl | hinterval
      · simp [ratOpenCombInterval]
      · exact ih (start + 1) hinterval

theorem ratOpenCombRegion_mem_bounds
    (w : ℕ) (interval : ℚ × ℚ)
    (hinterval : interval ∈ ratOpenCombRegion w) :
    0 ≤ interval.1 ∧ interval.2 ≤ 1 := by
  exact ratOpenCombTail_mem_bounds w 0 (w + 1) interval hinterval

theorem ratOpenPairRegion_mem_bounds
    (w₁ w₂ : ℕ) (interval : ℚ × ℚ)
    (hinterval : interval ∈ ratOpenPairRegion w₁ w₂) :
    0 ≤ interval.1 ∧ interval.2 ≤ 1 := by
  obtain ⟨leftSource, hleftSource, rightSource, hrightSource, rfl⟩ :=
    mem_strictInterMerge_sources
      (ratOpenCombRegion w₁) (ratOpenCombRegion w₂) hinterval
  have hleftBounds := ratOpenCombRegion_mem_bounds w₁ leftSource hleftSource
  have hrightBounds := ratOpenCombRegion_mem_bounds w₂ rightSource hrightSource
  unfold RatIntervals.clip
  exact ⟨hleftBounds.1.trans (le_max_left _ _),
    (min_le_left _ _).trans hleftBounds.2⟩

theorem norm_mem_lt
    (region : Region) (hnorm : Norm region)
    (interval : ℚ × ℚ) (hinterval : interval ∈ region) :
    interval.1 < interval.2 := by
  induction region with
  | nil => simp at hinterval
  | cons head tail ih =>
      rcases List.mem_cons.mp hinterval with rfl | hinterval
      · exact norm_head_lt hnorm
      · exact ih (norm_tail hnorm) hinterval

theorem length_eq_map_sub_sum
    (region : Region) (hnorm : Norm region) :
    RatIntervals.length region =
      (region.map fun interval => interval.2 - interval.1).sum := by
  induction region with
  | nil => simp [RatIntervals.length]
  | cons head tail ih =>
      have hhead := norm_head_lt hnorm
      unfold RatIntervals.length
      simp only [List.map_cons, List.sum_cons]
      rw [max_eq_right (sub_nonneg.mpr hhead.le)]
      have htail := ih (norm_tail hnorm)
      unfold RatIntervals.length at htail
      rw [htail]

theorem fin_sum_get_eq_list_map_sum
    {M : Type*} [AddCommMonoid M]
    (region : Region) (f : (ℚ × ℚ) → M) :
    (∑ component : Fin region.length, f (region.get component)) =
      (region.map f).sum := by
  calc
    (∑ component : Fin region.length, f (region.get component)) =
        (List.ofFn fun component => f (region.get component)).sum :=
      List.sum_ofFn.symm
    _ = (List.map f (List.ofFn region.get)).sum := by
      rw [List.map_ofFn]
      rfl
    _ = (region.map f).sum := by rw [List.ofFn_get]

theorem cast_length_eq_fin_sum_sub
    (region : Region) (hnorm : Norm region) :
    ((RatIntervals.length region : ℚ) : ℝ) =
      ∑ component : Fin region.length,
        ((((region.get component).2 : ℚ) : ℝ) -
          (((region.get component).1 : ℚ) : ℝ)) := by
  rw [length_eq_map_sub_sum region hnorm]
  push_cast
  rw [List.map_map]
  rw [← fin_sum_get_eq_list_map_sum region
    ((Rat.cast : ℚ → ℝ) ∘ (fun interval => interval.2 - interval.1))]
  simp only [Function.comp_apply]
  push_cast
  rfl

theorem not_strictMem_zero_ratOpenPairRegion (w₁ w₂ : ℕ) :
    ¬ strictMem 0 (ratOpenPairRegion w₁ w₂) := by
  rintro ⟨interval, hinterval, hinside⟩
  have hbounds := ratOpenPairRegion_mem_bounds w₁ w₂ interval hinterval
  exact (not_lt_of_ge hbounds.1) hinside.1

theorem strictMem_ratOpenPairRegion_natAbs_iff_dangers
    (firstSpeed secondSpeed : ℤ)
    (hfirst : firstSpeed ≠ 0) (hsecond : secondSpeed ≠ 0)
    (x : ℚ) (hx : (x : ℝ) ∈ Set.Ioo (0 : ℝ) 1) :
    strictMem x
        (ratOpenPairRegion firstSpeed.natAbs secondSpeed.natAbs) ↔
      (((x : ℚ) : ℝ) : UnitAddCircle) ∈
          LRCCommensuration.danger firstSpeed 0 (1 / 14) ∧
        (((x : ℚ) : ℝ) : UnitAddCircle) ∈
          LRCCommensuration.danger secondSpeed 0 (1 / 14) := by
  rw [strictMem_ratOpenPairRegion_iff_mem_real _ _
    (Int.natAbs_pos.mpr hfirst) (Int.natAbs_pos.mpr hsecond)]
  rw [Set.mem_inter_iff]
  rw [mem_openCombRegion_natAbs_iff_mem_danger firstSpeed hfirst (x : ℝ) hx]
  rw [mem_openCombRegion_natAbs_iff_mem_danger secondSpeed hsecond (x : ℝ) hx]

theorem strictMem_ratOpenPairRegion_iff_pairDanger
    (v : Fin 13 → ℤ) (first second : Fin 13)
    (hfirst : v first ≠ 0) (hsecond : v second ≠ 0)
    (x : ℚ) (hx : (x : ℝ) ∈ Set.Ioo (0 : ℝ) 1) :
    strictMem x (ratOpenPairRegion (v first).natAbs (v second).natAbs) ↔
      (((x : ℚ) : ℝ) : UnitAddCircle) ∈ pairDanger v {first, second} := by
  rw [strictMem_ratOpenPairRegion_natAbs_iff_dangers
    (v first) (v second) hfirst hsecond x hx]
  unfold pairDanger
  constructor
  · rintro ⟨hfirstDanger, hsecondDanger⟩ index hindex
    simp only [Finset.mem_insert, Finset.mem_singleton] at hindex
    rcases hindex with rfl | rfl
    · exact hfirstDanger
    · exact hsecondDanger
  · intro hall
    exact ⟨hall first (by simp), hall second (by simp)⟩

theorem strictMem_grid_iff_pairDanger
    (v : Fin 13 → ℤ) (first second : Fin 13)
    (hfirst : v first ≠ 0) (hsecond : v second ≠ 0)
    (q p : ℕ) (hq : 0 < q) (hp : 0 < p) (hpq : p < q) :
    strictMem ((p : ℚ) / q)
        (ratOpenPairRegion (v first).natAbs (v second).natAbs) ↔
      (((p : ℝ) / q : ℝ) : UnitAddCircle) ∈ pairDanger v {first, second} := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hpR : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  have hpqR : (p : ℝ) < (q : ℝ) := by exact_mod_cast hpq
  have hx : ((((p : ℚ) / (q : ℚ) : ℚ) : ℝ)) ∈ Set.Ioo (0 : ℝ) 1 := by
    constructor <;> push_cast
    · exact div_pos hpR hqR
    · exact (div_lt_one hqR).2 hpqR
  simpa using strictMem_ratOpenPairRegion_iff_pairDanger
    v first second hfirst hsecond ((p : ℚ) / q) hx

theorem zero_mem_pairDanger
    (v : Fin 13 → ℤ) (q : ℕ) (T : Finset (Fin 13)) (hq : 0 < q) :
    (((0 : ℝ) / q : ℝ) : UnitAddCircle) ∈ pairDanger v T := by
  unfold pairDanger
  intro index hindex
  simpa using (grid_mem_danger_iff_not_inBand v q 0 index hq).2 (by
      unfold LRC14Concrete.inBand
      simp
      omega)

theorem pairDangerGridCount_pair_eq
    (v : Fin 13 → ℤ) (first second : Fin 13)
    (hfirst : v first ≠ 0) (hsecond : v second ≠ 0)
    (q : ℕ) (hq : 0 < q) :
    pairDangerGridCount v q {first, second} =
      1 + naturalRegionGridCount q
        (ratOpenPairRegion (v first).natAbs (v second).natAbs) := by
  classical
  let region := ratOpenPairRegion (v first).natAbs (v second).natAbs
  have hfilter :
      (Finset.range q).filter (fun p : ℕ =>
        (((p : ℝ) / q : ℝ) : UnitAddCircle) ∈ pairDanger v {first, second}) =
        insert 0 ((Finset.range q).filter fun p : ℕ =>
          strictMem ((p : ℚ) / q) region) := by
    ext p
    simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_insert]
    constructor
    · rintro ⟨hpq, hdanger⟩
      rcases Nat.eq_zero_or_pos p with rfl | hp
      · exact Or.inl rfl
      · exact Or.inr ⟨hpq,
          (strictMem_grid_iff_pairDanger
            v first second hfirst hsecond q p hq hp hpq).2 hdanger⟩
    · rintro (rfl | ⟨hpq, hregion⟩)
      · refine ⟨hq, ?_⟩
        simpa using zero_mem_pairDanger v q {first, second} hq
      · have hp : 0 < p := by
          by_contra hp0
          have : p = 0 := by omega
          subst p
          exact not_strictMem_zero_ratOpenPairRegion
            (v first).natAbs (v second).natAbs (by simpa [region] using hregion)
        exact ⟨hpq, (strictMem_grid_iff_pairDanger
          v first second hfirst hsecond q p hq hp hpq).1 hregion⟩
  unfold pairDangerGridCount naturalRegionGridCount
  rw [hfilter, Finset.card_insert_of_notMem]
  · simp [region, Nat.add_comm]
  · simp only [Finset.mem_filter, Finset.mem_range, not_and]
    intro _
    simpa [region] using not_strictMem_zero_ratOpenPairRegion
      (v first).natAbs (v second).natAbs

theorem pairDanger_pair_eq_positivePairDanger
    (v : Fin 13 → ℤ) (first second : Fin 13) :
    pairDanger v {first, second} =
      positivePairDanger (v first).natAbs (v second).natAbs := by
  ext x
  unfold pairDanger positivePairDanger
  simp only [Set.mem_setOf_eq, Set.mem_inter_iff]
  constructor
  · intro hall
    exact ⟨by
      rw [danger_natAbs]
      exact hall first (by simp), by
      rw [danger_natAbs]
      exact hall second (by simp)⟩
  · rintro ⟨hfirst, hsecond⟩ index hindex
    simp only [Finset.mem_insert, Finset.mem_singleton] at hindex
    rcases hindex with rfl | rfl
    · rwa [danger_natAbs] at hfirst
    · rwa [danger_natAbs] at hsecond

theorem pairContinuumCorrelation_pair_eq_length
    (v : Fin 13 → ℤ) (first second : Fin 13)
    (hfirst : v first ≠ 0) (hsecond : v second ≠ 0) :
    pairContinuumCorrelation v {first, second} =
      (((RatIntervals.length
        (ratOpenPairRegion (v first).natAbs (v second).natAbs) : ℚ) : ℝ)) := by
  have hwfirst : 0 < (v first).natAbs := Int.natAbs_pos.mpr hfirst
  have hwsecond : 0 < (v second).natAbs := Int.natAbs_pos.mpr hsecond
  unfold pairContinuumCorrelation
  rw [pairDanger_pair_eq_positivePairDanger]
  rw [volume_positivePairDanger_eq_length _ _ hwfirst hwsecond]
  rw [ENNReal.toReal_ofReal]
  exact_mod_cast (RatIntervals.length_nonneg
    (ratOpenPairRegion (v first).natAbs (v second).natAbs))

theorem openPairIntervalLedger_pair
    (v : Fin 13 → ℤ) (first second : Fin 13) (hneq : first ≠ second)
    (hfirst : v first ≠ 0) (hsecond : v second ≠ 0)
    (q : ℕ) (hq : 0 < q) :
    Nonempty (OpenPairIntervalLedger v q {first, second}) := by
  let region := ratOpenPairRegion (v first).natAbs (v second).natAbs
  have hnorm : Norm region := by
    exact norm_ratOpenPairRegion _ _
      (Int.natAbs_pos.mpr hfirst) (Int.natAbs_pos.mpr hsecond)
  refine ⟨{
    componentCount := region.length
    left := fun component => ((region.get component).1 : ℝ)
    right := fun component => ((region.get component).2 : ℝ)
    left_nonneg := ?_
    right_le_one := ?_
    left_le_right := ?_
    componentCount_succ_le := ?_
    allGrid_count_eq := ?_
    continuum_eq := ?_
  }⟩
  · intro component
    have hbounds := ratOpenPairRegion_mem_bounds
      (v first).natAbs (v second).natAbs (region.get component) (by
        simp [region])
    exact_mod_cast hbounds.1
  · intro component
    have hbounds := ratOpenPairRegion_mem_bounds
      (v first).natAbs (v second).natAbs (region.get component) (by
        simp [region])
    exact_mod_cast hbounds.2
  · intro component
    have hlive := norm_mem_lt region hnorm (region.get component)
      (region.get_mem component)
    exact_mod_cast hlive.le
  · have hlength := length_ratOpenPairRegion_le
      (v first).natAbs (v second).natAbs
    have hpositiveFirst : 0 < (v first).natAbs := Int.natAbs_pos.mpr hfirst
    have hpositiveSecond : 0 < (v second).natAbs := Int.natAbs_pos.mpr hsecond
    have hbudget : region.length + 1 ≤
        2 * ((v first).natAbs + (v second).natAbs) := by
      dsimp [region]
      omega
    simpa [hneq] using hbudget
  · rw [allGridJointFail_eq_pairDanger_gridCount v q {first, second} hq]
    rw [pairDangerGridCount_pair_eq v first second hfirst hsecond q hq]
    rw [naturalRegionGridCount_eq_sum q region hnorm]
    rw [← fin_sum_get_eq_list_map_sum region
      (fun interval => openIntervalGridCount q (interval.1 : ℝ) (interval.2 : ℝ))]
  · rw [pairContinuumCorrelation_pair_eq_length v first second hfirst hsecond]
    exact cast_length_eq_fin_sum_sub region hnorm

/-- Every nonzero speed family has the strict-open pair ledgers required by
the complete-grid discrepancy theorem, at every positive modulus. -/
theorem openPairIntervalLedgersAt_of_nonzero
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (q : ℕ) (hq : 0 < q) :
    OpenPairIntervalLedgersAt v q := by
  intro support hsupport
  have hcard : support.card = 2 := (Finset.mem_powersetCard.mp hsupport).2
  obtain ⟨first, second, hneq, rfl⟩ := Finset.card_eq_two.mp hcard
  exact openPairIntervalLedger_pair v first second hneq
    (hv first) (hv second) q hq

theorem rawPairGridDiscrepancyAt_of_nonzero
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (q : ℕ) (hq : 0 < q) :
    RawPairGridDiscrepancyAt v q :=
  rawPairGridDiscrepancy_of_openIntervalLedgers v q hq
    (openPairIntervalLedgersAt_of_nonzero v hv q hq)

/-- The strict pair-grid socket no longer has a geometric premise: the only
remaining producer is the continuum floor. -/
theorem normalizedMass2_clean_534_gt_neg_target_of_continuum_floor'
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (hcontinuum :
      -(LRCPairRatioLayerArithmetic.negativePairTierBoundPathOnly : ℝ) ≤
        continuumMass2 v) :
    -(13 / 30 : ℚ) <
      LRCB5NormalizedBridge.normalizedMass2 v
        (LRCB5CleanModulus.cleanModulus v 534) := by
  apply normalizedMass2_clean_534_gt_neg_target_of_continuum_floor
    v hv hcontinuum
  exact rawPairGridDiscrepancyAt_of_nonzero v hv _
    (Nat.zero_lt_of_lt (LRCB5CleanModulus.one_lt_cleanModulus v 534 hv))

#print axioms rationalRegionGridCount_eq_sum
#print axioms naturalRegionGridCount_eq_sum
#print axioms ratOpenPairRegion_mem_bounds
#print axioms pairDangerGridCount_pair_eq
#print axioms pairContinuumCorrelation_pair_eq_length
#print axioms openPairIntervalLedger_pair
#print axioms openPairIntervalLedgersAt_of_nonzero
#print axioms rawPairGridDiscrepancyAt_of_nonzero
#print axioms normalizedMass2_clean_534_gt_neg_target_of_continuum_floor'

end
end LRCOpenPairLedger
end LonelyRunner
