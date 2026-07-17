/-
  TournamentH7.LRCRationalRegionProvider

  Concrete rational-interval analytic provider for THM-933. A normalized
  `RatIntervals.Region` is cast to a finite union of real half-open intervals,
  then periodized through `Int.fract`. The resulting survivor is measurable,
  its indicator is one-periodic, and its mass on one period is exactly the cast
  rational region length. These facts instantiate the centered-primitive and
  eta/q duality theorems from `LRCLocalDensityBlockGluing`. The final section
  specializes this provider to the concrete `rationalCircleSurvivor` and
  `rationalCircleChart`; their normalization and unit-chart premises are now
  unconditional. Exact deletion density is reduced further to a finite sum of
  rational chart-tooth clip lengths, so only evaluation of those explicit
  rational sums remains external.

  The exact-measure bridge needs precisely the rational measure discipline:
  `Norm region` for live, ordered, disjoint intervals, plus `regionInUnit region`
  to identify the listed chart with the fundamental period. No component-count
  topology is used.

  For the concrete moving circle charts, the unconditional normalization
  theorem discharges both measure premises. In particular, neither
  `BoundaryFaithfulRotation` nor cut-balance is an analytic assumption; those
  invariants remain confined to component-count transport. Zarankiewicz
  relation constraints can therefore narrow the later numeric region census
  without entering this provider interface.

  Assumption challenge: the periodicization quotient uses neither runners nor
  arcs as tournament vertices. It preserves survivor membership modulo one and
  exact period mass, while forgetting which integer translate represented a
  point. Tournament Analysis is therefore not a natural pairwise interface for
  this measure-provider lemma.

  No native_decide and no sorry.
-/

import Mathlib.MeasureTheory.Function.Floor
import TournamentH7.LRCLocalDensityBlockGluing

namespace LRC14.LocalDensityBlockGluing

open scoped BigOperators
open MeasureTheory Set
open LonelyRunner.RatIntervals

/-- Finite real union induced by a rational interval region. Each rational
half-open interval `[a,b)` is cast to the real interval `Ico a b`. -/
def rationalRegionSet : Region → Set ℝ
  | [] => ∅
  | interval :: region =>
      Ico (interval.1 : ℝ) (interval.2 : ℝ) ∪ rationalRegionSet region

/-- One-periodic real survivor induced by a rational region, obtained by
testing the fractional part against its fundamental-domain realization. -/
def periodicRationalRegion (region : Region) : Set ℝ :=
  Int.fract ⁻¹' rationalRegionSet region

/-- Every listed rational interval lies in the fundamental chart `[0,1]`. -/
def regionInUnit (region : Region) : Prop :=
  ∀ interval ∈ region, 0 ≤ interval.1 ∧ interval.2 ≤ 1

/-- Real density supplied by a rational region: the cast of its exact rational
length ledger. -/
noncomputable def rationalRegionDensity (region : Region) : ℝ :=
  (length region : ℝ)

/-- Membership in the cast real union is membership in one listed cast
half-open interval. -/
theorem mem_rationalRegionSet {region : Region} {x : ℝ} :
    x ∈ rationalRegionSet region ↔
      ∃ interval ∈ region, (interval.1 : ℝ) ≤ x ∧ x < (interval.2 : ℝ) := by
  induction region with
  | nil => simp [rationalRegionSet]
  | cons interval region inductionHypothesis =>
      simp [rationalRegionSet, inductionHypothesis]

/-- On rational points, the cast real union has exactly the original
`RatIntervals.mem` semantics. -/
theorem ratCast_mem_rationalRegionSet {region : Region} {x : ℚ} :
    (x : ℝ) ∈ rationalRegionSet region ↔
      LonelyRunner.RatIntervals.mem x region := by
  rw [mem_rationalRegionSet]
  unfold LonelyRunner.RatIntervals.mem
  constructor
  · rintro ⟨interval, hinterval, hleft, hright⟩
    exact ⟨interval, hinterval, by exact_mod_cast hleft, by exact_mod_cast hright⟩
  · rintro ⟨interval, hinterval, hleft, hright⟩
    exact ⟨interval, hinterval, by exact_mod_cast hleft, by exact_mod_cast hright⟩

/-- A finite union of cast rational intervals is measurable. -/
theorem measurableSet_rationalRegionSet (region : Region) :
    MeasurableSet (rationalRegionSet region) := by
  induction region with
  | nil => simp [rationalRegionSet]
  | cons interval region inductionHypothesis =>
      exact measurableSet_Ico.union inductionHypothesis

/-- The fractional-part preimage of a rational real region is measurable. -/
theorem measurableSet_periodicRationalRegion (region : Region) :
    MeasurableSet (periodicRationalRegion region) :=
  (measurableSet_rationalRegionSet region).preimage measurable_fract

/-- A unit-chart rational region casts into the real fundamental interval
`[0,1)`. -/
theorem rationalRegionSet_subset_Ico
    (region : Region) (hunit : regionInUnit region) :
    rationalRegionSet region ⊆ Ico (0 : ℝ) 1 := by
  intro x hx
  obtain ⟨interval, hinterval, hleft, hright⟩ := mem_rationalRegionSet.mp hx
  have hbounds := hunit interval hinterval
  constructor
  · exact le_trans (by exact_mod_cast hbounds.1) hleft
  · exact lt_of_lt_of_le hright (by exact_mod_cast hbounds.2)

/-- Under `Norm`, the first interval is disjoint from the cast union of the
tail. This is the sole disjointness input to the volume induction. -/
theorem head_disjoint_rationalRegionSet
    {interval : ℚ × ℚ} {region : Region} (hnorm : Norm (interval :: region)) :
    Disjoint (Ico (interval.1 : ℝ) (interval.2 : ℝ)) (rationalRegionSet region) := by
  apply Set.disjoint_left.mpr
  intro x hxhead hxtail
  obtain ⟨other, hother, hotherLeft, _hotherRight⟩ := mem_rationalRegionSet.mp hxtail
  have hordered : interval.2 ≤ other.1 := norm_head_le hnorm other hother
  have horderedReal : (interval.2 : ℝ) ≤ (other.1 : ℝ) := by exact_mod_cast hordered
  linarith [hxhead.2]

/-- Exact finite-union volume formula for a normalized rational region. -/
theorem volume_rationalRegionSet
    {region : Region} (hnorm : Norm region) :
    volume (rationalRegionSet region) = ENNReal.ofReal (rationalRegionDensity region) := by
  induction region with
  | nil => simp [rationalRegionSet, rationalRegionDensity, length]
  | cons interval region inductionHypothesis =>
      have hinterval : interval.1 < interval.2 := norm_head_lt hnorm
      have htail : Norm region := norm_tail hnorm
      have hintervalNonneg : (0 : ℚ) ≤ interval.2 - interval.1 := sub_nonneg.mpr hinterval.le
      have hdifferenceNonneg : (0 : ℝ) ≤ (interval.2 : ℝ) - (interval.1 : ℝ) := by
        exact_mod_cast hintervalNonneg
      have hdensityNonneg : 0 ≤ rationalRegionDensity region := by
        unfold rationalRegionDensity
        exact_mod_cast length_nonneg region
      rw [rationalRegionSet,
        measure_union (head_disjoint_rationalRegionSet hnorm)
          (measurableSet_rationalRegionSet region),
        Real.volume_Ico, inductionHypothesis htail,
        ← ENNReal.ofReal_add hdifferenceNonneg hdensityNonneg]
      congr 1
      simp only [rationalRegionDensity, length, List.map_cons, List.sum_cons]
      rw [max_eq_right hintervalNonneg]
      push_cast
      ring

/-- The indicator of the fractional-part periodization is one-periodic. -/
theorem periodicRationalRegion_indicator_periodic (region : Region) :
    Function.Periodic
      (fun t => (periodicRationalRegion region).indicator (fun _ => (1 : ℝ)) t) 1 := by
  intro point
  change (periodicRationalRegion region).indicator (fun _ => (1 : ℝ)) (point + 1) =
    (periodicRationalRegion region).indicator (fun _ => (1 : ℝ)) point
  have hmem : point + 1 ∈ periodicRationalRegion region ↔
      point ∈ periodicRationalRegion region := by
    simp [periodicRationalRegion, Int.fract_add_one]
  by_cases hx : point ∈ periodicRationalRegion region
  · rw [Set.indicator_of_mem hx, Set.indicator_of_mem (hmem.mpr hx)]
  · rw [Set.indicator_of_notMem hx,
      Set.indicator_of_notMem (fun h => hx (hmem.mp h))]

/-- Exact THM-933 provider density: on one real period, the periodized survivor
has Lebesgue measure equal to the cast rational region length. Half-open versus
half-closed endpoint conventions agree because endpoints are null. -/
theorem periodicRationalRegion_period_density
    (region : Region) (hnorm : Norm region) (hunit : regionInUnit region) :
    (volume (periodicRationalRegion region ∩ Ioc (0 : ℝ) 1)).toReal =
      rationalRegionDensity region := by
  have hinside : periodicRationalRegion region ∩ Ioo (0 : ℝ) 1 =
      rationalRegionSet region ∩ Ioo (0 : ℝ) 1 := by
    ext x
    simp only [periodicRationalRegion, Set.mem_inter_iff, Set.mem_preimage, Set.mem_Ioo]
    by_cases hx : 0 < x ∧ x < 1
    · have hfract : Int.fract x = x := Int.fract_eq_self.mpr ⟨hx.1.le, hx.2⟩
      simp [hx, hfract]
    · simp [hx]
  have hIocIoo : volume (periodicRationalRegion region ∩ Ioc (0 : ℝ) 1) =
      volume (periodicRationalRegion region ∩ Ioo (0 : ℝ) 1) :=
    measure_congr (ae_eq_set_inter (ae_eq_refl _) Ioo_ae_eq_Ioc.symm)
  have hbase : volume (rationalRegionSet region ∩ Ioo (0 : ℝ) 1) =
      volume (rationalRegionSet region) := by
    calc
      volume (rationalRegionSet region ∩ Ioo (0 : ℝ) 1) =
          volume (rationalRegionSet region ∩ Ico (0 : ℝ) 1) :=
        measure_congr (ae_eq_set_inter (ae_eq_refl _) Ioo_ae_eq_Ico)
      _ = volume (rationalRegionSet region) := by
        rw [inter_eq_left.mpr (rationalRegionSet_subset_Ico region hunit)]
  have hdensityNonneg : 0 ≤ rationalRegionDensity region := by
    unfold rationalRegionDensity
    exact_mod_cast length_nonneg region
  rw [hIocIoo, hinside, hbase, volume_rationalRegionSet hnorm,
    ENNReal.toReal_ofReal hdensityNonneg]

/-- Bundled analytic-provider hypotheses consumed by THM-933: measurability,
one-periodic indicator, and exact one-period density. -/
theorem rationalRegion_analytic_provider
    (region : Region) (hnorm : Norm region) (hunit : regionInUnit region) :
    MeasurableSet (periodicRationalRegion region) ∧
      Function.Periodic
        (fun t => (periodicRationalRegion region).indicator (fun _ => (1 : ℝ)) t) 1 ∧
      (volume (periodicRationalRegion region ∩ Ioc (0 : ℝ) 1)).toReal =
        rationalRegionDensity region :=
  ⟨measurableSet_periodicRationalRegion region,
    periodicRationalRegion_indicator_periodic region,
    periodicRationalRegion_period_density region hnorm hunit⟩

/-- Concrete centered-primitive identity for a rational-region survivor. -/
theorem rationalRegion_centeredPrimitive_interval_identity
    (region : Region) {a b : ℝ} (hab : a ≤ b) :
    (volume (periodicRationalRegion region ∩ Ioc a b)).toReal -
        rationalRegionDensity region * (b - a) =
      centeredPrimitive (periodicRationalRegion region) (rationalRegionDensity region) b -
        centeredPrimitive (periodicRationalRegion region) (rationalRegionDensity region) a :=
  centeredPrimitive_interval_identity
    (periodicRationalRegion region) (measurableSet_periodicRationalRegion region)
    (rationalRegionDensity region) hab

/-- Fixed-scale `eta` for a rational-region survivor is attained on `[0,1]`. -/
theorem rationalRegion_fixedScaleEta_attained
    (region : Region) {ell : ℝ} (hell : 0 < ell) :
    ∃ startPoint ∈ Icc (0 : ℝ) 1,
      fixedScaleEta (periodicRationalRegion region) ell =
          retainedWindow (periodicRationalRegion region) ell startPoint / ell ∧
      ∀ point ∈ Icc (0 : ℝ) 1,
        fixedScaleEta (periodicRationalRegion region) ell ≤
          retainedWindow (periodicRationalRegion region) ell point / ell :=
  fixedScaleEta_attained (periodicRationalRegion region)
    (measurableSet_periodicRationalRegion region) hell

/-- The eta/q extremal equality is attained at some scale `0 < ell ≤ 1` for a
normalized rational-region survivor. -/
theorem rationalRegion_eta_q_attained
    (region : Region) (hnorm : Norm region) (hunit : regionInUnit region) :
    ∃ ell ∈ Ioc (0 : ℝ) 1,
      ell * (rationalRegionDensity region -
        fixedScaleEta (periodicRationalRegion region) ell) =
          primitiveDiscrepancy (periodicRationalRegion region)
            (rationalRegionDensity region) :=
  exists_fixedScaleEta_deficit_eq_primitiveDiscrepancy
    (periodicRationalRegion region) (measurableSet_periodicRationalRegion region)
    (rationalRegionDensity region) (periodicRationalRegion_indicator_periodic region)
    (periodicRationalRegion_period_density region hnorm hunit)

/-- Full concrete eta/q duality for a normalized rational-region survivor. -/
theorem rationalRegion_eta_q_duality
    (region : Region) (hnorm : Norm region) (hunit : regionInUnit region) :
    primitiveDiscrepancy (periodicRationalRegion region) (rationalRegionDensity region) =
      fixedScaleDeficitSup (periodicRationalRegion region) (rationalRegionDensity region) :=
  primitiveDiscrepancy_eq_fixedScaleDeficitSup
    (periodicRationalRegion region) (measurableSet_periodicRationalRegion region)
    (rationalRegionDensity region) (periodicRationalRegion_indicator_periodic region)
    (periodicRationalRegion_period_density region hnorm hunit)

/-! ## Concrete rational-circle analytic providers -/

/-- One-periodic real set represented by the recursive rational circle
survivor. -/
def periodicRationalCircleSurvivor
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) : Set ℝ :=
  periodicRationalRegion (rationalCircleSurvivor shift tooth toothCount)

/-- One-periodic real set represented by the normalized rational circle
chart. -/
def periodicRationalCircleChart
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) : Set ℝ :=
  periodicRationalRegion (rationalCircleChart shift tooth toothCount)

/-- Exact real density ledger of the recursive rational circle survivor. -/
noncomputable def rationalCircleSurvivorDensity
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) : ℝ :=
  rationalRegionDensity (rationalCircleSurvivor shift tooth toothCount)

/-- Exact real density ledger of the normalized rational circle chart. -/
noncomputable def rationalCircleChartDensity
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) : ℝ :=
  rationalRegionDensity (rationalCircleChart shift tooth toothCount)

/-- The rational length ledger, hence its real density cast, is invariant
under interval-list permutation. -/
theorem rationalRegionDensity_eq_of_perm {left right : Region}
    (hperm : left.Perm right) :
    rationalRegionDensity left = rationalRegionDensity right := by
  unfold rationalRegionDensity LonelyRunner.RatIntervals.length
  congr 1
  exact (hperm.map (fun interval => max 0 (interval.2 - interval.1))).sum_eq

/-- A sorted rational circle translation preserves the exact density
ledger. -/
theorem rationalRegionDensity_sortedTranslateCirc
    (shift : ℚ) (region : Region) (hunit : RegionInUnit region) :
    rationalRegionDensity (sortedTranslateCirc shift region) =
      rationalRegionDensity region := by
  calc
    rationalRegionDensity (sortedTranslateCirc shift region) =
        rationalRegionDensity (translateCirc shift region) :=
      rationalRegionDensity_eq_of_perm (sortedTranslateCirc_perm shift region)
    _ = rationalRegionDensity region := by
      unfold rationalRegionDensity
      rw [length_translateCirc]
      intro interval hinterval
      have hbounds := hunit interval hinterval
      linarith

/-- Recharting changes coordinates and interval order, but not survivor
density. -/
theorem rationalCircleChartDensity_eq_survivorDensity
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    rationalCircleChartDensity shift tooth toothCount =
      rationalCircleSurvivorDensity shift tooth toothCount := by
  unfold rationalCircleChartDensity rationalCircleSurvivorDensity
  cases toothCount with
  | zero => rfl
  | succ toothCount =>
      simpa [rationalCircleChart, rationalCircleRechart] using
        rationalRegionDensity_sortedTranslateCirc
          (shift (toothCount + 1))
          (rationalCircleSurvivor shift tooth (toothCount + 1))
          (regionInUnit_rationalCircleSurvivor shift tooth (toothCount + 1))

/-! ### Exact anchored-deletion length and density ledger -/

/-- Removing degenerate interval pairs does not change the rational length
ledger. -/
theorem length_filterLive_eq (region : Region) :
    length (region.filter fun interval => decide (interval.1 < interval.2)) =
      length region := by
  induction region with
  | nil => rfl
  | cons interval region inductionHypothesis =>
      by_cases hinterval : interval.1 < interval.2
      · rw [List.filter_cons_of_pos (by simpa using hinterval)]
        unfold length at inductionHypothesis ⊢
        simp only [List.map_cons, List.sum_cons]
        rw [inductionHypothesis]
      · rw [List.filter_cons_of_neg (by simpa using hinterval)]
        unfold length at inductionHypothesis ⊢
        simp only [List.map_cons, List.sum_cons]
        rw [inductionHypothesis, max_eq_left (by
          push Not at hinterval
          linarith)]
        ring

/-- The computational filtered cut has the same exact rational length as the
unfiltered two-piece cut. -/
theorem length_cutF_eq_cut (interval toothInterval : ℚ × ℚ) :
    length (cutF interval toothInterval) =
      length (cut interval toothInterval) := by
  exact length_filterLive_eq _

/-- Exact filtered partition identity. For a live subtrahend, the retained
pieces plus the clipped overlap equal the original rational length. No `Norm`
or containment hypothesis is needed because this is an exact list ledger. -/
theorem length_diff1F_add_inter_exact
    (region : Region) {toothInterval : ℚ × ℚ}
    (htooth : toothInterval.1 ≤ toothInterval.2) :
    length (diff1F region toothInterval) +
        length (inter region [toothInterval]) =
      length region := by
  have hfiltered : length (diff1F region toothInterval) =
      length (diff1 region toothInterval) := by
    induction region with
    | nil => rfl
    | cons interval region inductionHypothesis =>
        unfold diff1F diff1 at inductionHypothesis ⊢
        simp only [List.flatMap_cons]
        rw [length_append, length_append, inductionHypothesis,
          length_cutF_eq_cut]
  rw [hfiltered]
  exact length_diff1_add_inter region htooth

/-- Universal exact anchored-tooth deletion ledger. The actual loss is the
overlap with `[0,tooth.width)`, not necessarily the full tooth width. -/
theorem length_deleteAnchoredTooth_add_overlap
    (region : Region) (tooth : AnchoredCircularTooth) :
    length (deleteAnchoredTooth region tooth) +
        length (inter region [((0 : ℚ), tooth.width)]) =
      length region := by
  exact length_diff1F_add_inter_exact region tooth.width_pos.le

/-- Intersection with one interval is the explicit sum of its per-source
clip lengths. -/
theorem length_inter_singleton_eq_clip_sum
    (region : Region) (interval : ℚ × ℚ) :
    length (inter region [interval]) =
      (region.map fun source =>
        max 0 (min source.2 interval.2 - max source.1 interval.1)).sum := by
  induction region with
  | nil => simp [inter, length]
  | cons source region inductionHypothesis =>
      rw [show inter (source :: region) [interval] =
        [clip source interval] ++ inter region [interval] by rfl]
      rw [length_append, inductionHypothesis]
      simp [length, clip]

/-- Rational amount removed by the next anchored tooth from the current
normalized circle chart. -/
def rationalCircleToothOverlapLength
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) : ℚ :=
  length (inter (rationalCircleChart shift tooth toothCount)
    [((0 : ℚ), (tooth toothCount).width)])

/-- Real cast of the exact rational next-tooth overlap ledger. -/
noncomputable def rationalCircleToothOverlapDensity
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) : ℝ :=
  (rationalCircleToothOverlapLength shift tooth toothCount : ℝ)

/-- Fully explicit rational clip-sum formula for the next-tooth overlap. Unit
chart containment removes the redundant `max interval.1 0` endpoint. -/
theorem rationalCircleToothOverlapLength_eq_clip_sum
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    rationalCircleToothOverlapLength shift tooth toothCount =
      ((rationalCircleChart shift tooth toothCount).map fun interval =>
        max 0 (min interval.2 (tooth toothCount).width - interval.1)).sum := by
  rw [rationalCircleToothOverlapLength, length_inter_singleton_eq_clip_sum]
  congr 1
  apply List.map_congr_left
  intro interval hinterval
  have hnonnegative :=
    (regionInUnit_rationalCircleChart shift tooth toothCount interval hinterval).1
  rw [max_eq_left hnonnegative]

/-- Coordinate recharting preserves the underlying exact rational length, not
only its real density cast. -/
theorem rationalCircleChart_length_eq_survivor_length
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    length (rationalCircleChart shift tooth toothCount) =
      length (rationalCircleSurvivor shift tooth toothCount) := by
  have hdensity :=
    rationalCircleChartDensity_eq_survivorDensity shift tooth toothCount
  unfold rationalCircleChartDensity rationalCircleSurvivorDensity
    rationalRegionDensity at hdensity
  exact_mod_cast hdensity

/-- Exact rational survivor recurrence: the next length plus the actual
chart-tooth overlap equals the current survivor length. -/
theorem rationalCircleSurvivor_length_succ_add_overlap
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    length (rationalCircleSurvivor shift tooth (toothCount + 1)) +
        rationalCircleToothOverlapLength shift tooth toothCount =
      length (rationalCircleSurvivor shift tooth toothCount) := by
  rw [rationalCircleSurvivor_succ]
  unfold rationalCircleToothOverlapLength
  rw [length_deleteAnchoredTooth_add_overlap]
  exact rationalCircleChart_length_eq_survivor_length shift tooth toothCount

/-- Subtraction form of the exact rational survivor recurrence. -/
theorem rationalCircleSurvivor_length_succ_eq_sub_overlap
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    length (rationalCircleSurvivor shift tooth (toothCount + 1)) =
      length (rationalCircleSurvivor shift tooth toothCount) -
        rationalCircleToothOverlapLength shift tooth toothCount := by
  linarith [rationalCircleSurvivor_length_succ_add_overlap
    shift tooth toothCount]

/-- Closed exact rational length formula: density evaluation is reduced to a
finite sum of explicit anchored chart-tooth overlaps. -/
theorem rationalCircleSurvivor_length_eq_one_sub_sum_overlap
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth) :
    ∀ toothCount,
      length (rationalCircleSurvivor shift tooth toothCount) =
        1 - ∑ count ∈ Finset.range toothCount,
          rationalCircleToothOverlapLength shift tooth count := by
  intro toothCount
  induction toothCount with
  | zero => simp [rationalCircleSurvivor, length]
  | succ toothCount inductionHypothesis =>
      rw [rationalCircleSurvivor_length_succ_eq_sub_overlap,
        inductionHypothesis, Finset.sum_range_succ]
      ring

/-- The next-tooth overlap is always nonnegative. -/
theorem rationalCircleToothOverlapLength_nonneg
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    0 ≤ rationalCircleToothOverlapLength shift tooth toothCount := by
  exact length_nonneg _

/-- A normalized chart can overlap an anchored tooth by at most the tooth's
full width. -/
theorem rationalCircleToothOverlapLength_le_width
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    rationalCircleToothOverlapLength shift tooth toothCount ≤
      (tooth toothCount).width := by
  unfold rationalCircleToothOverlapLength
  rw [length_inter_comm]
  have hbound := length_inter_le_left
    [((0 : ℚ), (tooth toothCount).width)]
    (norm_rationalCircleChart shift tooth toothCount)
  simpa [length, max_eq_right (tooth toothCount).width_pos.le] using hbound

/-- A concrete sufficient containment certificate for full-width loss: one chart
interval starts at the anchor and covers the whole tooth. -/
def RationalCircleAnchoredToothContained
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) : Prop :=
  ∃ interval ∈ rationalCircleChart shift tooth toothCount,
    interval.1 = 0 ∧ (tooth toothCount).width ≤ interval.2

/-- A concrete unit-chart disjointness certificate for zero loss: every chart
interval starts at or to the right of the tooth endpoint. -/
def RationalCircleAnchoredToothDisjoint
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) : Prop :=
  ∀ interval ∈ rationalCircleChart shift tooth toothCount,
    (tooth toothCount).width ≤ interval.1

/-- The containment certificate makes the actual overlap equal the full tooth
width. -/
theorem rationalCircleToothOverlapLength_eq_width_of_contained
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ)
    (hcontained : RationalCircleAnchoredToothContained shift tooth toothCount) :
    rationalCircleToothOverlapLength shift tooth toothCount =
      (tooth toothCount).width := by
  rw [rationalCircleToothOverlapLength_eq_clip_sum]
  apply le_antisymm
  · rw [← rationalCircleToothOverlapLength_eq_clip_sum]
    exact rationalCircleToothOverlapLength_le_width shift tooth toothCount
  · obtain ⟨interval, hinterval, hleft, hright⟩ := hcontained
    apply List.single_le_sum
    · intro value hvalue
      rw [List.mem_map] at hvalue
      obtain ⟨source, _hsource, rfl⟩ := hvalue
      exact le_max_left _ _
    · rw [List.mem_map]
      refine ⟨interval, hinterval, ?_⟩
      rw [min_eq_right hright, hleft, sub_zero,
        max_eq_right (tooth toothCount).width_pos.le]

/-- The disjointness certificate makes the actual overlap vanish. -/
theorem rationalCircleToothOverlapLength_eq_zero_of_disjoint
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ)
    (hdisjoint : RationalCircleAnchoredToothDisjoint shift tooth toothCount) :
    rationalCircleToothOverlapLength shift tooth toothCount = 0 := by
  rw [rationalCircleToothOverlapLength_eq_clip_sum]
  apply List.sum_eq_zero
  intro value hvalue
  rw [List.mem_map] at hvalue
  obtain ⟨interval, hinterval, rfl⟩ := hvalue
  apply max_eq_left
  have hleft := hdisjoint interval hinterval
  have hright := min_le_right interval.2 (tooth toothCount).width
  linarith

/-- Minimal algebraic certificate for the naive full-width recurrence: it is
exact precisely when the actual chart-tooth overlap has full width. -/
theorem rationalCircleSurvivor_length_succ_eq_sub_width_iff_overlap_eq_width
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    length (rationalCircleSurvivor shift tooth (toothCount + 1)) =
        length (rationalCircleSurvivor shift tooth toothCount) -
          (tooth toothCount).width ↔
      rationalCircleToothOverlapLength shift tooth toothCount =
        (tooth toothCount).width := by
  rw [rationalCircleSurvivor_length_succ_eq_sub_overlap]
  constructor <;> intro hypothesis <;> linarith

/-- Minimal algebraic certificate for zero density loss: deletion preserves
length precisely when the actual chart-tooth overlap is zero. -/
theorem rationalCircleSurvivor_length_succ_eq_iff_overlap_eq_zero
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    length (rationalCircleSurvivor shift tooth (toothCount + 1)) =
        length (rationalCircleSurvivor shift tooth toothCount) ↔
      rationalCircleToothOverlapLength shift tooth toothCount = 0 := by
  rw [rationalCircleSurvivor_length_succ_eq_sub_overlap]
  constructor <;> intro hypothesis <;> linarith

/-- Under full containment, the naive full-tooth-width recurrence is exact. -/
theorem rationalCircleSurvivor_length_succ_eq_sub_width_of_contained
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ)
    (hcontained : RationalCircleAnchoredToothContained shift tooth toothCount) :
    length (rationalCircleSurvivor shift tooth (toothCount + 1)) =
      length (rationalCircleSurvivor shift tooth toothCount) -
        (tooth toothCount).width := by
  rw [rationalCircleSurvivor_length_succ_eq_sub_overlap,
    rationalCircleToothOverlapLength_eq_width_of_contained
      shift tooth toothCount hcontained]

/-- Under chart-tooth disjointness, deletion leaves the rational length
unchanged. -/
theorem rationalCircleSurvivor_length_succ_eq_of_disjoint
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ)
    (hdisjoint : RationalCircleAnchoredToothDisjoint shift tooth toothCount) :
    length (rationalCircleSurvivor shift tooth (toothCount + 1)) =
      length (rationalCircleSurvivor shift tooth toothCount) := by
  rw [rationalCircleSurvivor_length_succ_eq_sub_overlap,
    rationalCircleToothOverlapLength_eq_zero_of_disjoint
      shift tooth toothCount hdisjoint, sub_zero]

/-- Exact real-density recurrence obtained by casting the rational deletion
ledger. -/
theorem rationalCircleSurvivorDensity_succ_add_overlap
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    rationalCircleSurvivorDensity shift tooth (toothCount + 1) +
        rationalCircleToothOverlapDensity shift tooth toothCount =
      rationalCircleSurvivorDensity shift tooth toothCount := by
  unfold rationalCircleSurvivorDensity rationalRegionDensity
    rationalCircleToothOverlapDensity
  exact_mod_cast rationalCircleSurvivor_length_succ_add_overlap
    shift tooth toothCount

/-- Subtraction form of the exact real-density recurrence. -/
theorem rationalCircleSurvivorDensity_succ_eq_sub_overlap
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    rationalCircleSurvivorDensity shift tooth (toothCount + 1) =
      rationalCircleSurvivorDensity shift tooth toothCount -
        rationalCircleToothOverlapDensity shift tooth toothCount := by
  linarith [rationalCircleSurvivorDensity_succ_add_overlap
    shift tooth toothCount]

/-- Closed exact real-density formula. The sole remaining numeric data are
the casts of the explicit rational clip sums. -/
theorem rationalCircleSurvivorDensity_eq_one_sub_sum_overlap
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth) :
    ∀ toothCount,
      rationalCircleSurvivorDensity shift tooth toothCount =
        1 - ∑ count ∈ Finset.range toothCount,
          rationalCircleToothOverlapDensity shift tooth count := by
  intro toothCount
  unfold rationalCircleSurvivorDensity rationalRegionDensity
    rationalCircleToothOverlapDensity
  rw [rationalCircleSurvivor_length_eq_one_sub_sum_overlap]
  push_cast
  rfl

/-- The exact overlap is nonnegative after casting to real density. -/
theorem rationalCircleToothOverlapDensity_nonneg
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    0 ≤ rationalCircleToothOverlapDensity shift tooth toothCount := by
  unfold rationalCircleToothOverlapDensity
  exact_mod_cast rationalCircleToothOverlapLength_nonneg shift tooth toothCount

/-- Successive exact survivor densities are monotone nonincreasing. -/
theorem rationalCircleSurvivorDensity_succ_le
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    rationalCircleSurvivorDensity shift tooth (toothCount + 1) ≤
      rationalCircleSurvivorDensity shift tooth toothCount := by
  linarith [rationalCircleSurvivorDensity_succ_add_overlap
      shift tooth toothCount,
    rationalCircleToothOverlapDensity_nonneg shift tooth toothCount]

/-- Full containment specializes the exact density recurrence to subtraction
of the full tooth width. -/
theorem rationalCircleSurvivorDensity_succ_eq_sub_width_of_contained
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ)
    (hcontained : RationalCircleAnchoredToothContained shift tooth toothCount) :
    rationalCircleSurvivorDensity shift tooth (toothCount + 1) =
      rationalCircleSurvivorDensity shift tooth toothCount -
        ((tooth toothCount).width : ℝ) := by
  unfold rationalCircleSurvivorDensity rationalRegionDensity
  exact_mod_cast rationalCircleSurvivor_length_succ_eq_sub_width_of_contained
    shift tooth toothCount hcontained

/-- Chart-tooth disjointness specializes the exact density recurrence to no
density loss. -/
theorem rationalCircleSurvivorDensity_succ_eq_of_disjoint
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ)
    (hdisjoint : RationalCircleAnchoredToothDisjoint shift tooth toothCount) :
    rationalCircleSurvivorDensity shift tooth (toothCount + 1) =
      rationalCircleSurvivorDensity shift tooth toothCount := by
  unfold rationalCircleSurvivorDensity rationalRegionDensity
  exact_mod_cast rationalCircleSurvivor_length_succ_eq_of_disjoint
    shift tooth toothCount hdisjoint

/-- The recursive survivor supplies all THM-933 analytic hypotheses without
an additional topology or normalization assumption. -/
theorem rationalCircleSurvivor_analytic_provider
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    MeasurableSet (periodicRationalCircleSurvivor shift tooth toothCount) ∧
      Function.Periodic
        (fun t =>
          (periodicRationalCircleSurvivor shift tooth toothCount).indicator
            (fun _ => (1 : ℝ)) t) 1 ∧
      (volume (periodicRationalCircleSurvivor shift tooth toothCount ∩
        Ioc (0 : ℝ) 1)).toReal =
          rationalCircleSurvivorDensity shift tooth toothCount := by
  simpa [periodicRationalCircleSurvivor, rationalCircleSurvivorDensity] using
    rationalRegion_analytic_provider
      (rationalCircleSurvivor shift tooth toothCount)
      (norm_rationalCircleSurvivor shift tooth toothCount)
      (regionInUnit_rationalCircleSurvivor shift tooth toothCount)

/-- The normalized moving chart supplies all THM-933 analytic hypotheses
without an additional topology or normalization assumption. -/
theorem rationalCircleChart_analytic_provider
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    MeasurableSet (periodicRationalCircleChart shift tooth toothCount) ∧
      Function.Periodic
        (fun t =>
          (periodicRationalCircleChart shift tooth toothCount).indicator
            (fun _ => (1 : ℝ)) t) 1 ∧
      (volume (periodicRationalCircleChart shift tooth toothCount ∩
        Ioc (0 : ℝ) 1)).toReal =
          rationalCircleChartDensity shift tooth toothCount := by
  simpa [periodicRationalCircleChart, rationalCircleChartDensity] using
    rationalRegion_analytic_provider
      (rationalCircleChart shift tooth toothCount)
      (norm_rationalCircleChart shift tooth toothCount)
      (regionInUnit_rationalCircleChart shift tooth toothCount)

/-- Concrete measurability of the recursive rational circle survivor. -/
theorem measurableSet_periodicRationalCircleSurvivor
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    MeasurableSet (periodicRationalCircleSurvivor shift tooth toothCount) :=
  (rationalCircleSurvivor_analytic_provider shift tooth toothCount).1

/-- Concrete measurability of the normalized rational circle chart. -/
theorem measurableSet_periodicRationalCircleChart
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    MeasurableSet (periodicRationalCircleChart shift tooth toothCount) :=
  (rationalCircleChart_analytic_provider shift tooth toothCount).1

/-- Concrete one-periodicity of the recursive survivor indicator. -/
theorem periodicRationalCircleSurvivor_indicator_periodic
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    Function.Periodic
      (fun t =>
        (periodicRationalCircleSurvivor shift tooth toothCount).indicator
          (fun _ => (1 : ℝ)) t) 1 :=
  (rationalCircleSurvivor_analytic_provider shift tooth toothCount).2.1

/-- Concrete one-periodicity of the normalized chart indicator. -/
theorem periodicRationalCircleChart_indicator_periodic
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    Function.Periodic
      (fun t => (periodicRationalCircleChart shift tooth toothCount).indicator
        (fun _ => (1 : ℝ)) t) 1 :=
  (rationalCircleChart_analytic_provider shift tooth toothCount).2.1

/-- Exact one-period mass of the recursive rational circle survivor. -/
theorem periodicRationalCircleSurvivor_period_density
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    (volume (periodicRationalCircleSurvivor shift tooth toothCount ∩
      Ioc (0 : ℝ) 1)).toReal =
        rationalCircleSurvivorDensity shift tooth toothCount :=
  (rationalCircleSurvivor_analytic_provider shift tooth toothCount).2.2

/-- Exact one-period mass of the normalized rational circle chart. -/
theorem periodicRationalCircleChart_period_density
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    (volume (periodicRationalCircleChart shift tooth toothCount ∩
      Ioc (0 : ℝ) 1)).toReal =
        rationalCircleChartDensity shift tooth toothCount :=
  (rationalCircleChart_analytic_provider shift tooth toothCount).2.2

/-- The chart period mass is the survivor density, so recharting introduces
no analytic debt. -/
theorem periodicRationalCircleChart_period_density_eq_survivorDensity
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    (volume (periodicRationalCircleChart shift tooth toothCount ∩
      Ioc (0 : ℝ) 1)).toReal =
        rationalCircleSurvivorDensity shift tooth toothCount := by
  rw [periodicRationalCircleChart_period_density,
    rationalCircleChartDensity_eq_survivorDensity]

/-- Exact one-period measure recurrence under the next anchored deletion. -/
theorem periodicRationalCircleSurvivor_period_density_succ_eq_sub_overlap
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    (volume (periodicRationalCircleSurvivor shift tooth (toothCount + 1) ∩
        Ioc (0 : ℝ) 1)).toReal =
      (volume (periodicRationalCircleSurvivor shift tooth toothCount ∩
        Ioc (0 : ℝ) 1)).toReal -
        rationalCircleToothOverlapDensity shift tooth toothCount := by
  rw [periodicRationalCircleSurvivor_period_density,
    periodicRationalCircleSurvivor_period_density]
  exact rationalCircleSurvivorDensity_succ_eq_sub_overlap
    shift tooth toothCount

/-- Closed exact one-period survivor mass formula in terms of the finite
rational overlap ledger. -/
theorem periodicRationalCircleSurvivor_period_density_eq_one_sub_sum_overlap
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    (volume (periodicRationalCircleSurvivor shift tooth toothCount ∩
        Ioc (0 : ℝ) 1)).toReal =
      1 - ∑ count ∈ Finset.range toothCount,
        rationalCircleToothOverlapDensity shift tooth count := by
  rw [periodicRationalCircleSurvivor_period_density,
    rationalCircleSurvivorDensity_eq_one_sub_sum_overlap]

/-- The moving chart has the same closed exact period-mass formula as its
survivor. -/
theorem periodicRationalCircleChart_period_density_eq_one_sub_sum_overlap
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    (volume (periodicRationalCircleChart shift tooth toothCount ∩
        Ioc (0 : ℝ) 1)).toReal =
      1 - ∑ count ∈ Finset.range toothCount,
        rationalCircleToothOverlapDensity shift tooth count := by
  rw [periodicRationalCircleChart_period_density_eq_survivorDensity,
    rationalCircleSurvivorDensity_eq_one_sub_sum_overlap]

/-- The survivor's centered primitive is one-periodic. -/
theorem rationalCircleSurvivor_centeredPrimitive_periodic
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    Function.Periodic
      (centeredPrimitive (periodicRationalCircleSurvivor shift tooth toothCount)
        (rationalCircleSurvivorDensity shift tooth toothCount)) 1 := by
  rcases rationalCircleSurvivor_analytic_provider shift tooth toothCount with
    ⟨hmeasurable, hperiodic, hdensity⟩
  exact centeredPrimitive_periodic _ hmeasurable _ hperiodic hdensity

/-- The moving chart's centered primitive is one-periodic. -/
theorem rationalCircleChart_centeredPrimitive_periodic
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    Function.Periodic
      (centeredPrimitive (periodicRationalCircleChart shift tooth toothCount)
        (rationalCircleChartDensity shift tooth toothCount)) 1 := by
  rcases rationalCircleChart_analytic_provider shift tooth toothCount with
    ⟨hmeasurable, hperiodic, hdensity⟩
  exact centeredPrimitive_periodic _ hmeasurable _ hperiodic hdensity

/-- Every length-one survivor window has exactly the symbolic survivor
density. -/
theorem rationalCircleSurvivor_retainedWindow_one_eq_density
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) (startPoint : ℝ) :
    retainedWindow (periodicRationalCircleSurvivor shift tooth toothCount)
        1 startPoint =
      rationalCircleSurvivorDensity shift tooth toothCount := by
  rcases rationalCircleSurvivor_analytic_provider shift tooth toothCount with
    ⟨hmeasurable, hperiodic, hdensity⟩
  exact retainedWindow_one_eq_density _ hmeasurable _ startPoint hperiodic hdensity

/-- Every length-one chart window has exactly the symbolic chart density. -/
theorem rationalCircleChart_retainedWindow_one_eq_density
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) (startPoint : ℝ) :
    retainedWindow (periodicRationalCircleChart shift tooth toothCount)
        1 startPoint =
      rationalCircleChartDensity shift tooth toothCount := by
  rcases rationalCircleChart_analytic_provider shift tooth toothCount with
    ⟨hmeasurable, hperiodic, hdensity⟩
  exact retainedWindow_one_eq_density _ hmeasurable _ startPoint hperiodic hdensity

/-- At full-period scale, survivor `eta` is its exact symbolic density. -/
theorem rationalCircleSurvivor_fixedScaleEta_one
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    fixedScaleEta (periodicRationalCircleSurvivor shift tooth toothCount) 1 =
      rationalCircleSurvivorDensity shift tooth toothCount := by
  rcases rationalCircleSurvivor_analytic_provider shift tooth toothCount with
    ⟨hmeasurable, hperiodic, hdensity⟩
  exact fixedScaleEta_one _ hmeasurable _ hperiodic hdensity

/-- At full-period scale, chart `eta` is its exact symbolic density. -/
theorem rationalCircleChart_fixedScaleEta_one
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    fixedScaleEta (periodicRationalCircleChart shift tooth toothCount) 1 =
      rationalCircleChartDensity shift tooth toothCount := by
  rcases rationalCircleChart_analytic_provider shift tooth toothCount with
    ⟨hmeasurable, hperiodic, hdensity⟩
  exact fixedScaleEta_one _ hmeasurable _ hperiodic hdensity

/-- Every positive-scale survivor deficit is bounded by its primitive
discrepancy. -/
theorem rationalCircleSurvivor_fixedScaleEta_deficit_le_primitiveDiscrepancy
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) {ell : ℝ} (hell : 0 < ell) :
    ell * (rationalCircleSurvivorDensity shift tooth toothCount -
      fixedScaleEta (periodicRationalCircleSurvivor shift tooth toothCount) ell) ≤
        primitiveDiscrepancy
          (periodicRationalCircleSurvivor shift tooth toothCount)
          (rationalCircleSurvivorDensity shift tooth toothCount) := by
  rcases rationalCircleSurvivor_analytic_provider shift tooth toothCount with
    ⟨hmeasurable, hperiodic, hdensity⟩
  exact fixedScaleEta_deficit_le_primitiveDiscrepancy
    _ hmeasurable _ hperiodic hdensity hell

/-- Every positive-scale chart deficit is bounded by its primitive
discrepancy. -/
theorem rationalCircleChart_fixedScaleEta_deficit_le_primitiveDiscrepancy
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) {ell : ℝ} (hell : 0 < ell) :
    ell * (rationalCircleChartDensity shift tooth toothCount -
      fixedScaleEta (periodicRationalCircleChart shift tooth toothCount) ell) ≤
        primitiveDiscrepancy
          (periodicRationalCircleChart shift tooth toothCount)
          (rationalCircleChartDensity shift tooth toothCount) := by
  rcases rationalCircleChart_analytic_provider shift tooth toothCount with
    ⟨hmeasurable, hperiodic, hdensity⟩
  exact fixedScaleEta_deficit_le_primitiveDiscrepancy
    _ hmeasurable _ hperiodic hdensity hell

/-- Survivor eta/q equality is attained at some scale `0 < ell ≤ 1`. -/
theorem rationalCircleSurvivor_eta_q_attained
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    ∃ ell ∈ Ioc (0 : ℝ) 1,
      ell * (rationalCircleSurvivorDensity shift tooth toothCount -
        fixedScaleEta (periodicRationalCircleSurvivor shift tooth toothCount) ell) =
          primitiveDiscrepancy
            (periodicRationalCircleSurvivor shift tooth toothCount)
            (rationalCircleSurvivorDensity shift tooth toothCount) := by
  simpa [periodicRationalCircleSurvivor, rationalCircleSurvivorDensity] using
    rationalRegion_eta_q_attained
      (rationalCircleSurvivor shift tooth toothCount)
      (norm_rationalCircleSurvivor shift tooth toothCount)
      (regionInUnit_rationalCircleSurvivor shift tooth toothCount)

/-- Chart eta/q equality is attained at some scale `0 < ell ≤ 1`. -/
theorem rationalCircleChart_eta_q_attained
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    ∃ ell ∈ Ioc (0 : ℝ) 1,
      ell * (rationalCircleChartDensity shift tooth toothCount -
        fixedScaleEta (periodicRationalCircleChart shift tooth toothCount) ell) =
          primitiveDiscrepancy
            (periodicRationalCircleChart shift tooth toothCount)
            (rationalCircleChartDensity shift tooth toothCount) := by
  simpa [periodicRationalCircleChart, rationalCircleChartDensity] using
    rationalRegion_eta_q_attained
      (rationalCircleChart shift tooth toothCount)
      (norm_rationalCircleChart shift tooth toothCount)
      (regionInUnit_rationalCircleChart shift tooth toothCount)

/-- Full attained eta/q duality for the recursive rational circle survivor. -/
theorem rationalCircleSurvivor_eta_q_duality
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    primitiveDiscrepancy
        (periodicRationalCircleSurvivor shift tooth toothCount)
        (rationalCircleSurvivorDensity shift tooth toothCount) =
      fixedScaleDeficitSup
        (periodicRationalCircleSurvivor shift tooth toothCount)
        (rationalCircleSurvivorDensity shift tooth toothCount) := by
  simpa [periodicRationalCircleSurvivor, rationalCircleSurvivorDensity] using
    rationalRegion_eta_q_duality
      (rationalCircleSurvivor shift tooth toothCount)
      (norm_rationalCircleSurvivor shift tooth toothCount)
      (regionInUnit_rationalCircleSurvivor shift tooth toothCount)

/-- Full attained eta/q duality for the normalized rational circle chart. -/
theorem rationalCircleChart_eta_q_duality
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    primitiveDiscrepancy
        (periodicRationalCircleChart shift tooth toothCount)
        (rationalCircleChartDensity shift tooth toothCount) =
      fixedScaleDeficitSup
        (periodicRationalCircleChart shift tooth toothCount)
        (rationalCircleChartDensity shift tooth toothCount) := by
  simpa [periodicRationalCircleChart, rationalCircleChartDensity] using
    rationalRegion_eta_q_duality
      (rationalCircleChart shift tooth toothCount)
      (norm_rationalCircleChart shift tooth toothCount)
      (regionInUnit_rationalCircleChart shift tooth toothCount)

/-! ## Axiom audit -/

#print axioms ratCast_mem_rationalRegionSet
#print axioms volume_rationalRegionSet
#print axioms periodicRationalRegion_indicator_periodic
#print axioms periodicRationalRegion_period_density
#print axioms rationalRegion_analytic_provider
#print axioms rationalRegion_centeredPrimitive_interval_identity
#print axioms rationalRegion_fixedScaleEta_attained
#print axioms rationalRegion_eta_q_attained
#print axioms rationalRegion_eta_q_duality
#print axioms rationalCircleChartDensity_eq_survivorDensity
#print axioms length_deleteAnchoredTooth_add_overlap
#print axioms rationalCircleToothOverlapLength_eq_clip_sum
#print axioms rationalCircleSurvivor_length_succ_add_overlap
#print axioms rationalCircleSurvivor_length_eq_one_sub_sum_overlap
#print axioms rationalCircleToothOverlapLength_eq_width_of_contained
#print axioms rationalCircleToothOverlapLength_eq_zero_of_disjoint
#print axioms rationalCircleSurvivorDensity_succ_add_overlap
#print axioms rationalCircleSurvivorDensity_eq_one_sub_sum_overlap
#print axioms periodicRationalCircleSurvivor_period_density_eq_one_sub_sum_overlap
#print axioms rationalCircleSurvivor_analytic_provider
#print axioms rationalCircleChart_analytic_provider
#print axioms rationalCircleSurvivor_centeredPrimitive_periodic
#print axioms rationalCircleChart_centeredPrimitive_periodic
#print axioms rationalCircleSurvivor_eta_q_attained
#print axioms rationalCircleChart_eta_q_attained
#print axioms rationalCircleSurvivor_eta_q_duality
#print axioms rationalCircleChart_eta_q_duality

end LRC14.LocalDensityBlockGluing
