/-
  TournamentH7.LRCRationalRegionProvider

  Concrete rational-interval analytic provider for THM-933. A normalized
  `RatIntervals.Region` is cast to a finite union of real half-open intervals,
  then periodized through `Int.fract`. The resulting survivor is measurable,
  its indicator is one-periodic, and its mass on one period is exactly the cast
  rational region length. These facts instantiate the centered-primitive and
  eta/q duality theorems from `LRCLocalDensityBlockGluing`.

  The exact-measure bridge needs precisely the rational measure discipline:
  `Norm region` for live, ordered, disjoint intervals, plus `regionInUnit region`
  to identify the listed chart with the fundamental period. No component-count
  topology is used.

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

end LRC14.LocalDensityBlockGluing
