/-
  TournamentH7.LRCLocalDensityBlockGluing

  Analytic and algebraic formalization core of THM-933, the sharp local-density
  block-gluing theorem (codex-2026-07-16-S21, HYP-7152).

  The geometric input to THM-933 says that every component J of an earlier
  survivor retains at least

      density * length(J) - discrepancy

  inside the next certified block.  This file proves the reusable pieces
  around that input:

  * the Lebesgue measure of a retained interval is exactly the increment of the
    centered indicator primitive, for both `Ioc` and `Icc` conventions;
  * the primitive and fixed-scale retained-window functions are continuous, so
    the primitive discrepancy `q` and fixed-scale density `eta` are attained;
  * for a one-periodic block of the advertised density, every fixed-scale
    deficit is at most the attained primitive discrepancy;
  * a bounded centered primitive gives the one-interval discrepancy inequality,
    with equality when an interval joins an upper extremum to a lower extremum;
  * summing the local inequality over components pays exactly card * q, and a
    component cap turns this into M * q;
  * exact rational interval subtraction proves the concrete cut-open-circle
    one-tooth split, and its atlas feeds the N-tooth component cap;
  * any sequence of nonnegative-density recurrence steps is bounded below by
    the recursive product-minus-weighted-debt ledger.

  It also kernel-checks the R = 7 LRC(14) arithmetic and the exact three-block
  certificate from the THM-933 referee.  No native_decide and no sorry.
-/

import Mathlib.MeasureTheory.Integral.IntervalIntegral.Periodic
import Mathlib.Tactic
import Mathlib.Topology.Order.Compact
import TournamentH7.LRCRegionDiff

open scoped BigOperators

namespace LRC14.LocalDensityBlockGluing

section PrimitiveBridge

variable {Point : Type*}

/-- The analytic core of THM-933 Lemma 1.  Once retained interval measure is
represented as a difference of a bounded centered primitive, its oscillation
is a valid additive discrepancy loss. -/
theorem primitive_interval_lower
    (primitive : Point → ℝ) (startPoint endPoint : Point)
    (density intervalLength retained lowerValue upperValue : ℝ)
    (hlower : ∀ point, lowerValue ≤ primitive point)
    (hupper : ∀ point, primitive point ≤ upperValue)
    (hidentity :
      retained - density * intervalLength =
        primitive endPoint - primitive startPoint) :
    density * intervalLength - (upperValue - lowerValue) ≤ retained := by
  have hstart := hupper startPoint
  have hend := hlower endPoint
  linarith

/-- The primitive loss is sharp when the oriented interval starts at an upper
extremum and ends at a lower extremum. -/
theorem primitive_interval_sharp
    (primitive : Point → ℝ) (startPoint endPoint : Point)
    (density intervalLength retained lowerValue upperValue : ℝ)
    (hstart : primitive startPoint = upperValue)
    (hend : primitive endPoint = lowerValue)
    (hidentity :
      retained - density * intervalLength =
        primitive endPoint - primitive startPoint) :
    retained = density * intervalLength - (upperValue - lowerValue) := by
  rw [hstart, hend] at hidentity
  linarith

/-- An attained fixed-scale deficit lies below the primitive oscillation.
Taking the supremum over scales and using `primitive_interval_sharp` is the
paper proof of `q = sup_ell ell * (delta - eta ell)`. -/
theorem fixedScale_deficit_le_discrepancy
    (primitive : Point → ℝ) (startPoint endPoint : Point)
    (density eta intervalLength retained lowerValue upperValue : ℝ)
    (hlower : ∀ point, lowerValue ≤ primitive point)
    (hupper : ∀ point, primitive point ≤ upperValue)
    (heta : retained = eta * intervalLength)
    (hidentity :
      retained - density * intervalLength =
        primitive endPoint - primitive startPoint) :
    intervalLength * (density - eta) ≤ upperValue - lowerValue := by
  have hstart := hupper startPoint
  have hend := hlower endPoint
  have hprimitive :
      primitive startPoint - primitive endPoint ≤ upperValue - lowerValue :=
    sub_le_sub hstart hend
  calc
    intervalLength * (density - eta)
        = density * intervalLength - eta * intervalLength := by ring
    _ = density * intervalLength - retained := by rw [heta]
    _ = primitive startPoint - primitive endPoint := by linarith
    _ ≤ upperValue - lowerValue := hprimitive

/-- At a fixed-scale minimizer whose arc joins primitive extrema, the
fixed-scale deficit equals the primitive discrepancy exactly. -/
theorem fixedScale_extremizer_eq_discrepancy
    (primitive : Point → ℝ) (startPoint endPoint : Point)
    (density eta intervalLength retained lowerValue upperValue : ℝ)
    (heta : retained = eta * intervalLength)
    (hstart : primitive startPoint = upperValue)
    (hend : primitive endPoint = lowerValue)
    (hidentity :
      retained - density * intervalLength =
        primitive endPoint - primitive startPoint) :
    intervalLength * (density - eta) = upperValue - lowerValue := by
  have hsharp := primitive_interval_sharp primitive startPoint endPoint
    density intervalLength retained lowerValue upperValue hstart hend hidentity
  calc
    intervalLength * (density - eta)
        = density * intervalLength - eta * intervalLength := by ring
    _ = density * intervalLength - retained := by rw [heta]
    _ = upperValue - lowerValue := by linarith

end PrimitiveBridge

section ConcretePrimitiveBridge

open MeasureTheory Set

/-- The centered indicator primitive from THM-933:
`H(u) = ∫₀ᵘ (1_S(t) - density) dt`. -/
noncomputable def centeredPrimitive (S : Set ℝ) (density u : ℝ) : ℝ :=
  ∫ t in (0 : ℝ)..u, S.indicator (fun _ => (1 : ℝ)) t - density

/-- Retained Lebesgue measure in the oriented positive-length window
`(startPoint, startPoint + ell]`. Endpoints are null, so this agrees with the
closed-interval convention used in THM-933. -/
noncomputable def retainedWindow (S : Set ℝ) (ell startPoint : ℝ) : ℝ :=
  (volume (S ∩ Ioc startPoint (startPoint + ell))).toReal

/-- Fixed-scale local density `eta(ell)`, defined as the infimum over one
period. The infimum is attained for measurable `S` and `ell > 0`; see
`fixedScaleEta_attained`. -/
noncomputable def fixedScaleEta (S : Set ℝ) (ell : ℝ) : ℝ :=
  sInf ((fun startPoint => retainedWindow S ell startPoint / ell) '' Icc (0 : ℝ) 1)

/-- Primitive discrepancy `q = max H - min H`, expressed using `sSup` and
`sInf` on one period. Both extrema are attained; see
`primitiveDiscrepancy_attained`. -/
noncomputable def primitiveDiscrepancy (S : Set ℝ) (density : ℝ) : ℝ :=
  sSup (centeredPrimitive S density '' Icc (0 : ℝ) 1) -
    sInf (centeredPrimitive S density '' Icc (0 : ℝ) 1)

/-- Supremum over the fixed-scale deficit ledger on scales `0 < ell ≤ 1`.
THM-933 identifies this with `primitiveDiscrepancy`; see
`primitiveDiscrepancy_eq_fixedScaleDeficitSup`. -/
noncomputable def fixedScaleDeficitSup (S : Set ℝ) (density : ℝ) : ℝ :=
  sSup ((fun ell => ell * (density - fixedScaleEta S ell)) '' Ioc (0 : ℝ) 1)

/-- A measurable indicator is integrable on every bounded real interval. -/
lemma indicator_intervalIntegrable (S : Set ℝ) (hS : MeasurableSet S) (a b : ℝ) :
    IntervalIntegrable (fun t => S.indicator (fun _ => (1 : ℝ)) t) volume a b := by
  have hone : IntervalIntegrable (fun _ : ℝ => (1 : ℝ)) volume a b :=
    intervalIntegrable_const
  constructor
  · exact hone.1.indicator hS
  · exact hone.2.indicator hS

/-- The centered indicator is integrable on every bounded real interval. -/
lemma centered_intervalIntegrable (S : Set ℝ) (hS : MeasurableSet S)
    (density a b : ℝ) :
    IntervalIntegrable (fun t => S.indicator (fun _ => (1 : ℝ)) t - density)
      volume a b :=
  (indicator_intervalIntegrable S hS a b).sub intervalIntegrable_const

/-- The interval integral of an indicator is the retained Lebesgue measure.
This is the concrete measure-theoretic input behind the THM-933 primitive
bridge. -/
theorem intervalIndicatorIntegral_eq_measure
    (S : Set ℝ) (hS : MeasurableSet S) {a b : ℝ} (hab : a ≤ b) :
    (∫ t in a..b, S.indicator (fun _ => (1 : ℝ)) t) =
      (volume (S ∩ Ioc a b)).toReal := by
  rw [intervalIntegral.integral_of_le hab]
  rw [← MeasureTheory.integral_indicator measurableSet_Ioc]
  rw [Set.indicator_indicator, inter_comm]
  exact MeasureTheory.integral_indicator_one (hS.inter measurableSet_Ioc)

/-- Concrete THM-933 primitive identity with the half-open interval convention:
`volume(S ∩ (a,b]) - density * (b-a) = H(b) - H(a)`. -/
theorem centeredPrimitive_interval_identity
    (S : Set ℝ) (hS : MeasurableSet S) (density : ℝ) {a b : ℝ} (hab : a ≤ b) :
    (volume (S ∩ Ioc a b)).toReal - density * (b - a) =
      centeredPrimitive S density b - centeredPrimitive S density a := by
  have hind := indicator_intervalIntegrable S hS a b
  have hcenter0b := centered_intervalIntegrable S hS density 0 b
  have hcenter0a := centered_intervalIntegrable S hS density 0 a
  calc
    (volume (S ∩ Ioc a b)).toReal - density * (b - a)
        = (∫ t in a..b, S.indicator (fun _ => (1 : ℝ)) t) -
            ∫ _t in a..b, density := by
              rw [intervalIndicatorIntegral_eq_measure S hS hab,
                intervalIntegral.integral_const, smul_eq_mul]
              ring
    _ = ∫ t in a..b, S.indicator (fun _ => (1 : ℝ)) t - density := by
      rw [intervalIntegral.integral_sub hind intervalIntegrable_const]
    _ = centeredPrimitive S density b - centeredPrimitive S density a := by
      rw [centeredPrimitive]
      exact (intervalIntegral.integral_interval_sub_left hcenter0b hcenter0a).symm

/-- Closed-interval form of `centeredPrimitive_interval_identity`. Lebesgue
measure ignores the left endpoint, so it is identical to the half-open form. -/
theorem centeredPrimitive_Icc_identity
    (S : Set ℝ) (hS : MeasurableSet S) (density : ℝ) {a b : ℝ} (hab : a ≤ b) :
    (volume (S ∩ Icc a b)).toReal - density * (b - a) =
      centeredPrimitive S density b - centeredPrimitive S density a := by
  have hnull : volume (S ∩ Ioc a b) = volume (S ∩ Icc a b) :=
    measure_congr (ae_eq_set_inter (ae_eq_refl S) Ioc_ae_eq_Icc)
  rw [← hnull]
  exact centeredPrimitive_interval_identity S hS density hab

/-- The centered primitive is continuous for every measurable block. -/
theorem continuous_centeredPrimitive
    (S : Set ℝ) (hS : MeasurableSet S) (density : ℝ) :
    Continuous (centeredPrimitive S density) := by
  exact intervalIntegral.continuous_primitive
    (fun a b => centered_intervalIntegrable S hS density a b) 0

/-- If the indicator is one-periodic and `density` is its mass on one period,
then the centered primitive is one-periodic as asserted in THM-933. -/
theorem centeredPrimitive_periodic
    (S : Set ℝ) (hS : MeasurableSet S) (density : ℝ)
    (hindicator : Function.Periodic
      (fun t => S.indicator (fun _ => (1 : ℝ)) t) 1)
    (hdensity : (volume (S ∩ Ioc (0 : ℝ) 1)).toReal = density) :
    Function.Periodic (centeredPrimitive S density) 1 := by
  let centeredIntegrand := fun t => S.indicator (fun _ => (1 : ℝ)) t - density
  have hperiodicIntegrand : Function.Periodic centeredIntegrand 1 := by
    intro point
    exact congrArg (fun value => value - density) (hindicator point)
  have hmean : (∫ t in (0 : ℝ)..1, centeredIntegrand t) = 0 := by
    have hidentity := centeredPrimitive_interval_identity S hS density
      (show (0 : ℝ) ≤ 1 by norm_num)
    simpa [centeredPrimitive, centeredIntegrand, hdensity] using hidentity.symm
  intro point
  change (∫ t in (0 : ℝ)..point + 1, centeredIntegrand t) =
    ∫ t in (0 : ℝ)..point, centeredIntegrand t
  calc
    (∫ t in (0 : ℝ)..point + 1, centeredIntegrand t) =
        (∫ t in (0 : ℝ)..point, centeredIntegrand t) +
          ∫ t in (0 : ℝ)..0 + 1, centeredIntegrand t :=
      hperiodicIntegrand.intervalIntegral_add_eq_add 0 point
        (fun a b => centered_intervalIntegrable S hS density a b)
    _ = ∫ t in (0 : ℝ)..point, centeredIntegrand t := by
      rw [zero_add, hmean, add_zero]

/-- A one-periodic block retains exactly its advertised one-period mass in
every length-one window. -/
theorem retainedWindow_one_eq_density
    (S : Set ℝ) (hS : MeasurableSet S) (density startPoint : ℝ)
    (hindicator : Function.Periodic
      (fun t => S.indicator (fun _ => (1 : ℝ)) t) 1)
    (hdensity : (volume (S ∩ Ioc (0 : ℝ) 1)).toReal = density) :
    retainedWindow S 1 startPoint = density := by
  rw [retainedWindow, ← intervalIndicatorIntegral_eq_measure S hS
    (show startPoint ≤ startPoint + 1 by linarith)]
  calc
    (∫ t in startPoint..startPoint + 1,
        S.indicator (fun _ => (1 : ℝ)) t) =
        ∫ t in (0 : ℝ)..0 + 1, S.indicator (fun _ => (1 : ℝ)) t :=
      hindicator.intervalIntegral_add_eq startPoint 0
    _ = density := by
      rw [zero_add, intervalIndicatorIntegral_eq_measure S hS (by norm_num), hdensity]

/-- A retained window is a difference of the uncentered primitive. -/
theorem retainedWindow_eq_primitive_sub
    (S : Set ℝ) (hS : MeasurableSet S) {ell : ℝ} (hell : 0 ≤ ell)
    (startPoint : ℝ) :
    retainedWindow S ell startPoint =
      centeredPrimitive S 0 (startPoint + ell) - centeredPrimitive S 0 startPoint := by
  rw [retainedWindow]
  have hidentity := centeredPrimitive_interval_identity S hS 0
    (show startPoint ≤ startPoint + ell by linarith)
  simpa using hidentity

/-- Retained measure in a fixed-length sliding window is continuous in its
starting point. This supplies the compactness premise for `eta`. -/
theorem continuous_retainedWindow
    (S : Set ℝ) (hS : MeasurableSet S) {ell : ℝ} (hell : 0 ≤ ell) :
    Continuous (retainedWindow S ell) := by
  have hprimitive := continuous_centeredPrimitive S hS 0
  have heq : retainedWindow S ell = fun startPoint =>
      centeredPrimitive S 0 (startPoint + ell) - centeredPrimitive S 0 startPoint := by
    funext startPoint
    exact retainedWindow_eq_primitive_sub S hS hell startPoint
  rw [heq]
  fun_prop

/-- The fixed-scale density has a minimizing start point on one period. -/
theorem exists_fixedScale_minimizer
    (S : Set ℝ) (hS : MeasurableSet S) {ell : ℝ} (hell : 0 < ell) :
    ∃ startPoint ∈ Icc (0 : ℝ) 1, ∀ point ∈ Icc (0 : ℝ) 1,
      retainedWindow S ell startPoint / ell ≤ retainedWindow S ell point / ell := by
  have hcontinuous : Continuous (fun point => retainedWindow S ell point / ell) :=
    (continuous_retainedWindow S hS hell.le).div_const ell
  obtain ⟨startPoint, hstartPoint, hmin⟩ :=
    isCompact_Icc.exists_isMinOn (show (Icc (0 : ℝ) 1).Nonempty by norm_num)
      hcontinuous.continuousOn
  exact ⟨startPoint, hstartPoint, fun point hpoint => hmin hpoint⟩

/-- `fixedScaleEta` is a genuine minimum, not merely an infimum. -/
theorem fixedScaleEta_attained
    (S : Set ℝ) (hS : MeasurableSet S) {ell : ℝ} (hell : 0 < ell) :
    ∃ startPoint ∈ Icc (0 : ℝ) 1,
      fixedScaleEta S ell = retainedWindow S ell startPoint / ell ∧
      ∀ point ∈ Icc (0 : ℝ) 1,
        fixedScaleEta S ell ≤ retainedWindow S ell point / ell := by
  let densityAt := fun point => retainedWindow S ell point / ell
  have hcontinuous : Continuous densityAt :=
    (continuous_retainedWindow S hS hell.le).div_const ell
  have hnonempty : (Icc (0 : ℝ) 1).Nonempty := by norm_num
  obtain ⟨startPoint, hstartPoint, hmin⟩ :=
    exists_fixedScale_minimizer S hS hell
  have hbddBelow : BddBelow (densityAt '' Icc (0 : ℝ) 1) :=
    isCompact_Icc.bddBelow_image hcontinuous.continuousOn
  have himageNonempty : (densityAt '' Icc (0 : ℝ) 1).Nonempty :=
    hnonempty.image densityAt
  have heta : fixedScaleEta S ell = densityAt startPoint := by
    unfold fixedScaleEta
    apply le_antisymm
    · exact csInf_le hbddBelow ⟨startPoint, hstartPoint, rfl⟩
    · exact le_csInf himageNonempty (by
        rintro value ⟨point, hpoint, rfl⟩
        exact hmin point hpoint)
  exact ⟨startPoint, hstartPoint, heta, fun point hpoint => heta.symm ▸ hmin point hpoint⟩

/-- At scale one, the fixed-scale minimum is exactly the period density. -/
theorem fixedScaleEta_one
    (S : Set ℝ) (hS : MeasurableSet S) (density : ℝ)
    (hindicator : Function.Periodic
      (fun t => S.indicator (fun _ => (1 : ℝ)) t) 1)
    (hdensity : (volume (S ∩ Ioc (0 : ℝ) 1)).toReal = density) :
    fixedScaleEta S 1 = density := by
  obtain ⟨startPoint, _hstartPoint, heta, _hminimum⟩ :=
    fixedScaleEta_attained S hS (show (0 : ℝ) < 1 by norm_num)
  calc
    fixedScaleEta S 1 = retainedWindow S 1 startPoint / 1 := heta
    _ = density := by rw [retainedWindow_one_eq_density S hS density startPoint
      hindicator hdensity, div_one]

/-- The centered primitive attains both its minimum and maximum on one period. -/
theorem exists_centeredPrimitive_extrema
    (S : Set ℝ) (hS : MeasurableSet S) (density : ℝ) :
    ∃ lowerPoint ∈ Icc (0 : ℝ) 1, ∃ upperPoint ∈ Icc (0 : ℝ) 1,
      (∀ point ∈ Icc (0 : ℝ) 1,
        centeredPrimitive S density lowerPoint ≤ centeredPrimitive S density point) ∧
      (∀ point ∈ Icc (0 : ℝ) 1,
        centeredPrimitive S density point ≤ centeredPrimitive S density upperPoint) := by
  have hcontinuous : ContinuousOn (centeredPrimitive S density) (Icc (0 : ℝ) 1) :=
    (continuous_centeredPrimitive S hS density).continuousOn
  have hnonempty : (Icc (0 : ℝ) 1).Nonempty := by norm_num
  obtain ⟨lowerPoint, hlowerPoint, hlower⟩ :=
    isCompact_Icc.exists_isMinOn hnonempty hcontinuous
  obtain ⟨upperPoint, hupperPoint, hupper⟩ :=
    isCompact_Icc.exists_isMaxOn hnonempty hcontinuous
  exact ⟨lowerPoint, hlowerPoint, upperPoint, hupperPoint,
    fun point hpoint => hlower hpoint, fun point hpoint => hupper hpoint⟩

/-- `primitiveDiscrepancy` is the difference of attained extrema. -/
theorem primitiveDiscrepancy_attained
    (S : Set ℝ) (hS : MeasurableSet S) (density : ℝ) :
    ∃ lowerPoint ∈ Icc (0 : ℝ) 1, ∃ upperPoint ∈ Icc (0 : ℝ) 1,
      primitiveDiscrepancy S density =
          centeredPrimitive S density upperPoint - centeredPrimitive S density lowerPoint ∧
      (∀ point ∈ Icc (0 : ℝ) 1,
        centeredPrimitive S density lowerPoint ≤ centeredPrimitive S density point) ∧
      (∀ point ∈ Icc (0 : ℝ) 1,
        centeredPrimitive S density point ≤ centeredPrimitive S density upperPoint) := by
  let primitive := centeredPrimitive S density
  have hcontinuous : ContinuousOn primitive (Icc (0 : ℝ) 1) :=
    (continuous_centeredPrimitive S hS density).continuousOn
  have hnonempty : (Icc (0 : ℝ) 1).Nonempty := by norm_num
  obtain ⟨lowerPoint, hlowerPoint, upperPoint, hupperPoint, hlower, hupper⟩ :=
    exists_centeredPrimitive_extrema S hS density
  have hbddBelow : BddBelow (primitive '' Icc (0 : ℝ) 1) :=
    isCompact_Icc.bddBelow_image hcontinuous
  have hbddAbove : BddAbove (primitive '' Icc (0 : ℝ) 1) :=
    isCompact_Icc.bddAbove_image hcontinuous
  have himageNonempty : (primitive '' Icc (0 : ℝ) 1).Nonempty :=
    hnonempty.image primitive
  have hinf : sInf (primitive '' Icc (0 : ℝ) 1) = primitive lowerPoint := by
    apply le_antisymm
    · exact csInf_le hbddBelow ⟨lowerPoint, hlowerPoint, rfl⟩
    · exact le_csInf himageNonempty (by
        rintro value ⟨point, hpoint, rfl⟩
        exact hlower point hpoint)
  have hsup : sSup (primitive '' Icc (0 : ℝ) 1) = primitive upperPoint := by
    apply le_antisymm
    · exact csSup_le himageNonempty (by
        rintro value ⟨point, hpoint, rfl⟩
        exact hupper point hpoint)
    · exact le_csSup hbddAbove ⟨upperPoint, hupperPoint, rfl⟩
  exact ⟨lowerPoint, hlowerPoint, upperPoint, hupperPoint, by
    simp only [primitiveDiscrepancy, primitive, hsup, hinf], hlower, hupper⟩

/-- Primitive discrepancy is nonnegative. -/
theorem primitiveDiscrepancy_nonneg
    (S : Set ℝ) (hS : MeasurableSet S) (density : ℝ) :
    0 ≤ primitiveDiscrepancy S density := by
  obtain ⟨lowerPoint, _hlowerPoint, upperPoint, hupperPoint, hq, hlower, _hupper⟩ :=
    primitiveDiscrepancy_attained S hS density
  rw [hq]
  exact sub_nonneg.mpr (hlower upperPoint hupperPoint)

/-- The attained primitive oscillation is nonnegative and bounds every
oriented difference with endpoints in one period. -/
theorem exists_attained_discrepancy
    (S : Set ℝ) (hS : MeasurableSet S) (density : ℝ) :
    ∃ lowerPoint ∈ Icc (0 : ℝ) 1, ∃ upperPoint ∈ Icc (0 : ℝ) 1,
      0 ≤ centeredPrimitive S density upperPoint - centeredPrimitive S density lowerPoint ∧
      ∀ startPoint ∈ Icc (0 : ℝ) 1, ∀ endPoint ∈ Icc (0 : ℝ) 1,
        centeredPrimitive S density startPoint - centeredPrimitive S density endPoint ≤
          centeredPrimitive S density upperPoint - centeredPrimitive S density lowerPoint := by
  obtain ⟨lowerPoint, hlowerPoint, upperPoint, hupperPoint, hlower, hupper⟩ :=
    exists_centeredPrimitive_extrema S hS density
  refine ⟨lowerPoint, hlowerPoint, upperPoint, hupperPoint, ?_, ?_⟩
  · exact sub_nonneg.mpr (hlower upperPoint hupperPoint)
  · intro startPoint hstartPoint endPoint hendPoint
    exact sub_le_sub (hupper startPoint hstartPoint) (hlower endPoint hendPoint)

/-- Under one-periodicity, the extrema found on `[0,1]` bound the primitive on
all of `ℝ`. -/
theorem exists_global_centeredPrimitive_extrema
    (S : Set ℝ) (hS : MeasurableSet S) (density : ℝ)
    (hindicator : Function.Periodic
      (fun t => S.indicator (fun _ => (1 : ℝ)) t) 1)
    (hdensity : (volume (S ∩ Ioc (0 : ℝ) 1)).toReal = density) :
    ∃ lowerPoint ∈ Icc (0 : ℝ) 1, ∃ upperPoint ∈ Icc (0 : ℝ) 1,
      (∀ point, centeredPrimitive S density lowerPoint ≤ centeredPrimitive S density point) ∧
      (∀ point, centeredPrimitive S density point ≤ centeredPrimitive S density upperPoint) := by
  obtain ⟨lowerPoint, hlowerPoint, upperPoint, hupperPoint, hlower, hupper⟩ :=
    exists_centeredPrimitive_extrema S hS density
  have hperiodic := centeredPrimitive_periodic S hS density hindicator hdensity
  have hreduce : ∀ point : ℝ,
      centeredPrimitive S density (Int.fract point) = centeredPrimitive S density point := by
    intro point
    simpa [Int.self_sub_floor] using
      (hperiodic.sub_int_mul_eq (x := point) ⌊point⌋)
  refine ⟨lowerPoint, hlowerPoint, upperPoint, hupperPoint, ?_, ?_⟩
  · intro point
    rw [← hreduce point]
    exact hlower (Int.fract point) ⟨Int.fract_nonneg point, (Int.fract_lt_one point).le⟩
  · intro point
    rw [← hreduce point]
    exact hupper (Int.fract point) ⟨Int.fract_nonneg point, (Int.fract_lt_one point).le⟩

/-- The rigorous upper half of THM-933's eta/q duality:
`ell * (density - eta(ell)) ≤ q` for every positive scale. No component-count
topology enters this statement or its proof. -/
theorem fixedScaleEta_deficit_le_primitiveDiscrepancy
    (S : Set ℝ) (hS : MeasurableSet S) (density : ℝ)
    (hindicator : Function.Periodic
      (fun t => S.indicator (fun _ => (1 : ℝ)) t) 1)
    (hdensity : (volume (S ∩ Ioc (0 : ℝ) 1)).toReal = density)
    {ell : ℝ} (hell : 0 < ell) :
    ell * (density - fixedScaleEta S ell) ≤ primitiveDiscrepancy S density := by
  obtain ⟨startPoint, _hstartPoint, heta, _hminimum⟩ :=
    fixedScaleEta_attained S hS hell
  obtain ⟨lowerPoint, _hlowerPoint, upperPoint, _hupperPoint, hq, hlower, hupper⟩ :=
    primitiveDiscrepancy_attained S hS density
  have hperiodic := centeredPrimitive_periodic S hS density hindicator hdensity
  have hreduce : ∀ point : ℝ,
      centeredPrimitive S density (Int.fract point) = centeredPrimitive S density point := by
    intro point
    simpa [Int.self_sub_floor] using
      (hperiodic.sub_int_mul_eq (x := point) ⌊point⌋)
  have hlowerGlobal : ∀ point, centeredPrimitive S density lowerPoint ≤
      centeredPrimitive S density point := by
    intro point
    rw [← hreduce point]
    exact hlower (Int.fract point) ⟨Int.fract_nonneg point, (Int.fract_lt_one point).le⟩
  have hupperGlobal : ∀ point, centeredPrimitive S density point ≤
      centeredPrimitive S density upperPoint := by
    intro point
    rw [← hreduce point]
    exact hupper (Int.fract point) ⟨Int.fract_nonneg point, (Int.fract_lt_one point).le⟩
  have hretained : retainedWindow S ell startPoint = fixedScaleEta S ell * ell := by
    rw [heta]
    exact (div_mul_cancel₀ _ hell.ne').symm
  have hidentity :
      retainedWindow S ell startPoint - density * ell =
        centeredPrimitive S density (startPoint + ell) -
          centeredPrimitive S density startPoint := by
    simpa [retainedWindow] using centeredPrimitive_interval_identity S hS density
      (show startPoint ≤ startPoint + ell by linarith)
  calc
    ell * (density - fixedScaleEta S ell) =
        density * ell - fixedScaleEta S ell * ell := by ring
    _ = density * ell - retainedWindow S ell startPoint := by rw [hretained]
    _ = centeredPrimitive S density startPoint -
        centeredPrimitive S density (startPoint + ell) := by linarith
    _ ≤ centeredPrimitive S density upperPoint -
        centeredPrimitive S density lowerPoint :=
      sub_le_sub (hupperGlobal startPoint) (hlowerGlobal (startPoint + ell))
    _ = primitiveDiscrepancy S density := hq.symm

/-- The reverse half of eta/q duality is attained: some scale
`0 < ell ≤ 1` has deficit exactly equal to the primitive discrepancy. The
positive-discrepancy case uses the oriented arc from a primitive maximum to a
primitive minimum; the zero case uses the full period. -/
theorem exists_fixedScaleEta_deficit_eq_primitiveDiscrepancy
    (S : Set ℝ) (hS : MeasurableSet S) (density : ℝ)
    (hindicator : Function.Periodic
      (fun t => S.indicator (fun _ => (1 : ℝ)) t) 1)
    (hdensity : (volume (S ∩ Ioc (0 : ℝ) 1)).toReal = density) :
    ∃ ell ∈ Ioc (0 : ℝ) 1,
      ell * (density - fixedScaleEta S ell) = primitiveDiscrepancy S density := by
  by_cases hqzero : primitiveDiscrepancy S density = 0
  · refine ⟨1, by norm_num, ?_⟩
    rw [fixedScaleEta_one S hS density hindicator hdensity, hqzero]
    ring
  · have hqpos : 0 < primitiveDiscrepancy S density :=
      lt_of_le_of_ne (primitiveDiscrepancy_nonneg S hS density) (Ne.symm hqzero)
    obtain ⟨lowerPoint, hlowerPoint, upperPoint, hupperPoint, hq, _hlower, _hupper⟩ :=
      primitiveDiscrepancy_attained S hS density
    have hdiffpos : 0 < centeredPrimitive S density upperPoint -
        centeredPrimitive S density lowerPoint := by
      rw [← hq]
      exact hqpos
    have hpointsNe : upperPoint ≠ lowerPoint := by
      intro hpoints
      rw [hpoints, sub_self] at hdiffpos
      exact lt_irrefl 0 hdiffpos
    have hperiodic := centeredPrimitive_periodic S hS density hindicator hdensity
    let ell := if upperPoint ≤ lowerPoint then
      lowerPoint - upperPoint else lowerPoint + 1 - upperPoint
    have hellpos : 0 < ell := by
      by_cases horder : upperPoint ≤ lowerPoint
      · simp only [ell, if_pos horder]
        exact sub_pos.mpr (lt_of_le_of_ne horder hpointsNe)
      · simp only [ell, if_neg horder]
        by_contra hnonpositive
        have hellnonpositive : lowerPoint + 1 - upperPoint ≤ 0 :=
          le_of_not_gt hnonpositive
        have hlowerZero : lowerPoint = 0 := by
          linarith [hlowerPoint.1, hupperPoint.2]
        have hupperOne : upperPoint = 1 := by
          linarith [hlowerPoint.1, hupperPoint.2]
        have hequalValues : centeredPrimitive S density upperPoint =
            centeredPrimitive S density lowerPoint := by
          rw [hlowerZero, hupperOne]
          simpa using hperiodic 0
        rw [hequalValues, sub_self] at hdiffpos
        exact lt_irrefl 0 hdiffpos
    have hellle : ell ≤ 1 := by
      by_cases horder : upperPoint ≤ lowerPoint
      · simp only [ell, if_pos horder]
        linarith [hlowerPoint.2, hupperPoint.1]
      · simp only [ell, if_neg horder]
        have hlowerUpper : lowerPoint < upperPoint := lt_of_not_ge horder
        linarith
    have hend : centeredPrimitive S density (upperPoint + ell) =
        centeredPrimitive S density lowerPoint := by
      by_cases horder : upperPoint ≤ lowerPoint
      · simp only [ell, if_pos horder]
        congr 1
        ring
      · simp only [ell, if_neg horder]
        have hadd : upperPoint + (lowerPoint + 1 - upperPoint) = lowerPoint + 1 := by ring
        rw [hadd]
        exact hperiodic lowerPoint
    have hidentity : retainedWindow S ell upperPoint - density * ell =
        centeredPrimitive S density (upperPoint + ell) -
          centeredPrimitive S density upperPoint := by
      simpa [retainedWindow] using centeredPrimitive_interval_identity S hS density
        (show upperPoint ≤ upperPoint + ell by linarith)
    have harcDeficit :
        ell * (density - retainedWindow S ell upperPoint / ell) =
          primitiveDiscrepancy S density := by
      have hcancel : retainedWindow S ell upperPoint / ell * ell =
          retainedWindow S ell upperPoint := div_mul_cancel₀ _ hellpos.ne'
      calc
        ell * (density - retainedWindow S ell upperPoint / ell) =
            density * ell - (retainedWindow S ell upperPoint / ell) * ell := by ring
        _ = density * ell - retainedWindow S ell upperPoint := by rw [hcancel]
        _ = centeredPrimitive S density upperPoint -
            centeredPrimitive S density (upperPoint + ell) := by linarith
        _ = centeredPrimitive S density upperPoint -
            centeredPrimitive S density lowerPoint := by rw [hend]
        _ = primitiveDiscrepancy S density := hq.symm
    obtain ⟨_minimumPoint, _hminimumPoint, _heta, hminimum⟩ :=
      fixedScaleEta_attained S hS hellpos
    have hreverse : primitiveDiscrepancy S density ≤
        ell * (density - fixedScaleEta S ell) := by
      rw [← harcDeficit]
      exact mul_le_mul_of_nonneg_left
        (sub_le_sub_left (hminimum upperPoint hupperPoint) density) hellpos.le
    have hforward := fixedScaleEta_deficit_le_primitiveDiscrepancy
      S hS density hindicator hdensity hellpos
    exact ⟨ell, ⟨hellpos, hellle⟩, le_antisymm hforward hreverse⟩

/-- Full attained eta/q duality from THM-933:
`q = sup_{0 < ell ≤ 1} ell * (density - eta(ell))`. -/
theorem primitiveDiscrepancy_eq_fixedScaleDeficitSup
    (S : Set ℝ) (hS : MeasurableSet S) (density : ℝ)
    (hindicator : Function.Periodic
      (fun t => S.indicator (fun _ => (1 : ℝ)) t) 1)
    (hdensity : (volume (S ∩ Ioc (0 : ℝ) 1)).toReal = density) :
    primitiveDiscrepancy S density = fixedScaleDeficitSup S density := by
  let deficit := fun ell => ell * (density - fixedScaleEta S ell)
  have hnonempty : (deficit '' Ioc (0 : ℝ) 1).Nonempty :=
    ⟨deficit 1, 1, by norm_num, rfl⟩
  have hbddAbove : BddAbove (deficit '' Ioc (0 : ℝ) 1) := by
    refine ⟨primitiveDiscrepancy S density, ?_⟩
    rintro value ⟨ell, hell, rfl⟩
    exact fixedScaleEta_deficit_le_primitiveDiscrepancy
      S hS density hindicator hdensity hell.1
  obtain ⟨ell, hell, heq⟩ :=
    exists_fixedScaleEta_deficit_eq_primitiveDiscrepancy
      S hS density hindicator hdensity
  change primitiveDiscrepancy S density = sSup (deficit '' Ioc (0 : ℝ) 1)
  apply le_antisymm
  · rw [← heq]
    exact le_csSup hbddAbove ⟨ell, hell, rfl⟩
  · exact csSup_le hnonempty (by
      rintro value ⟨scale, hscale, rfl⟩
      exact fixedScaleEta_deficit_le_primitiveDiscrepancy
        S hS density hindicator hdensity hscale.1)

end ConcretePrimitiveBridge

section ComponentSum

variable {ι : Type*}

/-- Summing a one-component local-density certificate pays one discrepancy
unit per component.  This is the algebraic content of THM-933 Lemma 2. -/
theorem local_to_component_sum
    (components : Finset ι) (length kept : ι → ℝ) (density discrepancy : ℝ)
    (hlocal : ∀ component ∈ components,
      density * length component - discrepancy ≤ kept component) :
    density * (∑ component ∈ components, length component)
        - (components.card : ℝ) * discrepancy
      ≤ ∑ component ∈ components, kept component := by
  calc
    density * (∑ component ∈ components, length component)
          - (components.card : ℝ) * discrepancy
        = ∑ component ∈ components,
            (density * length component - discrepancy) := by
              simp only [Finset.sum_sub_distrib, Finset.sum_const,
                nsmul_eq_mul, Finset.mul_sum]
    _ ≤ ∑ component ∈ components, kept component :=
      Finset.sum_le_sum fun component hcomponent => hlocal component hcomponent

/-- If the survivor has at most `complexity` components and discrepancy is
nonnegative, replace the exact cardinality debt by the advertised block
complexity debt. -/
theorem local_to_complexity_sum
    (components : Finset ι) (length kept : ι → ℝ)
    (density discrepancy complexity : ℝ)
    (hlocal : ∀ component ∈ components,
      density * length component - discrepancy ≤ kept component)
    (hcard : (components.card : ℝ) ≤ complexity)
    (hdiscrepancy : 0 ≤ discrepancy) :
    density * (∑ component ∈ components, length component)
        - complexity * discrepancy
      ≤ ∑ component ∈ components, kept component := by
  calc
    density * (∑ component ∈ components, length component)
          - complexity * discrepancy
        ≤ density * (∑ component ∈ components, length component)
          - (components.card : ℝ) * discrepancy := by
            exact sub_le_sub_left
              (mul_le_mul_of_nonneg_right hcard hdiscrepancy) _
    _ ≤ ∑ component ∈ components, kept component :=
      local_to_component_sum components length kept density discrepancy hlocal

/-- Fixed-scale sampling form (Opus S333 G1, after each component is tiled by
length-`ell` test intervals): summing the component losses pays
`eta * card * ell`.  THM-933's primitive form and this fixed-scale form are
dual certificate interfaces. -/
theorem fixedScale_sampling_sum
    (components : Finset ι) (length kept : ι → ℝ) (eta ell : ℝ)
    (hlocal : ∀ component ∈ components,
      eta * (length component - ell) ≤ kept component) :
      eta * ((∑ component ∈ components, length component)
        - (components.card : ℝ) * ell)
      ≤ ∑ component ∈ components, kept component := by
  have hsumConst :
      (∑ _component ∈ components, ell) = (components.card : ℝ) * ell := by
    rw [Finset.sum_const, nsmul_eq_mul]
  calc
    eta * ((∑ component ∈ components, length component)
          - (components.card : ℝ) * ell)
        = eta * ((∑ component ∈ components, length component)
            - ∑ _component ∈ components, ell) := by rw [hsumConst]
    _ = eta * (∑ component ∈ components,
          (length component - ell)) := by rw [Finset.sum_sub_distrib]
    _ = ∑ component ∈ components,
          eta * (length component - ell) := by rw [Finset.mul_sum]
    _ ≤ ∑ component ∈ components, kept component :=
      Finset.sum_le_sum fun component hcomponent => hlocal component hcomponent

/-- Opus's `eta * mass - loss` ledger is a valid weakening of the sharper
`eta * (mass - loss)` bound whenever `eta ≤ 1` and the loss is nonnegative. -/
theorem fixedScale_weaker_loss
    (mass boundaryLoss eta : ℝ) (heta : eta ≤ 1) (hloss : 0 ≤ boundaryLoss) :
    eta * mass - boundaryLoss ≤ eta * (mass - boundaryLoss) := by
  have hnonnegative : 0 ≤ (1 - eta) * boundaryLoss :=
    mul_nonneg (sub_nonneg.mpr heta) hloss
  calc
    eta * mass - boundaryLoss
        = eta * (mass - boundaryLoss) - (1 - eta) * boundaryLoss := by ring
    _ ≤ eta * (mass - boundaryLoss) := sub_le_self _ hnonnegative

end ComponentSum

section CircularToothSplit

open LonelyRunner.RatIntervals

/-- The number of cut-open interval pieces.  For a normalized maximal region,
this is its ordinary interval-component count. -/
def intervalComponentCount (region : Region) : ℕ := List.length region

/-- In a cut-open unit circle, the leftmost and rightmost pieces belong to one
circular component when both boundary sides occur.  The length guard prevents
the full-circle singleton `[(0, 1)]` from being incorrectly merged away. -/
def boundaryPiecesJoined (region : Region) : Bool :=
  decide (2 ≤ intervalComponentCount region) &&
    region.any (fun interval => interval.1 == 0) &&
    region.any (fun interval => interval.2 == 1)

/-- Component count of a cut-open circular interval region: merge the two
boundary pieces exactly when they are the two ends of one circular component. -/
def circularComponentCount (region : Region) : ℕ :=
  if boundaryPiecesJoined region then
    intervalComponentCount region - 1
  else
    intervalComponentCount region

theorem circularComponentCount_le_intervalComponentCount (region : Region) :
    circularComponentCount region ≤ intervalComponentCount region := by
  unfold circularComponentCount
  split <;> omega

/-- Cutting a circular component at one chart boundary creates at most one
extra interval piece. -/
theorem intervalComponentCount_le_circularComponentCount_add_one (region : Region) :
    intervalComponentCount region ≤ circularComponentCount region + 1 := by
  unfold circularComponentCount
  split <;> omega

/-- A connected circular tooth after rotating its initial endpoint to `0`.
The endpoint convention is half-open, matching `RatIntervals`. -/
structure AnchoredCircularTooth where
  width : ℚ
  width_pos : 0 < width
  width_le_one : width ≤ 1

/-- Delete an anchored circular tooth from a cut-open interval region. -/
def deleteAnchoredTooth (region : Region) (tooth : AnchoredCircularTooth) : Region :=
  diff1F region (0, tooth.width)

/-- `deleteAnchoredTooth` is exact set subtraction in the cut-open chart. -/
theorem mem_deleteAnchoredTooth {x : ℚ} {region : Region}
    {tooth : AnchoredCircularTooth} :
    mem x (deleteAnchoredTooth region tooth) ↔
      mem x region ∧ ¬ ((0 : ℚ) ≤ x ∧ x < tooth.width) := by
  simpa [deleteAnchoredTooth] using
    (mem_diff1F (x := x) (L := region) (q := ((0 : ℚ), tooth.width)))

/-- Subtracting an interval whose left endpoint is the chart boundary leaves
at most one live piece of each input interval. -/
theorem cutF_anchored_length_le_one
    (interval : ℚ × ℚ) (tooth : AnchoredCircularTooth)
    (hleft : 0 ≤ interval.1) :
    List.length (cutF interval (0, tooth.width)) ≤ 1 := by
  have hdead : ¬ interval.1 < min interval.2 0 := by
    exact not_lt.mpr (le_trans (min_le_right _ _) hleft)
  have hfilter := List.length_filter_le
    (fun other : ℚ × ℚ => decide (other.1 < other.2))
    [(max interval.1 tooth.width, interval.2)]
  simpa [cutF, cut, hdead] using hfilter

/-- Anchored tooth deletion cannot increase the cut-open piece count. -/
theorem intervalComponentCount_deleteAnchoredTooth_le
    (region : Region) (tooth : AnchoredCircularTooth)
    (hleft : ∀ interval ∈ region, 0 ≤ interval.1) :
    intervalComponentCount (deleteAnchoredTooth region tooth) ≤
      intervalComponentCount region := by
  induction region with
  | nil => simp [intervalComponentCount, deleteAnchoredTooth, diff1F]
  | cons interval region inductionHypothesis =>
      have hhead : 0 ≤ interval.1 := hleft interval (List.mem_cons_self ..)
      have htail : ∀ other ∈ region, 0 ≤ other.1 := by
        intro other hother
        exact hleft other (List.mem_cons_of_mem _ hother)
      calc
        intervalComponentCount
              (deleteAnchoredTooth (interval :: region) tooth)
            = List.length (cutF interval (0, tooth.width)) +
                intervalComponentCount (deleteAnchoredTooth region tooth) := by
                  simp [intervalComponentCount, deleteAnchoredTooth, diff1F]
        _ ≤ 1 + intervalComponentCount region :=
          Nat.add_le_add (cutF_anchored_length_le_one interval tooth hhead)
            (inductionHypothesis htail)
        _ = intervalComponentCount (interval :: region) := by
          simp [intervalComponentCount, Nat.add_comm]

/-- **Concrete one-circular-tooth split.**  Rotate the new connected tooth so
it starts at the cut.  The old circular survivor can gain one cut-open piece;
anchored subtraction gains none; closing the chart can only merge pieces.
Hence deleting one circular tooth raises component count by at most one. -/
theorem circularComponentCount_deleteAnchoredTooth_le_add_one
    (region : Region) (tooth : AnchoredCircularTooth)
    (hleft : ∀ interval ∈ region, 0 ≤ interval.1) :
    circularComponentCount (deleteAnchoredTooth region tooth) ≤
      circularComponentCount region + 1 := by
  calc
    circularComponentCount (deleteAnchoredTooth region tooth)
        ≤ intervalComponentCount (deleteAnchoredTooth region tooth) :=
      circularComponentCount_le_intervalComponentCount _
    _ ≤ intervalComponentCount region :=
      intervalComponentCount_deleteAnchoredTooth_le region tooth hleft
    _ ≤ circularComponentCount region + 1 :=
      intervalComponentCount_le_circularComponentCount_add_one region

/-- The complement of the first anchored tooth in the full circle has at most
one circular component. -/
theorem firstCircularTooth_componentCount_le_one (tooth : AnchoredCircularTooth) :
    circularComponentCount (deleteAnchoredTooth [(0, 1)] tooth) ≤ 1 := by
  calc
    circularComponentCount (deleteAnchoredTooth [(0, 1)] tooth)
        ≤ intervalComponentCount (deleteAnchoredTooth [(0, 1)] tooth) :=
      circularComponentCount_le_intervalComponentCount _
    _ ≤ intervalComponentCount [(0, 1)] := by
      apply intervalComponentCount_deleteAnchoredTooth_le
      simp
    _ = 1 := by simp [intervalComponentCount]

/-- A concrete cut-chart atlas for successive circular tooth deletions.
`chart n` is obtained by rotating the survivor after `n` deletions so that the
next tooth starts at `0`.  The only topology-facing fields are that this
recharting preserves circular component count and lists normalized unit-chart
pieces; deletion itself is the exact `RatIntervals.diff1F` operation. -/
structure CircularToothAtlas where
  survivor : ℕ → Region
  chart : ℕ → Region
  tooth : ℕ → AnchoredCircularTooth
  first_chart : chart 0 = [(0, 1)]
  chart_normalized : ∀ toothCount, Norm (chart toothCount)
  chart_in_unit : ∀ toothCount interval, interval ∈ chart toothCount →
    0 ≤ interval.1 ∧ interval.2 ≤ 1
  chart_count : ∀ toothCount,
    circularComponentCount (chart toothCount) =
      circularComponentCount (survivor toothCount)
  delete_eq : ∀ toothCount,
    survivor (toothCount + 1) =
      deleteAnchoredTooth (chart toothCount) (tooth toothCount)

/-! ### Concrete rational rotation charts -/

/-- Every listed interval is live. -/
def RegionLive (region : Region) : Prop :=
  ∀ interval ∈ region, interval.1 < interval.2

/-- Every listed interval lies in the cut-open unit chart. -/
def RegionInUnit (region : Region) : Prop :=
  ∀ interval ∈ region, 0 ≤ interval.1 ∧ interval.2 ≤ 1

theorem regionLive_cutF (interval toothInterval : ℚ × ℚ) :
    RegionLive (cutF interval toothInterval) := by
  intro piece hpiece
  exact decide_eq_true_eq.mp (List.mem_filter.mp hpiece).2

theorem regionLive_diff1F (region : Region) (toothInterval : ℚ × ℚ) :
    RegionLive (diff1F region toothInterval) := by
  intro piece hpiece
  unfold diff1F at hpiece
  rw [List.mem_flatMap] at hpiece
  obtain ⟨interval, _, hpiece⟩ := hpiece
  exact regionLive_cutF interval toothInterval piece hpiece

theorem regionLive_deleteAnchoredTooth
    (region : Region) (tooth : AnchoredCircularTooth) :
    RegionLive (deleteAnchoredTooth region tooth) := by
  exact regionLive_diff1F region (0, tooth.width)

/-- Every filtered cut piece stays inside its source interval. -/
theorem cutF_piece_bounds {interval toothInterval piece : ℚ × ℚ}
    (hpiece : piece ∈ cutF interval toothInterval) :
    interval.1 ≤ piece.1 ∧ piece.2 ≤ interval.2 := by
  unfold cutF cut at hpiece
  rw [List.mem_filter] at hpiece
  rcases List.mem_cons.mp hpiece.1 with rfl | hright
  · exact ⟨le_rfl, min_le_left _ _⟩
  · rcases List.mem_singleton.mp hright with rfl
    exact ⟨le_max_left _ _, le_rfl⟩

theorem regionInUnit_diff1F
    (region : Region) (toothInterval : ℚ × ℚ)
    (hunit : RegionInUnit region) :
    RegionInUnit (diff1F region toothInterval) := by
  intro piece hpiece
  unfold diff1F at hpiece
  rw [List.mem_flatMap] at hpiece
  obtain ⟨interval, hinterval, hpiece⟩ := hpiece
  have hsource := hunit interval hinterval
  have hbounds := cutF_piece_bounds hpiece
  exact ⟨le_trans hsource.1 hbounds.1, le_trans hbounds.2 hsource.2⟩

theorem regionInUnit_deleteAnchoredTooth
    (region : Region) (tooth : AnchoredCircularTooth)
    (hunit : RegionInUnit region) :
    RegionInUnit (deleteAnchoredTooth region tooth) := by
  exact regionInUnit_diff1F region (0, tooth.width) hunit

/-- A positive anchored tooth removes the chart origin: every surviving piece
starts strictly to its right. -/
theorem deleteAnchoredTooth_piece_start_pos
    (region : Region) (tooth : AnchoredCircularTooth)
    (hunit : RegionInUnit region) {piece : ℚ × ℚ}
    (hpiece : piece ∈ deleteAnchoredTooth region tooth) :
    0 < piece.1 := by
  unfold deleteAnchoredTooth diff1F at hpiece
  rw [List.mem_flatMap] at hpiece
  obtain ⟨interval, hinterval, hpiece⟩ := hpiece
  unfold cutF cut at hpiece
  rw [List.mem_filter] at hpiece
  obtain ⟨hpiece, hlive⟩ := hpiece
  rcases List.mem_cons.mp hpiece with rfl | hpiece
  · have hleft := (hunit interval hinterval).1
    have hlive' := decide_eq_true_eq.mp hlive
    have hright := min_le_right interval.2 0
    exfalso
    linarith
  · rcases List.mem_singleton.mp hpiece with rfl
    exact lt_of_lt_of_le tooth.width_pos (le_max_right _ _)

/-- `wrapOne` always lands inside the unit chart when its input width is at
most one. -/
theorem wrapOne_in_unit {interval piece : ℚ × ℚ}
    (hwidth : interval.2 - interval.1 ≤ 1)
    (hpiece : piece ∈ wrapOne interval) :
    0 ≤ piece.1 ∧ piece.2 ≤ 1 := by
  unfold wrapOne at hpiece
  set floorValue : ℤ := ⌊interval.1⌋ with hfloorValue
  have hfloor : (floorValue : ℚ) ≤ interval.1 := by
    simpa [hfloorValue] using Int.floor_le interval.1
  have hfloorNext : interval.1 < (floorValue : ℚ) + 1 := by
    rw [hfloorValue]
    exact Int.lt_floor_add_one interval.1
  by_cases hend : interval.2 - (floorValue : ℚ) ≤ 1
  · simp only [hend, if_true] at hpiece
    rcases List.mem_singleton.mp hpiece with rfl
    exact ⟨by linarith, hend⟩
  · simp only [hend, if_false] at hpiece
    rcases List.mem_cons.mp hpiece with rfl | hpiece
    · exact ⟨by linarith, le_rfl⟩
    · rcases List.mem_singleton.mp hpiece with rfl
      exact ⟨le_rfl, by linarith⟩

theorem regionInUnit_translateCirc
    (shift : ℚ) (region : Region) (hunit : RegionInUnit region) :
    RegionInUnit (translateCirc shift region) := by
  intro piece hpiece
  unfold translateCirc wrap translate at hpiece
  rw [List.mem_flatMap] at hpiece
  obtain ⟨translatedInterval, htranslatedInterval, hpiece⟩ := hpiece
  rw [List.mem_map] at htranslatedInterval
  obtain ⟨interval, hinterval, rfl⟩ := htranslatedInterval
  apply wrapOne_in_unit (piece := piece) (hpiece := hpiece)
  have hbounds := hunit interval hinterval
  linarith

theorem wrapOne_live {interval piece : ℚ × ℚ}
    (hlive : interval.1 < interval.2)
    (hpiece : piece ∈ wrapOne interval) :
    piece.1 < piece.2 := by
  unfold wrapOne at hpiece
  set floorValue : ℤ := ⌊interval.1⌋ with hfloorValue
  have hfloorNext : interval.1 < (floorValue : ℚ) + 1 := by
    rw [hfloorValue]
    exact Int.lt_floor_add_one interval.1
  by_cases hend : interval.2 - (floorValue : ℚ) ≤ 1
  · simp only [hend, if_true] at hpiece
    rcases List.mem_singleton.mp hpiece with rfl
    linarith
  · have hend' : 1 < interval.2 - (floorValue : ℚ) := lt_of_not_ge hend
    simp only [hend, if_false] at hpiece
    rcases List.mem_cons.mp hpiece with rfl | hpiece
    · linarith
    · rcases List.mem_singleton.mp hpiece with rfl
      linarith

theorem regionLive_translateCirc
    (shift : ℚ) (region : Region) (hlive : RegionLive region) :
    RegionLive (translateCirc shift region) := by
  intro piece hpiece
  unfold translateCirc wrap translate at hpiece
  rw [List.mem_flatMap] at hpiece
  obtain ⟨translatedInterval, htranslatedInterval, hpiece⟩ := hpiece
  rw [List.mem_map] at htranslatedInterval
  obtain ⟨interval, hinterval, rfl⟩ := htranslatedInterval
  apply wrapOne_live (hpiece := hpiece)
  simpa using hlive interval hinterval

/-- Sort a raw circle translation by left endpoint.  Sorting is semantic only:
it neither changes membership nor the boundary-aware component count. -/
def sortedTranslateCirc (shift : ℚ) (region : Region) : Region :=
  (translateCirc shift region).insertionSort
    (fun left right => left.1 ≤ right.1)

theorem sortedTranslateCirc_perm (shift : ℚ) (region : Region) :
    (sortedTranslateCirc shift region).Perm (translateCirc shift region) := by
  unfold sortedTranslateCirc
  exact List.perm_insertionSort _ _

theorem mem_eq_of_region_perm {left right : Region} (hperm : left.Perm right)
    (x : ℚ) : mem x left ↔ mem x right := by
  unfold mem
  constructor
  · rintro ⟨interval, hinterval, hx⟩
    exact ⟨interval, hperm.mem_iff.mp hinterval, hx⟩
  · rintro ⟨interval, hinterval, hx⟩
    exact ⟨interval, hperm.mem_iff.mpr hinterval, hx⟩

theorem intervalComponentCount_eq_of_perm {left right : Region}
    (hperm : left.Perm right) :
    intervalComponentCount left = intervalComponentCount right := by
  exact hperm.length_eq

theorem boundaryPiecesJoined_eq_of_perm {left right : Region}
    (hperm : left.Perm right) :
    boundaryPiecesJoined left = boundaryPiecesJoined right := by
  unfold boundaryPiecesJoined intervalComponentCount
  rw [hperm.length_eq, hperm.any_eq, hperm.any_eq]

theorem circularComponentCount_eq_of_perm {left right : Region}
    (hperm : left.Perm right) :
    circularComponentCount left = circularComponentCount right := by
  unfold circularComponentCount
  rw [boundaryPiecesJoined_eq_of_perm hperm,
    intervalComponentCount_eq_of_perm hperm]

theorem regionInUnit_sortedTranslateCirc
    (shift : ℚ) (region : Region) (hunit : RegionInUnit region) :
    RegionInUnit (sortedTranslateCirc shift region) := by
  intro interval hinterval
  exact regionInUnit_translateCirc shift region hunit interval
    ((sortedTranslateCirc_perm shift region).mem_iff.mp hinterval)

theorem regionLive_sortedTranslateCirc
    (shift : ℚ) (region : Region) (hlive : RegionLive region) :
    RegionLive (sortedTranslateCirc shift region) := by
  intro interval hinterval
  exact regionLive_translateCirc shift region hlive interval
    ((sortedTranslateCirc_perm shift region).mem_iff.mp hinterval)

/-- Linear pieces are separated when either one ends before the other starts. -/
def IntervalsSeparated (left right : ℚ × ℚ) : Prop :=
  left.2 ≤ right.1 ∨ right.2 ≤ left.1

def RegionSeparated (region : Region) : Prop :=
  region.Pairwise IntervalsSeparated

theorem intervalsSeparated_symmetric {left right : ℚ × ℚ} :
    IntervalsSeparated left right → IntervalsSeparated right left := by
  intro hseparated
  rcases hseparated with h | h
  · exact Or.inr h
  · exact Or.inl h

/-- Live, start-sorted, separated pieces form a `Norm` region. -/
theorem norm_of_live_pairwiseStart_separated : ∀ {region : Region},
    RegionLive region →
    region.Pairwise (fun left right => left.1 ≤ right.1) →
    RegionSeparated region → Norm region := by
  intro region
  induction region with
  | nil => intro _ _ _; trivial
  | cons interval region inductionHypothesis =>
      intro hlive hstart hseparated
      have hinterval : interval.1 < interval.2 :=
        hlive interval (List.mem_cons_self ..)
      cases region with
      | nil => exact hinterval
      | cons next remaining =>
          have htailLive : RegionLive (next :: remaining) := by
            intro piece hpiece
            exact hlive piece (List.mem_cons_of_mem _ hpiece)
          have htailStart := (List.pairwise_cons.mp hstart).2
          have htailSeparated := (List.pairwise_cons.mp hseparated).2
          have htailNorm := inductionHypothesis htailLive htailStart htailSeparated
          have hstartNext :=
            (List.pairwise_cons.mp hstart).1 next (List.mem_cons_self ..)
          have hseparatedNext :=
            (List.pairwise_cons.mp hseparated).1 next (List.mem_cons_self ..)
          have hnextLive : next.1 < next.2 :=
            htailLive next (List.mem_cons_self ..)
          have hendBefore : interval.2 ≤ next.1 := by
            rcases hseparatedNext with h | h
            · exact h
            · exfalso
              linarith
          exact ⟨hinterval, hendBefore, htailNorm⟩

theorem norm_sortedTranslateCirc_of_separated
    (shift : ℚ) (region : Region) (hlive : RegionLive region)
    (hseparated : RegionSeparated (translateCirc shift region)) :
    Norm (sortedTranslateCirc shift region) := by
  apply norm_of_live_pairwiseStart_separated
  · exact regionLive_sortedTranslateCirc shift region hlive
  · exact List.pairwise_insertionSort _ _
  · unfold RegionSeparated at hseparated ⊢
    exact hseparated.perm (sortedTranslateCirc_perm shift region).symm
      intervalsSeparated_symmetric

/-- Exact membership transport for a rational circle translation. -/
theorem mem_translateCirc_iff {x shift : ℚ} {region : Region}
    (hx0 : 0 ≤ x) (hx1 : x < 1)
    (hwidth : ∀ interval ∈ region, interval.2 - interval.1 ≤ 1) :
    mem x (translateCirc shift region) ↔
      ∃ integerShift : ℤ, mem (x + integerShift - shift) region := by
  unfold translateCirc
  have htranslatedWidth : ∀ interval ∈ translate shift region,
      interval.2 - interval.1 ≤ 1 := by
    intro interval hinterval
    unfold translate at hinterval
    rw [List.mem_map] at hinterval
    obtain ⟨source, hsource, rfl⟩ := hinterval
    simpa using hwidth source hsource
  rw [mem_wrap hx0 hx1 htranslatedWidth]
  constructor
  · rintro ⟨integerShift, hmember⟩
    exact ⟨integerShift, mem_translate.mp hmember⟩
  · rintro ⟨integerShift, hmember⟩
    exact ⟨integerShift, mem_translate.mpr hmember⟩

theorem mem_sortedTranslateCirc_iff {x shift : ℚ} {region : Region}
    (hx0 : 0 ≤ x) (hx1 : x < 1)
    (hwidth : ∀ interval ∈ region, interval.2 - interval.1 ≤ 1) :
    mem x (sortedTranslateCirc shift region) ↔
      ∃ integerShift : ℤ, mem (x + integerShift - shift) region := by
  rw [mem_eq_of_region_perm (sortedTranslateCirc_perm shift region)]
  exact mem_translateCirc_iff hx0 hx1 hwidth

/-- The one-bit correction converting cut-open pieces into circular
components. -/
def boundaryCorrection (region : Region) : ℕ :=
  if boundaryPiecesJoined region then 1 else 0

theorem circularComponentCount_add_boundaryCorrection (region : Region) :
    circularComponentCount region + boundaryCorrection region =
      intervalComponentCount region := by
  by_cases hjoined : boundaryPiecesJoined region = true
  · have htwo : 2 ≤ intervalComponentCount region := by
      have hjoined' := hjoined
      simp [boundaryPiecesJoined] at hjoined'
      exact hjoined'.1.1
    simp [circularComponentCount, boundaryCorrection, hjoined]
    omega
  · simp [circularComponentCount, boundaryCorrection, hjoined]

theorem boundaryPiecesJoined_deleteAnchoredTooth_eq_false
    (region : Region) (tooth : AnchoredCircularTooth)
    (hunit : RegionInUnit region) :
    boundaryPiecesJoined (deleteAnchoredTooth region tooth) = false := by
  have hnoLeft :
      (deleteAnchoredTooth region tooth).any
        (fun interval => interval.1 == 0) = false := by
    rw [List.any_eq_false]
    intro interval hinterval
    simp only [beq_iff_eq]
    exact ne_of_gt
      (deleteAnchoredTooth_piece_start_pos region tooth hunit hinterval)
  simp [boundaryPiecesJoined, hnoLeft]

/-- The exact combinatorial obligation for a rotation: a new cut either splits
one component or removes the old boundary merge, with those two effects
balancing. -/
def RotationCutBalance (shift : ℚ) (region : Region) : Prop :=
  intervalComponentCount (translateCirc shift region) +
      boundaryCorrection region =
    intervalComponentCount region +
      boundaryCorrection (translateCirc shift region)

/-- The complete topology certificate needed for one rational rechart.  The
separation field makes sorting a `Norm` presentation; `cutBalance` proves that
the new chart cut neither creates nor loses a circular component. -/
structure RotationTopologyCertificate (shift : ℚ) (region : Region) : Prop where
  separated : RegionSeparated (translateCirc shift region)
  cutBalance : RotationCutBalance shift region

/-- Positive-stage form of the topology certificate.  Since the preceding
anchored deletion removes the origin, the old boundary correction is zero;
the only count equation left says that translation creates one extra linear
piece exactly when the new chart has a boundary merge. -/
structure PositiveRotationTopologyCertificate
    (shift : ℚ) (region : Region) : Prop where
  separated : RegionSeparated (translateCirc shift region)
  pieceBalance :
    intervalComponentCount (translateCirc shift region) =
      intervalComponentCount region +
        boundaryCorrection (translateCirc shift region)

theorem rotationCutBalance_of_positiveCertificate
    (shift : ℚ) (region : Region)
    (hboundary : boundaryPiecesJoined region = false)
    (hcertificate : PositiveRotationTopologyCertificate shift region) :
    RotationCutBalance shift region := by
  unfold RotationCutBalance
  simp [boundaryCorrection, hboundary, hcertificate.pieceBalance]

theorem circularComponentCount_translateCirc_eq_of_cutBalance
    (shift : ℚ) (region : Region)
    (hbalance : RotationCutBalance shift region) :
    circularComponentCount (translateCirc shift region) =
      circularComponentCount region := by
  have htranslated :=
    circularComponentCount_add_boundaryCorrection (translateCirc shift region)
  have hsource := circularComponentCount_add_boundaryCorrection region
  unfold RotationCutBalance at hbalance
  omega

theorem circularComponentCount_sortedTranslateCirc_eq_of_cutBalance
    (shift : ℚ) (region : Region)
    (hbalance : RotationCutBalance shift region) :
    circularComponentCount (sortedTranslateCirc shift region) =
      circularComponentCount region := by
  calc
    circularComponentCount (sortedTranslateCirc shift region)
        = circularComponentCount (translateCirc shift region) :=
      circularComponentCount_eq_of_perm (sortedTranslateCirc_perm shift region)
    _ = circularComponentCount region :=
      circularComponentCount_translateCirc_eq_of_cutBalance shift region hbalance

/-- Stage zero uses the full circle in the coordinate system where the first
tooth is already anchored.  Later stages use sorted rational circle
translations. -/
def rationalCircleRechart (shift : ℕ → ℚ) (toothCount : ℕ)
    (region : Region) : Region :=
  match toothCount with
  | 0 => [(0, 1)]
  | count + 1 => sortedTranslateCirc (shift (count + 1)) region

/-- Actual rational circle survivors in their successive tooth-anchored
coordinate frames. -/
def rationalCircleSurvivor
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth) : ℕ → Region
  | 0 => [(0, 1)]
  | count + 1 =>
      deleteAnchoredTooth
        (rationalCircleRechart shift count
          (rationalCircleSurvivor shift tooth count))
        (tooth count)

/-- The normalized chart used immediately before the next tooth deletion. -/
def rationalCircleChart
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) : Region :=
  rationalCircleRechart shift toothCount
    (rationalCircleSurvivor shift tooth toothCount)

theorem rationalCircleChart_zero
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth) :
    rationalCircleChart shift tooth 0 = [(0, 1)] := rfl

theorem rationalCircleSurvivor_succ
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    rationalCircleSurvivor shift tooth (toothCount + 1) =
      deleteAnchoredTooth (rationalCircleChart shift tooth toothCount)
        (tooth toothCount) := rfl

theorem regionLive_rationalCircleSurvivor
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth) :
    ∀ toothCount, RegionLive (rationalCircleSurvivor shift tooth toothCount) := by
  intro toothCount
  induction toothCount with
  | zero => simp [RegionLive, rationalCircleSurvivor]
  | succ toothCount _ =>
      rw [rationalCircleSurvivor_succ]
      exact regionLive_deleteAnchoredTooth _ _

theorem regionInUnit_rationalCircleChart_of_survivor
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ)
    (hunit : RegionInUnit (rationalCircleSurvivor shift tooth toothCount)) :
    RegionInUnit (rationalCircleChart shift tooth toothCount) := by
  cases toothCount with
  | zero => simp [RegionInUnit, rationalCircleChart, rationalCircleRechart]
  | succ toothCount =>
      exact regionInUnit_sortedTranslateCirc _ _ hunit

theorem regionInUnit_rationalCircleSurvivor
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth) :
    ∀ toothCount, RegionInUnit (rationalCircleSurvivor shift tooth toothCount) := by
  intro toothCount
  induction toothCount with
  | zero => simp [RegionInUnit, rationalCircleSurvivor]
  | succ toothCount inductionHypothesis =>
      rw [rationalCircleSurvivor_succ]
      exact regionInUnit_deleteAnchoredTooth _ _
        (regionInUnit_rationalCircleChart_of_survivor shift tooth toothCount
          inductionHypothesis)

theorem boundaryPiecesJoined_rationalCircleSurvivor_eq_false
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) (hpositive : 1 ≤ toothCount) :
    boundaryPiecesJoined (rationalCircleSurvivor shift tooth toothCount) = false := by
  obtain ⟨previousCount, rfl⟩ := Nat.exists_eq_add_of_le hpositive
  rw [Nat.add_comm 1 previousCount, rationalCircleSurvivor_succ]
  apply boundaryPiecesJoined_deleteAnchoredTooth_eq_false
  exact regionInUnit_rationalCircleChart_of_survivor shift tooth previousCount
    (regionInUnit_rationalCircleSurvivor shift tooth previousCount)

theorem regionInUnit_rationalCircleChart
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    RegionInUnit (rationalCircleChart shift tooth toothCount) :=
  regionInUnit_rationalCircleChart_of_survivor shift tooth toothCount
    (regionInUnit_rationalCircleSurvivor shift tooth toothCount)

theorem rationalCircleChart_count_preserved
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (hbalance : ∀ toothCount, 1 ≤ toothCount →
      RotationCutBalance (shift toothCount)
        (rationalCircleSurvivor shift tooth toothCount)) :
    ∀ toothCount,
      circularComponentCount (rationalCircleChart shift tooth toothCount) =
        circularComponentCount
          (rationalCircleSurvivor shift tooth toothCount) := by
  intro toothCount
  cases toothCount with
  | zero => rfl
  | succ toothCount =>
      apply circularComponentCount_sortedTranslateCirc_eq_of_cutBalance
      exact hbalance (toothCount + 1) (by omega)

/-- **Concrete atlas constructor.**  The survivors, circle translations,
sorting, component-count transport, and exact deletions are now definitions or
theorems.  Consumers supply only the two genuine normalization facts for each
positive stage: sorted translated pieces are `Norm`, and the boundary
cut-balance law holds. -/
def rationalCircularToothAtlas
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (htopology : ∀ toothCount, 1 ≤ toothCount →
      RotationTopologyCertificate (shift toothCount)
        (rationalCircleSurvivor shift tooth toothCount)) :
    CircularToothAtlas where
  survivor := rationalCircleSurvivor shift tooth
  chart := rationalCircleChart shift tooth
  tooth := tooth
  first_chart := rationalCircleChart_zero shift tooth
  chart_normalized := by
    intro toothCount
    cases toothCount with
    | zero =>
        norm_num [rationalCircleChart, rationalCircleRechart,
          LonelyRunner.RatIntervals.Norm]
    | succ toothCount =>
        apply norm_sortedTranslateCirc_of_separated
        · exact regionLive_rationalCircleSurvivor shift tooth (toothCount + 1)
        · exact (htopology (toothCount + 1) (by omega)).separated
  chart_in_unit := by
    intro toothCount interval hinterval
    exact regionInUnit_rationalCircleChart shift tooth toothCount interval hinterval
  chart_count := rationalCircleChart_count_preserved shift tooth
    (fun toothCount htoothCount =>
      (htopology toothCount htoothCount).cutBalance)
  delete_eq := rationalCircleSurvivor_succ shift tooth

/-- Constructor using the reduced positive-stage certificate. -/
def rationalCircularToothAtlasOfPositiveCertificates
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (htopology : ∀ toothCount, 1 ≤ toothCount →
      PositiveRotationTopologyCertificate (shift toothCount)
        (rationalCircleSurvivor shift tooth toothCount)) :
    CircularToothAtlas :=
  rationalCircularToothAtlas shift tooth fun toothCount hpositive =>
    { separated := (htopology toothCount hpositive).separated
      cutBalance := rotationCutBalance_of_positiveCertificate
        (shift toothCount) (rationalCircleSurvivor shift tooth toothCount)
        (boundaryPiecesJoined_rationalCircleSurvivor_eq_false
          shift tooth toothCount hpositive)
        (htopology toothCount hpositive) }

/-- The recursive survivor really removes the next anchored tooth, pointwise. -/
theorem mem_rationalCircleSurvivor_succ
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) (x : ℚ) :
    mem x (rationalCircleSurvivor shift tooth (toothCount + 1)) ↔
      mem x (rationalCircleChart shift tooth toothCount) ∧
        ¬ ((0 : ℚ) ≤ x ∧ x < (tooth toothCount).width) := by
  rw [rationalCircleSurvivor_succ]
  exact mem_deleteAnchoredTooth

/-- At every positive stage, the sorted chart has the exact circle-rotation
semantics of `RatIntervals.translateCirc`. -/
theorem mem_rationalCircleChart_succ_iff
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) {x : ℚ} (hx0 : 0 ≤ x) (hx1 : x < 1) :
    mem x (rationalCircleChart shift tooth (toothCount + 1)) ↔
      ∃ integerShift : ℤ,
        mem (x + integerShift - shift (toothCount + 1))
          (rationalCircleSurvivor shift tooth (toothCount + 1)) := by
  apply mem_sortedTranslateCirc_iff hx0 hx1
  intro interval hinterval
  have hunit := regionInUnit_rationalCircleSurvivor shift tooth
    (toothCount + 1) interval hinterval
  linarith

end CircularToothSplit

section ComponentCap

/-- Abstract geometric ledger behind THM-933 Lemma 3.  The complement of one
deleted tooth has at most one component; if every additional tooth raises the
count by at most one, the complement of `toothCount` teeth has at most
`toothCount` components. -/
theorem component_count_le_tooth_count
    (componentCount : ℕ → ℕ)
    (hfirst : componentCount 1 ≤ 1)
    (hstep : ∀ toothCount, 1 ≤ toothCount →
      componentCount (toothCount + 1) ≤ componentCount toothCount + 1)
    (toothCount : ℕ) :
    1 ≤ toothCount → componentCount toothCount ≤ toothCount := by
  induction toothCount with
  | zero => omega
  | succ previousCount inductionHypothesis =>
      intro _
      by_cases hzero : previousCount = 0
      · subst previousCount
        simpa using hfirst
      · have hpositive : 1 ≤ previousCount := Nat.one_le_iff_ne_zero.mpr hzero
        calc
          componentCount (previousCount + 1)
              ≤ componentCount previousCount + 1 :=
                hstep previousCount hpositive
          _ ≤ previousCount + 1 :=
            Nat.add_le_add_right (inductionHypothesis hpositive) 1

/-- Concrete circular-interval instantiation of
`component_count_le_tooth_count`.  All count growth is discharged by
`circularComponentCount_deleteAnchoredTooth_le_add_one`; consumers only supply
the cut-chart atlas connecting their actual circle survivor to the rational
interval presentation. -/
theorem circular_component_count_le_tooth_count
    (atlas : CircularToothAtlas) (toothCount : ℕ) :
    1 ≤ toothCount →
      circularComponentCount (atlas.survivor toothCount) ≤ toothCount := by
  apply component_count_le_tooth_count
      (fun count => circularComponentCount (atlas.survivor count))
  · rw [atlas.delete_eq 0, atlas.first_chart]
    exact firstCircularTooth_componentCount_le_one (atlas.tooth 0)
  · intro previousCount _
    rw [atlas.delete_eq previousCount]
    calc
      circularComponentCount
            (deleteAnchoredTooth (atlas.chart previousCount)
              (atlas.tooth previousCount))
          ≤ circularComponentCount (atlas.chart previousCount) + 1 := by
            apply circularComponentCount_deleteAnchoredTooth_le_add_one
            intro interval hinterval
            exact (atlas.chart_in_unit previousCount interval hinterval).1
      _ = circularComponentCount (atlas.survivor previousCount) + 1 := by
        rw [atlas.chart_count]

/-- End-to-end component cap for the concrete rational survivor recursion. -/
theorem rational_circle_component_count_le_tooth_count
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (htopology : ∀ count, 1 ≤ count →
      PositiveRotationTopologyCertificate (shift count)
        (rationalCircleSurvivor shift tooth count))
    (toothCount : ℕ) :
    1 ≤ toothCount →
      circularComponentCount
          (rationalCircleSurvivor shift tooth toothCount) ≤ toothCount :=
  circular_component_count_le_tooth_count
    (rationalCircularToothAtlasOfPositiveCertificates shift tooth htopology)
    toothCount

end ComponentCap

section Recurrence

/-- Recursive lower ledger: multiply the previous lower bound by the next
local density and subtract the next boundary debt. -/
def lowerBound (initial : ℝ) (density debt : ℕ → ℝ) : ℕ → ℝ
  | 0 => initial
  | n + 1 => density n * lowerBound initial density debt n - debt n

/-- Recursive accumulated debt.  Earlier debts acquire every later density
factor, exactly as in the closed THM-933 formula. -/
def weightedDebt (density debt : ℕ → ℝ) : ℕ → ℝ
  | 0 => 0
  | n + 1 => density n * weightedDebt density debt n + debt n

/-- Product of all density factors used through stage `n`. -/
def densityProduct (density : ℕ → ℝ) (n : ℕ) : ℝ :=
  ∏ k ∈ Finset.range n, density k

/-- The recursive ledger is exactly product times the initial mass minus the
later-density-weighted debt. -/
theorem lowerBound_eq_product_sub_weightedDebt
    (initial : ℝ) (density debt : ℕ → ℝ) (n : ℕ) :
    lowerBound initial density debt n =
      densityProduct density n * initial - weightedDebt density debt n := by
  induction n with
  | zero => simp [lowerBound, densityProduct, weightedDebt]
  | succ n ih =>
      rw [lowerBound, weightedDebt, ih]
      simp only [densityProduct, Finset.prod_range_succ]
      ring

/-- Closed suffix-product form of the debt ledger.  The debt charged at stage
`k` is multiplied by every density factor from `k+1` through `n-1`. -/
theorem weightedDebt_eq_suffix_sum
    (density debt : ℕ → ℝ) (n : ℕ) :
    weightedDebt density debt n =
      ∑ k ∈ Finset.Ico 0 n,
        debt k * ∏ j ∈ Finset.Ico (k + 1) n, density j := by
  induction n with
  | zero => simp [weightedDebt]
  | succ n ih =>
      rw [weightedDebt, ih, Finset.sum_Ico_succ_top (Nat.zero_le n),
        Finset.Ico_self, Finset.prod_empty, mul_one]
      refine congr_arg (· + debt n) ?_
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro k hk
      rw [Finset.prod_Ico_succ_top (by
        have := Finset.mem_Ico.mp hk
        omega)]
      ring

/-- Fully unrolled product-minus-suffix-debt formula used in THM-933 (BG). -/
theorem lowerBound_eq_closed
    (initial : ℝ) (density debt : ℕ → ℝ) (n : ℕ) :
    lowerBound initial density debt n =
      densityProduct density n * initial
        - ∑ k ∈ Finset.Ico 0 n,
            debt k * ∏ j ∈ Finset.Ico (k + 1) n, density j := by
  rw [lowerBound_eq_product_sub_weightedDebt,
    weightedDebt_eq_suffix_sum]

/-- Soundness of the entire gluing recurrence.  Nonnegative density factors
allow every previously proved lower bound to be multiplied forward. -/
theorem lowerBound_le_actual
    (initial : ℝ) (density debt actual : ℕ → ℝ)
    (hinitial : initial ≤ actual 0)
    (hdensity : ∀ n, 0 ≤ density n)
    (hstep : ∀ n, density n * actual n - debt n ≤ actual (n + 1)) :
    ∀ n, lowerBound initial density debt n ≤ actual n := by
  intro n
  induction n with
  | zero => simpa [lowerBound] using hinitial
  | succ n ih =>
      calc
        lowerBound initial density debt (n + 1)
            = density n * lowerBound initial density debt n - debt n := rfl
        _ ≤ density n * actual n - debt n := by
          exact sub_le_sub_right (mul_le_mul_of_nonneg_left ih (hdensity n)) _
        _ ≤ actual (n + 1) := hstep n

/-- The explicit two-transition / three-block form consumed by THM-933:
`d₁d₂d₃ - (e₂d₃+e₃)`. -/
theorem three_block_gluing
    (d₁ d₂ d₃ e₂ e₃ w₁ w₂ w₃ : ℝ)
    (hd₂ : 0 ≤ d₂) (hd₃ : 0 ≤ d₃)
    (h₁ : d₁ ≤ w₁)
    (h₂ : d₂ * w₁ - e₂ ≤ w₂)
    (h₃ : d₃ * w₂ - e₃ ≤ w₃) :
    d₁ * d₂ * d₃ - (e₂ * d₃ + e₃) ≤ w₃ := by
  have h₂' : d₂ * d₁ - e₂ ≤ w₂ := by
    calc
      d₂ * d₁ - e₂ ≤ d₂ * w₁ - e₂ := by
        exact sub_le_sub_right (mul_le_mul_of_nonneg_left h₁ hd₂) _
      _ ≤ w₂ := h₂
  calc
    d₁ * d₂ * d₃ - (e₂ * d₃ + e₃) = d₃ * (d₂ * d₁ - e₂) - e₃ := by ring
    _ ≤ d₃ * w₂ - e₃ := by
      exact sub_le_sub_right (mul_le_mul_of_nonneg_left h₂' hd₃) _
    _ ≤ w₃ := h₃

end Recurrence

section ExactArithmetic

/-- Exact arithmetic behind the improved `R = 7` LRC(14) corollary. -/
theorem six_pow_twelve_gt_seven_pow_eleven :
    7 ^ 11 < 6 ^ 12 := by norm_num

/-- The explicit positive measure floor in THM-933 equation (10). -/
theorem ratio_seven_margin_pos :
    (0 : ℚ) < 199455593 / 13841287201 := by norm_num

/-- Kernel check of the three-block product-minus-debt identity. -/
theorem exact_three_block_ledger :
    (53 / 105 : ℚ) * (3 / 5) * (386 / 735)
        - (193 / 8575 + 193 / 6174)
      = 81253 / 771750 := by norm_num

/-- Kernel check of the stronger ledger using exact prefix component counts
`10` and `132` instead of tooth caps `15` and `375`. -/
theorem exact_three_block_component_ledger :
    (53 / 105 : ℚ) * (3 / 5) * (386 / 735)
        - (386 / 25725 + 4246 / 385875)
      = 7334 / 55125 := by norm_num

/-- The exact THM-933 block-gluing floor is positive. -/
theorem exact_three_block_margin_pos :
    (0 : ℚ) < 81253 / 771750 := by norm_num

/-- Exact arithmetic for the sharp fixed-scale Opus S333 7+6 witness. -/
theorem opus_fixedScale_two_block_ledger :
    (274025490881738650 / 1001359472502594621 : ℚ)
        * (12208485893419843 / 38882590600758450 - 1724 / 45000)
      = 60354211840721383388269695262412
        / 800043501647462161192289496375975 := by norm_num

/-- The fixed-scale 7+6 composition margin is strictly positive. -/
theorem opus_fixedScale_two_block_margin_pos :
    (0 : ℚ) < 60354211840721383388269695262412
      / 800043501647462161192289496375975 := by norm_num

/-- The direct exact safe measure dominates the proved block-gluing floor. -/
theorem exact_measure_dominates_ledger :
    (81253 / 771750 : ℚ) ≤ 55063 / 330750 := by norm_num

end ExactArithmetic

/-! ## Axiom audit -/

#print axioms local_to_component_sum
#print axioms local_to_complexity_sum
#print axioms fixedScale_sampling_sum
#print axioms primitive_interval_lower
#print axioms primitive_interval_sharp
#print axioms fixedScale_deficit_le_discrepancy
#print axioms fixedScale_extremizer_eq_discrepancy
#print axioms intervalIndicatorIntegral_eq_measure
#print axioms centeredPrimitive_interval_identity
#print axioms centeredPrimitive_Icc_identity
#print axioms continuous_centeredPrimitive
#print axioms centeredPrimitive_periodic
#print axioms retainedWindow_one_eq_density
#print axioms continuous_retainedWindow
#print axioms fixedScaleEta_attained
#print axioms fixedScaleEta_one
#print axioms primitiveDiscrepancy_attained
#print axioms primitiveDiscrepancy_nonneg
#print axioms exists_global_centeredPrimitive_extrema
#print axioms fixedScaleEta_deficit_le_primitiveDiscrepancy
#print axioms exists_fixedScaleEta_deficit_eq_primitiveDiscrepancy
#print axioms primitiveDiscrepancy_eq_fixedScaleDeficitSup
#print axioms fixedScale_weaker_loss
#print axioms mem_deleteAnchoredTooth
#print axioms circularComponentCount_deleteAnchoredTooth_le_add_one
#print axioms mem_rationalCircleChart_succ_iff
#print axioms rationalCircleChart_count_preserved
#print axioms component_count_le_tooth_count
#print axioms circular_component_count_le_tooth_count
#print axioms rational_circle_component_count_le_tooth_count
#print axioms lowerBound_eq_product_sub_weightedDebt
#print axioms weightedDebt_eq_suffix_sum
#print axioms lowerBound_eq_closed
#print axioms lowerBound_le_actual
#print axioms three_block_gluing
#print axioms six_pow_twelve_gt_seven_pow_eleven
#print axioms ratio_seven_margin_pos
#print axioms exact_three_block_ledger
#print axioms exact_three_block_component_ledger
#print axioms exact_three_block_margin_pos
#print axioms opus_fixedScale_two_block_ledger
#print axioms opus_fixedScale_two_block_margin_pos
#print axioms exact_measure_dominates_ledger

end LRC14.LocalDensityBlockGluing
