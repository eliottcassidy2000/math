import TournamentH7.LRCB5NormalizedBridge
import TournamentH7.LRCCommensuration
import TournamentH7.LRCPairContinuumBridge

namespace LonelyRunner
namespace LRCB5PairGridBridge

open Finset MeasureTheory
open LRC14Concrete LRCB5NormalizedBridge
open LRCB5CleanModulus LRCPairContinuumBridge LRCPairRatioLayerArithmetic

noncomputable section

/-- The joint-failure count on the complete `q`-point grid, including `p = 0`. -/
def allGridJointFail (v : Fin 13 → ℤ) (q : ℕ) (T : Finset (Fin 13)) : ℕ :=
  ((Finset.range q).filter fun p => ∀ i ∈ T, ¬ inBand v q p i).card

/-- The complete-grid and nonzero-grid counts differ by the always-bad point `p = 0`. -/
theorem allGridJointFail_eq_jointFail_add_one
    (v : Fin 13 → ℤ) (q : ℕ) (T : Finset (Fin 13))
    (hq : 0 < q) :
    allGridJointFail v q T = jointFail v q T + 1 := by
  have hrange : Finset.range q = insert 0 (Finset.Ioo 0 q) := by
    ext p
    simp only [Finset.mem_range, Finset.mem_insert, Finset.mem_Ioo]
    omega
  have hzero : ∀ i ∈ T, ¬ inBand v q 0 i := by
    intro i hi
    unfold inBand
    simp only [Int.ofNat_zero, mul_zero, Int.zero_emod]
    omega
  unfold allGridJointFail jointFail
  rw [hrange]
  rw [Finset.filter_insert, if_pos hzero]
  rw [Finset.card_insert_of_notMem]
  simp

/-- The nonzero-grid density of a pair-support failure event. -/
def pairGridDensity (v : Fin 13 → ℤ) (q : ℕ) (T : Finset (Fin 13)) : ℝ :=
  (jointFail v q T : ℝ) / ((q : ℝ) - 1)

/-- The complete-grid density of a pair-support failure event. -/
def allPairGridDensity (v : Fin 13 → ℤ) (q : ℕ) (T : Finset (Fin 13)) : ℝ :=
  (allGridJointFail v q T : ℝ) / q

/-- The continuous joint danger event on the parallel-class circle. -/
def pairDanger (v : Fin 13 → ℤ) (T : Finset (Fin 13)) : Set UnitAddCircle :=
  {x | ∀ i ∈ T,
    x ∈ LRCCommensuration.danger (v i) (0 : UnitAddCircle) (1 / 14)}

/-- Continuous pair correlation on one period. -/
def pairContinuumCorrelation (v : Fin 13 → ℤ) (T : Finset (Fin 13)) : ℝ :=
  (volume (pairDanger v T)).toReal

/-- The continuum pair mass `C(v)`: sum of pair correlations minus equilibrium. -/
def continuumMass2 (v : Fin 13 → ℤ) : ℝ :=
  ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
    (pairContinuumCorrelation v T - 1 / 49)

/-- Total speed mass in the real normalization used by the discrepancy bound. -/
def speedL1 (v : Fin 13 → ℤ) : ℝ :=
  ∑ i, |(v i : ℝ)|

/-- The endpoint budget attached to one pair support. -/
def pairCutWeight (v : Fin 13 → ℤ) (T : Finset (Fin 13)) : ℝ :=
  2 * ∑ i ∈ T, |(v i : ℝ)|

theorem pairContinuumCorrelation_nonneg
    (v : Fin 13 → ℤ) (T : Finset (Fin 13)) :
    0 ≤ pairContinuumCorrelation v T := by
  exact ENNReal.toReal_nonneg

theorem pairContinuumCorrelation_le_one
    (v : Fin 13 → ℤ) (T : Finset (Fin 13)) :
    pairContinuumCorrelation v T ≤ 1 := by
  unfold pairContinuumCorrelation
  have hmeasure : volume (pairDanger v T) ≤ volume (Set.univ : Set UnitAddCircle) :=
    measure_mono (Set.subset_univ _)
  have hreal := ENNReal.toReal_mono (measure_ne_top volume (Set.univ : Set UnitAddCircle)) hmeasure
  simpa using hreal

/-- Exact normalized pair ledger in real coordinates. -/
theorem normalizedMass2_real_eq_pairGridLedger
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q) :
    (normalizedMass2 v q : ℝ) =
      ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
        (pairGridDensity v q T - 1 / 49) := by
  have hdenQ : ((q : ℚ) - 1) ≠ 0 := by
    have : (1 : ℚ) < q := by exact_mod_cast hq
    linarith
  have hrat : normalizedMass2 v q =
      ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
        ((jointFail v q T : ℚ) / ((q : ℚ) - 1) - 1 / 49) := by
    unfold normalizedMass2 normalizedAggregateDeviation aggregateDeviation
    simp_rw [Finset.sum_div]
    apply Finset.sum_congr rfl
    intro T hT
    unfold LRC14Concrete.deviation
    rw [(Finset.mem_powersetCard.mp hT).2]
    norm_num
    field_simp [hdenQ]
  rw [hrat]
  push_cast
  rfl

/-- Each runner lies in exactly twelve two-element supports. -/
theorem pair_incidence_sum (f : Fin 13 → ℝ) :
    (∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
      ∑ i ∈ T, f i) = 12 * ∑ i, f i := by
  classical
  let pairs := (Finset.univ : Finset (Fin 13)).powersetCard 2
  have hcard : ∀ i : Fin 13, (pairs.filter fun T => i ∈ T).card = 12 := by
    intro i
    have h := Finset.card_filter_powersetCard_subset
      ({i} : Finset (Fin 13)) (Finset.univ : Finset (Fin 13)) 2
      (by simp) (by simp)
    simpa [pairs, Finset.singleton_subset_iff] using h
  calc
    (∑ T ∈ pairs, ∑ i ∈ T, f i) =
        ∑ T ∈ pairs, ∑ i, if i ∈ T then f i else 0 := by
          apply Finset.sum_congr rfl
          intro T hT
          rw [← Finset.sum_filter]
          simp
    _ = ∑ i, ∑ T ∈ pairs, if i ∈ T then f i else 0 := Finset.sum_comm
    _ = ∑ i, 12 * f i := by
          apply Finset.sum_congr rfl
          intro i hi
          rw [← Finset.sum_filter, Finset.sum_const, hcard i, nsmul_eq_mul]
          norm_num
    _ = 12 * ∑ i, f i := by rw [Finset.mul_sum]

theorem pair_support_card :
    ((Finset.univ : Finset (Fin 13)).powersetCard 2).card = 78 := by
  rw [Finset.card_powersetCard]
  norm_num [Nat.choose]

/-- Algebraic cost of deleting the zero grid point and changing `q` to `q-1`. -/
theorem remove_zero_normalization_discrepancy
    (q allCount nonzeroCount : ℕ) (correlation endpointCost : ℝ)
    (hq : 1 < q) (hcount : allCount = nonzeroCount + 1)
    (hcorrelation0 : 0 ≤ correlation) (hcorrelation1 : correlation ≤ 1)
    (hendpointCost : 0 ≤ endpointCost)
    (hraw : |(allCount : ℝ) / q - correlation| ≤ endpointCost / q) :
    |(nonzeroCount : ℝ) / ((q : ℝ) - 1) - correlation| ≤
      (endpointCost + 1) / ((q : ℝ) - 1) := by
  have hqR : (0 : ℝ) < q := by positivity
  have hqOne : (1 : ℝ) < q := by exact_mod_cast hq
  have hden : (0 : ℝ) < (q : ℝ) - 1 := by linarith
  have hid : (nonzeroCount : ℝ) / ((q : ℝ) - 1) - correlation =
      ((q : ℝ) / ((q : ℝ) - 1)) * ((allCount : ℝ) / q - correlation) +
        (correlation - 1) / ((q : ℝ) - 1) := by
    rw [hcount]
    push_cast
    field_simp [ne_of_gt hqR, ne_of_gt hden]
    ring
  rw [hid]
  calc
    |((q : ℝ) / ((q : ℝ) - 1)) * ((allCount : ℝ) / q - correlation) +
        (correlation - 1) / ((q : ℝ) - 1)|
        ≤ |((q : ℝ) / ((q : ℝ) - 1)) *
            ((allCount : ℝ) / q - correlation)| +
          |(correlation - 1) / ((q : ℝ) - 1)| := abs_add_le _ _
    _ = ((q : ℝ) / ((q : ℝ) - 1)) *
          |(allCount : ℝ) / q - correlation| +
          (1 - correlation) / ((q : ℝ) - 1) := by
            rw [abs_mul, abs_of_pos (div_pos hqR hden), abs_div]
            rw [abs_of_pos hden, abs_of_nonpos (sub_nonpos.mpr hcorrelation1)]
            ring
    _ ≤ ((q : ℝ) / ((q : ℝ) - 1)) * (endpointCost / q) +
          (1 - correlation) / ((q : ℝ) - 1) := by
            gcongr
    _ ≤ ((q : ℝ) / ((q : ℝ) - 1)) * (endpointCost / q) +
          1 / ((q : ℝ) - 1) := by
            have : 1 - correlation ≤ 1 := by linarith
            gcongr
    _ = (endpointCost + 1) / ((q : ℝ) - 1) := by
          field_simp [ne_of_gt hqR, ne_of_gt hden]

/-- The one remaining geometric obligation: complete-grid discrepancy is bounded by
the number of circle cut points. -/
def RawPairGridDiscrepancyAt (v : Fin 13 → ℤ) (q : ℕ) : Prop :=
  ∀ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
    |allPairGridDensity v q T - pairContinuumCorrelation v T| ≤
      pairCutWeight v T / q

/-- The raw complete-grid estimate gives the exact nonzero-grid pair estimate. -/
theorem pairGridDiscrepancy_of_raw
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q)
    (hraw : RawPairGridDiscrepancyAt v q) :
    ∀ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
      |pairGridDensity v q T - pairContinuumCorrelation v T| ≤
        (pairCutWeight v T + 1) / ((q : ℝ) - 1) := by
  intro T hT
  apply remove_zero_normalization_discrepancy q
    (allGridJointFail v q T) (jointFail v q T)
    (pairContinuumCorrelation v T) (pairCutWeight v T) hq
  · exact allGridJointFail_eq_jointFail_add_one v q T (by omega)
  · exact pairContinuumCorrelation_nonneg v T
  · exact pairContinuumCorrelation_le_one v T
  · unfold pairCutWeight
    positivity
  · exact hraw T hT

/-- Integer count assigned to one half-open interval by its two ceiling cuts. -/
def ceilIntervalCount (q : ℕ) (left right : ℝ) : ℕ :=
  (⌈right * q⌉ - ⌈left * q⌉).toNat

/-- A single half-open component costs at most `2/q` in complete-grid
discrepancy.  The sharp ceiling calculation actually gives `< 1/q`; the
factor two deliberately absorbs endpoint conventions. -/
theorem ceilIntervalCount_discrepancy
    (q : ℕ) (hq : 0 < q) (left right : ℝ) (hlr : left ≤ right) :
    |(ceilIntervalCount q left right : ℝ) / q - (right - left)| ≤ 2 / q := by
  let leftCut : ℤ := ⌈left * q⌉
  let rightCut : ℤ := ⌈right * q⌉
  have hqR : (0 : ℝ) < q := by exact_mod_cast hq
  have hcuts : leftCut ≤ rightCut := by
    apply Int.ceil_mono
    exact mul_le_mul_of_nonneg_right hlr (Nat.cast_nonneg q)
  have hdiff : 0 ≤ rightCut - leftCut := sub_nonneg.mpr hcuts
  have hcount : (ceilIntervalCount q left right : ℝ) =
      (rightCut : ℝ) - (leftCut : ℝ) := by
    unfold ceilIntervalCount leftCut rightCut
    have hnat := Int.toNat_of_nonneg hdiff
    have hnatR := congrArg (fun value : ℤ => (value : ℝ)) hnat
    push_cast at hnatR
    exact hnatR
  have hleftLo : left * q ≤ (leftCut : ℝ) := by
    exact Int.le_ceil _
  have hleftHi : (leftCut : ℝ) < left * q + 1 := by
    exact Int.ceil_lt_add_one _
  have hrightLo : right * q ≤ (rightCut : ℝ) := by
    exact Int.le_ceil _
  have hrightHi : (rightCut : ℝ) < right * q + 1 := by
    exact Int.ceil_lt_add_one _
  let error : ℝ := (rightCut : ℝ) - (leftCut : ℝ) - q * (right - left)
  have herror : |error| ≤ 2 := by
    rw [abs_le]
    dsimp [error]
    constructor <;> linarith
  have hid : (ceilIntervalCount q left right : ℝ) / q - (right - left) =
      error / q := by
    rw [hcount]
    dsimp [error]
    field_simp [ne_of_gt hqR]
  rw [hid, abs_div, abs_of_pos hqR]
  exact div_le_div_of_nonneg_right herror hqR.le

/-- Actual complete-grid count in a strict real interval. -/
def openIntervalGridCount (q : ℕ) (left right : ℝ) : ℕ :=
  ((Finset.Ico (0 : ℤ) (q : ℤ)).filter fun (point : ℤ) =>
    left < (point : ℝ) / q ∧ (point : ℝ) / q < right).card

/-- A strict interval inside `[0,1]` has complete-grid discrepancy at most
`1/q`.  This is the endpoint lemma needed for the open LRC danger arcs. -/
theorem openIntervalGridCount_discrepancy
    (q : ℕ) (hq : 0 < q) (left right : ℝ)
    (hleft : 0 ≤ left) (hright : right ≤ 1) (hlr : left ≤ right) :
    |(openIntervalGridCount q left right : ℝ) / q - (right - left)| ≤ 1 / q := by
  let leftCut : ℤ := ⌊left * q⌋
  let rightCut : ℤ := ⌈right * q⌉
  have hqR : (0 : ℝ) < q := by exact_mod_cast hq
  have hset :
      (Finset.Ico (0 : ℤ) (q : ℤ)).filter (fun (point : ℤ) =>
        left < (point : ℝ) / q ∧ (point : ℝ) / q < right) =
        Finset.Ioo leftCut rightCut := by
    ext point
    simp only [Finset.mem_filter, Finset.mem_Ico, Finset.mem_Ioo]
    constructor
    · rintro ⟨hpoint, hleftPoint, hpointRight⟩
      have hleftScaled : left * q < (point : ℝ) := by
        exact (lt_div_iff₀ hqR).mp hleftPoint
      have hrightScaled : (point : ℝ) < right * q := by
        exact (div_lt_iff₀ hqR).mp hpointRight
      exact ⟨Int.floor_lt.mpr hleftScaled, Int.lt_ceil.mpr hrightScaled⟩
    · rintro ⟨hleftPoint, hpointRight⟩
      have hleftScaled : left * q < (point : ℝ) := Int.floor_lt.mp hleftPoint
      have hrightScaled : (point : ℝ) < right * q := Int.lt_ceil.mp hpointRight
      have hleftCutNonneg : 0 ≤ leftCut := by
        apply Int.floor_nonneg.mpr
        positivity
      have hrightCutLe : rightCut ≤ (q : ℤ) := by
        apply Int.ceil_le.mpr
        push_cast
        exact mul_le_of_le_one_left (Nat.cast_nonneg q) hright
      refine ⟨⟨by omega, by omega⟩, ?_, ?_⟩
      · exact (lt_div_iff₀ hqR).mpr hleftScaled
      · exact (div_lt_iff₀ hqR).mpr hrightScaled
  unfold openIntervalGridCount
  rw [hset, Int.card_Ioo]
  let difference : ℤ := rightCut - leftCut - 1
  by_cases hdifference : 0 ≤ difference
  · have hcount : ((difference.toNat : ℕ) : ℝ) = (difference : ℝ) := by
      have hnat := Int.toNat_of_nonneg hdifference
      have hnatR := congrArg (fun value : ℤ => (value : ℝ)) hnat
      push_cast at hnatR
      exact hnatR
    have hleftLo : (leftCut : ℝ) ≤ left * q := by
      exact Int.floor_le _
    have hleftHi : left * q < (leftCut : ℝ) + 1 := by
      exact Int.lt_floor_add_one _
    have hrightLo : right * q ≤ (rightCut : ℝ) := by
      exact Int.le_ceil _
    have hrightHi : (rightCut : ℝ) < right * q + 1 := by
      exact Int.ceil_lt_add_one _
    let error : ℝ :=
      (rightCut : ℝ) - (leftCut : ℝ) - 1 - q * (right - left)
    have herror : |error| ≤ 1 := by
      rw [abs_le]
      dsimp [error]
      constructor <;> linarith
    have hid : (difference.toNat : ℝ) / q - (right - left) = error / q := by
      rw [hcount]
      dsimp [difference, error]
      push_cast
      field_simp [ne_of_gt hqR]
    rw [hid, abs_div, abs_of_pos hqR]
    exact div_le_div_of_nonneg_right herror hqR.le
  · have hdifferenceNeg : difference < 0 := lt_of_not_ge hdifference
    have hcuts : rightCut ≤ leftCut := by
      dsimp [difference] at hdifferenceNeg
      omega
    have hleftLo : (leftCut : ℝ) ≤ left * q := Int.floor_le _
    have hrightLo : right * q ≤ (rightCut : ℝ) := Int.le_ceil _
    have heqScaled : left * q = right * q := by
      have hlrScaled := mul_le_mul_of_nonneg_right hlr (Nat.cast_nonneg q)
      have hcutsR : (rightCut : ℝ) ≤ (leftCut : ℝ) := by exact_mod_cast hcuts
      linarith
    have heq : left = right := by
      apply mul_right_cancel₀ (ne_of_gt hqR)
      simpa [mul_comm] using heqScaled
    have htoNat : difference.toNat = 0 := Int.toNat_of_nonpos (le_of_lt hdifferenceNeg)
    rw [htoNat, heq]
    norm_num

/-- Exact interval-ledger interface for one pair.  Constructing this ledger
is precisely the remaining component theorem: at most `|vᵢ|+|vⱼ|`
half-open components, with their grid count and Lebesgue mass both exposed. -/
structure PairIntervalLedger
    (v : Fin 13 → ℤ) (q : ℕ) (T : Finset (Fin 13)) where
  componentCount : ℕ
  left : Fin componentCount → ℝ
  right : Fin componentCount → ℝ
  left_le_right : ∀ component, left component ≤ right component
  componentCount_le : componentCount ≤ ∑ i ∈ T, (v i).natAbs
  allGrid_count_eq :
    allGridJointFail v q T =
      ∑ component, ceilIntervalCount q (left component) (right component)
  continuum_eq :
    pairContinuumCorrelation v T =
      ∑ component, (right component - left component)

/-- A component ledger proves the raw pair discrepancy. -/
theorem rawPairGridDiscrepancy_of_intervalLedger
    (v : Fin 13 → ℤ) (q : ℕ) (T : Finset (Fin 13)) (hq : 0 < q)
    (ledger : PairIntervalLedger v q T) :
    |allPairGridDensity v q T - pairContinuumCorrelation v T| ≤
      pairCutWeight v T / q := by
  have hqR : (0 : ℝ) < q := by exact_mod_cast hq
  have hid : allPairGridDensity v q T - pairContinuumCorrelation v T =
      ∑ component,
        ((ceilIntervalCount q (ledger.left component) (ledger.right component) : ℝ) / q -
          (ledger.right component - ledger.left component)) := by
    unfold allPairGridDensity
    rw [ledger.allGrid_count_eq, ledger.continuum_eq]
    push_cast
    rw [Finset.sum_div]
    simp_rw [Finset.sum_sub_distrib]
  rw [hid]
  calc
    |∑ component,
        ((ceilIntervalCount q (ledger.left component) (ledger.right component) : ℝ) / q -
          (ledger.right component - ledger.left component))|
        ≤ ∑ component,
            |(ceilIntervalCount q (ledger.left component) (ledger.right component) : ℝ) / q -
              (ledger.right component - ledger.left component)| :=
          Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _component : Fin ledger.componentCount, (2 : ℝ) / (q : ℝ) := by
          apply Finset.sum_le_sum
          intro component hcomponent
          exact ceilIntervalCount_discrepancy q hq
            (ledger.left component) (ledger.right component)
            (ledger.left_le_right component)
    _ = (ledger.componentCount : ℝ) * ((2 : ℝ) / (q : ℝ)) := by
          rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    _ ≤ (∑ i ∈ T, |(v i : ℝ)|) * ((2 : ℝ) / (q : ℝ)) := by
          have hmass : ((∑ i ∈ T, (v i).natAbs : ℕ) : ℝ) =
              ∑ i ∈ T, |(v i : ℝ)| := by
            push_cast
            simp
          have hcomponent : (ledger.componentCount : ℝ) ≤
              ∑ i ∈ T, |(v i : ℝ)| := by
            rw [← hmass]
            exact_mod_cast ledger.componentCount_le
          gcongr
    _ = pairCutWeight v T / q := by
          unfold pairCutWeight
          field_simp [ne_of_gt hqR]

/-- The exact remaining geometric payload, isolated pair by pair. -/
def PairIntervalLedgersAt (v : Fin 13 → ℤ) (q : ℕ) : Prop :=
  ∀ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
    Nonempty (PairIntervalLedger v q T)

theorem rawPairGridDiscrepancy_of_intervalLedgers
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    (hledgers : PairIntervalLedgersAt v q) :
    RawPairGridDiscrepancyAt v q := by
  intro T hT
  exact rawPairGridDiscrepancy_of_intervalLedger v q T hq (hledgers T hT).some

/-- Strict-open linear component ledger.  Cutting circular components at zero

can at most double their number, so the honest cap is
`2 * (|vᵢ|+|vⱼ|)`; the sharp open-interval lemma charges only `1/q` each. -/
structure OpenPairIntervalLedger
    (v : Fin 13 → ℤ) (q : ℕ) (T : Finset (Fin 13)) where
  componentCount : ℕ
  left : Fin componentCount → ℝ
  right : Fin componentCount → ℝ
  left_nonneg : ∀ component, 0 ≤ left component
  right_le_one : ∀ component, right component ≤ 1
  left_le_right : ∀ component, left component ≤ right component
  componentCount_le : componentCount ≤ 2 * ∑ i ∈ T, (v i).natAbs
  allGrid_count_eq :
    allGridJointFail v q T =
      ∑ component, openIntervalGridCount q (left component) (right component)
  continuum_eq :
    pairContinuumCorrelation v T =
      ∑ component, (right component - left component)

/-- The strongest elementary raw endpoint lemma in this scratch file: an
exact strict-open component ledger directly yields the desired
`2(|vᵢ|+|vⱼ|)/q` complete-grid discrepancy. -/
theorem rawPairGridDiscrepancy_of_openIntervalLedger
    (v : Fin 13 → ℤ) (q : ℕ) (T : Finset (Fin 13)) (hq : 0 < q)
    (ledger : OpenPairIntervalLedger v q T) :
    |allPairGridDensity v q T - pairContinuumCorrelation v T| ≤
      pairCutWeight v T / q := by
  have hqR : (0 : ℝ) < q := by exact_mod_cast hq
  have hid : allPairGridDensity v q T - pairContinuumCorrelation v T =
      ∑ component,
        ((openIntervalGridCount q (ledger.left component) (ledger.right component) : ℝ) / q -
          (ledger.right component - ledger.left component)) := by
    unfold allPairGridDensity
    rw [ledger.allGrid_count_eq, ledger.continuum_eq]
    push_cast
    rw [Finset.sum_div]
    simp_rw [Finset.sum_sub_distrib]
  rw [hid]
  calc
    |∑ component,
        ((openIntervalGridCount q (ledger.left component) (ledger.right component) : ℝ) / q -
          (ledger.right component - ledger.left component))|
        ≤ ∑ component,
            |(openIntervalGridCount q (ledger.left component) (ledger.right component) : ℝ) / q -
              (ledger.right component - ledger.left component)| :=
          Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _component : Fin ledger.componentCount, (1 : ℝ) / (q : ℝ) := by
          apply Finset.sum_le_sum
          intro component hcomponent
          exact openIntervalGridCount_discrepancy q hq
            (ledger.left component) (ledger.right component)
            (ledger.left_nonneg component) (ledger.right_le_one component)
            (ledger.left_le_right component)
    _ = (ledger.componentCount : ℝ) * ((1 : ℝ) / (q : ℝ)) := by
          rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    _ ≤ (2 * ∑ i ∈ T, |(v i : ℝ)|) * ((1 : ℝ) / (q : ℝ)) := by
          have hmass : ((∑ i ∈ T, (v i).natAbs : ℕ) : ℝ) =
              ∑ i ∈ T, |(v i : ℝ)| := by
            push_cast
            simp
          have hcomponent : (ledger.componentCount : ℝ) ≤
              2 * ∑ i ∈ T, |(v i : ℝ)| := by
            rw [← hmass]
            exact_mod_cast ledger.componentCount_le
          gcongr
    _ = pairCutWeight v T / q := by
          unfold pairCutWeight
          field_simp [ne_of_gt hqR]

/-- Named remaining geometric obligation in its strict-open form. -/
def OpenPairIntervalLedgersAt (v : Fin 13 → ℤ) (q : ℕ) : Prop :=
  ∀ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
    Nonempty (OpenPairIntervalLedger v q T)

theorem rawPairGridDiscrepancy_of_openIntervalLedgers
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    (hledgers : OpenPairIntervalLedgersAt v q) :
    RawPairGridDiscrepancyAt v q := by
  intro T hT
  exact rawPairGridDiscrepancy_of_openIntervalLedger v q T hq (hledgers T hT).some

/-- The full elementary aggregation bridge.  Once the concrete cut-point discrepancy
is supplied pairwise, `normalizedMass2` differs from the continuum mass by at most
`(24 * sumAbs + 78)/(q-1)`. -/
theorem normalizedMass2_sub_continuumMass2_le
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q)
    (hraw : RawPairGridDiscrepancyAt v q) :
    |(normalizedMass2 v q : ℝ) - continuumMass2 v| ≤
      (24 * speedL1 v + 78) / ((q : ℝ) - 1) := by
  let pairs := (Finset.univ : Finset (Fin 13)).powersetCard 2
  have hpair := pairGridDiscrepancy_of_raw v q hq hraw
  have hid : (normalizedMass2 v q : ℝ) - continuumMass2 v =
      ∑ T ∈ pairs, (pairGridDensity v q T - pairContinuumCorrelation v T) := by
    rw [normalizedMass2_real_eq_pairGridLedger v q hq]
    unfold continuumMass2
    rw [← Finset.sum_sub_distrib]
    apply Finset.sum_congr rfl
    intro T hT
    ring
  rw [hid]
  calc
    |∑ T ∈ pairs, (pairGridDensity v q T - pairContinuumCorrelation v T)|
        ≤ ∑ T ∈ pairs,
            |pairGridDensity v q T - pairContinuumCorrelation v T| :=
          Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ T ∈ pairs, (pairCutWeight v T + 1) / ((q : ℝ) - 1) := by
          apply Finset.sum_le_sum
          intro T hT
          exact hpair T hT
    _ = (24 * speedL1 v + 78) / ((q : ℝ) - 1) := by
          rw [← Finset.sum_div]
          congr 1
          rw [Finset.sum_add_distrib]
          unfold pairCutWeight
          rw [← Finset.mul_sum]
          rw [pair_incidence_sum]
          rw [Finset.sum_const, pair_support_card, nsmul_eq_mul]
          unfold speedL1
          ring


/-- End-to-end component-ledger form of the target estimate. -/
theorem normalizedMass2_sub_continuumMass2_le_of_openIntervalLedgers
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q)
    (hledgers : OpenPairIntervalLedgersAt v q) :
    |(normalizedMass2 v q : ℝ) - continuumMass2 v| ≤
      (24 * speedL1 v + 78) / ((q : ℝ) - 1) :=
  normalizedMass2_sub_continuumMass2_le v q hq
    (rawPairGridDiscrepancy_of_openIntervalLedgers v q (by omega) hledgers)


/-- Real and rational presentations of the aggregate error budget agree. -/
theorem pairGridErrorBudget_cast (v : Fin 13 → ℤ) (q : ℕ) :
    (pairGridErrorBudget v q : ℝ) =
      (24 * speedL1 v + 78) / ((q : ℝ) - 1) := by
  unfold pairGridErrorBudget speedL1
  push_cast
  congr 3
  apply Finset.sum_congr rfl
  intro i _
  simp

/-- THM-954's continuum floor and the raw complete-grid discrepancy imply the
strict normalized pair socket at the explicit clean height `534`. -/
theorem normalizedMass2_clean_534_gt_neg_target_of_continuum_floor
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hcontinuum :
      -(negativePairTierBoundPathOnly : ℝ) ≤ continuumMass2 v)
    (hraw : RawPairGridDiscrepancyAt v (cleanModulus v 534)) :
    -(13 / 30 : ℚ) < normalizedMass2 v (cleanModulus v 534) := by
  have hq := one_lt_cleanModulus v 534 hv
  have hdiscrepancy := normalizedMass2_sub_continuumMass2_le
    v (cleanModulus v 534) hq hraw
  rw [← pairGridErrorBudget_cast] at hdiscrepancy
  have hlower :
      -(pairGridErrorBudget v (cleanModulus v 534) : ℝ) ≤
        (normalizedMass2 v (cleanModulus v 534) : ℝ) - continuumMass2 v :=
    neg_le_of_abs_le hdiscrepancy
  have hmarginQ := cleanModulus_534_error_lt_path_margin v hv
  have hmargin :
      (pairGridErrorBudget v (cleanModulus v 534) : ℝ) <
        (13 / 30 - negativePairTierBoundPathOnly : ℚ) := by
    exact_mod_cast hmarginQ
  have hreal :
      -(13 / 30 : ℝ) <
        (normalizedMass2 v (cleanModulus v 534) : ℝ) := by
    push_cast at hmargin
    nlinarith
  have hcast :
      ((-(13 / 30 : ℚ) : ℚ) : ℝ) <
        ((normalizedMass2 v (cleanModulus v 534) : ℚ) : ℝ) := by
    norm_num at hreal ⊢
    exact hreal
  exact (Rat.cast_lt (K := ℝ)).mp hcast

/-- Half-open component ledgers are a concrete sufficient payload for the
strict clean-modulus pair socket. -/
theorem normalizedMass2_clean_534_gt_neg_target_of_intervalLedgers
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hcontinuum :
      -(negativePairTierBoundPathOnly : ℝ) ≤ continuumMass2 v)
    (hledgers : PairIntervalLedgersAt v (cleanModulus v 534)) :
    -(13 / 30 : ℚ) < normalizedMass2 v (cleanModulus v 534) := by
  apply normalizedMass2_clean_534_gt_neg_target_of_continuum_floor
    v hv hcontinuum
  exact rawPairGridDiscrepancy_of_intervalLedgers
    v (cleanModulus v 534)
      (Nat.zero_lt_of_lt (one_lt_cleanModulus v 534 hv)) hledgers

/-- Strict-open component ledgers match the actual danger-arc endpoint
convention and supply the same clean-modulus pair socket. -/
theorem normalizedMass2_clean_534_gt_neg_target_of_openIntervalLedgers
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hcontinuum :
      -(negativePairTierBoundPathOnly : ℝ) ≤ continuumMass2 v)
    (hledgers : OpenPairIntervalLedgersAt v (cleanModulus v 534)) :
    -(13 / 30 : ℚ) < normalizedMass2 v (cleanModulus v 534) := by
  apply normalizedMass2_clean_534_gt_neg_target_of_continuum_floor
    v hv hcontinuum
  exact rawPairGridDiscrepancy_of_openIntervalLedgers
    v (cleanModulus v 534)
      (Nat.zero_lt_of_lt (one_lt_cleanModulus v 534 hv)) hledgers

#print axioms ceilIntervalCount_discrepancy
#print axioms openIntervalGridCount_discrepancy
#print axioms rawPairGridDiscrepancy_of_intervalLedger
#print axioms rawPairGridDiscrepancy_of_openIntervalLedger
#print axioms normalizedMass2_sub_continuumMass2_le
#print axioms normalizedMass2_sub_continuumMass2_le_of_openIntervalLedgers
#print axioms normalizedMass2_clean_534_gt_neg_target_of_openIntervalLedgers

/-! ## Axiom audit -/

#print axioms allGridJointFail_eq_jointFail_add_one
#print axioms normalizedMass2_real_eq_pairGridLedger
#print axioms pair_incidence_sum
#print axioms remove_zero_normalization_discrepancy
#print axioms pairGridDiscrepancy_of_raw
#print axioms normalizedMass2_sub_continuumMass2_le

end
end LRCB5PairGridBridge
end LonelyRunner
