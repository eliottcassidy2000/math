import TournamentH7.LRCB5WrappedToothOverlap
import TournamentH7.LRCB5DifferenceFibers
import TournamentH7.LRCB5PairOverlapSum
import TournamentH7.LRCB5MergeLength

/-!
Exact finite reindexing of rational open-comb intersections by cyclic
difference classes.  The geometric tooth-pair observable is transported from
`Fin first × Fin second` to `ZMod (g * p * q)`; for coprime reduced speeds
`p,q`, every difference fiber has cardinality `g`.  Thus the circle's
parallel classes contribute a uniform factor `g`, closing the scaled overlap
ledger required by the Bernoulli covariance formula.

The quotient preserves the complete summed pair-overlap length but deliberately
discards the labels of individual tooth pairs.  This challenges the default
runner-vertex viewpoint: the effective vertices here are tooth pairs and cyclic
residue classes.  Tournament Analysis is not natural for this equality proof,
because no comparison switch or tie orientation is introduced; imposing one
would discard the additive fiber invariant used by the reindexing.
-/

namespace LonelyRunner
namespace LRCB5CombReindexing

open RatIntervals LRCRationalOpenComb
open LRCB5WrappedToothOverlap LRCB5DifferenceFibers
open LRCB5PairOverlapSum LRCB5MergeLength
open LRC14.ClosedBudget

noncomputable section

def rawToothStart (w k : ℕ) : ℚ :=
  ((k : ℚ) - 1 / 14) / w

def toothLengthQ (w : ℕ) : ℚ :=
  1 / (7 * w)

def circularTooth (w k : ℕ) : Region :=
  circleInterval (Int.fract (rawToothStart w k)) (toothLengthQ w)

def circularCombRegion (w : ℕ) : Region :=
  (List.range w).flatMap (circularTooth w)

theorem fract_rawToothStart_zero
    (w : ℕ) (hw : 0 < w) :
    Int.fract (rawToothStart w 0) = 1 - 1 / (14 * w : ℚ) := by
  have hwQ : (0 : ℚ) < w := by exact_mod_cast hw
  let small : ℚ := 1 / (14 * w)
  have hsmall0 : 0 < small := by dsimp [small]; positivity
  have hsmall1 : small < 1 := by
    dsimp [small]
    rw [div_lt_one (by positivity : (0 : ℚ) < 14 * w)]
    have hwOne : (1 : ℚ) ≤ w := by exact_mod_cast hw
    linarith
  have hfract : Int.fract small = small :=
    Int.fract_eq_self.mpr ⟨hsmall0.le, hsmall1⟩
  have hneg := Int.fract_neg (show Int.fract small ≠ 0 by
    rw [hfract]
    exact hsmall0.ne')
  rw [hfract] at hneg
  rw [show rawToothStart w 0 = -small by
    dsimp [rawToothStart, small]
    field_simp
    ring]
  exact hneg

theorem fract_rawToothStart_internal
    (w k : ℕ) (hw : 0 < w) (hk0 : 0 < k) (hkw : k < w) :
    Int.fract (rawToothStart w k) = rawToothStart w k := by
  apply Int.fract_eq_self.mpr
  have hwQ : (0 : ℚ) < w := by exact_mod_cast hw
  have hk0Q : (1 : ℚ) ≤ k := by exact_mod_cast hk0
  have hkwQ : (k : ℚ) < w := by exact_mod_cast hkw
  unfold rawToothStart
  constructor
  · exact div_nonneg (by linarith [show (0 : ℚ) < 1 / 14 by norm_num]) hwQ.le
  · rw [div_lt_one hwQ]
    linarith

theorem ratOpenCombInterval_zero
    (w : ℕ) (hw : 0 < w) :
    ratOpenCombInterval w 0 = (0, 1 / (14 * w : ℚ)) := by
  have hwQ : (0 : ℚ) < w := by exact_mod_cast hw
  unfold ratOpenCombInterval
  norm_num only [Nat.cast_zero]
  apply Prod.ext
  · simp only [Prod.fst]
    have hnonpos : (-(1 / 14 : ℚ) / w) ≤ 0 :=
      div_nonpos_of_nonpos_of_nonneg (by norm_num) hwQ.le
    rw [max_eq_left hnonpos]
  · simp only [Prod.snd]
    have hle : ((1 / 14 : ℚ) / w) ≤ 1 := by
      rw [div_le_one hwQ]
      have hwOne : (1 : ℚ) ≤ w := by exact_mod_cast hw
      linarith
    rw [min_eq_right hle]
    field_simp

theorem ratOpenCombInterval_boundary
    (w : ℕ) (hw : 0 < w) :
    ratOpenCombInterval w w = (1 - 1 / (14 * w : ℚ), 1) := by
  have hwQ : (0 : ℚ) < w := by exact_mod_cast hw
  unfold ratOpenCombInterval
  apply Prod.ext
  · simp only [Prod.fst]
    have hnonneg : 0 ≤ (((w : ℚ) - 1 / 14) / w) := by
      apply div_nonneg
      · have hwOne : (1 : ℚ) ≤ w := by exact_mod_cast hw
        linarith
      · exact hwQ.le
    rw [max_eq_right hnonneg]
    field_simp
  · simp only [Prod.snd]
    have hone : 1 ≤ (((w : ℚ) + 1 / 14) / w) := by
      rw [one_le_div hwQ]
      linarith
    rw [min_eq_left hone]

theorem circularTooth_zero_eq_boundary
    (w : ℕ) (hw : 0 < w) :
    circularTooth w 0 =
      [ratOpenCombInterval w w, ratOpenCombInterval w 0] := by
  have hwQ : (0 : ℚ) < w := by exact_mod_cast hw
  rw [circularTooth, fract_rawToothStart_zero w hw,
    ratOpenCombInterval_zero w hw, ratOpenCombInterval_boundary w hw]
  unfold circleInterval toothLengthQ
  have hwrap : ¬1 - 1 / (14 * w : ℚ) + 1 / (7 * w : ℚ) ≤ 1 := by
    intro h
    have hpositive : (0 : ℚ) < 1 / (14 * w) := by positivity
    field_simp at h
    linarith
  rw [if_neg hwrap]
  rw [List.cons.injEq]
  constructor
  · rfl
  · rw [List.cons.injEq]
    constructor
    · apply Prod.ext
      · rfl
      · simp only [Prod.snd]
        field_simp
        ring
    · rfl

theorem circularTooth_internal_eq_interval
    (w k : ℕ) (hw : 0 < w) (hk0 : 0 < k) (hkw : k < w) :
    circularTooth w k = [ratOpenCombInterval w k] := by
  have hwQ : (0 : ℚ) < w := by exact_mod_cast hw
  have hk0Q : (1 : ℚ) ≤ k := by exact_mod_cast hk0
  have hkwQ : (k : ℚ) < w := by exact_mod_cast hkw
  have hkwMarginQ : (k : ℚ) + 1 ≤ w := by exact_mod_cast (show k + 1 ≤ w by omega)
  rw [circularTooth, fract_rawToothStart_internal w k hw hk0 hkw]
  unfold circleInterval toothLengthQ rawToothStart ratOpenCombInterval
  have hnowrap :
      ((k : ℚ) - 1 / 14) / w + 1 / (7 * w : ℚ) ≤ 1 := by
    field_simp
    linarith
  rw [if_pos hnowrap]
  congr 1
  apply Prod.ext
  · simp only [Prod.fst]
    rw [max_eq_right]
    exact div_nonneg (by linarith) hwQ.le
  · simp only [Prod.snd]
    rw [min_eq_right]
    · field_simp
      ring
    · rw [div_le_one hwQ]
      linarith

theorem ratOpenCombTail_succ_eq_append
    (w start count : ℕ) :
    ratOpenCombTail w start (count + 1) =
      ratOpenCombTail w start count ++
        [ratOpenCombInterval w (start + count)] := by
  induction count generalizing start with
  | zero => simp [ratOpenCombTail]
  | succ count inductionHypothesis =>
      have hcons := congrArg
        (fun tail => ratOpenCombInterval w start :: tail)
        (inductionHypothesis (start + 1))
      simpa [ratOpenCombTail, Nat.add_comm, Nat.add_left_comm,
        Nat.add_assoc] using hcons

theorem ratOpenCombRegion_eq_tail_append_boundary (w : ℕ) :
    ratOpenCombRegion w =
      ratOpenCombTail w 0 w ++ [ratOpenCombInterval w w] := by
  unfold ratOpenCombRegion
  simpa using ratOpenCombTail_succ_eq_append w 0 w

theorem circularInternalFlatMap_eq_tail
    (w start count : ℕ) (hw : 0 < w) (hstart : 0 < start)
    (hbound : start + count ≤ w) :
    (List.range' start count).flatMap (circularTooth w) =
      ratOpenCombTail w start count := by
  induction count generalizing start with
  | zero => simp [ratOpenCombTail]
  | succ count inductionHypothesis =>
      rw [List.range'_succ, List.flatMap_cons,
        circularTooth_internal_eq_interval w start hw hstart (by omega),
        ratOpenCombTail,
        inductionHypothesis (start + 1) (by omega) (by omega)]
      rfl

theorem circularCombRegion_eq_boundary_cons_tail
    (w : ℕ) (hw : 0 < w) :
    circularCombRegion w =
      ratOpenCombInterval w w :: ratOpenCombTail w 0 w := by
  obtain ⟨count, rfl⟩ := Nat.exists_eq_succ_of_ne_zero hw.ne'
  unfold circularCombRegion
  rw [List.range_eq_range', List.range'_succ, List.flatMap_cons,
    circularTooth_zero_eq_boundary (count + 1) (by omega),
    circularInternalFlatMap_eq_tail (count + 1) 1 count
      (by omega) (by omega) (by omega),
    ratOpenCombTail]
  rfl

theorem circularCombRegion_perm_ratOpenCombRegion
    (w : ℕ) (hw : 0 < w) :
    (circularCombRegion w).Perm (ratOpenCombRegion w) := by
  rw [circularCombRegion_eq_boundary_cons_tail w hw,
    ratOpenCombRegion_eq_tail_append_boundary]
  simpa using
    (List.perm_append_comm (l₁ := [ratOpenCombInterval w w])
      (l₂ := ratOpenCombTail w 0 w))

theorem length_eq_of_perm {left right : Region} (hperm : left.Perm right) :
    RatIntervals.length left = RatIntervals.length right := by
  unfold RatIntervals.length
  exact (hperm.map fun interval => max 0 (interval.2 - interval.1)).sum_eq

theorem length_inter_eq_of_perm_left
    {first first' second : Region} (hperm : first.Perm first') :
    RatIntervals.length (inter first second) =
      RatIntervals.length (inter first' second) := by
  unfold inter
  apply length_eq_of_perm
  exact hperm.flatMap fun _ _ => List.Perm.refl _

theorem length_inter_eq_of_perm_right
    {first second second' : Region} (hperm : second.Perm second') :
    RatIntervals.length (inter first second) =
      RatIntervals.length (inter first second') := by
  rw [RatIntervals.length_inter_comm first second,
    RatIntervals.length_inter_comm first second']
  exact length_inter_eq_of_perm_left hperm

theorem length_inter_ratOpenCombRegion_eq_circularCombRegion
    (first second : ℕ) (hfirst : 0 < first) (hsecond : 0 < second) :
    RatIntervals.length
        (inter (ratOpenCombRegion first) (ratOpenCombRegion second)) =
      RatIntervals.length
        (inter (circularCombRegion first) (circularCombRegion second)) := by
  calc
    RatIntervals.length
        (inter (ratOpenCombRegion first) (ratOpenCombRegion second)) =
        RatIntervals.length
          (inter (circularCombRegion first) (ratOpenCombRegion second)) :=
      length_inter_eq_of_perm_left
        (circularCombRegion_perm_ratOpenCombRegion first hfirst).symm
    _ = RatIntervals.length
          (inter (circularCombRegion first) (circularCombRegion second)) :=
      length_inter_eq_of_perm_right
        (circularCombRegion_perm_ratOpenCombRegion second hsecond).symm

theorem fract_sub_eq_fract_fract_sub_fractQ (first second : ℚ) :
    Int.fract (first - second) =
      Int.fract (Int.fract first - Int.fract second) := by
  rw [Int.fract_eq_fract]
  refine ⟨⌊first⌋ - ⌊second⌋, ?_⟩
  simp only [Int.fract]
  rw [Int.cast_sub]
  ring

theorem toothLengthQ_nonneg (w : ℕ) :
    0 ≤ toothLengthQ w := by
  unfold toothLengthQ
  positivity

theorem toothLengthQ_le_one (w : ℕ) (hw : 0 < w) :
    toothLengthQ w ≤ 1 := by
  have hdenNat : 1 ≤ 7 * w := by omega
  have hden : (1 : ℚ) ≤ 7 * w := by exact_mod_cast hdenNat
  unfold toothLengthQ
  exact (div_le_one (by positivity)).2 hden

theorem length_inter_circularTooth_eq_pairOverlapQ
    (first second firstIndex secondIndex : ℕ)
    (hfirst : 0 < first) (hsecond : 0 < second) :
    RatIntervals.length
        (inter (circularTooth first firstIndex)
          (circularTooth second secondIndex)) =
      pairOverlapQ (toothLengthQ first) (toothLengthQ second)
        (Int.fract
          (rawToothStart first firstIndex -
            rawToothStart second secondIndex)) := by
  unfold circularTooth
  rw [length_inter_circleInterval_eq_pairOverlapQ
      (Int.fract (rawToothStart first firstIndex))
      (Int.fract (rawToothStart second secondIndex))
      (toothLengthQ first) (toothLengthQ second)
      (Int.fract_nonneg _) (Int.fract_lt_one _)
      (Int.fract_nonneg _) (Int.fract_lt_one _)
      (toothLengthQ_nonneg first) (toothLengthQ_le_one first hfirst)
      (toothLengthQ_nonneg second) (toothLengthQ_le_one second hsecond)]
  rw [← fract_sub_eq_fract_fract_sub_fractQ]

theorem length_inter_flatMap_right
    {index : Type} (left : Region) (indices : List index)
    (region : index → Region) :
    RatIntervals.length (inter left (indices.flatMap region)) =
      (indices.map fun value =>
        RatIntervals.length (inter left (region value))).sum := by
  induction indices with
  | nil => simp [RatIntervals.length_inter_nil]
  | cons head tail inductionHypothesis =>
      rw [List.flatMap_cons, RatIntervals.length_inter_append_right,
        inductionHypothesis, List.map_cons, List.sum_cons]

theorem length_inter_circularCombRegion_eq_doubleSum
    (first second : ℕ) (hfirst : 0 < first) (hsecond : 0 < second) :
    RatIntervals.length
        (inter (circularCombRegion first) (circularCombRegion second)) =
      ((List.range second).map fun secondIndex =>
        ((List.range first).map fun firstIndex =>
          pairOverlapQ (toothLengthQ first) (toothLengthQ second)
            (Int.fract
              (rawToothStart first firstIndex -
                rawToothStart second secondIndex))).sum).sum := by
  unfold circularCombRegion
  rw [length_inter_flatMap_right]
  apply congrArg List.sum
  apply List.map_congr_left
  intro secondIndex hsecondIndex
  rw [RatIntervals.length_inter_comm, length_inter_flatMap_right]
  apply congrArg List.sum
  apply List.map_congr_left
  intro firstIndex hfirstIndex
  rw [RatIntervals.length_inter_comm,
    length_inter_circularTooth_eq_pairOverlapQ
      first second firstIndex secondIndex hfirst hsecond]

theorem length_inter_ratOpenCombRegion_eq_doubleSum
    (first second : ℕ) (hfirst : 0 < first) (hsecond : 0 < second) :
    RatIntervals.length
        (inter (ratOpenCombRegion first) (ratOpenCombRegion second)) =
      ((List.range second).map fun secondIndex =>
        ((List.range first).map fun firstIndex =>
          pairOverlapQ (toothLengthQ first) (toothLengthQ second)
            (Int.fract
              (rawToothStart first firstIndex -
                rawToothStart second secondIndex))).sum).sum := by
  rw [length_inter_ratOpenCombRegion_eq_circularCombRegion
    first second hfirst hsecond]
  exact length_inter_circularCombRegion_eq_doubleSum
    first second hfirst hsecond

theorem cast_pairOverlapQ (length₁ length₂ shift : ℚ) :
    ((pairOverlapQ length₁ length₂ shift : ℚ) : ℝ) =
      pairOverlap (length₁ : ℝ) (length₂ : ℝ) (shift : ℝ) := by
  unfold pairOverlapQ pairOverlap
  push_cast
  rfl

theorem cast_list_sum (values : List ℚ) :
    ((values.sum : ℚ) : ℝ) =
      (values.map fun value => (value : ℝ)).sum := by
  induction values with
  | nil => simp
  | cons head tail inductionHypothesis =>
      simp [inductionHypothesis]

theorem list_sum_range_eq_fin_sum
    (count : ℕ) (weight : ℕ → ℝ) :
    ((List.range count).map weight).sum =
      ∑ index : Fin count, weight index := by
  have hlist :
      ((List.range count).map weight).sum =
        ∑ index ∈ Finset.range count, weight index := by
    induction count with
    | zero => simp
    | succ count inductionHypothesis =>
        rw [List.range_succ, List.map_append, List.sum_append,
          Finset.sum_range_succ, inductionHypothesis]
        simp
  rw [hlist]
  exact (Fin.sum_univ_eq_sum_range weight count).symm

def realToothOverlap
    (first second firstIndex secondIndex : ℕ) : ℝ :=
  pairOverlap
    (1 / (7 * (first : ℝ)))
    (1 / (7 * (second : ℝ)))
    (Int.fract
      ((((firstIndex : ℕ) : ℝ) - 1 / 14) / first -
        (((secondIndex : ℕ) : ℝ) - 1 / 14) / second))

theorem cast_pairOverlapQ_rawToothStart
    (first second firstIndex secondIndex : ℕ) :
    ((pairOverlapQ (toothLengthQ first) (toothLengthQ second)
        (Int.fract
          (rawToothStart first firstIndex -
            rawToothStart second secondIndex)) : ℚ) : ℝ) =
      realToothOverlap first second firstIndex secondIndex := by
  rw [cast_pairOverlapQ, Rat.cast_fract]
  unfold toothLengthQ rawToothStart realToothOverlap
  push_cast
  rfl

def rationalToothOverlap
    (first second firstIndex secondIndex : ℕ) : ℚ :=
  pairOverlapQ (toothLengthQ first) (toothLengthQ second)
    (Int.fract
      (rawToothStart first firstIndex - rawToothStart second secondIndex))

def rationalOverlapRow (first second secondIndex : ℕ) : ℚ :=
  ((List.range first).map fun firstIndex =>
    rationalToothOverlap first second firstIndex secondIndex).sum

def rationalOverlapLedger (first second : ℕ) : ℚ :=
  ((List.range second).map fun secondIndex =>
    rationalOverlapRow first second secondIndex).sum

def realOverlapRow (first second secondIndex : ℕ) : ℝ :=
  ((List.range first).map fun firstIndex =>
    realToothOverlap first second firstIndex secondIndex).sum

theorem cast_rationalToothOverlap
    (first second firstIndex secondIndex : ℕ) :
    ((rationalToothOverlap first second firstIndex secondIndex : ℚ) : ℝ) =
      realToothOverlap first second firstIndex secondIndex := by
  exact cast_pairOverlapQ_rawToothStart
    first second firstIndex secondIndex

theorem cast_rationalOverlapRow
    (first second secondIndex : ℕ) :
    ((rationalOverlapRow first second secondIndex : ℚ) : ℝ) =
      realOverlapRow first second secondIndex := by
  unfold rationalOverlapRow realOverlapRow
  rw [cast_list_sum]
  simp only [bind_pure_comp, List.map_eq_map, List.map_map,
    Function.comp_apply]
  apply congrArg List.sum
  apply List.map_congr_left
  intro firstIndex hfirstIndex
  exact cast_rationalToothOverlap
    first second firstIndex secondIndex

theorem cast_rationalOverlapLedger_eq_finSum
    (first second : ℕ) :
    ((rationalOverlapLedger first second : ℚ) : ℝ) =
      ∑ secondIndex : Fin second, ∑ firstIndex : Fin first,
        realToothOverlap first second firstIndex secondIndex := by
  unfold rationalOverlapLedger
  rw [cast_list_sum]
  simp only [bind_pure_comp, List.map_eq_map, List.map_map,
    Function.comp_apply]
  calc
    ((List.range second).map fun secondIndex =>
        ((rationalOverlapRow first second secondIndex : ℚ) : ℝ)).sum =
        ((List.range second).map fun secondIndex =>
          realOverlapRow first second secondIndex).sum := by
      apply congrArg List.sum
      apply List.map_congr_left
      intro secondIndex hsecondIndex
      exact cast_rationalOverlapRow first second secondIndex
    _ = ∑ secondIndex : Fin second,
          realOverlapRow first second secondIndex :=
      list_sum_range_eq_fin_sum second _
    _ = ∑ secondIndex : Fin second, ∑ firstIndex : Fin first,
          realToothOverlap first second firstIndex secondIndex := by
      apply Fintype.sum_congr
      intro secondIndex
      unfold realOverlapRow
      exact list_sum_range_eq_fin_sum first _

theorem cast_length_inter_ratOpenCombRegion_eq_finProduct
    (first second : ℕ) (hfirst : 0 < first) (hsecond : 0 < second) :
    ((RatIntervals.length
        (inter (ratOpenCombRegion first) (ratOpenCombRegion second)) : ℚ) : ℝ) =
      ∑ pair : Fin first × Fin second,
        realToothOverlap first second pair.1 pair.2 := by
  have hlength := length_inter_ratOpenCombRegion_eq_doubleSum
    first second hfirst hsecond
  have hledger :
      RatIntervals.length
          (inter (ratOpenCombRegion first) (ratOpenCombRegion second)) =
        rationalOverlapLedger first second := by
    simpa [rationalOverlapLedger, rationalOverlapRow,
      rationalToothOverlap] using hlength
  have hcast := congrArg (fun value : ℚ => (value : ℝ)) hledger
  calc
    ((RatIntervals.length
        (inter (ratOpenCombRegion first) (ratOpenCombRegion second)) : ℚ) : ℝ) =
        ((rationalOverlapLedger first second : ℚ) : ℝ) := hcast
    _ =
        ∑ secondIndex : Fin second, ∑ firstIndex : Fin first,
          realToothOverlap first second firstIndex secondIndex :=
      cast_rationalOverlapLedger_eq_finSum first second
    _ = ∑ firstIndex : Fin first, ∑ secondIndex : Fin second,
          realToothOverlap first second firstIndex secondIndex := by
      rw [Finset.sum_comm]
    _ = ∑ pair : Fin first × Fin second,
          realToothOverlap first second pair.1 pair.2 := by
      rw [Fintype.sum_prod_type]

def zmodToothOverlap
    (first second : ℕ) (pair : ZMod first × ZMod second) : ℝ :=
  realToothOverlap first second pair.1.val pair.2.val

theorem val_finEquiv
    (modulus : ℕ) [NeZero modulus] (index : Fin modulus) :
    ((ZMod.finEquiv modulus) index).val = index.val := by
  cases modulus with
  | zero => exact (NeZero.ne 0 rfl).elim
  | succ modulus => rfl

theorem finProduct_sum_eq_zmodProduct
    (first second : ℕ) [NeZero first] [NeZero second] :
    (∑ pair : Fin first × Fin second,
        realToothOverlap first second pair.1 pair.2) =
      ∑ pair : ZMod first × ZMod second,
        zmodToothOverlap first second pair := by
  let equivalence : Fin first × Fin second ≃
      ZMod first × ZMod second :=
    Equiv.prodCongr (ZMod.finEquiv first).toEquiv
      (ZMod.finEquiv second).toEquiv
  apply Fintype.sum_equiv equivalence
    (fun pair : Fin first × Fin second =>
      realToothOverlap first second pair.1 pair.2)
    (zmodToothOverlap first second)
  intro pair
  unfold zmodToothOverlap
  rw [show (equivalence pair).1.val = pair.1.val by
      exact val_finEquiv first pair.1,
    show (equivalence pair).2.val = pair.2.val by
      exact val_finEquiv second pair.2]

def differenceOverlapWeight
    (g p q : ℕ) (residue : ZMod (g * p * q)) : ℝ :=
  pairOverlap
    (1 / (7 * (g * p : ℕ) : ℝ))
    (1 / (7 * (g * q : ℕ) : ℝ))
    (Int.fract
      (1 / (14 * (g * q : ℕ) : ℝ) -
        1 / (14 * (g * p : ℕ) : ℝ) +
        residue.val / (g * p * q : ℕ)))

theorem fract_rawDifference_eq_differenceHom
    (g p q : ℕ) (hg : 0 < g) (hp : 0 < p) (hq : 0 < q)
    (pair : ZMod (g * p) × ZMod (g * q)) :
    Int.fract
        ((((pair.1.val : ℕ) : ℝ) - 1 / 14) / (g * p : ℕ) -
          (((pair.2.val : ℕ) : ℝ) - 1 / 14) / (g * q : ℕ)) =
      Int.fract
        (1 / (14 * (g * q : ℕ) : ℝ) -
          1 / (14 * (g * p : ℕ) : ℝ) +
          (differenceHom g p q pair).val / (g * p * q : ℕ)) := by
  letI : NeZero g := ⟨hg.ne'⟩
  letI : NeZero p := ⟨hp.ne'⟩
  letI : NeZero q := ⟨hq.ne'⟩
  letI : NeZero (g * p) := ⟨(Nat.mul_pos hg hp).ne'⟩
  letI : NeZero (g * q) := ⟨(Nat.mul_pos hg hq).ne'⟩
  letI : NeZero (g * p * q) :=
    ⟨(Nat.mul_pos (Nat.mul_pos hg hp) hq).ne'⟩
  let difference : ℤ :=
    (q : ℤ) * pair.1.val - (p : ℤ) * pair.2.val
  let residue : ℤ := (differenceHom g p q pair).val
  have hpair : pair =
      (((pair.1.val : ℕ) : ZMod (g * p)),
        ((pair.2.val : ℕ) : ZMod (g * q))) := by
    apply Prod.ext
    · exact (ZMod.natCast_zmod_val pair.1).symm
    · exact (ZMod.natCast_zmod_val pair.2).symm
  have hdifference :
      differenceHom g p q pair =
        (difference : ZMod (g * p * q)) := by
    rw [hpair]
    simpa [difference] using
      differenceHom_intCast g p q
        (pair.1.val : ℤ) (pair.2.val : ℤ)
  have hresidue :
      differenceHom g p q pair =
        (residue : ZMod (g * p * q)) := by
    simpa [residue] using
      (ZMod.natCast_zmod_val (differenceHom g p q pair)).symm
  have hequal :
      (difference : ZMod (g * p * q)) =
        (residue : ZMod (g * p * q)) := hdifference.symm.trans hresidue
  have hdivides :
      (g * p * q : ℤ) ∣ residue - difference :=
    (ZMod.intCast_eq_intCast_iff_dvd_sub
      difference residue (g * p * q)).mp hequal
  obtain ⟨multiple, hmultiple⟩ := hdivides
  have hmultipleR :
      (residue : ℝ) - (difference : ℝ) =
        (g * p * q : ℕ) * (multiple : ℝ) := by
    exact_mod_cast hmultiple
  have hperiodR : (0 : ℝ) < (g * p * q : ℕ) := by
    exact_mod_cast Nat.mul_pos (Nat.mul_pos hg hp) hq
  apply Int.fract_eq_fract.mpr
  refine ⟨-multiple, ?_⟩
  have halgebra :
      ((((pair.1.val : ℕ) : ℝ) - 1 / 14) / (g * p : ℕ) -
          (((pair.2.val : ℕ) : ℝ) - 1 / 14) / (g * q : ℕ)) -
        (1 / (14 * (g * q : ℕ) : ℝ) -
          1 / (14 * (g * p : ℕ) : ℝ) +
          (differenceHom g p q pair).val / (g * p * q : ℕ)) =
        (((difference : ℝ) - residue) / (g * p * q : ℕ)) := by
    dsimp [difference, residue]
    push_cast
    field_simp
    ring
  rw [halgebra]
  rw [Int.cast_neg]
  apply (div_eq_iff hperiodR.ne').2
  push_cast at hmultipleR ⊢
  linarith

theorem zmodToothOverlap_eq_differenceWeight
    (g p q : ℕ) (hg : 0 < g) (hp : 0 < p) (hq : 0 < q)
    (pair : ZMod (g * p) × ZMod (g * q)) :
    zmodToothOverlap (g * p) (g * q) pair =
      differenceOverlapWeight g p q (differenceHom g p q pair) := by
  unfold zmodToothOverlap realToothOverlap differenceOverlapWeight
  rw [fract_rawDifference_eq_differenceHom g p q hg hp hq pair]

theorem sum_differenceOverlapWeight_eq_gridPairOverlap
    (g p q : ℕ) [NeZero (g * p * q)]
    (hg : 0 < g) (hp : 0 < p) (hq : 0 < q) :
    (∑ residue : ZMod (g * p * q),
        differenceOverlapWeight g p q residue) =
      gridPairOverlap (g * p * q)
        (1 / (7 * (g * p : ℕ) : ℝ))
        (1 / (7 * (g * q : ℕ) : ℝ))
        (1 / (14 * (g * q : ℕ) : ℝ) -
          1 / (14 * (g * p : ℕ) : ℝ)) := by
  let equivalence : Fin (g * p * q) ≃ ZMod (g * p * q) :=
    (ZMod.finEquiv (g * p * q)).toEquiv
  calc
    (∑ residue : ZMod (g * p * q),
        differenceOverlapWeight g p q residue) =
        ∑ index : Fin (g * p * q),
          differenceOverlapWeight g p q (equivalence index) := by
      exact (Equiv.sum_comp equivalence
        (differenceOverlapWeight g p q)).symm
    _ = ∑ index : Fin (g * p * q),
          pairOverlap
            (1 / (7 * (g * p : ℕ) : ℝ))
            (1 / (7 * (g * q : ℕ) : ℝ))
            (Int.fract
              (1 / (14 * (g * q : ℕ) : ℝ) -
                1 / (14 * (g * p : ℕ) : ℝ) +
                (index : ℕ) / (g * p * q : ℕ))) := by
      apply Fintype.sum_congr
      intro index
      unfold differenceOverlapWeight
      rw [show (equivalence index).val = index.val by
        exact val_finEquiv (g * p * q) index]
    _ = gridPairOverlap (g * p * q)
          (1 / (7 * (g * p : ℕ) : ℝ))
          (1 / (7 * (g * q : ℕ) : ℝ))
          (1 / (14 * (g * q : ℕ) : ℝ) -
            1 / (14 * (g * p : ℕ) : ℝ)) := by
      unfold gridPairOverlap
      exact Fin.sum_univ_eq_sum_range
        (fun index : ℕ =>
          pairOverlap
            (1 / (7 * (g * p : ℕ) : ℝ))
            (1 / (7 * (g * q : ℕ) : ℝ))
            (Int.fract
              (1 / (14 * (g * q : ℕ) : ℝ) -
                1 / (14 * (g * p : ℕ) : ℝ) +
                (index : ℝ) / (g * p * q : ℕ))))
        (g * p * q)

theorem cast_length_inter_gp_gq_eq_scaledPairOverlapLedger
    (g p q : ℕ) (hg : 0 < g) (hp : 0 < p) (hq : 0 < q)
    (hcoprime : Nat.Coprime p q) :
    ((RatIntervals.length
        (inter (ratOpenCombRegion (g * p))
          (ratOpenCombRegion (g * q))) : ℚ) : ℝ) =
      scaledPairOverlapLedger g p q := by
  letI : NeZero g := ⟨hg.ne'⟩
  letI : NeZero p := ⟨hp.ne'⟩
  letI : NeZero q := ⟨hq.ne'⟩
  letI : NeZero (g * p) := ⟨(Nat.mul_pos hg hp).ne'⟩
  letI : NeZero (g * q) := ⟨(Nat.mul_pos hg hq).ne'⟩
  calc
    ((RatIntervals.length
        (inter (ratOpenCombRegion (g * p))
          (ratOpenCombRegion (g * q))) : ℚ) : ℝ) =
        ∑ pair : Fin (g * p) × Fin (g * q),
          realToothOverlap (g * p) (g * q) pair.1 pair.2 :=
      cast_length_inter_ratOpenCombRegion_eq_finProduct
        (g * p) (g * q) (Nat.mul_pos hg hp) (Nat.mul_pos hg hq)
    _ = ∑ pair : ZMod (g * p) × ZMod (g * q),
          zmodToothOverlap (g * p) (g * q) pair :=
      finProduct_sum_eq_zmodProduct (g * p) (g * q)
    _ = ∑ pair : ZMod (g * p) × ZMod (g * q),
          differenceOverlapWeight g p q (differenceHom g p q pair) := by
      apply Fintype.sum_congr
      intro pair
      exact zmodToothOverlap_eq_differenceWeight g p q hg hp hq pair
    _ = (g : ℝ) * ∑ residue : ZMod (g * p * q),
          differenceOverlapWeight g p q residue :=
      sum_weight_differenceHom g p q hg hp hq hcoprime
        (differenceOverlapWeight g p q)
    _ = scaledPairOverlapLedger g p q := by
      rw [sum_differenceOverlapWeight_eq_gridPairOverlap g p q hg hp hq]
      rfl

theorem cast_length_ratOpenPairRegion_gp_gq_eq_scaledPairOverlapLedger
    (g p q : ℕ) (hg : 0 < g) (hp : 0 < p) (hq : 0 < q)
    (hcoprime : Nat.Coprime p q) :
    ((RatIntervals.length
        (ratOpenPairRegion (g * p) (g * q)) : ℚ) : ℝ) =
      scaledPairOverlapLedger g p q := by
  unfold ratOpenPairRegion
  rw [length_strictInterMerge_eq_inter
    (ratOpenCombRegion (g * p)) (ratOpenCombRegion (g * q))
    (norm_ratOpenCombRegion (g * p) (Nat.mul_pos hg hp))
    (norm_ratOpenCombRegion (g * q) (Nat.mul_pos hg hq))]
  exact cast_length_inter_gp_gq_eq_scaledPairOverlapLedger
    g p q hg hp hq hcoprime

theorem cast_length_ratOpenPairRegion_eq_scaledPairOverlapLedger
    (first second : ℕ) (hfirst : 0 < first) (hsecond : 0 < second) :
    ((RatIntervals.length
        (ratOpenPairRegion first second) : ℚ) : ℝ) =
      scaledPairOverlapLedger (Nat.gcd first second)
        (first / Nat.gcd first second)
        (second / Nat.gcd first second) := by
  let g := Nat.gcd first second
  let p := first / g
  let q := second / g
  have hg : 0 < g := by
    dsimp [g]
    exact Nat.gcd_pos_of_pos_left second hfirst
  have hp : 0 < p := by
    dsimp [p]
    exact Nat.div_pos (Nat.gcd_le_left second hfirst) hg
  have hq : 0 < q := by
    dsimp [q]
    exact Nat.div_pos (Nat.gcd_le_right first hsecond) hg
  have hcoprime : Nat.Coprime p q := by
    dsimp [p, q, g]
    exact Nat.coprime_div_gcd_div_gcd hg
  have hfirstFactor : first = g * p := by
    dsimp [g, p]
    exact (Nat.mul_div_cancel' (Nat.gcd_dvd_left first second)).symm
  have hsecondFactor : second = g * q := by
    dsimp [g, q]
    exact (Nat.mul_div_cancel' (Nat.gcd_dvd_right first second)).symm
  change ((RatIntervals.length
      (ratOpenPairRegion first second) : ℚ) : ℝ) =
    scaledPairOverlapLedger g p q
  calc
    ((RatIntervals.length
        (ratOpenPairRegion first second) : ℚ) : ℝ) =
        ((RatIntervals.length
          (ratOpenPairRegion (g * p) (g * q)) : ℚ) : ℝ) := by
      rw [← hfirstFactor, ← hsecondFactor]
    _ = scaledPairOverlapLedger g p q :=
      cast_length_ratOpenPairRegion_gp_gq_eq_scaledPairOverlapLedger
        g p q hg hp hq hcoprime

#print axioms circularCombRegion_perm_ratOpenCombRegion
#print axioms length_inter_ratOpenCombRegion_eq_doubleSum
#print axioms fract_rawDifference_eq_differenceHom
#print axioms cast_length_inter_gp_gq_eq_scaledPairOverlapLedger
#print axioms cast_length_ratOpenPairRegion_eq_scaledPairOverlapLedger

end
end LRCB5CombReindexing
end LonelyRunner


