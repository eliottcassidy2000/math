import TournamentH7.LRCDenseCoreRelationTrap
import TournamentH7.LRCEndgameUniformThreePhase
import TournamentH7.LRCPairTowerGapTwoProducer
import TournamentH7.LRCRealRelationLock
import TournamentH7.LRCSelectedWitnessTwoFourFourFrequency

/-!
# Exact selected-witness relation router

The zero-frequency residues of the `q333`, `q244`, and `q4488` parallel-class
circles all produce the same arithmetic object: a three-term relation on the
original selected speeds with coefficients in `{1, -1}`.  Such a relation has
below-mass at most two, so the dense-core factor-three ladder localizes its top
index to at most one step above the last dense pair.  This is sharper for
support three than the general unit-relation trap.

Only the exact zero branches are routed.  A small nonzero scalar frequency is
not rounded to zero; the explicit `q4488` example below satisfies all residue
types and the joint zero/small residual with one genuinely nonzero frequency.

Relation-lock guardrail: `real_relation_lock3` can transfer the relation to
integer witness labels only when the three danger inequalities hold at one
common real phase.  The matching-wall witnesses naturally live at different
shifted branch phases `(u+c)/q`, so the parallel-class producer alone does not
supply that hypothesis.  Likewise, one support-three relation is not an
analytic B5 relation-mass certificate.

Tournament / Zarankiewicz audit: the useful vertices here are the three
relation-support positions (or the four/eight branch residues before taking
the quotient), not runners.  The quotient to a signed relation preserves the
zero predicate and the top support, but destroys branch chronology, witness
phase alignment, and the distinction between the two q4 choices.  Thus a
`K_{2,2}` pattern on the parallel-class circle suggests the two q4/q8
rectangles, but no Zarankiewicz edge count is used as a vanishing theorem.
-/

namespace LonelyRunner
namespace LRCSelectedWitnessRelationRouter

open LRC14Grand
open LRCPairTowerGapTwoFrequency
open LRCPairTowerGapTwoProducer

/-- A three-term integer relation whose three displayed coefficients are
units.  Coordinates need not be distinct; applications record distinct
selected runners separately. -/
def SignedUnitThreeRelation (x y z : ℤ) : Prop :=
  ∃ coefficientX coefficientY coefficientZ : ℤ,
    (coefficientX = 1 ∨ coefficientX = -1) ∧
    (coefficientY = 1 ∨ coefficientY = -1) ∧
    (coefficientZ = 1 ∨ coefficientZ = -1) ∧
    coefficientX * x + coefficientY * y + coefficientZ * z = 0

theorem intSign_eq_one_or_neg_one_of_ne_zero {value : ℤ}
    (hvalue : value ≠ 0) : value.sign = 1 ∨ value.sign = -1 := by
  rcases lt_or_gt_of_ne hvalue with hnegative | hpositive
  · right
    exact Int.sign_eq_neg_one_iff_neg.mpr hnegative
  · left
    exact Int.sign_eq_one_iff_pos.mpr hpositive

/-- Orienting coefficients by the signs of three nonzero speeds transports a
signed unit relation to their absolute values without changing its support. -/
theorem signedUnitThreeRelation_abs_of_ne_zero
    {x y z : ℤ} (hx : x ≠ 0) (hy : y ≠ 0) (hz : z ≠ 0)
    (relation : SignedUnitThreeRelation x y z) :
    SignedUnitThreeRelation |x| |y| |z| := by
  obtain ⟨coefficientX, coefficientY, coefficientZ,
      hcoefficientX, hcoefficientY, hcoefficientZ, hrelation⟩ := relation
  have hsignX := intSign_eq_one_or_neg_one_of_ne_zero hx
  have hsignY := intSign_eq_one_or_neg_one_of_ne_zero hy
  have hsignZ := intSign_eq_one_or_neg_one_of_ne_zero hz
  refine ⟨x.sign * coefficientX, y.sign * coefficientY,
    z.sign * coefficientZ, ?_, ?_, ?_, ?_⟩
  · rcases hsignX with hsignX | hsignX <;>
      rcases hcoefficientX with hcoefficientX | hcoefficientX <;>
      simp [hsignX, hcoefficientX]
  · rcases hsignY with hsignY | hsignY <;>
      rcases hcoefficientY with hcoefficientY | hcoefficientY <;>
      simp [hsignY, hcoefficientY]
  · rcases hsignZ with hsignZ | hsignZ <;>
      rcases hcoefficientZ with hcoefficientZ | hcoefficientZ <;>
      simp [hsignZ, hcoefficientZ]
  · calc
      (x.sign * coefficientX) * |x| +
          (y.sign * coefficientY) * |y| +
          (z.sign * coefficientZ) * |z| =
        coefficientX * (x.sign * |x|) +
          coefficientY * (y.sign * |y|) +
          coefficientZ * (z.sign * |z|) := by ring
      _ = coefficientX * x + coefficientY * y + coefficientZ * z := by
        rw [Int.sign_mul_abs, Int.sign_mul_abs, Int.sign_mul_abs]
      _ = 0 := hrelation

/-- The normalized `q333` zero branch is a signed unit support-three relation
on the original speeds. -/
theorem signedUnitThreeRelation_of_threePhaseFrequency_eq_zero
    (speed sign : Fin 3 → ℤ)
    (hsign : ∀ index, sign index = 1 ∨ sign index = -1)
    (hresidue : ∀ index, (3 : ℤ) ∣ sign index * speed index - 1)
    (hzero : threePhaseFrequency (fun index => sign index * speed index) = 0) :
    SignedUnitThreeRelation (speed 0) (speed 1) (speed 2) := by
  refine ⟨sign 0, sign 1, sign 2, hsign 0, hsign 1, hsign 2, ?_⟩
  simpa using
    (threePhaseFrequency_eq_zero_iff_three_term_relation
      (fun index => sign index * speed index) hresidue).mp hzero

/-- The exact identity published by the normalized `q244` wall turns its zero
frequency branch into a signed unit support-three relation. -/
theorem signedUnitThreeRelation_of_twoFourFourFrequency_eq_zero
    (deltaTwo deltaFourA deltaFourB modulus : ℤ)
    (signTwo signFourA signFourB
      numeratorTwo numeratorFourA numeratorFourB : ℤ)
    (hsignTwo : signTwo = 1 ∨ signTwo = -1)
    (hsignFourA : signFourA = 1 ∨ signFourA = -1)
    (hsignFourB : signFourB = 1 ∨ signFourB = -1)
    (hidentity :
      signTwo * deltaTwo + signFourA * deltaFourA +
          signFourB * deltaFourB =
        modulus * twoFourFourPhaseFrequency
          numeratorTwo numeratorFourA numeratorFourB)
    (hzero :
      twoFourFourPhaseFrequency numeratorTwo numeratorFourA numeratorFourB = 0) :
    SignedUnitThreeRelation deltaTwo deltaFourA deltaFourB := by
  refine ⟨signTwo, signFourA, signFourB,
    hsignTwo, hsignFourA, hsignFourB, ?_⟩
  rw [hzero, mul_zero] at hidentity
  exact hidentity

/-- A signed reduced normal form gives the exact scale equation needed to
return from primitive numerators to the original selected speed. -/
theorem signedReducedNormalForm_scale
    (delta modulus denominator numerator : ℤ)
    (normalForm :
      SignedReducedNormalForm delta modulus denominator numerator) :
    denominator * (normalForm.sign * delta) = modulus * numerator := by
  rcases normalForm with
    ⟨divisor, sign, hdivisor, hdelta, hmodulus, hsign, hcoprime⟩
  subst delta
  subst modulus
  rcases hsign with rfl | rfl <;> dsimp <;> ring

/-- The zero scalar relation from either q4 side of the `q4488` rectangle is
a signed unit relation on that q4 speed and the two q8 speeds. -/
theorem signedUnitThreeRelation_of_fourEightEight_relation_eq_zero
    (deltaFour deltaEightA deltaEightB modulus
      numeratorFour numeratorEightA numeratorEightB : ℤ)
    (normalFour :
      SignedReducedNormalForm deltaFour modulus 4 numeratorFour)
    (normalEightA :
      SignedReducedNormalForm deltaEightA modulus 8 numeratorEightA)
    (normalEightB :
      SignedReducedNormalForm deltaEightB modulus 8 numeratorEightB)
    (hzero : -2 * numeratorFour + numeratorEightA + numeratorEightB = 0) :
    SignedUnitThreeRelation deltaFour deltaEightA deltaEightB := by
  have hsignFour : -normalFour.sign = 1 ∨ -normalFour.sign = -1 := by
    rcases normalFour.sign_eq with hsign | hsign <;> omega
  refine ⟨-normalFour.sign, normalEightA.sign, normalEightB.sign,
    hsignFour, normalEightA.sign_eq, normalEightB.sign_eq, ?_⟩
  have hscaleFour := signedReducedNormalForm_scale
    deltaFour modulus 4 numeratorFour normalFour
  have hscaleEightA := signedReducedNormalForm_scale
    deltaEightA modulus 8 numeratorEightA normalEightA
  have hscaleEightB := signedReducedNormalForm_scale
    deltaEightB modulus 8 numeratorEightB normalEightB
  have hmultiple :
      8 * ((-normalFour.sign) * deltaFour +
          normalEightA.sign * deltaEightA +
          normalEightB.sign * deltaEightB) =
        modulus *
          (-2 * numeratorFour + numeratorEightA + numeratorEightB) := by
    calc
      8 * ((-normalFour.sign) * deltaFour +
          normalEightA.sign * deltaEightA +
          normalEightB.sign * deltaEightB) =
          -2 * (4 * (normalFour.sign * deltaFour)) +
            8 * (normalEightA.sign * deltaEightA) +
            8 * (normalEightB.sign * deltaEightB) := by ring
      _ = -2 * (modulus * numeratorFour) +
            modulus * numeratorEightA + modulus * numeratorEightB := by
              rw [hscaleFour, hscaleEightA, hscaleEightB]
      _ = modulus *
          (-2 * numeratorFour + numeratorEightA + numeratorEightB) := by ring
  rw [hzero, mul_zero] at hmultiple
  omega

/-- A support-three unit relation cannot have its displayed top coordinate
two or more steps above the last dense pair.  This is the support-three
specialization of `no_low_mass_relation_above_pair`: two lower unit terms have
mass at most two, while one ladder step already multiplies by three. -/
theorem signedUnitThreeRelation_top_le_lastDensePair_add_one
    (weight : Fin 13 → ℤ) (hpositive : ∀ index, 0 < weight index)
    (hmonotone : Monotone weight) (lastDensePair : Fin 12)
    (hladder : ∀ index : Fin 12, lastDensePair < index →
      3 * weight index.castSucc ≤ weight index.succ)
    (first second top : Fin 13)
    (hfirstTop : first < top) (hsecondTop : second < top)
    (hrelation :
      SignedUnitThreeRelation (weight first) (weight second) (weight top)) :
    (top : ℕ) ≤ (lastDensePair : ℕ) + 1 := by
  obtain ⟨coefficientFirst, coefficientSecond, coefficientTop,
      hcoefficientFirst, hcoefficientSecond, hcoefficientTop,
      hrelation⟩ := hrelation
  by_contra htop
  have htopTwo : (lastDensePair : ℕ) + 2 ≤ (top : ℕ) := by omega
  let predecessor : Fin 13 := ⟨(top : ℕ) - 1, by omega⟩
  have hstep : 3 * weight predecessor ≤ weight top := by
    have h := ladder_top_step weight lastDensePair hladder
      (top : ℕ) htopTwo top.isLt
    have htopEq : (⟨(top : ℕ), top.isLt⟩ : Fin 13) = top := Fin.ext rfl
    rw [htopEq] at h
    exact h
  have hfirstPredecessor : first ≤ predecessor := by
    show (first : ℕ) ≤ (top : ℕ) - 1
    omega
  have hsecondPredecessor : second ≤ predecessor := by
    show (second : ℕ) ≤ (top : ℕ) - 1
    omega
  have hfirstWeight : weight first ≤ weight predecessor :=
    hmonotone hfirstPredecessor
  have hsecondWeight : weight second ≤ weight predecessor :=
    hmonotone hsecondPredecessor
  have hpredecessorPositive : 0 < weight predecessor :=
    hpositive predecessor
  have hfirstPositive : 0 < weight first := hpositive first
  have hsecondPositive : 0 < weight second := hpositive second
  rcases hcoefficientFirst with rfl | rfl <;>
    rcases hcoefficientSecond with rfl | rfl <;>
    rcases hcoefficientTop with rfl | rfl <;>
    simp only [one_mul, neg_one_mul] at hrelation <;> omega

/-- Direct selected-speed router.  After a sorting permutation, an exact
signed support-three relation on original speeds can only use a top position
at most one step above the last dense pair. -/
theorem selectedSpeedRelation_top_le_lastDensePair_add_one
    (speed : Fin 13 → ℤ) (hspeed : ∀ index, speed index ≠ 0)
    (order : Equiv.Perm (Fin 13))
    (hmonotone : Monotone (fun index => |speed (order index)|))
    (lastDensePair : Fin 12)
    (hladder : ∀ index : Fin 12, lastDensePair < index →
      3 * |speed (order index.castSucc)| ≤ |speed (order index.succ)|)
    (first second top : Fin 13)
    (hfirstTop : first < top) (hsecondTop : second < top)
    (hrelation : SignedUnitThreeRelation
      (speed (order first)) (speed (order second)) (speed (order top))) :
    (top : ℕ) ≤ (lastDensePair : ℕ) + 1 := by
  apply signedUnitThreeRelation_top_le_lastDensePair_add_one
    (fun index => |speed (order index)|)
    (fun index => abs_pos.mpr (hspeed (order index)))
    hmonotone lastDensePair hladder first second top hfirstTop hsecondTop
  exact signedUnitThreeRelation_abs_of_ne_zero
    (hspeed (order first)) (hspeed (order second)) (hspeed (order top)) hrelation

/-- At a genuinely common phase, the exact support-three speed relation locks
the three integer witness labels as well.  Shifted parallel-class wall phases
must not be supplied to this theorem as though they were common. -/
theorem SignedUnitThreeRelation.lock_common_phase
    {speedX speedY speedZ : ℤ}
    (relation : SignedUnitThreeRelation speedX speedY speedZ)
    (witnessX witnessY witnessZ : ℤ) (phase : ℝ)
    (hbadX : |(speedX : ℝ) * phase - witnessX| < 1 / 14)
    (hbadY : |(speedY : ℝ) * phase - witnessY| < 1 / 14)
    (hbadZ : |(speedZ : ℝ) * phase - witnessZ| < 1 / 14) :
    ∃ coefficientX coefficientY coefficientZ : ℤ,
      (coefficientX = 1 ∨ coefficientX = -1) ∧
      (coefficientY = 1 ∨ coefficientY = -1) ∧
      (coefficientZ = 1 ∨ coefficientZ = -1) ∧
      coefficientX * witnessX + coefficientY * witnessY +
        coefficientZ * witnessZ = 0 := by
  obtain ⟨coefficientX, coefficientY, coefficientZ,
      hcoefficientX, hcoefficientY, hcoefficientZ, hspeedRelation⟩ := relation
  refine ⟨coefficientX, coefficientY, coefficientZ,
    hcoefficientX, hcoefficientY, hcoefficientZ, ?_⟩
  have habsX : |coefficientX| = 1 := by
    rcases hcoefficientX with rfl | rfl <;> norm_num
  have habsY : |coefficientY| = 1 := by
    rcases hcoefficientY with rfl | rfl <;> norm_num
  have habsZ : |coefficientZ| = 1 := by
    rcases hcoefficientZ with rfl | rfl <;> norm_num
  exact LRCRealRelationLock.real_relation_lock3
    speedX speedY speedZ witnessX witnessY witnessZ
    coefficientX coefficientY coefficientZ phase hspeedRelation
    (by rw [habsX, habsY, habsZ]; norm_num) hbadX hbadY hbadZ

/-- Small nonzero frequencies survive the exact residue geometry.  Here both
q4 numerators and both q8 numerators have common residue one, the first q4
frequency is `-1`, the second is zero, and the full joint residual holds at
`B = 1`. -/
theorem q4488_small_nonzero_frequency_guardrail :
    (4 : ℤ) ∣ 5 - 1 ∧
    (4 : ℤ) ∣ 1 - 1 ∧
    (8 : ℤ) ∣ 1 - 1 ∧
    fourEightEightPhaseFrequency 5 1 1 = -1 ∧
    fourEightEightPhaseFrequency 1 1 1 = 0 ∧
    TwoFourEightEightFrequencyEscapeOrRelation 5 1 1 1 1 := by
  norm_num [fourEightEightPhaseFrequency,
    TwoFourEightEightFrequencyEscapeOrRelation]

#print axioms signedUnitThreeRelation_of_threePhaseFrequency_eq_zero
#print axioms signedUnitThreeRelation_of_twoFourFourFrequency_eq_zero
#print axioms signedUnitThreeRelation_of_fourEightEight_relation_eq_zero
#print axioms signedUnitThreeRelation_top_le_lastDensePair_add_one
#print axioms selectedSpeedRelation_top_le_lastDensePair_add_one
#print axioms SignedUnitThreeRelation.lock_common_phase
#print axioms q4488_small_nonzero_frequency_guardrail

end LRCSelectedWitnessRelationRouter
end LonelyRunner
