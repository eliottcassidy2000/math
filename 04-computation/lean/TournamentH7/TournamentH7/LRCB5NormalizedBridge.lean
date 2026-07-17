/-
  TournamentH7.LRCB5NormalizedBridge

  THM-940 already expresses the concrete integer B5 as an alternating sum of
  subset deviations.  Dividing its support-k aggregates by q-1 and applying a
  triangular Möbius transform produces the exact THM-935 coefficients.  The
  singleton aggregate need not vanish: the proved unit-bijection formula makes
  it nonpositive, so its negative sign can only increase B5.

  Thus the concrete consumer needs only two quantitative inputs:

    normalizedMass2 >= -13/30,
    (24/49) normalizedMass3 + (2/7) normalizedMass4 + normalizedMass5
      < 7712/84035.

  Tournament-analysis audit: the natural vertices here are subset supports,
  not runners or arcs.  Pair data forms a symmetric support-incidence relation;
  orienting ties by an index gauge would add no Fourier-sign information.  The
  Möbius quotient preserves the complete alternating B5 ledger but destroys
  the individual phase locations inside each deviation.  The challenged
  assumption is that singleton deviations must be zero: their favorable
  nonpositive sign is sufficient and strictly more general.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCB5RelationBudget
import TournamentH7.LRCB5SubsetExpansion
import TournamentH7.LRCB5CleanModulus
import TournamentH7.LRCDeviationSingles

namespace LonelyRunner
namespace LRCB5NormalizedBridge

open Finset
open LRCB5RelationBudget
open LRCB5CleanModulus

noncomputable section

/-- THM-940's total deviation on the support-`k` stratum. -/
def aggregateDeviation (v : Fin 13 → ℤ) (q k : ℕ) : ℚ :=
  ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard k,
    LRC14Concrete.deviation v q T

/-- Aggregate deviation per nonzero multiplier. -/
def normalizedAggregateDeviation (v : Fin 13 → ℤ) (q k : ℕ) : ℚ :=
  aggregateDeviation v q k / ((q : ℚ) - 1)

/-- Exact-support-two normalized mass. -/
def normalizedMass2 (v : Fin 13 → ℤ) (q : ℕ) : ℚ :=
  normalizedAggregateDeviation v q 2

/-- Exact-support-three normalized mass after binomial Möbius inversion. -/
def normalizedMass3 (v : Fin 13 → ℤ) (q : ℕ) : ℚ :=
  normalizedAggregateDeviation v q 3 -
    (11 / 7) * normalizedAggregateDeviation v q 2

/-- Exact-support-four normalized mass after binomial Möbius inversion. -/
def normalizedMass4 (v : Fin 13 → ℤ) (q : ℕ) : ℚ :=
  normalizedAggregateDeviation v q 4 -
    (10 / 7) * normalizedAggregateDeviation v q 3 +
    (55 / 49) * normalizedAggregateDeviation v q 2

/-- Exact-support-five normalized mass after binomial Möbius inversion. -/
def normalizedMass5 (v : Fin 13 → ℤ) (q : ℕ) : ℚ :=
  normalizedAggregateDeviation v q 5 -
    (9 / 7) * normalizedAggregateDeviation v q 4 +
    (45 / 49) * normalizedAggregateDeviation v q 3 -
    (165 / 343) * normalizedAggregateDeviation v q 2

/-- Rational version of the signed relation model. -/
def relationModelQ (mass2 mass3 mass4 mass5 : ℚ) : ℚ :=
  2052 / 16807 + (24 / 343) * mass2 - (24 / 49) * mass3 -
    (2 / 7) * mass4 - mass5

/-- The signed higher-support combination which can lower B5. -/
def normalizedHigherMobiusSigned (v : Fin 13 → ℤ) (q : ℕ) : ℚ :=
  (24 / 49) * normalizedMass3 v q +
    (2 / 7) * normalizedMass4 v q + normalizedMass5 v q

/-- Normalized second factorial moment of the band-depth histogram. -/
def normalizedPairDepthMoment (v : Fin 13 → ℤ) (q : ℕ) : ℚ :=
  (LRC14Concrete.momentS v q 2 : ℚ) / ((q : ℚ) - 1)

theorem aggregateDeviation_zero (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 ≤ q) :
    aggregateDeviation v q 0 = 0 := by
  simp [aggregateDeviation, LRC14Concrete.deviation_empty v q hq]

/-- Aggregate deviations are factorial moments minus their equilibrium
values. -/
theorem aggregateDeviation_eq_moment_sub_equilibrium
    (v : Fin 13 → ℤ) (q k : ℕ) :
    aggregateDeviation v q k =
      (LRC14Concrete.momentS v q k : ℚ) -
        (Nat.choose 13 k : ℚ) * (((q : ℚ) - 1) / 7 ^ k) := by
  rw [LRC14Concrete.momentS_eq_sum_jointFail]
  simp only [aggregateDeviation, LRC14Concrete.deviation,
    Finset.sum_sub_distrib]
  push_cast
  congr 1
  calc
    ∑ support ∈ (Finset.univ : Finset (Fin 13)).powersetCard k,
        ((q : ℚ) - 1) / 7 ^ support.card =
      ∑ _support ∈ (Finset.univ : Finset (Fin 13)).powersetCard k,
        ((q : ℚ) - 1) / 7 ^ k := by
          apply Finset.sum_congr rfl
          intro support hsupport
          rw [(Finset.mem_powersetCard.mp hsupport).2]
    _ = (Nat.choose 13 k : ℚ) * (((q : ℚ) - 1) / 7 ^ k) := by
      rw [Finset.sum_const, Finset.card_powersetCard, nsmul_eq_mul]
      simp

/-- Exact pair-mass reading in the depth histogram. -/
theorem normalizedMass2_eq_pairDepthMoment
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q) :
    normalizedMass2 v q = normalizedPairDepthMoment v q - 78 / 49 := by
  rw [normalizedMass2, normalizedAggregateDeviation,
    aggregateDeviation_eq_moment_sub_equilibrium]
  have hden : ((q : ℚ) - 1) ≠ 0 := by
    have : (1 : ℚ) < q := by exact_mod_cast hq
    linarith
  norm_num [normalizedPairDepthMoment, Nat.choose]
  field_simp [hden]

/-- The pair-tail target is exactly the normalized second factorial-moment
floor `1703/1470`. -/
theorem normalizedMass2_lower_iff_pairDepthMoment
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q) :
    -(13 / 30 : ℚ) ≤ normalizedMass2 v q ↔
      1703 / 1470 ≤ normalizedPairDepthMoment v q := by
  rw [normalizedMass2_eq_pairDepthMoment v q hq]
  constructor <;> intro h <;> linarith

theorem pair_floor_gap_to_horizon_thirty :
    (1703 / 1470 : ℚ) - 6 / 7 = 443 / 1470 := by
  norm_num

/-- The pair-moment excess beyond the convex baseline.  Depth zero contributes
one unit; depth `d ≥ 1` contributes `choose (d-1) 2`. -/
def pairExcessTerm (depth : ℕ) : ℕ :=
  if depth = 0 then 1 else (depth - 1).choose 2

def pairExcessMoment (v : Fin 13 → ℤ) (q : ℕ) : ℕ :=
  ∑ p ∈ Finset.Ioo 0 q,
    pairExcessTerm (LRC14Concrete.bandCount v q p)

def normalizedPairExcessMoment (v : Fin 13 → ℤ) (q : ℕ) : ℚ :=
  (pairExcessMoment v q : ℚ) / ((q : ℚ) - 1)

/-- The depth-at-least-three part of pair excess.  The term vanishes exactly
at depths zero, one, and two. -/
def shiftedPairDepthMoment (v : Fin 13 → ℤ) (q : ℕ) : ℕ :=
  ∑ p ∈ Finset.Ioo 0 q,
    (LRC14Concrete.bandCount v q p - 1).choose 2

theorem pairExcessMoment_eq_live_add_shifted
    (v : Fin 13 → ℤ) (q : ℕ) :
    pairExcessMoment v q =
      LRC14Concrete.liveCount v q + shiftedPairDepthMoment v q := by
  unfold pairExcessMoment shiftedPairDepthMoment LRC14Concrete.liveCount
  calc
    ∑ p ∈ Finset.Ioo 0 q,
        pairExcessTerm (LRC14Concrete.bandCount v q p) =
      ∑ p ∈ Finset.Ioo 0 q,
        ((if LRC14Concrete.bandCount v q p = 0 then 1 else 0) +
          (LRC14Concrete.bandCount v q p - 1).choose 2) := by
            apply Finset.sum_congr rfl
            intro p _
            by_cases hzero : LRC14Concrete.bandCount v q p = 0
            · simp [pairExcessTerm, hzero]
            · simp [pairExcessTerm, hzero]
    _ = (∑ p ∈ Finset.Ioo 0 q,
          if LRC14Concrete.bandCount v q p = 0 then 1 else 0) +
        ∑ p ∈ Finset.Ioo 0 q,
          (LRC14Concrete.bandCount v q p - 1).choose 2 := by
            rw [Finset.sum_add_distrib]
    _ = ((Finset.Ioo 0 q).filter fun p =>
          LRC14Concrete.bandCount v q p = 0).card +
        ∑ p ∈ Finset.Ioo 0 q,
          (LRC14Concrete.bandCount v q p - 1).choose 2 := by
            congr 1
            simp only [Finset.card_filter]

/-- Exact pointwise identity summed over multipliers:
`S₂ + (q-1) = S₁ + excess`. -/
theorem momentS_two_add_nonzero_eq_momentS_one_add_pairExcess
    (v : Fin 13 → ℤ) (q : ℕ) :
    LRC14Concrete.momentS v q 2 + (q - 1) =
      LRC14Concrete.momentS v q 1 + pairExcessMoment v q := by
  unfold LRC14Concrete.momentS pairExcessMoment
  have hones : q - 1 = ∑ _p ∈ Finset.Ioo 0 q, 1 := by
    simp [Nat.card_Ioo]
  rw [hones, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro p _
  have hdepth := LRC14Concrete.bandCount_le_thirteen v q p
  interval_cases LRC14Concrete.bandCount v q p <;>
    norm_num [pairExcessTerm, Nat.choose]

/-- Exact normalized coverage identity requested by the pair budget:
the pair moment equals the first moment minus one, plus live mass and the
depth-at-least-three shifted moment. -/
theorem normalizedPairDepthMoment_eq_first_sub_one_add_coverage_excess
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q) :
    normalizedPairDepthMoment v q =
      (LRC14Concrete.momentS v q 1 : ℚ) / ((q : ℚ) - 1) - 1 +
        ((LRC14Concrete.liveCount v q : ℚ) +
          (shiftedPairDepthMoment v q : ℚ)) / ((q : ℚ) - 1) := by
  have hid := momentS_two_add_nonzero_eq_momentS_one_add_pairExcess v q
  rw [pairExcessMoment_eq_live_add_shifted] at hid
  have hidQ : (LRC14Concrete.momentS v q 2 : ℚ) + ((q : ℚ) - 1) =
      (LRC14Concrete.momentS v q 1 : ℚ) +
        ((LRC14Concrete.liveCount v q : ℚ) +
          (shiftedPairDepthMoment v q : ℚ)) := by
    have hq1 : 1 ≤ q := le_of_lt hq
    exact_mod_cast hid
  have hden : ((q : ℚ) - 1) ≠ 0 := by
    have : (1 : ℚ) < q := by exact_mod_cast hq
    linarith
  dsimp [normalizedPairDepthMoment]
  field_simp [hden]
  linarith

/-- With exact singleton equilibrium, the horizon-thirty pair floor is
equivalent to an excess-overlap floor of exactly `443/1470`. -/
theorem pairDepthMoment_target_iff_excess_of_first_moment
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q)
    (hfirst : (LRC14Concrete.momentS v q 1 : ℚ) / ((q : ℚ) - 1) = 13 / 7) :
    1703 / 1470 ≤ normalizedPairDepthMoment v q ↔
      443 / 1470 ≤ normalizedPairExcessMoment v q := by
  have hid := momentS_two_add_nonzero_eq_momentS_one_add_pairExcess v q
  have hidQ : (LRC14Concrete.momentS v q 2 : ℚ) + ((q : ℚ) - 1) =
      (LRC14Concrete.momentS v q 1 : ℚ) + (pairExcessMoment v q : ℚ) := by
    have hq1 : 1 ≤ q := le_of_lt hq
    exact_mod_cast hid
  have hden : (0 : ℚ) < (q : ℚ) - 1 := by
    have : (1 : ℚ) < q := by exact_mod_cast hq
    linarith
  have hfirst' : (LRC14Concrete.momentS v q 1 : ℚ) =
      (13 / 7) * ((q : ℚ) - 1) := by
    rw [div_eq_iff (ne_of_gt hden)] at hfirst
    linarith
  dsimp [normalizedPairDepthMoment, normalizedPairExcessMoment]
  constructor <;> intro h
  · apply (le_div_iff₀ hden).2
    apply (le_div_iff₀ hden).1 at h
    nlinarith
  · apply (le_div_iff₀ hden).2
    apply (le_div_iff₀ hden).1 at h
    nlinarith

/-- Direct live/depth form of the exact `443/1470` shortfall. -/
theorem pairDepthMoment_target_iff_live_add_shifted_of_first_moment
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q)
    (hfirst : (LRC14Concrete.momentS v q 1 : ℚ) / ((q : ℚ) - 1) = 13 / 7) :
    1703 / 1470 ≤ normalizedPairDepthMoment v q ↔
      443 / 1470 ≤
        ((LRC14Concrete.liveCount v q : ℚ) +
          (shiftedPairDepthMoment v q : ℚ)) / ((q : ℚ) - 1) := by
  rw [pairDepthMoment_target_iff_excess_of_first_moment v q hq hfirst]
  simp only [normalizedPairExcessMoment,
    pairExcessMoment_eq_live_add_shifted, Nat.cast_add]

/-- The support-one aggregate is the sum of the singleton deviations. -/
theorem aggregateDeviation_one_eq_sum_singletons
    (v : Fin 13 → ℤ) (q : ℕ) :
    aggregateDeviation v q 1 =
      ∑ i, LRC14Concrete.deviation v q {i} := by
  rw [aggregateDeviation, Finset.powersetCard_one]
  simp

/-- At a modulus congruent to one modulo fourteen, every coprime singleton
deviation vanishes exactly. -/
theorem deviation_singleton_eq_zero_of_mod_fourteen_eq_one
    (v : Fin 13 → ℤ) (q : ℕ) (i : Fin 13)
    (hq : 14 ≤ q) (hmod : q % 14 = 1)
    (hgcd : Int.gcd (v i) (q : ℤ) = 1) :
    LRC14Concrete.deviation v q {i} = 0 := by
  have hq0 : 0 < q := by omega
  unfold LRC14Concrete.deviation
  rw [LRC14Concrete.jointFail_singleton_eq v q i hq0 hgcd,
    Finset.card_singleton, LRC14Concrete.bandSize_eq q hq]
  have hqform : q = 14 * (q / 14) + 1 := by omega
  have hband : 13 * q / 14 - (q + 13) / 14 + 1 = 12 * (q / 14) := by
    omega
  rw [hband]
  have hsub : q - 1 - 12 * (q / 14) = 2 * (q / 14) := by omega
  rw [hsub]
  have hqformQ : (q : ℚ) = 14 * ((q / 14 : ℕ) : ℚ) + 1 := by
    exact_mod_cast hqform
  push_cast
  rw [hqformQ]
  ring

/-- Exact singleton equilibrium at every coprime `1 mod 14` modulus. -/
theorem normalizedFirstDepthMoment_eq_thirteen_sevenths
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 14 ≤ q) (hmod : q % 14 = 1)
    (hcoprime : ∀ i, Int.gcd (v i) (q : ℤ) = 1) :
    (LRC14Concrete.momentS v q 1 : ℚ) / ((q : ℚ) - 1) = 13 / 7 := by
  have haggregate : aggregateDeviation v q 1 = 0 := by
    rw [aggregateDeviation_one_eq_sum_singletons]
    simp [deviation_singleton_eq_zero_of_mod_fourteen_eq_one
      v q _ hq hmod (hcoprime _)]
  have hmoment := aggregateDeviation_eq_moment_sub_equilibrium v q 1
  rw [haggregate] at hmoment
  have hden : ((q : ℚ) - 1) ≠ 0 := by
    have : (1 : ℚ) < q := by exact_mod_cast (show 1 < q by omega)
    linarith
  norm_num [Nat.choose] at hmoment ⊢
  field_simp [hden]
  linarith

/-- The clean modulus has exact singleton equilibrium. -/
theorem normalizedFirstDepthMoment_at_cleanModulus_eq_thirteen_sevenths
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (height : ℕ) :
    (LRC14Concrete.momentS v (cleanModulus v height) 1 : ℚ) /
        ((cleanModulus v height : ℚ) - 1) = 13 / 7 := by
  apply normalizedFirstDepthMoment_eq_thirteen_sevenths
  · exact fourteen_le_cleanModulus v height hv
  · exact cleanModulus_mod_fourteen v height
  · intro i
    exact int_gcd_speed_cleanModulus_eq_one v height i

/-- At the clean modulus, the pair floor is exactly a `443/1470` lower bound
for live mass plus depth-at-least-three pair excess. -/
theorem pairDepthMoment_target_at_cleanModulus_iff_coverage_excess
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (height : ℕ) :
    1703 / 1470 ≤ normalizedPairDepthMoment v (cleanModulus v height) ↔
      443 / 1470 ≤
        ((LRC14Concrete.liveCount v (cleanModulus v height) : ℚ) +
          (shiftedPairDepthMoment v (cleanModulus v height) : ℚ)) /
            ((cleanModulus v height : ℚ) - 1) := by
  exact pairDepthMoment_target_iff_live_add_shifted_of_first_moment
    v (cleanModulus v height) (one_lt_cleanModulus v height hv)
      (normalizedFirstDepthMoment_at_cleanModulus_eq_thirteen_sevenths
        v hv height)

/-- Nonpositive singleton deviations give a nonpositive normalized singleton
aggregate. -/
theorem normalizedAggregateDeviation_one_nonpos_of_singletons
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q)
    (hsingle : ∀ i, LRC14Concrete.deviation v q {i} ≤ 0) :
    normalizedAggregateDeviation v q 1 ≤ 0 := by
  have hsum : (∑ i, LRC14Concrete.deviation v q {i}) ≤ 0 :=
    Finset.sum_nonpos fun i _ => hsingle i
  have hden : (0 : ℚ) < (q : ℚ) - 1 := by
    have : (1 : ℚ) < q := by exact_mod_cast hq
    linarith
  rw [normalizedAggregateDeviation, aggregateDeviation_one_eq_sum_singletons]
  exact div_nonpos_of_nonpos_of_nonneg hsum hden.le

/-- THM-940 divided by `q-1`, with every support stratum exposed. -/
theorem normalized_B5_eq_equilibrium_add_aggregates
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q) :
    (LRC14Concrete.B5 v q : ℚ) / ((q : ℚ) - 1) =
      2052 / 16807 - normalizedAggregateDeviation v q 1 +
        normalizedAggregateDeviation v q 2 -
        normalizedAggregateDeviation v q 3 +
        normalizedAggregateDeviation v q 4 -
        normalizedAggregateDeviation v q 5 := by
  have hq1 : 1 ≤ q := le_of_lt hq
  have hqQ : ((q : ℚ) - 1) ≠ 0 := by
    have : (1 : ℚ) < q := by exact_mod_cast hq
    linarith
  have hledger :
      (∑ k ∈ range 6, (-1 : ℚ) ^ k *
        ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard k,
          LRC14Concrete.deviation v q T) =
        -aggregateDeviation v q 1 + aggregateDeviation v q 2 -
          aggregateDeviation v q 3 + aggregateDeviation v q 4 -
          aggregateDeviation v q 5 := by
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, pow_zero, pow_one,
      one_mul, neg_mul, Even.neg_one_pow (by decide : Even 2),
      Odd.neg_one_pow (by decide : Odd 3),
      Even.neg_one_pow (by decide : Even 4),
      Odd.neg_one_pow (by decide : Odd 5)]
    rw [show (∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 0,
        LRC14Concrete.deviation v q T) = 0 by
          simpa [aggregateDeviation] using aggregateDeviation_zero v q hq1]
    simp only [zero_add, aggregateDeviation]
    ring
  rw [LRC14Concrete.B5_eq_equilibrium_add_deviation, hledger]
  simp only [normalizedAggregateDeviation]
  field_simp [hqQ]
  ring

/-- The triangular transforms reproduce the exact THM-935 coefficients. -/
theorem aggregate_relationModel_identity
    (aggregate2 aggregate3 aggregate4 aggregate5 : ℚ) :
    2052 / 16807 + aggregate2 - aggregate3 + aggregate4 - aggregate5 =
      relationModelQ aggregate2
        (aggregate3 - (11 / 7) * aggregate2)
        (aggregate4 - (10 / 7) * aggregate3 + (55 / 49) * aggregate2)
        (aggregate5 - (9 / 7) * aggregate4 +
          (45 / 49) * aggregate3 - (165 / 343) * aggregate2) := by
  norm_num [relationModelQ]
  ring

/-- Favorable singleton deviations make the relation model a lower bound for
the normalized concrete B5 count. -/
theorem relationModelQ_le_normalized_B5_of_singletons_nonpos
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q)
    (hsingle : ∀ i, LRC14Concrete.deviation v q {i} ≤ 0) :
    relationModelQ (normalizedMass2 v q) (normalizedMass3 v q)
        (normalizedMass4 v q) (normalizedMass5 v q) ≤
      (LRC14Concrete.B5 v q : ℚ) / ((q : ℚ) - 1) := by
  have hA1 :=
    normalizedAggregateDeviation_one_nonpos_of_singletons v q hq hsingle
  rw [normalized_B5_eq_equilibrium_add_aggregates v q hq]
  have hid := aggregate_relationModel_identity
    (normalizedAggregateDeviation v q 2)
    (normalizedAggregateDeviation v q 3)
    (normalizedAggregateDeviation v q 4)
    (normalizedAggregateDeviation v q 5)
  simp only [normalizedMass2, normalizedMass3, normalizedMass4,
    normalizedMass5] at ⊢
  rw [← hid]
  linarith

/-- In raw THM-940 coordinates the harmful higher-support combination collapses
to one alternating expression. -/
theorem normalizedHigherMobiusSigned_expanded
    (v : Fin 13 → ℤ) (q : ℕ) :
    normalizedHigherMobiusSigned v q =
      normalizedAggregateDeviation v q 5 -
        normalizedAggregateDeviation v q 4 +
        normalizedAggregateDeviation v q 3 -
        (319 / 343) * normalizedAggregateDeviation v q 2 := by
  norm_num [normalizedHigherMobiusSigned, normalizedMass2, normalizedMass3,
    normalizedMass4, normalizedMass5]
  ring

/-- The harmful higher-support polynomial evaluated at one band depth. -/
def harmfulDepthTerm (depth : ℕ) : ℚ :=
  (depth.choose 5 : ℚ) - (depth.choose 4 : ℚ) +
    (depth.choose 3 : ℚ) - (319 / 343) * (depth.choose 2 : ℚ)

/-- Total harmful depth moment over the nonzero multipliers. -/
def harmfulDepthMoment (v : Fin 13 → ℤ) (q : ℕ) : ℚ :=
  ∑ p ∈ Finset.Ioo 0 q,
    harmfulDepthTerm (LRC14Concrete.bandCount v q p)

/-- Depths at most six contribute credit, never harmful debt. -/
theorem harmfulDepthTerm_nonpos_of_le_six
    (depth : ℕ) (hdepth : depth ≤ 6) :
    harmfulDepthTerm depth ≤ 0 := by
  interval_cases depth <;> norm_num [harmfulDepthTerm,
    show (4 : ℕ).choose 2 = 6 by decide,
    show (4 : ℕ).choose 3 = 4 by decide,
    show (5 : ℕ).choose 2 = 10 by decide,
    show (5 : ℕ).choose 3 = 10 by decide,
    show (5 : ℕ).choose 4 = 5 by decide,
    show (6 : ℕ).choose 2 = 15 by decide,
    show (6 : ℕ).choose 3 = 20 by decide,
    show (6 : ℕ).choose 4 = 15 by decide,
    show (6 : ℕ).choose 5 = 6 by decide]

/-- On the thirteen-runner range, depths at least seven contribute strictly
positive harmful debt. -/
theorem harmfulDepthTerm_pos_of_seven_le_of_le_thirteen
    (depth : ℕ) (hlower : 7 ≤ depth) (hupper : depth ≤ 13) :
    0 < harmfulDepthTerm depth := by
  interval_cases depth <;> norm_num [harmfulDepthTerm,
    show (7 : ℕ).choose 2 = 21 by decide,
    show (7 : ℕ).choose 3 = 35 by decide,
    show (7 : ℕ).choose 4 = 35 by decide,
    show (7 : ℕ).choose 5 = 21 by decide,
    show (8 : ℕ).choose 2 = 28 by decide,
    show (8 : ℕ).choose 3 = 56 by decide,
    show (8 : ℕ).choose 4 = 70 by decide,
    show (8 : ℕ).choose 5 = 56 by decide,
    show (9 : ℕ).choose 2 = 36 by decide,
    show (9 : ℕ).choose 3 = 84 by decide,
    show (9 : ℕ).choose 4 = 126 by decide,
    show (9 : ℕ).choose 5 = 126 by decide,
    show (10 : ℕ).choose 2 = 45 by decide,
    show (10 : ℕ).choose 3 = 120 by decide,
    show (10 : ℕ).choose 4 = 210 by decide,
    show (10 : ℕ).choose 5 = 252 by decide,
    show (11 : ℕ).choose 2 = 55 by decide,
    show (11 : ℕ).choose 3 = 165 by decide,
    show (11 : ℕ).choose 4 = 330 by decide,
    show (11 : ℕ).choose 5 = 462 by decide,
    show (12 : ℕ).choose 2 = 66 by decide,
    show (12 : ℕ).choose 3 = 220 by decide,
    show (12 : ℕ).choose 4 = 495 by decide,
    show (12 : ℕ).choose 5 = 792 by decide,
    show (13 : ℕ).choose 2 = 78 by decide,
    show (13 : ℕ).choose 3 = 286 by decide,
    show (13 : ℕ).choose 4 = 715 by decide,
    show (13 : ℕ).choose 5 = 1287 by decide]

theorem harmfulDepthMoment_eq_moment_combination
    (v : Fin 13 → ℤ) (q : ℕ) :
    harmfulDepthMoment v q =
      (LRC14Concrete.momentS v q 5 : ℚ) -
        (LRC14Concrete.momentS v q 4 : ℚ) +
        (LRC14Concrete.momentS v q 3 : ℚ) -
        (319 / 343) * (LRC14Concrete.momentS v q 2 : ℚ) := by
  simp only [harmfulDepthMoment, harmfulDepthTerm, LRC14Concrete.momentS,
    Nat.cast_sum]
  rw [Finset.sum_sub_distrib, Finset.sum_add_distrib,
    Finset.sum_sub_distrib, ← Finset.mul_sum]

/-- Exact leverage/depth-spectrum form of the signed higher mass. -/
theorem normalizedHigherMobiusSigned_eq_depthMoment
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q) :
    normalizedHigherMobiusSigned v q =
      harmfulDepthMoment v q / ((q : ℚ) - 1) + 14586 / 16807 := by
  rw [normalizedHigherMobiusSigned_expanded,
    harmfulDepthMoment_eq_moment_combination]
  simp only [normalizedAggregateDeviation]
  rw [aggregateDeviation_eq_moment_sub_equilibrium,
    aggregateDeviation_eq_moment_sub_equilibrium,
    aggregateDeviation_eq_moment_sub_equilibrium,
    aggregateDeviation_eq_moment_sub_equilibrium]
  have hden : ((q : ℚ) - 1) ≠ 0 := by
    have : (1 : ℚ) < q := by exact_mod_cast hq
    linarith
  norm_num [Nat.choose]
  field_simp [hden]
  ring

/-- The signed remaining budget is exactly a negative average-depth-moment
target. -/
theorem signed_remaining_iff_depthMoment
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q) :
    normalizedHigherMobiusSigned v q < 7712 / 84035 ↔
      harmfulDepthMoment v q / ((q : ℚ) - 1) < -(65218 / 84035) := by
  rw [normalizedHigherMobiusSigned_eq_depthMoment v q hq]
  constructor <;> intro h <;> linarith

/-- The sharp rational horizon-thirty budget. -/
theorem normalized_relationModelQ_pos_of_signed_remaining_budget
    (v : Fin 13 → ℤ) (q : ℕ)
    (hpair : -(13 / 30 : ℚ) ≤ normalizedMass2 v q)
    (hhigher : normalizedHigherMobiusSigned v q < 7712 / 84035) :
    0 < relationModelQ (normalizedMass2 v q) (normalizedMass3 v q)
      (normalizedMass4 v q) (normalizedMass5 v q) := by
  dsimp [relationModelQ, normalizedHigherMobiusSigned] at ⊢ hhigher
  nlinarith

/-- End-to-end THM-940 consumer with singleton sign, pair lower tail, and the
signed higher-support remainder. -/
theorem B5_pos_of_normalized_signed_budget
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 < q)
    (hsingle : ∀ i, LRC14Concrete.deviation v q {i} ≤ 0)
    (hpair : -(13 / 30 : ℚ) ≤ normalizedMass2 v q)
    (hhigher : normalizedHigherMobiusSigned v q < 7712 / 84035) :
    0 < LRC14Concrete.B5 v q := by
  have hmodel := normalized_relationModelQ_pos_of_signed_remaining_budget
    v q hpair hhigher
  have hlower := relationModelQ_le_normalized_B5_of_singletons_nonpos
    v q hq hsingle
  have hnorm :
      (0 : ℚ) < (LRC14Concrete.B5 v q : ℚ) / ((q : ℚ) - 1) :=
    lt_of_lt_of_le hmodel hlower
  have hden : (0 : ℚ) < (q : ℚ) - 1 := by
    have : (1 : ℚ) < q := by exact_mod_cast hq
    linarith
  have hB5Q : (0 : ℚ) < (LRC14Concrete.B5 v q : ℚ) := by
    rcases div_pos_iff.mp hnorm with hpositive | hnegative
    · exact hpositive.1
    · exact absurd hnegative.2 (not_lt_of_ge hden.le)
  exact_mod_cast hB5Q

/-- S36 discharges the singleton premise automatically at any modulus at least
fourteen which is coprime to every speed. -/
theorem B5_pos_of_coprime_normalized_signed_budget
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 14 ≤ q)
    (hcoprime : ∀ i, Int.gcd (v i) (q : ℤ) = 1)
    (hpair : -(13 / 30 : ℚ) ≤ normalizedMass2 v q)
    (hhigher : normalizedHigherMobiusSigned v q < 7712 / 84035) :
    0 < LRC14Concrete.B5 v q := by
  apply B5_pos_of_normalized_signed_budget v q (by omega)
  · intro i
    exact (LRC14Concrete.deviation_singleton_bounds v q i hq
      (hcoprime i)).2
  · exact hpair
  · exact hhigher

/-- End-to-end version exposing the true depth-seven crux instead of the
Möbius-mass coordinates. -/
theorem B5_pos_of_coprime_pair_and_depth_budget
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 14 ≤ q)
    (hcoprime : ∀ i, Int.gcd (v i) (q : ℤ) = 1)
    (hpair : -(13 / 30 : ℚ) ≤ normalizedMass2 v q)
    (hdepth : harmfulDepthMoment v q / ((q : ℚ) - 1) <
      -(65218 / 84035)) :
    0 < LRC14Concrete.B5 v q := by
  apply B5_pos_of_coprime_normalized_signed_budget v q hq hcoprime hpair
  exact (signed_remaining_iff_depthMoment v q (by omega)).2 hdepth

/-- Canonical clean-modulus consumer.  Once a coefficient height is chosen,
coprimality and the singleton rung are automatic; only the pair and depth
budgets at that explicit modulus remain. -/
theorem B5_pos_at_cleanModulus_of_pair_and_depth_budget
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (height : ℕ)
    (hpair : -(13 / 30 : ℚ) ≤
      normalizedMass2 v (cleanModulus v height))
    (hdepth :
      harmfulDepthMoment v (cleanModulus v height) /
          (((cleanModulus v height : ℕ) : ℚ) - 1) <
        -(65218 / 84035)) :
    0 < LRC14Concrete.B5 v (cleanModulus v height) := by
  apply B5_pos_of_coprime_pair_and_depth_budget
    v (cleanModulus v height) (fourteen_le_cleanModulus v height hv)
  · intro i
    exact int_gcd_speed_cleanModulus_eq_one v height i
  · exact hpair
  · exact hdepth

/-- Histogram-only clean-modulus consumer.  Both remaining inputs are now
moments of `bandCount`: the second factorial moment floor and the signed
depth-seven tail polynomial. -/
theorem B5_pos_at_cleanModulus_of_depthMoment_budgets
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (height : ℕ)
    (hpair : 1703 / 1470 ≤
      normalizedPairDepthMoment v (cleanModulus v height))
    (hdepth :
      harmfulDepthMoment v (cleanModulus v height) /
          (((cleanModulus v height : ℕ) : ℚ) - 1) <
        -(65218 / 84035)) :
    0 < LRC14Concrete.B5 v (cleanModulus v height) := by
  apply B5_pos_at_cleanModulus_of_pair_and_depth_budget v hv height
  · exact (normalizedMass2_lower_iff_pairDepthMoment
      v (cleanModulus v height)
      (one_lt_cleanModulus v height hv)).2 hpair
  · exact hdepth

/-- Final coverage-language clean-modulus consumer.  The pair side is no
longer an opaque moment: it is exactly live mass plus depth-at-least-three
pair excess above the sharp `443/1470` floor. -/
theorem B5_pos_at_cleanModulus_of_coverage_excess_and_depth_budget
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (height : ℕ)
    (hcoverage : 443 / 1470 ≤
      ((LRC14Concrete.liveCount v (cleanModulus v height) : ℚ) +
        (shiftedPairDepthMoment v (cleanModulus v height) : ℚ)) /
          ((cleanModulus v height : ℚ) - 1))
    (hdepth :
      harmfulDepthMoment v (cleanModulus v height) /
          (((cleanModulus v height : ℕ) : ℚ) - 1) <
        -(65218 / 84035)) :
    0 < LRC14Concrete.B5 v (cleanModulus v height) := by
  apply B5_pos_at_cleanModulus_of_depthMoment_budgets v hv height
  · exact (pairDepthMoment_target_at_cleanModulus_iff_coverage_excess
      v hv height).2 hcoverage
  · exact hdepth

/-! ## Axiom audit -/

#print axioms normalized_B5_eq_equilibrium_add_aggregates
#print axioms aggregate_relationModel_identity
#print axioms normalizedMass2_lower_iff_pairDepthMoment
#print axioms relationModelQ_le_normalized_B5_of_singletons_nonpos
#print axioms normalizedHigherMobiusSigned_expanded
#print axioms harmfulDepthTerm_nonpos_of_le_six
#print axioms harmfulDepthTerm_pos_of_seven_le_of_le_thirteen
#print axioms normalizedHigherMobiusSigned_eq_depthMoment
#print axioms signed_remaining_iff_depthMoment
#print axioms B5_pos_of_normalized_signed_budget
#print axioms B5_pos_of_coprime_normalized_signed_budget
#print axioms B5_pos_of_coprime_pair_and_depth_budget
#print axioms B5_pos_at_cleanModulus_of_pair_and_depth_budget
#print axioms B5_pos_at_cleanModulus_of_depthMoment_budgets
#print axioms B5_pos_at_cleanModulus_of_coverage_excess_and_depth_budget

end
end LRCB5NormalizedBridge
end LonelyRunner
