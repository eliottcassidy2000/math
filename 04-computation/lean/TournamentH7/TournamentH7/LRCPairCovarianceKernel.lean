import TournamentH7.LRCClosedBudget
import TournamentH7.LRCB5PairGridBridge

/-!
Rational residue-14 Bernoulli kernel for the LRC(14) pair covariance.  All
casting, negative-part, and support bookkeeping is discharged here; the only
producer-facing analytic premise is the exact two-runner danger correlation
identity stated by `correlation_lower_of_pairKernel_identity`.
-/

namespace LonelyRunner
namespace LRCPairCovarianceKernel

open LRC14.ClosedBudget
open LRCB5PairGridBridge

/-- `B₂` at the fractional residue `k/14`, kept rational. -/
def bernoulliResidue14 (k : ℤ) : ℚ :=
  let residue : ℚ := (k % 14 : ℤ)
  residue ^ 2 / 196 - residue / 14 + 1 / 6

theorem B2R_int_div_fourteen (k : ℤ) :
    B2R ((k : ℝ) / 14) = (bernoulliResidue14 k : ℝ) := by
  unfold B2R bernoulliResidue14
  have hfract : Int.fract ((k : ℝ) / 14) =
      ((k % (14 : ℕ) : ℤ) : ℝ) / 14 :=
    Int.fract_div_intCast_eq_div_intCast_mod
  rw [hfract]
  push_cast
  ring

def reducedFirst (first second : ℤ) : ℕ :=
  first.natAbs / Nat.gcd first.natAbs second.natAbs

def reducedSecond (first second : ℤ) : ℕ :=
  second.natAbs / Nat.gcd first.natAbs second.natAbs

/-- Exact rational primitive-ratio covariance kernel. -/
def pairKernel (first second : ℤ) : ℚ :=
  let reducedFirst := reducedFirst first second
  let reducedSecond := reducedSecond first second
  (bernoulliResidue14 ((reducedFirst : ℤ) - reducedSecond) -
      bernoulliResidue14 ((reducedFirst : ℤ) + reducedSecond)) /
    ((reducedFirst : ℚ) * reducedSecond)

theorem pairKernel_cast (first second : ℤ) :
    (pairKernel first second : ℝ) =
      (B2R ((((reducedFirst first second : ℕ) : ℝ) -
          (reducedSecond first second : ℕ)) / 14) -
        B2R ((((reducedFirst first second : ℕ) : ℝ) +
          (reducedSecond first second : ℕ)) / 14)) /
      ((reducedFirst first second : ℕ) *
        (reducedSecond first second : ℕ)) := by
  rw [show (((reducedFirst first second : ℕ) : ℝ) -
        (reducedSecond first second : ℕ)) / 14 =
      (((reducedFirst first second : ℕ) : ℤ) -
        (reducedSecond first second : ℕ) : ℤ) / 14 by push_cast; ring]
  rw [show (((reducedFirst first second : ℕ) : ℝ) +
        (reducedSecond first second : ℕ)) / 14 =
      (((reducedFirst first second : ℕ) : ℤ) +
        (reducedSecond first second : ℕ) : ℤ) / 14 by push_cast; ring]
  rw [B2R_int_div_fourteen, B2R_int_div_fourteen]
  unfold pairKernel
  push_cast
  rfl

/-- The nonnegative rational edge weight used by the ratio-layer theorem. -/
def pairDeficit (first second : ℤ) : ℚ :=
  max 0 (-pairKernel first second)

theorem neg_pairDeficit_le_pairKernel (first second : ℤ) :
    -(pairDeficit first second : ℝ) ≤ (pairKernel first second : ℝ) := by
  have hq : -pairDeficit first second ≤ pairKernel first second := by
    unfold pairDeficit
    simpa using neg_le_neg (le_max_right (0 : ℚ) (-pairKernel first second))
  exact_mod_cast hq

abbrev PairSupport := Finset (Fin 13)

def pairSupports : Finset PairSupport :=
  (Finset.univ : Finset (Fin 13)).powersetCard 2

/-- Once the single two-runner danger-set identity is available, the rational
negative part gives exactly the termwise minorant required by the continuum
floor adapter. -/
theorem correlation_lower_of_pairKernel_identity
    (v : Fin 13 → ℤ) (weight : PairSupport → ℚ)
    (hweight : ∀ first second, first ≠ second →
      weight {first, second} = pairDeficit (v first) (v second))
    (hpair : ∀ first second, first ≠ second →
      pairContinuumCorrelation v {first, second} - 1 / 49 =
        (pairKernel (v first) (v second) : ℝ)) :
    ∀ support ∈ pairSupports,
      -(weight support : ℝ) ≤
        pairContinuumCorrelation v support - 1 / 49 := by
  intro support hsupport
  rw [pairSupports, Finset.mem_powersetCard] at hsupport
  obtain ⟨_, hcard⟩ := hsupport
  obtain ⟨first, second, hne, rfl⟩ := Finset.card_eq_two.mp hcard
  rw [hweight first second hne, hpair first second hne]
  exact neg_pairDeficit_le_pairKernel (v first) (v second)

#print axioms B2R_int_div_fourteen
#print axioms pairKernel_cast
#print axioms neg_pairDeficit_le_pairKernel
#print axioms correlation_lower_of_pairKernel_identity

end LRCPairCovarianceKernel
end LonelyRunner
