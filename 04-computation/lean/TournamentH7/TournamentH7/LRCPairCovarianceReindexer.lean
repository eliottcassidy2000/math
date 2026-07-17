import TournamentH7.LRCOpenPairLedger
import TournamentH7.LRCB5PairOverlapSum
import TournamentH7.LRCPairCovarianceKernel

/-!
Adapter from the single remaining finite comb reindexing identity to the exact
circle pair covariance formula.  All measure, merge, and Bernoulli arithmetic
has already been discharged before this interface.
-/

namespace LonelyRunner
namespace LRCPairCovarianceReindexer

open RatIntervals LRCRationalOpenComb
open LRCOpenPairLedger LRCB5PairOverlapSum
open LRCB5PairGridBridge LRCPairCovarianceKernel

/-- The sole finite geometric identity left in the pair covariance producer. -/
def PairOverlapReindexing (first second : ℤ) : Prop :=
  ((length (ratOpenPairRegion first.natAbs second.natAbs) : ℚ) : ℝ) =
    scaledPairOverlapLedger (Nat.gcd first.natAbs second.natAbs)
      (first.natAbs / Nat.gcd first.natAbs second.natAbs)
      (second.natAbs / Nat.gcd first.natAbs second.natAbs)

/-- Once the explicit finite clip-sum reindexer is supplied, the circle pair
correlation is exactly the rational Bernoulli kernel. -/
theorem pairContinuumCorrelation_sub_eq_pairKernel_of_reindexing
    (v : Fin 13 → ℤ) (first second : Fin 13)
    (hfirst : v first ≠ 0) (hsecond : v second ≠ 0)
    (hreindex : PairOverlapReindexing (v first) (v second)) :
    pairContinuumCorrelation v {first, second} - 1 / 49 =
      (pairKernel (v first) (v second) : ℝ) := by
  rw [pairContinuumCorrelation_pair_eq_length v first second hfirst hsecond]
  rw [length_ratOpenPairRegion_sub_eq_B2_of_scaledLedger
    (v first).natAbs (v second).natAbs
    (Int.natAbs_pos.mpr hfirst) (Int.natAbs_pos.mpr hsecond) hreindex]
  rw [pairKernel_cast]
  rfl

/-- Family-level form matching the producer premise of
`correlation_lower_of_pairKernel_identity`. -/
def PairOverlapReindexingFor (v : Fin 13 → ℤ) : Prop :=
  ∀ first second, first ≠ second →
    PairOverlapReindexing (v first) (v second)

theorem pairKernel_identity_of_reindexing
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (hreindex : PairOverlapReindexingFor v) :
    ∀ first second, first ≠ second →
      pairContinuumCorrelation v {first, second} - 1 / 49 =
        (pairKernel (v first) (v second) : ℝ) := by
  intro first second hneq
  exact pairContinuumCorrelation_sub_eq_pairKernel_of_reindexing
    v first second (hv first) (hv second) (hreindex first second hneq)

#print axioms pairContinuumCorrelation_sub_eq_pairKernel_of_reindexing
#print axioms pairKernel_identity_of_reindexing

end LRCPairCovarianceReindexer
end LonelyRunner
