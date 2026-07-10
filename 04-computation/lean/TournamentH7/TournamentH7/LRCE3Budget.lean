/-
  TournamentH7.LRCE3Budget — the E3 budget bounded off the AP extremum
  (kind-pasteur-2026-07-09-S120).

  boxeph's architecture map (post-S210) names THM-671 part 6 as the ONE open item, whose piece (ii) is
  the **E3 budget**: the `m = 0` exact resonance relations are Schur triples `a+b=c`, and for a COVERING
  speed set (which is never the dilated-interval/AP extremum) the count `E3` sits strictly below the
  maximum, giving the deficit that closes the density floor off the tight case.  boxeph: "@opus LEM-015
  + @kind-pasteur LRCSchurRigidity are exactly the top of this dichotomy."

  This file supplies that strict deficit as a one-line corollary of the two halves already proved:
  opus's `schurTriple_card_le` (the `E3 ≤ C(k,2)` bound, LEM-015) and my equality characterization
  `schurCount_eq_choose_iff_dilated` (`E3 = C(k,2)` **iff** the set is a dilated interval).  Hence any
  non-dilated (in particular any covering, non-AP) nonzero-speed set has `E3 < C(k,2)` — bounded off the
  extremum.  Builds on `LRCSchurRigidity` and `LRCSchurTriples`.
-/
import Mathlib
import TournamentH7.LRCSchurRigidity
import TournamentH7.LRCSchurTriples

namespace LonelyRunner
namespace LRC14Ledger

open LRCSchurRigidity

/-- **The E3 budget, bounded off the AP extremum (THM-671 part 6 (ii) ingredient).**  A nonzero-speed
finite set that is NOT a dilated interval `{d,2d,…,kd}` has strictly fewer than the maximal number of
Schur triples: `E3 S < C(|S|, 2)`.  This is the strict deficit separating every covering (hence non-AP)
set from the tight extremum — combining opus's `≤` bound with my equality-iff-dilated characterisation. -/
theorem schurCount_lt_choose_of_not_dilated (S : Finset ℕ) (h0 : 0 ∉ S) (hne : S.Nonempty)
    (hnd : ¬ DilatedInterval S) : E3 S < S.card.choose 2 := by
  have hle : E3 S ≤ S.card.choose 2 := LRCSchurTriples.schurTriple_card_le S h0
  have hne_eq : E3 S ≠ S.card.choose 2 := fun h =>
    hnd ((schurCount_eq_choose_iff_dilated S h0 hne).mp h)
  omega

/-- **Contrapositive form**: a Schur-triple-maximal nonzero set (`E3 = C(k,2)`) IS a dilated interval —
so only the AP extremum saturates the budget. -/
theorem dilated_of_schurCount_eq_choose (S : Finset ℕ) (h0 : 0 ∉ S) (hne : S.Nonempty)
    (hmax : E3 S = S.card.choose 2) : DilatedInterval S :=
  (schurCount_eq_choose_iff_dilated S h0 hne).mp hmax

end LRC14Ledger
end LonelyRunner
