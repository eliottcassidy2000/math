import Mathlib

/-!
# The q=14 beat-puncture integer core (THM-1215)

This module formalizes the finite arithmetic contradiction used by THM-1215.
For a consecutive block of `N` beat numerators, the THM-1192 cap at reduced
period 14 is

`U(N,14) = N / 14 + min (N % 14) 1 = ceil(N/14)`.

If `N >= 6`, then `N > 5 U(N,14)`.  Consequently the phase-free necessary
law `N - U <= 4 U` cannot hold.

The geometric and modular bridge is deliberately explicit in the hypotheses
of `q14_stalk_bridge_impossible`: `hN` is supplied by the slow-gap length
`q >= 7a`; `hdef` and `hcomp` are supplied by the five gcd identities
`gcd(c,q)=d` with `q=14d`; and `hlaw` is THM-1192's phase-free necessary law.
This file does not silently assume the slow-gap cover, identify its integer
block, or prove the gcd reduction from the original speeds.
-/

namespace LRC14.BeatPunctureQ14

/-- The natural-number ceiling of `N / 14`. -/
def ceilDiv14 (N : ℕ) : ℕ :=
  (N + 13) / 14

/-- THM-1192's full-period-plus-remainder cap specialized to `Q=14`, where
the strict dangerous window contains exactly one residue per period. -/
def period14Cap (N : ℕ) : ℕ :=
  N / 14 + min (N % 14) 1

/-- The THM-1192 period cap at `Q=14` is exactly `ceil(N/14)`. -/
theorem period14Cap_eq_ceilDiv14 (N : ℕ) :
    period14Cap N = ceilDiv14 N := by
  unfold period14Cap ceilDiv14
  omega

/-- Sharp integer core: from the first possible slow-gap block size onward,
five copies of the q=14 cap are strictly smaller than the block. -/
theorem five_ceilDiv14_lt_of_six_le (N : ℕ) (hN : 6 ≤ N) :
    5 * ceilDiv14 N < N := by
  unfold ceilDiv14
  omega

/-- Equivalent version stated directly with THM-1192's cap formula. -/
theorem five_period14Cap_lt_of_six_le (N : ℕ) (hN : 6 ≤ N) :
    5 * period14Cap N < N := by
  rw [period14Cap_eq_ceilDiv14]
  exact five_ceilDiv14_lt_of_six_le N hN

/-- The q=14 specialization of THM-1192's phase-free law is false for every
block of at least six consecutive beat numerators. -/
theorem not_phaseFreeLaw_of_six_le (N : ℕ) (hN : 6 ≤ N) :
    ¬(N - period14Cap N ≤ 4 * period14Cap N) := by
  intro hlaw
  have hstrict := five_period14Cap_lt_of_six_le N hN
  omega

/-- Typed consumer for the geometric/gcd bridge.  The defining comb has cap
`U`; the four complementary combs have total capacity at most `4U`; and the
THM-1192 necessary law asks the surviving defining-pair block to fit inside
that total.  These explicit provider hypotheses are inconsistent once the
geometric block supplier gives `N >= 6`. -/
theorem q14_stalk_bridge_impossible
    (N definingCap complementaryTotal : ℕ)
    (hN : 6 ≤ N)
    (hdef : definingCap = period14Cap N)
    (hcomp : complementaryTotal ≤ 4 * period14Cap N)
    (hlaw : N - definingCap ≤ complementaryTotal) : False := by
  subst definingCap
  have hphase : N - period14Cap N ≤ 4 * period14Cap N :=
    hlaw.trans hcomp
  exact not_phaseFreeLaw_of_six_le N hN hphase

#print axioms period14Cap_eq_ceilDiv14
#print axioms five_ceilDiv14_lt_of_six_le
#print axioms not_phaseFreeLaw_of_six_le
#print axioms q14_stalk_bridge_impossible

end LRC14.BeatPunctureQ14
