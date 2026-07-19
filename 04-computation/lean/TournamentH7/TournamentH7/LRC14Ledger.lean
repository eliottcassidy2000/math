/- LRC14Ledger.lean — mac-mini-2026-07-16-S129.
   THE ASSEMBLY FILE of the LRC(14) formalization ladder.
   Imports every built rung, states the LRC(14) target as a Prop (no axiom, no
   sorry), and records the mapping from the target to the canon chain.

   BUILT RUNGS (oleans, kernel-checked):
   · FragmentationCount — THM-883's kernel (arc count; fragmentation measure bound)
   · KillerBudget       — the budget composition + the good-set floor
   · TrivialLoneliness  — END-TO-END existence: `trivial_LRC` proves the Lonely
                          Runner statement at the trivial constant 1/(6k)
   · TieSplitWalk       — THM-866's kernel (F3 flip arithmetic; pigeonhole)

   THE TARGET (stated below, not yet formally discharged): `LRC14` — thirteen
   positive speeds admit a 1/14-lonely time. The informal proof is the canon
   chain (covering: THM-922 sign-off + THM-883 rigidity + bands + low-M;
   S3 residual: THM-527 reformulation + THM-924 glue + THM-925/926 floor
   c₀ = 17/84; level-5 wall: THM-896/897 + THM-926 T2). The formal gap between
   `LRC13speeds_at_gap_26` and `LRC14` is the remaining constant improvement
   1/26 → 1/14; each canon page is finite/rational and slots
   into this namespace as a further rung. -/
import TournamentH7.FragmentationCount
import TournamentH7.KillerBudget
import TournamentH7.TrivialLoneliness
import TournamentH7.TieSplitWalk
import TournamentH7.UnitBudgetEndpoint
import TournamentH7.LRCMod19LedgerBridge

open MeasureTheory Set

namespace LRC14

/-- **Uniqueness-axis rung (n=12 gap, problem (C)).**  Re-export of the mod-19 antipodal-spread
    constraint wired into this ledger via `LRCMod19LedgerBridge`: any family in the mod-19 gap
    regime (margin `< 2/19` at every rational time `b/19`, implied by `M(C) < 2/19`) with no speed
    divisible by 19 has residues covering every antipodal unit-pair of `ℤ/19`.  A proved necessary
    condition on any family populating the uniqueness gap, on the same (C)-axis as `LRCLadderD1`'s
    `2/25` reach bound (see `antipodal_cover_of_margin`). -/
theorem gap_regime_mod19_spread {k : ℕ} [Nonempty (Fin k)] (v : Fin k → ℤ)
    (hunit : ∀ i, ¬ ((19 : ℤ) ∣ v i))
    (hgap : ∀ b : ℤ, TournamentH7.LRCWitness.margin v ((b : ℝ) / 19) < 2 / 19)
    (u : ZMod 19) (hu : u ≠ 0) :
    ∃ i, (v i : ZMod 19) = u ∨ (v i : ZMod 19) = -u :=
  antipodal_cover_of_margin v hunit hgap u hu

/-- **The LRC(14) target statement.** Thirteen distinct positive integer speeds
    admit a time `t` with `‖w·t‖ ≥ 1/14` for every speed — stated via integer
    distance, matching `trivial_LRC`'s conclusion shape. -/
def LRC14 : Prop :=
  ∀ W : Finset ℕ, (∀ w ∈ W, 0 < w) → W.card = 13 →
    ∃ t ∈ Icc (0 : ℝ) 1, ∀ w ∈ W, ∀ a : ℤ, (1 : ℝ)/14 ≤ |(w : ℝ) * t - a|

/-- The ladder's current unconditional instance: the same statement at the
    trivial constant — for thirteen speeds, loneliness at gap `1/79` (any
    `lam < 1/78` works; `1/79` keeps the arithmetic decidable by `norm_num`). -/
theorem LRC13speeds_at_trivial_gap :
    ∀ W : Finset ℕ, (∀ w ∈ W, 0 < w) → W.card = 13 →
      ∃ t ∈ Icc (0 : ℝ) 1, ∀ w ∈ W, ∀ a : ℤ, (1 : ℝ)/79 ≤ |(w : ℝ) * t - a| := by
  intro W hW hcard
  refine trivial_LRC W hW (1/79) (by norm_num) (by norm_num) ?_
  rw [hcard]
  norm_num

end LRC14

namespace LRC14

/-- **Thirteen speeds at gap 1/27** — a strict-budget rational anchor retained
    for computations.  `LRC13speeds_at_gap_26` is the stronger compactness
    endpoint; the remaining program is `1/26 → 1/14`. -/
theorem LRC13speeds_at_gap_27 :
    ∀ W : Finset ℕ, (∀ w ∈ W, 0 < w) → W.card = 13 →
      ∃ t ∈ Icc (0 : ℝ) 1, ∀ w ∈ W, ∀ a : ℤ, (1 : ℝ)/27 ≤ |(w : ℝ) * t - a| := by
  intro W hW hcard
  refine exists_lonely_sharp W hW (1/27) (by norm_num) (by norm_num) ?_
  rw [hcard]
  norm_num

end LRC14
