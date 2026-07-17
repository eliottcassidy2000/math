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
   `trivial_LRC` and `LRC14` is exactly the certified chain's constant
   improvement 1/(6·13) → 1/14; each canon page is finite/rational and slots
   into this namespace as a further rung. -/
import TournamentH7.FragmentationCount
import TournamentH7.KillerBudget
import TournamentH7.TrivialLoneliness
import TournamentH7.TieSplitWalk

open Set

namespace LRC14

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
