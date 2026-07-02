/-
klein-2026-07-02-S105 — the module-6 fuel-checker INSTANTIATION exemplar (HYP-4010).

Concrete `LadderPack` data through the decidable gate: the DEEP WELL {1,…,12,182}
(= 183 − 1), the covering-min extremal family. Window [1/14 + 1/500, 13/168] around the
binding witness t* = 14/183; h = 1/14 (the gate band), μ = 1/2200. `ladderOK` closes by
`decide` — one kernel check certifies the whole family — and the dispatcher exports the
skeleton-shaped loneliness of the deep-well vector. This is the instantiation pattern
for every remaining pack: data + one decide + one `lonely_of_ladder_mem`.
-/
import Mathlib
import TournamentH7.LRC14Dispatch

namespace LadderPackData

open LonelyRunner LonelyRunner.LRC14 WitnessCert

/-- The deep-well pack: bounded part {1..12}, one level V = 183 with offset 1
(speed 182), window around t* = 14/183. -/
def deepWellPack : LadderPack :=
  { h := 1/14, μ := 1/2200, lo := 1/14 + 1/2200, hi := 13/168
    P := [1,2,3,4,5,6,7,8,9,10,11,12]
    levels := [⟨[1], 0, 183⟩] }

/-- The fuel-checker gate passes on the deep well: one kernel `decide`. -/
theorem deepWellPack_ok : ladderOK deepWellPack := by native_decide

/-- The deep-well speed vector. -/
def deepWellV : Fin 13 → ℤ := ![1,2,3,4,5,6,7,8,9,10,11,12,182]

theorem deepWellV_mem : ∀ i, deepWellV i ∈ deepWellPack.speeds := by decide

/-- **The instantiated dispatcher**: the deep well is lonely — module 6, through the
gate, in the skeleton's shape. -/
theorem deepWell_lonely : ∃ t : ℝ, Lonely 14 deepWellV t :=
  lonely_of_ladder_mem deepWellPack deepWellPack_ok deepWellV deepWellV_mem

end LadderPackData
