/-
  TournamentH7.LRCMultiFoldRows — THE CONSECUTIVE MULTI-FOLD TOWER FLOOR ROWS
  (kind-pasteur-2026-07-05-S17, HYP-4217; the kps side of the HYP-4212 lift
  domination co-authorship with mac-mini).

  The consecutive multi-fold towers `D_l = {1..12−l} ∪ {14r : r = 13−l..12}`
  (the deep well `D_1 = {1..11,168}` iterated: the top `l` residues lifted to
  their `14r` rungs).  CORRECTION LOGGED THIS SESSION (under HYP-4177): the
  closed-form law `M(D_l) = 14/(14(13−l)+1)` is FALSE at `l = 4, 5` — exact
  enumeration gives `M(D_4) = 17/155` and `M(D_5) = 19/155` (both BELOW the
  claimed values; the binding pair migrates to the `154+1` grid mid-ladder;
  `l = 1, 2, 3, 6` match the law).  WHAT THE DOMINATION ASSEMBLY NEEDS is only
  the FLOOR, and it survives everywhere:

      `M(D_l) ≥ 2/25`  for all `l = 1..6`  (the dichotomy's β, HYP-4212's
      corrected floor),

  certified here by five `rational_point_margin` rows (HYP-4102 atom, all
  hypotheses by `decide`) plus `≥ 2/25` corollaries:

    l=2: (k,s,μ) = (14,155,14)   l=3: (14,141,14)   l=4: (17,155,17)
    l=5: (19,155,19)             l=6: (14, 99,14)
    (l=1 is `ladder_margin_r12` in `LRCLiftFloorRows` — the deep well.)

  Kernel-pure (`decide` on integer mod conditions only; no `native_decide`).
-/
import TournamentH7.LRCHarmonicGate
import TournamentH7.LRCLiftFloorRows

namespace LonelyRunner
namespace MultiFoldRows

open HarmonicGate

/-- The consecutive multi-fold towers `D_l`, l = 2..6. -/
def tower2 : Fin 12 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 154, 168]
def tower3 : Fin 12 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 9, 140, 154, 168]
def tower4 : Fin 12 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 126, 140, 154, 168]
def tower5 : Fin 12 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 112, 126, 140, 154, 168]
def tower6 : Fin 12 → ℤ := ![1, 2, 3, 4, 5, 6, 98, 112, 126, 140, 154, 168]

/-- Tower l=2: margin ≥ 14/155 at t = 14/155 (the law holds here). -/
theorem tower2_margin :
    ∀ i, ∀ m : ℤ, (14 : ℝ) / 155 ≤ |(tower2 i : ℝ) * ((14 : ℝ) / 155) - m| :=
  rational_point_margin tower2 14 155 14 (by norm_num) (by decide)

/-- Tower l=3: margin ≥ 14/141 at t = 14/141. -/
theorem tower3_margin :
    ∀ i, ∀ m : ℤ, (14 : ℝ) / 141 ≤ |(tower3 i : ℝ) * ((14 : ℝ) / 141) - m| :=
  rational_point_margin tower3 14 141 14 (by norm_num) (by decide)

/-- Tower l=4: margin ≥ 17/155 at t = 17/155 — the CORRECTED witness (the
`14/127` law value is false here; the binding pair migrates to the 154+1 grid). -/
theorem tower4_margin :
    ∀ i, ∀ m : ℤ, (17 : ℝ) / 155 ≤ |(tower4 i : ℝ) * ((17 : ℝ) / 155) - m| :=
  rational_point_margin tower4 17 155 17 (by norm_num) (by decide)

/-- Tower l=5: margin ≥ 19/155 at t = 19/155 — the CORRECTED witness. -/
theorem tower5_margin :
    ∀ i, ∀ m : ℤ, (19 : ℝ) / 155 ≤ |(tower5 i : ℝ) * ((19 : ℝ) / 155) - m| :=
  rational_point_margin tower5 19 155 19 (by norm_num) (by decide)

/-- Tower l=6: margin ≥ 14/99 at t = 14/99 (the law resumes). -/
theorem tower6_margin :
    ∀ i, ∀ m : ℤ, (14 : ℝ) / 99 ≤ |(tower6 i : ℝ) * ((14 : ℝ) / 99) - m| :=
  rational_point_margin tower6 14 99 14 (by norm_num) (by decide)

/-! ## The 2/25 floor corollaries (what HYP-4212's domination assembly consumes) -/

theorem tower2_floor :
    ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(tower2 i : ℝ) * ((14 : ℝ) / 155) - m| :=
  fun i m => le_trans (by norm_num) (tower2_margin i m)

theorem tower3_floor :
    ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(tower3 i : ℝ) * ((14 : ℝ) / 141) - m| :=
  fun i m => le_trans (by norm_num) (tower3_margin i m)

theorem tower4_floor :
    ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(tower4 i : ℝ) * ((17 : ℝ) / 155) - m| :=
  fun i m => le_trans (by norm_num) (tower4_margin i m)

theorem tower5_floor :
    ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(tower5 i : ℝ) * ((19 : ℝ) / 155) - m| :=
  fun i m => le_trans (by norm_num) (tower5_margin i m)

theorem tower6_floor :
    ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(tower6 i : ℝ) * ((14 : ℝ) / 99) - m| :=
  fun i m => le_trans (by norm_num) (tower6_margin i m)

#print axioms tower4_margin
#print axioms tower5_floor

end MultiFoldRows
end LonelyRunner
