/-
  TournamentH7.LRCLiftFloorRows — THE MULTI-LIFT FLOOR CERTIFICATE ROWS
  (mac-mini-2026-07-05-S53, HYP-4109; certificates computed in S52, HYP-4103).

  Seven decidable instances of kps-S2's `rational_point_margin` atom (HYP-4102),
  formalizing the lift-stratum floor witnesses of the hdich/TightLooseDichotomy
  analysis:

  * `block46_margin` — the DOUBLE BLOCK LIFT `{1..12}\{4,6} ∪ {17,19}` has margin
    ≥ 2/25 at `t = 6/25`.  This is THE multi-lift floor family (S52: exact
    M = 2/25, attained): the lift stratum meets the global n=12 second value.
    Binders: runner 8 at −2, runner 17 at +2 (mod 25).

  * `ladder_margin_r7 … _r12` — the 14r LADDER `W_r = {1..12}\{r} ∪ {14r}`,
    r = 7..12, margin ≥ 14/(13(r+1)) at `t = a_r/(13(r+1))` where
    `(13−r)·a_r ≡ −14 (mod 13(r+1))` (the S52 witness congruence; binding pair
    (13−r, 14r) at ∓14).  Row r = 12 is the n=13 deep well `{1..11,168}` at
    margin 14/169 — the single-lift floor.  (S52 exact values: these margins
    are EXACT, i.e. M(W_r) = 14/(13(r+1)); this file certifies the ≥ direction,
    which is what the dichotomy's loose branch consumes.)

  All hypotheses discharged by `decide` on 12 integer mod conditions each.
  a_r values (mirror representatives): r=7: 15, r=8: 44, r=9: 29, r=10: 43,
  r=11: 71, r=12: 14.
-/
import TournamentH7.LRCHarmonicGate

namespace LonelyRunner
namespace LiftFloorRows

open HarmonicGate

/-- The double block lift `{1..12}\{4,6} ∪ {17,19}`: lift 4→17, 6→19 (k=1 both). -/
def block46 : Fin 12 → ℤ := ![1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 17, 19]

/-- The 14r ladder families `{1..12}\{r} ∪ {14r}`, r = 7..12. -/
def ladder7  : Fin 12 → ℤ := ![1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 98]
def ladder8  : Fin 12 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 112]
def ladder9  : Fin 12 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 126]
def ladder10 : Fin 12 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 140]
def ladder11 : Fin 12 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 154]
def ladder12 : Fin 12 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 168]

/-- The multi-lift floor row: margin ≥ 2/25 at t = 6/25 for the {4,6}-block lift.
    (S52: M = 2/25 exactly — the lift stratum attains the n=12 second value.) -/
theorem block46_margin :
    ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(block46 i : ℝ) * ((6 : ℝ) / 25) - m| :=
  rational_point_margin block46 6 25 2 (by norm_num)
    (by decide)

/-- Ladder row r=7: margin ≥ 14/104 = 7/52 at t = 15/104 for {1..6,8..12,98}. -/
theorem ladder_margin_r7 :
    ∀ i, ∀ m : ℤ, (14 : ℝ) / 104 ≤ |(ladder7 i : ℝ) * ((15 : ℝ) / 104) - m| :=
  rational_point_margin ladder7 15 104 14 (by norm_num) (by decide)

/-- Ladder row r=8: margin ≥ 14/117 at t = 44/117. -/
theorem ladder_margin_r8 :
    ∀ i, ∀ m : ℤ, (14 : ℝ) / 117 ≤ |(ladder8 i : ℝ) * ((44 : ℝ) / 117) - m| :=
  rational_point_margin ladder8 44 117 14 (by norm_num) (by decide)

/-- Ladder row r=9: margin ≥ 14/130 = 7/65 at t = 29/130. -/
theorem ladder_margin_r9 :
    ∀ i, ∀ m : ℤ, (14 : ℝ) / 130 ≤ |(ladder9 i : ℝ) * ((29 : ℝ) / 130) - m| :=
  rational_point_margin ladder9 29 130 14 (by norm_num) (by decide)

/-- Ladder row r=10: margin ≥ 14/143 at t = 43/143. -/
theorem ladder_margin_r10 :
    ∀ i, ∀ m : ℤ, (14 : ℝ) / 143 ≤ |(ladder10 i : ℝ) * ((43 : ℝ) / 143) - m| :=
  rational_point_margin ladder10 43 143 14 (by norm_num) (by decide)

/-- Ladder row r=11: margin ≥ 14/156 = 7/78 at t = 71/156. -/
theorem ladder_margin_r11 :
    ∀ i, ∀ m : ℤ, (14 : ℝ) / 156 ≤ |(ladder11 i : ℝ) * ((71 : ℝ) / 156) - m| :=
  rational_point_margin ladder11 71 156 14 (by norm_num) (by decide)

/-- Ladder row r=12 — the n=13 DEEP WELL {1..11,168}: margin ≥ 14/169 at
    t = 14/169.  The single-lift floor (S51/S77/MISTAKE-104). -/
theorem ladder_margin_r12 :
    ∀ i, ∀ m : ℤ, (14 : ℝ) / 169 ≤ |(ladder12 i : ℝ) * ((14 : ℝ) / 169) - m| :=
  rational_point_margin ladder12 14 169 14 (by norm_num) (by decide)

end LiftFloorRows
end LonelyRunner
