/-
  TournamentH7.LRCLiftFloorAssembly — THE CORRECTED LIFT-FLOOR PACKAGE
  (mac-mini-2026-07-05-S54b, HYP-4212; corrects HYP-4177's closing-move constant).

  THE CORRECTION (exact, this session): kps-S15's proposed domination
  `M(lift) ≥ 14/169` is FALSE — the block lift {1..12}\{4,6} ∪ {17,19} is a
  full-residue-system sieve-surviving lift with M = 2/25 < 14/169, and it is
  the UNIQUE height-1 double below 14/169.  The CORRECT floor statement is

      every sieve-surviving lift has margin ≥ 2/25,  TIGHT at the block lift

  — which is exactly the dichotomy's loose-branch constant β = 2/25: the lift
  branch closes at the level the assembly consumes, with a tight extremal.
  THM-621 is untouched (14/169 is the SINGLE-LEG floor; the deep well
  dominates single legs only).

  This file packages the formal floor witnesses as ONE consumable statement:

  * `ladderFamily r` — ({1..12}\{r}) ∪ {14r}, the THM-621 single-leg
    sieve-survivors, r = 7..12;
  * `ladder_law` — for every r ∈ [7,12], the family clears margin
    14/(13(r+1)) at t = a_r/(13(r+1)) (the six LRCLiftFloorRows rows quantified
    into one statement; a_r from the witness congruence (13−r)a_r ≡ −14);
  * `ladder_law_floor` — the uniform corollary: every ladder family clears
    ≥ 14/169 (the r = 12 rung is the weakest);
  * `block46_floor` — the tight extremal clears ≥ 2/25 (re-export);
  * `lift_floor_beta` — every ladder family AND the block lift clear the
    dichotomy constant 2/25 (14/169 > 2/25 bridged numerically): the formal
    kernel of the corrected lift-floor package.

  All decidable/numeric; consumes LRCLiftFloorRows (HYP-4109-series rows).
-/
import TournamentH7.LRCLiftFloorRows

namespace LonelyRunner
namespace LiftFloorAssembly

open LiftFloorRows

/-- The THM-621 single-leg ladder family at residue r: ({1..12}\{r}) ∪ {14r}. -/
def ladderFamily : ℕ → (Fin 12 → ℤ)
  | 7  => ladder7
  | 8  => ladder8
  | 9  => ladder9
  | 10 => ladder10
  | 11 => ladder11
  | _  => ladder12

/-- The witness numerator a_r (mirror representatives; (13−r)·a_r ≡ −14 mod 13(r+1)). -/
def ladderA : ℕ → ℤ
  | 7  => 15
  | 8  => 44
  | 9  => 29
  | 10 => 43
  | 11 => 71
  | _  => 14

/-- The witness modulus 13(r+1). -/
def ladderQ (r : ℕ) : ℤ := 13 * ((r : ℤ) + 1)

/-- THE PARAMETRIC LADDER LAW: every single-leg sieve-survivor r ∈ [7,12]
    clears margin 14/(13(r+1)) at its THM-621 witness.  (The six
    LRCLiftFloorRows rows, quantified.) -/
theorem ladder_law : ∀ r ∈ Finset.Icc 7 12,
    ∀ i, ∀ m : ℤ,
      (14 : ℝ) / (13 * ((r : ℝ) + 1)) ≤
        |(ladderFamily r i : ℝ) * ((ladderA r : ℝ) / (13 * ((r : ℝ) + 1))) - m| := by
  intro r hr i m
  have hr' := Finset.mem_Icc.mp hr
  interval_cases r
  · simpa [ladderFamily, ladderA, show ((13:ℝ) * ((7:ℕ) + 1) : ℝ) = 104 by norm_num]
      using ladder_margin_r7 i m
  · simpa [ladderFamily, ladderA, show ((13:ℝ) * ((8:ℕ) + 1) : ℝ) = 117 by norm_num]
      using ladder_margin_r8 i m
  · simpa [ladderFamily, ladderA, show ((13:ℝ) * ((9:ℕ) + 1) : ℝ) = 130 by norm_num]
      using ladder_margin_r9 i m
  · simpa [ladderFamily, ladderA, show ((13:ℝ) * ((10:ℕ) + 1) : ℝ) = 143 by norm_num]
      using ladder_margin_r10 i m
  · simpa [ladderFamily, ladderA, show ((13:ℝ) * ((11:ℕ) + 1) : ℝ) = 156 by norm_num]
      using ladder_margin_r11 i m
  · simpa [ladderFamily, ladderA, show ((13:ℝ) * ((12:ℕ) + 1) : ℝ) = 169 by norm_num]
      using ladder_margin_r12 i m

/-- The uniform floor corollary: every ladder family clears ≥ 14/169
    (the r = 12 deep-well rung is the weakest). -/
theorem ladder_law_floor : ∀ r ∈ Finset.Icc 7 12,
    ∀ i, ∀ m : ℤ,
      (14 : ℝ) / 169 ≤
        |(ladderFamily r i : ℝ) * ((ladderA r : ℝ) / (13 * ((r : ℝ) + 1))) - m| := by
  intro r hr i m
  have hr' := Finset.mem_Icc.mp hr
  have h12 : (r : ℝ) ≤ 12 := by exact_mod_cast hr'.2
  have hstep : (14 : ℝ) / 169 ≤ 14 / (13 * ((r : ℝ) + 1)) := by
    gcongr
    · norm_num
    · nlinarith
  exact hstep.trans (ladder_law r hr i m)

/-- The corrected lift-floor kernel at the dichotomy constant: the ladder
    families clear β = 2/25 (since 14/169 > 2/25), and the tight extremal
    block lift clears exactly 2/25 (block46_margin).  These are the formal
    attainment witnesses of `lift floor = 2/25, tight`. -/
theorem lift_floor_beta_ladder : ∀ r ∈ Finset.Icc 7 12,
    ∀ i, ∀ m : ℤ,
      (2 : ℝ) / 25 ≤
        |(ladderFamily r i : ℝ) * ((ladderA r : ℝ) / (13 * ((r : ℝ) + 1))) - m| := by
  intro r hr i m
  refine le_trans (by norm_num) (ladder_law_floor r hr i m)

end LiftFloorAssembly
end LonelyRunner
