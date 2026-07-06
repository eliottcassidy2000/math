/-
  TournamentH7.LRCSingleLift — THE SINGLE-LIFT RIGIDITY: every single-13-lift of
  the AP {1,…,12} is loose at 2/25 (kind-pasteur-2026-07-06-S28).

  The density floor (the genuinely-open residual of (G)) is the rigidity of the
  AP's roots-of-unity fiber: the AP is the UNIQUE covering family at 2/25, and
  every perturbation jumps M out of the gap.  Its BASE CASE is the single-13-lift:
  replace one AP runner `j+1` by `j+1+13` (still a full residue system mod 13 —
  the "residue-pinned" side of the S23 split).  These are the simplest non-AP
  near-AP families; mac-mini's density-floor probe (S17) finds the safe-minimising
  perturbations exactly here (M = 1/12, 2/23).

  This file formalises that ALL 12 single-lifts are LOOSE — each carries an
  explicit `2/25`-margin point — via the `rational_point_margin` atom with the
  computed witnesses (lrc single-lift scan):

    runner 1→14 : t=1/16   2→15 : 9/19   3→16 : 6/17   4→17 : 5/19
    5→18 : 4/19  6→19 : 4/23 (=2/23, the tightest)   7→20 : 1/7   8→21 : 1/8
    9→22 : 1/9   10→23 : 1/10   11→24 : 1/11   12→25 : 1/12

  Each margin exceeds `2/25`, so the single-lift fiber never enters the gap — the
  formal base case of the near-AP rigidity.  Kernel-pure; `decide` on the finite
  residue conditions, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCHarmonicGate

namespace LonelyRunner
namespace SingleLift

open HarmonicGate
open scoped Classical

/-- The AP `{1,…,12}` with runner `j` (0-indexed) lifted by `13`. -/
def apLift (j : Fin 12) (i : Fin 12) : ℤ :=
  ((i : ℕ) + 1 : ℤ) + (if i = j then 13 else 0)

/-- **Loose certificate from a rational witness**: if at `t = k/s` every runner's
residue lands in `[μ, s−μ]` and `2/25 ≤ μ/s`, the family is loose at `2/25`. -/
theorem loose_cert (v : Fin 12 → ℤ) (k s μ : ℤ) (hs : 0 < s)
    (hμ : (2 : ℝ) / 25 ≤ (μ : ℝ) / s)
    (hc : ∀ i : Fin 12, μ ≤ (v i * k) % s ∧ (v i * k) % s ≤ s - μ) :
    ∃ t : ℝ, ∀ i : Fin 12, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(v i : ℝ) * t - m| :=
  ⟨(k : ℝ) / s, fun i m => le_trans hμ (rational_point_margin v k s μ hs hc i m)⟩

/-- **THE SINGLE-LIFT RIGIDITY**: every single-13-lift of the AP is loose at
`2/25` — the formal base case of the density floor (the near-AP fiber rigidity). -/
theorem ap_single_lift_loose (j : Fin 12) :
    ∃ t : ℝ, ∀ i : Fin 12, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(apLift j i : ℝ) * t - m| := by
  fin_cases j
  · exact loose_cert _ 1 16 2 (by norm_num) (by norm_num) (by decide)
  · exact loose_cert _ 9 19 2 (by norm_num) (by norm_num) (by decide)
  · exact loose_cert _ 6 17 2 (by norm_num) (by norm_num) (by decide)
  · exact loose_cert _ 5 19 2 (by norm_num) (by norm_num) (by decide)
  · exact loose_cert _ 4 19 2 (by norm_num) (by norm_num) (by decide)
  · exact loose_cert _ 4 23 2 (by norm_num) (by norm_num) (by decide)
  · exact loose_cert _ 1 7 1 (by norm_num) (by norm_num) (by decide)
  · exact loose_cert _ 1 8 1 (by norm_num) (by norm_num) (by decide)
  · exact loose_cert _ 1 9 1 (by norm_num) (by norm_num) (by decide)
  · exact loose_cert _ 1 10 1 (by norm_num) (by norm_num) (by decide)
  · exact loose_cert _ 1 11 1 (by norm_num) (by norm_num) (by decide)
  · exact loose_cert _ 1 12 1 (by norm_num) (by norm_num) (by decide)

#print axioms ap_single_lift_loose

end SingleLift
end LonelyRunner
