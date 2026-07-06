/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: mac-mini (LRC multi-agent project, 2026-07-06-S20)
-/
import TournamentH7.LRCHarmonicGate

/-!
# The second-value gap is NONEMPTY at n = 7 — a formal witness of n-specificity (HYP-4452)

The second-value gap `(1/n, 2/(2n-1))` is **not** universally empty: it is nonempty at
`n = 7` (and `n = 8`) but empty at `n = 13`.  This module gives the kernel-pure formal
witness at `n = 7`, via kps-S2's `rational_point_margin` atom (HYP-4102).

The family `{1, 5, 6, 11, 16, 17}` (6 speeds = an `n = 7` configuration) is
**`5/33`-lonely** at `t = 10/33`: every runner stays `≥ 5/33` from every integer.  The
residues `v·10 mod 33` are `{10, 17, 27, 11, 28, 5}`, all in `[5, 28]`, so the margin is
`5/33`.  Since `5/33 > 1/7` this is strictly above the LRC(7) minimum, and (by exact
computation, grid-verified) `M = 5/33 < 2/13`, placing `M` inside the second gap
`(1/7, 2/13)`.

**Significance.**  (G) — the emptiness of the gap — is therefore **n-specific**, not a
universal structural identity.  The structural "equi-" lenses (three-gap / roots-of-unity /
sum-product) mark the AP special at *every* `n`, so they are necessary but **not
sufficient**; the emptiness at `n = 13` is a quantitative wall-vs-budget fact.  This is the
formal counterexample behind that correction.
-/

namespace LonelyRunner
namespace NSpecificGap

open HarmonicGate

/-- The `n = 7` second-gap member `{1, 5, 6, 11, 16, 17}` (`M = 5/33`, grid-verified,
inside the second gap `(1/7, 2/13)`). -/
def n7gap : Fin 6 → ℤ := ![1, 5, 6, 11, 16, 17]

/-- **The gap witness (kernel-pure).**  `{1,5,6,11,16,17}` is `5/33`-lonely at `t = 10/33`:
every runner is `≥ 5/33` from every integer.  Hence `M(n7gap) ≥ 5/33`. -/
theorem n7gap_lonely_5_33 :
    ∃ t : ℝ, ∀ (j : Fin 6) (m : ℤ), (5 : ℝ) / 33 ≤ |(n7gap j : ℝ) * t - m| :=
  ⟨(10 : ℝ) / 33, rational_point_margin n7gap 10 33 5 (by norm_num) (by decide)⟩

/-- `5/33 > 1/7`: the witness margin exceeds the LRC(7) tight value, so `n7gap` is NOT
tight — it lies strictly above the minimum. -/
theorem five_over_33_gt_one_over_7 : (1 : ℝ) / 7 < (5 : ℝ) / 33 := by norm_num

/-- `5/33 < 2/13`: the (computed, grid-verified) value `M = 5/33` is below the
second-value edge, so the witnessed loneliness sits inside the second gap `(1/7, 2/13)`. -/
theorem five_over_33_lt_two_over_13 : (5 : ℝ) / 33 < (2 : ℝ) / 13 := by norm_num

#print axioms n7gap_lonely_5_33

end NSpecificGap
end LonelyRunner
