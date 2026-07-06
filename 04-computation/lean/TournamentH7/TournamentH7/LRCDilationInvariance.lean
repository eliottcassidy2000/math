/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S110)
-/
import TournamentH7.LRCWitnessAttainment

/-!
# Dilation invariance of the loneliness minimax (HYP-4426)

The loneliness constant `M(v) = ⨆_t margin v t` depends only on the *ratios* of the speeds:
`M(c·v) = M(v)` for any nonzero integer dilation `c`.  Scaling every speed by `c` merely
reparametrizes time by `c` (`margin_const_mul`), and `t ↦ c·t` is a surjection of `ℝ`, so the
supremum is unchanged.

This is the formal **"up to dilation" primitive** the project has used informally throughout.
Its role in the (G) endgame is the *normalization* step of the indirect (compactness) route to
the height bound: since `M` is dilation-invariant, a putative sequence of integer gap-members
at heights `H_k → ∞` may be rescaled to a bounded window `w_i = v_i / H_k ∈ [1, ρ]`, whose
limit is a bounded-ratio real family with `M ∈ [1/13, 2/25]`.  The height coordinate is
*gauge* — only the ratio pattern (equivalently the height-independent difference-core; S110)
carries loneliness.  So "bound the height" is really "control the ratio pattern," and dilation
invariance is what makes that reduction legitimate.
-/

namespace LonelyRunner
namespace Dilation

open TournamentH7.LRCWitness

variable {k : ℕ} [Nonempty (Fin k)]

/-- **Pointwise dilation.**  Scaling every speed by `c` reparametrizes time by `c`:
`margin (c·v) t = margin v (c·t)`. -/
lemma margin_const_mul (c : ℤ) (v : Fin k → ℤ) (t : ℝ) :
    margin (fun i => c * v i) t = margin v ((c : ℝ) * t) := by
  unfold margin
  have hfun : (fun i => distZ ((↑(c * v i) : ℝ) * t))
      = (fun i => distZ ((v i : ℝ) * ((c : ℝ) * t))) := by
    funext i; congr 1; push_cast; ring
  rw [hfun]

/-- **Dilation invariance of the minimax constant.**  `M(c·v) = M(v)` for `c ≠ 0`: the
loneliness constant `⨆_t margin` depends only on the speed *ratios*.  This is the formal
"up to dilation" primitive — the normalization that lets the descent/compactness argument
assume a gap-member family has bounded height (the height is a gauge coordinate). -/
theorem iSup_margin_const_mul (c : ℤ) (hc : c ≠ 0) (v : Fin k → ℤ) :
    (⨆ t : ℝ, margin (fun i => c * v i) t) = ⨆ t : ℝ, margin v t := by
  have hcR : (c : ℝ) ≠ 0 := by exact_mod_cast hc
  have hsurj : Function.Surjective (fun t : ℝ => (c : ℝ) * t) := by
    intro s; exact ⟨s / c, by field_simp⟩
  simp_rw [margin_const_mul c v]
  exact hsurj.iSup_comp (margin v)

#print axioms iSup_margin_const_mul

end Dilation
end LonelyRunner
