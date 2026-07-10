/-
  TournamentH7.LRCP2Eval  (boxeph-2026-07-09-S10)

  Concrete evaluation of the periodized Bernoulli-2 basis function — the
  single enabler for every decide-shaped consumer of the aliasing surface
  (per-(E,V) window certificates, the THM-666 clamp instantiation, the D_m
  product-clamp checks): with `P₂` reduced to the explicit quadratic
  `fract(x)² − fract(x) + 1/6`, the real part of `plW` at rational grid
  points becomes a finite explicit rational expression, and the domination
  hypothesis of `exists_good_of_mean_pos` becomes `norm_num`-checkable
  per grid point from cell-engine `(α, a, w₀)` data.

    * `bernoulliFun_two` — `B₂(x) = x² − x + 1/6` (the general-`x`
      evaluation, by the ZetaValues example's own tactic recipe).
    * `P2_eval` — `P₂(x) = fract(x)² − fract(x) + 1/6`.
    * `P2_grid` — at in-range grid points `j/V` the fract is the identity:
      `P₂(j/V) = (j/V)² − j/V + 1/6`.

  Kernel-pure target: no `sorry`, no `native_decide` (audited below).
-/
import Mathlib
import TournamentH7.LRCPLFourier

open Complex Finset

namespace LonelyRunner
namespace Aliasing

noncomputable section

/-- **The Bernoulli-2 polynomial, evaluated:** `B₂(x) = x² − x + 1/6`. -/
theorem bernoulliFun_two (x : ℝ) : bernoulliFun 2 x = x ^ 2 - x + 1 / 6 := by
  rw [bernoulliFun]
  simp_rw [Polynomial.bernoulli, Finset.sum_range_succ, Finset.sum_range_zero,
    Polynomial.map_add, Polynomial.map_zero, Polynomial.map_monomial,
    Polynomial.eval_add, Polynomial.eval_zero, Polynomial.eval_monomial]
  rw [_root_.bernoulli_one, bernoulli_eq_bernoulli'_of_ne_one zero_ne_one, bernoulli'_zero,
    bernoulli_eq_bernoulli'_of_ne_one (by decide : (2 : ℕ) ≠ 1), bernoulli'_two]
  norm_num
  ring

/-- **P₂, evaluated:** `P₂(x) = fract(x)² − fract(x) + 1/6` (as a complex
number via the real coercion). -/
theorem P2_eval (x : ℝ) :
    P2 x = ((Int.fract x ^ 2 - Int.fract x + 1 / 6 : ℝ) : ℂ) := by
  rw [P2, bernoulliFun_two]

/-- **P₂ at in-range grid points:** for `j < V`,
`P₂(j/V) = (j/V)² − j/V + 1/6` — the fract is the identity on `[0,1)`,
so the value is an explicit rational in `j` and `V`. -/
theorem P2_grid (V : ℕ) (hV : 0 < V) (j : ℕ) (hj : j < V) :
    P2 ((j : ℝ) / (V : ℝ))
      = ((((j : ℝ) / (V : ℝ)) ^ 2 - (j : ℝ) / (V : ℝ) + 1 / 6 : ℝ) : ℂ) := by
  rw [P2_eval]
  have hfr : Int.fract ((j : ℝ) / (V : ℝ)) = (j : ℝ) / (V : ℝ) := by
    apply Int.fract_eq_self.mpr
    have hVR : (0 : ℝ) < (V : ℝ) := Nat.cast_pos.mpr hV
    constructor
    · positivity
    · rw [div_lt_one hVR]
      exact_mod_cast hj
  rw [hfr]

end

end Aliasing
end LonelyRunner

-- kernel-purity audit
#print axioms LonelyRunner.Aliasing.bernoulliFun_two
#print axioms LonelyRunner.Aliasing.P2_eval
#print axioms LonelyRunner.Aliasing.P2_grid
