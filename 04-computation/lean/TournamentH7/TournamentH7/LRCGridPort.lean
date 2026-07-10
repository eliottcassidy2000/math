/-
  TournamentH7.LRCGridPort  (boxeph-2026-07-09-S8, HYP-5818)

  The REVERSE-TRIANGLE payoff layer of the aliasing program: converting
  `thm665_full`'s deviation bound into grid-point EXISTENCE certificates —
  the abstract Lean core of THM-666 (the clamped grid-port).

  The move: `thm665_full` bounds the DEVIATION `‖E_grid[W] − ∫W‖ ≤
  TV(W′)/(12V²)`; the proof needs EXISTENCE of good ruler points.  The
  reverse triangle inequality `Re(E_grid[W]) ≥ ∫W − ‖E_grid[W] − ∫W‖`
  turns the deviation bound into a mean floor; pointwise domination
  `Re W(j/V) ≤ 1[P j]` turns the mean floor into a count floor; and a
  positive count is a witness.  This is the metric analogue of Rédei's
  parity move: existence from `count ≥ mean·V − deviation·V > 0`
  instead of `count ≡ 1 (mod 2)`.

  Instantiating `P j := "maxgap at j/V ≥ θ″"` with the clamp's
  `(α, a, w₀)` Bernoulli-basis data (the clamp is continuous piecewise
  linear; the cell engine emits the data) gives THM-666:
  `gridfrac{j : maxgap(j/V) ≥ θ″} ≥ μ(θ′) − TV(C′)/(12V²)`; and
  `exists_good_of_mean_pos` is the a-priori good-period existence on the
  V-ruler for `V > V₀ = √(TV/(12·∫C))`.

  Kernel-pure target: no `sorry`, no `native_decide` (audited below).
-/
import Mathlib
import TournamentH7.LRCPLFourier

set_option maxHeartbeats 1000000

open Complex Finset

namespace LonelyRunner
namespace Aliasing

noncomputable section

/-- **The reverse-triangle mean floor.**  For real mean data `w₀ : ℝ`, the
real part of the grid average is at least `w₀ − TV(W′)/(12V²)`:
`Re z ≥ w₀ − ‖z − w₀‖`, composed with `thm665_full`. -/
theorem grid_mean_re_ge (k : ℕ) (α : Fin k → ℝ) (a : Fin k → ℝ) (w0 : ℝ)
    (V : ℕ) (hV : 0 < V) :
    w0 - (2 * ∑ i, |α i|) / (12 * (V : ℝ) ^ 2)
      ≤ ((V : ℂ)⁻¹ * ∑ j ∈ range V, plW k α a (w0 : ℂ) ((j : ℝ) / (V : ℝ))).re := by
  have h := thm665_full k α a (w0 : ℂ) V hV
  set z : ℂ := (V : ℂ)⁻¹ * ∑ j ∈ range V, plW k α a (w0 : ℂ) ((j : ℝ) / (V : ℝ)) with hz
  have h1 : |(z - (w0 : ℂ)).re| ≤ ‖z - (w0 : ℂ)‖ := Complex.abs_re_le_norm _
  have h2 : (z - (w0 : ℂ)).re = z.re - w0 := by
    simp
  rw [h2] at h1
  have h3 : |z.re - w0| ≤ (2 * ∑ i, |α i|) / (12 * (V : ℝ) ^ 2) := le_trans h1 h
  have h4 := (abs_le.mp h3).1
  linarith

/-- **The grid-count port (THM-666's abstract core).**  If `Re W` is
dominated on the grid by the indicator of `P`, the number of grid points
satisfying `P` is at least `V·(w₀ − TV(W′)/(12V²))`. -/
theorem good_count_ge (k : ℕ) (α : Fin k → ℝ) (a : Fin k → ℝ) (w0 : ℝ)
    (V : ℕ) (hV : 0 < V) (P : ℕ → Prop) [DecidablePred P]
    (hdom : ∀ j ∈ range V, (plW k α a (w0 : ℂ) ((j : ℝ) / (V : ℝ))).re
      ≤ if P j then (1 : ℝ) else 0) :
    (w0 - (2 * ∑ i, |α i|) / (12 * (V : ℝ) ^ 2)) * (V : ℝ)
      ≤ (((range V).filter P).card : ℝ) := by
  have hVR : (0 : ℝ) < (V : ℝ) := Nat.cast_pos.mpr hV
  have hmean := grid_mean_re_ge k α a w0 V hV
  -- Re of the average = (1/V) · Σ of the real parts
  have hre : ((V : ℂ)⁻¹ * ∑ j ∈ range V, plW k α a (w0 : ℂ) ((j : ℝ) / (V : ℝ))).re
      = ((V : ℝ))⁻¹ * ∑ j ∈ range V, (plW k α a (w0 : ℂ) ((j : ℝ) / (V : ℝ))).re := by
    rw [show ((V : ℂ))⁻¹ = ((((V : ℝ))⁻¹ : ℝ) : ℂ) by push_cast; ring]
    rw [Complex.mul_re, Complex.re_sum]
    simp
  rw [hre] at hmean
  -- the summed real parts are at most the count
  have hcount : ∑ j ∈ range V, (plW k α a (w0 : ℂ) ((j : ℝ) / (V : ℝ))).re
      ≤ (((range V).filter P).card : ℝ) := by
    calc ∑ j ∈ range V, (plW k α a (w0 : ℂ) ((j : ℝ) / (V : ℝ))).re
        ≤ ∑ j ∈ range V, (if P j then (1 : ℝ) else 0) := Finset.sum_le_sum hdom
      _ = (((range V).filter P).card : ℝ) := by
          rw [Finset.sum_boole]
  -- multiply the mean floor by V and chain
  have h2 : (w0 - (2 * ∑ i, |α i|) / (12 * (V : ℝ) ^ 2)) * (V : ℝ)
      ≤ ∑ j ∈ range V, (plW k α a (w0 : ℂ) ((j : ℝ) / (V : ℝ))).re := by
    have h1 := mul_le_mul_of_nonneg_right hmean hVR.le
    rwa [inv_mul_eq_div, div_mul_cancel₀ _ hVR.ne'] at h1
  linarith

/-- **A-priori existence on the ruler.**  If the mean beats the aliasing
deviation — `TV(W′)/(12V²) < w₀`, i.e. `V > V₀ = √(TV/(12·w₀))` — then some
grid point satisfies `P`.  The metric counterpart of Rédei's parity
argument: a witness from `count > 0` by mean-minus-deviation. -/
theorem exists_good_of_mean_pos (k : ℕ) (α : Fin k → ℝ) (a : Fin k → ℝ)
    (w0 : ℝ) (V : ℕ) (hV : 0 < V) (P : ℕ → Prop) [DecidablePred P]
    (hdom : ∀ j ∈ range V, (plW k α a (w0 : ℂ) ((j : ℝ) / (V : ℝ))).re
      ≤ if P j then (1 : ℝ) else 0)
    (hpos : (2 * ∑ i, |α i|) / (12 * (V : ℝ) ^ 2) < w0) :
    ∃ j, j < V ∧ P j := by
  have hVR : (0 : ℝ) < (V : ℝ) := Nat.cast_pos.mpr hV
  have h := good_count_ge k α a w0 V hV P hdom
  have hcard : (0 : ℝ) < (((range V).filter P).card : ℝ) := by
    have hlow : (0 : ℝ) < (w0 - (2 * ∑ i, |α i|) / (12 * (V : ℝ) ^ 2)) * (V : ℝ) :=
      mul_pos (by linarith) hVR
    linarith
  have hne : ((range V).filter P).Nonempty := by
    rw [← Finset.card_pos]
    exact_mod_cast hcard
  obtain ⟨j, hj⟩ := hne
  rw [Finset.mem_filter, Finset.mem_range] at hj
  exact ⟨j, hj.1, hj.2⟩

end

end Aliasing
end LonelyRunner

-- kernel-purity audit
#print axioms LonelyRunner.Aliasing.grid_mean_re_ge
#print axioms LonelyRunner.Aliasing.good_count_ge
#print axioms LonelyRunner.Aliasing.exists_good_of_mean_pos
