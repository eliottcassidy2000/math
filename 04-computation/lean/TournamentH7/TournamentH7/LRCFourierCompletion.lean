/-
  TournamentH7.LRCFourierCompletion — LEM-022's Fourier completion, Stage A
  (death-star-2026-07-09-S13, HYP-5890).

  The last ℂ piece of the t2 hyperbola lemma (LEM-022, death-star-S9; heart in Lean S10,
  dyadic assembly S12).  STAGE A here: **the interval exponential-sum bound**

  > `‖Σ_{r < len} exp(2πi·h·(lo + r)/q)‖ ≤ q/(2·d)`

  under the sine witness `2·(d/q) ≤ |sin(π·h/q)|` (supplied at the call site from
  `d = cdist (h : ZMod q)` via Jordan's inequality and 1-periodicity).

  Mechanism: geometric sum with ratio `z = exp(2πi·h/q) ≠ 1`; `‖z^len − 1‖ ≤ 2`;
  `‖z − 1‖ = 2·|sin(π·h/q)|` (the factoring `exp(iθ) − 1 = exp(iθ/2)·(2i·sin(θ/2))`);
  Jordan's `Real.mul_le_sin` packaged as `two_mul_le_sin_pi_mul`.

  STAGE B (documented, next session): self-derived orthogonality
  `Σ_{x < q} exp(2πi·h·x/q) = q·1{q ∣ h}` (same geometric machinery), the completion identity
  `C_w = b²/q + (1/q)·Σ_{h≠0} B̂(h)·conj(B̂(w⁻¹h))`, and the assembly with
  `hyperbola_box_count` (S10) + `harmonic_ratio_sum_mul_le` (S12), noting `P(w⁻¹) = P(w)`.

  Kernel-pure target: no `sorry`, no `native_decide`.
-/
import Mathlib

namespace LonelyRunner
namespace FourierCompletion

open Complex Real

/-- **Jordan's inequality on the unit scale**: for `0 ≤ t ≤ 1/2`, `2·t ≤ sin (π·t)`.
(`Real.mul_le_sin` at `x = π·t`.) -/
theorem two_mul_le_sin_pi_mul {t : ℝ} (h0 : 0 ≤ t) (h2 : t ≤ 1 / 2) :
    2 * t ≤ Real.sin (π * t) := by
  have hx0 : 0 ≤ π * t := by positivity
  have hx2 : π * t ≤ π / 2 := by
    have hπ : (0 : ℝ) < π := Real.pi_pos
    nlinarith
  have := Real.mul_le_sin hx0 hx2
  calc 2 * t = 2 / π * (π * t) := by
        field_simp
      _ ≤ Real.sin (π * t) := this

/-- **`‖exp(iθ) − 1‖ = 2·|sin(θ/2)|`** — the standard factoring
`exp(iθ) − 1 = exp(iθ/2)·(2i·sin(θ/2))`. -/
theorem norm_exp_I_sub_one (θ : ℝ) :
    ‖Complex.exp ((θ : ℂ) * Complex.I) - 1‖ = 2 * |Real.sin (θ / 2)| := by
  set E : ℂ := Complex.exp ((θ / 2 : ℝ) * Complex.I) with hE
  set E' : ℂ := Complex.exp (-((θ / 2 : ℝ) : ℂ) * Complex.I) with hE'
  have hEE' : E * E' = 1 := by
    rw [hE, hE', ← Complex.exp_add]
    ring_nf
    exact Complex.exp_zero
  have h2Isin : (2 : ℂ) * Complex.I * ((Real.sin (θ / 2) : ℝ) : ℂ) = E - E' := by
    rw [Complex.ofReal_sin, Complex.sin, hE, hE']
    have hI2 : (Complex.I : ℂ) * Complex.I = -1 := Complex.I_mul_I
    push_cast
    field_simp
    ring_nf
    rw [Complex.I_sq]
    ring
  have hsplit : ((θ : ℂ)) * Complex.I = ((θ / 2 : ℝ) : ℂ) * Complex.I
      + ((θ / 2 : ℝ) : ℂ) * Complex.I := by
    push_cast
    ring
  have hfactor : Complex.exp ((θ : ℂ) * Complex.I) - 1
      = E * ((2 : ℂ) * Complex.I * ((Real.sin (θ / 2) : ℝ) : ℂ)) := by
    rw [h2Isin]
    calc Complex.exp ((θ : ℂ) * Complex.I) - 1
        = E * E - 1 := by rw [hsplit, Complex.exp_add, hE]
      _ = E * E - E * E' := by rw [hEE']
      _ = E * (E - E') := by ring
  rw [hfactor, norm_mul]
  have hnE : ‖E‖ = 1 := by
    rw [hE, Complex.norm_exp]
    simp
  rw [hnE, one_mul]
  rw [show ((2 : ℂ) * Complex.I * ((Real.sin (θ / 2) : ℝ) : ℂ))
      = ((2 * Real.sin (θ / 2) : ℝ) : ℂ) * Complex.I from by push_cast; ring]
  rw [norm_mul, Complex.norm_I, mul_one, Complex.norm_real, Real.norm_eq_abs, abs_mul]
  norm_num

/-- **Stage A: the interval exponential-sum bound.**  Under the sine witness
`hsin : 2·(d/q) ≤ |sin(π·h/q)|` with `d > 0`:

  `‖Σ_{r < len} exp(2π·(h·(lo + r)/q)·i)‖ ≤ q/(2·d)`. -/
theorem norm_expSum_le (q : ℕ) (hq : 0 < q) (h : ℤ) (lo : ℤ) (len : ℕ) (d : ℝ)
    (hd : 0 < d)
    (hsin : 2 * (d / q) ≤ |Real.sin (π * h / q)|) :
    ‖∑ r ∈ Finset.range len,
        Complex.exp ((2 * π * ((h : ℝ) * ((lo : ℝ) + r) / q) : ℝ) * Complex.I)‖
      ≤ (q : ℝ) / (2 * d) := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  set z : ℂ := Complex.exp ((2 * π * ((h : ℝ) / q) : ℝ) * Complex.I) with hz
  -- each term = (unit constant) · z^r
  have hterm : ∀ r : ℕ, Complex.exp ((2 * π * ((h : ℝ) * ((lo : ℝ) + r) / q) : ℝ) * Complex.I)
      = Complex.exp ((2 * π * ((h : ℝ) * (lo : ℝ) / q) : ℝ) * Complex.I) * z ^ r := by
    intro r
    rw [hz, ← Complex.exp_nat_mul, ← Complex.exp_add]
    congr 1
    push_cast
    field_simp
  rw [Finset.sum_congr rfl (fun r _ => hterm r), ← Finset.mul_sum, norm_mul]
  have hone : ‖Complex.exp ((2 * π * ((h : ℝ) * (lo : ℝ) / q) : ℝ) * Complex.I)‖ = 1 := by
    rw [Complex.norm_exp]
    simp
  rw [hone, one_mul]
  -- ‖z − 1‖ = 2|sin(π h/q)| ≥ 4d/q > 0, so z ≠ 1
  have hzsub : ‖z - 1‖ = 2 * |Real.sin (π * h / q)| := by
    rw [hz]
    have harg : (2 * π * ((h : ℝ) / q) : ℝ) = 2 * (π * h / q) := by ring
    rw [harg, norm_exp_I_sub_one,
      show (2 * (π * (h : ℝ) / q)) / 2 = π * (h : ℝ) / q from by ring]
  have hlow : 4 * d / q ≤ ‖z - 1‖ := by
    rw [hzsub]
    calc 4 * d / q = 2 * (2 * (d / q)) := by ring
      _ ≤ 2 * |Real.sin (π * h / q)| := by linarith
  have hzsub_pos : 0 < ‖z - 1‖ := lt_of_lt_of_le (by positivity) hlow
  have hzne : z ≠ 1 := by
    intro hcon
    rw [hcon] at hzsub_pos
    simp at hzsub_pos
  rw [geom_sum_eq hzne, norm_div]
  have hnum : ‖z ^ len - 1‖ ≤ 2 := by
    calc ‖z ^ len - 1‖ ≤ ‖z ^ len‖ + ‖(1 : ℂ)‖ := norm_sub_le _ _
      _ = 1 + 1 := by
          rw [norm_pow, hz, Complex.norm_exp]
          simp
      _ = 2 := by norm_num
  calc ‖z ^ len - 1‖ / ‖z - 1‖
      ≤ 2 / (4 * d / q) := div_le_div₀ (by norm_num) hnum (by positivity) hlow
    _ = (q : ℝ) / (2 * d) := by
        field_simp
        ring

/-! ## Axiom audit -/
#print axioms two_mul_le_sin_pi_mul
#print axioms norm_exp_I_sub_one
#print axioms norm_expSum_le

end FourierCompletion
end LonelyRunner
