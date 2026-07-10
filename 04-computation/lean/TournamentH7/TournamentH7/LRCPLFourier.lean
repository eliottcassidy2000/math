/-
  TournamentH7.LRCPLFourier  (boxeph-2026-07-09-S7, HYP-5778 completion)

  The piecewise-linear Fourier input for THM-665 — via the PERIODIZED
  BERNOULLI BASIS, with NO integration by parts.

  Key identity: every continuous 1-periodic piecewise-linear W with
  W′-jumps j_i at breakpoints a_i is EXACTLY
      W x = w₀ + Σ_i α_i · P₂(x − a_i),   α_i = −j_i/2,  Σ_i α_i = 0,
  where P₂ = periodized Bernoulli-2 (B₂ ∘ fract).  Mathlib's ZetaValues
  supplies P₂'s pointwise Fourier expansion
  (`hasSum_one_div_pow_mul_fourier_mul_bernoulliFun`, k = 2):
      B₂(fract x) = (1/2π²)·Σ'_n (1/n²)·e(2πi n x),
  so the coefficients of W are
      c_n(W) = Σ_i α_i e(−2πi n a_i)/(2π²n²)   (n ≠ 0),   c_0 = w₀,
  giving ‖c_n‖ ≤ (Σ|α_i|)/(2π²n²) = TV(W′)/(4π²n²)  [TV(W′) = 2Σ|α_i|]
  — precisely the hypotheses of `LRCAliasingBound.aliasing_tail_bound`.

  Main result `thm665_full`: for W in this presentation,
      ‖(1/V)·Σ_{j<V} W(j/V) − w₀‖ ≤ (2·Σ|α_i|)/(12·V²),
  i.e. |E_grid[W](V) − ∫W| ≤ TV(W′)/(12V²) — THM-665, complete.

  Kernel-pure target: no `sorry`, no `native_decide` (audited below).
-/
import Mathlib
import TournamentH7.LRCAliasingBound

set_option maxHeartbeats 1000000

open Complex Finset

namespace LonelyRunner
namespace Aliasing

noncomputable section

/-- `e(2πi n x)` only sees `x` mod 1. -/
theorem exp_two_pi_int_fract (n : ℤ) (x : ℝ) :
    Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * (Int.fract x : ℝ))
      = Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * x) := by
  have hfr : ((Int.fract x : ℝ) : ℂ) = (x : ℂ) - ((⌊x⌋ : ℤ) : ℂ) := by
    rw [Int.fract]
    push_cast
    ring
  rw [hfr, show 2 * (Real.pi : ℂ) * Complex.I * n * ((x : ℂ) - ((⌊x⌋ : ℤ) : ℂ))
      = 2 * (Real.pi : ℂ) * Complex.I * n * x
        - ((n * ⌊x⌋ : ℤ) : ℂ) * (2 * (Real.pi : ℂ) * Complex.I) by push_cast; ring,
    Complex.exp_sub, Complex.exp_int_mul_two_pi_mul_I, div_one]

/-- The periodized Bernoulli-2 function, as a map `ℝ → ℂ`. -/
def P2 (x : ℝ) : ℂ := ((bernoulliFun 2 (Int.fract x) : ℝ) : ℂ)

/-- **P₂'s pointwise Fourier expansion** (from Mathlib's
`hasSum_one_div_pow_mul_fourier_mul_bernoulliFun` at `k = 2`):
`HasSum (fun n => e(2πi n x)/(2π²n²)) (P₂ x)` for every real `x`. -/
theorem hasSum_P2 (x : ℝ) :
    HasSum (fun n : ℤ => 1 / (2 * (Real.pi : ℂ) ^ 2 * (n : ℂ) ^ 2)
      * Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * x)) (P2 x) := by
  have hx : Int.fract x ∈ Set.Icc (0 : ℝ) 1 :=
    ⟨Int.fract_nonneg x, (Int.fract_lt_one x).le⟩
  have h := hasSum_one_div_pow_mul_fourier_mul_bernoulliFun (k := 2) le_rfl hx
  have h2 := h.mul_left (1 / (2 * (Real.pi : ℂ) ^ 2))
  have hfun : (fun n : ℤ => 1 / (2 * (Real.pi : ℂ) ^ 2)
      * (1 / (n : ℂ) ^ 2 * fourier n ((Int.fract x : ℝ) : AddCircle (1 : ℝ))))
      = fun n : ℤ => 1 / (2 * (Real.pi : ℂ) ^ 2 * (n : ℂ) ^ 2)
      * Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * x) := by
    funext n
    rw [fourier_coe_apply, Complex.ofReal_one, div_one, exp_two_pi_int_fract n x]
    ring
  have hval : 1 / (2 * (Real.pi : ℂ) ^ 2)
      * (-(2 * (Real.pi : ℂ) * Complex.I) ^ 2 / (Nat.factorial 2)
        * (bernoulliFun 2 (Int.fract x) : ℝ)) = P2 x := by
    have hI : (-(2 * (Real.pi : ℂ) * Complex.I) ^ 2 : ℂ) = 4 * (Real.pi : ℂ) ^ 2 := by
      rw [mul_pow, mul_pow, Complex.I_sq]
      ring
    rw [hI, P2]
    have hpi : (Real.pi : ℂ) ≠ 0 := by
      exact_mod_cast Real.pi_ne_zero
    rw [show (Nat.factorial 2 : ℂ) = 2 by norm_num [Nat.factorial]]
    field_simp
    ring
  rw [hfun] at h2
  rwa [hval] at h2

/-! ### The piecewise-linear class in the Bernoulli basis -/

/-- A continuous periodic PL function presented in the periodized-Bernoulli
basis: `W x = w₀ + Σ_i α_i · P₂(x − a_i)`.  (Every continuous 1-periodic PL
function has this form with `α_i = −(jump of W′ at a_i)/2`, `Σα_i = 0`,
`TV(W′) = 2·Σ|α_i|`; the cell engine emits this data directly.) -/
def plW (k : ℕ) (α : Fin k → ℝ) (a : Fin k → ℝ) (w0 : ℂ) : ℝ → ℂ :=
  fun x => w0 + ∑ i, (α i : ℂ) * P2 (x - a i)

/-- The Fourier coefficients of `plW`. -/
def plC (k : ℕ) (α : Fin k → ℝ) (a : Fin k → ℝ) (w0 : ℂ) : ℤ → ℂ :=
  fun n => if n = 0 then w0
    else (∑ i, (α i : ℂ)
      * Complex.exp (-(2 * (Real.pi : ℂ) * Complex.I * n * a i)))
      / (2 * (Real.pi : ℂ) ^ 2 * (n : ℂ) ^ 2)

/-- The pointwise Fourier representation of `plW`. -/
theorem plW_hasSum (k : ℕ) (α : Fin k → ℝ) (a : Fin k → ℝ) (w0 : ℂ) (x : ℝ) :
    HasSum (fun n : ℤ => plC k α a w0 n
      * Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * x))
      (plW k α a w0 x) := by
  -- the α-combination part
  have hterm : ∀ i : Fin k, HasSum
      (fun n : ℤ => (α i : ℂ) * (1 / (2 * (Real.pi : ℂ) ^ 2 * (n : ℂ) ^ 2)
        * Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * ((x - a i : ℝ) : ℂ))))
      ((α i : ℂ) * P2 (x - a i)) :=
    fun i => (hasSum_P2 (x - a i)).mul_left ((α i : ℂ))
  have hsum : HasSum
      (fun n : ℤ => ∑ i, (α i : ℂ) * (1 / (2 * (Real.pi : ℂ) ^ 2 * (n : ℂ) ^ 2)
        * Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * ((x - a i : ℝ) : ℂ))))
      (∑ i, (α i : ℂ) * P2 (x - a i)) :=
    hasSum_sum fun i _ => hterm i
  -- the constant part, concentrated at n = 0
  have hconst : HasSum (fun n : ℤ => if n = 0 then w0 else 0) w0 :=
    hasSum_ite_eq 0 w0
  have htot := hconst.add hsum
  -- identify the total series with plC · exp
  have hfun : (fun n : ℤ => (if n = 0 then w0 else 0)
      + ∑ i, (α i : ℂ) * (1 / (2 * (Real.pi : ℂ) ^ 2 * (n : ℂ) ^ 2)
        * Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * ((x - a i : ℝ) : ℂ))))
      = fun n : ℤ => plC k α a w0 n
      * Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * x) := by
    funext n
    by_cases hn : n = 0
    · subst hn
      simp [plC]
    · rw [if_neg hn, plC, if_neg hn, zero_add]
      rw [Finset.sum_div, Finset.sum_mul]
      congr 1
      funext i
      rw [show 2 * (Real.pi : ℂ) * Complex.I * n * ((x - a i : ℝ) : ℂ)
          = 2 * (Real.pi : ℂ) * Complex.I * n * x
            + -(2 * (Real.pi : ℂ) * Complex.I * n * (a i : ℝ)) by push_cast; ring,
        Complex.exp_add]
      ring
    -- (note: the `congr 1; funext i` above works since both sides are Fin k sums)
  rw [hfun] at htot
  exact htot

/-- The coefficient bound: `‖c_n‖ ≤ (2Σ|α_i|)/(4π²n²)` for `n ≠ 0` —
the THM-665 BV constant with `TV(W′) = 2Σ|α_i|`. -/
theorem plC_norm_le (k : ℕ) (α : Fin k → ℝ) (a : Fin k → ℝ) (w0 : ℂ)
    (n : ℤ) (hn : n ≠ 0) :
    ‖plC k α a w0 n‖ ≤ (2 * ∑ i, |α i|) / (4 * Real.pi ^ 2 * (n : ℝ) ^ 2) := by
  rw [plC, if_neg hn, norm_div]
  have hnum : ‖∑ i, (α i : ℂ)
      * Complex.exp (-(2 * (Real.pi : ℂ) * Complex.I * n * a i))‖
      ≤ ∑ i, |α i| := by
    refine le_trans (norm_sum_le _ _) (le_of_eq (Finset.sum_congr rfl fun i _ => ?_))
    rw [norm_mul]
    have hexp : ‖Complex.exp (-(2 * (Real.pi : ℂ) * Complex.I * n * a i))‖ = 1 := by
      rw [show -(2 * (Real.pi : ℂ) * Complex.I * n * a i)
          = ((-(2 * Real.pi * n * a i) : ℝ) : ℂ) * Complex.I by push_cast; ring]
      exact Complex.norm_exp_ofReal_mul_I _
    rw [hexp, mul_one, Complex.norm_real, Real.norm_eq_abs]
  have hden : ‖(2 * (Real.pi : ℂ) ^ 2 * (n : ℂ) ^ 2)‖ = 2 * Real.pi ^ 2 * (n : ℝ) ^ 2 := by
    rw [show (2 * (Real.pi : ℂ) ^ 2 * (n : ℂ) ^ 2)
        = ((2 * Real.pi ^ 2 * (n : ℝ) ^ 2 : ℝ) : ℂ) by push_cast; ring,
      Complex.norm_real, Real.norm_eq_abs, abs_of_pos]
    have hnR : ((n : ℝ)) ≠ 0 := Int.cast_ne_zero.mpr hn
    positivity
  rw [hden]
  have hnR : ((n : ℝ)) ≠ 0 := Int.cast_ne_zero.mpr hn
  have hpi : Real.pi ≠ 0 := Real.pi_ne_zero
  calc ‖∑ i, (α i : ℂ) * Complex.exp (-(2 * (Real.pi : ℂ) * Complex.I * n * a i))‖
      / (2 * Real.pi ^ 2 * (n : ℝ) ^ 2)
      ≤ (∑ i, |α i|) / (2 * Real.pi ^ 2 * (n : ℝ) ^ 2) := by
        first
        | (gcongr; exact hnum)
        | gcongr
    _ = (2 * ∑ i, |α i|) / (4 * Real.pi ^ 2 * (n : ℝ) ^ 2) := by
        field_simp
        ring

/-- Absolute summability of the coefficients. -/
theorem plC_summable (k : ℕ) (α : Fin k → ℝ) (a : Fin k → ℝ) (w0 : ℂ) :
    Summable fun n : ℤ => ‖plC k α a w0 n‖ := by
  have hmaj : Summable fun n : ℤ =>
      (2 * ∑ i, |α i|) / (4 * Real.pi ^ 2) * (1 / (n : ℝ) ^ 2)
        + (if n = 0 then ‖w0‖ else 0) := by
    apply Summable.add
    · exact ((Real.summable_one_div_int_pow.mpr one_lt_two).mul_left _)
    · exact (hasSum_ite_eq 0 ‖w0‖).summable
  apply Summable.of_nonneg_of_le (fun n => norm_nonneg _) _ hmaj
  intro n
  by_cases hn : n = 0
  · subst hn
    simp [plC]
  · rw [if_neg hn, add_zero]
    have h1 := plC_norm_le k α a w0 n hn
    have heq : (2 * ∑ i, |α i|) / (4 * Real.pi ^ 2 * (n : ℝ) ^ 2)
        = (2 * ∑ i, |α i|) / (4 * Real.pi ^ 2) * (1 / (n : ℝ) ^ 2) := by
      have hnR : ((n : ℝ)) ≠ 0 := Int.cast_ne_zero.mpr hn
      have hpi : Real.pi ≠ 0 := Real.pi_ne_zero
      field_simp
    rwa [heq] at h1

/-- **THM-665, complete.**  For a continuous 1-periodic piecewise-linear `W`
presented in the periodized-Bernoulli basis (`W x = w₀ + Σ_i α_i P₂(x−a_i)`,
mean `w₀`, `TV(W′) = 2Σ|α_i|`), the `V`-grid average satisfies

  `‖E_grid[W](V) − ∫₀¹W‖ ≤ TV(W′)/(12·V²)`. -/
theorem thm665_full (k : ℕ) (α : Fin k → ℝ) (a : Fin k → ℝ) (w0 : ℂ)
    (V : ℕ) (hV : 0 < V) :
    ‖(V : ℂ)⁻¹ * ∑ j ∈ range V, plW k α a w0 ((j : ℝ) / (V : ℝ)) - w0‖
      ≤ (2 * ∑ i, |α i|) / (12 * (V : ℝ) ^ 2) := by
  have hC : (0 : ℝ) ≤ 2 * ∑ i, |α i| := by positivity
  have h := aliasing_tail_bound (plW k α a w0) (plC k α a w0)
    (2 * ∑ i, |α i|) V hV hC (plC_summable k α a w0)
    (plW_hasSum k α a w0) (plC_norm_le k α a w0)
  simpa [plC] using h

end

end Aliasing
end LonelyRunner

-- kernel-purity audit
#print axioms LonelyRunner.Aliasing.hasSum_P2
#print axioms LonelyRunner.Aliasing.plW_hasSum
#print axioms LonelyRunner.Aliasing.plC_norm_le
#print axioms LonelyRunner.Aliasing.thm665_full
