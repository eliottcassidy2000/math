/-
  TournamentH7.LRCAliasingBound  (boxeph-2026-07-09-S6, HYP-5778)

  THM-665 (monad-explorer-S1), the sharp aliasing bound — the analytic core
  in Lean, with the piecewise-linear Fourier input left as named hypotheses.

  Paper statement: for the 1-periodic continuous piecewise-linear uncovered
  measure `W`,  |E_grid[W](V) − ∫₀¹ W| ≤ TV(W′)/(12·V²).  Here (abstract
  coefficients `c : ℤ → ℂ`, `c 0` = the mean):

    * `grid_char_sum` — grid orthogonality `Σ_{j<V} e(2πi n j/V) = V·1[V∣n]`
      (finite geometric sum).
    * `grid_aliasing` — the Poisson-aliasing identity
      `(1/V)·Σ_{j<V} W(j/V) = Σ'_m c(m·V)`, conditional on the pointwise
      absolutely-summable Fourier representation of `W`.
    * `aliasing_tail_bound` — with the BV decay `‖c n‖ ≤ C/(4π²n²)` (n ≠ 0),
      `‖grid average − c 0‖ ≤ C/(12·V²)`: the SHARP constant
      `1/12 = 2·(π²/6)/(4π²)` via the Basel sum `hasSum_zeta_two`.

  Remaining for full THM-665 (named, feeds the hypotheses here): for W
  piecewise linear with breakpoint data, (a) the Fourier representation with
  absolute summability, and (b) `‖c n‖ ≤ TV(W′)/(4π²n²)` — W′ is a STEP
  function, so (b) is two discrete summations by parts over finitely many
  jumps (finite explicit exponential integrals; no BV theory).

  Kernel-pure target: no `sorry`, no `native_decide` (audited below).
-/
import Mathlib
import TournamentH7.LonelyRunner

open Complex Finset

namespace LonelyRunner
namespace Aliasing

noncomputable section

/-- The grid character `e(2πi·n·j/V)`. -/
def gchar (V : ℕ) (n : ℤ) (j : ℕ) : ℂ :=
  Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * j / V)

/-- **Grid orthogonality.**  `Σ_{j<V} e(2πi n j/V)` is `V` if `V ∣ n`, else `0`. -/
theorem grid_char_sum (V : ℕ) (hV : 0 < V) (n : ℤ) :
    ∑ j ∈ range V, gchar V n j = if (V : ℤ) ∣ n then (V : ℂ) else 0 := by
  have hVne : (V : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hV.ne'
  have hw : (2 : ℂ) * (Real.pi : ℂ) * Complex.I ≠ 0 := by
    simp [Real.pi_ne_zero, Complex.I_ne_zero]
  set z : ℂ := Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n / V) with hz
  have hpow : ∀ j : ℕ, gchar V n j = z ^ j := by
    intro j
    rw [hz, ← Complex.exp_nat_mul]
    show Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * j / V) = _
    congr 1
    push_cast
    ring
  by_cases hdvd : (V : ℤ) ∣ n
  · obtain ⟨m, rfl⟩ := hdvd
    have hz1 : z = 1 := by
      rw [hz]
      have harg : 2 * (Real.pi : ℂ) * Complex.I * ((V * m : ℤ) : ℂ) / V
          = (m : ℂ) * (2 * (Real.pi : ℂ) * Complex.I) := by
        first
        | (push_cast; field_simp; ring)
        | (push_cast; field_simp)
        | push_cast
      rw [harg]
      exact Complex.exp_int_mul_two_pi_mul_I m
    have hsum : ∑ j ∈ range V, gchar V (V * m) j = (V : ℂ) := by
      calc ∑ j ∈ range V, gchar V (V * m) j
          = ∑ j ∈ range V, z ^ j := sum_congr rfl fun j _ => hpow j
        _ = ∑ j ∈ range V, 1 := by rw [hz1]; simp
        _ = (V : ℂ) := by simp
    rw [hsum, if_pos ⟨m, rfl⟩]
  · have hzne : z ≠ 1 := by
      intro h1
      apply hdvd
      rw [hz, Complex.exp_eq_one_iff] at h1
      obtain ⟨k, hk⟩ := h1
      rw [div_eq_iff hVne] at hk
      have h2 : (2 : ℂ) * (Real.pi : ℂ) * Complex.I * (n : ℂ)
          = (2 : ℂ) * (Real.pi : ℂ) * Complex.I * ((k * V : ℤ) : ℂ) := by
        first
        | (push_cast; linear_combination hk)
        | push_cast
      have h3 : (n : ℂ) = ((k * V : ℤ) : ℂ) := mul_left_cancel₀ hw h2
      have h4 : n = k * V := by exact_mod_cast h3
      exact ⟨k, by linarith⟩
    have hzV : z ^ V = 1 := by
      rw [hz, ← Complex.exp_nat_mul]
      have harg : (V : ℂ) * (2 * (Real.pi : ℂ) * Complex.I * n / V)
          = (n : ℂ) * (2 * (Real.pi : ℂ) * Complex.I) := by
        first
        | (field_simp; ring)
        | field_simp
      rw [harg]
      exact Complex.exp_int_mul_two_pi_mul_I n
    calc ∑ j ∈ range V, gchar V n j
        = ∑ j ∈ range V, z ^ j := sum_congr rfl fun j _ => hpow j
      _ = (z ^ V - 1) / (z - 1) := geom_sum_eq hzne V
      _ = 0 := by rw [hzV]; simp
      _ = if (V : ℤ) ∣ n then (V : ℂ) else 0 := by rw [if_neg hdvd]

/-- **The Poisson-aliasing identity (THM-665 step (i), conditional).**  If `W`
has the pointwise absolutely-summable Fourier representation, its `V`-grid
average equals the aliased coefficient sum `Σ'_m c(m·V)`. -/
theorem grid_aliasing (W : ℝ → ℂ) (c : ℤ → ℂ) (V : ℕ) (hV : 0 < V)
    (habs : Summable fun n : ℤ => ‖c n‖)
    (hrep : ∀ x : ℝ, HasSum
      (fun n : ℤ => c n * Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * x)) (W x)) :
    (V : ℂ)⁻¹ * ∑ j ∈ range V, W ((j : ℝ) / (V : ℝ)) = ∑' m : ℤ, c (m * V) := by
  have hVne : (V : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hV.ne'
  have hVZne : (V : ℤ) ≠ 0 := by exact_mod_cast hV.ne'
  have hinj : Function.Injective (fun m : ℤ => m * V) := fun a b hab =>
    mul_right_cancel₀ hVZne hab
  -- grid values as tsums
  have hgrid : ∀ j : ℕ, W ((j : ℝ) / (V : ℝ)) = ∑' n : ℤ, c n * gchar V n j := by
    intro j
    rw [← (hrep ((j : ℝ) / (V : ℝ))).tsum_eq]
    congr 1
    funext n
    show c n * Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * (((j : ℝ) / (V : ℝ) : ℝ) : ℂ))
        = c n * Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * j / V)
    congr 2
    push_cast
    ring
  have hnorm : ∀ (n : ℤ) (j : ℕ), ‖gchar V n j‖ = 1 := by
    intro n j
    have harg : 2 * (Real.pi : ℂ) * Complex.I * n * j / V
        = ((2 * Real.pi * n * j / V : ℝ) : ℂ) * Complex.I := by
      push_cast
      ring
    show ‖Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * j / V)‖ = 1
    rw [harg, Complex.norm_exp_ofReal_mul_I]
  have hsummable : ∀ j ∈ range V, Summable fun n : ℤ => c n * gchar V n j := by
    intro j _
    apply Summable.of_norm
    apply Summable.of_nonneg_of_le (fun n => norm_nonneg _) _ habs
    intro n
    rw [norm_mul, hnorm n j, mul_one]
  have hsupp : Function.support (fun n : ℤ => if (V : ℤ) ∣ n then c n else 0)
      ⊆ Set.range (fun m : ℤ => m * V) := by
    intro n hn
    by_cases hdvd : (V : ℤ) ∣ n
    · obtain ⟨m, rfl⟩ := hdvd
      exact ⟨m, by ring⟩
    · simp [Function.mem_support, hdvd] at hn
  have key : ∑ j ∈ range V, W ((j : ℝ) / (V : ℝ)) = (V : ℂ) * ∑' m : ℤ, c (m * V) := by
    calc ∑ j ∈ range V, W ((j : ℝ) / (V : ℝ))
        = ∑ j ∈ range V, ∑' n : ℤ, c n * gchar V n j :=
          sum_congr rfl fun j _ => hgrid j
      _ = ∑' n : ℤ, ∑ j ∈ range V, c n * gchar V n j :=
          (Summable.tsum_finsetSum hsummable).symm
      _ = ∑' n : ℤ, c n * ∑ j ∈ range V, gchar V n j := by
          congr 1
          funext n
          rw [mul_sum]
      _ = ∑' n : ℤ, (if (V : ℤ) ∣ n then (V : ℂ) * c n else 0) := by
          congr 1
          funext n
          rw [grid_char_sum V hV n]
          split_ifs <;> ring
      _ = (V : ℂ) * ∑' n : ℤ, (if (V : ℤ) ∣ n then c n else 0) := by
          rw [← tsum_mul_left]
          congr 1
          funext n
          split_ifs <;> ring
      _ = (V : ℂ) * ∑' m : ℤ, c (m * V) := by
          congr 1
          rw [← Function.Injective.tsum_eq hinj hsupp]
          congr 1
          funext m
          show (if (V : ℤ) ∣ m * V then c (m * V) else 0) = c (m * V)
          rw [if_pos (dvd_mul_left (V : ℤ) m)]
  rw [key, ← mul_assoc, inv_mul_cancel₀ hVne, one_mul]

/-- **The sharp tail bound (THM-665 step (iii)).**  Under the representation
hypotheses, BV decay `‖c n‖ ≤ C/(4π²n²)` (n ≠ 0) gives
`‖grid average − c 0‖ ≤ C/(12·V²)` — the constant `1/12 = 2·(π²/6)/(4π²)`
via the Basel value (`hasSum_zeta_two`). -/
theorem aliasing_tail_bound (W : ℝ → ℂ) (c : ℤ → ℂ) (C : ℝ) (V : ℕ)
    (hV : 0 < V) (hC : 0 ≤ C)
    (habs : Summable fun n : ℤ => ‖c n‖)
    (hrep : ∀ x : ℝ, HasSum
      (fun n : ℤ => c n * Complex.exp (2 * (Real.pi : ℂ) * Complex.I * n * x)) (W x))
    (hcoef : ∀ n : ℤ, n ≠ 0 → ‖c n‖ ≤ C / (4 * Real.pi ^ 2 * (n : ℝ) ^ 2)) :
    ‖(V : ℂ)⁻¹ * ∑ j ∈ range V, W ((j : ℝ) / (V : ℝ)) - c 0‖
      ≤ C / (12 * (V : ℝ) ^ 2) := by
  have hVZne : (V : ℤ) ≠ 0 := by exact_mod_cast hV.ne'
  have hVRne : (V : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hV.ne'
  have hinj : Function.Injective (fun m : ℤ => m * V) := fun a b hab =>
    mul_right_cancel₀ hVZne hab
  rw [grid_aliasing W c V hV habs hrep]
  -- summability of the aliased sequence
  have hsubnorm : Summable fun m : ℤ => ‖c (m * V)‖ := habs.comp_injective hinj
  have hsub : Summable fun m : ℤ => c (m * V) := hsubnorm.of_norm
  -- split off m = 0 (c (0·V) = c 0)
  rw [hsub.tsum_eq_add_tsum_ite 0]
  simp only [zero_mul]
  rw [add_sub_cancel_left]
  -- the scale K and the integer Basel majorant
  set f : ℤ → ℝ := fun m => if m = 0 then (0 : ℝ) else 1 / (m : ℝ) ^ 2 with hf
  set K : ℝ := C / (4 * Real.pi ^ 2 * (V : ℝ) ^ 2) with hK
  have hKnn : 0 ≤ K := by
    rw [hK]
    positivity
  -- Basel over ℕ, both branches
  have hbasel : HasSum (fun n : ℕ => (1 : ℝ) / (n : ℝ) ^ 2) (Real.pi ^ 2 / 6) :=
    hasSum_zeta_two
  have hpos : HasSum (fun n : ℕ => f (n : ℤ)) (Real.pi ^ 2 / 6) := by
    have heq : (fun n : ℕ => f (n : ℤ)) = fun n : ℕ => (1 : ℝ) / (n : ℝ) ^ 2 := by
      funext n
      rcases Nat.eq_zero_or_pos n with h | h
      · subst h
        simp [hf]
      · rw [hf]
        simp only []
        rw [if_neg (by exact_mod_cast h.ne')]
        push_cast
        ring
    rw [heq]
    exact hbasel
  have hshift : HasSum (fun n : ℕ => (1 : ℝ) / ((n : ℝ) + 1) ^ 2) (Real.pi ^ 2 / 6) := by
    have h1 : HasSum (fun n : ℕ => (1 : ℝ) / ((n + 1 : ℕ) : ℝ) ^ 2) (Real.pi ^ 2 / 6) := by
      apply (hasSum_nat_add_iff (f := fun n : ℕ => (1 : ℝ) / (n : ℝ) ^ 2) 1).mpr
      simpa using hbasel
    have heq : (fun n : ℕ => (1 : ℝ) / ((n + 1 : ℕ) : ℝ) ^ 2)
        = fun n : ℕ => (1 : ℝ) / ((n : ℝ) + 1) ^ 2 := by
      funext n
      push_cast
      ring
    rwa [heq] at h1
  have hneg : HasSum (fun n : ℕ => f (-((n : ℤ) + 1))) (Real.pi ^ 2 / 6) := by
    have heq : (fun n : ℕ => f (-((n : ℤ) + 1)))
        = fun n : ℕ => (1 : ℝ) / ((n : ℝ) + 1) ^ 2 := by
      funext n
      rw [hf]
      simp only []
      rw [if_neg (by omega)]
      push_cast
      ring
    rw [heq]
    exact hshift
  have hint : HasSum f (Real.pi ^ 2 / 3) := by
    have h := HasSum.of_nat_of_neg_add_one hpos hneg
    have h36 : Real.pi ^ 2 / 6 + Real.pi ^ 2 / 6 = Real.pi ^ 2 / 3 := by ring
    rwa [h36] at h
  have hmaj : HasSum (fun m : ℤ => K * f m) (K * (Real.pi ^ 2 / 3)) := hint.mul_left K
  -- pointwise bound of the split tail by the majorant
  have hb : ∀ m : ℤ, ‖(if m = 0 then (0 : ℂ) else c (m * V))‖ ≤ K * f m := by
    intro m
    by_cases hm : m = 0
    · simp [hm, hf]
    · rw [if_neg hm, hf]
      simp only []
      rw [if_neg hm]
      have hmV : (m * V : ℤ) ≠ 0 := mul_ne_zero hm hVZne
      have h1 : ‖c (m * V)‖ ≤ C / (4 * Real.pi ^ 2 * ((m * V : ℤ) : ℝ) ^ 2) :=
        hcoef _ hmV
      have hmR : ((m : ℝ)) ≠ 0 := Int.cast_ne_zero.mpr hm
      have h2 : C / (4 * Real.pi ^ 2 * ((m * V : ℤ) : ℝ) ^ 2)
          = K * (1 / (m : ℝ) ^ 2) := by
        rw [hK]
        first
        | (push_cast; field_simp; ring)
        | (push_cast; field_simp)
      rw [h2] at h1
      exact h1
  -- summability of the split tail norms
  have hsummite : Summable fun m : ℤ => ‖(if m = 0 then (0 : ℂ) else c (m * V))‖ := by
    apply Summable.of_nonneg_of_le (fun m => norm_nonneg _) _ hsubnorm
    intro m
    by_cases hm : m = 0 <;> simp [hm]
  -- assemble
  calc ‖∑' m : ℤ, (if m = 0 then (0 : ℂ) else c (m * V))‖
      ≤ ∑' m : ℤ, ‖(if m = 0 then (0 : ℂ) else c (m * V))‖ :=
        norm_tsum_le_tsum_norm hsummite
    _ ≤ ∑' m : ℤ, K * f m :=
        hsummite.tsum_le_tsum hb hmaj.summable
    _ = K * (Real.pi ^ 2 / 3) := hmaj.tsum_eq
    _ = C / (12 * (V : ℝ) ^ 2) := by
        rw [hK]
        have hpi : Real.pi ≠ 0 := Real.pi_ne_zero
        field_simp
        ring

end

end Aliasing
end LonelyRunner

-- kernel-purity audit
#print axioms LonelyRunner.Aliasing.grid_char_sum
#print axioms LonelyRunner.Aliasing.grid_aliasing
#print axioms LonelyRunner.Aliasing.aliasing_tail_bound
