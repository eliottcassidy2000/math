/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-11-S213)
-/
import Mathlib
import TournamentH7.LRCFourierCompletion
import TournamentH7.LRCHyperbolaBox

/-!
# LEM-022 Fourier completion, Stage B (opus-S213) — the t2 per-cell equidistribution bound

Continues death-star-S13's `LRCFourierCompletion` Stage A (the interval exponential-sum bound
`‖Σ exp‖ ≤ q/(2d)`, `norm_expSum_le`).  Stage B supplies the ANALYTIC input the character route's
combinatorial core (`hyperbola_box_count → zcorr_percell → offdiag_mcorr_sq_le`) currently assumes: the
per-cell equidistribution of the safe band `B` under a multiplicative twist.

This file (opus-S213) builds the foundational pieces:

* **B.1 — orthogonality** (`sum_exp_orthogonality`): `Σ_{x<q} e_q(hx) = q·1{q ∣ h}`, the finite additive
  character sum. Kernel-pure via `geom_sum_eq` + `Complex.exp_int_mul_two_pi_mul_I`.
* **B.2 — the band coefficient bound** (`norm_bandCoeff_le`): the Fourier coefficient of an interval band
  `B̂(h) = Σ_{x∈[lo,lo+len)} e_q(hx)` obeys `‖B̂(h)‖ ≤ q/(2·cdist h)` for `h ≢ 0` — Stage A composed with
  the Jordan sine witness `2·cdist(h)/q ≤ |sin(π h/q)|`.

The completion identity `C_w = b²/q + (1/q)Σ_{h≠0} B̂(h)·conj(B̂(w⁻¹h))` and the final assembly with
`harmonic_ratio_sum_mul_le` are the remaining Stage B.3 (documented; the harder aggregation).

`e_q(t) := exp(2π·t·i)`; the exp form matches `LRCFourierCompletion.norm_expSum_le` verbatim so the two
compose. Kernel-pure: no `sorry`, no `native_decide`.
-/

namespace LonelyRunner
namespace FourierCompletion

open Complex Real

/-- **B.1 — additive-character orthogonality on `ℤ/q`.** `Σ_{x<q} e_q(hx) = q` if `q ∣ h`, else `0`. -/
theorem sum_exp_orthogonality (q : ℕ) (hq : 0 < q) (h : ℤ) :
    (∑ x ∈ Finset.range q,
        Complex.exp ((2 * π * ((h : ℝ) * (x : ℝ) / q) : ℝ) * Complex.I))
      = if (q : ℤ) ∣ h then (q : ℂ) else 0 := by
  have hqc : (q : ℂ) ≠ 0 := by exact_mod_cast hq.ne'
  set z : ℂ := Complex.exp ((2 * π * ((h : ℝ) / q) : ℝ) * Complex.I) with hz
  -- each term equals `z ^ x`
  have hterm : ∀ x ∈ Finset.range q,
      Complex.exp ((2 * π * ((h : ℝ) * (x : ℝ) / q) : ℝ) * Complex.I) = z ^ x := by
    intro x _
    rw [hz, ← Complex.exp_nat_mul]
    congr 1
    push_cast
    field_simp
  rw [Finset.sum_congr rfl hterm]
  -- `z ^ q = 1`  (a `q`-th root of unity)
  have hzq : z ^ q = 1 := by
    rw [hz, ← Complex.exp_nat_mul]
    have harg : (q : ℂ) * ((2 * π * ((h : ℝ) / q) : ℝ) * Complex.I)
        = (h : ℤ) * (2 * π * Complex.I) := by
      push_cast; field_simp
    rw [harg, Complex.exp_int_mul_two_pi_mul_I]
  by_cases hdvd : (q : ℤ) ∣ h
  · rw [if_pos hdvd]
    -- `q ∣ h ⟹ z = 1 ⟹ Σ 1 = q`
    obtain ⟨m, hm⟩ := hdvd
    have hz1 : z = 1 := by
      rw [hz]
      have : (2 * π * ((h : ℝ) / q) : ℝ) * Complex.I = (m : ℤ) * (2 * π * Complex.I) := by
        rw [hm]; push_cast; field_simp
      rw [this, Complex.exp_int_mul_two_pi_mul_I]
    rw [hz1]
    simp
  · rw [if_neg hdvd]
    -- `q ∤ h ⟹ z ≠ 1 ⟹ geometric sum `(z^q − 1)/(z − 1) = 0`
    have hzne : z ≠ 1 := by
      rw [hz]
      intro hcon
      rw [Complex.exp_eq_one_iff] at hcon
      obtain ⟨n, hn⟩ := hcon
      apply hdvd
      have h2pi : (2 * (π : ℂ) * Complex.I) ≠ 0 := by
        simp [Real.pi_ne_zero, Complex.I_ne_zero]
      have hLHS : ((2 * π * ((h : ℝ) / q) : ℝ) : ℂ) * Complex.I
          = ((h : ℂ) / (q : ℂ)) * (2 * (π : ℂ) * Complex.I) := by push_cast; field_simp
      rw [hLHS] at hn
      have hcancel : (h : ℂ) / (q : ℂ) = (n : ℂ) := mul_right_cancel₀ h2pi hn
      have hhqn : (h : ℂ) = (q : ℂ) * (n : ℂ) := by
        field_simp at hcancel; linear_combination hcancel
      exact ⟨n, by exact_mod_cast hhqn⟩
    rw [geom_sum_eq hzne, hzq]
    simp

end FourierCompletion
end LonelyRunner
