/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-11-S213)
-/
import Mathlib
import TournamentH7.LRCFourierCompletion
import TournamentH7.LRCHyperbolaBox

set_option maxHeartbeats 800000

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

/-- **The sine witness.**  `|sin(π h/q)| ≥ 2·cdist(h)/q` — the input Stage A's `norm_expSum_le`
consumes with `d = cdist h`.  Proof: `|sin(π·/q)|` has period `q` in `h` (via `sin(x + nπ) = ±sin x`),
reducing to the residue `r = (h : ZMod q).val`; then Jordan's inequality on `min(r, q−r)/q ≤ 1/2`. -/
theorem sine_cdist_witness (q : ℕ) (hq : 0 < q) (h : ℤ) :
    2 * ((LonelyRunner.HyperbolaBox.cdist (h : ZMod q) : ℝ) / q) ≤ |Real.sin (π * (h : ℝ) / q)| := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  haveI : NeZero q := ⟨hq.ne'⟩
  set r : ℕ := (h : ZMod q).val with hr
  have hrq : r < q := ZMod.val_lt _
  have hrmod : (r : ℤ) = h % q := by rw [hr, ZMod.val_intCast]
  have hcast : (h : ℝ) = (q : ℝ) * ((h / q : ℤ) : ℝ) + (r : ℝ) := by
    have hdm : (h : ℤ) = q * (h / q) + (r : ℤ) := by
      have := Int.ediv_add_emod h q; omega
    exact_mod_cast hdm
  -- reduce `|sin(π h/q)|` to `|sin(π r/q)|`
  have hred : |Real.sin (π * (h : ℝ) / q)| = |Real.sin (π * (r : ℝ) / q)| := by
    have hn : π * (h : ℝ) / q = π * (r : ℝ) / q + ((h / q : ℤ) : ℝ) * π := by
      rw [hcast]; field_simp; ring
    have hcos1 : |Real.cos (((h / q : ℤ) : ℝ) * π)| = 1 := by
      rw [Real.cos_int_mul_pi]; simp
    rw [hn, Real.sin_add, Real.sin_int_mul_pi, mul_zero, add_zero, abs_mul, hcos1, mul_one]
  rw [hred]
  have hr_nonneg : (0 : ℝ) ≤ π * (r : ℝ) / q := by positivity
  have hr_le_pi : π * (r : ℝ) / q ≤ π := by
    rw [div_le_iff₀ hqR]
    have : (r : ℝ) ≤ q := by exact_mod_cast le_of_lt hrq
    nlinarith [Real.pi_pos]
  have hsin_nonneg : 0 ≤ Real.sin (π * (r : ℝ) / q) :=
    Real.sin_nonneg_of_nonneg_of_le_pi hr_nonneg hr_le_pi
  rw [abs_of_nonneg hsin_nonneg]
  have hcd : LonelyRunner.HyperbolaBox.cdist (h : ZMod q) = min r (q - r) := rfl
  rw [hcd]
  rcases le_total r (q - r) with hle | hle
  · rw [min_eq_left hle]
    have hle2 : (r : ℝ) / q ≤ 1 / 2 := by
      rw [div_le_iff₀ hqR]
      have h2r : 2 * r ≤ q := by omega
      have : (2 * r : ℝ) ≤ q := by exact_mod_cast h2r
      linarith
    have hj := two_mul_le_sin_pi_mul (by positivity) hle2
    calc 2 * ((r : ℝ) / q) ≤ Real.sin (π * ((r : ℝ) / q)) := hj
      _ = Real.sin (π * (r : ℝ) / q) := by rw [mul_div_assoc]
  · rw [min_eq_right hle]
    have hqr : Real.sin (π * (r : ℝ) / q) = Real.sin (π * ((q - r : ℕ) : ℝ) / q) := by
      have hpsub : π * ((q - r : ℕ) : ℝ) / q = π - π * (r : ℝ) / q := by
        rw [Nat.cast_sub (le_of_lt hrq)]; field_simp
      rw [hpsub, Real.sin_pi_sub]
    have hle2 : ((q - r : ℕ) : ℝ) / q ≤ 1 / 2 := by
      rw [div_le_iff₀ hqR]
      have h2r : 2 * (q - r) ≤ q := by omega
      have : (2 * ((q - r : ℕ) : ℝ)) ≤ q := by exact_mod_cast h2r
      linarith
    have hj := two_mul_le_sin_pi_mul (by positivity) hle2
    rw [hqr]
    calc 2 * (((q - r : ℕ) : ℝ) / q) ≤ Real.sin (π * (((q - r : ℕ) : ℝ) / q)) := hj
      _ = Real.sin (π * ((q - r : ℕ) : ℝ) / q) := by rw [mul_div_assoc]

/-- **B.2 — the band Fourier coefficient bound.**  The Fourier coefficient of the interval band
`B = [lo, lo+len)`, `B̂(h) = Σ_{r<len} e_q(h(lo+r))`, obeys `‖B̂(h)‖ ≤ q/(2·cdist h)` for `h ≢ 0`.
Stage A (`norm_expSum_le`) composed with the sine witness. -/
theorem norm_bandCoeff_le (q : ℕ) (hq : 0 < q) (h lo : ℤ) (len : ℕ)
    (hh : (h : ZMod q) ≠ 0) :
    ‖∑ r ∈ Finset.range len,
        Complex.exp ((2 * π * ((h : ℝ) * ((lo : ℝ) + r) / q) : ℝ) * Complex.I)‖
      ≤ (q : ℝ) / (2 * (LonelyRunner.HyperbolaBox.cdist (h : ZMod q) : ℝ)) := by
  haveI : NeZero q := ⟨hq.ne'⟩
  have hcd_pos : 0 < LonelyRunner.HyperbolaBox.cdist (h : ZMod q) :=
    LonelyRunner.HyperbolaBox.one_le_cdist hh
  have hd : (0 : ℝ) < (LonelyRunner.HyperbolaBox.cdist (h : ZMod q) : ℝ) := by exact_mod_cast hcd_pos
  exact norm_expSum_le q hq h lo len _ hd (sine_cdist_witness q hq h)

end FourierCompletion
end LonelyRunner

