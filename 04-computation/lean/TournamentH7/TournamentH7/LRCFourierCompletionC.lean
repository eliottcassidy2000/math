/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-11-S214)
-/
import Mathlib
import TournamentH7.LRCFourierCompletionB

set_option maxHeartbeats 1200000

/-!
# LEM-022 Fourier completion, Stage B.3 (opus-S214) — the completion identity

The final analytic piece of the t2 hyperbola bound (LEM-022, death-star-S9).  Builds on Stage B.1/B.2
(`LRCFourierCompletionB`, opus-S213): the additive orthogonality `Σ_{x<q} e_q(hx) = q·1{q∣h}` and the band
coefficient bound `‖B̂(h)‖ ≤ q/(2·cdist h)`.

**The completion identity** (canon LEM-022, "Completion" step): for a band `B ⊆ [0,q)` and a twist `w`,
the multiplicative pair-correlation `C_w := #{s ∈ B : (w·s) mod q ∈ B}` satisfies

> `(C_w : ℂ) = (1/q)·Σ_{k<q} B̂(k)·conj(B̂(w·k))`,   `B̂(j) := Σ_{r∈B} e_q(−j·r)`.

Everything is over the **integers on `range q`** (not `ZMod q` characters) — `e_q` is `q`-periodic, so the
twist `w·k` needs no reduction, and the orthogonality is exactly B.1. Derivation: expand both coefficients,
swap the `k`-sum inside, collapse `Σ_{k<q} e_q(k·(ws−r)) = q·1{q ∣ ws−r}` (B.1), and count the surviving
pairs `(r,s)` with `ws ≡ r`.

Kernel-pure: no `sorry`, no `native_decide`.
-/

namespace LonelyRunner
namespace FourierCompletion

open Complex Real

/-- `e_q(n) := exp(2π·n·i/q)`, integer argument — the summand of B.1/B.2 packaged. -/
noncomputable def eInt (q : ℕ) (n : ℤ) : ℂ :=
  Complex.exp ((2 * π * (n : ℝ) / q : ℝ) * Complex.I)

/-- `e_q` is an additive character: `e_q(a)·e_q(b) = e_q(a+b)`. -/
theorem eInt_add (q : ℕ) (a b : ℤ) : eInt q a * eInt q b = eInt q (a + b) := by
  rw [eInt, eInt, eInt, ← Complex.exp_add]; congr 1; push_cast; ring

/-- `conj(e_q(n)) = e_q(−n)`. -/
theorem eInt_conj (q : ℕ) (n : ℤ) : (starRingEnd ℂ) (eInt q n) = eInt q (-n) := by
  rw [eInt, eInt, ← Complex.exp_conj]
  congr 1
  simp only [map_mul, Complex.conj_I, map_div₀, Complex.conj_ofReal]
  push_cast
  ring

/-- **B.1 in `eInt` form.** `Σ_{k<q} e_q(k·h) = q·1{q ∣ h}`. -/
theorem eInt_orthog (q : ℕ) (hq : 0 < q) (h : ℤ) :
    (∑ k ∈ Finset.range q, eInt q ((k : ℤ) * h)) = if (q : ℤ) ∣ h then (q : ℂ) else 0 := by
  rw [← sum_exp_orthogonality q hq h]
  refine Finset.sum_congr rfl (fun x _ => ?_)
  rw [eInt]; congr 1; push_cast; ring

/-- The band DFT coefficient `B̂(j) = Σ_{r∈B} e_q(−j·r)`. -/
noncomputable def bandDFT (q : ℕ) (B : Finset ℕ) (j : ℤ) : ℂ :=
  ∑ r ∈ B, eInt q (-(j * (r : ℤ)))

/-- The multiplicative pair-correlation `C_w = #{s ∈ B : (w·s) mod q ∈ B}`. -/
def corrCount (q : ℕ) (B : Finset ℕ) (w : ℕ) : ℕ :=
  (B.filter (fun s => (w * s) % q ∈ B)).card

/-- for `r < q`, the off-line congruence `q ∣ w·s − r` pins `r` to the residue `(w·s) mod q`. -/
theorem dvd_sub_iff_eq_emod (q : ℕ) (hq : 0 < q) (w s r : ℕ) (hr : r < q) :
    (q : ℤ) ∣ ((w : ℤ) * (s : ℤ) - (r : ℤ)) ↔ r = (w * s) % q := by
  rw [Int.dvd_iff_emod_eq_zero, ← Int.emod_eq_emod_iff_emod_sub_eq_zero,
    show (w : ℤ) * (s : ℤ) = ((w * s : ℕ) : ℤ) from by push_cast; ring,
    ← Int.natCast_mod,
    Int.emod_eq_of_lt (show (0 : ℤ) ≤ (r : ℤ) from by positivity)
      (show (r : ℤ) < (q : ℤ) from by exact_mod_cast hr),
    eq_comm]
  exact Int.natCast_inj

/-- **B.3 — the completion identity.**  For a band `B ⊆ [0,q)` and a twist `w`, the multiplicative
pair-correlation equals the diagonal-suppressed Fourier form:
`(C_w : ℂ) = (1/q)·Σ_{k<q} B̂(k)·conj(B̂(w·k))`. -/
theorem completion_identity (q : ℕ) (hq : 0 < q) (B : Finset ℕ) (hB : B ⊆ Finset.range q) (w : ℕ) :
    (corrCount q B w : ℂ)
      = (1 / (q : ℂ)) * ∑ k ∈ Finset.range q,
          bandDFT q B (k : ℤ) * (starRingEnd ℂ) (bandDFT q B ((w : ℤ) * (k : ℤ))) := by
  have hqc : (q : ℂ) ≠ 0 := by exact_mod_cast hq.ne'
  -- expand each summand to `Σ_r Σ_s e_q(k·(ws−r))`
  have hsummand : ∀ k : ℕ, bandDFT q B (k : ℤ) * (starRingEnd ℂ) (bandDFT q B ((w : ℤ) * (k : ℤ)))
      = ∑ r ∈ B, ∑ s ∈ B, eInt q ((k : ℤ) * ((w : ℤ) * (s : ℤ) - (r : ℤ))) := by
    intro k
    rw [bandDFT, bandDFT, map_sum, Finset.sum_mul_sum]
    refine Finset.sum_congr rfl (fun r _ => Finset.sum_congr rfl (fun s _ => ?_))
    rw [eInt_conj, neg_neg, eInt_add]
    congr 1
    push_cast; ring
  -- the whole sum, reordered to `Σ_r Σ_s Σ_k`, then collapse the `k`-sum by orthogonality
  have hcollapsed : (∑ k ∈ Finset.range q,
        bandDFT q B (k : ℤ) * (starRingEnd ℂ) (bandDFT q B ((w : ℤ) * (k : ℤ))))
      = ∑ r ∈ B, ∑ s ∈ B, (if (q : ℤ) ∣ ((w : ℤ) * (s : ℤ) - (r : ℤ)) then (q : ℂ) else 0) := by
    rw [Finset.sum_congr rfl (fun k _ => hsummand k)]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl (fun r _ => ?_)
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl (fun s _ => ?_)
    exact eInt_orthog q hq ((w : ℤ) * (s : ℤ) - (r : ℤ))
  rw [hcollapsed]
  -- count: for each `s`, `Σ_r [q ∣ ws−r]·q = q·1{ws mod q ∈ B}`; sum over `s` gives `q·C_w`
  have hswap : (∑ r ∈ B, ∑ s ∈ B, (if (q : ℤ) ∣ ((w : ℤ) * (s : ℤ) - (r : ℤ)) then (q : ℂ) else 0))
      = ∑ s ∈ B, ∑ r ∈ B, (if (q : ℤ) ∣ ((w : ℤ) * (s : ℤ) - (r : ℤ)) then (q : ℂ) else 0) :=
    Finset.sum_comm
  rw [hswap]
  have hper_s : ∀ s ∈ B,
      (∑ r ∈ B, (if (q : ℤ) ∣ ((w : ℤ) * (s : ℤ) - (r : ℤ)) then (q : ℂ) else 0))
        = if (w * s) % q ∈ B then (q : ℂ) else 0 := by
    intro s _
    rw [Finset.sum_ite, Finset.sum_const_zero, add_zero, Finset.sum_const,
      Finset.filter_congr
        (fun (r : ℕ) hr => dvd_sub_iff_eq_emod q hq w s r (Finset.mem_range.mp (hB hr))),
      Finset.filter_eq']
    by_cases h : (w * s) % q ∈ B <;> simp [h]
  rw [Finset.sum_congr rfl hper_s]
  -- `Σ_s 1{ws mod q ∈ B}·q = q·C_w`
  rw [Finset.sum_ite, Finset.sum_const_zero, add_zero, Finset.sum_const, nsmul_eq_mul]
  show (((B.filter (fun s => (w * s) % q ∈ B)).card : ℕ) : ℂ)
    = (1 / (q : ℂ)) * ((B.filter (fun s => (w * s) % q ∈ B)).card * (q : ℂ))
  field_simp

end FourierCompletion
end LonelyRunner

