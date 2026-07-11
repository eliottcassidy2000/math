/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-11-S214)
-/
import Mathlib
import TournamentH7.LRCFourierCompletionB
import TournamentH7.LRCFourierAggregation

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

/-- `e_q` is `q`-periodic: `q ∣ a − b ⟹ e_q(a) = e_q(b)`. -/
theorem eInt_periodic (q : ℕ) (hq : 0 < q) (a b : ℤ) (h : (q : ℤ) ∣ (a - b)) :
    eInt q a = eInt q b := by
  obtain ⟨t, ht⟩ := h
  have hqcC : (q : ℂ) ≠ 0 := by exact_mod_cast hq.ne'
  rw [eInt, eInt,
    show (2 * π * (a : ℝ) / q : ℝ) * Complex.I
      = (2 * π * (b : ℝ) / q : ℝ) * Complex.I
        + (2 * π * ((a - b : ℤ) : ℝ) / q : ℝ) * Complex.I from by push_cast; ring,
    Complex.exp_add,
    show (2 * π * ((a - b : ℤ) : ℝ) / q : ℝ) * Complex.I = (t : ℤ) * (2 * π * Complex.I) from by
      rw [ht]; push_cast; rw [mul_div_assoc, mul_div_cancel_left₀ _ hqcC]; ring,
    Complex.exp_int_mul_two_pi_mul_I, mul_one]

/-- The band DFT coefficient `B̂(j) = Σ_{r∈B} e_q(−j·r)`. -/
noncomputable def bandDFT (q : ℕ) (B : Finset ℕ) (j : ℤ) : ℂ :=
  ∑ r ∈ B, eInt q (-(j * (r : ℤ)))

/-- `B̂` is `q`-periodic in its argument. -/
theorem bandDFT_periodic (q : ℕ) (hq : 0 < q) (B : Finset ℕ) (a b : ℤ) (h : (q : ℤ) ∣ (a - b)) :
    bandDFT q B a = bandDFT q B b := by
  rw [bandDFT, bandDFT]
  refine Finset.sum_congr rfl (fun r _ => ?_)
  refine eInt_periodic q hq _ _ ?_
  have hd : (q : ℤ) ∣ (b - a) := by
    rw [show b - a = -(a - b) from by ring]; exact dvd_neg.mpr h
  rw [show -(a * (r : ℤ)) - -(b * (r : ℤ)) = (b - a) * (r : ℤ) from by ring]
  exact Dvd.dvd.mul_right hd (r : ℤ)

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

/-- `B̂(0) = |B|` — the DC coefficient counts the band. -/
theorem bandDFT_zero (q : ℕ) (B : Finset ℕ) : bandDFT q B 0 = (B.card : ℂ) := by
  rw [bandDFT]
  have he0 : eInt q 0 = 1 := by rw [eInt]; norm_num
  simp only [zero_mul, neg_zero, he0, Finset.sum_const, nsmul_eq_mul, mul_one]

/-- **B.3 — the completion difference bound.**  Splitting off the `k = 0` main term `b²/q`
(`b = |B|`) and applying the triangle inequality:
`‖C_w − b²/q‖ ≤ (1/q)·Σ_{k≠0} ‖B̂(k)‖·‖B̂(wk)‖`.  The `k = 0` term of the identity is `B̂(0)² = b²`. -/
theorem completion_diff_bound (q : ℕ) (hq : 0 < q) (B : Finset ℕ) (hB : B ⊆ Finset.range q) (w : ℕ) :
    ‖(corrCount q B w : ℂ) - (B.card : ℂ) ^ 2 / q‖
      ≤ (1 / q) * ∑ k ∈ (Finset.range q).filter (fun k => k ≠ 0),
          ‖bandDFT q B (k : ℤ)‖ * ‖bandDFT q B ((w : ℤ) * (k : ℤ))‖ := by
  have hqc : (q : ℂ) ≠ 0 := by exact_mod_cast hq.ne'
  have hqR : (0 : ℝ) < q := by exact_mod_cast hq
  have h0mem : (0 : ℕ) ∈ Finset.range q := Finset.mem_range.mpr hq
  set F : ℕ → ℂ := fun k => bandDFT q B (k : ℤ) * (starRingEnd ℂ) (bandDFT q B ((w : ℤ) * (k : ℤ)))
    with hF
  -- the `k = 0` summand is `b²`
  have hF0 : F 0 = (B.card : ℂ) ^ 2 := by
    rw [hF]; simp only [Nat.cast_zero, mul_zero, bandDFT_zero, map_natCast]; ring
  -- split off `k = 0`
  have hsplit : (∑ k ∈ Finset.range q, F k)
      = (B.card : ℂ) ^ 2 + ∑ k ∈ (Finset.range q).filter (fun k => k ≠ 0), F k := by
    rw [← Finset.add_sum_erase (Finset.range q) F h0mem, hF0]
    congr 1
    apply Finset.sum_congr _ (fun _ _ => rfl)
    ext k; simp [Finset.mem_erase, Finset.mem_filter, and_comm]
  -- the difference is the off-diagonal sum, scaled
  have hdiff : (corrCount q B w : ℂ) - (B.card : ℂ) ^ 2 / q
      = (1 / q) * ∑ k ∈ (Finset.range q).filter (fun k => k ≠ 0), F k := by
    rw [completion_identity q hq B hB w, ← hF, hsplit]
    field_simp
    ring
  rw [hdiff, norm_mul]
  have h1q : ‖(1 / (q : ℂ))‖ = 1 / q := by
    rw [norm_div, norm_one, Complex.norm_natCast]
  rw [h1q]
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  refine le_trans (norm_sum_le _ _) ?_
  apply Finset.sum_le_sum
  intro k _
  rw [hF, norm_mul, RCLike.norm_conj]

open LonelyRunner.HyperbolaBox in
/-- **B.3 — the final LEM-022 bound.**  Combining the completion identity with kps's off-diagonal
aggregation (`offDiag_bandSum_le_closed`): for a band `B ⊆ [0,q)`, a unit twist `w`, and the ratio-lattice
floor `P`, if the band coefficients obey opus's B.2 bound then
`‖C_w − b²/q‖ ≤ 5q(log₂q+1)²/P` — the t2 hyperbola bound, kernel-pure end to end.

The seam: my identity is over `range q` integers, kps's aggregation over `ZMod q`; the bridge reindexes
`k ↔ h.val` and matches the twist `B̂(w·k) = B̂((w·h).val)` by `bandDFT_periodic`. -/
theorem completion_final {q : ℕ} [NeZero q] (B : Finset ℕ) (hB : B ⊆ Finset.range q)
    (w : ℕ) (hw : IsUnit ((w : ZMod q)))
    (P : ℕ) (hP : 0 < P)
    (hPmin : ∀ h : ZMod q, h ≠ 0 → P ≤ cdist h * cdist ((w : ZMod q) * h))
    (hbc : ∀ h : ZMod q, h ≠ 0 → ‖bandDFT q B ((h.val : ℤ))‖ ≤ (q : ℝ) / (2 * (cdist h : ℝ))) :
    ‖(corrCount q B w : ℂ) - (B.card : ℂ) ^ 2 / q‖
      ≤ 5 * (q : ℝ) * ((Nat.log 2 q : ℝ) + 1) ^ 2 / P := by
  have hq : 0 < q := Nat.pos_of_ne_zero (NeZero.ne q)
  have hqc : (q : ℂ) ≠ 0 := by exact_mod_cast hq.ne'
  have hqR : (0 : ℝ) < q := by exact_mod_cast hq
  have hPR : (0 : ℝ) < P := by exact_mod_cast hP
  set bc : ZMod q → ℂ := fun h => bandDFT q B ((h.val : ℤ)) with hbcdef
  set F : ℕ → ℂ := fun k => bandDFT q B (k : ℤ) * (starRingEnd ℂ) (bandDFT q B ((w : ℤ) * (k : ℤ)))
    with hFdef
  have h0mem : (0 : ℕ) ∈ Finset.range q := Finset.mem_range.mpr hq
  -- k = 0 split
  have hF0 : F 0 = (B.card : ℂ) ^ 2 := by
    rw [hFdef]; simp only [Nat.cast_zero, mul_zero, bandDFT_zero, map_natCast]; ring
  have hsplit : (∑ k ∈ Finset.range q, F k)
      = (B.card : ℂ) ^ 2 + ∑ k ∈ (Finset.range q).filter (fun k => k ≠ 0), F k := by
    rw [← Finset.add_sum_erase (Finset.range q) F h0mem, hF0]
    congr 1
    apply Finset.sum_congr _ (fun _ _ => rfl)
    ext k; simp [Finset.mem_erase, Finset.mem_filter, and_comm]
  -- reindex `k ≠ 0` (range q) to `h ≠ 0` (ZMod q), matching the twist by periodicity
  have hoffeq : (∑ k ∈ (Finset.range q).filter (fun k => k ≠ 0), F k)
      = ∑ h ∈ (Finset.univ : Finset (ZMod q)).filter (fun h => h ≠ 0),
          bc h * (starRingEnd ℂ) (bc ((w : ZMod q) * h)) := by
    refine Finset.sum_nbij' (fun k => (k : ZMod q)) (fun h => h.val) ?_ ?_ ?_ ?_ ?_
    · intro k hk
      rw [Finset.mem_filter, Finset.mem_range] at hk
      simp only [Finset.mem_filter, Finset.mem_univ, true_and]
      intro hz
      exact hk.2 (Nat.eq_zero_of_dvd_of_lt ((CharP.cast_eq_zero_iff (ZMod q) q k).mp hz) hk.1)
    · intro h hh
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hh
      simp only [Finset.mem_filter, Finset.mem_range]
      refine ⟨ZMod.val_lt h, ?_⟩
      intro hv
      apply hh
      have hnc := ZMod.natCast_zmod_val h
      rw [hv, Nat.cast_zero] at hnc
      exact hnc.symm
    · intro k hk
      rw [Finset.mem_filter, Finset.mem_range] at hk
      exact ZMod.val_natCast_of_lt hk.1
    · intro h _; exact ZMod.natCast_zmod_val h
    · intro k hk
      rw [Finset.mem_filter, Finset.mem_range] at hk
      simp only [hFdef, hbcdef, ZMod.val_natCast_of_lt hk.1]
      congr 2
      refine bandDFT_periodic q hq B _ _ ?_
      have hval : (((w : ZMod q) * (k : ZMod q)).val : ℤ) = (w : ℤ) * (k : ℤ) % q := by
        rw [show (w : ZMod q) * (k : ZMod q) = ((w * k : ℕ) : ZMod q) from by push_cast; ring,
          ZMod.val_natCast, Int.natCast_mod, Nat.cast_mul]
      rw [hval]
      exact ⟨(w : ℤ) * (k : ℤ) / q, by rw [Int.emod_def]; ring⟩
  -- assemble the bridge equality
  have hbridge : (corrCount q B w : ℂ) - (B.card : ℂ) ^ 2 / q
      = (1 / q) * ∑ h ∈ (Finset.univ : Finset (ZMod q)).filter (fun h => h ≠ 0),
          bc h * (starRingEnd ℂ) (bc ((w : ZMod q) * h)) := by
    rw [completion_identity q hq B hB w, ← hFdef, hsplit, ← hoffeq]
    field_simp
    ring
  rw [hbridge, norm_mul]
  have h1q : ‖(1 / (q : ℂ))‖ = 1 / q := by rw [norm_div, norm_one, Complex.norm_natCast]
  rw [h1q]
  have hkps := offDiag_bandSum_le_closed (w : ZMod q) hw bc P hP hPmin hbc
  calc (1 / (q : ℝ)) * ‖∑ h ∈ (Finset.univ : Finset (ZMod q)).filter (fun h => h ≠ 0),
          bc h * (starRingEnd ℂ) (bc ((w : ZMod q) * h))‖
      ≤ (1 / (q : ℝ)) * (5 * (q : ℝ) ^ 2 * ((Nat.log 2 q : ℝ) + 1) ^ 2 / P) :=
        mul_le_mul_of_nonneg_left hkps (by positivity)
    _ = 5 * (q : ℝ) * ((Nat.log 2 q : ℝ) + 1) ^ 2 / P := by field_simp

open LonelyRunner.HyperbolaBox in
/-- **The band bridge.**  For the interval band `Icc lo hi`, `B̂(j)` is exactly B.2's interval exponential
sum (via `sum_Ico_eq_sum_range` reindex and `e_q` periodicity), so `‖B̂(j)‖ ≤ q/(2·cdist j)` for `j ≢ 0`
(`cdist_neg` reconciles the `e_q(−jr)` sign). -/
theorem norm_bandDFT_Icc_le (q : ℕ) (hq : 0 < q) (lo hi : ℕ) (j : ℤ) (hj : ((j : ℤ) : ZMod q) ≠ 0) :
    ‖bandDFT q (Finset.Icc lo hi) j‖ ≤ (q : ℝ) / (2 * (cdist ((j : ℤ) : ZMod q) : ℝ)) := by
  haveI : NeZero q := ⟨hq.ne'⟩
  have hnj : (((-j : ℤ)) : ZMod q) ≠ 0 := by rw [Int.cast_neg]; exact neg_ne_zero.mpr hj
  have hsum : bandDFT q (Finset.Icc lo hi) j
      = ∑ r ∈ Finset.range (hi + 1 - lo),
          Complex.exp ((2 * π * (((-j : ℤ) : ℝ) * ((lo : ℝ) + r) / q) : ℝ) * Complex.I) := by
    rw [bandDFT,
      show Finset.Icc lo hi = Finset.Ico lo (hi + 1) from by
        ext x; simp only [Finset.mem_Icc, Finset.mem_Ico]; omega,
      Finset.sum_Ico_eq_sum_range]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [eInt]; congr 1; push_cast; ring
  rw [hsum]
  calc ‖∑ r ∈ Finset.range (hi + 1 - lo),
          Complex.exp ((2 * π * (((-j : ℤ) : ℝ) * ((lo : ℝ) + r) / q) : ℝ) * Complex.I)‖
      ≤ (q : ℝ) / (2 * (cdist (((-j : ℤ)) : ZMod q) : ℝ)) := norm_bandCoeff_le q hq (-j) lo _ hnj
    _ = (q : ℝ) / (2 * (cdist ((j : ℤ) : ZMod q) : ℝ)) := by rw [Int.cast_neg, cdist_neg]

open LonelyRunner.HyperbolaBox in
/-- **LEM-022, the t2 hyperbola bound — for the actual band, unconditional.**  For the interval band
`B = Icc lo hi ⊆ [0,q)`, a unit twist `w`, and the ratio-lattice floor `P`, the multiplicative
pair-correlation obeys `‖C_w − b²/q‖ ≤ 5q(log₂q+1)²/P` — the whole Fourier-completion node of the OffLine
gate (THM-680), kernel-pure: the completion identity + kps's aggregation + opus's B.2, composed. -/
theorem completion_band {q : ℕ} [NeZero q] (lo hi : ℕ) (hhi : hi < q)
    (w : ℕ) (hw : IsUnit ((w : ZMod q)))
    (P : ℕ) (hP : 0 < P)
    (hPmin : ∀ h : ZMod q, h ≠ 0 → P ≤ cdist h * cdist ((w : ZMod q) * h)) :
    ‖(corrCount q (Finset.Icc lo hi) w : ℂ) - ((Finset.Icc lo hi).card : ℂ) ^ 2 / q‖
      ≤ 5 * (q : ℝ) * ((Nat.log 2 q : ℝ) + 1) ^ 2 / P := by
  have hq : 0 < q := Nat.pos_of_ne_zero (NeZero.ne q)
  have hsub : Finset.Icc lo hi ⊆ Finset.range q := by
    intro r hr; rw [Finset.mem_Icc] at hr; rw [Finset.mem_range]; omega
  refine completion_final (Finset.Icc lo hi) hsub w hw P hP hPmin ?_
  intro h hh
  have hval : (((h.val : ℤ) : ZMod q)) = h := by
    rw [Int.cast_natCast, ZMod.natCast_zmod_val]
  have hb := norm_bandDFT_Icc_le q hq lo hi (h.val : ℤ) (by rw [hval]; exact hh)
  rwa [hval] at hb

end FourierCompletion
end LonelyRunner




