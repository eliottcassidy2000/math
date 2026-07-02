/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-02-S47)
-/
import TournamentH7.LonelyRunner

/-!
# The kernel witness gate: census rows without `native_decide`

**The mathlib-acceptance upgrade.**  mathlib does not accept `native_decide` (it adds the
compiler to the trusted base); 40 corpus files currently carry it.  The obstruction was never
the arithmetic — kernel `decide` on ℤ/ℕ literals is GMP-accelerated — but ℚ normalization
(`Nat.gcd` well-founded recursion) inside `Rat` operations.

This gate re-expresses the rational-witness check in PURE INTEGER arithmetic: for witness
`num/den` and speed `s`, with `r := (s·num) % den`,

    den ≤ 14 · min r (den − r)

is a literal ℤ comparison, kernel-decidable.  `lonely_of_kernelWitness` (kernel-pure) turns a
passing 13-tuple into `Lonely 14 v (num/den)`.  Census rows checked this way carry ONLY the
standard axioms — no `Lean.ofReduceBool` anywhere.  Migration: mechanical re-emit of the
witness packs with `by decide`; sample migrated rows from `LRCWindowPack1` below.
-/

namespace LonelyRunner
namespace KernelGate

/-- One speed clears the `1/14` band at witness `num/den`, in pure integer form. -/
def speedOK (s num : ℤ) (den : ℕ) : Prop :=
  (den : ℤ) ≤ 14 * min ((s * num) % (den : ℤ)) ((den : ℤ) - (s * num) % (den : ℤ))

instance (s num : ℤ) (den : ℕ) : Decidable (speedOK s num den) := by
  unfold speedOK; infer_instance

/-- The integer distance step: `|p − m·den| ≥ min (p % den) (den − p % den)` for every `m`. -/
theorem int_dist_ge {p : ℤ} {den : ℕ} (hden : 0 < den) (m : ℤ) :
    min (p % (den : ℤ)) ((den : ℤ) - p % (den : ℤ)) ≤ |p - m * den| := by
  have hdenZ : (0 : ℤ) < (den : ℤ) := by exact_mod_cast hden
  have hr0 : 0 ≤ p % (den : ℤ) := Int.emod_nonneg p hdenZ.ne'
  have hrlt : p % (den : ℤ) < (den : ℤ) := Int.emod_lt_of_pos p hdenZ
  set r : ℤ := p % (den : ℤ) with hr
  set q : ℤ := p / (den : ℤ) with hq
  have hpq : p = (den : ℤ) * q + r := by rw [hq, hr]; exact (Int.ediv_add_emod p _).symm
  have hval : p - m * den = (den : ℤ) * (q - m) + r := by rw [hpq]; ring
  rcases eq_or_ne q m with rfl | hne
  · simp only [hval, sub_self, mul_zero, zero_add]
    rw [abs_of_nonneg hr0]
    exact min_le_left _ _
  · have h1 : 1 ≤ |q - m| := by
      rcases lt_or_gt_of_ne (sub_ne_zero.mpr hne) with h | h
      · rw [abs_of_neg h]; omega
      · rw [abs_of_pos h]; omega
    have h2 : (den : ℤ) ≤ |(den : ℤ) * (q - m)| := by
      rw [abs_mul, abs_of_nonneg hdenZ.le]
      nlinarith
    have h3 : |(den : ℤ) * (q - m)| - r ≤ |(den : ℤ) * (q - m) + r| := by
      rcases abs_cases ((den : ℤ) * (q - m) + r) with ⟨he, _⟩ | ⟨he, _⟩ <;>
        rcases abs_cases ((den : ℤ) * (q - m)) with ⟨he2, _⟩ | ⟨he2, _⟩ <;> omega
    calc min r ((den : ℤ) - r) ≤ (den : ℤ) - r := min_le_right _ _
      _ ≤ |(den : ℤ) * (q - m)| - r := by linarith
      _ ≤ |(den : ℤ) * (q - m) + r| := h3
      _ = |p - m * den| := by rw [hval]

/-- **The kernel gate (soundness, kernel-pure)**: a 13-tuple passing the integer check at
`num/den` is `1/14`-lonely at the rational time.  Rows through this gate carry only the
standard axioms. -/
theorem lonely_of_kernelWitness {v : Fin 13 → ℤ} {num : ℤ} {den : ℕ} (hden : 0 < den)
    (h : ∀ i, speedOK (v i) num den) :
    Lonely 14 v ((((num : ℚ) / (den : ℚ)) : ℚ) : ℝ) := by
  intro i m
  have hdenZ : (0 : ℤ) < (den : ℤ) := by exact_mod_cast hden
  have hdenR : (0 : ℝ) < (den : ℝ) := by exact_mod_cast hden
  have hkey := int_dist_ge (p := v i * num) hden m
  have hOK := h i
  unfold speedOK at hOK
  have hval : (v i : ℝ) * (((num : ℚ) / (den : ℚ) : ℚ) : ℝ) - (m : ℝ)
      = ((v i * num - m * den : ℤ) : ℝ) / (den : ℝ) := by
    push_cast
    field_simp
  rw [hval, abs_div, abs_of_pos hdenR, le_div_iff₀ hdenR]
  have hZ : (den : ℤ) ≤ 14 * |v i * num - m * den| := by
    calc (den : ℤ)
        ≤ 14 * min ((v i * num) % (den : ℤ)) ((den : ℤ) - (v i * num) % (den : ℤ)) := hOK
      _ ≤ 14 * |v i * num - m * den| := by nlinarith [hkey]
  have hZR : ((den : ℤ) : ℝ) ≤ ((14 * |v i * num - m * den| : ℤ) : ℝ) := by exact_mod_cast hZ
  push_cast at hZR
  push_cast
  linarith

/-! ### Migrated sample rows (kernel `decide` — NO `native_decide`) -/

/-- Pack row `{1,…,13}` at witness `1/14`, kernel-checked. -/
theorem kernelRow_AP : Lonely 14 (![1,2,3,4,5,6,7,8,9,10,11,12,13] : Fin 13 → ℤ)
    ((((1 : ℚ) / (14 : ℚ)) : ℚ) : ℝ) :=
  lonely_of_kernelWitness (by norm_num) (by decide)

/-- Pack row `{1,…,12,14}` at its witness, kernel-checked. -/
theorem kernelRow_2 : Lonely 14 (![1,2,3,4,5,6,7,8,9,10,11,12,14] : Fin 13 → ℤ)
    ((((1 : ℚ) / (13 : ℚ)) : ℚ) : ℝ) :=
  lonely_of_kernelWitness (by norm_num) (by decide)

end KernelGate
end LonelyRunner
