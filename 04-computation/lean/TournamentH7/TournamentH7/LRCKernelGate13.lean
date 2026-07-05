/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-05-S77)
-/
import TournamentH7.LRCKernelGate

/-!
# The strict kernel gate at band 1/13 (hdich lift-rigidity rows)

The 13-band twin of `LRCKernelGate` (S47), STRICT: a 12-tuple passing the integer check
`den < 13 · min r (den − r)` at witness `num/den` satisfies `1/13 < |v·t − m|` for every
runner and integer — "strictly loose", exactly what the lift-rigidity rows must certify.
Kernel `decide` throughout; standard axioms only.
-/

namespace LonelyRunner
namespace KernelGate13

/-- One speed clears the `1/13` band STRICTLY at witness `num/den`, integer form. -/
def speedOK13 (s num : ℤ) (den : ℕ) : Prop :=
  (den : ℤ) < 13 * min ((s * num) % (den : ℤ)) ((den : ℤ) - (s * num) % (den : ℤ))

instance (s num : ℤ) (den : ℕ) : Decidable (speedOK13 s num den) := by
  unfold speedOK13; infer_instance

/-- Strict loneliness at band `1/13`. -/
def StrictLonely13 {ι : Type*} (v : ι → ℤ) (t : ℝ) : Prop :=
  ∀ i, ∀ m : ℤ, (1 : ℝ) / 13 < |(v i : ℝ) * t - m|

/-- **The strict 13-gate (soundness, kernel-pure)**: a 12-tuple passing the strict integer
check at `num/den` is strictly `1/13`-lonely at the rational time. -/
theorem strictLonely13_of_kernelWitness {v : Fin 12 → ℤ} {num : ℤ} {den : ℕ} (hden : 0 < den)
    (h : ∀ i, speedOK13 (v i) num den) :
    ∃ t : ℝ, StrictLonely13 v t := by
  refine ⟨(((num : ℚ) / (den : ℚ) : ℚ) : ℝ), fun i m => ?_⟩
  have hdenZ : (0 : ℤ) < (den : ℤ) := by exact_mod_cast hden
  have hdenR : (0 : ℝ) < (den : ℝ) := by exact_mod_cast hden
  have hkey := KernelGate.int_dist_ge (p := v i * num) hden m
  have hOK := h i
  unfold speedOK13 at hOK
  have hval : (v i : ℝ) * (((num : ℚ) / (den : ℚ) : ℚ) : ℝ) - (m : ℝ)
      = ((v i * num - m * den : ℤ) : ℝ) / (den : ℝ) := by
    push_cast
    field_simp
  rw [hval, abs_div, abs_of_pos hdenR, lt_div_iff₀ hdenR]
  have hZ : (den : ℤ) < 13 * |v i * num - m * den| := by
    calc (den : ℤ)
        < 13 * min ((v i * num) % (den : ℤ)) ((den : ℤ) - (v i * num) % (den : ℤ)) := hOK
      _ ≤ 13 * |v i * num - m * den| := by nlinarith [hkey]
  have hZR : ((den : ℤ) : ℝ) < ((13 * |v i * num - m * den| : ℤ) : ℝ) := by exact_mod_cast hZ
  push_cast at hZR
  push_cast
  linarith

end KernelGate13
end LonelyRunner
