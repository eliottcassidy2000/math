/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-10-S205)
-/
import Mathlib
import TournamentH7.LRC14WindowWiring

/-!
# The primitivity peel — LRC(14) reduces to primitive families

**Why this matters (opus-S205).** The `ResidualObligation` domain, as stated, ADMITS NON-PRIMITIVE
DILATES. Machine-checked witness: `v = 2·[1,2,3,4,5,6,7,8,9,11,12,13,20] =
[2,4,6,8,10,12,14,16,18,22,24,26,40]` satisfies every residual clause (covering, scale gap `> 13`,
compressed, distinct `|vᵢ|`, some `|vᵢ| ≥ 23`, divisor-closed, no nontrivial common residue) with
`gcd = 2`, while its core has `Vmax = 20 ≤ 22` (already handled by the window census).

Since `α ↦ c·α` is measure-preserving on the circle, the safe-set measure satisfies `μ(c·w) = μ(w)`
EXACTLY. Cores range over the window census, which contains near-APs with `μ` arbitrarily small (the tight
AP has `μ = 0`). Hence

> `inf μ = 0` over the residual as stated — **no uniform measure floor `μ₀ > 0` exists**,

so (via klein's THM-685 Kronecker transfer, where a floor `μ ≥ μ₀` buys liveness at `q > Σv/μ₀`) any proof
routing through a UNIFORM measure floor is unreachable. The obligation is still *true*; the target was
ill-posed.

**The fix, here.** Dilates are trivial: `c·w` is lonely at `t/c` whenever `w` is lonely at `t`
(`lonely_scale`). So LRC(14) reduces to PRIMITIVE families, and the residual may be restated with
`tupleGcd v = 1`, after which the dilates are gone and a uniform floor becomes well-posed (empirically
`min μ ≈ 0.0094` over the primitive residual, vs iid `(6/7)^13 ≈ 0.1348`).

`lrc14_of_primitive : (∀ primitive w, lonely) → LRC14Statement`.

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace LRC14
namespace PrimitivePeel

open Finset

/-- **Dilates are free.** If the core `w` has a lonely instant, so does `c·w` (at `t/c`). -/
theorem lonely_of_dilate (w : Fin 13 → ℤ) (c : ℤ) (hc : c ≠ 0)
    (h : ∃ t : ℝ, Lonely 14 w t) : ∃ t : ℝ, Lonely 14 (fun i => c * w i) t := by
  obtain ⟨t, ht⟩ := h
  exact ⟨t / (c : ℝ), (lonely_scale 14 w t c hc).mpr ht⟩

/-- **The primitivity peel.** LRC(14) follows from LRC(14) on PRIMITIVE families (`tupleGcd = 1`). -/
theorem lrc14_of_primitive
    (hprim : ∀ w : Fin 13 → ℤ, (∀ i, w i ≠ 0) → tupleGcd w = 1 → ∃ t : ℝ, Lonely 14 w t) :
    LRC14Statement := by
  intro v hv
  have hcpos : 0 < tupleGcd v := tupleGcd_pos v hv
  have hcne : ((tupleGcd v : ℕ) : ℤ) ≠ 0 := by exact_mod_cast hcpos.ne'
  set c : ℤ := ((tupleGcd v : ℕ) : ℤ) with hcdef
  set w : Fin 13 → ℤ := fun i => v i / c with hwdef
  have hvi : ∀ i, v i = c * w i := fun i => (Int.mul_ediv_cancel' (tupleGcd_dvd v i)).symm
  have hwne : ∀ i, w i ≠ 0 := by
    intro i hwi
    exact hv i (by rw [hvi i, hwi, mul_zero])
  have habs : ∀ i, (v i).natAbs = tupleGcd v * (w i).natAbs := by
    intro i
    rw [hvi i, Int.natAbs_mul, hcdef, Int.natAbs_natCast]
  have hwg : tupleGcd w = 1 := by
    have h1 : tupleGcd v = tupleGcd v * tupleGcd w := by
      conv_lhs => rw [tupleGcd]
      rw [show (fun i => (v i).natAbs) = (fun i => tupleGcd v * (w i).natAbs) from funext habs]
      rw [Finset.gcd_mul_left]
      simp [tupleGcd]
    have h2 : tupleGcd v * 1 = tupleGcd v * tupleGcd w := by simpa using h1
    exact (Nat.eq_of_mul_eq_mul_left hcpos h2).symm
  obtain ⟨t, ht⟩ := hprim w hwne hwg
  refine ⟨t / (c : ℝ), ?_⟩
  have hL : Lonely 14 (fun i => c * w i) (t / (c : ℝ)) := (lonely_scale 14 w t c hcne).mpr ht
  have hveq : v = fun i => c * w i := funext hvi
  rw [hveq]; exact hL

end PrimitivePeel
end LRC14
end LonelyRunner
