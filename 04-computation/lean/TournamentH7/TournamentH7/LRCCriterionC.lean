/-
  TournamentH7.LRCCriterionC — the algebraic core of THM-527 Part A (klein-2026-07-09-S204).

  Part A (`thm527_partA_density_pos_implies_reach`) is the deep node: `0 < ρ*(shapeOf v) ⟹ Mreach v ≥
  1/14`, the slow-fast change of variables + criterion C + equidistribution `ρ_K → ρ*`.  This file
  formalizes **criterion C's algebraic heart** — the reason a fast phase clearing the cluster teeth
  produces loneliness — and thereby REDUCES Part A to the pure realization/equidistribution step.

  THE IDENTITY.  With the co-offset convention `e_i = Vmax − v_i` (THM-527: `e = 0` for the observer
  `Vmax`), the runner position is the fast phase minus the tooth, modulo 1:

      ‖v_i·τ‖ = nearInt(v_i·τ) = nearInt( frac(Vmax·τ) − frac(e_i·τ) ).

  So if at some real time `τ` the fast phase `φ = frac(Vmax·τ)` clears every tooth `frac(e_i·τ)` by
  `> 1/14`, then EVERY runner is `> 1/14` from the origin — `minReach v τ ≥ 1/14`, hence `Mreach ≥ 1/14`.
  This is exactly kps-S31 `GapReach`'s `nearInt(φ − c) > 1/14` clearance, now identified with the
  concrete `minReach`.  The one thing left for Part A is the REALIZATION: that a good-period gap is
  achieved by `φ = frac(Vmax·τ)` at a `τ` whose teeth are the gap's cluster — the equidistribution
  `ρ_K → ρ*` (with the `O(1/Vmax)` finite-ruler correction).  Self-contained on `LRCMreachConcrete`.
-/
import Mathlib
import TournamentH7.LRCMreachConcrete

namespace LonelyRunner
namespace LRC14Concrete

open Set

/-- `nearInt` is `ℤ`-periodic (depends only on the fractional part). -/
lemma nearInt_add_int (x : ℝ) (n : ℤ) : nearInt (x + (n : ℝ)) = nearInt x := by
  unfold nearInt; rw [Int.fract_add_intCast]

/-- **Criterion-C identity.**  For the co-offset `e = Vmax − v`, the runner phase `v·τ` equals — modulo
`1`, hence under `nearInt` — the fast phase `frac(Vmax·τ)` minus the tooth `frac(e·τ)`. -/
theorem nearInt_speed_eq (v Vmax e : ℤ) (τ : ℝ) (he : e = Vmax - v) :
    nearInt ((v : ℝ) * τ) = nearInt (Int.fract ((Vmax : ℝ) * τ) - Int.fract ((e : ℝ) * τ)) := by
  have h1 : (v : ℝ) * τ = (Vmax : ℝ) * τ - (e : ℝ) * τ := by rw [he]; push_cast; ring
  have h2 : (Vmax : ℝ) * τ - (e : ℝ) * τ
      = (Int.fract ((Vmax : ℝ) * τ) - Int.fract ((e : ℝ) * τ))
        + (((⌊(Vmax : ℝ) * τ⌋ - ⌊(e : ℝ) * τ⌋ : ℤ)) : ℝ) := by
    simp only [Int.fract]; push_cast; ring
  rw [h1, h2, nearInt_add_int]

/-- **Criterion C (the reach core).**  If the co-offsets are `e_i = Vmax − v_i`, and at time `τ` the
fast phase `frac(Vmax·τ)` clears every tooth `frac(e_i·τ)` by `≥ 1/14`, then every runner is `≥ 1/14`
from the origin: `minReach v τ ≥ 1/14`. -/
theorem minReach_ge_of_fastphase_clears (v : Fin 13 → ℤ) (Vmax : ℤ) (e : Fin 13 → ℤ) (τ : ℝ)
    (he : ∀ i, e i = Vmax - v i)
    (hclear : ∀ i, (1 : ℝ) / 14 ≤
      nearInt (Int.fract ((Vmax : ℝ) * τ) - Int.fract ((e i : ℝ) * τ))) :
    (1 : ℝ) / 14 ≤ minReach v τ := by
  unfold minReach
  apply le_ciInf
  intro i
  rw [nearInt_speed_eq (v i) Vmax (e i) τ (he i)]
  exact hclear i

/-- **Criterion C ⟹ the reach floor `Mreach ≥ 1/14`.**  A witness time `τ ∈ [0,1]` at which the fast
phase clears all teeth bounds the global `Mreach` from below.  This is THM-527 Part A **reduced to the
realization**: supply a `τ` whose fast phase lands in the cluster's good-period gap (the equidistribution
`ρ_K → ρ*`), and the reach `M(S) ≥ 1/14` follows. -/
theorem Mreach_ge_of_fastphase_clears (v : Fin 13 → ℤ) (Vmax : ℤ) (e : Fin 13 → ℤ) (τ : ℝ)
    (hτ : τ ∈ Icc (0 : ℝ) 1) (he : ∀ i, e i = Vmax - v i)
    (hclear : ∀ i, (1 : ℝ) / 14 ≤
      nearInt (Int.fract ((Vmax : ℝ) * τ) - Int.fract ((e i : ℝ) * τ))) :
    (1 : ℝ) / 14 ≤ Mreach v := by
  refine le_trans (minReach_ge_of_fastphase_clears v Vmax e τ he hclear)
    (le_csSup ?_ (Set.mem_image_of_mem _ hτ))
  refine ⟨1 / 2, ?_⟩
  rintro y ⟨s, -, rfl⟩
  exact minReach_le_half v s

/-! ## Axiom audit -/
#print axioms nearInt_speed_eq
#print axioms minReach_ge_of_fastphase_clears
#print axioms Mreach_ge_of_fastphase_clears

end LRC14Concrete
end LonelyRunner
