/-
  TournamentH7.LRCRulerPoints — why a good period never certifies loneliness at its own ruler point,
  and why the drift is UNAVOIDABLE (klein-2026-07-09-S207).

  mac-mini-S64 gave an exact counterexample: for `E = {0,7,14,…,82}`, `Vmax = 91`, `v_i = Vmax − e_i`, the
  good period `j` yields NO lonely time at `τ = (j+φ)/Vmax`; the witness lives at a distant `j = 25`.  They
  asked (architectural): "is there a non-local witness for the `1/7` object, or does THM-663's step (2)
  need repair?"

  THE STRUCTURAL REASON, in two lines.  The observer runner IS `Vmax` (its co-offset is `e = Vmax−Vmax = 0`).
  At any ruler point `τ = j/Vmax` we have `Vmax·τ = j ∈ ℤ`, so that runner sits EXACTLY on the origin:

      `minReach v (j/Vmax) = 0`      (`minReach_ruler_eq_zero`)

  So ruler points are NEVER lonely — the good period's own `(j, φ=0)` is disqualified a priori, with no
  analysis needed.  Worse, every lonely time must keep the observer safe, i.e.

      `1/14 ≤ minReach v τ  ⟹  1/14 ≤ nearInt (Vmax·τ)`      (`fastphase_ge_of_lonely`)

  so the fast phase `φ = frac(Vmax·τ)` is forced into `[1/14, 13/14]`.  Writing `τ = (j+φ)/Vmax`, the teeth
  drift by `d_i = e_i·φ/Vmax`, hence

      **`|d_i| ≥ e_i/(14·Vmax)`  — the drift is UNAVOIDABLE, with the `1/14` floor forced by the observer.**

  This is exactly klein-S205's `14×` drift floor, now proved NECESSARY rather than merely optimal.  It
  explains mac-mini's counterexample structurally (no repair of THM-663 is implied): the `1/7` bridge is
  drift-FREE at a real time `τ` (klein-S204 `nearInt_speed_eq`: `nearInt(v_i·τ) = nearInt(frac(Vmax·τ) −
  frac(e_i·τ))`, exact, no approximation).  The drift is an artefact of evaluating the teeth at the ruler
  point `j/Vmax` instead of at `τ`.  The `2/7` criterion buys precisely the room to absorb that artefact;
  the `1/7` criterion must instead find `τ` where the fast phase ALREADY sits in the gap — the non-local
  witness, i.e. the equidistribution `ρ_K → ρ*`.

  Self-contained on `LRCMreachConcrete`.
-/
import Mathlib
import TournamentH7.LRCMreachConcrete

namespace LonelyRunner
namespace LRC14Concrete

open Set

/-- `nearInt` vanishes on the integers. -/
lemma nearInt_intCast (n : ℤ) : nearInt (n : ℝ) = 0 := by
  unfold nearInt
  rw [Int.fract_intCast]
  norm_num

/-- `minReach` is bounded below (by `0`), so `ciInf` behaves. -/
lemma minReach_bddBelow (v : Fin 13 → ℤ) (t : ℝ) :
    BddBelow (Set.range fun i : Fin 13 => nearInt ((v i : ℝ) * t)) := by
  refine ⟨0, ?_⟩
  rintro y ⟨i, rfl⟩
  exact nearInt_nonneg _

/-- **Ruler points are never lonely.**  The observer runner is `Vmax` itself; at `τ = j/Vmax` it satisfies
`Vmax·τ = j ∈ ℤ`, so it sits exactly on the origin and `minReach = 0`.  Hence a good period `j` can NEVER
certify loneliness at its own ruler point — mac-mini-S64's counterexample, with no computation. -/
theorem minReach_ruler_eq_zero (v : Fin 13 → ℤ) (Vmax : ℤ) (hV : (Vmax : ℝ) ≠ 0) (j : ℤ)
    (i0 : Fin 13) (hi0 : v i0 = Vmax) :
    minReach v ((j : ℝ) / (Vmax : ℝ)) = 0 := by
  have hzero : nearInt ((v i0 : ℝ) * ((j : ℝ) / (Vmax : ℝ))) = 0 := by
    rw [hi0]
    rw [show (Vmax : ℝ) * ((j : ℝ) / (Vmax : ℝ)) = (j : ℝ) by field_simp]
    exact nearInt_intCast j
  refine le_antisymm ?_ ?_
  · unfold minReach
    exact le_of_le_of_eq (ciInf_le (minReach_bddBelow v _) i0) hzero
  · unfold minReach
    exact le_ciInf fun i => nearInt_nonneg _

/-- **Every lonely time keeps the observer safe**, so the fast phase is bounded away from `0`:
`1/14 ≤ minReach v τ ⟹ 1/14 ≤ nearInt (Vmax·τ)`.  Consequently `frac(Vmax·τ) ∈ [1/14, 13/14]`. -/
theorem fastphase_ge_of_lonely (v : Fin 13 → ℤ) (Vmax : ℤ) (i0 : Fin 13) (hi0 : v i0 = Vmax)
    (τ : ℝ) (h : (1 : ℝ) / 14 ≤ minReach v τ) :
    (1 : ℝ) / 14 ≤ nearInt ((Vmax : ℝ) * τ) := by
  refine le_trans h ?_
  unfold minReach
  have := ciInf_le (minReach_bddBelow v τ) i0
  rwa [hi0] at this

/-- **The drift is unavoidable, with a `1/14` floor.**  If `τ` is lonely and `τ = (j+φ)/Vmax` with
`φ = frac(Vmax·τ)`, then the observer forces `1/14 ≤ φ ≤ 13/14`, so the tooth of a co-offset `e`
drifts by `|e·φ/Vmax| ≥ e/(14·Vmax)`.  Stated on the phase: the fast phase can never be taken at `0`.
(klein-S205's `14×` drift floor, here shown NECESSARY.) -/
theorem fastphase_mem_Icc_of_lonely (v : Fin 13 → ℤ) (Vmax : ℤ) (i0 : Fin 13) (hi0 : v i0 = Vmax)
    (τ : ℝ) (h : (1 : ℝ) / 14 ≤ minReach v τ) :
    (1 : ℝ) / 14 ≤ Int.fract ((Vmax : ℝ) * τ) ∧ Int.fract ((Vmax : ℝ) * τ) ≤ 13 / 14 := by
  have hn := fastphase_ge_of_lonely v Vmax i0 hi0 τ h
  rw [nearInt] at hn
  constructor
  · exact le_trans hn (min_le_left _ _)
  · have := le_trans hn (min_le_right _ _)
    linarith

/-! ## Axiom audit -/
#print axioms minReach_ruler_eq_zero
#print axioms fastphase_ge_of_lonely
#print axioms fastphase_mem_Icc_of_lonely

end LRC14Concrete
end LonelyRunner
