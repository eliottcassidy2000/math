/-
  TournamentH7.LRCReachWitness — the reach side of hpartA, in usable form (kind-pasteur-2026-07-09-S99b).

  `hpartA` (the skeleton's Part A) is `0 < witnessG2 (shapeOf v) → 1/14 ≤ Mreach v`.  Its final step
  is: a positive covering witness produces a LONELY INSTANT — a time `t` at which every runner
  `v i · t` sits `≥ 1/14` from the origin (the observer) — and any lonely instant forces
  `Mreach v ≥ 1/14`.  `LRCMreachConcrete` proves the compactness/attainment core; `LRC14Assembly`
  has `Mreach_ge_of_witness` for `minReach v t` on `[0,1]`.  This file pins the CONCRETE TARGET the
  covering→reach bridge must hit, in the shape it naturally produces it:

    * `le_minReach_iff` — `c ≤ minReach v t ↔ ∀ i, c ≤ nearInt (v i · t)` (the observer clears every
      runner iff it clears the nearest);  a clean reusable unfolding of the `⨅`.
    * `Mreach_ge_of_lonely_instant` — `(∃ t, ∀ i, 1/14 ≤ nearInt (v i · t)) → 1/14 ≤ Mreach v`, at
      ANY real `t` (via the global-`sSup` form, no `[0,1]` bookkeeping).

  Together with `lonely_of_Mreach_ge` (LRCMreachConcrete) this is the exact reach target for hpartA.
-/
import Mathlib
import TournamentH7.LRCMreachConcrete

namespace LonelyRunner
namespace LRC14Concrete

/-- **The observer clears every runner iff it clears the nearest.**  `c ≤ minReach v t` exactly when
every residue `nearInt (v i · t)` is `≥ c` — the definitional unfolding of the finite infimum. -/
theorem le_minReach_iff (v : Fin 13 → ℤ) (t c : ℝ) :
    c ≤ minReach v t ↔ ∀ i, c ≤ nearInt ((v i : ℝ) * t) := by
  have hbdd : BddBelow (Set.range fun i => nearInt ((v i : ℝ) * t)) := by
    refine ⟨0, ?_⟩
    rintro _ ⟨i, rfl⟩
    exact nearInt_nonneg _
  simpa [minReach] using le_ciInf_iff hbdd

/-- **A lonely instant forces the reach.**  If at some time `t` every runner `v i · t` is `≥ 1/14`
from the origin, then `Mreach v ≥ 1/14`.  This is the reach conclusion of Part A, stated at an
arbitrary real `t` (the global supremum sees every time by `1`-periodicity). -/
theorem Mreach_ge_of_lonely_instant (v : Fin 13 → ℤ)
    (h : ∃ t : ℝ, ∀ i, (1 : ℝ) / 14 ≤ nearInt ((v i : ℝ) * t)) :
    (1 : ℝ) / 14 ≤ Mreach v := by
  obtain ⟨t, ht⟩ := h
  have hmin : (1 : ℝ) / 14 ≤ minReach v t := (le_minReach_iff v t _).mpr ht
  have hbddA : BddAbove (Set.range (minReach v)) :=
    ⟨1 / 2, by rintro _ ⟨s, rfl⟩; exact minReach_le_half v s⟩
  rw [Mreach_eq_global_sSup]
  exact le_trans hmin (le_csSup hbddA ⟨t, rfl⟩)

end LRC14Concrete
end LonelyRunner
