/-
  TournamentH7.LRCHlargeRoute  (opus-2026-07-03-S54)

  THE CASE-ROUTING SKELETON for the large-magnitude branch `hlarge` of the two-sided architecture
  (`lrc14_of_magnitude_split`, kps).  The renormalization-depth architecture (HYP-4041) routes every
  covering family into four complementary closers:
     (1) global ratio <= 13            => `spread13_lonely` (kps)              -- PROVED
     (2) scale gap, <= 6 far runners   => `lonely_of_simul_peel` (kps)        -- PROVED
     (3) scale gap, near-equal cluster => `scale_separation` (opus THM-608)   -- PROVED
     (3') resonant 13-comb cluster     => `scale_separation_phase` (opus)     -- PROVED
     (4) bounded magnitude             => bounded-denominator census          -- finite
  All the ENGINES are landed, kernel-pure.  What was missing is the ROUTING that dispatches a family to
  the right engine.  This file lands route (1) unconditionally and reduces `hlarge` to the single GAP
  obligation (ratio > 13 => a scale gap to peel), so the remaining surface is exactly the
  peel/tower dispatch on families that are NOT all-comparable.

  Route (1) is the base of the renormalization tower: it is the case with NO scale to separate.
-/
import TournamentH7.LRCSpread13
import TournamentH7.LRC14CertRoute

namespace LonelyRunner
namespace HlargeRoute

open LonelyRunner
open LonelyRunner.LRC14

/-- **Route (1): the all-comparable base case.**  If every pair of speeds is within a factor 13
(`∀ i j, |v i| ≤ 13·|v j|`, i.e. `max ≤ 13·min`), the family is `Lonely 14` by `spread13_lonely` at
`t = 1/(min+max)`; otherwise there is a scale GAP and loneliness is deferred to `hgap`.  This discharges
architecture route (1) and isolates routes (2)–(4) as the single `hgap` input. -/
theorem lonely14_of_ratio13_or_gap
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hgap : (¬ ∀ i j : Fin 13, |v i| ≤ 13 * |v j|) → ∃ t : ℝ, Lonely 14 v t) :
    ∃ t : ℝ, Lonely 14 v t := by
  by_cases hr : ∀ i j : Fin 13, |v i| ≤ 13 * |v j|
  · obtain ⟨imax, -, hmax⟩ :=
      Finset.exists_max_image (Finset.univ) (fun i => |v i|) ⟨0, Finset.mem_univ 0⟩
    obtain ⟨jmin, -, hmin⟩ :=
      Finset.exists_min_image (Finset.univ) (fun i => |v i|) ⟨0, Finset.mem_univ 0⟩
    exact ⟨_, spread13_lonely v (|v jmin|) (|v imax|)
      (abs_pos.mpr (hv jmin))
      (fun i => hmin i (Finset.mem_univ i))
      (fun i => hmax i (Finset.mem_univ i))
      (hr imax jmin)⟩
  · exact hgap hr

/-- **`hlarge` reduced to the GAP obligation.**  The large-magnitude branch of the two-sided architecture
follows from the single obligation: covering families that are NOT all-comparable (a scale gap, ratio > 13)
are lonely.  The all-comparable slice is discharged internally by route (1).  So the entire remaining content
of `hlarge` is the peel/tower dispatch on gap families — the structural decomposition. -/
theorem hlarge_of_gap (M : ℤ)
    (hgap : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v → (∃ i, M < |v i|) →
      (¬ ∀ i j : Fin 13, |v i| ≤ 13 * |v j|) → ∃ t : ℝ, Lonely 14 v t) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v → (∃ i, M < |v i|) →
      ∃ t : ℝ, Lonely 14 v t := by
  intro v hv hcov hfar
  exact lonely14_of_ratio13_or_gap v hv (fun hnr => hgap v hv hcov hfar hnr)

end HlargeRoute
end LonelyRunner
