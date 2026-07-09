/-
  TournamentH7.LRCGoodPeriodReach — wiring the good-period dichotomy to the reach core
  (klein-2026-07-09-S203).

  The good-period leg (LEM-012 ∪ LEM-013, assembled by kps-S99 `LRCGoodPeriodDispatch`) produces a
  concrete `HasGoodPeriod E Vmax` (klein's `LRCGoodPeriodMaxgap`: `∃ j, 7·maxCircGap > Vmax`).  kps-S31
  `LRCGapReach` turns a free real gap `> 1/7` into a `nearInt` loneliness margin `> 1/14`.  THIS FILE
  supplies the missing LINK between them: a good period's residue gap IS such a free real gap.

  Chain (bridge to `M(S) ≥ 1/14`):
    HasGoodPeriod E Vmax                         [dichotomy: LEM-012/013, kps-S99]
      → ∃ free real interval of teeth, len > 1/7  [THIS FILE: `reachMargin_of_residueGap`]
      → ∃ φ, nearInt(φ − tooth) > 1/14 ∀ teeth    [GapReach.exists_nearInt_margin_of_gap]
      → ∃ real τ, minReach v τ ≥ 1/14             [RULER EMBEDDING — deferred Part-A, the single
                                                    remaining hypothesis, shared with the density route]
      → 1/14 ≤ Mreach v                           [Assembly.Mreach_ge_of_witness]

  Everything except the ruler embedding is proven here / cited sorry-free.  Self-contained.
-/
import Mathlib
import TournamentH7.LRCGoodPeriodMaxgap
import TournamentH7.LRCGapReach

namespace LRC14

open LonelyRunner

/-- The **teeth** of a cluster at ruler dilation `j`: the residues `e·j mod Vmax`, normalized to
`[0,1)`, as a `Finset ℝ`. -/
noncomputable def teeth (E : List ℕ) (Vmax j : ℕ) : Finset ℝ :=
  (E.map (fun e => ((e * j % Vmax : ℕ) : ℝ) / (Vmax : ℝ))).toFinset

/-- **The reach-margin from a residue gap.**  Suppose the normalized open interval
`(a, a + g)` with `g > 1/7` contains no integer translate of any tooth.  Then a single fast phase `φ`
clears every tooth by more than `1/14` — the loneliness margin.  (A direct repackaging of kps-S31's
`GapReach.exists_nearInt_margin_of_gap` at the concrete teeth of a good period.) -/
theorem reachMargin_of_residueGap (E : List ℕ) (Vmax j : ℕ) (a g : ℝ) (hg : 1 / 7 < g)
    (hfree : ∀ c ∈ teeth E Vmax j, ∀ n : ℤ, (c + (n : ℝ)) ∉ Set.Ioo a (a + g)) :
    ∃ φ : ℝ, ∀ c ∈ teeth E Vmax j, 1 / 14 < LRC14Concrete.nearInt (φ - c) :=
  GapReach.exists_nearInt_margin_of_gap (teeth E Vmax j) a g hg hfree

/-- **The good-period → reach bridge, modulo the ruler embedding.**  Given a good period (hence,
via the free-gap link `hlink`, a `1/14` margin from every tooth) and the RULER EMBEDDING `hembed`
(the deferred Part-A content: a fast phase clearing all teeth is realized at a real time `τ` with
`minReach v τ ≥ 1/14`), the reach floor `1/14 ≤ Mreach v` holds.  The ONE remaining hypothesis
`hembed` is exactly THM-527 Part A's `Vmax`-ruler embedding, shared with the density-floor route. -/
theorem mreach_ge_of_goodPeriod
    (E : List ℕ) (Vmax : ℕ) (v : Fin 13 → ℤ)
    (hgp : HasGoodPeriod E Vmax)
    (hlink : HasGoodPeriod E Vmax →
      ∃ (j : ℕ) (a g : ℝ), 1 / 7 < g ∧
        ∀ c ∈ teeth E Vmax j, ∀ n : ℤ, (c + (n : ℝ)) ∉ Set.Ioo a (a + g))
    (hembed : (∃ (j : ℕ) (φ : ℝ), ∀ c ∈ teeth E Vmax j,
        1 / 14 < LRC14Concrete.nearInt (φ - c)) →
      ∃ τ : ℝ, τ ∈ Set.Icc (0 : ℝ) 1 ∧ (1 : ℝ) / 14 ≤ LRC14Concrete.minReach v τ) :
    (1 : ℝ) / 14 ≤ LRC14Concrete.Mreach v := by
  obtain ⟨j, a, g, hg, hfree⟩ := hlink hgp
  have hmargin : ∃ (j : ℕ) (φ : ℝ), ∀ c ∈ teeth E Vmax j,
      1 / 14 < LRC14Concrete.nearInt (φ - c) := by
    obtain ⟨φ, hφ⟩ := reachMargin_of_residueGap E Vmax j a g hg hfree
    exact ⟨j, φ, hφ⟩
  obtain ⟨τ, hτmem, hτ⟩ := hembed hmargin
  -- inline `Mreach_ge_of_witness`: a single witness time bounds the sup `Mreach = sSup (minReach '' [0,1])`
  refine le_trans hτ (le_csSup ?_ (Set.mem_image_of_mem _ hτmem))
  refine ⟨1 / 2, ?_⟩
  rintro y ⟨s, -, rfl⟩
  exact LRC14Concrete.minReach_le_half v s

/-! ## Axiom audit -/
#print axioms reachMargin_of_residueGap
#print axioms mreach_ge_of_goodPeriod

end LRC14
