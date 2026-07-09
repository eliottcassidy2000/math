/-
  TournamentH7.LRCTeethGap — building blocks for klein-S203's `hlink` (kind-pasteur-2026-07-09-S101).

  `LRCGoodPeriodReach.mreach_ge_of_goodPeriod` leaves two links: `hlink` (a good period gives a free
  residue gap `> 1/7`) and `hembed` (the ruler embedding, deferred Part-A).  This file discharges the
  GEOMETRIC core of `hlink`'s conclusion: the "no integer translate of a tooth lands in the gap"
  requirement collapses to "no tooth lands in the gap" whenever the gap is a NON-WRAPPING subinterval
  `(a, a+g) ⊆ [0,1]` — because every tooth lies in `[0,1)`, so the only integer translate that could
  reach a subinterval of `(0,1)` is the trivial one (`n = 0`).

  Provides:
    * `teeth_subset_Ico` — every tooth lies in `[0,1)` (residue `< Vmax`);
    * `free_translate_of_free_subInterval` — the integer-translate reduction for a non-wrapping gap.

  Together these reduce `hlink` (on the non-wrapping branch) to producing a tooth-free subinterval of
  length `> 1/7` from `maxCircGap` — the remaining piece being the `mergeSort` argmax extraction.
  Self-contained on `LRCGoodPeriodReach`.
-/
import Mathlib
import TournamentH7.LRCGoodPeriodReach

namespace LRC14

open LonelyRunner

/-- **Every tooth lies in `[0,1)`.**  The teeth are residues `e·j mod Vmax` divided by `Vmax`, so
they sit in `[0,1)` whenever `Vmax > 0`. -/
theorem teeth_subset_Ico (E : List ℕ) (Vmax j : ℕ) (hV : 0 < Vmax) :
    ∀ c ∈ teeth E Vmax j, c ∈ Set.Ico (0 : ℝ) 1 := by
  intro c hc
  simp only [teeth, List.mem_toFinset, List.mem_map] at hc
  obtain ⟨e, -, rfl⟩ := hc
  have hlt : e * j % Vmax < Vmax := Nat.mod_lt _ hV
  have hVpos : (0 : ℝ) < (Vmax : ℝ) := by exact_mod_cast hV
  constructor
  · positivity
  · rw [div_lt_one hVpos]; exact_mod_cast hlt

/-- **Integer-translate reduction (non-wrapping gap).**  If every point of a finite set `S` lies in
`[0,1)`, and the open interval `(a, a+g)` is a subinterval of `[0,1]` (`0 ≤ a`, `a+g ≤ 1`) that
contains no point of `S`, then it contains no INTEGER TRANSLATE of any point of `S` either: for every
`c ∈ S` and `n : ℤ`, `c + n ∉ (a, a+g)`.  (Any translate reaching `(a,a+g) ⊆ (0,1)` forces `n = 0`.) -/
theorem free_translate_of_free_subInterval (S : Finset ℝ) (a g : ℝ)
    (hSIco : ∀ c ∈ S, c ∈ Set.Ico (0 : ℝ) 1)
    (ha : 0 ≤ a) (hag : a + g ≤ 1)
    (hfree : ∀ c ∈ S, c ∉ Set.Ioo a (a + g)) :
    ∀ c ∈ S, ∀ n : ℤ, (c + (n : ℝ)) ∉ Set.Ioo a (a + g) := by
  intro c hc n hmem
  obtain ⟨hlo, hhi⟩ := hmem
  obtain ⟨hc0, hc1⟩ := hSIco c hc
  -- `c + n ∈ (a, a+g) ⊆ (0,1)` and `c ∈ [0,1)` force `n = 0`.
  have hpos : (0 : ℝ) < c + (n : ℝ) := lt_of_le_of_lt ha hlo
  have hlt1 : c + (n : ℝ) < 1 := lt_of_lt_of_le hhi hag
  have hnlt : (n : ℝ) < 1 := by linarith
  have hngt : (-1 : ℝ) < (n : ℝ) := by linarith
  have hn0 : n = 0 := by
    have h1 : n < 1 := by exact_mod_cast hnlt
    have h2 : (-1 : ℤ) < n := by exact_mod_cast hngt
    omega
  subst hn0
  simp only [Int.cast_zero, add_zero] at hlo hhi
  exact hfree c hc ⟨hlo, hhi⟩

end LRC14
