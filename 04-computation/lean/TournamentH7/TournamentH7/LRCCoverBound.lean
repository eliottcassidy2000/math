/-
  TournamentH7.LRCCoverBound -- elementary cores of the wide cover bound `hp0cap`
  (`p0(E) ≤ cap_k`), kind-pasteur-2026-06-22-S31, HYP-2839.

  hp0cap = `p0(E) = meas(coverSet E) ≤ cap_k` for wide `E` is the sector route's
  deep node.  Its full proof needs the analytic decorrelation/Tornheim tail.  This
  file formalizes the ELEMENTARY, sorry-free cores (the analogue of `LRCGapReach`
  for hpartA):

  * **coverage monotonicity** -- `coverSet` is monotone in the speed list, so `p0`
    is monotone; the wide bound's worst case is the full speed set.
  * **the small-k pigeonhole** -- six disjoint inner sectors need six distinct
    speeds, so `|E| < 6 ⟹ coverSet E = ∅ ⟹ p0(E) = 0` (the analogue of the
    hnu1 max-gap pigeonhole).  Closes hp0cap trivially for `k < 6`.

  The binding cases `k = 8..12` remain the analytic residual (`p0 ≤ p0_decorr ≤
  Q(k-1) < cap`, the resonance/Tornheim tail), isolated in `LRCWideBoundReduction`.
-/

import TournamentH7.LRCDenseCovers

namespace LonelyRunner
namespace DenseCovers

open MeasureTheory

/-! ## Coverage monotonicity -/

/-- `coverSet` is monotone in the speed list: more speeds ⟹ easier to hit every
sector. -/
theorem coverSet_mono {E E' : List ℤ} (h : ∀ e ∈ E, e ∈ E') :
    coverSet E ⊆ coverSet E' := by
  intro x hx j hj1 hj6
  obtain ⟨e, heE, hb⟩ := hx j hj1 hj6
  exact ⟨e, h e heE, hb⟩

/-- `p0 = meas(coverSet)` is monotone in the speed list. -/
theorem slowμ_coverSet_mono {E E' : List ℤ} (h : ∀ e ∈ E, e ∈ E') :
    slowμ (coverSet E) ≤ slowμ (coverSet E') :=
  measure_mono (coverSet_mono h)

/-! ## The small-k pigeonhole: `|E| < 6 ⟹ p0 = 0` -/

/-- **Six disjoint inner sectors need six distinct speeds.**  If at some slow time
`x` the phases `{frac(e·x) : e ∈ E}` hit all six inner sectors, then `E` has at
least six distinct elements.  (Different sectors are witnessed by different speeds,
since one phase lies in one half-open sector.) -/
theorem six_le_card_of_coverSet_mem {E : List ℤ} {x : ℝ} (hx : x ∈ coverSet E) :
    6 ≤ E.toFinset.card := by
  classical
  -- a witnessing speed for inner sector `m = (j:ℕ)+1 ∈ {1,…,6}`, indexed by `j : Fin 6`
  have hw : ∀ j : Fin 6, ∃ e ∈ E,
      (((j : ℕ) + 1 : ℕ) : ℝ) / 7 ≤ Int.fract ((e : ℝ) * x) ∧
        Int.fract ((e : ℝ) * x) < ((((j : ℕ) + 1 : ℕ) : ℝ) + 1) / 7 := by
    intro j
    exact hx ((j : ℕ) + 1) (by omega) (by have := j.isLt; omega)
  choose w hwE hwlo hwhi using hw
  -- `w` is injective: a shared speed would put one phase in two disjoint sectors
  have hinj : Function.Injective w := by
    intro a b hab
    have hla := hwlo a
    have hhb := hwhi b
    have hlb := hwlo b
    have hha := hwhi a
    rw [hab] at hla hha
    -- now `frac (w b · x)` lies in both `[(a+1)/7,(a+2)/7)` and `[(b+1)/7,(b+2)/7)`
    have hab1 : (((a : ℕ) + 1 : ℕ) : ℝ) < (((b : ℕ) + 1 : ℕ) : ℝ) + 1 := by
      have := lt_of_le_of_lt hla hhb; linarith
    have hba1 : (((b : ℕ) + 1 : ℕ) : ℝ) < (((a : ℕ) + 1 : ℕ) : ℝ) + 1 := by
      have := lt_of_le_of_lt hlb hha; linarith
    have hab2 : (a : ℕ) + 1 < (b : ℕ) + 1 + 1 := by exact_mod_cast hab1
    have hba2 : (b : ℕ) + 1 < (a : ℕ) + 1 + 1 := by exact_mod_cast hba1
    exact Fin.ext (by omega)
  -- the six witness speeds sit inside `E.toFinset`
  have hsub : Finset.image w Finset.univ ⊆ E.toFinset := by
    intro e he
    simp only [Finset.mem_image, Finset.mem_univ, true_and] at he
    obtain ⟨j, rfl⟩ := he
    exact List.mem_toFinset.mpr (hwE j)
  calc (6 : ℕ) = (Finset.univ : Finset (Fin 6)).card := by simp
    _ = (Finset.image w Finset.univ).card :=
        (Finset.card_image_of_injective _ hinj).symm
    _ ≤ E.toFinset.card := Finset.card_le_card hsub

/-- **Small-k vanishing.**  If `E` has fewer than six distinct speeds, no slow time
covers all six inner sectors: `coverSet E = ∅`. -/
theorem coverSet_eq_empty_of_card_lt_six {E : List ℤ} (h : E.toFinset.card < 6) :
    coverSet E = ∅ := by
  ext x
  simp only [Set.mem_empty_iff_false, iff_false]
  intro hx
  exact absurd (six_le_card_of_coverSet_mem hx) (by omega)

/-- **hp0cap for `k < 6`:** `p0(E) = 0 ≤ cap`.  The cover bound is trivial below
six speeds. -/
theorem slowμ_coverSet_eq_zero_of_card_lt_six {E : List ℤ} (h : E.toFinset.card < 6) :
    slowμ (coverSet E) = 0 := by
  rw [coverSet_eq_empty_of_card_lt_six h, measure_empty]

/-! ## Axiom audit -/

#print axioms coverSet_mono
#print axioms six_le_card_of_coverSet_mem
#print axioms slowμ_coverSet_eq_zero_of_card_lt_six

end DenseCovers
end LonelyRunner
