/-
  TournamentH7.LRC14WindowWiring — the hwindow NORMALIZATION WIRING
  (kind-pasteur-2026-07-02-S15, HYP-3972).  Single-writer satellite.

  klein-S107's certificate names two finite computations between the corpus and an
  unconditional `LRC14Statement`; one is `hwindow` — loneliness of EVERY bounded family
  `|v i| ≤ W`.  The window packs (opus S46/S47 WindowPack1, klein S108 WindowPack2, …)
  certify exactly the PRIMITIVE, COVERING, SORTED-POSITIVE tuples.  This file proves that
  class SUFFICES:

      hwindow_of_normalized_census :
        (every monotone positive primitive covering tuple ≤ W is lonely)  →
        (every nonzero family with |v i| ≤ W is lonely).

  The reduction chain, all kernel-pure and previously machine-checked: signs
  (`lonely_abs_iff`), relabeling (`Tuple.sort` + `lonely_comp_equiv`), integer scale
  (Nat-gcd primitivization + `lonely_exists_of_scale`), and the sieve dichotomy at the
  primitive level (`sieve_one_div`): a primitivized tuple that stopped covering misses
  some `q ∈ {2,…,14}` and is lonely at `1/q` directly.

  With this, the packs' per-band rows literally SUM to `hwindow W` (one finite
  membership case-split per band), and the final surface becomes
  `lrc14_of_dispatchComplete_and_census` — both remaining hypotheses census-shaped.
-/
import TournamentH7.LRC14Dispatch
import Mathlib.Data.Fin.Tuple.Sort
import Mathlib.Algebra.GCDMonoid.Finset

namespace LonelyRunner
namespace LRC14

/-- The Nat-gcd of a 13-tuple's absolute values. -/
def tupleGcd (w : Fin 13 → ℤ) : ℕ := Finset.univ.gcd fun i => (w i).natAbs

theorem tupleGcd_dvd (w : Fin 13 → ℤ) (i : Fin 13) : (tupleGcd w : ℤ) ∣ w i := by
  have h : tupleGcd w ∣ (w i).natAbs := Finset.gcd_dvd (Finset.mem_univ i)
  exact dvd_trans (Int.natCast_dvd_natCast.mpr h) (Int.natAbs_dvd.mpr dvd_rfl)

theorem tupleGcd_pos (w : Fin 13 → ℤ) (hw : ∀ i, w i ≠ 0) : 0 < tupleGcd w := by
  rcases Nat.eq_zero_or_pos (tupleGcd w) with h0 | hpos
  · exfalso
    have := tupleGcd_dvd w 0
    rw [h0] at this
    push_cast at this
    exact hw 0 (zero_dvd_iff.mp this)
  · exact hpos

/-- The primitivized tuple: divide out the gcd. -/
def primPart (w : Fin 13 → ℤ) : Fin 13 → ℤ := fun i => w i / (tupleGcd w : ℤ)

theorem primPart_mul (w : Fin 13 → ℤ) (i : Fin 13) :
    w i = (tupleGcd w : ℤ) * primPart w i := by
  unfold primPart
  exact (Int.mul_ediv_cancel' (tupleGcd_dvd w i)).symm

theorem primPart_pos (w : Fin 13 → ℤ) (hw : ∀ i, 0 < w i) (i : Fin 13) :
    0 < primPart w i := by
  have hne : ∀ j, w j ≠ 0 := fun j => (hw j).ne'
  have hg : 0 < (tupleGcd w : ℤ) := by exact_mod_cast tupleGcd_pos w hne
  have hmul := primPart_mul w i
  by_contra hle
  have hle2 : primPart w i ≤ 0 := le_of_not_gt hle
  have hnp : (tupleGcd w : ℤ) * primPart w i ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos hg.le hle2
  linarith [hw i]

theorem primPart_gcd_eq_one (w : Fin 13 → ℤ) (hw : ∀ i, w i ≠ 0) :
    tupleGcd (primPart w) = 1 := by
  have hg : 0 < tupleGcd w := tupleGcd_pos w hw
  have hkey : tupleGcd w * tupleGcd (primPart w) ∣ tupleGcd w := by
    apply Finset.dvd_gcd
    intro i _
    have h1 : (tupleGcd (primPart w) : ℤ) ∣ primPart w i := tupleGcd_dvd (primPart w) i
    have h2 := primPart_mul w i
    have h3 : ((tupleGcd w * tupleGcd (primPart w) : ℕ) : ℤ) ∣ w i := by
      push_cast
      rw [h2]
      exact mul_dvd_mul_left _ h1
    have h4 : (tupleGcd w * tupleGcd (primPart w)) ∣ (w i).natAbs := by
      rw [← Int.natAbs_natCast (tupleGcd w * tupleGcd (primPart w))]
      exact Int.natAbs_dvd_natAbs.mpr h3
    exact h4
  have hle := Nat.le_of_dvd hg hkey
  have hpos : 0 < tupleGcd (primPart w) := by
    apply tupleGcd_pos
    intro i h0
    have := primPart_mul w i
    rw [h0, mul_zero] at this
    exact hw i this
  nlinarith [hle, hpos, hg]

/-- **THE WINDOW WIRING THEOREM**: the packs' census class (monotone, positive,
primitive, covering, bounded) suffices for the full `hwindow` signature. -/
theorem hwindow_of_normalized_census (W : ℤ)
    (hcensus : ∀ u : Fin 13 → ℤ, (∀ i, 0 < u i) → Monotone u → (∀ i, u i ≤ W) →
      CoveringFamily u → tupleGcd u = 1 → ∃ t : ℝ, Lonely 14 u t) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → (∀ i, |v i| ≤ W) → ∃ t : ℝ, Lonely 14 v t := by
  intro v hv hW
  -- Step 1: absolute values
  set va : Fin 13 → ℤ := fun i => |v i| with hva
  have hva_pos : ∀ i, 0 < va i := fun i => abs_pos.mpr (hv i)
  -- Step 2: sort
  set σ : Equiv.Perm (Fin 13) := Tuple.sort va with hσ
  set w : Fin 13 → ℤ := va ∘ σ with hw
  have hw_mono : Monotone w := Tuple.monotone_sort va
  have hw_pos : ∀ i, 0 < w i := fun i => hva_pos (σ i)
  have hw_W : ∀ i, w i ≤ W := fun i => by
    have := hW (σ i)
    simpa [hw, hva] using this
  have hw_ne : ∀ i, w i ≠ 0 := fun i => (hw_pos i).ne'
  -- Step 3: primitivize
  set g : ℕ := tupleGcd w with hg
  have hg_pos : 0 < g := tupleGcd_pos w hw_ne
  set u : Fin 13 → ℤ := primPart w with hu
  have hu_pos : ∀ i, 0 < u i := primPart_pos w hw_pos
  have hu_mul : ∀ i, w i = (g : ℤ) * u i := primPart_mul w
  have hu_mono : Monotone u := by
    intro i j hij
    have hgZ : (0 : ℤ) < (g : ℤ) := by exact_mod_cast hg_pos
    have := hw_mono hij
    rw [hu_mul i, hu_mul j] at this
    exact le_of_mul_le_mul_left this hgZ
  have hu_W : ∀ i, u i ≤ W := fun i => by
    have h1 : u i ≤ w i := by
      have hgZ : (1 : ℤ) ≤ (g : ℤ) := by exact_mod_cast hg_pos
      have := hu_mul i
      nlinarith [hu_pos i]
    exact le_trans h1 (hw_W i)
  have hu_prim : tupleGcd u = 1 := primPart_gcd_eq_one w hw_ne
  -- Step 4: sieve-or-census dichotomy at the primitive level
  have hu_lonely : ∃ t : ℝ, Lonely 14 u t := by
    by_cases hcov : CoveringFamily u
    · exact hcensus u hu_pos hu_mono hu_W hcov hu_prim
    · unfold CoveringFamily at hcov
      push Not at hcov
      obtain ⟨q, hq2, hq14, hdiv⟩ := hcov
      exact ⟨(1 : ℝ) / q, sieve_one_div 14 q u hq14 (by omega) hdiv⟩
  -- Transport back: scale, then relabeling, then signs
  have hw_lonely : ∃ t : ℝ, Lonely 14 w t := by
    have := lonely_exists_of_scale 14 u (g : ℤ) (by exact_mod_cast hg_pos.ne') hu_lonely
    obtain ⟨t, ht⟩ := this
    refine ⟨t, ?_⟩
    have hfun : w = fun i => (g : ℤ) * u i := funext hu_mul
    rw [hfun]
    exact ht
  obtain ⟨t, ht⟩ := hw_lonely
  have hva_lonely : Lonely 14 va t := (lonely_comp_equiv 14 va t σ).mp ht
  exact ⟨t, (lonely_abs_iff 14 v t).mp hva_lonely⟩

/-- **THE FINAL SURFACE, both hypotheses census-shaped**: LRC(14) from the THM-602
dispatch completeness plus the normalized bounded census.  Everything else in both
reductions is machine-checked, kernel-pure glue. -/
theorem lrc14_of_dispatchComplete_and_census (W : ℤ)
    (hcomp : DispatchComplete W)
    (hcensus : ∀ u : Fin 13 → ℤ, (∀ i, 0 < u i) → Monotone u → (∀ i, u i ≤ W) →
      CoveringFamily u → tupleGcd u = 1 → ∃ t : ℝ, Lonely 14 u t) :
    LRC14Statement :=
  lrc14_of_dispatchComplete W hcomp (hwindow_of_normalized_census W hcensus)

end LRC14
end LonelyRunner
