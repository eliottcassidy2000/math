/-
  TournamentH7.LRCSelectedWitnessFrequencyEndpoint

  The exact unit-budget endpoint applies to *distinct absolute harmonic
  quotient frequencies*, rather than to the ten off-selected runner labels.
  Thus a three-detuned decomposition with at most seven such frequencies has
  a harmonic-good phase without invoking the LRC(<=13) citation.  Equivalently,
  failure of harmonic clearance forces at least eight absolute frequencies.

  Tournament / alternate-vertex audit.  The vertices here are equivalence
  classes of off-selected harmonic quotients under `q ~ -q`, not runners,
  arcs, or wall-crossing events.  The pairwise observable is equality of
  absolute quotient frequency; the sign switch is the gauge `q |-> q.natAbs`.
  This quotient preserves every distance-to-integer constraint and destroys
  the sign and multiplicity of its runner labels.  Equality gives a tie graph,
  not a binary orientation, so a tie Hamiltonian path would be merely an
  arbitrary ordering of frequency classes; tournament cycle/SCC/edge-flip
  fingerprints carry no additional information in this reduction.  The
  challenged assumption is therefore that tournament vertices should be
  runners: here the proof-relevant vertices are frequency classes.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCDetunedOverlap
import TournamentH7.UnitBudgetEndpoint

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

noncomputable section

/-- Distinct absolute values of the harmonic quotient frequencies away from
the three selected indices.  Passing to `natAbs` is lossless for the
distance-to-integer predicate, while removing repetitions is exactly what lets
the endpoint theorem see fewer than ten frequencies. -/
def threeDetunedAbsQuotientFrequencies
    (v : Fin 13 → ℤ) (g : ℤ) (i₁ i₂ i₃ : Fin 13) : Finset ℕ :=
  (Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))).image
    (fun j => (v j / g).natAbs)

private theorem offSelected_card_eq_ten
    (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃) :
    (Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))).card = 10 := by
  have hselectedCard : ({i₁, i₂, i₃} : Finset (Fin 13)).card = 3 :=
    Finset.card_eq_three.mpr ⟨i₁, i₂, i₃, h12, h13, h23, rfl⟩
  rw [Finset.card_sdiff, Finset.inter_univ, hselectedCard, Finset.card_univ,
    Fintype.card_fin]

private theorem quotient_ne_zero_of_offSelected
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ)
    (i₁ i₂ i₃ j : Fin 13)
    (hdvd : ∀ k, k ≠ i₁ → k ≠ i₂ → k ≠ i₃ → g ∣ v k)
    (hj : j ∈ Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))) :
    v j / g ≠ 0 := by
  rw [Finset.mem_sdiff, Finset.mem_insert, Finset.mem_insert,
    Finset.mem_singleton] at hj
  push Not at hj
  intro hzero
  have hfactor : v j = g * (v j / g) :=
    (Int.mul_ediv_cancel' (hdvd j hj.2.1 hj.2.2.1 hj.2.2.2)).symm
  apply hv j
  rw [hfactor, hzero, mul_zero]

private theorem absQuotientFrequencies_nonempty
    (v : Fin 13 → ℤ) (g : ℤ) (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃) :
    (threeDetunedAbsQuotientFrequencies v g i₁ i₂ i₃).Nonempty := by
  have hcard := offSelected_card_eq_ten i₁ i₂ i₃ h12 h13 h23
  obtain ⟨j, hj⟩ := Finset.card_pos.mp (by omega :
    0 < (Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))).card)
  refine ⟨(v j / g).natAbs, ?_⟩
  exact Finset.mem_image.mpr ⟨j, hj, rfl⟩

/-- **Seven-frequency harmonic endpoint.**  If the ten off-selected harmonic
quotients occupy at most seven absolute frequency classes, compactness at the
exact unit-budget endpoint supplies a phase clearing every class by `1/14`.
This result is independent of `LRCUpTo13`. -/
theorem exists_threeDetunedHarmonicGoodAt_of_absQuotient_card_le_seven
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (hcard : (threeDetunedAbsQuotientFrequencies v g i₁ i₂ i₃).card ≤ 7) :
    ∃ u : ℝ, ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u := by
  let W := threeDetunedAbsQuotientFrequencies v g i₁ i₂ i₃
  have hWnonempty : W.Nonempty :=
    absQuotientFrequencies_nonempty v g i₁ i₂ i₃ h12 h13 h23
  have hWcard : 0 < W.card := Finset.card_pos.mpr hWnonempty
  have hcardW : W.card ≤ 7 := by simpa [W] using hcard
  have hWpos : ∀ w ∈ W, 0 < w := by
    intro w hw
    rw [show W = threeDetunedAbsQuotientFrequencies v g i₁ i₂ i₃ from rfl,
      threeDetunedAbsQuotientFrequencies, Finset.mem_image] at hw
    obtain ⟨j, hj, rfl⟩ := hw
    exact Int.natAbs_pos.mpr
      (quotient_ne_zero_of_offSelected v hv g i₁ i₂ i₃ j hdvd hj)
  obtain ⟨u, _, hu⟩ := LRC14.exists_lonely_unit_endpoint W hWpos hWcard
  refine ⟨u, ?_⟩
  intro j hj1 hj2 hj3 integer
  have hjOff : j ∈ Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13)) := by
    rw [Finset.mem_sdiff, Finset.mem_insert, Finset.mem_insert,
      Finset.mem_singleton]
    exact ⟨Finset.mem_univ j, by
      push Not
      exact ⟨hj1, hj2, hj3⟩⟩
  have hw : (v j / g).natAbs ∈ W := by
    exact Finset.mem_image.mpr ⟨j, hjOff, rfl⟩
  have hgap : (1 : ℝ) / 14 ≤ 1 / (2 * W.card) := by
    apply one_div_le_one_div_of_le (by positivity)
    have hcardNat : 2 * W.card ≤ 14 := by omega
    exact_mod_cast hcardNat
  have hq : v j / g ≠ 0 :=
    quotient_ne_zero_of_offSelected v hv g i₁ i₂ i₃ j hdvd hjOff
  rcases lt_or_gt_of_ne hq with hneg | hpos
  · have hqReal : ((v j / g : ℤ) : ℝ) = -(((v j / g).natAbs : ℕ) : ℝ) := by
      exact_mod_cast congrArg (fun z : ℤ => (z : ℝ))
        (by omega : v j / g = -((v j / g).natAbs : ℤ))
    calc
      (1 : ℝ) / 14 ≤ 1 / (2 * W.card) := hgap
      _ ≤ |((v j / g).natAbs : ℝ) * u - ((-integer : ℤ) : ℝ)| :=
        hu (v j / g).natAbs hw (-integer)
      _ = |((v j / g : ℤ) : ℝ) * u - (integer : ℝ)| := by
        rw [Int.cast_neg, hqReal]
        calc
          |((v j / g).natAbs : ℝ) * u - -(integer : ℝ)| =
              |-(-((v j / g).natAbs : ℝ) * u - (integer : ℝ))| := by
                congr 1
                ring
          _ = |-((v j / g).natAbs : ℝ) * u - (integer : ℝ)| := abs_neg _
  · have hqInt : (((v j / g).natAbs : ℕ) : ℤ) = v j / g :=
      Int.natAbs_of_nonneg hpos.le
    have hqReal : (((v j / g).natAbs : ℕ) : ℝ) = ((v j / g : ℤ) : ℝ) := by
      have hcast : (((((v j / g).natAbs : ℕ) : ℤ)) : ℝ) =
          (((v j / g).natAbs : ℕ) : ℝ) := Int.cast_natCast _
      rw [← hcast, hqInt]
    simpa only [hqReal] using hgap.trans (hu (v j / g).natAbs hw integer)

/-- Contrapositive structural residue: a decomposition with no harmonic-good
phase must contain at least eight distinct absolute harmonic quotient
frequencies. -/
theorem eight_absQuotient_frequencies_of_no_threeDetunedHarmonicGoodAt
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (hno : ¬ ∃ u : ℝ, ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u) :
    8 ≤ (threeDetunedAbsQuotientFrequencies v g i₁ i₂ i₃).card := by
  by_contra hlt
  have hcard :
      (threeDetunedAbsQuotientFrequencies v g i₁ i₂ i₃).card ≤ 7 := by
    omega
  exact hno (exists_threeDetunedHarmonicGoodAt_of_absQuotient_card_le_seven
    v hv g i₁ i₂ i₃ h12 h13 h23 hdvd hcard)

/-! ## Axiom audit -/

#print axioms exists_threeDetunedHarmonicGoodAt_of_absQuotient_card_le_seven
#print axioms eight_absQuotient_frequencies_of_no_threeDetunedHarmonicGoodAt

end

end LRC14Grand
end LonelyRunner
