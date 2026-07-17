import TournamentH7.LRCSelectedWitnessGapTwoEscape
import TournamentH7.LRCPairTowerGapTwoProducer

/-!
# The exact residual of the `(4,4,8,8)` scalar escape

The two possible q4 anchors give an affine pair of scalar frequencies.  Their
difference is exactly one quarter of the normalized q4-numerator difference.
Consequently both frequencies can vanish only on the repeated-normalized-row
diagonal, and a q4 gap of size `60B` forces one of the two sharp large-frequency
gates `15B <= 2|F|`.

`LRCPairTowerGapTwoProducer` now manufactures both phase-matching walls from
the intrinsic `(4,4,8,8)` failure partition.  This module composes that producer
with the large-frequency consumer in `LRCSelectedWitnessGapTwoEscape` and
identifies the exact surviving case: both signed relations are below threshold
and the two q4 numerators lie in the strict `60B` pencil.

Tournament-analysis audit: the vertices are the two q4 anchor choices, and the
observable is their signed scalar wall frequency.  Thresholding `2|F|` at
`15B` gives the useful binary relation.  The two-vertex Hamiltonian path is the
q4-numerator order.  Quotienting further to runners loses the exact affine
difference, while quotienting to branch residues loses the frequency size.
-/

namespace LonelyRunner
namespace LRCPairTowerGapTwoFrequency

open LonelyRunner
open scoped Classical

noncomputable section

/-- The common residue normalization makes the q4488 scalar numerator exactly
divisible by eight. -/
theorem eight_dvd_fourEightEightPhaseNumerator
    (a₄ a₈a a₈b residue : ℤ)
    (hres₄ : (4 : ℤ) ∣ a₄ - residue)
    (hres₈a : (8 : ℤ) ∣ a₈a - residue)
    (hres₈b : (8 : ℤ) ∣ a₈b - residue) :
    (8 : ℤ) ∣ -2 * a₄ + a₈a + a₈b := by
  obtain ⟨k₄, hk₄⟩ := hres₄
  obtain ⟨k₈a, hk₈a⟩ := hres₈a
  obtain ⟨k₈b, hk₈b⟩ := hres₈b
  refine ⟨-k₄ + k₈a + k₈b, ?_⟩
  omega

/-- Exact cleared-denominator form of the q4488 phase frequency. -/
theorem eight_mul_fourEightEightPhaseFrequency_eq
    (a₄ a₈a a₈b residue : ℤ)
    (hres₄ : (4 : ℤ) ∣ a₄ - residue)
    (hres₈a : (8 : ℤ) ∣ a₈a - residue)
    (hres₈b : (8 : ℤ) ∣ a₈b - residue) :
    8 * fourEightEightPhaseFrequency a₄ a₈a a₈b =
      -2 * a₄ + a₈a + a₈b := by
  have hdiv := eight_dvd_fourEightEightPhaseNumerator
    a₄ a₈a a₈b residue hres₄ hres₈a hres₈b
  have hcancel := Int.ediv_mul_cancel hdiv
  simpa [fourEightEightPhaseFrequency, mul_comm] using hcancel

/-- The two q4-anchor frequencies form an affine pencil whose exact step is
one quarter of the q4-numerator step. -/
theorem four_mul_fourEightEightPhaseFrequency_sub_eq
    (a₄a a₄b a₈a a₈b residue : ℤ)
    (hres₄a : (4 : ℤ) ∣ a₄a - residue)
    (hres₄b : (4 : ℤ) ∣ a₄b - residue)
    (hres₈a : (8 : ℤ) ∣ a₈a - residue)
    (hres₈b : (8 : ℤ) ∣ a₈b - residue) :
    4 * (fourEightEightPhaseFrequency a₄a a₈a a₈b -
      fourEightEightPhaseFrequency a₄b a₈a a₈b) =
        a₄b - a₄a := by
  have ha := eight_mul_fourEightEightPhaseFrequency_eq
    a₄a a₈a a₈b residue hres₄a hres₈a hres₈b
  have hb := eight_mul_fourEightEightPhaseFrequency_eq
    a₄b a₈a a₈b residue hres₄b hres₈a hres₈b
  omega

/-- Both q4-derived frequencies vanish only on the diagonal where the two
normalized q4 numerators coincide. -/
theorem qFour_eq_of_both_fourEightEightPhaseFrequencies_zero
    (a₄a a₄b a₈a a₈b residue : ℤ)
    (hres₄a : (4 : ℤ) ∣ a₄a - residue)
    (hres₄b : (4 : ℤ) ∣ a₄b - residue)
    (hres₈a : (8 : ℤ) ∣ a₈a - residue)
    (hres₈b : (8 : ℤ) ∣ a₈b - residue)
    (hzeroA : fourEightEightPhaseFrequency a₄a a₈a a₈b = 0)
    (hzeroB : fourEightEightPhaseFrequency a₄b a₈a a₈b = 0) :
    a₄a = a₄b := by
  have hstep := four_mul_fourEightEightPhaseFrequency_sub_eq
    a₄a a₄b a₈a a₈b residue
      hres₄a hres₄b hres₈a hres₈b
  rw [hzeroA, hzeroB] at hstep
  omega

/-- Distinct normalized q4 rows force at least one nonzero wall frequency. -/
theorem one_fourEightEightPhaseFrequency_ne_zero_of_qFour_ne
    (a₄a a₄b a₈a a₈b residue : ℤ)
    (hres₄a : (4 : ℤ) ∣ a₄a - residue)
    (hres₄b : (4 : ℤ) ∣ a₄b - residue)
    (hres₈a : (8 : ℤ) ∣ a₈a - residue)
    (hres₈b : (8 : ℤ) ∣ a₈b - residue)
    (hne : a₄a ≠ a₄b) :
    fourEightEightPhaseFrequency a₄a a₈a a₈b ≠ 0 ∨
      fourEightEightPhaseFrequency a₄b a₈a a₈b ≠ 0 := by
  by_cases hzeroA : fourEightEightPhaseFrequency a₄a a₈a a₈b = 0
  · right
    intro hzeroB
    exact hne (qFour_eq_of_both_fourEightEightPhaseFrequencies_zero
      a₄a a₄b a₈a a₈b residue
        hres₄a hres₄b hres₈a hres₈b hzeroA hzeroB)
  · exact Or.inl hzeroA

/-- A normalized q4 gap of `60B` forces one of the two sharp q4488 escape
gates.  Thus the unresolved large-frequency complement is a close affine
pencil, not two unrelated small frequencies. -/
theorem one_fourEightEightPhaseFrequency_large_of_qFour_gap
    (a₄a a₄b a₈a a₈b residue : ℤ) (B : ℝ)
    (hres₄a : (4 : ℤ) ∣ a₄a - residue)
    (hres₄b : (4 : ℤ) ∣ a₄b - residue)
    (hres₈a : (8 : ℤ) ∣ a₈a - residue)
    (hres₈b : (8 : ℤ) ∣ a₈b - residue)
    (hgap : 60 * B ≤ |(a₄b : ℝ) - a₄a|) :
    15 * B ≤
        2 * |(fourEightEightPhaseFrequency a₄a a₈a a₈b : ℝ)| ∨
      15 * B ≤
        2 * |(fourEightEightPhaseFrequency a₄b a₈a a₈b : ℝ)| := by
  by_contra hlarge
  push Not at hlarge
  rcases hlarge with ⟨hsmallA, hsmallB⟩
  have hstepZ := four_mul_fourEightEightPhaseFrequency_sub_eq
    a₄a a₄b a₈a a₈b residue
      hres₄a hres₄b hres₈a hres₈b
  have hstepR :
      (4 : ℝ) *
          ((fourEightEightPhaseFrequency a₄a a₈a a₈b : ℝ) -
            (fourEightEightPhaseFrequency a₄b a₈a a₈b : ℝ)) =
        (a₄b : ℝ) - a₄a := by
    exact_mod_cast hstepZ
  have htriangle := abs_sub
    (fourEightEightPhaseFrequency a₄a a₈a a₈b : ℝ)
    (fourEightEightPhaseFrequency a₄b a₈a a₈b : ℝ)
  have habsStep :
      |(a₄b : ℝ) - a₄a| =
        4 * |(fourEightEightPhaseFrequency a₄a a₈a a₈b : ℝ) -
          (fourEightEightPhaseFrequency a₄b a₈a a₈b : ℝ)| := by
    rw [← hstepR, abs_mul]
    norm_num
  rw [habsStep] at hgap
  nlinarith

#print axioms eight_dvd_fourEightEightPhaseNumerator
#print axioms eight_mul_fourEightEightPhaseFrequency_eq
#print axioms four_mul_fourEightEightPhaseFrequency_sub_eq
#print axioms qFour_eq_of_both_fourEightEightPhaseFrequencies_zero
#print axioms one_fourEightEightPhaseFrequency_ne_zero_of_qFour_ne
#print axioms one_fourEightEightPhaseFrequency_large_of_qFour_gap

end
end LRCPairTowerGapTwoFrequency
end LonelyRunner

namespace LonelyRunner
namespace LRCPairTowerGapTwo

open LonelyRunner LRCPairTowerGapTwoFrequency LRCPairTowerGapTwoProducer
open scoped Classical

noncomputable section

/-- Phase-uniform matching-wall interface supplied by the intrinsic q4488
failure producer.  Naming it keeps the affine-pencil composition readable. -/
def FourEightEightMatchingWall
    (v : Fin 13 → ℤ) (g : ℤ) (i₁ i₂ i₃ i₄ : Fin 13)
    (a₄ a₈a a₈b : ℤ) : Prop :=
  ∀ u : ℝ,
    ¬ HasFourDetunedGoodBranch (v i₁) (v i₂) (v i₃) (v i₄) g u →
    ∃ c₄ c₈a c₈b : ℤ,
      (8 : ℤ) ∣ -2 * c₄ + c₈a + c₈b ∧
      (∃ n : ℤ,
        |(a₄ : ℝ) * ((u + (c₄ : ℝ)) / 4) - n| < 1 / 14) ∧
      (∃ n : ℤ,
        |(a₈a : ℝ) * ((u + (c₈a : ℝ)) / 8) - n| < 1 / 14) ∧
      (∃ n : ℤ,
        |(a₈b : ℝ) * ((u + (c₈b : ℝ)) / 8) - n| < 1 / 14)

/-- Two matching walls plus a `60B` separation of their q4 anchors close the
q4488 selected phase.  This is the direct composition of the affine-pencil
identity with the sharp `LRC(9)` interval consumer. -/
theorem selectedWitness_of_cite_and_two_matchingWalls_of_qFour_gap
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (i₁ i₂ i₃ i₄ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h14 : i₁ ≠ i₄)
    (h23 : i₂ ≠ i₃) (h24 : i₂ ≠ i₄) (h34 : i₃ ≠ i₄)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ → g ∣ v j)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ →
      |((v j / g : ℤ) : ℝ)| ≤ B)
    (a₄a a₄b a₈a a₈b residue : ℤ)
    (hres₄a : (4 : ℤ) ∣ a₄a - residue)
    (hres₄b : (4 : ℤ) ∣ a₄b - residue)
    (hres₈a : (8 : ℤ) ∣ a₈a - residue)
    (hres₈b : (8 : ℤ) ∣ a₈b - residue)
    (hmatchingA : FourEightEightMatchingWall
      v g i₁ i₂ i₃ i₄ a₄a a₈a a₈b)
    (hmatchingB : FourEightEightMatchingWall
      v g i₁ i₂ i₃ i₄ a₄b a₈a a₈b)
    (hgap : 60 * B ≤ |(a₄b : ℝ) - a₄a|) :
    ∃ u : ℝ,
      FourDetunedHarmonicGoodAt v g i₁ i₂ i₃ i₄ u ∧
      HasFourDetunedGoodBranch (v i₁) (v i₂) (v i₃) (v i₄) g u := by
  rcases one_fourEightEightPhaseFrequency_large_of_qFour_gap
      a₄a a₄b a₈a a₈b residue B
        hres₄a hres₄b hres₈a hres₈b hgap with hlargeA | hlargeB
  · have hnonzeroA : fourEightEightPhaseFrequency a₄a a₈a a₈b ≠ 0 := by
      intro hzero
      rw [hzero] at hlargeA
      norm_num at hlargeA
      linarith
    exact selectedWitness_of_cite_and_fourEightEight_matchingWall
      cite v hv g i₁ i₂ i₃ i₄ h12 h13 h14 h23 h24 h34 hdvd
        B hB0 hB a₄a a₈a a₈b residue
        hres₄a hres₈a hres₈b hnonzeroA hlargeA hmatchingA
  · have hnonzeroB : fourEightEightPhaseFrequency a₄b a₈a a₈b ≠ 0 := by
      intro hzero
      rw [hzero] at hlargeB
      norm_num at hlargeB
      linarith
    exact selectedWitness_of_cite_and_fourEightEight_matchingWall
      cite v hv g i₁ i₂ i₃ i₄ h12 h13 h14 h23 h24 h34 hdvd
        B hB0 hB a₄b a₈a a₈b residue
        hres₄b hres₈a hres₈b hnonzeroB hlargeB hmatchingB

/-- **Producer-consumer closure for the q4488 branch.**  Starting from only
the four denominator equations and one failing anchor, the unconditional
producer supplies common signed primitive numerators and both matching walls.
Either the sharp interval consumer gives a selected witness, or the exact
residual is a pair of zero/small signed relations whose q4 anchors lie in the
strict `60B` affine pencil.

Thus every q4488 configuration with q4 separation at least `60B` is closed;
the displayed close pencil is the complete residual of this route. -/
theorem selectedWitness_or_close_q4488_relation_pencil_of_failure
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (i₁ i₂ i₃ i₄ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h14 : i₁ ≠ i₄)
    (h23 : i₂ ≠ i₃) (h24 : i₂ ≠ i₄) (h34 : i₃ ≠ i₄)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ → g ∣ v j)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ →
      |((v j / g : ℤ) : ℝ)| ≤ B)
    (hg : 2 ≤ g)
    (hq4a : g / (Int.gcd (v i₁) g : ℤ) = 4)
    (hq4b : g / (Int.gcd (v i₂) g : ℤ) = 4)
    (hq8a : g / (Int.gcd (v i₃) g : ℤ) = 8)
    (hq8b : g / (Int.gcd (v i₄) g : ℤ) = 8)
    (anchor : ℝ)
    (hfailure : ¬ HasFourDetunedGoodBranch
      (v i₁) (v i₂) (v i₃) (v i₄) g anchor) :
    (∃ u : ℝ,
      FourDetunedHarmonicGoodAt v g i₁ i₂ i₃ i₄ u ∧
      HasFourDetunedGoodBranch (v i₁) (v i₂) (v i₃) (v i₄) g u) ∨
    ∃ a₄a a₄b a₈a a₈b residue : ℤ,
      (4 : ℤ) ∣ a₄a - residue ∧
      (4 : ℤ) ∣ a₄b - residue ∧
      (8 : ℤ) ∣ a₈a - residue ∧
      (8 : ℤ) ∣ a₈b - residue ∧
      a₄a ≠ a₄b ∧
      |(a₄b : ℝ) - a₄a| < 60 * B ∧
      (((-2 * a₄a + a₈a + a₈b = 0) ∨
          (fourEightEightPhaseFrequency a₄a a₈a a₈b ≠ 0 ∧
            2 * |(fourEightEightPhaseFrequency a₄a a₈a a₈b : ℝ)| <
              15 * B)) ∧
        ((-2 * a₄b + a₈a + a₈b = 0) ∨
          (fourEightEightPhaseFrequency a₄b a₈a a₈b ≠ 0 ∧
            2 * |(fourEightEightPhaseFrequency a₄b a₈a a₈b : ℝ)| <
              15 * B)) ∧
        (fourEightEightPhaseFrequency a₄a a₈a a₈b ≠ 0 ∨
          fourEightEightPhaseFrequency a₄b a₈a a₈b ≠ 0)) := by
  obtain ⟨a₄a, a₄b, a₈a, a₈b, residue,
      _hnorm4a, _hnorm4b, _hnorm8a, _hnorm8b,
      hres₄a, hres₄b, hres₈a, hres₈b, hdistinct,
      hmatching, hsplit⟩ :=
    exists_unconditional_q4488_failure_producer
      (v i₁) (v i₂) (v i₃) (v i₄) g anchor hg
        hq4a hq4b hq8a hq8b B hB0 hfailure
  have hmatchingA : FourEightEightMatchingWall
      v g i₁ i₂ i₃ i₄ a₄a a₈a a₈b := by
    intro u hfail
    exact (hmatching u hfail).1
  have hmatchingB : FourEightEightMatchingWall
      v g i₁ i₂ i₃ i₄ a₄b a₈a a₈b := by
    intro u hfail
    exact (hmatching u hfail).2
  rcases hsplit with hlargeA | hlargeB | hresidual
  · left
    exact selectedWitness_of_cite_and_fourEightEight_matchingWall
      cite v hv g i₁ i₂ i₃ i₄ h12 h13 h14 h23 h24 h34 hdvd
        B hB0 hB a₄a a₈a a₈b residue
        hres₄a hres₈a hres₈b hlargeA.1 hlargeA.2 hmatchingA
  · left
    exact selectedWitness_of_cite_and_fourEightEight_matchingWall
      cite v hv g i₁ i₂ i₃ i₄ h12 h13 h14 h23 h24 h34 hdvd
        B hB0 hB a₄b a₈a a₈b residue
        hres₄b hres₈a hres₈b hlargeB.1 hlargeB.2 hmatchingB
  · right
    have hsmallA :
        2 * |(fourEightEightPhaseFrequency a₄a a₈a a₈b : ℝ)| <
          15 * B := by
      rcases hresidual.1 with hzeroA | hsmallA
      · have hfrequencyZero :
            fourEightEightPhaseFrequency a₄a a₈a a₈b = 0 := by
          have hcleared := eight_mul_fourEightEightPhaseFrequency_eq
            a₄a a₈a a₈b residue hres₄a hres₈a hres₈b
          omega
        rw [hfrequencyZero]
        norm_num
        linarith
      · exact hsmallA.2
    have hsmallB :
        2 * |(fourEightEightPhaseFrequency a₄b a₈a a₈b : ℝ)| <
          15 * B := by
      rcases hresidual.2.1 with hzeroB | hsmallB
      · have hfrequencyZero :
            fourEightEightPhaseFrequency a₄b a₈a a₈b = 0 := by
          have hcleared := eight_mul_fourEightEightPhaseFrequency_eq
            a₄b a₈a a₈b residue hres₄b hres₈a hres₈b
          omega
        rw [hfrequencyZero]
        norm_num
        linarith
      · exact hsmallB.2
    have hclose : |(a₄b : ℝ) - a₄a| < 60 * B := by
      apply lt_of_not_ge
      intro hwide
      rcases one_fourEightEightPhaseFrequency_large_of_qFour_gap
          a₄a a₄b a₈a a₈b residue B
            hres₄a hres₄b hres₈a hres₈b hwide with
        hlargeA | hlargeB
      · linarith
      · linarith
    exact ⟨a₄a, a₄b, a₈a, a₈b, residue,
      hres₄a, hres₄b, hres₈a, hres₈b,
      hdistinct, hclose, hresidual⟩

#print axioms selectedWitness_of_cite_and_two_matchingWalls_of_qFour_gap
#print axioms selectedWitness_or_close_q4488_relation_pencil_of_failure

end
end LRCPairTowerGapTwo
end LonelyRunner
