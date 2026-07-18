/-
  TournamentH7.LRCPairTowerFreshBranch

  Exact fixed-phase classification of the first pair-tower lift.  Once at
  least one coordinate crosses the wall from `g` to `2g`, a branch clears the
  inherited pair and every fresh row if and only if neither obstruction named
  by `PairTowerParallelObstruction` occurs.  Thus the existing obstruction is
  not merely sufficient at a fixed phase: it is exact.

  This theorem does not choose the phase and does not prove
  `ManyLiftFailurePhaseSelector` or `NonterminatingPairTowerSupply`.  Harmonic
  safety of the doubled core remains a separate condition.  The q4488
  small-frequency residue and the dense-core B5 payment are also independent.

  Tournament-analysis audit: the faithful vertices are fresh divisibility-wall
  indices together with branch residues.  Pair membership in
  `detunedBadBranches` is the observable; overlap versus q22 opposition is the
  binary switch, and integer order on `Ico 0 (2*g)` is the tie Hamiltonian path.
  This carrier preserves the common fixed phase and simultaneous branch
  clearing.  A runner tournament loses bad-row equality and the q244 common
  branch.  The challenged assumption is that absence of the named obstruction
  was only a one-way selector hypothesis; the theorem below proves the converse
  at each fixed phase while deliberately forgetting phase motion and core
  harmonic safety.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCPairTowerReduction

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

/-- At one fixed phase, one branch clears both inherited q-four rows and every
fresh q-two row exposed between `g` and `2g`. -/
def PairTowerFreshBranchAt
    (v : Fin 13 → ℤ) (g : ℤ) (first second : Fin 13) (u : ℝ) : Prop :=
  ∃ branch : ℤ,
    branch ∈ Finset.Ico 0 (2 * g) ∧
    branch ∉ detunedBadBranches (v first) (2 * g) u ∧
    branch ∉ detunedBadBranches (v second) (2 * g) u ∧
    ∀ j : Fin 13, (g ∣ v j ∧ ¬ 2 * g ∣ v j) →
      branch ∉ detunedBadBranches (v j) (2 * g) u

/-- Exact fixed-phase pair-tower classifier.  With a nonempty fresh layer, a
simultaneous clearing branch exists exactly when there is neither a q22
opposition among fresh rows nor a q244 failure against the inherited pair. -/
theorem pairTowerFreshBranchAt_iff_not_parallelObstruction
    (v : Fin 13 → ℤ) (g : ℤ) (first second : Fin 13) (u : ℝ)
    (hg : 2 ≤ g) (hfresh : 0 < liftFailureCard v g) :
    PairTowerFreshBranchAt v g first second u ↔
      ¬ PairTowerParallelObstruction v g first second u := by
  constructor
  · rintro ⟨branch, hbranch, hfirstGood, hsecondGood, hallFreshGood⟩
    rintro (h22 | h244)
    · obtain ⟨i, j, hi, hj, _hij, hopposition⟩ := h22
      have hqI := odd_harmonic_lifts_to_q_two
        (v i) g hg hi.1 hi.2
      have hqJ := odd_harmonic_lifts_to_q_two
        (v j) g hg hj.1 hj.2
      have hgood :
          HasThreeDetunedGoodBranch
            (v i) (v j) (v first) (2 * g) u :=
        ⟨branch, hbranch, hallFreshGood i hi, hallFreshGood j hj,
          hfirstGood⟩
      exact (not_hasThreeDetunedGoodBranch_two_two_three_of_opposition
        (v i) (v j) (v first) (2 * g) u (by omega)
          hqI hqJ hopposition) hgood
    · obtain ⟨i, hi, hfail⟩ := h244
      exact hfail
        ⟨branch, hbranch, hallFreshGood i hi, hfirstGood, hsecondGood⟩
  · intro hnoObstruction
    have hfreshNonempty :
        (Finset.univ.filter fun i =>
          g ∣ v i ∧ ¬ 2 * g ∣ v i).Nonempty := by
      rw [← Finset.card_pos]
      simpa [liftFailureCard] using hfresh
    obtain ⟨defaultIndex, hdefaultMem⟩ := hfreshNonempty
    have hdefault :
        g ∣ v defaultIndex ∧ ¬ 2 * g ∣ v defaultIndex := by
      simpa using hdefaultMem
    obtain ⟨fresh, hfreshRow, hmode⟩ :
        ∃ i : Fin 13, (g ∣ v i ∧ ¬ 2 * g ∣ v i) ∧
          ((detunedBadBranches (v i) (2 * g) u).Nonempty ∨
            ∀ j : Fin 13, (g ∣ v j ∧ ¬ 2 * g ∣ v j) →
              ¬ (detunedBadBranches (v j) (2 * g) u).Nonempty) := by
      by_cases hactive : ∃ i : Fin 13,
          (g ∣ v i ∧ ¬ 2 * g ∣ v i) ∧
            (detunedBadBranches (v i) (2 * g) u).Nonempty
      · obtain ⟨i, hi, hrow⟩ := hactive
        exact ⟨i, hi, Or.inl hrow⟩
      · refine ⟨defaultIndex, hdefault, Or.inr ?_⟩
        intro j hj hrow
        exact hactive ⟨j, hj, hrow⟩
    have hgood :
        HasThreeDetunedGoodBranch
          (v fresh) (v first) (v second) (2 * g) u := by
      by_contra hfail
      exact hnoObstruction (Or.inr ⟨fresh, hfreshRow, hfail⟩)
    obtain ⟨branch, hbranch, hfreshGood, hfirstGood, hsecondGood⟩ := hgood
    have hallFreshGood : ∀ j : Fin 13,
        (g ∣ v j ∧ ¬ 2 * g ∣ v j) →
        branch ∉ detunedBadBranches (v j) (2 * g) u := by
      intro j hj
      rcases hmode with hfreshActive | hallEmpty
      · by_cases hjActive :
          (detunedBadBranches (v j) (2 * g) u).Nonempty
        · by_cases hjFresh : j = fresh
          · subst j
            exact hfreshGood
          · have hnotOpposition :
                ¬ TwoTwoThreePhaseOpposition
                  (v fresh) (v j) (2 * g) u := by
              intro hopposition
              exact hnoObstruction (Or.inl
                ⟨fresh, j, hfreshRow, hj, Ne.symm hjFresh, hopposition⟩)
            have hnotDisjoint :
                ¬ Disjoint
                  (detunedBadBranches (v fresh) (2 * g) u)
                  (detunedBadBranches (v j) (2 * g) u) := by
              intro hdisjoint
              exact hnotOpposition ⟨hfreshActive, hjActive, hdisjoint⟩
            have hinter :
                (detunedBadBranches (v fresh) (2 * g) u ∩
                  detunedBadBranches (v j) (2 * g) u).Nonempty := by
              obtain ⟨branch', hleft, hright⟩ :=
                Finset.not_disjoint_iff.mp hnotDisjoint
              exact ⟨branch', Finset.mem_inter.mpr ⟨hleft, hright⟩⟩
            have hqFresh := odd_harmonic_lifts_to_q_two
              (v fresh) g hg hfreshRow.1 hfreshRow.2
            have hqJ := odd_harmonic_lifts_to_q_two
              (v j) g hg hj.1 hj.2
            have hrowsEq :=
              detunedBadBranches_eq_of_overlap_same_reducedDenominator
                (v fresh) (v j) (2 * g) 2 u (by omega) (by norm_num)
                  hqFresh hqJ hinter
            rw [hrowsEq] at hfreshGood
            exact hfreshGood
        · intro hjBranch
          exact hjActive ⟨branch, hjBranch⟩
      · intro hjBranch
        exact hallEmpty j hj ⟨branch, hjBranch⟩
    exact ⟨branch, hbranch, hfirstGood, hsecondGood, hallFreshGood⟩

#print axioms pairTowerFreshBranchAt_iff_not_parallelObstruction

end LRC14Grand
end LonelyRunner
