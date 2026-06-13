/-
  TournamentH7.StaircaseConnectivity — connecting concrete tilings to THM-330

  This module packages the bridge from the concrete `StTiling` coordinate
  model to the abstract base-path tournament theorem:

    * `StTiling.toTournament` has the base path;
    * a good cut of a tiling is exactly a crossing-upward cut of the induced
      tournament;
    * therefore the top good-cut bucket is exactly strong connectivity.
-/

import TournamentH7.GoodCuts
import TournamentH7.StaircaseModel

namespace Tournament

variable {n : ℕ}

noncomputable section

/-- Good cuts in the concrete tiling model are exactly crossing-upward cuts
    in the induced base-path tournament. -/
theorem StTiling.isGoodCut_iff_crossesUpward_toTournament
    (b : StTiling n) {k : ℕ} :
    StTiling.IsGoodCut b k ↔ CrossesUpward b.toTournament k := by
  constructor
  · rintro ⟨t, ht, hcross⟩
    refine ⟨t.hi, t.lo, hcross.2, hcross.1, t.gap, ?_⟩
    calc
      b.toTournament.arc t.lo t.hi = b.arc t.lo t.hi := rfl
      _ = b t := StTiling.arc_tile_up b t
      _ = true := ht
  · rintro ⟨i, j, hi_ge, hj_lt, hgap, harc⟩
    let t : StTile n := ⟨i, j, hgap⟩
    have harc' : b.arc j i = true := by
      simpa using harc
    have ht : b t = true := by
      rw [← StTiling.arc_nonconsec_up (b := b) (lo := j) (hi := i) hgap]
      exact harc'
    exact ⟨t, ht, ⟨hj_lt, hi_ge⟩⟩

/-- Set form: all legal cuts are good iff every THM-330 cut has a
    crossing-upward arc in the induced tournament. -/
theorem StTiling.all_goodCuts_iff_all_crossesUpward_toTournament
    (b : StTiling n) :
    (∀ k, k ∈ cutSet n → StTiling.IsGoodCut b k) ↔
      ∀ k, 1 ≤ k → k < n → CrossesUpward b.toTournament k := by
  constructor
  · intro h k hk1 hkn
    have hk : k ∈ cutSet n := by
      unfold cutSet
      simp
      omega
    exact (StTiling.isGoodCut_iff_crossesUpward_toTournament b).mp (h k hk)
  · intro h k hk
    have hk1 : 1 ≤ k := by
      unfold cutSet at hk
      simp at hk
      exact hk.1
    have hkn : k < n := by
      unfold cutSet at hk
      simp at hk
      omega
    exact (StTiling.isGoodCut_iff_crossesUpward_toTournament b).mpr
      (h k hk1 hkn)

/-- The full cut set is good exactly when the induced tournament is strongly
    connected. -/
theorem StTiling.goodCuts_eq_cutSet_iff_toTournament_stronglyConnected
    (b : StTiling n) :
    b.goodCuts = cutSet n ↔ IsStronglyConnected b.toTournament := by
  rw [StTiling.goodCuts_eq_cutSet_iff]
  rw [thm330_SC_iff_all_cuts_crossing b.toTournament
    (StTiling.toTournament_hasBasePath b)]
  exact StTiling.all_goodCuts_iff_all_crossesUpward_toTournament b

/-- Top good-cut bucket iff strong connectivity of the induced tournament. -/
theorem StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected
    (b : StTiling n) :
    b.goodCutCount = n - 1 ↔ IsStronglyConnected b.toTournament := by
  rw [StTiling.goodCutCount_eq_top_iff]
  exact StTiling.goodCuts_eq_cutSet_iff_toTournament_stronglyConnected b

theorem StTiling.toTournament_stronglyConnected_iff_goodCutCount_eq_top
    (b : StTiling n) :
    IsStronglyConnected b.toTournament ↔ b.goodCutCount = n - 1 :=
  (StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected b).symm

theorem StTiling.toTournament_not_stronglyConnected_iff_goodCutCount_ne_top
    (b : StTiling n) :
    ¬ IsStronglyConnected b.toTournament ↔ b.goodCutCount ≠ n - 1 := by
  rw [← StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected b]

/-- The all-up tiling is strongly connected once there are legal
    non-consecutive tiles crossing every cut. -/
theorem StTiling.allUp_toTournament_stronglyConnected (hn : 3 ≤ n) :
    IsStronglyConnected (StTiling.allUp n).toTournament := by
  exact (StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected
    (StTiling.allUp n)).mp (StTiling.goodCutCount_allUp hn)

/-- The complement of any all-down tiling is strongly connected for `n ≥ 3`. -/
theorem StTiling.complement_of_allDown_toTournament_stronglyConnected
    {b : StTiling n} (hn : 3 ≤ n)
    (h : ∀ t : StTile n, b t = false) :
    IsStronglyConnected b.complement.toTournament := by
  exact (StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected
    b.complement).mp (StTiling.goodCutCount_complement_of_all_down hn h)

/-- The all-down tiling is not strongly connected as soon as the base path has
    a genuine source/sink separation. -/
theorem StTiling.allDown_toTournament_not_stronglyConnected (hn : 2 ≤ n) :
    ¬ IsStronglyConnected (StTiling.allDown n).toTournament := by
  rw [StTiling.toTournament_not_stronglyConnected_iff_goodCutCount_ne_top]
  have hzero : (StTiling.allDown n).goodCutCount = 0 := by
    rw [StTiling.goodCutCount_eq_zero_iff_all_down]
    intro t
    simp
  rw [hzero]
  omega

theorem StTiling.exists_stronglyConnected_toTournament (hn : 3 ≤ n) :
    ∃ b : StTiling n, IsStronglyConnected b.toTournament :=
  ⟨StTiling.allUp n, StTiling.allUp_toTournament_stronglyConnected hn⟩

theorem StTiling.exists_not_stronglyConnected_toTournament (hn : 2 ≤ n) :
    ∃ b : StTiling n, ¬ IsStronglyConnected b.toTournament :=
  ⟨StTiling.allDown n, StTiling.allDown_toTournament_not_stronglyConnected hn⟩

end

end Tournament
