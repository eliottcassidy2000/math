/-
  TournamentH7.ApexBridge — Bridge from tile apex to tournament SC

  Connects:
   • The staircase apex tile (n-1, 0) with bit = true (upward).
   • The induced base-path tournament being strongly connected.
-/

import TournamentH7.StaircaseConnectivity
import TournamentH7.GoodCuts

namespace Tournament

variable {n : ℕ}

/-! ### The apex tile -/

/-- The apex tile (n-1, 0) for n ≥ 3. -/
noncomputable def apexTile (n : ℕ) (hn : 3 ≤ n) : StTile n where
  hi := ⟨n - 1, by omega⟩
  lo := ⟨0, by omega⟩
  gap := by show 0 + 2 ≤ n - 1; omega

@[simp] lemma apexTile_hi (n : ℕ) (hn : 3 ≤ n) :
    (apexTile n hn).hi.val = n - 1 := rfl

@[simp] lemma apexTile_lo (n : ℕ) (hn : 3 ≤ n) :
    (apexTile n hn).lo.val = 0 := rfl

/-- The apex tile crosses every legal base-path cut. -/
theorem apexTile_cutInterval_eq_cutSet (n : ℕ) (hn : 3 ≤ n) :
    (apexTile n hn).cutInterval = cutSet n := by
  ext k
  rw [StTile.mem_cutInterval]
  constructor
  · exact And.left
  · intro hk
    refine ⟨hk, ?_⟩
    unfold StTile.crossesCut
    unfold cutSet at hk
    simp at hk
    constructor
    · rw [apexTile_lo]
      omega
    · rw [apexTile_hi]
      omega

/-! ### The "singleUp apex" tiling is SC

    Setting only the apex tile to UP (and others to DOWN) gives a tiling
    whose induced tournament is strongly connected, because the apex
    crosses every cut. -/

/-- A tiling with the apex tile UP has the apex arc 0 → (n-1) in its
    induced tournament. -/
theorem singleUp_apex_has_apex_arc (n : ℕ) (hn : 3 ≤ n) :
    (StTiling.singleUp (apexTile n hn)).arc
      ⟨0, by omega⟩ ⟨n - 1, by omega⟩ = true := by
  have h_apex_gap : (0 : ℕ) + 2 ≤ n - 1 := by omega
  have h_eq := StTiling.arc_nonconsec_up (b := StTiling.singleUp (apexTile n hn))
    (lo := ⟨0, by omega⟩) (hi := ⟨n - 1, by omega⟩) h_apex_gap
  rw [h_eq]
  show (StTiling.singleUp (apexTile n hn)) _ = true
  -- The tile is { hi := ⟨n-1,_⟩, lo := ⟨0,_⟩, gap := _ } = apexTile n hn (up to gap proof).
  -- Apply StTiling.singleUp_eq_true_iff with s = this tile, t = apexTile n hn.
  rw [StTiling.singleUp_eq_true_iff]
  -- Need to show: { hi := ⟨n-1, _⟩, lo := ⟨0, _⟩, gap := h_apex_gap } = apexTile n hn.
  apply StTile.ext <;> rfl

/-- The tiling with only the apex tile UP induces a strongly connected
    tournament.  PROVED. -/
theorem singleUp_apex_toTournament_SC (n : ℕ) (hn : 3 ≤ n) :
    IsStronglyConnected (StTiling.singleUp (apexTile n hn)).toTournament := by
  -- Apply apex_implies_SC.
  have h_apex := singleUp_apex_has_apex_arc n hn
  have h_bp : HasBasePath (StTiling.singleUp (apexTile n hn)).toTournament :=
    StTiling.toTournament_hasBasePath _
  exact apex_implies_SC _ h_bp hn (by omega) (by omega) h_apex

/-- The single-up apex tiling has every legal cut good. -/
theorem singleUp_apex_goodCuts_eq_cutSet (n : ℕ) (hn : 3 ≤ n) :
    (StTiling.singleUp (apexTile n hn)).goodCuts = cutSet n := by
  rw [StTiling.goodCuts_singleUp_eq_cutInterval]
  exact apexTile_cutInterval_eq_cutSet n hn

/-- The single-up apex tiling lies in the top good-cut bucket. -/
theorem singleUp_apex_goodCutCount_top (n : ℕ) (hn : 3 ≤ n) :
    (StTiling.singleUp (apexTile n hn)).goodCutCount = n - 1 := by
  unfold StTiling.goodCutCount
  rw [singleUp_apex_goodCuts_eq_cutSet n hn, cutSet_card]

end Tournament
