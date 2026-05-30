/-
  TournamentH7.StaircaseBucketTransport -- concrete staircase transport checksums

  This module specializes the abstract Boolean-cube bucket transport layer to
  the concrete staircase tiling cube `StTiling n = StTile n -> Bool`.
-/

import TournamentH7.BucketBalance
import TournamentH7.ApexBridge

namespace Tournament

variable {n : ℕ}

namespace StTiling

noncomputable section

/-! ### Concrete staircase mask families -/

/-- A staircase mask is nonzero exactly when at least one tile bit is true. -/
theorem isNonzeroMask_iff_exists_true (u : StTiling n) :
    BucketBalance.IsNonzeroMask u ↔ ∃ t : StTile n, u t = true :=
  Iff.rfl

/-- A one-tile flip is a nonzero Boolean mask. -/
theorem singleUp_isNonzeroMask (t : StTile n) :
    BucketBalance.IsNonzeroMask (StTiling.singleUp t) := by
  exact ⟨t, StTiling.singleUp_apply_self t⟩

/-- For `n >= 3`, the complement mask is nonzero. -/
theorem allUp_isNonzeroMask_of_three_le (hn : 3 ≤ n) :
    BucketBalance.IsNonzeroMask (StTiling.allUp n) := by
  exact ⟨apexTile n hn, rfl⟩

/-- The finite family of all nonzero staircase xor masks. -/
def nonzeroMasks (n : ℕ) : Finset (StTiling n) := by
  classical
  exact Finset.univ.filter fun u => BucketBalance.IsNonzeroMask u

theorem mem_nonzeroMasks {u : StTiling n} :
    u ∈ nonzeroMasks n ↔ BucketBalance.IsNonzeroMask u := by
  classical
  simp [nonzeroMasks]

theorem nonzeroMasks_all_nonzero (n : ℕ) :
    ∀ u, u ∈ nonzeroMasks n -> BucketBalance.IsNonzeroMask u := by
  intro u hu
  exact (mem_nonzeroMasks (n := n)).mp hu

/-- The finite family of single-tile staircase xor masks. -/
def singleTileMasks (n : ℕ) : Finset (StTiling n) := by
  classical
  exact Finset.univ.image StTiling.singleUp

theorem mem_singleTileMasks_nonzero {u : StTiling n}
    (hu : u ∈ singleTileMasks n) : BucketBalance.IsNonzeroMask u := by
  classical
  rcases Finset.mem_image.mp hu with ⟨t, _ht, rfl⟩
  exact singleUp_isNonzeroMask t

/-- The singleton family containing the complement mask. -/
def complementMask (n : ℕ) : Finset (StTiling n) :=
  {StTiling.allUp n}

theorem complementMask_all_nonzero (hn : 3 ≤ n) :
    ∀ u, u ∈ complementMask n -> BucketBalance.IsNonzeroMask u := by
  intro u hu
  have hu_eq : u = StTiling.allUp n := by
    simpa [complementMask] using hu
  subst u
  exact allUp_isNonzeroMask_of_three_le hn

/-! ### Concrete unordered and transport-row balances -/

theorem unordered_balance_masks
    {beta : Type*} [DecidableEq beta]
    (q : StTiling n -> beta) (moves : Finset (StTiling n)) (b : beta)
    (hmoves : ∀ u, u ∈ moves -> BucketBalance.IsNonzeroMask u) :
    2 * BucketBalance.internalLineCount q BucketBalance.xorMask moves b +
        (BucketBalance.crossHalf q BucketBalance.xorMask moves b).card =
      (BucketBalance.fiber q b).card * moves.card := by
  exact BucketBalance.unordered_balance_boolCube_masks q moves b hmoves

theorem transport_row_balance_masks
    {beta : Type*} [Fintype beta] [DecidableEq beta]
    (q : StTiling n -> beta) (moves : Finset (StTiling n)) (b : beta)
    (hmoves : ∀ u, u ∈ moves -> BucketBalance.IsNonzeroMask u) :
    2 * BucketBalance.internalLineCount q BucketBalance.xorMask moves b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf q BucketBalance.xorMask moves b c).card) =
      (BucketBalance.fiber q b).card * moves.card := by
  exact BucketBalance.transport_row_balance_boolCube_masks q moves b hmoves

theorem unordered_balance_allNonzeroMasks
    {beta : Type*} [DecidableEq beta] (q : StTiling n -> beta) (b : beta) :
    2 * BucketBalance.internalLineCount q BucketBalance.xorMask (nonzeroMasks n) b +
        (BucketBalance.crossHalf q BucketBalance.xorMask (nonzeroMasks n) b).card =
      (BucketBalance.fiber q b).card * (nonzeroMasks n).card := by
  exact unordered_balance_masks q (nonzeroMasks n) b (nonzeroMasks_all_nonzero n)

theorem transport_row_balance_allNonzeroMasks
    {beta : Type*} [Fintype beta] [DecidableEq beta]
    (q : StTiling n -> beta) (b : beta) :
    2 * BucketBalance.internalLineCount q BucketBalance.xorMask (nonzeroMasks n) b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf q BucketBalance.xorMask (nonzeroMasks n) b c).card) =
      (BucketBalance.fiber q b).card * (nonzeroMasks n).card := by
  exact transport_row_balance_masks q (nonzeroMasks n) b (nonzeroMasks_all_nonzero n)

theorem unordered_balance_singleTileMasks
    {beta : Type*} [DecidableEq beta] (q : StTiling n -> beta) (b : beta) :
    2 * BucketBalance.internalLineCount q BucketBalance.xorMask (singleTileMasks n) b +
        (BucketBalance.crossHalf q BucketBalance.xorMask (singleTileMasks n) b).card =
      (BucketBalance.fiber q b).card * (singleTileMasks n).card := by
  refine unordered_balance_masks q (singleTileMasks n) b ?_
  intro u hu
  exact mem_singleTileMasks_nonzero hu

theorem transport_row_balance_singleTileMasks
    {beta : Type*} [Fintype beta] [DecidableEq beta]
    (q : StTiling n -> beta) (b : beta) :
    2 * BucketBalance.internalLineCount q BucketBalance.xorMask (singleTileMasks n) b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf q BucketBalance.xorMask (singleTileMasks n) b c).card) =
      (BucketBalance.fiber q b).card * (singleTileMasks n).card := by
  refine transport_row_balance_masks q (singleTileMasks n) b ?_
  intro u hu
  exact mem_singleTileMasks_nonzero hu

theorem unordered_balance_complementMask
    {beta : Type*} [DecidableEq beta]
    (q : StTiling n -> beta) (b : beta) (hn : 3 ≤ n) :
    2 * BucketBalance.internalLineCount q BucketBalance.xorMask (complementMask n) b +
        (BucketBalance.crossHalf q BucketBalance.xorMask (complementMask n) b).card =
      (BucketBalance.fiber q b).card * (complementMask n).card := by
  exact unordered_balance_masks q (complementMask n) b (complementMask_all_nonzero hn)

theorem transport_row_balance_complementMask
    {beta : Type*} [Fintype beta] [DecidableEq beta]
    (q : StTiling n -> beta) (b : beta) (hn : 3 ≤ n) :
    2 * BucketBalance.internalLineCount q BucketBalance.xorMask (complementMask n) b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf q BucketBalance.xorMask (complementMask n) b c).card) =
      (BucketBalance.fiber q b).card * (complementMask n).card := by
  exact transport_row_balance_masks q (complementMask n) b (complementMask_all_nonzero hn)

/-! ### The concrete good-cut quotient -/

/-- Good-cut count as a finite bucket target. -/
def goodCutBucket (u : StTiling n) : Fin (n + 1) :=
  ⟨u.goodCutCount, by
    have h := StTiling.goodCutCount_le_n_minus_one u
    omega⟩

/-- The top good-cut bucket, represented in `Fin (n+1)`. -/
def topGoodCutBucket (n : ℕ) : Fin (n + 1) :=
  ⟨n - 1, by omega⟩

theorem goodCutBucket_eq_zero_iff_all_down (u : StTiling n) :
    goodCutBucket u = (0 : Fin (n + 1)) ↔ ∀ t : StTile n, u t = false := by
  constructor
  · intro h
    have hval : u.goodCutCount = 0 := congrArg Fin.val h
    exact (StTiling.goodCutCount_eq_zero_iff_all_down u).mp hval
  · intro h
    apply Fin.ext
    exact (StTiling.goodCutCount_eq_zero_iff_all_down u).mpr h

theorem goodCutBucket_eq_top_iff_toTournament_SC (u : StTiling n) :
    goodCutBucket u = topGoodCutBucket n ↔ IsStronglyConnected u.toTournament := by
  constructor
  · intro h
    have hval : u.goodCutCount = n - 1 := congrArg Fin.val h
    exact (StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected u).mp hval
  · intro h
    apply Fin.ext
    exact (StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected u).mpr h

/-! ### Good-cut quotient gaps -/

theorem goodCutBucket_image_iff (hn : 3 ≤ n) (b : Fin (n + 1)) :
    (∃ u : StTiling n, goodCutBucket u = b) ↔
      b.val = 0 ∨ (2 ≤ b.val ∧ b.val ≤ n - 1) := by
  constructor
  · rintro ⟨u, hu⟩
    have hval : u.goodCutCount = b.val := congrArg Fin.val hu
    exact (StTiling.goodCutCount_spectrum (n := n) (r := b.val) hn).mp
      ⟨u, hval⟩
  · intro hb
    rcases (StTiling.goodCutCount_spectrum (n := n) (r := b.val) hn).mpr hb
      with ⟨u, hu⟩
    refine ⟨u, ?_⟩
    apply Fin.ext
    exact hu

theorem goodCutBucket_ne_one (hn : 1 ≤ n) (u : StTiling n) :
    goodCutBucket u ≠ (⟨1, by omega⟩ : Fin (n + 1)) := by
  intro h
  have hval : u.goodCutCount = 1 := congrArg Fin.val h
  exact StTiling.goodCutCount_ne_one u hval

theorem goodCutBucket_ne_overTop (hn : 1 ≤ n) (u : StTiling n) :
    goodCutBucket u ≠ (⟨n, by omega⟩ : Fin (n + 1)) := by
  intro h
  have hval : u.goodCutCount = n := congrArg Fin.val h
  have hle := StTiling.goodCutCount_le_n_minus_one u
  omega

theorem goodCutBucket_fiber_one_eq_empty (hn : 1 ≤ n) :
    BucketBalance.fiber goodCutBucket (⟨1, by omega⟩ : Fin (n + 1)) = ∅ := by
  classical
  ext u
  simp [BucketBalance.fiber, goodCutBucket_ne_one hn u]

theorem goodCutBucket_fiber_overTop_eq_empty (hn : 1 ≤ n) :
    BucketBalance.fiber goodCutBucket (⟨n, by omega⟩ : Fin (n + 1)) = ∅ := by
  classical
  ext u
  simp [BucketBalance.fiber, goodCutBucket_ne_overTop hn u]

theorem transport_row_balance_goodCutBucket_allNonzeroMasks (b : Fin (n + 1)) :
    2 * BucketBalance.internalLineCount goodCutBucket BucketBalance.xorMask
        (nonzeroMasks n) b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf goodCutBucket BucketBalance.xorMask
            (nonzeroMasks n) b c).card) =
      (BucketBalance.fiber goodCutBucket b).card * (nonzeroMasks n).card := by
  exact transport_row_balance_allNonzeroMasks goodCutBucket b

theorem transport_row_balance_goodCutBucket_singleTileMasks (b : Fin (n + 1)) :
    2 * BucketBalance.internalLineCount goodCutBucket BucketBalance.xorMask
        (singleTileMasks n) b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf goodCutBucket BucketBalance.xorMask
            (singleTileMasks n) b c).card) =
      (BucketBalance.fiber goodCutBucket b).card * (singleTileMasks n).card := by
  exact transport_row_balance_singleTileMasks goodCutBucket b

theorem transport_row_balance_goodCutBucket_complementMask
    (b : Fin (n + 1)) (hn : 3 ≤ n) :
    2 * BucketBalance.internalLineCount goodCutBucket BucketBalance.xorMask
        (complementMask n) b +
        (∑ c ∈ (Finset.univ.erase b),
          (BucketBalance.transportHalf goodCutBucket BucketBalance.xorMask
            (complementMask n) b c).card) =
      (BucketBalance.fiber goodCutBucket b).card * (complementMask n).card := by
  exact transport_row_balance_complementMask goodCutBucket b hn

end

end StTiling
end Tournament
