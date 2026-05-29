/-
  TournamentH7.GoodCuts — good-cut buckets for staircase tilings

  This module formalises the first layer of THM-336 in the concrete
  `StTiling` model:

   • a tile crosses a cut when its interval spans that cut;
   • a cut is good for a tiling when some upward tile crosses it;
   • the good-cut bucket `goodCutCount` is never equal to 1;
   • bucket 0 is exactly the all-down/transitive tiling;
   • grid reflection preserves good-cut bucket size.

  The point is to make the "bucket constraints" available to Lean without
  depending on any project axioms.
-/

import TournamentH7.StaircaseTileModel
import Mathlib.Tactic

namespace Tournament

variable {n : ℕ}

noncomputable section

/-! ### Cuts and crossing -/

/-- The legal cuts in the base-path staircase: `1, ..., n-1`. -/
def cutSet (n : ℕ) : Finset ℕ := Finset.Icc 1 (n - 1)

/-- A tile `(hi, lo)` crosses cut `k` iff `lo < k ≤ hi`. -/
def StTile.crossesCut (t : StTile n) (k : ℕ) : Prop :=
  t.lo.val < k ∧ k ≤ t.hi.val

namespace StTile

/-- Value of the high endpoint after grid reflection. -/
@[simp] lemma reflect_hi_val (t : StTile n) :
    t.reflect.hi.val = n - 1 - t.lo.val := rfl

/-- Value of the low endpoint after grid reflection. -/
@[simp] lemma reflect_lo_val (t : StTile n) :
    t.reflect.lo.val = n - 1 - t.hi.val := rfl

/-- The first cut crossed by a tile. -/
lemma crossesCut_lo_succ (t : StTile n) :
    t.crossesCut (t.lo.val + 1) := by
  have h_gap := t.gap
  unfold crossesCut
  constructor <;> omega

/-- The second cut crossed by a tile. This is where non-consecutiveness enters. -/
lemma crossesCut_lo_succ_succ (t : StTile n) :
    t.crossesCut (t.lo.val + 2) := by
  have h_gap := t.gap
  unfold crossesCut
  constructor <;> omega

lemma lo_succ_mem_cutSet (t : StTile n) :
    t.lo.val + 1 ∈ cutSet n := by
  unfold cutSet
  have h_hi : t.hi.val < n := t.hi.is_lt
  have h_gap := t.gap
  simp
  omega

lemma lo_succ_succ_mem_cutSet (t : StTile n) :
    t.lo.val + 2 ∈ cutSet n := by
  unfold cutSet
  have h_hi : t.hi.val < n := t.hi.is_lt
  have h_gap := t.gap
  simp
  omega

lemma cutSet_reflect_mem {k : ℕ} (hk : k ∈ cutSet n) :
    n - k ∈ cutSet n := by
  unfold cutSet at hk ⊢
  simp at hk ⊢
  omega

lemma reflect_reflect_cut {k : ℕ} (hk : k ∈ cutSet n) :
    n - (n - k) = k := by
  unfold cutSet at hk
  simp at hk
  omega

/-- Grid-reflection sends cut `k` to cut `n-k`.
    The equivalence is arithmetically true for every natural `k`; outside
    the legal cut range both sides are false. -/
lemma reflect_crossesCut {t : StTile n} {k : ℕ} :
    t.reflect.crossesCut (n - k) ↔ t.crossesCut k := by
  simp [StTile.crossesCut]
  constructor <;> intro h <;> constructor <;> omega

end StTile

/-! ### Good-cut buckets -/

/-- A good cut is one crossed by at least one upward tile. -/
def StTiling.IsGoodCut (b : StTiling n) (k : ℕ) : Prop :=
  ∃ t : StTile n, b t = true ∧ t.crossesCut k

/-- The finite set of good cuts of a tiling. -/
def StTiling.goodCuts (b : StTiling n) : Finset ℕ := by
  classical
  exact (cutSet n).filter (fun k => StTiling.IsGoodCut b k)

/-- The good-cut bucket index of a tiling. -/
def StTiling.goodCutCount (b : StTiling n) : ℕ :=
  b.goodCuts.card

lemma StTiling.mem_goodCuts {b : StTiling n} {k : ℕ} :
    k ∈ b.goodCuts ↔ k ∈ cutSet n ∧ StTiling.IsGoodCut b k := by
  unfold StTiling.goodCuts
  simp

lemma StTiling.goodCut_of_upward_tile {b : StTiling n} {t : StTile n}
    (ht : b t = true) :
    t.lo.val + 1 ∈ b.goodCuts ∧ t.lo.val + 2 ∈ b.goodCuts := by
  constructor
  · rw [StTiling.mem_goodCuts]
    exact ⟨t.lo_succ_mem_cutSet, ⟨t, ht, t.crossesCut_lo_succ⟩⟩
  · rw [StTiling.mem_goodCuts]
    exact ⟨t.lo_succ_succ_mem_cutSet, ⟨t, ht, t.crossesCut_lo_succ_succ⟩⟩

/-- Bucket 0: if all tiles are down, there are no good cuts. -/
theorem StTiling.goodCuts_empty_of_all_down {b : StTiling n}
    (h : ∀ t : StTile n, b t = false) :
    b.goodCuts = ∅ := by
  classical
  ext k
  rw [StTiling.mem_goodCuts]
  simp [StTiling.IsGoodCut, h]

/-- Bucket 0 is exactly the all-down tiling. -/
theorem StTiling.all_down_of_goodCuts_empty {b : StTiling n}
    (h : b.goodCuts = ∅) :
    ∀ t : StTile n, b t = false := by
  intro t
  cases ht : b t with
  | false => rfl
  | true =>
      have hg := StTiling.goodCut_of_upward_tile (b := b) (t := t) ht
      have : t.lo.val + 1 ∈ (∅ : Finset ℕ) := by
        rw [← h]
        exact hg.1
      simp at this

theorem StTiling.goodCuts_empty_iff_all_down (b : StTiling n) :
    b.goodCuts = ∅ ↔ ∀ t : StTile n, b t = false :=
  ⟨StTiling.all_down_of_goodCuts_empty, StTiling.goodCuts_empty_of_all_down⟩

/-- THM-336, Lean core: no tiling lives in good-cut bucket 1. -/
theorem StTiling.goodCutCount_ne_one (b : StTiling n) :
    b.goodCutCount ≠ 1 := by
  unfold StTiling.goodCutCount
  intro h
  rcases Finset.card_eq_one.mp h with ⟨k, hk⟩
  have hk_mem : k ∈ b.goodCuts := by
    rw [hk]
    simp
  rw [StTiling.mem_goodCuts] at hk_mem
  rcases hk_mem with ⟨_, t, ht, _⟩
  have htwo := StTiling.goodCut_of_upward_tile (b := b) (t := t) ht
  have h1 : t.lo.val + 1 = k := by
    have := htwo.1
    rw [hk] at this
    simpa using this
  have h2 : t.lo.val + 2 = k := by
    have := htwo.2
    rw [hk] at this
    simpa using this
  omega

/-! ### Grid-reflection preserves buckets -/

private lemma StTiling.reflect_goodCut_forward {b : StTiling n} {k : ℕ}
    (hk : k ∈ cutSet n) :
    k ∈ b.goodCuts → n - k ∈ b.reflect.goodCuts := by
  intro hgood
  rw [StTiling.mem_goodCuts] at hgood ⊢
  rcases hgood with ⟨_, t, ht, hcross⟩
  refine ⟨StTile.cutSet_reflect_mem hk, ⟨t.reflect, ?_, ?_⟩⟩
  · unfold StTiling.reflect
    rw [StTile.reflect_reflect]
    exact ht
  · rw [StTile.reflect_crossesCut]
    exact hcross

/-- Membership form: grid-reflection sends the good cut `k` to `n-k`. -/
theorem StTiling.mem_goodCuts_reflect {b : StTiling n} {k : ℕ}
    (hk : k ∈ cutSet n) :
    n - k ∈ b.reflect.goodCuts ↔ k ∈ b.goodCuts := by
  constructor
  · intro h
    have hk' : n - k ∈ cutSet n := StTile.cutSet_reflect_mem hk
    have h2 := StTiling.reflect_goodCut_forward (b := b.reflect) (k := n - k) hk' h
    rw [StTiling.reflect_reflect] at h2
    have hkk : n - (n - k) = k := StTile.reflect_reflect_cut hk
    rwa [hkk] at h2
  · exact StTiling.reflect_goodCut_forward (b := b) (k := k) hk

/-- Grid-reflection preserves the good-cut bucket index. -/
theorem StTiling.goodCutCount_reflect (b : StTiling n) :
    b.reflect.goodCutCount = b.goodCutCount := by
  unfold StTiling.goodCutCount
  symm
  apply Finset.card_bij (fun k _ => n - k)
  · intro k hk
    have hk_cut : k ∈ cutSet n := (StTiling.mem_goodCuts.mp hk).1
    exact (StTiling.mem_goodCuts_reflect (b := b) hk_cut).2 hk
  · intro a ha c hc h
    have ha_cut : a ∈ cutSet n := (StTiling.mem_goodCuts.mp ha).1
    have hc_cut : c ∈ cutSet n := (StTiling.mem_goodCuts.mp hc).1
    unfold cutSet at ha_cut hc_cut
    simp at ha_cut hc_cut
    omega
  · intro k hk
    have hk_cut : k ∈ cutSet n := (StTiling.mem_goodCuts.mp hk).1
    refine ⟨n - k, ?_, ?_⟩
    · have hnk_cut : n - k ∈ cutSet n := StTile.cutSet_reflect_mem hk_cut
      exact (StTiling.mem_goodCuts_reflect (b := b) hnk_cut).1 (by
        rwa [StTile.reflect_reflect_cut hk_cut])
    · exact StTile.reflect_reflect_cut hk_cut

end

end Tournament
