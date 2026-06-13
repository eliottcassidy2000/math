/-
  TournamentH7.GoodCuts — good-cut buckets for staircase tilings

  This module formalises the first layer of THM-336 in the concrete
  `StTiling` model:

   • a tile crosses a cut when its interval spans that cut;
   • a cut is good for a tiling when some upward tile crosses it;
   • the good-cut bucket `goodCutCount` is either 0 or at least 2;
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

@[simp] lemma cutSet_card (n : ℕ) :
    (cutSet n).card = n - 1 := by
  unfold cutSet
  simp

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

lemma exists_crossesCut_of_mem_cutSet (hn : 3 ≤ n) {k : ℕ}
    (hk : k ∈ cutSet n) :
    ∃ t : StTile n, t.crossesCut k := by
  unfold cutSet at hk
  simp at hk
  by_cases hk1 : k = 1
  · subst hk1
    have h2n : 2 < n := by omega
    have h0n : 0 < n := by omega
    let t : StTile n := ⟨⟨2, h2n⟩, ⟨0, h0n⟩, by norm_num⟩
    refine ⟨t, ?_⟩
    unfold StTile.crossesCut
    simp [t]
  · have hk_ge : 2 ≤ k := by omega
    have hk_lt : k < n := by omega
    have hlo_lt : k - 2 < n := by omega
    have hgap : (k - 2) + 2 ≤ k := by omega
    let t : StTile n := ⟨⟨k, hk_lt⟩, ⟨k - 2, hlo_lt⟩, hgap⟩
    refine ⟨t, ?_⟩
    unfold StTile.crossesCut
    simp [t]
    omega

lemma reflect_reflect_cut {k : ℕ} (hk : k ∈ cutSet n) :
    n - (n - k) = k := by
  unfold cutSet at hk
  simp at hk
  omega

lemma crossesCut_mem_cutSet {t : StTile n} {k : ℕ}
    (h : t.crossesCut k) :
    k ∈ cutSet n := by
  unfold cutSet
  have h_hi : t.hi.val < n := t.hi.is_lt
  unfold crossesCut at h
  simp
  omega

/-- The interval of cuts crossed by a tile. -/
def cutInterval (t : StTile n) : Finset ℕ := by
  classical
  exact (cutSet n).filter (fun k => t.crossesCut k)

lemma mem_cutInterval {t : StTile n} {k : ℕ} :
    k ∈ t.cutInterval ↔ k ∈ cutSet n ∧ t.crossesCut k := by
  unfold cutInterval
  classical
  simp

/-- Grid-reflection sends cut `k` to cut `n-k`.
    The equivalence is arithmetically true for every natural `k`; outside
    the legal cut range both sides are false. -/
lemma reflect_crossesCut {t : StTile n} {k : ℕ} :
    t.reflect.crossesCut (n - k) ↔ t.crossesCut k := by
  simp [StTile.crossesCut]
  constructor <;> intro h <;> constructor <;> omega

end StTile

/-! ### Good-cut buckets -/

/-- The all-down tiling. -/
def StTiling.allDown (n : ℕ) : StTiling n := fun _ => false

/-- The all-up tiling. -/
def StTiling.allUp (n : ℕ) : StTiling n := fun _ => true

@[simp] lemma StTiling.allDown_apply (t : StTile n) :
    StTiling.allDown n t = false := rfl

@[simp] lemma StTiling.allUp_apply (t : StTile n) :
    StTiling.allUp n t = true := rfl

@[simp] lemma StTiling.complement_allDown :
    (StTiling.allDown n).complement = StTiling.allUp n := by
  funext t
  simp [StTiling.complement, StTiling.allDown, StTiling.allUp]

@[simp] lemma StTiling.complement_allUp :
    (StTiling.allUp n).complement = StTiling.allDown n := by
  funext t
  simp [StTiling.complement, StTiling.allDown, StTiling.allUp]

/-- The tiling with exactly one named tile upward. -/
def StTiling.singleUp (t : StTile n) : StTiling n :=
  fun s => decide (s.hi = t.hi ∧ s.lo = t.lo)

@[simp] lemma StTiling.singleUp_apply_self (t : StTile n) :
    StTiling.singleUp t t = true := by
  simp [StTiling.singleUp]

lemma StTiling.singleUp_eq_true_iff {s t : StTile n} :
    StTiling.singleUp t s = true ↔ s = t := by
  constructor
  · intro hs
    have hpair : s.hi = t.hi ∧ s.lo = t.lo := by
      simpa [StTiling.singleUp] using hs
    exact StTile.ext hpair.1 hpair.2
  · intro h
    subst h
    simp [StTiling.singleUp]

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

/-- A cut is good exactly when it lies in the cut interval of some upward tile. -/
theorem StTiling.isGoodCut_iff_exists_upward_tile_interval
    {b : StTiling n} {k : ℕ} :
    StTiling.IsGoodCut b k ↔
      ∃ t : StTile n, b t = true ∧ k ∈ t.cutInterval := by
  constructor
  · rintro ⟨t, ht, hcross⟩
    exact ⟨t, ht, (StTile.mem_cutInterval).2
      ⟨StTile.crossesCut_mem_cutSet hcross, hcross⟩⟩
  · rintro ⟨t, ht, hk⟩
    exact ⟨t, ht, (StTile.mem_cutInterval.mp hk).2⟩

theorem StTiling.mem_goodCuts_iff_exists_upward_tile_interval
    {b : StTiling n} {k : ℕ} :
    k ∈ b.goodCuts ↔
      ∃ t : StTile n, b t = true ∧ k ∈ t.cutInterval := by
  constructor
  · intro hk
    exact StTiling.isGoodCut_iff_exists_upward_tile_interval.mp
      (StTiling.mem_goodCuts.mp hk).2
  · rintro ⟨t, ht, hk⟩
    rw [StTiling.mem_goodCuts]
    exact ⟨(StTile.mem_cutInterval.mp hk).1,
      StTiling.isGoodCut_iff_exists_upward_tile_interval.mpr ⟨t, ht, hk⟩⟩

lemma StTiling.cutInterval_subset_goodCuts_of_upward_tile
    {b : StTiling n} {t : StTile n} (ht : b t = true) :
    t.cutInterval ⊆ b.goodCuts := by
  intro k hk
  rw [StTiling.mem_goodCuts_iff_exists_upward_tile_interval]
  exact ⟨t, ht, hk⟩

lemma StTiling.goodCuts_subset_cutSet (b : StTiling n) :
    b.goodCuts ⊆ cutSet n := by
  intro k hk
  exact (StTiling.mem_goodCuts.mp hk).1

lemma StTiling.goodCutCount_le_cutSet_card (b : StTiling n) :
    b.goodCutCount ≤ (cutSet n).card := by
  unfold StTiling.goodCutCount
  exact Finset.card_le_card (StTiling.goodCuts_subset_cutSet b)

lemma StTiling.goodCutCount_le_n_minus_one (b : StTiling n) :
    b.goodCutCount ≤ n - 1 := by
  simpa using StTiling.goodCutCount_le_cutSet_card (b := b)

lemma StTiling.goodCuts_mono {b c : StTiling n}
    (h : ∀ t : StTile n, b t = true → c t = true) :
    b.goodCuts ⊆ c.goodCuts := by
  intro k hk
  rw [StTiling.mem_goodCuts] at hk ⊢
  rcases hk with ⟨hk_cut, t, ht, hcross⟩
  exact ⟨hk_cut, ⟨t, h t ht, hcross⟩⟩

lemma StTiling.goodCutCount_mono {b c : StTiling n}
    (h : ∀ t : StTile n, b t = true → c t = true) :
    b.goodCutCount ≤ c.goodCutCount := by
  unfold StTiling.goodCutCount
  exact Finset.card_le_card (StTiling.goodCuts_mono h)

lemma StTiling.goodCut_of_upward_tile {b : StTiling n} {t : StTile n}
    (ht : b t = true) :
    t.lo.val + 1 ∈ b.goodCuts ∧ t.lo.val + 2 ∈ b.goodCuts := by
  constructor
  · rw [StTiling.mem_goodCuts]
    exact ⟨t.lo_succ_mem_cutSet, ⟨t, ht, t.crossesCut_lo_succ⟩⟩
  · rw [StTiling.mem_goodCuts]
    exact ⟨t.lo_succ_succ_mem_cutSet, ⟨t, ht, t.crossesCut_lo_succ_succ⟩⟩

/-- Any upward tile makes the good-cut set nonempty. -/
theorem StTiling.goodCuts_nonempty_of_upward_tile {b : StTiling n} {t : StTile n}
    (ht : b t = true) :
    b.goodCuts.Nonempty := by
  exact ⟨t.lo.val + 1, (StTiling.goodCut_of_upward_tile (b := b) (t := t) ht).1⟩

/-- A nonempty good-cut set witnesses at least one upward tile. -/
theorem StTiling.exists_upward_tile_of_goodCuts_nonempty {b : StTiling n}
    (h : b.goodCuts.Nonempty) :
    ∃ t : StTile n, b t = true := by
  rcases h with ⟨k, hk⟩
  rw [StTiling.mem_goodCuts] at hk
  rcases hk with ⟨_, t, ht, _⟩
  exact ⟨t, ht⟩

/-- Positive good-cut support is equivalent to the existence of an upward tile. -/
theorem StTiling.goodCuts_nonempty_iff_exists_upward_tile (b : StTiling n) :
    b.goodCuts.Nonempty ↔ ∃ t : StTile n, b t = true := by
  constructor
  · exact StTiling.exists_upward_tile_of_goodCuts_nonempty
  · rintro ⟨t, ht⟩
    exact StTiling.goodCuts_nonempty_of_upward_tile (b := b) (t := t) ht

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

/-- Cardinality form of bucket 0. -/
theorem StTiling.goodCutCount_eq_zero_iff_all_down (b : StTiling n) :
    b.goodCutCount = 0 ↔ ∀ t : StTile n, b t = false := by
  unfold StTiling.goodCutCount
  rw [Finset.card_eq_zero]
  exact StTiling.goodCuts_empty_iff_all_down b

/-- Positive good-cut count is equivalent to the existence of an upward tile. -/
theorem StTiling.goodCutCount_pos_iff_exists_upward_tile (b : StTiling n) :
    0 < b.goodCutCount ↔ ∃ t : StTile n, b t = true := by
  unfold StTiling.goodCutCount
  rw [Finset.card_pos]
  exact StTiling.goodCuts_nonempty_iff_exists_upward_tile b

/-- A tiling has positive good-cut count iff it is not the all-down tiling. -/
theorem StTiling.goodCutCount_pos_iff_not_all_down (b : StTiling n) :
    0 < b.goodCutCount ↔ ¬ ∀ t : StTile n, b t = false := by
  rw [StTiling.goodCutCount_pos_iff_exists_upward_tile]
  constructor
  · rintro ⟨t, ht⟩ hall
    rw [hall t] at ht
    contradiction
  · intro hnot
    by_contra hnone
    apply hnot
    intro t
    cases ht : b t with
    | false => rfl
    | true => exact False.elim (hnone ⟨t, ht⟩)

/-- A single upward tile forces at least two distinct good cuts. -/
theorem StTiling.two_le_goodCutCount_of_upward_tile {b : StTiling n} {t : StTile n}
    (ht : b t = true) :
    2 ≤ b.goodCutCount := by
  classical
  unfold StTiling.goodCutCount
  have htwo := StTiling.goodCut_of_upward_tile (b := b) (t := t) ht
  let pair : Finset ℕ := {t.lo.val + 1, t.lo.val + 2}
  have hsub : pair ⊆ b.goodCuts := by
    intro k hk
    simp [pair] at hk
    rcases hk with hk | hk
    · rw [hk]
      exact htwo.1
    · rw [hk]
      exact htwo.2
  have hcard : pair.card = 2 := by
    simp [pair]
  have hle := Finset.card_le_card hsub
  omega

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

/-- THM-336 strengthened: every tiling has either zero or at least two good cuts. -/
theorem StTiling.goodCutCount_eq_zero_or_two_le (b : StTiling n) :
    b.goodCutCount = 0 ∨ 2 ≤ b.goodCutCount := by
  by_cases h : ∃ t : StTile n, b t = true
  · rcases h with ⟨t, ht⟩
    exact Or.inr (StTiling.two_le_goodCutCount_of_upward_tile (b := b) (t := t) ht)
  · left
    rw [StTiling.goodCutCount_eq_zero_iff_all_down]
    intro t
    cases ht : b t with
    | false => rfl
    | true => exact False.elim (h ⟨t, ht⟩)

/-- Set-cardinality form: the good-cut set is empty or has at least two elements. -/
theorem StTiling.goodCuts_empty_or_two_le_card (b : StTiling n) :
    b.goodCuts = ∅ ∨ 2 ≤ b.goodCuts.card := by
  have h := StTiling.goodCutCount_eq_zero_or_two_le b
  unfold StTiling.goodCutCount at h
  rcases h with h0 | h2
  · exact Or.inl (Finset.card_eq_zero.mp h0)
  · exact Or.inr h2

theorem StTiling.goodCutCount_bucket_bounds (b : StTiling n) :
    b.goodCutCount = 0 ∨
      (2 ≤ b.goodCutCount ∧ b.goodCutCount ≤ n - 1) := by
  rcases StTiling.goodCutCount_eq_zero_or_two_le b with h0 | h2
  · exact Or.inl h0
  · exact Or.inr ⟨h2, StTiling.goodCutCount_le_n_minus_one b⟩

/-- The top bucket is exactly the full cut set. -/
theorem StTiling.goodCuts_eq_cutSet_iff (b : StTiling n) :
    b.goodCuts = cutSet n ↔
      ∀ k, k ∈ cutSet n → StTiling.IsGoodCut b k := by
  constructor
  · intro h k hk
    have hk_good : k ∈ b.goodCuts := by
      rw [h]
      exact hk
    exact (StTiling.mem_goodCuts.mp hk_good).2
  · intro h
    apply Finset.Subset.antisymm
    · exact StTiling.goodCuts_subset_cutSet b
    · intro k hk
      rw [StTiling.mem_goodCuts]
      exact ⟨hk, h k hk⟩

theorem StTiling.goodCutCount_eq_cutSet_card_iff (b : StTiling n) :
    b.goodCutCount = (cutSet n).card ↔ b.goodCuts = cutSet n := by
  unfold StTiling.goodCutCount
  constructor
  · intro h
    apply Finset.eq_of_subset_of_card_le (StTiling.goodCuts_subset_cutSet b)
    omega
  · intro h
    rw [h]

theorem StTiling.goodCutCount_eq_top_iff (b : StTiling n) :
    b.goodCutCount = n - 1 ↔ b.goodCuts = cutSet n := by
  rw [← cutSet_card n]
  exact StTiling.goodCutCount_eq_cutSet_card_iff b

theorem StTiling.goodCutCount_eq_top_iff_all_cuts_good (b : StTiling n) :
    b.goodCutCount = n - 1 ↔
      ∀ k, k ∈ cutSet n → StTiling.IsGoodCut b k := by
  rw [StTiling.goodCutCount_eq_top_iff,
    StTiling.goodCuts_eq_cutSet_iff]

theorem StTiling.goodCuts_singleUp_eq_cutInterval (t : StTile n) :
    (StTiling.singleUp t).goodCuts = t.cutInterval := by
  ext k
  rw [StTiling.mem_goodCuts, StTile.mem_cutInterval]
  constructor
  · rintro ⟨hk, s, hs, hcross⟩
    have hst : s = t := StTiling.singleUp_eq_true_iff.mp hs
    subst s
    exact ⟨hk, hcross⟩
  · rintro ⟨hk, hcross⟩
    exact ⟨hk, t, by simp, hcross⟩

theorem StTiling.exists_goodCutCount_eq_of_allowed
    (hn : 3 ≤ n) {r : ℕ} (hr2 : 2 ≤ r) (hrn : r ≤ n - 1) :
    ∃ b : StTiling n, b.goodCutCount = r := by
  have hr_lt : r < n := by omega
  have h0_lt : 0 < n := by omega
  let t : StTile n := ⟨⟨r, hr_lt⟩, ⟨0, h0_lt⟩, by omega⟩
  refine ⟨StTiling.singleUp t, ?_⟩
  unfold StTiling.goodCutCount
  rw [StTiling.goodCuts_singleUp_eq_cutInterval]
  have hinterval : t.cutInterval = Finset.Icc 1 r := by
    ext k
    rw [StTile.mem_cutInterval]
    unfold cutSet StTile.crossesCut
    simp [t]
    constructor <;> intro h <;> omega
  rw [hinterval]
  simp

theorem StTiling.exists_goodCutCount_eq_zero (n : ℕ) :
    ∃ b : StTiling n, b.goodCutCount = 0 := by
  refine ⟨StTiling.allDown n, ?_⟩
  rw [StTiling.goodCutCount_eq_zero_iff_all_down]
  intro t
  simp

theorem StTiling.goodCutCount_spectrum (hn : 3 ≤ n) {r : ℕ} :
    (∃ b : StTiling n, b.goodCutCount = r) ↔
      r = 0 ∨ (2 ≤ r ∧ r ≤ n - 1) := by
  constructor
  · rintro ⟨b, hb⟩
    rcases StTiling.goodCutCount_bucket_bounds b with h0 | hpos
    · left
      omega
    · right
      omega
  · intro h
    rcases h with h0 | ⟨hr2, hrn⟩
    · subst h0
      exact StTiling.exists_goodCutCount_eq_zero n
    · exact StTiling.exists_goodCutCount_eq_of_allowed hn hr2 hrn

theorem StTiling.allUp_all_cuts_good (hn : 3 ≤ n) :
    ∀ k, k ∈ cutSet n → StTiling.IsGoodCut (StTiling.allUp n) k := by
  intro k hk
  rcases StTile.exists_crossesCut_of_mem_cutSet (n := n) hn hk with ⟨t, hcross⟩
  exact ⟨t, by simp, hcross⟩

theorem StTiling.goodCuts_allUp_eq_cutSet (hn : 3 ≤ n) :
    (StTiling.allUp n).goodCuts = cutSet n := by
  exact (StTiling.goodCuts_eq_cutSet_iff (StTiling.allUp n)).2
    (StTiling.allUp_all_cuts_good hn)

theorem StTiling.goodCutCount_allUp (hn : 3 ≤ n) :
    (StTiling.allUp n).goodCutCount = n - 1 := by
  exact (StTiling.goodCutCount_eq_top_iff (StTiling.allUp n)).2
    (StTiling.goodCuts_allUp_eq_cutSet hn)

theorem StTiling.goodCuts_complement_of_all_down {b : StTiling n}
    (hn : 3 ≤ n) (h : ∀ t : StTile n, b t = false) :
    b.complement.goodCuts = cutSet n := by
  have hb : b.complement = StTiling.allUp n := by
    funext t
    simp [StTiling.complement, StTiling.allUp, h t]
  rw [hb]
  exact StTiling.goodCuts_allUp_eq_cutSet hn

theorem StTiling.goodCutCount_complement_of_all_down {b : StTiling n}
    (hn : 3 ≤ n) (h : ∀ t : StTile n, b t = false) :
    b.complement.goodCutCount = n - 1 := by
  rw [StTiling.goodCutCount_eq_top_iff]
  exact StTiling.goodCuts_complement_of_all_down hn h

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
