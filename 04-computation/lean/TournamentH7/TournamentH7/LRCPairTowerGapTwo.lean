import TournamentH7.LRCDetunedD3

/-!
# The `(4,4,8,8)` valuation-gap-two equality wall

At the second lift of a two-minimum pair tower, the two fresh rows have
reduced denominator four and the inherited pair has reduced denominator
eight.  All four bad-row bounds are exactly one quarter of the branch circle.
This module proves the sharp local-density statement: one pair collision saves
one branch and produces a common good branch.  Consequently failure is an
exact pairwise-disjoint parallel partition of the branch circle.

The remaining crux is therefore phase chronology: choose a harmonic-good
phase at which this partition normal form does not occur.
-/

namespace LonelyRunner
namespace LRCPairTowerGapTwo

open LonelyRunner
open scoped Classical

noncomputable section

def detunedBadBranches (δ g : ℤ) (u : ℝ) : Finset ℤ :=
  (Finset.Ico (0 : ℤ) g).filter fun c =>
    ∃ n : ℤ, |(δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n| < 1 / 14

theorem detunedBadBranches_subset_Ico (δ g : ℤ) (u : ℝ) :
    detunedBadBranches δ g u ⊆ Finset.Ico (0 : ℤ) g := by
  intro c hc
  exact (Finset.mem_filter.mp hc).1

theorem detunedBadBranches_card_le (δ g : ℤ) (hg : 1 ≤ g) (u : ℝ) :
    (detunedBadBranches δ g u).card ≤ DetunedD3.badCount δ g := by
  simpa [detunedBadBranches, DetunedD3.badCount] using
    LRCIntervalCount.bad_count_le δ g hg u

theorem four_mul_badCount_eq_of_reducedDenominator_four
    (δ g : ℤ) (hg : 0 < g)
    (hq4 : g / (Int.gcd δ g : ℤ) = 4) :
    4 * DetunedD3.badCount δ g = g.toNat := by
  have hdvd : ((Int.gcd δ g : ℤ)) ∣ g := Int.gcd_dvd_right δ g
  have htoNat : (g / (Int.gcd δ g : ℤ)).toNat = 4 := by
    rw [hq4]
    rfl
  have hbad : DetunedD3.badCount δ g = Int.gcd δ g := by
    rw [DetunedD3.badCount, htoNat]
    norm_num
  have hfactor := Int.mul_ediv_cancel' hdvd
  rw [hq4] at hfactor
  rw [hbad]
  omega

theorem four_mul_badCount_eq_of_reducedDenominator_eight
    (δ g : ℤ) (hg : 0 < g)
    (hq8 : g / (Int.gcd δ g : ℤ) = 8) :
    4 * DetunedD3.badCount δ g = g.toNat := by
  have hdvd : ((Int.gcd δ g : ℤ)) ∣ g := Int.gcd_dvd_right δ g
  have htoNat : (g / (Int.gcd δ g : ℤ)).toNat = 8 := by
    rw [hq8]
    rfl
  have hbad : DetunedD3.badCount δ g = 2 * Int.gcd δ g := by
    rw [DetunedD3.badCount, htoNat]
    norm_num [Nat.mul_comm]
  have hfactor := Int.mul_ediv_cancel' hdvd
  rw [hq8] at hfactor
  rw [hbad]
  omega

def HasFourDetunedGoodBranch
    (δ₁ δ₂ δ₃ δ₄ g : ℤ) (u : ℝ) : Prop :=
  ∃ c ∈ Finset.Ico (0 : ℤ) g,
    c ∉ detunedBadBranches δ₁ g u ∧
    c ∉ detunedBadBranches δ₂ g u ∧
    c ∉ detunedBadBranches δ₃ g u ∧
    c ∉ detunedBadBranches δ₄ g u

def FourBadRowOverlap
    (δ₁ δ₂ δ₃ δ₄ g : ℤ) (u : ℝ) : Prop :=
  (detunedBadBranches δ₁ g u ∩ detunedBadBranches δ₂ g u).Nonempty ∨
  (detunedBadBranches δ₁ g u ∩ detunedBadBranches δ₃ g u).Nonempty ∨
  (detunedBadBranches δ₁ g u ∩ detunedBadBranches δ₄ g u).Nonempty ∨
  (detunedBadBranches δ₂ g u ∩ detunedBadBranches δ₃ g u).Nonempty ∨
  (detunedBadBranches δ₂ g u ∩ detunedBadBranches δ₄ g u).Nonempty ∨
  (detunedBadBranches δ₃ g u ∩ detunedBadBranches δ₄ g u).Nonempty

theorem card_four_union_lt_of_pair_overlap
    {α : Type*} [DecidableEq α]
    (first second third fourth : Finset α)
    (hoverlap :
      (first ∩ second).Nonempty ∨
      (first ∩ third).Nonempty ∨
      (first ∩ fourth).Nonempty ∨
      (second ∩ third).Nonempty ∨
      (second ∩ fourth).Nonempty ∨
      (third ∩ fourth).Nonempty) :
    (first ∪ second ∪ third ∪ fourth).card <
      first.card + second.card + third.card + fourth.card := by
  have pair_saving (left right rest₁ rest₂ : Finset α)
      (hinter : (left ∩ right).Nonempty) :
      (left ∪ right ∪ rest₁ ∪ rest₂).card <
        left.card + right.card + rest₁.card + rest₂.card := by
    have hinterPos : 0 < (left ∩ right).card := hinter.card_pos
    have hunionPair := Finset.card_union_add_card_inter left right
    have hpair : (left ∪ right).card < left.card + right.card := by omega
    calc
      (left ∪ right ∪ rest₁ ∪ rest₂).card
          ≤ (left ∪ right ∪ rest₁).card + rest₂.card :=
        Finset.card_union_le _ _
      _ ≤ ((left ∪ right).card + rest₁.card) + rest₂.card := by
        exact Nat.add_le_add_right (Finset.card_union_le _ _) _
      _ < (left.card + right.card + rest₁.card) + rest₂.card := by
        omega
  rcases hoverlap with h12 | h13 | h14 | h23 | h24 | h34
  · exact pair_saving first second third fourth h12
  · simpa [Finset.union_assoc, Finset.union_left_comm, Finset.union_comm,
      Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using
      pair_saving first third second fourth h13
  · simpa [Finset.union_assoc, Finset.union_left_comm, Finset.union_comm,
      Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using
      pair_saving first fourth second third h14
  · simpa [Finset.union_assoc, Finset.union_left_comm, Finset.union_comm,
      Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using
      pair_saving second third first fourth h23
  · simpa [Finset.union_assoc, Finset.union_left_comm, Finset.union_comm,
      Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using
      pair_saving second fourth first third h24
  · simpa [Finset.union_assoc, Finset.union_left_comm, Finset.union_comm,
      Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using
      pair_saving third fourth first second h34

theorem hasFourDetunedGoodBranch_of_pairOverlap
    (δ₁ δ₂ δ₃ δ₄ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hbudget :
      DetunedD3.badCount δ₁ g + DetunedD3.badCount δ₂ g +
        DetunedD3.badCount δ₃ g + DetunedD3.badCount δ₄ g ≤ g.toNat)
    (hoverlap : FourBadRowOverlap δ₁ δ₂ δ₃ δ₄ g u) :
    HasFourDetunedGoodBranch δ₁ δ₂ δ₃ δ₄ g u := by
  let first := detunedBadBranches δ₁ g u
  let second := detunedBadBranches δ₂ g u
  let third := detunedBadBranches δ₃ g u
  let fourth := detunedBadBranches δ₄ g u
  let branches := Finset.Ico (0 : ℤ) g
  have hfirst : first.card ≤ DetunedD3.badCount δ₁ g :=
    detunedBadBranches_card_le δ₁ g hg u
  have hsecond : second.card ≤ DetunedD3.badCount δ₂ g :=
    detunedBadBranches_card_le δ₂ g hg u
  have hthird : third.card ≤ DetunedD3.badCount δ₃ g :=
    detunedBadBranches_card_le δ₃ g hg u
  have hfourth : fourth.card ≤ DetunedD3.badCount δ₄ g :=
    detunedBadBranches_card_le δ₄ g hg u
  have hbranches : branches.card = g.toNat := by
    dsimp [branches]
    rw [Int.card_Ico]
    congr 1
    omega
  have hunion_lt :
      (first ∪ second ∪ third ∪ fourth).card < branches.card := by
    calc
      (first ∪ second ∪ third ∪ fourth).card <
          first.card + second.card + third.card + fourth.card :=
        card_four_union_lt_of_pair_overlap first second third fourth hoverlap
      _ ≤ DetunedD3.badCount δ₁ g + DetunedD3.badCount δ₂ g +
          DetunedD3.badCount δ₃ g + DetunedD3.badCount δ₄ g := by omega
      _ ≤ g.toNat := hbudget
      _ = branches.card := hbranches.symm
  have hnotsub :
      ¬ branches ⊆ first ∪ second ∪ third ∪ fourth := by
    intro hsub
    exact (Nat.not_le_of_lt hunion_lt) (Finset.card_le_card hsub)
  rw [Finset.not_subset] at hnotsub
  obtain ⟨c, hcBranches, hcUnion⟩ := hnotsub
  simp only [Finset.mem_union, not_or] at hcUnion
  exact ⟨c, hcBranches, hcUnion.1.1.1, hcUnion.1.1.2,
    hcUnion.1.2, hcUnion.2⟩

structure FourRowParallelPartition
    {α : Type*} [DecidableEq α]
    (branches first second third fourth : Finset α) : Prop where
  cover : first ∪ second ∪ third ∪ fourth = branches
  disjoint₁₂ : Disjoint first second
  disjoint₁₃ : Disjoint first third
  disjoint₁₄ : Disjoint first fourth
  disjoint₂₃ : Disjoint second third
  disjoint₂₄ : Disjoint second fourth
  disjoint₃₄ : Disjoint third fourth

theorem fourRowParallelPartition_of_failure
    (δ₁ δ₂ δ₃ δ₄ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hbudget :
      DetunedD3.badCount δ₁ g + DetunedD3.badCount δ₂ g +
        DetunedD3.badCount δ₃ g + DetunedD3.badCount δ₄ g ≤ g.toNat)
    (hfailure : ¬ HasFourDetunedGoodBranch δ₁ δ₂ δ₃ δ₄ g u) :
    FourRowParallelPartition (Finset.Ico (0 : ℤ) g)
      (detunedBadBranches δ₁ g u) (detunedBadBranches δ₂ g u)
      (detunedBadBranches δ₃ g u) (detunedBadBranches δ₄ g u) := by
  let first := detunedBadBranches δ₁ g u
  let second := detunedBadBranches δ₂ g u
  let third := detunedBadBranches δ₃ g u
  let fourth := detunedBadBranches δ₄ g u
  let branches := Finset.Ico (0 : ℤ) g
  have no_pair (hoverlap : FourBadRowOverlap δ₁ δ₂ δ₃ δ₄ g u) : False :=
    hfailure (hasFourDetunedGoodBranch_of_pairOverlap
      δ₁ δ₂ δ₃ δ₄ g u hg hbudget hoverlap)
  have hcover : first ∪ second ∪ third ∪ fourth = branches := by
    apply Finset.Subset.antisymm
    · intro branch hbranch
      simp only [Finset.mem_union] at hbranch
      rcases hbranch with ((hfirst | hsecond) | hthird) | hfourth
      · exact detunedBadBranches_subset_Ico δ₁ g u hfirst
      · exact detunedBadBranches_subset_Ico δ₂ g u hsecond
      · exact detunedBadBranches_subset_Ico δ₃ g u hthird
      · exact detunedBadBranches_subset_Ico δ₄ g u hfourth
    · intro branch hbranch
      by_contra houtside
      apply hfailure
      simp only [Finset.mem_union, not_or] at houtside
      exact ⟨branch, hbranch, houtside.1.1.1, houtside.1.1.2,
        houtside.1.2, houtside.2⟩
  have h12 : Disjoint first second := by
    rw [Finset.disjoint_left]
    intro branch hfirst hsecond
    exact no_pair (Or.inl ⟨branch, Finset.mem_inter.mpr ⟨hfirst, hsecond⟩⟩)
  have h13 : Disjoint first third := by
    rw [Finset.disjoint_left]
    intro branch hfirst hthird
    exact no_pair (Or.inr (Or.inl
      ⟨branch, Finset.mem_inter.mpr ⟨hfirst, hthird⟩⟩))
  have h14 : Disjoint first fourth := by
    rw [Finset.disjoint_left]
    intro branch hfirst hfourth
    exact no_pair (Or.inr (Or.inr (Or.inl
      ⟨branch, Finset.mem_inter.mpr ⟨hfirst, hfourth⟩⟩)))
  have h23 : Disjoint second third := by
    rw [Finset.disjoint_left]
    intro branch hsecond hthird
    exact no_pair (Or.inr (Or.inr (Or.inr (Or.inl
      ⟨branch, Finset.mem_inter.mpr ⟨hsecond, hthird⟩⟩))))
  have h24 : Disjoint second fourth := by
    rw [Finset.disjoint_left]
    intro branch hsecond hfourth
    exact no_pair (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl
      ⟨branch, Finset.mem_inter.mpr ⟨hsecond, hfourth⟩⟩)))))
  have h34 : Disjoint third fourth := by
    rw [Finset.disjoint_left]
    intro branch hthird hfourth
    exact no_pair (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr
      ⟨branch, Finset.mem_inter.mpr ⟨hthird, hfourth⟩⟩)))))
  exact ⟨hcover, h12, h13, h14, h23, h24, h34⟩

/-- The exact-cover half of a parallel partition already forces failure;
neither a cardinality budget nor pairwise disjointness is needed in this
direction. -/
theorem failure_of_fourRowParallelPartition
    (δ₁ δ₂ δ₃ δ₄ g : ℤ) (u : ℝ)
    (hpartition :
      FourRowParallelPartition (Finset.Ico (0 : ℤ) g)
        (detunedBadBranches δ₁ g u) (detunedBadBranches δ₂ g u)
        (detunedBadBranches δ₃ g u) (detunedBadBranches δ₄ g u)) :
    ¬ HasFourDetunedGoodBranch δ₁ δ₂ δ₃ δ₄ g u := by
  rintro ⟨branch, hbranch, hfirst, hsecond, hthird, hfourth⟩
  have hbad :
      branch ∈ detunedBadBranches δ₁ g u ∪ detunedBadBranches δ₂ g u ∪
        detunedBadBranches δ₃ g u ∪ detunedBadBranches δ₄ g u := by
    rw [hpartition.cover]
    exact hbranch
  simp only [Finset.mem_union] at hbad
  rcases hbad with ((hbad | hbad) | hbad) | hbad
  · exact hfirst hbad
  · exact hsecond hbad
  · exact hthird hbad
  · exact hfourth hbad

/-- **Static gap-two failure wall.**  Under the sharp four-row cardinality
budget, failure of local-density branch gluing is equivalent to the four
actual bad rows forming their exact pairwise-disjoint parallel partition of
the branch circle. -/
theorem failure_iff_fourRowParallelPartition
    (δ₁ δ₂ δ₃ δ₄ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hbudget :
      DetunedD3.badCount δ₁ g + DetunedD3.badCount δ₂ g +
        DetunedD3.badCount δ₃ g + DetunedD3.badCount δ₄ g ≤ g.toNat) :
    (¬ HasFourDetunedGoodBranch δ₁ δ₂ δ₃ δ₄ g u) ↔
      FourRowParallelPartition (Finset.Ico (0 : ℤ) g)
        (detunedBadBranches δ₁ g u) (detunedBadBranches δ₂ g u)
        (detunedBadBranches δ₃ g u) (detunedBadBranches δ₄ g u) := by
  constructor
  · exact fourRowParallelPartition_of_failure δ₁ δ₂ δ₃ δ₄ g u hg hbudget
  · exact failure_of_fourRowParallelPartition δ₁ δ₂ δ₃ δ₄ g u

theorem four_four_eight_eight_failure_parallel_partition
    (δ₄a δ₄b δ₈a δ₈b g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4)
    (hq8a : g / (Int.gcd δ₈a g : ℤ) = 8)
    (hq8b : g / (Int.gcd δ₈b g : ℤ) = 8)
    (hfailure : ¬ HasFourDetunedGoodBranch δ₄a δ₄b δ₈a δ₈b g u) :
    FourRowParallelPartition (Finset.Ico (0 : ℤ) g)
        (detunedBadBranches δ₄a g u) (detunedBadBranches δ₄b g u)
        (detunedBadBranches δ₈a g u) (detunedBadBranches δ₈b g u) ∧
      4 * (detunedBadBranches δ₄a g u).card = g.toNat ∧
      4 * (detunedBadBranches δ₄b g u).card = g.toNat ∧
      4 * (detunedBadBranches δ₈a g u).card = g.toNat ∧
      4 * (detunedBadBranches δ₈b g u).card = g.toNat := by
  have hg0 : 0 < g := by omega
  have hbad4a := four_mul_badCount_eq_of_reducedDenominator_four
    δ₄a g hg0 hq4a
  have hbad4b := four_mul_badCount_eq_of_reducedDenominator_four
    δ₄b g hg0 hq4b
  have hbad8a := four_mul_badCount_eq_of_reducedDenominator_eight
    δ₈a g hg0 hq8a
  have hbad8b := four_mul_badCount_eq_of_reducedDenominator_eight
    δ₈b g hg0 hq8b
  have hbudget :
      DetunedD3.badCount δ₄a g + DetunedD3.badCount δ₄b g +
        DetunedD3.badCount δ₈a g + DetunedD3.badCount δ₈b g ≤ g.toNat := by
    omega
  have hpartition := fourRowParallelPartition_of_failure
    δ₄a δ₄b δ₈a δ₈b g u (by omega) hbudget hfailure
  have hcard4a := detunedBadBranches_card_le δ₄a g (by omega) u
  have hcard4b := detunedBadBranches_card_le δ₄b g (by omega) u
  have hcard8a := detunedBadBranches_card_le δ₈a g (by omega) u
  have hcard8b := detunedBadBranches_card_le δ₈b g (by omega) u
  have hbranches : (Finset.Ico (0 : ℤ) g).card = g.toNat := by
    rw [Int.card_Ico]
    congr 1
    omega
  have hunionCard :
      (detunedBadBranches δ₄a g u ∪ detunedBadBranches δ₄b g u ∪
        detunedBadBranches δ₈a g u ∪ detunedBadBranches δ₈b g u).card =
        g.toNat := by
    rw [hpartition.cover, hbranches]
  have hunionLe :
      (detunedBadBranches δ₄a g u ∪ detunedBadBranches δ₄b g u ∪
        detunedBadBranches δ₈a g u ∪ detunedBadBranches δ₈b g u).card ≤
      (detunedBadBranches δ₄a g u).card +
        (detunedBadBranches δ₄b g u).card +
        (detunedBadBranches δ₈a g u).card +
        (detunedBadBranches δ₈b g u).card := by
    calc
      _ ≤ (detunedBadBranches δ₄a g u ∪
          detunedBadBranches δ₄b g u ∪ detunedBadBranches δ₈a g u).card +
          (detunedBadBranches δ₈b g u).card := Finset.card_union_le _ _
      _ ≤ ((detunedBadBranches δ₄a g u ∪
          detunedBadBranches δ₄b g u).card +
          (detunedBadBranches δ₈a g u).card) +
          (detunedBadBranches δ₈b g u).card := by
        exact Nat.add_le_add_right (Finset.card_union_le _ _) _
      _ ≤ (((detunedBadBranches δ₄a g u).card +
          (detunedBadBranches δ₄b g u).card) +
          (detunedBadBranches δ₈a g u).card) +
          (detunedBadBranches δ₈b g u).card := by
        exact Nat.add_le_add_right
          (Nat.add_le_add_right (Finset.card_union_le _ _) _) _
  refine ⟨hpartition, ?_, ?_, ?_, ?_⟩ <;> omega

def FourDetunedHarmonicGoodAt
    (v : Fin 13 → ℤ) (g : ℤ)
    (i₁ i₂ i₃ i₄ : Fin 13) (u : ℝ) : Prop :=
  ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ → ∀ n : ℤ,
    (1 : ℝ) / 14 ≤ |((v j / g : ℤ) : ℝ) * u - n|

theorem HasFourDetunedGoodBranch.clearances
    {δ₁ δ₂ δ₃ δ₄ g : ℤ} {u : ℝ}
    (hgood : HasFourDetunedGoodBranch δ₁ δ₂ δ₃ δ₄ g u) :
    ∃ c : ℤ,
      (∀ n : ℤ, (1 : ℝ) / 14 ≤
        |(δ₁ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|) ∧
      (∀ n : ℤ, (1 : ℝ) / 14 ≤
        |(δ₂ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|) ∧
      (∀ n : ℤ, (1 : ℝ) / 14 ≤
        |(δ₃ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|) ∧
      (∀ n : ℤ, (1 : ℝ) / 14 ≤
        |(δ₄ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|) := by
  obtain ⟨c, hcIco, hc1, hc2, hc3, hc4⟩ := hgood
  refine ⟨c, ?_, ?_, ?_, ?_⟩ <;> intro n
  · exact not_lt.mp fun hlt => hc1 (by
      rw [detunedBadBranches, Finset.mem_filter]
      exact ⟨hcIco, n, hlt⟩)
  · exact not_lt.mp fun hlt => hc2 (by
      rw [detunedBadBranches, Finset.mem_filter]
      exact ⟨hcIco, n, hlt⟩)
  · exact not_lt.mp fun hlt => hc3 (by
      rw [detunedBadBranches, Finset.mem_filter]
      exact ⟨hcIco, n, hlt⟩)
  · exact not_lt.mp fun hlt => hc4 (by
      rw [detunedBadBranches, Finset.mem_filter]
      exact ⟨hcIco, n, hlt⟩)

theorem lonely14_of_four_detuned_selectedWitness
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (i₁ i₂ i₃ i₄ : Fin 13)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ → g ∣ v j)
    (hselect : ∃ u : ℝ,
      FourDetunedHarmonicGoodAt v g i₁ i₂ i₃ i₄ u ∧
      HasFourDetunedGoodBranch (v i₁) (v i₂) (v i₃) (v i₄) g u) :
    ∃ t : ℝ, Lonely 14 v t := by
  obtain ⟨u, hharm, hgood⟩ := hselect
  obtain ⟨c, hc1, hc2, hc3, hc4⟩ := hgood.clearances
  refine ⟨(u + (c : ℝ)) / (g : ℝ), fun i n => ?_⟩
  rcases eq_or_ne i i₁ with rfl | h1
  · exact hc1 n
  rcases eq_or_ne i i₂ with rfl | h2
  · exact hc2 n
  rcases eq_or_ne i i₃ with rfl | h3
  · exact hc3 n
  rcases eq_or_ne i i₄ with rfl | h4
  · exact hc4 n
  have hdiv := hdvd i h1 h2 h3 h4
  have hvi : (v i : ℝ) = (g : ℝ) * ((v i / g : ℤ) : ℝ) := by
    have : v i = g * (v i / g) := (Int.mul_ediv_cancel' hdiv).symm
    exact_mod_cast this
  have hval :
      (v i : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n =
        ((v i / g : ℤ) : ℝ) * u -
          (((n - (v i / g) * c : ℤ)) : ℝ) := by
    rw [hvi]
    push_cast
    field_simp
    ring
  rw [hval]
  exact hharm i h1 h2 h3 h4 (n - (v i / g) * c)

#print axioms hasFourDetunedGoodBranch_of_pairOverlap
#print axioms fourRowParallelPartition_of_failure
#print axioms failure_of_fourRowParallelPartition
#print axioms failure_iff_fourRowParallelPartition
#print axioms four_four_eight_eight_failure_parallel_partition
#print axioms lonely14_of_four_detuned_selectedWitness

end
end LRCPairTowerGapTwo
end LonelyRunner
