/-
  TournamentH7.LRCDetunedOverlap

  Pair-intersection sharpening for the saturated three-detuned branch count.
  A shared bad branch saves one unit in the three-set union bound.  Thus the
  uniform q=3 residual can fail only when its three bad rows partition the
  parallel-class circle exactly and each row attains its full g/3 count.

  Tournament-analysis audit: the vertices here are detuned bad rows, not
  runners or arcs.  The pairwise observable is nonempty row intersection; at
  fixed (g,u), thresholding its cardinality at one gives the binary collision
  relation.  This relation is symmetric rather than a tournament orientation,
  so an oriented tie Hamiltonian path would require an arbitrary index gauge
  and would discard which shared branch supplies the saving.  We retain the
  exact intersection witness instead.  The quotient preserves the good-branch
  predicate and destroys cyclic placement inside each bad row.  The challenged
  assumption is that saturated row degrees alone can obstruct clearing: they
  can do so only in the exact pairwise-disjoint partition normal form below.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCEndgameParameterDischargeTwoThree

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

noncomputable section

/-- Branch classes at which one detuned coordinate fails the `1/14`
clearance test. -/
def detunedBadBranches (δ g : ℤ) (u : ℝ) : Finset ℤ :=
  (Finset.Ico (0 : ℤ) g).filter fun c =>
    ∃ n : ℤ, |(δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n| < 1 / 14

/-- A branch in `[0,g)` which clears all three detuned coordinates. -/
def HasThreeDetunedGoodBranch (δ₁ δ₂ δ₃ g : ℤ) (u : ℝ) : Prop :=
  ∃ c ∈ Finset.Ico (0 : ℤ) g,
    c ∉ detunedBadBranches δ₁ g u ∧
    c ∉ detunedBadBranches δ₂ g u ∧
    c ∉ detunedBadBranches δ₃ g u

theorem detunedBadBranches_subset_Ico (δ g : ℤ) (u : ℝ) :
    detunedBadBranches δ g u ⊆ Finset.Ico (0 : ℤ) g := by
  intro c hc
  exact (Finset.mem_filter.mp hc).1

theorem detunedBadBranches_card_le (δ g : ℤ) (hg : 1 ≤ g) (u : ℝ) :
    (detunedBadBranches δ g u).card ≤ DetunedD3.badCount δ g := by
  simpa [detunedBadBranches, DetunedD3.badCount] using
    LRCIntervalCount.bad_count_le δ g hg u

/-- Any pairwise collision saves at least one unit from the naive three-row
degree sum. -/
theorem card_three_union_lt_of_pair_overlap {α : Type*} [DecidableEq α]
    (first second third : Finset α)
    (hoverlap :
      (first ∩ second).Nonempty ∨
      (first ∩ third).Nonempty ∨
      (second ∩ third).Nonempty) :
    (first ∪ second ∪ third).card <
      first.card + second.card + third.card := by
  have pair_saving (left right rest : Finset α)
      (hinter : (left ∩ right).Nonempty) :
      (left ∪ right ∪ rest).card < left.card + right.card + rest.card := by
    have hinterPos : 0 < (left ∩ right).card := hinter.card_pos
    have hunionPair := Finset.card_union_add_card_inter left right
    have hpair : (left ∪ right).card < left.card + right.card := by omega
    exact lt_of_le_of_lt (Finset.card_union_le (left ∪ right) rest)
      (Nat.add_lt_add_right hpair rest.card)
  rcases hoverlap with h12 | h13 | h23
  · exact pair_saving first second third h12
  · simpa [Finset.union_assoc, Finset.union_left_comm, Finset.union_comm,
      Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using
      pair_saving first third second h13
  · simpa [Finset.union_assoc, Finset.union_left_comm, Finset.union_comm,
      Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using
      pair_saving second third first h23

/-- A pair intersection upgrades the weak saturated degree budget `≤ g` to a
strict union bound, hence produces a common good branch. -/
theorem hasThreeDetunedGoodBranch_of_pairOverlap
    (δ₁ δ₂ δ₃ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hbudget : DetunedD3.badCount δ₁ g + DetunedD3.badCount δ₂ g +
      DetunedD3.badCount δ₃ g ≤ g.toNat)
    (hoverlap :
      (detunedBadBranches δ₁ g u ∩ detunedBadBranches δ₂ g u).Nonempty ∨
      (detunedBadBranches δ₁ g u ∩ detunedBadBranches δ₃ g u).Nonempty ∨
      (detunedBadBranches δ₂ g u ∩ detunedBadBranches δ₃ g u).Nonempty) :
    HasThreeDetunedGoodBranch δ₁ δ₂ δ₃ g u := by
  let first := detunedBadBranches δ₁ g u
  let second := detunedBadBranches δ₂ g u
  let third := detunedBadBranches δ₃ g u
  let branches := Finset.Ico (0 : ℤ) g
  have hfirst : first.card ≤ DetunedD3.badCount δ₁ g :=
    detunedBadBranches_card_le δ₁ g hg u
  have hsecond : second.card ≤ DetunedD3.badCount δ₂ g :=
    detunedBadBranches_card_le δ₂ g hg u
  have hthird : third.card ≤ DetunedD3.badCount δ₃ g :=
    detunedBadBranches_card_le δ₃ g hg u
  have hbranches : branches.card = g.toNat := by
    dsimp [branches]
    rw [Int.card_Ico]
    congr 1
    omega
  have hunion_lt : (first ∪ second ∪ third).card < branches.card := by
    calc
      (first ∪ second ∪ third).card < first.card + second.card + third.card :=
        card_three_union_lt_of_pair_overlap first second third hoverlap
      _ ≤ DetunedD3.badCount δ₁ g + DetunedD3.badCount δ₂ g +
          DetunedD3.badCount δ₃ g := by omega
      _ ≤ g.toNat := hbudget
      _ = branches.card := hbranches.symm
  have hnotsub : ¬ branches ⊆ first ∪ second ∪ third := by
    intro hsub
    exact (Nat.not_le_of_lt hunion_lt) (Finset.card_le_card hsub)
  rw [Finset.not_subset] at hnotsub
  obtain ⟨c, hcBranches, hcUnion⟩ := hnotsub
  simp only [Finset.mem_union, not_or] at hcUnion
  exact ⟨c, hcBranches, hcUnion.1.1, hcUnion.1.2, hcUnion.2⟩

/-- Nonmembership in the bad-row representation is exactly the clearance
statement consumed by the three-detuned construction. -/
theorem HasThreeDetunedGoodBranch.clearances
    {δ₁ δ₂ δ₃ g : ℤ} {u : ℝ}
    (hgood : HasThreeDetunedGoodBranch δ₁ δ₂ δ₃ g u) :
    ∃ c : ℤ,
      (∀ n : ℤ, (1 : ℝ) / 14 ≤
        |(δ₁ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|) ∧
      (∀ n : ℤ, (1 : ℝ) / 14 ≤
        |(δ₂ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|) ∧
      (∀ n : ℤ, (1 : ℝ) / 14 ≤
        |(δ₃ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|) := by
  obtain ⟨c, hcIco, hc1, hc2, hc3⟩ := hgood
  refine ⟨c, ?_, ?_, ?_⟩
  · intro n
    exact not_lt.mp fun hlt => hc1 (by
      rw [detunedBadBranches, Finset.mem_filter]
      exact ⟨hcIco, n, hlt⟩)
  · intro n
    exact not_lt.mp fun hlt => hc2 (by
      rw [detunedBadBranches, Finset.mem_filter]
      exact ⟨hcIco, n, hlt⟩)
  · intro n
    exact not_lt.mp fun hlt => hc3 (by
      rw [detunedBadBranches, Finset.mem_filter]
      exact ⟨hcIco, n, hlt⟩)

/-- A phase-by-phase pair collision supplies the instance-clearing interface
consumed by the harmonic LRC reduction. -/
theorem threeDetunedInstanceClearing_of_pairOverlap
    (δ₁ δ₂ δ₃ g : ℤ) (hg : 1 ≤ g)
    (hbudget : DetunedD3.badCount δ₁ g + DetunedD3.badCount δ₂ g +
      DetunedD3.badCount δ₃ g ≤ g.toNat)
    (hoverlap : ∀ u : ℝ,
      (detunedBadBranches δ₁ g u ∩ detunedBadBranches δ₂ g u).Nonempty ∨
      (detunedBadBranches δ₁ g u ∩ detunedBadBranches δ₃ g u).Nonempty ∨
      (detunedBadBranches δ₂ g u ∩ detunedBadBranches δ₃ g u).Nonempty) :
    DetunedD3.ThreeDetunedInstanceClearing δ₁ δ₂ δ₃ g := by
  intro u
  exact (hasThreeDetunedGoodBranch_of_pairOverlap
    δ₁ δ₂ δ₃ g u hg hbudget (hoverlap u)).clearances

/-- Direct LRC consumer for the Zarankiewicz overlap saving.  The LRC(10)
citation handles the harmonic coordinates, while a pairwise bad-row collision
at every harmonic phase supplies the common branch. -/
theorem lonely14_of_three_detuned_pairOverlap (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₁ i₂ i₃ : Fin 13) (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (hbudget : DetunedD3.badCount (v i₁) g + DetunedD3.badCount (v i₂) g +
      DetunedD3.badCount (v i₃) g ≤ g.toNat)
    (hoverlap : ∀ u : ℝ,
      (detunedBadBranches (v i₁) g u ∩
        detunedBadBranches (v i₂) g u).Nonempty ∨
      (detunedBadBranches (v i₁) g u ∩
        detunedBadBranches (v i₃) g u).Nonempty ∨
      (detunedBadBranches (v i₂) g u ∩
        detunedBadBranches (v i₃) g u).Nonempty) :
    ∃ t : ℝ, Lonely 14 v t :=
  DetunedD3.lonely14_of_three_detuned_instance cite v hv g hg
    i₁ i₂ i₃ h12 h13 h23 hdvd
    (threeDetunedInstanceClearing_of_pairOverlap
      (v i₁) (v i₂) (v i₃) g (by omega) hbudget hoverlap)

/-- Exact obstruction certificate for the uniform `(3,3,3)` row-degree
pattern: the three bad rows partition all branch classes, attain size `g/3`,
and are pairwise disjoint. -/
structure UniformThreeBadPartition (δ₁ δ₂ δ₃ g : ℤ) (u : ℝ) : Prop where
  cover :
    detunedBadBranches δ₁ g u ∪ detunedBadBranches δ₂ g u ∪
      detunedBadBranches δ₃ g u = Finset.Ico (0 : ℤ) g
  first_card : 3 * (detunedBadBranches δ₁ g u).card = g.toNat
  second_card : 3 * (detunedBadBranches δ₂ g u).card = g.toNat
  third_card : 3 * (detunedBadBranches δ₃ g u).card = g.toNat
  first_second_disjoint :
    Disjoint (detunedBadBranches δ₁ g u) (detunedBadBranches δ₂ g u)
  first_third_disjoint :
    Disjoint (detunedBadBranches δ₁ g u) (detunedBadBranches δ₃ g u)
  second_third_disjoint :
    Disjoint (detunedBadBranches δ₂ g u) (detunedBadBranches δ₃ g u)

/-- If no common clearing branch exists in the uniform `q=3` pattern, every
inequality in the degree/union argument is forced to equality. -/
theorem uniformThreeBadPartition_of_noGoodBranch
    (δ₁ δ₂ δ₃ g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq1 : g / (Int.gcd δ₁ g : ℤ) = 3)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 3)
    (hq3 : g / (Int.gcd δ₃ g : ℤ) = 3)
    (hno : ¬ HasThreeDetunedGoodBranch δ₁ δ₂ δ₃ g u) :
    UniformThreeBadPartition δ₁ δ₂ δ₃ g u := by
  let first := detunedBadBranches δ₁ g u
  let second := detunedBadBranches δ₂ g u
  let third := detunedBadBranches δ₃ g u
  let branches := Finset.Ico (0 : ℤ) g
  have hg0 : 0 < g := by omega
  have hbad1 := three_mul_badCount_eq δ₁ g hg0 hq1
  have hbad2 := three_mul_badCount_eq δ₂ g hg0 hq2
  have hbad3 := three_mul_badCount_eq δ₃ g hg0 hq3
  have hbadTotal : DetunedD3.badCount δ₁ g + DetunedD3.badCount δ₂ g +
      DetunedD3.badCount δ₃ g = g.toNat := by omega
  have hfirst : first.card ≤ DetunedD3.badCount δ₁ g :=
    detunedBadBranches_card_le δ₁ g hg u
  have hsecond : second.card ≤ DetunedD3.badCount δ₂ g :=
    detunedBadBranches_card_le δ₂ g hg u
  have hthird : third.card ≤ DetunedD3.badCount δ₃ g :=
    detunedBadBranches_card_le δ₃ g hg u
  have hbranches : branches.card = g.toNat := by
    dsimp [branches]
    rw [Int.card_Ico]
    congr 1
    omega
  have hcover : first ∪ second ∪ third = branches := by
    apply Finset.Subset.antisymm
    · exact Finset.union_subset
        (Finset.union_subset
          (detunedBadBranches_subset_Ico δ₁ g u)
          (detunedBadBranches_subset_Ico δ₂ g u))
        (detunedBadBranches_subset_Ico δ₃ g u)
    · intro c hc
      by_contra hcUnion
      apply hno
      simp only [Finset.mem_union, not_or] at hcUnion
      exact ⟨c, hc, hcUnion.1.1, hcUnion.1.2, hcUnion.2⟩
  have hsumLower : g.toNat ≤ first.card + second.card + third.card := by
    calc
      g.toNat = branches.card := hbranches.symm
      _ = (first ∪ second ∪ third).card := congrArg Finset.card hcover.symm
      _ ≤ (first ∪ second).card + third.card := Finset.card_union_le _ _
      _ ≤ first.card + second.card + third.card :=
        Nat.add_le_add_right (Finset.card_union_le first second) third.card
  have hsumCards : first.card + second.card + third.card = g.toNat := by
    omega
  have hfirstCard : 3 * first.card = g.toNat := by omega
  have hsecondCard : 3 * second.card = g.toNat := by omega
  have hthirdCard : 3 * third.card = g.toNat := by omega
  have hbudget : DetunedD3.badCount δ₁ g + DetunedD3.badCount δ₂ g +
      DetunedD3.badCount δ₃ g ≤ g.toNat := hbadTotal.le
  have hdisjoint12 : Disjoint first second := by
    rw [Finset.disjoint_left]
    intro c hcFirst hcSecond
    apply hno
    apply hasThreeDetunedGoodBranch_of_pairOverlap δ₁ δ₂ δ₃ g u hg hbudget
    exact Or.inl ⟨c, by simpa [first, second] using And.intro hcFirst hcSecond⟩
  have hdisjoint13 : Disjoint first third := by
    rw [Finset.disjoint_left]
    intro c hcFirst hcThird
    apply hno
    apply hasThreeDetunedGoodBranch_of_pairOverlap δ₁ δ₂ δ₃ g u hg hbudget
    exact Or.inr (Or.inl
      ⟨c, by simpa [first, third] using And.intro hcFirst hcThird⟩)
  have hdisjoint23 : Disjoint second third := by
    rw [Finset.disjoint_left]
    intro c hcSecond hcThird
    apply hno
    apply hasThreeDetunedGoodBranch_of_pairOverlap δ₁ δ₂ δ₃ g u hg hbudget
    exact Or.inr (Or.inr
      ⟨c, by simpa [second, third] using And.intro hcSecond hcThird⟩)
  exact ⟨hcover, hfirstCard, hsecondCard, hthirdCard,
    hdisjoint12, hdisjoint13, hdisjoint23⟩

/-! ## Axiom audit -/

#print axioms card_three_union_lt_of_pair_overlap
#print axioms hasThreeDetunedGoodBranch_of_pairOverlap
#print axioms HasThreeDetunedGoodBranch.clearances
#print axioms threeDetunedInstanceClearing_of_pairOverlap
#print axioms lonely14_of_three_detuned_pairOverlap
#print axioms uniformThreeBadPartition_of_noGoodBranch

end
end LRC14Grand
end LonelyRunner
