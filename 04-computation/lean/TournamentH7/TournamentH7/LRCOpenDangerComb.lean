import TournamentH7.LRCCommensuration

/-!
Finite strict-open danger combs on the fundamental interval `(0,1)` and the
line-order collapse of their bipartite overlap graph.  These are shared by the
pair-grid ledger and the sorted-tooth cluster-gap route.
-/

namespace LonelyRunner
namespace LRCOpenDangerComb

open Finset MeasureTheory

noncomputable section

/-- Edges in the bipartite overlap graph of two finite real-interval classes. -/
noncomputable def intervalOverlapEdges
    {α β : Type*} [Fintype α] [Fintype β]
    (leftA rightA : α → ℝ) (leftB rightB : β → ℝ) : Finset (α × β) := by
  classical
  exact (Finset.univ.product Finset.univ).filter fun edge =>
    max (leftA edge.1) (leftB edge.2) < min (rightA edge.1) (rightB edge.2)

/-- **Parallel-class Zarankiewicz collapse on the line.**  If the intervals
inside each of two classes are pairwise disjoint, their bipartite overlap
graph has at most `|α|+|β|` edges.  The proof maps an overlap to the class
vertex with the later left endpoint.  Two edges cannot have the same owner,
because their opposite-class intervals would then meet just to the right of
that endpoint.  Thus the generic `K₂,₂` Zarankiewicz scale collapses to an
additive bound for ordered interval classes. -/
theorem card_intervalOverlapEdges_le_add
    {α β : Type*} [Fintype α] [Fintype β]
    (leftA rightA : α → ℝ) (leftB rightB : β → ℝ)
    (hA : ∀ a a', a ≠ a' →
      Disjoint (Set.Ioo (leftA a) (rightA a)) (Set.Ioo (leftA a') (rightA a')))
    (hB : ∀ b b', b ≠ b' →
      Disjoint (Set.Ioo (leftB b) (rightB b)) (Set.Ioo (leftB b') (rightB b'))) :
    (intervalOverlapEdges leftA rightA leftB rightB).card ≤
      Fintype.card α + Fintype.card β := by
  classical
  let owner : α × β → α ⊕ β := fun edge =>
    if leftB edge.2 ≤ leftA edge.1 then Sum.inl edge.1 else Sum.inr edge.2
  have hinjective : Set.InjOn owner
      ↑(intervalOverlapEdges leftA rightA leftB rightB) := by
    rintro ⟨a, b⟩ hedge ⟨a', b'⟩ hedge' howner
    change (a, b) ∈ intervalOverlapEdges leftA rightA leftB rightB at hedge
    change (a', b') ∈ intervalOverlapEdges leftA rightA leftB rightB at hedge'
    simp [intervalOverlapEdges] at hedge hedge'
    by_cases hab : leftB b ≤ leftA a
    · by_cases hab' : leftB b' ≤ leftA a'
      · simp only [owner, hab, hab', if_true, Sum.inl.injEq] at howner
        subst a'
        congr 1
        by_contra hbb'
        have haRightB : leftA a < rightB b := hedge.2.1
        have haRightB' : leftA a < rightB b' := hedge'.2.1
        let point : ℝ :=
          (leftA a + min (rightB b) (rightB b')) / 2
        have hpointB : point ∈ Set.Ioo (leftB b) (rightB b) := by
          constructor
          · dsimp [point]
            have : leftA a < min (rightB b) (rightB b') := lt_min haRightB haRightB'
            linarith
          · dsimp [point]
            have hmin : min (rightB b) (rightB b') ≤ rightB b := min_le_left _ _
            have : leftA a < min (rightB b) (rightB b') := lt_min haRightB haRightB'
            linarith
        have hpointB' : point ∈ Set.Ioo (leftB b') (rightB b') := by
          constructor
          · dsimp [point]
            have : leftA a < min (rightB b) (rightB b') := lt_min haRightB haRightB'
            linarith
          · dsimp [point]
            have hmin : min (rightB b) (rightB b') ≤ rightB b' := min_le_right _ _
            have : leftA a < min (rightB b) (rightB b') := lt_min haRightB haRightB'
            linarith
        have hdisjoint := hB b b' hbb'
        rw [Set.disjoint_left] at hdisjoint
        exact hdisjoint hpointB hpointB'
      · simp only [owner, hab, hab', if_true, if_false] at howner
        cases howner
    · by_cases hab' : leftB b' ≤ leftA a'
      · simp only [owner, hab, hab', if_true, if_false] at howner
        cases howner
      · simp only [owner, hab, hab', if_false, Sum.inr.injEq] at howner
        subst b'
        congr 1
        by_contra haa'
        have hbRightA : leftB b < rightA a := hedge.1.2
        have hbRightA' : leftB b < rightA a' := hedge'.1.2
        have hleftA : leftA a < leftB b := lt_of_not_ge hab
        have hleftA' : leftA a' < leftB b := lt_of_not_ge hab'
        let point : ℝ :=
          (leftB b + min (rightA a) (rightA a')) / 2
        have hpointA : point ∈ Set.Ioo (leftA a) (rightA a) := by
          constructor
          · dsimp [point]
            have : leftB b < min (rightA a) (rightA a') := lt_min hbRightA hbRightA'
            linarith
          · dsimp [point]
            have hmin : min (rightA a) (rightA a') ≤ rightA a := min_le_left _ _
            have : leftB b < min (rightA a) (rightA a') := lt_min hbRightA hbRightA'
            linarith
        have hpointA' : point ∈ Set.Ioo (leftA a') (rightA a') := by
          constructor
          · dsimp [point]
            have : leftB b < min (rightA a) (rightA a') := lt_min hbRightA hbRightA'
            linarith
          · dsimp [point]
            have hmin : min (rightA a) (rightA a') ≤ rightA a' := min_le_right _ _
            have : leftB b < min (rightA a) (rightA a') := lt_min hbRightA hbRightA'
            linarith
        have hdisjoint := hA a a' haa'
        rw [Set.disjoint_left] at hdisjoint
        exact hdisjoint hpointA hpointA'
  calc
    (intervalOverlapEdges leftA rightA leftB rightB).card =
        ((intervalOverlapEdges leftA rightA leftB rightB).image owner).card := by
      symm
      exact Finset.card_image_iff.mpr hinjective
    _ ≤ (Finset.univ : Finset (α ⊕ β)).card := by
      apply Finset.card_le_card
      exact Finset.subset_univ _
    _ = Fintype.card α + Fintype.card β := by simp

/-- Left endpoint of the `k`-th positive-speed danger tooth after clipping to
the open fundamental interval `(0,1)`. -/
def openCombLeft (w : ℕ) (k : Fin (w + 1)) : ℝ :=
  max 0 ((((k : ℕ) : ℝ) - (1 : ℝ) / 14) / (w : ℝ))

/-- Right endpoint of the `k`-th positive-speed danger tooth after clipping to
the open fundamental interval `(0,1)`. -/
def openCombRight (w : ℕ) (k : Fin (w + 1)) : ℝ :=
  min 1 ((((k : ℕ) : ℝ) + (1 : ℝ) / 14) / (w : ℝ))

/-- Finite strict comb carried by one positive-speed danger set on `(0,1)`.
There are `w+1` clipped teeth: the two end teeth are the two linear pieces of
the single circular component through zero. -/
def openCombRegion (w : ℕ) : Set ℝ :=
  ⋃ k : Fin (w + 1), Set.Ioo (openCombLeft w k) (openCombRight w k)

/-- Reversing a speed only reverses its circle position, so its centered
danger comb depends on the absolute speed. -/
theorem danger_natAbs (speed : ℤ) (radius : ℝ) :
    LRCCommensuration.danger (speed.natAbs : ℤ) 0 radius =
      LRCCommensuration.danger speed 0 radius := by
  ext x
  simp only [LRCCommensuration.danger, LRCCommensuration.runnerMap,
    Set.mem_preimage, Metric.mem_ball, add_zero, dist_zero_right]
  by_cases hspeed : 0 ≤ speed
  · have hcast : (speed.natAbs : ℤ) = speed := by
      rw [Int.natCast_natAbs, abs_of_nonneg hspeed]
    rw [hcast]
  · have hspeedNeg : speed < 0 := lt_of_not_ge hspeed
    have hcast : (speed.natAbs : ℤ) = -speed := by
      rw [Int.natCast_natAbs, abs_of_neg hspeedNeg]
    rw [hcast, neg_smul, norm_neg]

/-- **Concrete normalized comb carrier.**  On the open fundamental interval,
the finite clipped comb is exactly the pullback of the positive-speed danger
set on the unit circle. -/
theorem mem_openCombRegion_iff_mem_danger
    (w : ℕ) (hw : 0 < w) (x : ℝ) (hx : x ∈ Set.Ioo (0 : ℝ) 1) :
    x ∈ openCombRegion w ↔
      (x : UnitAddCircle) ∈
        LRCCommensuration.danger (w : ℤ) (0 : UnitAddCircle) (1 / 14) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  have hpoint :
      (w : ℤ) • (x : UnitAddCircle) = (((w : ℝ) * x : ℝ) : UnitAddCircle) := by
    rw [← AddCircle.coe_zsmul]
    congr 1
    simp [zsmul_eq_mul]
  have hdanger :
      (x : UnitAddCircle) ∈
          LRCCommensuration.danger (w : ℤ) (0 : UnitAddCircle) (1 / 14) ↔
        |(w : ℝ) * x - round ((w : ℝ) * x)| < (1 : ℝ) / 14 := by
    simp only [LRCCommensuration.danger, LRCCommensuration.runnerMap,
      Set.mem_preimage, Metric.mem_ball, add_zero, dist_zero_right]
    rw [hpoint, UnitAddCircle.norm_eq]
  have hwitness :
      (x : UnitAddCircle) ∈
          LRCCommensuration.danger (w : ℤ) (0 : UnitAddCircle) (1 / 14) ↔
        ∃ a : ℤ, |(w : ℝ) * x - (a : ℝ)| < (1 : ℝ) / 14 := by
    rw [hdanger]
    constructor
    · intro h
      exact ⟨round ((w : ℝ) * x), h⟩
    · rintro ⟨a, ha⟩
      exact lt_of_le_of_lt (round_le ((w : ℝ) * x) a) ha
  rw [hwitness]
  constructor
  · intro hcomb
    obtain ⟨k, hk⟩ := Set.mem_iUnion.mp hcomb
    refine ⟨((k : ℕ) : ℤ), ?_⟩
    rw [abs_lt]
    have hlower :
        ((((k : ℕ) : ℝ) - (1 : ℝ) / 14) / (w : ℝ)) < x :=
      lt_of_le_of_lt (le_max_right _ _) hk.1
    have hupper :
        x < ((((k : ℕ) : ℝ) + (1 : ℝ) / 14) / (w : ℝ)) :=
      lt_of_lt_of_le hk.2 (min_le_right _ _)
    constructor
    · have := (div_lt_iff₀ hwR).mp hlower
      push_cast
      nlinarith
    · have := (lt_div_iff₀ hwR).mp hupper
      push_cast
      nlinarith
  · rintro ⟨a, ha⟩
    rw [abs_lt] at ha
    have hprodPos : (0 : ℝ) < (w : ℝ) * x := mul_pos hwR hx.1
    have hprodLt : (w : ℝ) * x < (w : ℝ) := by
      nlinarith [mul_lt_mul_of_pos_left hx.2 hwR]
    have haLowerR : (-1 : ℝ) < (a : ℝ) := by
      nlinarith [ha.2]
    have haUpperR : (a : ℝ) < (w : ℝ) + 1 := by
      nlinarith [ha.1]
    have haLowerZ : (-1 : ℤ) < a := by exact_mod_cast haLowerR
    have haUpperZ : a < (w : ℤ) + 1 := by exact_mod_cast haUpperR
    have ha0 : (0 : ℤ) ≤ a := by omega
    have haNatLt : a.toNat < w + 1 := by omega
    let k : Fin (w + 1) := ⟨a.toNat, haNatLt⟩
    have hkCast : (((k : ℕ) : ℝ)) = (a : ℝ) := by
      dsimp [k]
      have haCast : (a.toNat : ℤ) = a := Int.toNat_of_nonneg ha0
      exact_mod_cast haCast
    apply Set.mem_iUnion.mpr
    refine ⟨k, ?_, ?_⟩
    · unfold openCombLeft
      apply max_lt
      · exact hx.1
      · rw [hkCast]
        rw [div_lt_iff₀ hwR]
        nlinarith [ha.1]
    · unfold openCombRight
      apply lt_min
      · exact hx.2
      · rw [hkCast]
        rw [lt_div_iff₀ hwR]
        nlinarith [ha.2]

/-- Absolute-speed form used for arbitrary nonzero integer runner speeds. -/
theorem mem_openCombRegion_natAbs_iff_mem_danger
    (speed : ℤ) (hspeed : speed ≠ 0) (x : ℝ) (hx : x ∈ Set.Ioo (0 : ℝ) 1) :
    x ∈ openCombRegion speed.natAbs ↔
      (x : UnitAddCircle) ∈
        LRCCommensuration.danger speed (0 : UnitAddCircle) (1 / 14) := by
  rw [← danger_natAbs speed (1 / 14)]
  exact mem_openCombRegion_iff_mem_danger speed.natAbs
    (Int.natAbs_pos.mpr hspeed) x hx

#print axioms card_intervalOverlapEdges_le_add
#print axioms mem_openCombRegion_iff_mem_danger
#print axioms mem_openCombRegion_natAbs_iff_mem_danger

end
end LRCOpenDangerComb
end LonelyRunner
