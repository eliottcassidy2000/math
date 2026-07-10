/-
  TournamentH7.LRCUnionVolume — brick (iii): the finite-union volume floor
  (mac-mini-2026-07-09-S65 cont.20).

  The last measure-theory brick of the witness-floor route: a SORTED-DISJOINT list of
  subintervals of [0,1) lying inside S forces `slowμ S ≥ ofReal (Σ lengths)`.  Composed with
  `Ico_subset_safeSet_of_bounds` (cont.18: each engine interval sits inside `safeSet P` by a
  checkable band condition), this discharges the FULL m_P floors of hsmall3/hlarge from the
  engine's merged interval lists (cont.16) — no exact set-identity needed: the ⊆ direction
  plus this floor is everything a lower bound requires.
-/
import Mathlib
import TournamentH7.LRCFourteenSkeleton

namespace LonelyRunner
namespace LRC14
namespace UnionVolume

open DenseCovers MeasureTheory

/-- Finite unions of `Ico`s over a list are measurable. -/
theorem measurable_listUnion_Ico (l : List (ℝ × ℝ)) :
    MeasurableSet (⋃ p ∈ l, Set.Ico p.1 p.2) := by
  induction l with
  | nil => simp
  | cons q t ih =>
      simp only [List.mem_cons, Set.iUnion_iUnion_eq_or_left]
      exact measurableSet_Ico.union ih

/-- **Brick (iii), the volume floor.**  A sorted-disjoint list of subintervals of `[0,1)`
inside `S` forces `slowμ S ≥ ofReal (Σ lengths)`. -/
theorem slowμ_ge_sum_of_sorted_Ico_subset {S : Set ℝ} (l : List (ℝ × ℝ))
    (hsub : ∀ p ∈ l, Set.Ico p.1 p.2 ⊆ S)
    (hin : ∀ p ∈ l, 0 ≤ p.1 ∧ p.1 ≤ p.2 ∧ p.2 ≤ 1)
    (hsorted : l.Pairwise (fun p q => p.2 ≤ q.1)) :
    ENNReal.ofReal ((l.map (fun p => p.2 - p.1)).sum) ≤ slowμ S := by
  have key : ∀ (m : List (ℝ × ℝ)),
      (∀ p ∈ m, 0 ≤ p.1 ∧ p.1 ≤ p.2 ∧ p.2 ≤ 1) →
      m.Pairwise (fun p q => p.2 ≤ q.1) →
      ENNReal.ofReal ((m.map (fun p => p.2 - p.1)).sum)
        ≤ slowμ (⋃ p ∈ m, Set.Ico p.1 p.2) := by
    intro m
    induction m with
    | nil => intro _ _; simp
    | cons q t ih =>
        intro hin' hpw'
        rw [List.pairwise_cons] at hpw'
        have hq := hin' q (by simp)
        have hint : ∀ p ∈ t, 0 ≤ p.1 ∧ p.1 ≤ p.2 ∧ p.2 ≤ 1 :=
          fun p hp => hin' p (by simp [hp])
        have hdisj : Disjoint (Set.Ico q.1 q.2) (⋃ p ∈ t, Set.Ico p.1 p.2) := by
          rw [Set.disjoint_iUnion₂_right]
          intro p hp
          rw [Set.disjoint_left]
          intro x hxq hxp
          exact absurd (lt_of_lt_of_le hxq.2 (hpw'.1 p hp)) (not_lt.mpr hxp.1)
        have hIco : slowμ (Set.Ico q.1 q.2) = ENNReal.ofReal (q.2 - q.1) := by
          unfold DenseCovers.slowμ
          rw [Measure.restrict_apply' measurableSet_Ico,
            Set.inter_eq_self_of_subset_left (Set.Ico_subset_Ico hq.1 hq.2.2)]
          exact Real.volume_Ico
        simp only [List.mem_cons, Set.iUnion_iUnion_eq_or_left, List.map_cons, List.sum_cons]
        rw [measure_union hdisj (measurable_listUnion_Ico t), hIco]
        calc ENNReal.ofReal ((q.2 - q.1) + (t.map (fun p => p.2 - p.1)).sum)
            ≤ ENNReal.ofReal (q.2 - q.1)
              + ENNReal.ofReal ((t.map (fun p => p.2 - p.1)).sum) := ENNReal.ofReal_add_le
          _ ≤ ENNReal.ofReal (q.2 - q.1) + slowμ (⋃ p ∈ t, Set.Ico p.1 p.2) := by
              gcongr
              exact ih hint hpw'.2
  refine le_trans (key l hin hsorted) (measure_mono ?_)
  exact Set.iUnion₂_subset hsub

end UnionVolume
end LRC14
end LonelyRunner
