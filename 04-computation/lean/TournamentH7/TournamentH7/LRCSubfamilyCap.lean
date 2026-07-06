/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S115)
-/
import TournamentH7.LRCWitnessAttainment

/-!
# The subfamily cap: adding runners can only lower the margin (HYP-4476)

A new input for the bounded case.  Computationally (opus-S115), the loneliness constant `M` of a
*defected dilated AP* `c·{1,…,12}` with a fixed defect pattern is **height-independent** — it is
pinned to a rung (`1/12, 1/11, 2/19, …`, all `≥ 2/25`, or `1/13` for the pure AP), *identically*
across scales `c = 1, 2, …, 30`.  Increasing the height does not drift `M` toward `1/13`; the
family inherits its sub-AP subfamily's rung.  This is exactly why the bounded-height census
(kps-S32: 149k generalized APs, none in the gap) is **height-robust** — there is no large-height
escape into the window.

The provable core of that phenomenon is the **subfamily cap**: at any time, the margin of a
family is `≤` the margin of any subfamily.  Reindexing the speeds through any nonempty
`e : Fin m → Fin k` (a subfamily selector) can only *raise* the pointwise margin, because a min
over fewer terms is larger.  Hence `M(S) ≤ M(S')` for every subfamily `S' ⊆ S`: a family's
loneliness constant is capped by every one of its subfamilies' constants.  In particular a family
containing a dilated `m`-AP has `M ≤ 1/(m+1)` — the height-independent rung it lands on.
-/

namespace LonelyRunner
namespace SubfamilyCap

open TournamentH7.LRCWitness

variable {k m : ℕ} [Nonempty (Fin k)] [Nonempty (Fin m)]

/-- **Subfamily cap (pointwise).**  Reindexing the speeds through any `e : Fin m → Fin k` raises
the margin: `margin v t ≤ margin (v ∘ e) t`.  (A min over the sub-collection `{v (e i)}` is `≥`
the min over all of `{v j}`.) -/
theorem margin_le_comp (v : Fin k → ℤ) (e : Fin m → Fin k) (t : ℝ) :
    margin v t ≤ margin (v ∘ e) t := by
  unfold margin
  refine Finset.le_inf' Finset.univ_nonempty _ (fun i _ => ?_)
  exact Finset.inf'_le _ (Finset.mem_univ (e i))

/-- Every margin value is `≤ 1/2` (each `distZ` term is). -/
theorem margin_le_half (v : Fin k → ℤ) (t : ℝ) : margin v t ≤ 1 / 2 := by
  have half : ∀ y : ℝ, distZ y ≤ 1 / 2 := by
    intro y
    calc distZ y ≤ dist y ((round y : ℤ) : ℝ) :=
          Metric.infDist_le_dist_of_mem ⟨round y, rfl⟩
      _ = |y - (round y : ℝ)| := by rw [Real.dist_eq]
      _ ≤ 1 / 2 := abs_sub_round y
  obtain ⟨i⟩ := (inferInstance : Nonempty (Fin k))
  have hle : margin v t ≤ distZ ((v i : ℝ) * t) := by
    unfold margin; exact Finset.inf'_le _ (Finset.mem_univ i)
  exact hle.trans (half ((v i : ℝ) * t))

/-- **Subfamily cap (loneliness constant).**  `M(S) ≤ M(S')` for every subfamily selector
`e : Fin m → Fin k`: adding runners can only lower the loneliness constant `⨆_t margin`.  So a
family's `M` is bounded above by every subfamily's `M` — the mechanism by which a defected
dilated AP inherits its sub-AP's height-independent rung. -/
theorem iSup_margin_le_comp (v : Fin k → ℤ) (e : Fin m → Fin k) :
    (⨆ t : ℝ, margin v t) ≤ ⨆ t : ℝ, margin (v ∘ e) t := by
  have hbdd : BddAbove (Set.range fun t : ℝ => margin (v ∘ e) t) :=
    ⟨1 / 2, by rintro _ ⟨t, rfl⟩; exact margin_le_half (v ∘ e) t⟩
  exact ciSup_mono hbdd (fun t => margin_le_comp v e t)

#print axioms margin_le_comp
#print axioms iSup_margin_le_comp

end SubfamilyCap
end LonelyRunner
