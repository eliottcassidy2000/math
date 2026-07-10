/-
  TournamentH7.LRCWitnessG2Discharge — the FIRST DISCHARGED `hfloor` INSTANCES against the
  skeleton's now-CONCRETE `witnessG2` (death-star-2026-07-09-S4).

  With the de-opaquing of `LRCFourteenSkeleton.witnessG2` (:= `slowμ(goodSet E ∩ safeSet P)`)
  and `shapeOf` (absolute speeds ≤ 13 → `P`; the rest → co-offsets against `Vmax`), the
  witness-floor hypothesis `hfloor : witnessMP ≤ witnessG2 (shapeOf v)` becomes a PROVABLE
  statement instead of an assertion about an opaque token.  This file proves it on the first
  regime the certificate stack covers:

    · `witnessG2_pure_cluster` — with `P = []` the safe event is everything, so
      `witnessG2 ([], E) = slowμ(goodSet E)`;
    · `hfloor_pure_cluster_diam75` — for a nonempty co-offset list of diameter ≤ 75,
      `witnessMP ≤ witnessG2 ([], E)` (via `LRCGoodSetBridge` ← AP76 certificate ← Farey roof);
    · `hfloor_of_large_speeds` — the same against the skeleton's own `shapeOf`: any 13-family
      with all `|v i| > 13` and pairwise speed spread ≤ 75 satisfies the `hfloor` instance
      `witnessMP ≤ witnessG2 (shapeOf v)`.

  This is the first time the density-floor certificates (death-star-S2 AP20 → boxeph AP30/44 →
  opus-S145 AP76 → mac-mini TailDiameter → death-star-S3 GoodSetBridge) reach the skeleton's own
  vocabulary end-to-end.  Remaining for the full `hfloor`: the small-`P` (mixed) shapes — the
  Bonferroni leg (`LRCWitnessBonferroni` glue + Lemma B's `G_P` measure + Lemma A or the per-k
  moment floors) — and the diameter > 75 tail (LEM-005/006, exact arithmetic, unformalized).

  Kernel-pure target: no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCFourteenSkeleton
import TournamentH7.LRCGoodSetBridge

namespace LonelyRunner
namespace LRC14

open Set

/-- Every member bounds the `foldr max` from below. -/
theorem le_foldr_max {l : List ℤ} {a b : ℤ} (h : a ∈ l) : a ≤ l.foldr max b := by
  induction l with
  | nil => cases h
  | cons x xs ih =>
    rcases List.mem_cons.mp h with rfl | h
    · exact le_max_left _ _
    · exact le_trans (ih h) (le_max_right _ _)

/-- `foldr max` is bounded by any common upper bound of the seed and the members. -/
theorem foldr_max_le {l : List ℤ} {b c : ℤ} (hb : b ≤ c) (h : ∀ a ∈ l, a ≤ c) :
    l.foldr max b ≤ c := by
  induction l with
  | nil => exact hb
  | cons x xs ih =>
    exact max_le (h x List.mem_cons_self)
      (ih fun a ha => h a (List.mem_cons_of_mem x ha))

/-- **The empty small part imposes no constraint**: `safeSet [] = univ`. -/
theorem safeSet_nil : DenseCovers.safeSet ([] : List ℤ) = Set.univ := by
  ext x
  simp [DenseCovers.safeSet]

/-- **Pure clusters reduce to the bare GOOD measure**:
`witnessG2 ([], E) = slowμ(goodSet E)`. -/
theorem witnessG2_pure_cluster (E : List ℤ) :
    witnessG2 (([] : List ℤ), E)
      = (DenseCovers.slowμ (TournamentH7.GoodSet.goodSet E)).toReal := by
  unfold witnessG2
  rw [show ((([] : List ℤ), E) : Shape).1 = ([] : List ℤ) from rfl,
    show ((([] : List ℤ), E) : Shape).2 = E from rfl, safeSet_nil, Set.inter_univ]

/-- **The first discharged `hfloor` instance (pure cluster, diameter ≤ 75)**:
`witnessMP ≤ witnessG2 ([], E)` for any nonempty co-offset list within a window of
width ≤ 75.  Consumes the whole certificate chain: Farey roof → AP76 → TailDiameter
→ GoodSetBridge. -/
theorem hfloor_pure_cluster_diam75 (E : List ℤ) (hE : E.toFinset.Nonempty) (m D : ℤ)
    (hD75 : D ≤ 75) (hEd : ∀ e ∈ E.toFinset, e - m ∈ Finset.Icc (0 : ℤ) D) :
    (witnessMP : ℝ) ≤ witnessG2 (([] : List ℤ), E) := by
  rw [witnessG2_pure_cluster]
  have h := GoodSetBridge.slowμ_goodSet_toReal_ge_mP_diam75 E hE m D hD75 hEd
  have hcast : (witnessMP : ℝ) = (14249 : ℝ) / 252252 := by
    norm_num [witnessMP]
  rw [hcast]
  exact h

/-- **The `hfloor` instance against the skeleton's own `shapeOf`**: a 13-family whose
absolute speeds all exceed 13 (so the small part is empty) and whose pairwise speed
spread is ≤ 75 (so the co-offset cluster has diameter ≤ 75) satisfies the witness
floor `witnessMP ≤ witnessG2 (shapeOf v)` — the exact hypothesis shape
`lrc14_from_witness_floor` consumes, discharged on this regime. -/
theorem hfloor_of_large_speeds (v : Fin 13 → ℤ)
    (hlarge : ∀ i, 13 < |v i|)
    (hdiam : ∀ i j, |v i| ≤ |v j| + 75) :
    (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) := by
  classical
  set speeds : List ℤ := List.ofFn fun i => |v i| with hspeeds
  set Vmax : ℤ := speeds.foldr max 0 with hVmax
  -- the small part is empty: no absolute speed is ≤ 13
  have hPnil : speeds.filter (fun a => a ≤ 13) = [] := by
    rw [List.filter_eq_nil_iff]
    intro a ha
    rcases List.mem_ofFn.mp (hspeeds ▸ ha) with ⟨i, rfl⟩
    simpa using not_le.mpr (hlarge i)
  -- every absolute speed survives the cluster filter
  have hEfull : speeds.filter (fun a => 13 < a) = speeds := by
    rw [List.filter_eq_self]
    intro a ha
    rcases List.mem_ofFn.mp (hspeeds ▸ ha) with ⟨i, rfl⟩
    simpa using hlarge i
  -- the skeleton's shape of `v`, concretely
  have hshape : shapeOf v = (([] : List ℤ), speeds.map (fun a => Vmax - a)) := by
    simp only [shapeOf]
    rw [← hspeeds, ← hVmax, hPnil, hEfull]
  rw [hshape]
  -- the co-offset list: nonempty, and inside [0, 75]
  have hmem0 : |v 0| ∈ speeds := hspeeds ▸ List.mem_ofFn.mpr ⟨0, rfl⟩
  apply hfloor_pure_cluster_diam75 _ ?_ 0 75 le_rfl ?_
  · -- nonempty: the co-offset of runner 0 is a member
    exact ⟨Vmax - |v 0|, List.mem_toFinset.mpr (List.mem_map.mpr ⟨|v 0|, hmem0, rfl⟩)⟩
  · -- bounds: every co-offset Vmax − |v i| lies in [0, 75]
    intro e he
    rcases List.mem_map.mp (List.mem_toFinset.mp he) with ⟨a, ha, rfl⟩
    rcases List.mem_ofFn.mp (hspeeds ▸ ha) with ⟨i, rfl⟩
    have hup : |v i| ≤ Vmax := le_foldr_max (hspeeds ▸ List.mem_ofFn.mpr ⟨i, rfl⟩)
    have hdown : Vmax ≤ |v i| + 75 := by
      apply foldr_max_le
      · have := abs_nonneg (v i); linarith
      · intro a ha
        rcases List.mem_ofFn.mp (hspeeds ▸ ha) with ⟨j, rfl⟩
        exact hdiam j i
    simp only [Finset.mem_Icc, sub_zero]
    omega

/-! ## Axiom audit -/
#print axioms safeSet_nil
#print axioms witnessG2_pure_cluster
#print axioms hfloor_pure_cluster_diam75
#print axioms hfloor_of_large_speeds

end LRC14
end LonelyRunner
