/-
  TournamentH7.LRCClusterGapBrick — L1 OF THE CLUSTER-GAP ROUTE (THM-955)
  (death-star-2026-07-17-S45, on opus-S338's staged draft
  `04-computation/lean-drafts/LRCClusterGap.lean`; opus: "whoever takes L1
  gets the last elementary brick").

  `sorted_gap_pigeonhole`: removing `N` open rational intervals of total
  clipped length ≤ `B` from `[a, b]` leaves a closed subinterval of width
  ≥ (b − a − B)/(N + 1) avoiding every tooth.

  TWO DELIVERABLES:

  1. **The draft statement is FALSE without positivity** — teeth may overhang
     the window, and when `B ≥ b − a` a single tooth `(a−1, b+1)` covers
     `[a, b]` entirely: no avoiding pair `(a', b')` with `a ≤ a'`, `b' ≤ b`,
     `0 ≤ b' − a'` exists.  `gap_pigeonhole_positivity_necessary` proves the
     counterexample in-kernel.  The fix: hypothesis `0 < b − a − B` (which
     L3's width formula supplies whenever its bound is useful).

  2. **The corrected L1, proved** — not by the drafted mergeSort-partition
     plan but by a leaner *sorted-head induction*: sort teeth by left
     endpoint; the head tooth has minimal `t.1`, so `[a, c]` (with `c` the
     head's clamped left end) meets NO tooth.  Either `[a, c]` is already
     wide enough, or `c − a < (b−a−B)/(N+2)` and we excise the head's clamped
     span `[c, d]`, recursing on `[d, b]` with budget `B − (d − c)` and one
     fewer tooth; the mediant inequality
     (b−d−B')/(N+1) = (b−c−B)/(N+1) > (b−a−B)/(N+2)
     keeps the sharp constant.  The unsorted wrapper transfers along
     `List.perm_insertionSort` (mass and length are permutation-invariant).

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib

namespace LonelyRunner
namespace LRC14

/-- **The positivity hypothesis is necessary**: on `[0, 1]` with the single
overhanging tooth `(−1, 2)` (clipped mass `1 ≤ B = 1`), no avoiding pair
exists — the hypothesis-free draft form of L1 is false. -/
theorem gap_pigeonhole_positivity_necessary :
    ¬ (∃ a' b' : ℚ, (0 : ℚ) ≤ a' ∧ b' ≤ 1 ∧
        ((1 : ℚ) - 0 - 1) / (([((-1 : ℚ), (2 : ℚ))].length : ℚ) + 1) ≤ b' - a' ∧
        ∀ t ∈ [((-1 : ℚ), (2 : ℚ))], b' ≤ t.1 ∨ t.2 ≤ a') := by
  rintro ⟨a', b', ha, hb, hw, havoid⟩
  have h := havoid ((-1 : ℚ), (2 : ℚ)) (by simp)
  simp only at h
  norm_num at hw
  rcases h with h | h <;> linarith

/-- Auxiliary: the gap pigeonhole for teeth sorted by left endpoint, by
structural induction generalizing the window's left end and the budget. -/
theorem sorted_gap_aux (b : ℚ) (teeth : List (ℚ × ℚ)) :
    ∀ a B : ℚ, a ≤ b → 0 ≤ B →
      ((teeth.map fun t => max 0 (min t.2 b - max t.1 a)).sum ≤ B) →
      0 < b - a - B →
      teeth.Pairwise (fun s t : ℚ × ℚ => s.1 ≤ t.1) →
      ∃ a' b' : ℚ, a ≤ a' ∧ b' ≤ b ∧
        (b - a - B) / (teeth.length + 1) ≤ b' - a' ∧
        ∀ t ∈ teeth, b' ≤ t.1 ∨ t.2 ≤ a' := by
  induction teeth with
  | nil =>
      intro a B hab hB _hmass hpos _hsort
      refine ⟨a, b, le_rfl, le_rfl, ?_, by simp⟩
      simp only [List.length_nil, Nat.cast_zero, zero_add, div_one]
      linarith
  | cons t rest ih =>
      intro a B hab hB hmass hpos hsort
      obtain ⟨hhead, hrestsort⟩ := List.pairwise_cons.mp hsort
      simp only [List.map_cons, List.sum_cons] at hmass
      -- the head's clamped span [c, d] inside [a, b]
      have hca : a ≤ min (max t.1 a) b := le_min (le_max_right t.1 a) hab
      have hcb : min (max t.1 a) b ≤ b := min_le_right _ b
      set c : ℚ := min (max t.1 a) b with hc_def
      have hcd : c ≤ max (min t.2 b) c := le_max_right _ c
      have hdb : max (min t.2 b) c ≤ b := max_le (min_le_right t.2 b) hcb
      set d : ℚ := max (min t.2 b) c with hd_def
      -- the head's clipped mass dominates its clamped span
      have hrest_nonneg : 0 ≤ (rest.map fun t' => max 0 (min t'.2 b - max t'.1 a)).sum := by
        apply List.sum_nonneg
        intro x hx
        rcases List.mem_map.mp hx with ⟨t', _, rfl⟩
        exact le_max_left 0 _
      have hheadB : max 0 (min t.2 b - max t.1 a) ≤ B := by linarith
      have hdc_head : d - c ≤ max 0 (min t.2 b - max t.1 a) := by
        rcases le_total (min t.2 b) c with h | h
        · rw [hd_def, max_eq_right h]
          simp
        · rw [hd_def, max_eq_left h]
          rcases le_total (max t.1 a) b with h2 | h2
          · rw [hc_def, min_eq_left h2]
            exact le_max_right _ _
          · rw [hc_def, min_eq_right h2]
            have : min t.2 b - b ≤ 0 := by
              have := min_le_right t.2 b
              linarith
            exact le_trans this (le_max_left 0 _)
      -- normalize the goal's length cast
      simp only [List.length_cons, Nat.cast_add, Nat.cast_one]
      have hden : ((rest.length : ℚ) + 1) + 1 = (rest.length : ℚ) + 2 := by ring
      rw [hden]
      have hdenpos : (0 : ℚ) < (rest.length : ℚ) + 2 := by positivity
      rcases le_or_gt ((b - a - B) / ((rest.length : ℚ) + 2)) (c - a) with hwide | hnarrow
      · -- the tooth-free left segment [a, c] is already wide enough
        refine ⟨a, c, le_rfl, hcb, hwide, ?_⟩
        have hcpos : a < c := by
          have hWpos : 0 < (b - a - B) / ((rest.length : ℚ) + 2) := div_pos hpos hdenpos
          linarith
        have hct : c ≤ t.1 := by
          have hcm : c ≤ max t.1 a := by rw [hc_def]; exact min_le_left _ _
          rcases max_cases t.1 a with ⟨he, _⟩ | ⟨he, _⟩
          · rw [he] at hcm; exact hcm
          · rw [he] at hcm; linarith [hcpos]
        intro t' ht'
        rcases List.mem_cons.mp ht' with rfl | ht'
        · exact Or.inl hct
        · exact Or.inl (le_trans hct (hhead t' ht'))
      · -- excise the head's span; recurse on [d, b] with the reduced budget
        have hy : (c - a) * ((rest.length : ℚ) + 2) < b - a - B :=
          (lt_div_iff₀ hdenpos).mp hnarrow
        have hkey : 0 < b - c - B := by nlinarith [hpos, hdenpos, hca]
        have had : a ≤ d := le_trans hca hcd
        have hB₂ : 0 ≤ B - (d - c) := by linarith
        have hmass₂ : ((rest.map fun t' => max 0 (min t'.2 b - max t'.1 d)).sum)
            ≤ B - (d - c) := by
          have hmono : ((rest.map fun t' => max 0 (min t'.2 b - max t'.1 d)).sum)
              ≤ ((rest.map fun t' => max 0 (min t'.2 b - max t'.1 a)).sum) := by
            apply List.sum_le_sum
            intro t' _
            apply max_le_max le_rfl
            apply sub_le_sub_left
            exact max_le_max le_rfl had
          linarith
        have hpos₂ : 0 < b - d - (B - (d - c)) := by linarith [hkey]
        obtain ⟨a', b', ha'd, hb'b, hwidth₂, havoid₂⟩ :=
          ih d (B - (d - c)) hdb hB₂ hmass₂ hpos₂ hrestsort
        have hnum : b - d - (B - (d - c)) = b - c - B := by ring
        rw [hnum] at hwidth₂
        refine ⟨a', b', le_trans had ha'd, hb'b, ?_, ?_⟩
        · have hstep : (b - a - B) / ((rest.length : ℚ) + 2)
              ≤ (b - c - B) / ((rest.length : ℚ) + 1) := by
            rw [div_le_div_iff₀ (by positivity) (by positivity)]
            nlinarith [hy]
          exact le_trans hstep hwidth₂
        · intro t' ht'
          rcases List.mem_cons.mp ht' with rfl | ht'
          · -- the head tooth: its clipped right end sits left of a'
            right
            rcases le_total t'.2 b with h2b | hb2
            · have hmd : min t'.2 b ≤ d := by rw [hd_def]; exact le_max_left _ _
              rw [min_eq_left h2b] at hmd
              exact le_trans hmd ha'd
            · -- t'.2 > b would force the clipped tooth to reach b — but then
              -- its mass eats the whole remaining budget, contradicting hkey
              exfalso
              have hbc : b - c ≤ max 0 (min t'.2 b - max t'.1 a) := by
                rcases le_total (max t'.1 a) b with h | h
                · have hmm : min t'.2 b = b := min_eq_right hb2
                  have hcc : c = max t'.1 a := by rw [hc_def, min_eq_left h]
                  have hbr : b - max t'.1 a = min t'.2 b - max t'.1 a := by rw [hmm]
                  rw [hcc, hbr]
                  exact le_max_right _ _
                · have hcc : c = b := by rw [hc_def, min_eq_right h]
                  rw [hcc]
                  simp
              linarith [hheadB]
          · exact havoid₂ t' ht'

/-- **L1, THE SORTED-GAP PIGEONHOLE** (opus-S338's draft shape + the necessary
positivity): removing `N` open rational intervals of total clipped length
≤ `B < b − a` from `[a, b]` leaves a closed subinterval of width
≥ (b − a − B)/(N + 1) avoiding every tooth. -/
theorem sorted_gap_pigeonhole (a b B : ℚ) (hab : a ≤ b)
    (teeth : List (ℚ × ℚ)) (hB : 0 ≤ B)
    (hmass : (teeth.map fun t => max 0 (min t.2 b - max t.1 a)).sum ≤ B)
    (hpos : 0 < b - a - B) :
    ∃ a' b' : ℚ, a ≤ a' ∧ b' ≤ b ∧
      (b - a - B) / (teeth.length + 1) ≤ b' - a' ∧
      ∀ t ∈ teeth, b' ≤ t.1 ∨ t.2 ≤ a' := by
  classical
  haveI : DecidableRel (fun s t : ℚ × ℚ => s.1 ≤ t.1) :=
    fun s t => inferInstanceAs (Decidable (s.1 ≤ t.1))
  haveI : Std.Total (fun s t : ℚ × ℚ => s.1 ≤ t.1) := ⟨fun s t => le_total s.1 t.1⟩
  haveI : IsTrans (ℚ × ℚ) (fun s t : ℚ × ℚ => s.1 ≤ t.1) :=
    ⟨fun _ _ _ h1 h2 => le_trans h1 h2⟩
  have hperm := List.perm_insertionSort (fun s t : ℚ × ℚ => s.1 ≤ t.1) teeth
  have hsorted := List.pairwise_insertionSort (fun s t : ℚ × ℚ => s.1 ≤ t.1) teeth
  have hmass' : (((teeth.insertionSort (fun s t : ℚ × ℚ => s.1 ≤ t.1)).map
      fun t => max 0 (min t.2 b - max t.1 a)).sum ≤ B) := by
    rw [(hperm.map _).sum_eq]
    exact hmass
  obtain ⟨a', b', h1, h2, h3, h4⟩ :=
    sorted_gap_aux b (teeth.insertionSort (fun s t : ℚ × ℚ => s.1 ≤ t.1))
      a B hab hB hmass' hpos hsorted
  refine ⟨a', b', h1, h2, ?_, ?_⟩
  · rwa [hperm.length_eq] at h3
  · intro t ht
    exact h4 t (hperm.mem_iff.mpr ht)

/-! ## Axiom audit -/
#print axioms gap_pigeonhole_positivity_necessary
#print axioms sorted_gap_aux
#print axioms sorted_gap_pigeonhole

end LRC14
end LonelyRunner
