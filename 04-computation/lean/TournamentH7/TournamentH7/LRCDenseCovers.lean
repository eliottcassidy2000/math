/-
  TournamentH7.LRCDenseCovers -- the elementary pointwise inclusion at the heart
  of the witness/p0 UNIFICATION (kind-pasteur-2026-06-22-S30, HYP-2832).

  This is the combinatorial core of `D(E) <= p0(E)`:

    if a finite set of circle-phases `S ⊆ [0,1)` (with the anchor `0 ∈ S`) has
    every circular gap `≤ 1/7` ("1/7-dense"), then every INNER sector
    `[j/7, (j+1)/7)`, `j = 1,…,6`, contains a phase.

  Pointwise this is `{x : phases 1/7-dense} ⊆ {x : all inner sectors hit}`, hence
  (by measure monotonicity, once the measures are defined) `D(E) ≤ p0(E)` — the
  inclusion that makes the 1/7-witness floor a corollary of the p0 wide bound.

  The proof is the half-open-sector argument: if sector `S_j` held no phase, the
  gap straddling it would exceed `1/7`; the boundary case `= 1/7` is excluded
  because a phase at `j/7` lies IN `S_j` (left-closed), and the wrap-around case
  is excluded by the anchor `0 ∈ S`.  Sorry-free, no measure theory.
-/

import Mathlib.Tactic
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic

namespace LonelyRunner
namespace DenseCovers

open Finset MeasureTheory

/-- `S` is **1/7-dense**: every adjacent circular gap is `≤ 1/7`.  First clause:
each consecutive pair `a < b` with no phase strictly between has `b - a ≤ 1/7`.
Second clause (wrap): the maximal phase `b` has `1 - b ≤ 1/7` (the wrap gap from
`b` back to the anchor `0`). -/
def Dense17 (S : Finset ℝ) : Prop :=
  (∀ a ∈ S, ∀ b ∈ S, a < b → (∀ c ∈ S, a < c → b ≤ c) → b - a ≤ 1 / 7) ∧
  (∀ b ∈ S, (∀ c ∈ S, c ≤ b) → 1 - b ≤ 1 / 7)

/-- **The pointwise inclusion `D ⊆ p0`.**  If `S ⊆ [0,1)` is 1/7-dense and the
anchor `0 ∈ S`, then every inner sector `[j/7, (j+1)/7)` (`1 ≤ j ≤ 6`) contains a
phase of `S`. -/
theorem inner_sector_covered
    (S : Finset ℝ) (_hsub : ∀ s ∈ S, 0 ≤ s ∧ s < 1) (h0 : (0 : ℝ) ∈ S)
    (hdense : Dense17 S) (j : ℕ) (hj1 : 1 ≤ j) (hj6 : j ≤ 6) :
    ∃ s ∈ S, (j : ℝ) / 7 ≤ s ∧ s < ((j : ℝ) + 1) / 7 := by
  obtain ⟨hadj, hwrap⟩ := hdense
  -- abbreviations for the two sector endpoints
  set jr : ℝ := (j : ℝ) / 7 with hjr
  set jr1 : ℝ := ((j : ℝ) + 1) / 7 with hjr1
  have hj1R : (1 : ℝ) ≤ (j : ℝ) := by exact_mod_cast hj1
  have hj6R : (j : ℝ) ≤ 6 := by exact_mod_cast hj6
  have hjr_pos : 0 < jr := by rw [hjr]; positivity
  have hjr_lt_jr1 : jr < jr1 := by rw [hjr, hjr1]; linarith
  -- L = phases below the right endpoint jr1; nonempty since 0 ∈ L (0 < jr1)
  have h0_lt_jr1 : (0 : ℝ) < jr1 := lt_trans hjr_pos hjr_lt_jr1
  have h0L : (0 : ℝ) ∈ S.filter (· < jr1) := by
    rw [mem_filter]; exact ⟨h0, h0_lt_jr1⟩
  have hLne : (S.filter (· < jr1)).Nonempty := ⟨0, h0L⟩
  set a : ℝ := (S.filter (· < jr1)).max' hLne with ha
  have haL : a ∈ S.filter (· < jr1) := max'_mem _ hLne
  have haS : a ∈ S := (mem_filter.mp haL).1
  have ha_lt_jr1 : a < jr1 := (mem_filter.mp haL).2
  -- it suffices to show jr ≤ a (then a is in the sector)
  refine ⟨a, haS, ?_, ha_lt_jr1⟩
  by_contra hlt
  push Not at hlt          -- hlt : a < jr
  -- CASE on whether some phase exceeds a
  by_cases hU : ∃ c ∈ S, a < c
  · -- let b be the least phase > a; it is the adjacent successor of a
    obtain ⟨c0, hc0S, hc0⟩ := hU
    -- the successor set U = {c ∈ S : a < c} is nonempty
    have hUne : (S.filter (a < ·)).Nonempty := ⟨c0, by rw [mem_filter]; exact ⟨hc0S, hc0⟩⟩
    set b : ℝ := (S.filter (a < ·)).min' hUne with hb
    have hbU : b ∈ S.filter (a < ·) := min'_mem _ hUne
    have hbS : b ∈ S := (mem_filter.mp hbU).1
    have ha_lt_b : a < b := (mem_filter.mp hbU).2
    -- b is the immediate successor: any phase > a is ≥ b
    have hsucc : ∀ c ∈ S, a < c → b ≤ c := by
      intro c hcS hac
      exact min'_le _ c (by rw [mem_filter]; exact ⟨hcS, hac⟩)
    -- adjacent gap bound: b - a ≤ 1/7
    have hgap : b - a ≤ 1 / 7 := hadj a haS b hbS ha_lt_b hsucc
    -- b ≥ jr1 : since a is the max of {< jr1}, anything > a is ≥ jr1
    have hb_ge : jr1 ≤ b := by
      by_contra hbb
      push Not at hbb       -- b < jr1
      have hbL : b ∈ S.filter (· < jr1) := by rw [mem_filter]; exact ⟨hbS, hbb⟩
      have : b ≤ a := le_max' _ b hbL
      linarith
    -- but a < jr and b ≥ jr1 = jr + 1/7 give b - a > 1/7, contradiction
    have hjr1_eq : jr1 = jr + 1 / 7 := by rw [hjr, hjr1]; ring
    linarith
  · -- no phase exceeds a, so a is the maximum; the wrap clause applies
    have hamax : ∀ c ∈ S, c ≤ a := by
      intro c hcS
      by_contra hc
      push Not at hc            -- hc : a < c
      exact hU ⟨c, hcS, hc⟩
    have hwrapa : 1 - a ≤ 1 / 7 := hwrap a haS hamax
    -- so a ≥ 6/7 ≥ jr, contradicting a < jr
    have hjr_le : jr ≤ 6 / 7 := by rw [hjr]; linarith
    linarith

/-- **The packaged `D ⊆ p0` pointwise inclusion.**  A 1/7-dense phase set (with the
anchor `0 ∈ S`) hits *every* inner sector simultaneously — exactly the `p0`/`measS7`
cover event.  Pointwise this is `{1/7-dense} ⊆ {all 6 inner sectors hit}`; the
measure form `D(E) ≤ p0(E)` follows by monotonicity once the two events are given
their Lebesgue measures.  This is the elementary inclusion that makes the
1/7-witness floor a corollary of the p0 wide bound (HYP-2832, the UNIFICATION). -/
theorem dense_covers_all_inner
    (S : Finset ℝ) (hsub : ∀ s ∈ S, 0 ≤ s ∧ s < 1) (h0 : (0 : ℝ) ∈ S)
    (hdense : Dense17 S) :
    ∀ j : ℕ, 1 ≤ j → j ≤ 6 → ∃ s ∈ S, (j : ℝ) / 7 ≤ s ∧ s < ((j : ℝ) + 1) / 7 := by
  intro j hj1 hj6
  exact inner_sector_covered S hsub h0 hdense j hj1 hj6

/-! ## The measure-level inclusion `D(E) ≤ p0(E)`

Lifting the pointwise inclusion to slow-time measure.  `E : List ℤ` is the cluster
co-offsets; the phase set at slow time `x` is `{frac(e x) : e ∈ E}`.  The dense set
`{x : phases 1/7-dense}` is contained in the cover set `{x : all inner sectors hit}`,
so `volume(dense) ≤ volume(cover)` — the measure form of `D(E) ≤ p0(E)`.  Only the
*inclusion* is needed (`measure_mono` is monotone on arbitrary sets), so no
measurability of the dense set is required. -/

/-- The phase finset of `E` at slow time `x`: `{frac(e x) : e ∈ E}`. -/
noncomputable def phaseFinset (E : List ℤ) (x : ℝ) : Finset ℝ :=
  (E.map (fun e => Int.fract ((e : ℝ) * x))).toFinset

/-- The **cover set** `{x : the phases hit every inner sector}` (`= measS7` support). -/
def coverSet (E : List ℤ) : Set ℝ :=
  {x | ∀ j : ℕ, 1 ≤ j → j ≤ 6 →
      ∃ e ∈ E, (j : ℝ) / 7 ≤ Int.fract ((e : ℝ) * x) ∧
        Int.fract ((e : ℝ) * x) < ((j : ℝ) + 1) / 7}

/-- The **dense set** `{x : the phases are 1/7-dense}`. -/
def denseSet (E : List ℤ) : Set ℝ := {x | Dense17 (phaseFinset E x)}

/-- Membership in the phase finset is "is the fractional part of some `e*x`". -/
theorem mem_phaseFinset {E : List ℤ} {x s : ℝ} :
    s ∈ phaseFinset E x ↔ ∃ e ∈ E, Int.fract ((e : ℝ) * x) = s := by
  simp [phaseFinset]

/-- The anchor `0 ∈ E` puts the phase `0` in every phase finset. -/
theorem zero_mem_phaseFinset {E : List ℤ} (hE : (0 : ℤ) ∈ E) (x : ℝ) :
    (0 : ℝ) ∈ phaseFinset E x := by
  rw [mem_phaseFinset]
  refine ⟨0, hE, ?_⟩
  rw [Int.cast_zero, zero_mul, Int.fract_zero]

/-- Every phase lies in `[0,1)`. -/
theorem phaseFinset_mem_Ico {E : List ℤ} {x s : ℝ} (hs : s ∈ phaseFinset E x) :
    0 ≤ s ∧ s < 1 := by
  rw [mem_phaseFinset] at hs
  obtain ⟨e, _, he⟩ := hs
  exact ⟨he ▸ Int.fract_nonneg _, he ▸ Int.fract_lt_one _⟩

/-- **Pointwise `D ⊆ p0`.**  If the phases are 1/7-dense at `x` (and `0 ∈ E`), every
inner sector is hit at `x`. -/
theorem denseSet_subset_coverSet (E : List ℤ) (hE : (0 : ℤ) ∈ E) :
    denseSet E ⊆ coverSet E := by
  intro x hx j hj1 hj6
  obtain ⟨s, hsS, hs1, hs2⟩ :=
    inner_sector_covered (phaseFinset E x)
      (fun s hs => phaseFinset_mem_Ico hs) (zero_mem_phaseFinset hE x) hx j hj1 hj6
  obtain ⟨e, heE, he⟩ := mem_phaseFinset.mp hsS
  exact ⟨e, heE, he ▸ hs1, he ▸ hs2⟩

/-- **The measure-level inclusion `D(E) ≤ p0(E)`.**  Restricted to ONE period
`[0,1)` (the sets are 1-periodic, so the unrestricted `volume` is infinite and
vacuous): `volume(denseSet ∩ [0,1)) ≤ volume(coverSet ∩ [0,1))`.  This is the
genuine slow-time measure form of `D(E) ≤ p0(E)` (`D, p0 ∈ [0,1]`), the elementary
inclusion behind the witness/p0 unification (HYP-2832).  Proved by monotonicity of
Lebesgue measure on the pointwise inclusion (no measurability of `denseSet` needed). -/
theorem volume_denseSet_le_coverSet (E : List ℤ) (hE : (0 : ℤ) ∈ E) :
    volume (denseSet E ∩ Set.Ico (0 : ℝ) 1) ≤ volume (coverSet E ∩ Set.Ico (0 : ℝ) 1) :=
  measure_mono (Set.inter_subset_inter_left _ (denseSet_subset_coverSet E hE))

/-- Sanity: the restricted cover measure is a genuine probability mass in `[0,1]`
(`p0 ≤ 1`), so the inclusion above is a non-vacuous bound. -/
theorem volume_coverSet_inter_le_one (E : List ℤ) :
    volume (coverSet E ∩ Set.Ico (0 : ℝ) 1) ≤ 1 := by
  calc volume (coverSet E ∩ Set.Ico (0 : ℝ) 1)
      ≤ volume (Set.Ico (0 : ℝ) 1) := measure_mono (Set.inter_subset_right)
    _ = 1 := by rw [Real.volume_Ico]; norm_num

/-! ## Axiom audit -/

#print axioms inner_sector_covered
#print axioms dense_covers_all_inner
#print axioms denseSet_subset_coverSet
#print axioms volume_denseSet_le_coverSet

end DenseCovers
end LonelyRunner
