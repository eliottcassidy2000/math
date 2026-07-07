/-
  TournamentH7.LRCFareyRoofBridge — bridging the Farey roof (opus-S135, THM-637) to
  the tail-diameter good set (mac-mini-S42 / monad-S2 HYP-4817).  boxeph-2026-07-07.

  `LRCFareyRoof.zero_gap_empty` proves that on the open Farey-k cell `(p/q, p'/q')`
  no config point of the AP `{1,…,k}` lies strictly inside the open interval
  `(q'x − p', qx − p)`, which contains `0` and has length
  `roof := (qx − p) + (p' − q'x)`.

  `LRCTailDiameter.Good θ E = {x | ∃ a, ∀ e ∈ E, θ < fract(e·x − a)}` is the good set
  (some closed θ-arc is empty of the orbit); its measure `muGood` is the density
  `μ_θ`, and the whole diameter route is GREEN conditional on the AP₇₆ certificate
  `muGood (1/7) {0..75} ≥ 2314528732/40290957525`.

  THIS FILE supplies the missing pointwise link:

    roof(x) > θ   ⟹   x ∈ Good θ (Finset.Icc 1 k).

  Consequence: `{x in a Farey-k cell : roof > θ} ⊆ Good θ (AP_k)`, so
  `muGood θ (AP_k) ≥ meas{ roof > θ }` — reducing the certificate to a PURE
  real-superlevel (Farey-sum) measure computation, with no orbit reasoning left.

  Proof idea: put the θ-arc's left end at `a := (q'x − p') + (roof − θ)/2`, the
  midpoint offset that makes the closed θ-arc `[a, a+θ]` sit STRICTLY inside the
  empty roof-interval.  Any orbit point in `(a, a+θ]` would then land strictly
  inside `(q'x − p', qx − p)`, contradicting `zero_gap_empty`.

  Kernel-pure: no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LRCFareyRoof
import TournamentH7.LRCTailDiameter

namespace LonelyRunner
namespace FareyRoofBridge

open TailDiameter
open MeasureTheory
open scoped ENNReal

/-- **The roof→good bridge.**  On the open Farey-`k` cell `(p/q, p'/q')` (encoded in
cleared form `p < q·x`, `q'·x < p'`, determinant `q·p' − p·q' = 1`, `k < q + q'`),
if the roof `(q·x − p) + (p' − q'·x)` exceeds `θ ≥ 0`, then `x` is in the AP-`k`
good set: some closed θ-arc is empty of the orbit `{frac(i·x) : 1 ≤ i ≤ k}`. -/
theorem good_of_roof_gt {p q p' q' k : ℤ} {x θ : ℝ}
    (hq : 0 < q) (hq' : 0 < q') (hdet : q * p' - p * q' = 1) (hsum : k < q + q')
    (hx : (p : ℝ) < q * x) (hx' : (q' : ℝ) * x < p')
    (hroof : θ < ((q : ℝ) * x - p) + ((p' : ℝ) - q' * x)) :
    x ∈ Good θ (Finset.Icc (1 : ℤ) k) := by
  -- roof interval endpoints  a0 = q'x − p' < 0 < b0 = qx − p
  set a0 : ℝ := (q' : ℝ) * x - p' with ha0
  set b0 : ℝ := (q : ℝ) * x - p with hb0
  have hroof' : θ < b0 - a0 := by rw [ha0, hb0]; linarith [hroof]
  -- witness: left end of the θ-arc, pushed in by half the slack
  refine ⟨a0 + (b0 - a0 - θ) / 2, ?_⟩
  intro e he
  rcases Finset.mem_Icc.mp he with ⟨he1, hek⟩
  set a : ℝ := a0 + (b0 - a0 - θ) / 2 with ha
  by_contra hcon
  have hcon' : Int.fract ((e : ℝ) * x - a) ≤ θ := not_lt.mp hcon
  have hf0 : (0 : ℝ) ≤ Int.fract ((e : ℝ) * x - a) := Int.fract_nonneg _
  -- n := −⌊e·x − a⌋ realizes  e·x + n = a + fract(e·x − a)
  set n : ℤ := -⌊(e : ℝ) * x - a⌋ with hn
  have hfl : ((⌊(e : ℝ) * x - a⌋ : ℤ) : ℝ) + Int.fract ((e : ℝ) * x - a)
      = (e : ℝ) * x - a := Int.floor_add_fract _
  have hval : (e : ℝ) * x + (n : ℝ) = a + Int.fract ((e : ℝ) * x - a) := by
    rw [hn]; push_cast; linarith [hfl]
  -- a + fract ∈ (a0, b0): left end from slack>0 and fract≥0; right end from fract≤θ
  have hlo : a0 < (e : ℝ) * x + (n : ℝ) := by rw [hval, ha]; nlinarith [hf0, hroof']
  have hhi : (e : ℝ) * x + (n : ℝ) < b0 := by rw [hval, ha]; nlinarith [hcon', hroof']
  exact FareyRoof.zero_gap_empty hq hq' hdet hsum he1 hek hx hx' n ⟨hlo, hhi⟩

/-- Set form: within a fixed Farey-`k` cell, the roof-superlevel set is contained
in the AP-`k` good set.  (The cell hypotheses are carried as a predicate on `x`.) -/
theorem roof_superlevel_subset_good
    {p q p' q' k : ℤ} {θ : ℝ}
    (hq : 0 < q) (hq' : 0 < q') (hdet : q * p' - p * q' = 1) (hsum : k < q + q') :
    {x : ℝ | (p : ℝ) < q * x ∧ (q' : ℝ) * x < p' ∧
        θ < ((q : ℝ) * x - p) + ((p' : ℝ) - q' * x)} ⊆ Good θ (Finset.Icc (1 : ℤ) k) := by
  rintro x ⟨hx, hx', hroof⟩
  exact good_of_roof_gt hq hq' hdet hsum hx hx' hroof

/-- **Bridge → measure atom.**  If the open interval `(c, d) ⊆ [0,1]` lies inside the
Farey-`k` cell and the roof exceeds `θ` throughout it, then that interval's length is a
lower bound for the good-set measure `μ_θ(AP_k)`.  Summing this atom over the finitely
many roof-superlevel intervals (only cells adjacent to a `q ≤ 6` node contribute, since
the roof node values are `1/q`) reduces the `AP₇₆` certificate to a pure interval-length
arithmetic — exactly the Farey-cell sum the fleet computed numerically. -/
theorem muGood_ge_of_cell_interval {p q p' q' k : ℤ} {θ c d : ℝ}
    (hq : 0 < q) (hq' : 0 < q') (hdet : q * p' - p * q' = 1) (hsum : k < q + q')
    (hc0 : 0 ≤ c) (hd1 : d ≤ 1)
    (hcell : ∀ x ∈ Set.Ioo c d, (p : ℝ) < q * x ∧ (q' : ℝ) * x < p')
    (hroof : ∀ x ∈ Set.Ioo c d, θ < ((q : ℝ) * x - p) + ((p' : ℝ) - q' * x)) :
    ENNReal.ofReal (d - c) ≤ muGood θ (Finset.Icc (1 : ℤ) k) := by
  have hsub : Set.Ioo c d ⊆ Good θ (Finset.Icc (1 : ℤ) k) ∩ Set.Icc (0 : ℝ) 1 := by
    intro x hx
    refine ⟨?_, ?_⟩
    · obtain ⟨hx1, hx2⟩ := hcell x hx
      exact good_of_roof_gt hq hq' hdet hsum hx1 hx2 (hroof x hx)
    · rcases hx with ⟨hcx, hxd⟩
      exact ⟨le_of_lt (lt_of_le_of_lt hc0 hcx), le_of_lt (lt_of_lt_of_le hxd hd1)⟩
  calc ENNReal.ofReal (d - c)
      = volume (Set.Ioo c d) := (Real.volume_Ioo).symm
    _ ≤ volume (Good θ (Finset.Icc (1 : ℤ) k) ∩ Set.Icc (0 : ℝ) 1) := measure_mono hsub
    _ = muGood θ (Finset.Icc (1 : ℤ) k) := rfl

/-- **The aggregator.**  A finite family of pairwise-disjoint open intervals, each contained
in `Good θ E ∩ [0,1]`, contributes the SUM of its lengths as a lower bound for `μ_θ(E)`.
Instantiated with the roof-superlevel intervals of the contributing Farey cells (each
membership discharged by `good_of_roof_gt`), this turns the `AP₇₆` certificate into a single
decidable interval-length sum. -/
theorem muGood_ge_sum_intervals {θ : ℝ} {E : Finset ℤ} {ι : Type*} (s : Finset ι)
    (c d : ι → ℝ)
    (hsub : ∀ i ∈ s, Set.Ioo (c i) (d i) ⊆ Good θ E ∩ Set.Icc (0 : ℝ) 1)
    (hdisj : (↑s : Set ι).PairwiseDisjoint (fun i => Set.Ioo (c i) (d i))) :
    ∑ i ∈ s, ENNReal.ofReal (d i - c i) ≤ muGood θ E := by
  have hunion : (⋃ i ∈ s, Set.Ioo (c i) (d i)) ⊆ Good θ E ∩ Set.Icc (0 : ℝ) 1 := by
    simpa only [Set.iUnion_subset_iff] using hsub
  have hmeas : volume (⋃ i ∈ s, Set.Ioo (c i) (d i)) = ∑ i ∈ s, volume (Set.Ioo (c i) (d i)) :=
    measure_biUnion_finset hdisj (fun i _ => measurableSet_Ioo)
  calc ∑ i ∈ s, ENNReal.ofReal (d i - c i)
      = ∑ i ∈ s, volume (Set.Ioo (c i) (d i)) :=
        Finset.sum_congr rfl (fun i _ => (Real.volume_Ioo).symm)
    _ = volume (⋃ i ∈ s, Set.Ioo (c i) (d i)) := hmeas.symm
    _ ≤ volume (Good θ E ∩ Set.Icc (0 : ℝ) 1) := measure_mono hunion
    _ = muGood θ E := rfl

end FareyRoofBridge
end LonelyRunner
