/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-02-S44)
-/
import TournamentH7.RatIntervals
import TournamentH7.LRCSevenTranslate

/-!
# The far-element rate lemma over `Region` (the last unformalized statement, reduced)

klein's HYP-4001(b) rate lemma — the engine of the completeness peel — is, in `Region` terms:
clipping a region `A` against a fast comb removes at most density `2r` plus an end-effect
`≤ 4r/v` per component of `A`.  Per component `p = (a, b)`:

    Σ_{k < v} toothclip_k  ≤  2r·(b − a) + 4r/v.

**Proof architecture (this file):** each tooth `((k+φ−r)/v, (k+φ+r)/v)` is CONTAINED in the
cell `[c + k/v, c + (k+1)/v)` of the shifted uniform chain (`c = (φ−r)/v`), since `2r ≤ 1`.
Hence `toothclip_k ≤ min(2r/v, cellclip_k)`.  On FULL cells (`cellclip = 1/v`) the min equals
`2r · cellclip`; kps's `clip_chain_sum` evaluates `Σ cellclip` to the covered length `≤ b − a`.
The single genuinely new arithmetic fact is the **two-partial-cells bound** — an interval
meets at most two cells of a uniform chain partially (floor-uniqueness at each endpoint):

    Σ_{k < v} (min(2r/v, cellclip_k) − 2r·cellclip_k)  ≤  2 · (2r/v).             (hpartial)

This file proves the rate lemma FROM `hpartial` (sorry-free glue, kernel-pure) and
machine-checks `hpartial` instances (the pattern for its full proof or per-use discharge).
The rate lemma thus moves from "unformalized" to "one named finite-flavored sub-lemma".
-/

namespace LonelyRunner
namespace RatIntervals

/-- The clip length of interval `(a,b)` against the `k`-th tooth of `comb v r φ`. -/
def toothClip (v : ℕ) (r φ a b : ℚ) (k : ℕ) : ℚ :=
  max 0 (min b (((k : ℚ) + φ + r) / v) - max a (((k : ℚ) + φ - r) / v))

/-- The clip length of `(a,b)` against the `k`-th SHIFTED CHAIN cell
`[(φ−r)/v + k/v, (φ−r)/v + (k+1)/v)`. -/
def cellClip (v : ℕ) (r φ a b : ℚ) (k : ℕ) : ℚ :=
  max 0 (min b ((φ - r) / v + ((k : ℚ) + 1) / v) - max a ((φ - r) / v + (k : ℚ) / v))

/-- Tooth `k` is contained in shifted cell `k` (width `2r ≤ 1`), so its clip is dominated:
`toothClip ≤ min (2r/v) cellClip`. -/
theorem toothClip_le {v : ℕ} (hv : 0 < v) {r φ a b : ℚ} (hr : 0 ≤ r) (h2r : 2 * r ≤ 1)
    (k : ℕ) :
    toothClip v r φ a b k ≤ min (2 * r / v) (cellClip v r φ a b k) := by
  have hvQ : (0 : ℚ) < (v : ℚ) := by exact_mod_cast hv
  unfold toothClip cellClip
  apply le_min
  · -- width bound: the tooth has length 2r/v
    have h1 : min b (((k : ℚ) + φ + r) / v) - max a (((k : ℚ) + φ - r) / v)
        ≤ ((k : ℚ) + φ + r) / v - ((k : ℚ) + φ - r) / v :=
      sub_le_sub (min_le_right _ _) (le_max_right _ _)
    have h2 : ((k : ℚ) + φ + r) / v - ((k : ℚ) + φ - r) / v = 2 * r / v := by
      field_simp; ring
    calc max 0 (min b (((k : ℚ) + φ + r) / v) - max a (((k : ℚ) + φ - r) / v))
        ≤ max 0 (2 * r / v) := by
          apply max_le (le_max_left _ _)
          exact le_trans (h2 ▸ h1) (le_max_right _ _)
      _ = 2 * r / v := max_eq_right (by positivity)
  · -- containment bound: tooth k ⊆ cell k
    apply max_le (le_max_left _ _)
    have hlo : (φ - r) / v + (k : ℚ) / v ≤ ((k : ℚ) + φ - r) / v :=
      le_of_eq (by field_simp; ring)
    have hhi : ((k : ℚ) + φ + r) / v ≤ (φ - r) / v + ((k : ℚ) + 1) / v := by
      have e : (φ - r) / v + ((k : ℚ) + 1) / v - ((k : ℚ) + φ + r) / v = (1 - 2 * r) / v := by
        field_simp; ring
      have hpos : (0 : ℚ) ≤ (1 - 2 * r) / v := div_nonneg (by linarith) hvQ.le
      linarith
    have := sub_le_sub (min_le_min (le_refl b) hhi) (max_le_max (le_refl a) hlo)
    exact le_trans this (le_max_right _ _)

/-- **The far-element rate lemma over `Region`, from the two-partial-cells bound.**
Given `hpartial` (the single remaining arithmetic sub-lemma: partial cells contribute at most
two end-corrections), one component of a region loses at most density `2r` plus `4r/v` to a
speed-`v` comb. -/
theorem rate_lemma_component {v : ℕ} (hv : 0 < v) {r φ a b : ℚ} (hr : 0 ≤ r)
    (h2r : 2 * r ≤ 1) (hab : a ≤ b)
    (hcover : ((List.range v).map fun k => cellClip v r φ a b k).sum ≤ b - a)
    (hpartial : ((List.range v).map fun k =>
        min (2 * r / v) (cellClip v r φ a b k) - 2 * r * cellClip v r φ a b k).sum
      ≤ 2 * (2 * r / v)) :
    ((List.range v).map fun k => toothClip v r φ a b k).sum
      ≤ 2 * r * (b - a) + 4 * r / v := by
  have h1 : ((List.range v).map fun k => toothClip v r φ a b k).sum
      ≤ ((List.range v).map fun k => min (2 * r / v) (cellClip v r φ a b k)).sum := by
    apply List.sum_le_sum
    intro k hk
    simp only [List.mem_map] at *
    exact toothClip_le hv hr h2r k
  have h2 : ((List.range v).map fun k => min (2 * r / v) (cellClip v r φ a b k)).sum
      = ((List.range v).map fun k =>
          min (2 * r / v) (cellClip v r φ a b k) - 2 * r * cellClip v r φ a b k).sum
        + 2 * r * ((List.range v).map fun k => cellClip v r φ a b k).sum := by
    induction (List.range v) with
    | nil => simp
    | cons x L ih => simp only [List.map_cons, List.sum_cons, ih]; ring
  have h3 : 2 * r * ((List.range v).map fun k => cellClip v r φ a b k).sum
      ≤ 2 * r * (b - a) := by
    apply mul_le_mul_of_nonneg_left hcover (by positivity)
  calc ((List.range v).map fun k => toothClip v r φ a b k).sum
      ≤ _ := h1
    _ = _ := h2
    _ ≤ 2 * (2 * r / v) + 2 * r * (b - a) := add_le_add hpartial h3
    _ = 2 * r * (b - a) + 4 * r / v := by ring

/-! ### Dischargeability of the two hypotheses (machine-checked instances)

Both `hcover` and `hpartial` are decidable per instance.  `hcover` is moreover kps's
`clip_chain_sum` in inequality form (the chain-clip sum is EXACTLY the covered length, which
is `≤ b − a`); `hpartial`'s general proof is the floor-uniqueness two-endpoint argument.
Instances below establish the discharge pattern at the working parameters. -/

theorem hcover_instance :
    ((List.range 50).map fun k => cellClip 50 (1/14) (399/4000) (7/20) (3/8) k).sum
      ≤ 3/8 - 7/20 := by native_decide

theorem hpartial_instance :
    ((List.range 50).map fun k =>
        min (2 * (1/14) / 50) (cellClip 50 (1/14) (399/4000) (7/20) (3/8) k)
          - 2 * (1/14) * cellClip 50 (1/14) (399/4000) (7/20) (3/8) k).sum
      ≤ 2 * (2 * (1/14) / 50) := by native_decide

/-- The assembled instance: the rate bound at the S36 level-1 parameters. -/
theorem rate_lemma_instance :
    ((List.range 50).map fun k => toothClip 50 (1/14) (399/4000) (7/20) (3/8) k).sum
      ≤ 2 * (1/14) * (3/8 - 7/20) + 4 * (1/14) / 50 :=
  rate_lemma_component (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    hcover_instance hpartial_instance

end RatIntervals
end LonelyRunner
