/-
  TournamentH7.LRCThreeGapSampling -- the equally-spaced SAMPLING bound of the
  apex-ruler three-gap lemma (THM-565, kind-pasteur-2026-06-22-S34, THREAD 1).

  Part (2) of THM-565.  The apex `V` samples the slow coordinate `x` at the `V`
  equally-spaced points `x_j = (j + a)/V`, `j = 0,…,V-1`.  For the good set
  `G = ⊔_{i<m} A_i` a finite union of `m` arcs (THM-565 part (1)), the number of
  good samples is within `m` of `V·meas(G)`:

      #{ j : x_j ∈ G }  ∈  [ V·meas(G) - m , V·meas(G) + m ].

  The content is the per-interval lattice-point count: an interval of length `L`
  contains a number of `1/V`-spaced points that differs from `V·L` by `< 1`.
  This module proves that per-interval fact and the `m`-interval sum, sorry-free,
  as the elementary combinatorial core THM-527 Part A's equidistribution feeds on.

  We work with the half-open count `N(lo,hi) = ⌈V·hi - a⌉ - ⌈V·lo - a⌉`
  (= #{ j ∈ ℤ : lo ≤ (j+a)/V < hi }) and show `|N - V·(hi-lo)| ≤ 1`; summing the
  signed errors over `m` arcs gives the `±m` band.
-/

import Mathlib.Tactic

namespace LonelyRunner
namespace ThreeGapSampling

open scoped BigOperators

/-- The integer count of lattice points `(j+a)/V ∈ [lo, hi)`, written via ceilings:
`#{ j : lo ≤ (j+a)/V < hi } = ⌈V*hi - a⌉ - ⌈V*lo - a⌉` for `V > 0`. -/
noncomputable def latticeCount (V : ℕ) (a lo hi : ℝ) : ℤ :=
  ⌈(V : ℝ) * hi - a⌉ - ⌈(V : ℝ) * lo - a⌉

/-- **Per-interval lattice-count error bound.**  For `V ≥ 1`, the lattice count of an
interval `[lo,hi)` differs from `V·(hi-lo)` by at most `1` in absolute value.
This is the heart of part (2): equally-spaced points sample an interval up to a
single endpoint's worth of error. -/
theorem latticeCount_sub_abs_le_one (V : ℕ) (a lo hi : ℝ) (hV : 1 ≤ V) :
    |(latticeCount V a lo hi : ℝ) - (V : ℝ) * (hi - lo)| ≤ 1 := by
  have hVR : (1 : ℝ) ≤ (V : ℝ) := by exact_mod_cast hV
  set Hh := (V : ℝ) * hi - a with hHh
  set Hl := (V : ℝ) * lo - a with hHl
  have hdiff : Hh - Hl = (V : ℝ) * (hi - lo) := by rw [hHh, hHl]; ring
  -- ceil bounds:  x ≤ ⌈x⌉ < x + 1
  have hh1 : Hh ≤ (⌈Hh⌉ : ℝ) := Int.le_ceil Hh
  have hh2 : (⌈Hh⌉ : ℝ) < Hh + 1 := Int.ceil_lt_add_one Hh
  have hl1 : Hl ≤ (⌈Hl⌉ : ℝ) := Int.le_ceil Hl
  have hl2 : (⌈Hl⌉ : ℝ) < Hl + 1 := Int.ceil_lt_add_one Hl
  have hcount : (latticeCount V a lo hi : ℝ) = (⌈Hh⌉ : ℝ) - (⌈Hl⌉ : ℝ) := by
    rw [latticeCount]; push_cast; rw [hHh, hHl]
  rw [hcount]
  rw [abs_le]
  constructor
  · -- (⌈Hh⌉ - ⌈Hl⌉) - (Hh - Hl) ≥ -1  ⟸  ⌈Hh⌉ ≥ Hh and ⌈Hl⌉ < Hl + 1
    have : (⌈Hh⌉ : ℝ) - (⌈Hl⌉ : ℝ) ≥ Hh - (Hl + 1) := by linarith
    rw [hdiff] at *
    linarith
  · -- (⌈Hh⌉ - ⌈Hl⌉) - (Hh - Hl) ≤ 1  ⟸  ⌈Hh⌉ < Hh + 1 and ⌈Hl⌉ ≥ Hl
    have : (⌈Hh⌉ : ℝ) - (⌈Hl⌉ : ℝ) ≤ (Hh + 1) - Hl := by linarith
    rw [hdiff] at *
    linarith

/-- **The `m`-interval sampling band.**  If `G` is partitioned into `m` arcs with
lengths `L i` and the good-sample count is `Σ_i N_i` with each `N_i` within `1` of
`V·L i`, then the total count is within `m` of `V·Σ L i = V·meas(G)`.  Stated
abstractly over the per-arc counts and lengths. -/
theorem sum_count_sub_abs_le_card
    {ι : Type*} (I : Finset ι) (V : ℕ) (a : ℝ) (lo hi : ι → ℝ) (hV : 1 ≤ V) :
    |(∑ i ∈ I, (latticeCount V a (lo i) (hi i) : ℝ))
        - (V : ℝ) * (∑ i ∈ I, (hi i - lo i))| ≤ (I.card : ℝ) := by
  have hsplit :
      (∑ i ∈ I, (latticeCount V a (lo i) (hi i) : ℝ))
        - (V : ℝ) * (∑ i ∈ I, (hi i - lo i))
      = ∑ i ∈ I, ((latticeCount V a (lo i) (hi i) : ℝ) - (V : ℝ) * (hi i - lo i)) := by
    rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
  rw [hsplit]
  calc
    |∑ i ∈ I, ((latticeCount V a (lo i) (hi i) : ℝ) - (V : ℝ) * (hi i - lo i))|
        ≤ ∑ i ∈ I, |(latticeCount V a (lo i) (hi i) : ℝ) - (V : ℝ) * (hi i - lo i)| :=
          Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _i ∈ I, (1 : ℝ) := by
          apply Finset.sum_le_sum
          intro i _
          exact latticeCount_sub_abs_le_one V a (lo i) (hi i) hV
    _ = (I.card : ℝ) := by rw [Finset.sum_const, nsmul_eq_mul, mul_one]

/-- **The lower half = THM-565 part (2)'s `#good ≥ V·meas(G) - m`.**  From the band,
the good count is at least `V·meas(G) - m`.  This is the form used downstream: a
positive `V·c - m` forces `#good ≥ 1`. -/
theorem count_ge_measure_sub_card
    {ι : Type*} (I : Finset ι) (V : ℕ) (a : ℝ) (lo hi : ι → ℝ) (hV : 1 ≤ V) :
    (V : ℝ) * (∑ i ∈ I, (hi i - lo i)) - (I.card : ℝ)
      ≤ ∑ i ∈ I, (latticeCount V a (lo i) (hi i) : ℝ) := by
  have h := sum_count_sub_abs_le_card I V a lo hi hV
  have := (abs_le.mp h).1
  linarith

/-- **Positivity threshold (the Diophantine `V > m/c`).**  If `V·meas(G) > m`
(i.e. `V > m / meas(G) = V*`), then the good count is strictly positive, so at least
one apex period is good — the finite-`V` witness exists. -/
theorem count_pos_of_measure_gt_card
    {ι : Type*} (I : Finset ι) (V : ℕ) (a : ℝ) (lo hi : ι → ℝ) (hV : 1 ≤ V)
    (hgt : (I.card : ℝ) < (V : ℝ) * (∑ i ∈ I, (hi i - lo i))) :
    0 < ∑ i ∈ I, (latticeCount V a (lo i) (hi i) : ℝ) := by
  have h := count_ge_measure_sub_card I V a lo hi hV
  linarith

/-! ## Axiom audit -/

#print axioms latticeCount_sub_abs_le_one
#print axioms sum_count_sub_abs_le_card
#print axioms count_ge_measure_sub_card
#print axioms count_pos_of_measure_gt_card

end ThreeGapSampling
end LonelyRunner
