/-
  TournamentH7.LRCArcCount — the ARC-COUNT PIGEONHOLE for good-period existence
  (opus-2026-07-09-S169).

  The finite-`Vmax` glue of the LRC(14) covering capstone (mac-mini-S58/S61, opus-S167/S168).
  A *good period* is an integer `j` with `j/V ∈ Good θ E` (some θ-arc empty of the orbit at the
  rational phase `j/V`).  The dissociated branch closes by EXISTENCE of a good period, obtained
  from the pigeonhole

      (#arcs of `Good θ E ∩ [0,1)`)  <  ρ*·V        ⟹        a good period exists,

  where `ρ* = vol(Good θ E ∩ [0,1))`.  Concretely: if the good set contains `N` intervals of
  total length `> N/V`, one of them is longer than `1/V`, hence catches a grid point `j/V`.

  This file proves that pigeonhole in full (kernel-pure).  The two inputs are supplied elsewhere:
    • `ρ* ≥ D3 ≥ bar > 0` — the density floor (THM-661 / klein finite `D3_c` table), and
    • `#arcs ≤ c·spread` with `c < ρ*` — the bounded-arc-count bound (mac-mini-S58; the trivial
      `O(k² spread)` bound is `~10³×` too loose, so the small-constant a-priori bound remains the
      one open analytic item; here it enters as the hypothesis `N/V < total-length`).

  The heart — `exists_gridpoint_Ico` (a long interval catches a grid point) and `exists_long_Ico`
  (of `N` intervals summing to `> N/V`, one exceeds `1/V`) — is elementary and unconditional.
-/
import TournamentH7.LRCTailDiameter

namespace LonelyRunner
namespace TailDiameter

open scoped ENNReal
open Set

/-! ### 1. A long interval catches a grid point

If `Ico a b` is longer than `1/V`, then some grid point `j/V` (`j : ℤ`) lies in it.
Take `j = ⌈a·V⌉`: then `a·V ≤ j < a·V + 1 < b·V`, so `a ≤ j/V < b`. -/

theorem exists_gridpoint_Ico {V : ℕ} (hV : 0 < V) {a b : ℝ}
    (hlen : 1 / (V : ℝ) < b - a) : ∃ j : ℤ, (j : ℝ) / V ∈ Ico a b := by
  have hVR : (0 : ℝ) < V := by exact_mod_cast hV
  refine ⟨⌈a * V⌉, ?_, ?_⟩
  · -- `a ≤ ⌈a·V⌉ / V`
    rw [le_div_iff₀ hVR]
    exact Int.le_ceil (a * V)
  · -- `⌈a·V⌉ / V < b`
    rw [div_lt_iff₀ hVR]
    -- from `1/V < b - a` : `V·a + 1 < V·b`
    have h1 : (1 : ℝ) < V * (b - a) := by
      have := (div_lt_iff₀ hVR).1 hlen        -- `1 < (b - a) * V`
      nlinarith [this]
    have hceil : (⌈a * V⌉ : ℝ) < a * V + 1 := Int.ceil_lt_add_one (a * V)
    nlinarith [hceil, h1]

/-! ### 2. Pigeonhole on interval lengths

Of `N` intervals whose lengths sum to more than `N/V`, at least one exceeds `1/V`. -/

theorem exists_long_Ico {V : ℕ} (hV : 0 < V) {N : ℕ} (a b : Fin N → ℝ)
    (hsum : (N : ℝ) / V < ∑ i, (b i - a i)) : ∃ i, 1 / (V : ℝ) < b i - a i := by
  by_contra h
  push_neg at h                                  -- `∀ i, b i - a i ≤ 1/V`
  have hle : (∑ i, (b i - a i)) ≤ ∑ _i : Fin N, (1 / (V : ℝ)) :=
    Finset.sum_le_sum (fun i _ => h i)
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hle
  -- `∑ ≤ N · (1/V) = N/V`, contradicting `N/V < ∑`
  rw [mul_one_div] at hle
  linarith [hsum, hle]

/-! ### 3. The arc-count pigeonhole: a good period exists

If the good set `Good θ E` contains `N` intervals of total length `> N/V`, then some grid
point `j/V` is good — a *good period*.  (Take `Good θ E ∩ [0,1)` to be the union of its `N`
maximal arcs; then `∑ lengths = ρ*` and the hypothesis is exactly `N < ρ*·V`.) -/

theorem good_period_of_arccount {θ : ℝ} {E : Finset ℤ} {V : ℕ} (hV : 0 < V) {N : ℕ}
    (a b : Fin N → ℝ)
    (hcover : (⋃ i, Ico (a i) (b i)) ⊆ Good θ E)
    (hmeas : (N : ℝ) / V < ∑ i, (b i - a i)) :
    ∃ j : ℤ, ((j : ℝ) / V) ∈ Good θ E := by
  obtain ⟨i, hi⟩ := exists_long_Ico hV a b hmeas
  obtain ⟨j, hj⟩ := exists_gridpoint_Ico hV hi
  exact ⟨j, hcover (mem_iUnion.2 ⟨i, hj⟩)⟩

/-- Contrapositive reading: **no good period ⟹ the good arcs are few** (`ρ*·V ≤ #arcs`).
This is the exact obstruction — realized only by the tight complete-residue AP at `V` prime
(opus-S164), the LRC extremal.  Any speed set beating `#arcs < ρ*·V` has a good period. -/
theorem arccount_le_of_no_good_period {θ : ℝ} {E : Finset ℤ} {V : ℕ} (hV : 0 < V) {N : ℕ}
    (a b : Fin N → ℝ)
    (hcover : (⋃ i, Ico (a i) (b i)) ⊆ Good θ E)
    (hno : ∀ j : ℤ, ((j : ℝ) / V) ∉ Good θ E) :
    (∑ i, (b i - a i)) ≤ (N : ℝ) / V := by
  by_contra h
  push_neg at h
  obtain ⟨j, hj⟩ := good_period_of_arccount hV a b hcover h
  exact hno j hj

end TailDiameter
end LonelyRunner
