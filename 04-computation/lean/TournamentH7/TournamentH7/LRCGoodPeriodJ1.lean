/-
  TournamentH7.LRCGoodPeriodJ1 — LEM-010(i): the `j = 1` wraparound good-period lemma
  (mac-mini-2026-07-08-S59).

  Part of closing THM-527-A (the finite-`Vmax` glue of the LRC(14) covering case).  A "good
  period" `j ∈ {1,…,Vmax-1}` is one where the cluster phases `{frac(e_i · j / Vmax)}` leave a
  circular gap `> 1/7` (so the observer runner `Vmax` can sit `> 1/14` from every cluster
  runner, forcing `M(S) ≥ 1/14`).

  This file proves the ELEMENTARY half (LEM-010(i)): if the cluster's co-offsets `E` all lie in
  `[0, spread]` with `spread < 6·Vmax/7`, then at `j = 1` the phases `{e/Vmax}` all lie in
  `[0, spread/Vmax]`, so the open arc `(spread/Vmax, 1)` is empty and has length
  `1 - spread/Vmax > 1/7` — a good period.  (The observer's own free window.)

  The other half (Dirichlet, `Vmax > 3^{k-1}`) and the sharp `j* = O(k)` bound are LEM-010(ii),(iii).
-/
import Mathlib

namespace LRC14

/-- **LEM-010(i) — the `j = 1` wraparound good period.**
If every co-offset `e ∈ E` satisfies `e ≤ spread` and `7 * spread < 6 * Vmax` (the cluster is
compressed below `6/7` of the ruler), then there is a gap of length `> 1/7` open past every
phase `e / Vmax`: concretely, `gapLen = 1 - spread/Vmax` exceeds `1/7`, and every phase
`e/Vmax` sits at or before the start of that gap (`e/Vmax + gapLen ≤ 1`).  A good period exists
at `j = 1`. -/
theorem good_period_j1_wraparound
    (Vmax spread : ℕ) (E : Finset ℕ)
    (hV : 0 < Vmax)
    (hE : ∀ e ∈ E, e ≤ spread)
    (hsmall : 7 * spread < 6 * Vmax) :
    ∃ gapLen : ℚ, (1 : ℚ) / 7 < gapLen ∧
      ∀ e ∈ E, (e : ℚ) / Vmax + gapLen ≤ 1 := by
  have hVQ : (0 : ℚ) < Vmax := by exact_mod_cast hV
  have hsmallQ : (7 : ℚ) * (spread : ℚ) < 6 * (Vmax : ℚ) := by exact_mod_cast hsmall
  refine ⟨1 - (spread : ℚ) / Vmax, ?_, ?_⟩
  · -- 1/7 < 1 - spread/Vmax  ⟺  spread/Vmax < 6/7  ⟺  7·spread < 6·Vmax
    have h : (spread : ℚ) / Vmax < 6 / 7 := by
      rw [div_lt_iff₀ hVQ]; linarith [hsmallQ]
    linarith
  · intro e he
    have heQ : (e : ℚ) ≤ (spread : ℚ) := by exact_mod_cast hE e he
    have hle : (e : ℚ) / Vmax ≤ (spread : ℚ) / Vmax := by gcongr
    linarith

/-- **LEM-010(i), non-strict — the knife-edge wraparound good period (mac-mini-S64).**
The LRC(14) loneliness criterion is `M(S) ≥ 1/14`, i.e. a gap `≥ 1/7` (NON-strict: equality
`gap = 1/7` gives `M = 1/14` exactly, which satisfies the conjecture).  So the good-period
hypothesis is the non-strict `7 * spread ≤ 6 * Vmax` (the closed compression `spread ≤ 6·Vmax/7`),
delivering `gapLen ≥ 1/7`.  This covers the **wraparound-boundary knife-edge** `spread = 6·Vmax/7`
(e.g. the 7-structured covering set `{0,7,10,14,18,20,21,26,28,35,36,37,42}` at `Vmax = 49`, where
`spread = 42 = 6·49/7` and `gapLen = 1 − 6/7 = 1/7` exactly) — the extremal case that the strict
lemma `good_period_j1_wraparound` (`<`, `>`) EXCLUDES.  Verified (`lrc14_nonstrict_knife_edge`):
over 7-structured dissociated `k=13` covering sets on the resonant grid `7 ∣ Vmax`, the loneliness
margin `maxgap·7 − Vmax` is `≥ 0` everywhere (no counterexample) with equality achieved exactly at
this boundary — the genuine `n = 14` knife-edge. -/
theorem good_period_j1_wraparound_nonstrict
    (Vmax spread : ℕ) (E : Finset ℕ)
    (hV : 0 < Vmax)
    (hE : ∀ e ∈ E, e ≤ spread)
    (hsmall : 7 * spread ≤ 6 * Vmax) :
    ∃ gapLen : ℚ, (1 : ℚ) / 7 ≤ gapLen ∧
      ∀ e ∈ E, (e : ℚ) / Vmax + gapLen ≤ 1 := by
  have hVQ : (0 : ℚ) < Vmax := by exact_mod_cast hV
  have hsmallQ : (7 : ℚ) * (spread : ℚ) ≤ 6 * (Vmax : ℚ) := by exact_mod_cast hsmall
  refine ⟨1 - (spread : ℚ) / Vmax, ?_, ?_⟩
  · -- 1/7 ≤ 1 - spread/Vmax  ⟺  spread/Vmax ≤ 6/7  ⟺  7·spread ≤ 6·Vmax
    have h : (spread : ℚ) / Vmax ≤ 6 / 7 := by
      rw [div_le_iff₀ hVQ]; linarith [hsmallQ]
    linarith
  · intro e he
    have heQ : (e : ℚ) ≤ (spread : ℚ) := by exact_mod_cast hE e he
    have hle : (e : ℚ) / Vmax ≤ (spread : ℚ) / Vmax := by gcongr
    linarith

/-- **The good-period core (reusable).** If a finite set `P` of phases all lies in an interval
`[lo, hi]` of length `< 6/7`, then the complementary circular arc `(hi, lo+1)` — of length
`1 − (hi − lo) > 1/7` — is empty: a gap `> 1/7`, i.e. a good period.  This is the shared engine of
LEM-010: `j=1` (phases in `[0, spread/Vmax]`), the Dirichlet pigeonhole (phases in a `2/3`-arc),
and the AP lemma (a `k`-term AP of small step spans `< 6/7`).  `gapLen := 1 − (hi − lo)` witnesses
the gap; every phase `p` satisfies `p + gapLen ≤ 1 + lo` (it sits at or before the arc's start). -/
theorem good_gap_of_phases_in_interval
    (P : Finset ℚ) (lo hi : ℚ)
    (hbound : ∀ p ∈ P, lo ≤ p ∧ p ≤ hi)
    (hlen : hi - lo < 6 / 7) :
    ∃ gapLen : ℚ, (1 : ℚ) / 7 < gapLen ∧ ∀ p ∈ P, p + gapLen ≤ 1 + lo := by
  refine ⟨1 - (hi - lo), by linarith, ?_⟩
  intro p hp
  have := (hbound p hp).2
  linarith

/-- **The capstone reduction (opus-S165).** Let `W j := ` the uncovered measure at period `j`
(`W j > 0 ⟺ j` is a good period), which is nonnegative. Then the partial sum `S_N = Σ_{j=1}^N W j`
is positive **iff** some `j ∈ {1,…,N}` is a good period.  This is why the whole finite-`Vmax` glue
reduces to the single inequality `S_N > 0` (equivalently `r_N < 1`): a sum of nonnegatives is
positive exactly when one summand is. -/
theorem goodPeriod_iff_partialSum_pos
    (W : ℕ → ℚ) (N : ℕ) (hW : ∀ j, 0 ≤ W j) :
    0 < ∑ j ∈ Finset.Icc 1 N, W j ↔ ∃ j ∈ Finset.Icc 1 N, 0 < W j := by
  constructor
  · intro hsum
    by_contra hcon
    push_neg at hcon
    have : ∑ j ∈ Finset.Icc 1 N, W j ≤ 0 := Finset.sum_nonpos (fun j hj => hcon j hj)
    linarith
  · rintro ⟨j, hj, hpos⟩
    exact Finset.sum_pos' (fun i _ => hW i) ⟨j, hj, hpos⟩

/-- **The gap-split pigeonhole — the combinatorial core of LEM-012 (near-AP branch, klein-S196).**
After the sub-AP is Dirichlet-clustered into a small arc, its complement is one gap of length
`> (m+1)/7`, and the `m = k − L` stray points cut that gap into `m+1` sub-gaps `g_0,…,g_m ≥ 0`
summing to the gap's length `> (m+1)·(1/7)`.  Then the LARGEST sub-gap exceeds `1/7` — a good
period.  (A sum of `m+1` nonnegatives exceeding `(m+1)/7` must have a term `> 1/7`.)  This is why
`L ≥ k−5` (`m ≤ 5`) is exactly the range where the pigeonhole bites. -/
theorem gap_split_pigeonhole (m : ℕ) (g : Fin (m + 1) → ℝ)
    (hsum : ((m : ℝ) + 1) * (1 / 7) < ∑ i, g i) :
    ∃ i, (1 : ℝ) / 7 < g i := by
  by_contra h
  push_neg at h
  have hle : ∑ i, g i ≤ ((m : ℝ) + 1) * (1 / 7) := by
    calc ∑ i, g i ≤ ∑ _i : Fin (m + 1), (1 / 7 : ℝ) := Finset.sum_le_sum (fun i _ => h i)
      _ = ((m : ℝ) + 1) * (1 / 7) := by
          rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin]; ring
  linarith

/-- **The AP clustering good period — the pure-AP case of LEM-012 (exact-AP, mac-mini-S59).**
For an arithmetic progression of `k` phases `{0, α, 2α, …, (k−1)α}` with a small positive step
`α ≥ 0`, if the span `(k−1)·α < 6/7` then all `k` phases lie in `[0, (k−1)α] ⊆ [0, 6/7)`, so the
complement arc — of length `1 − (k−1)α > 1/7` — is an empty gap: a good period.  This is the
`m = 0` case of `gap_split_pigeonhole` and the Dirichlet-consuming step: after Dirichlet gives a
dilation `j ≤ ⌈7(k−1)/6⌉` with `‖jd/Vmax‖ = α < 6/(7(k−1))`, this lemma delivers the gap.
(Mathlib's `Real.exists_nat_abs_mul_sub_round_le` supplies the Dirichlet `α`.) -/
theorem ap_clustered_good_period (k : ℕ) (α : ℚ) (hα0 : 0 ≤ α)
    (hspan : ((k : ℚ) - 1) * α < 6 / 7) :
    ∃ gapLen : ℚ, (1 : ℚ) / 7 < gapLen ∧ ∀ i : ℕ, i < k → (i : ℚ) * α + gapLen ≤ 1 := by
  refine ⟨1 - ((k : ℚ) - 1) * α, by linarith, ?_⟩
  intro i hi
  have hik : (i : ℚ) ≤ (k : ℚ) - 1 := by
    have : (i : ℚ) + 1 ≤ (k : ℚ) := by exact_mod_cast hi
    linarith
  have : (i : ℚ) * α ≤ ((k : ℚ) - 1) * α := mul_le_mul_of_nonneg_right hik hα0
  linarith

end LRC14
