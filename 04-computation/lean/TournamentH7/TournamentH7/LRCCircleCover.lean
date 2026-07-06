/-
  TournamentH7.LRCCircleCover — THE CIRCLE-COVERING FLOOR AND THE (A)-WINDOW
  REDUCTION (kind-pasteur-2026-07-06-S20d, HYP-4247 follow-on).

  The torus-split rung (LRCTorusSplit) killed the (A) window for ≤ 6 lifted by
  a DENSITY argument (`2ρl < 1`).  Numerics (lrc_seven_comb_covering_kps_S20b,
  lrc_covering_threshold_kps_S20c) show something stronger: DISTINCT-frequency
  combs of radius `2/25` leave a positive uncovered floor for `l ≤ 9`
  (uncovered ≥ 0.025 even for consecutive `{1..9}`).  CORRECTION (S20f, aligns
  with mac-mini S5's self-correction): at `l ≥ 10` CONSECUTIVE frequencies
  `{1..l}` TILE at `2/25` (uncovered → 0: 0.0078 at l=10, 0.0064 at l=11) — so
  the distinct-frequency floor FAILS for `l ≥ 10`.  My earlier "dead through
  l=14" was a RANDOM-search artifact (random sets avoid the tiling; consecutive
  tiles).  So the fixed-slice covering-impossibility closes the (A) window only
  for `l ≤ 9`; `l = 10, 11` join the census/phase-orbit bucket (the clear point
  must come from the `t`-variation — `torus_forced_rectangle` — not a single
  base-clear slice).

  This file makes the reduction precise:

  * `circle_clear_of_density`: `2ρ·|S| < 1` ⟹ a `ρ`-clear point exists on the
    circle — the density lane, PROVED (via `tooth_visit_count_rho`, evaluated
    at an uncovered GRID point, no measure theory).

  * `CircleClearFloor ρ l`: the NAMED obligation — `l` distinct-frequency combs
    at radius `ρ` always leave a clear point.  `circleClearFloor_of_le6`
    discharges it for `l ≤ 6` (density).  It is TRUE (numerically) for
    `7 ≤ l ≤ 9` (the genuine Newman-shaped distinctness lemma, opus's `φ > 0`
    lane) but FALSE for `l ≥ 10` (consecutive tiles) — so `torus_A_window_empty`
    below closes the distinct-freq (A) window only up to `|L| ≤ 9`.

  * `torus_A_window_empty`: base citation + `CircleClearFloor (2/25) |L|` ⟹ a
    proper coupled 2-torus system with DISTINCT lifted frequencies has a
    `2/25`-clear point, i.e. its value is outside the open gap `(1/13, 2/25)`.
    With `circleClearFloor_of_le6` this is UNCONDITIONAL for `|L| ≤ 6`; the
    DISTINCT-frequency (A) window reduces to the named floor at `7 ≤ l ≤ 9`
    (true, Newman-shaped).  `|L| = 10, 11` do NOT reduce this way (the floor is
    false there — consecutive tiles), and fall to the `t`-variation / census.

  SCOPE (honest, complementary to mac-mini HYP-4292): the `hdistinct`
  hypothesis restricts to lifted runners of DISTINCT θ-frequency — the
  SINGLE-lift-class / k-stratification stratum (k-multiples have distinct
  quotients).  The GENERAL 7-spread residual is MULTI-class: the 12 pair-
  vectors fall into ≥ 3 direction-classes, runners PARALLEL within a class
  share a frequency (repeated frequencies, which CAN measure-cover T²).  That
  multi-class stratum is mac-mini HYP-4292 — census-clean at infimum EXACTLY
  `1/6 > 2/25` (the 5-5-2 minimizer), reduced to one ≥3-class covering lemma.
  Together: distinct-freq `l ≤ 9` (this file, floor ≥ 0.025) + `l = 10,11`
  (t-variation / census, since consecutive tiles) + multi-class (HYP-4292,
  infimum 1/6) = the full (A) residual.

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCTorusSplit

namespace LonelyRunner
namespace CircleCover

open TeethR ClusterGcd ClusterGcdSharp TorusSplit
open scoped Classical

/-! ## the density lane on the circle -/

/-- **CIRCLE DENSITY FLOOR**: `|S|` combs of radius `ρ` (frequencies `r i > 0`,
shifts `s i`) with total density `2ρ|S| < 1` cannot cover the circle — a
`ρ`-clear `θ` exists.  Proved by counting an integer grid, evaluated at an
uncovered grid point (no measure theory). -/
theorem circle_clear_of_density (S : Finset (Fin 12)) (r : Fin 12 → ℤ)
    (s : Fin 12 → ℝ) (ρ : ℝ) (hρ0 : 0 < ρ) (hρh : ρ ≤ 1/2)
    (hr : ∀ i ∈ S, 0 < r i) (hdens : 2 * ρ * S.card < 1) :
    ∃ θ : ℝ, ∀ i ∈ S, ∀ m : ℤ, ρ ≤ |(r i : ℝ) * θ + s i - m| := by
  classical
  rcases S.eq_empty_or_nonempty with hSe | hSne
  · -- empty S: the goal is vacuous
    exact ⟨0, by rw [hSe]; intro i hi; simp at hi⟩
  by_contra hcon
  push_neg at hcon
  -- hcon : ∀ θ, ∃ i ∈ S, ∃ m, |r i θ + s i - m| < ρ  (some comb covers θ)
  have hcover : ∀ θ : ℝ, ∃ i, i ∈ S ∧ ∃ m : ℤ, |(r i : ℝ) * θ + s i - m| < ρ := by
    intro θ
    obtain ⟨i, hiS, m, hm⟩ := hcon θ
    exact ⟨i, hiS, m, hm⟩
  set R : ℤ := ∑ i ∈ S, r i with hR
  -- density gives the gap and the large grid
  have hgap : 0 < 1 - 2 * ρ * (S.card : ℝ) := by linarith
  have hR0 : 0 < R := Finset.sum_pos (fun i hi => hr i hi) hSne
  have hRR : (0 : ℝ) < (R : ℝ) := by exact_mod_cast hR0
  obtain ⟨D, hD⟩ := exists_nat_gt (3 * (R : ℝ) / (1 - 2 * ρ * (S.card : ℝ)))
  have hDpos0 : (0 : ℝ) < 3 * (R : ℝ) / (1 - 2 * ρ * (S.card : ℝ)) :=
    div_pos (by linarith) hgap
  have hD0 : 0 < D := by
    have : (0 : ℝ) < (D : ℝ) := lt_trans hDpos0 hD
    exact_mod_cast this
  -- classify each grid point θ = j/D to a covering comb
  have hchoice : ∀ j : ℕ, ∃ i : Fin 12, j ∈ Finset.range D →
      (i ∈ S ∧ ∃ m : ℤ, |(r i : ℝ) * ((j : ℝ) / D) + s i - m| < ρ) := by
    intro j
    by_cases hj : j ∈ Finset.range D
    · obtain ⟨i, hiS, hm⟩ := hcover ((j : ℝ) / D)
      exact ⟨i, fun _ => ⟨hiS, hm⟩⟩
    · exact ⟨0, fun hc => absurd hc hj⟩
  choose F hFspec using hchoice
  have hFmem : ∀ j ∈ Finset.range D, F j ∈ S := fun j hj => (hFspec j hj).1
  have hsplit : (Finset.range D).card
      = ∑ i ∈ S, ((Finset.range D).filter (fun j => F j = i)).card :=
    Finset.card_eq_sum_card_fiberwise hFmem
  -- per fiber: inside the ρ-visit filter of r i
  have hfiber : ∀ i ∈ S, (((Finset.range D).filter (fun j => F j = i)).card : ℝ)
      ≤ 2 * ρ * D + 3 * ((r i : ℤ) : ℝ) := by
    intro i hi
    have hsub : (Finset.range D).filter (fun j => F j = i)
        ⊆ (Finset.range D).filter
          (fun j : ℕ => ∃ m : ℤ, |s i + (r i : ℝ) * ((j : ℝ) / D) - m| < ρ) := by
      intro j hj
      rw [Finset.mem_filter] at hj ⊢
      obtain ⟨hjr, hjF⟩ := hj
      have hm := (hFspec j hjr).2
      rw [hjF] at hm
      obtain ⟨m, hm⟩ := hm
      refine ⟨hjr, m, ?_⟩
      have harg : s i + (r i : ℝ) * ((j : ℝ) / D) - m
          = (r i : ℝ) * ((j : ℝ) / D) + s i - m := by ring
      rw [harg]; exact hm
    have hcount := tooth_visit_count_rho (r i) (hr i hi) D hD0 (s i) ρ hρ0 hρh
    calc (((Finset.range D).filter (fun j => F j = i)).card : ℝ)
        ≤ (((Finset.range D).filter
          (fun j : ℕ => ∃ m : ℤ, |s i + (r i : ℝ) * ((j : ℝ) / D) - m| < ρ)).card : ℝ) := by
          exact_mod_cast Finset.card_le_card hsub
      _ ≤ 2 * ρ * D + 3 * ((r i : ℤ) : ℝ) := hcount
  -- assemble: D ≤ 2ρ|S|D + 3R
  have hsum : ((D : ℕ) : ℝ) ≤ 2 * ρ * D * S.card + 3 * ((R : ℤ) : ℝ) := by
    have hcards : ((D : ℕ) : ℝ)
        = ∑ i ∈ S, (((Finset.range D).filter (fun j => F j = i)).card : ℝ) := by
      have h0 : ((Finset.range D).card : ℝ)
          = ∑ i ∈ S, (((Finset.range D).filter (fun j => F j = i)).card : ℝ) := by
        exact_mod_cast congrArg (Nat.cast : ℕ → ℝ) hsplit
      rw [← h0, Finset.card_range]
    conv_lhs => rw [hcards]
    have hRcast : ((R : ℤ) : ℝ) = ∑ i ∈ S, ((r i : ℤ) : ℝ) := by
      rw [hR]; push_cast; ring
    calc (∑ i ∈ S, (((Finset.range D).filter (fun j => F j = i)).card : ℝ))
        ≤ ∑ i ∈ S, (2 * ρ * (D : ℝ) + 3 * ((r i : ℤ) : ℝ)) := Finset.sum_le_sum hfiber
      _ = 2 * ρ * D * S.card + 3 * ((R : ℤ) : ℝ) := by
          rw [Finset.sum_add_distrib, Finset.sum_const, ← Finset.mul_sum, hRcast]
          ring
  have hDbig : 3 * (R : ℝ) < (1 - 2 * ρ * (S.card : ℝ)) * D := by
    rw [div_lt_iff₀ hgap] at hD; linarith [hD]
  nlinarith [hsum, hDbig]

/-! ## the named covering-floor obligation -/

/-- **THE CIRCLE CLEAR FLOOR** at radius `ρ` for `l` combs: any `l`
DISTINCT-frequency combs (positive frequencies, arbitrary shifts) leave a
`ρ`-clear point.  For `ρ = 2/25` this is:
  · PROVED for `l ≤ 6` (`circleClearFloor_of_le6`, density);
  · the Newman-shaped distinctness lemma for `7 ≤ l ≤ 11` (numerically
    confirmed: uncovered floor ≥ 0.06). -/
def CircleClearFloor (ρ : ℝ) (l : ℕ) : Prop :=
  ∀ (S : Finset (Fin 12)) (r : Fin 12 → ℤ) (s : Fin 12 → ℝ),
    S.card = l → (∀ i ∈ S, 0 < r i) →
    (∀ i ∈ S, ∀ j ∈ S, r i = r j → i = j) →
    ∃ θ : ℝ, ∀ i ∈ S, ∀ m : ℤ, ρ ≤ |(r i : ℝ) * θ + s i - m|

/-- The floor is UNCONDITIONAL for `l ≤ 6` — the density lane (distinctness
unused here; it becomes essential only at `l ≥ 7`, where `2ρl ≥ 1`). -/
theorem circleClearFloor_of_le6 (l : ℕ) (hl : l ≤ 6) :
    CircleClearFloor (2/25) l := by
  intro S r s hcard hr _
  apply circle_clear_of_density S r s (2/25) (by norm_num) (by norm_num) hr
  rw [hcard]
  have hl6 : (l : ℝ) ≤ 6 := by exact_mod_cast hl
  linarith

/-! ## the (A)-window reduction -/

/-- **THE (A) WINDOW REDUCTION**: a PROPER coupled 2-torus system (base
`w i` on `t`, distinct lifted frequencies `r i > 0` with couplings `a i`,
`L ≠ univ` so the base is nonempty) whose lifted count meets the circle clear
floor has a `2/25`-clear point `(t, θ)` — hence its value is `≥ 2/25`, outside
the open gap `(1/13, 2/25)`.

With `circleClearFloor_of_le6` this is UNCONDITIONAL for `|L| ≤ 6`; the whole
(A) window (mac-mini HYP-4262's coupled-2-torus crux) reduces to the single
named floor at `7 ≤ |L| ≤ 11`. -/
theorem torus_A_window_empty (cite : LRCUpTo13)
    (L : Finset (Fin 12)) (hL : L.Nonempty) (hLne : L ≠ Finset.univ)
    (w r a : Fin 12 → ℤ)
    (hw : ∀ i, i ∉ L → w i ≠ 0) (hr : ∀ i, i ∈ L → 0 < r i)
    (hdistinct : ∀ i ∈ L, ∀ j ∈ L, r i = r j → i = j)
    (floor : CircleClearFloor (2/25) L.card) :
    ∃ t θ : ℝ, (∀ i, i ∉ L → ∀ m : ℤ, 2/25 ≤ |(w i : ℝ) * t - m|) ∧
               (∀ i, i ∈ L → ∀ m : ℤ, 2/25 ≤ |(r i : ℝ) * θ + (a i : ℝ) * t - m|) := by
  classical
  -- the base is nonempty (L ≠ univ), cite it
  set T : Finset (Fin 12) := Finset.univ \ L with hT
  have hTmem : ∀ i, i ∈ T ↔ i ∉ L := by
    intro i; rw [hT, Finset.mem_sdiff]; simp only [Finset.mem_univ, true_and]
  have hTne : T.Nonempty := by
    rw [hT]
    rw [Finset.sdiff_nonempty]
    intro huniv
    exact hLne (Finset.eq_univ_of_forall (fun i => huniv (Finset.mem_univ i)))
  have hTcard : T.card ≤ 12 := by
    calc T.card ≤ Finset.univ.card := Finset.card_le_card (Finset.subset_univ T)
      _ = 12 := by rw [Finset.card_univ, Fintype.card_fin]
  have hT11 : T.card ≤ 11 := by
    have hL1 : 1 ≤ L.card := Finset.card_pos.mpr hL
    have : T.card = 12 - L.card := by
      rw [hT, Finset.card_univ_diff, Fintype.card_fin]
    omega
  obtain ⟨t₀, hmargin⟩ := cite_margin_gen cite
    (fun i => if i ∈ L then 1 else w i) T hTcard
    (fun i hi => by
      show (if i ∈ L then 1 else w i) ≠ 0
      rw [if_neg ((hTmem i).mp hi)]
      exact hw i ((hTmem i).mp hi))
  have hbase : ∀ i, i ∉ L → ∀ m : ℤ, (2:ℝ)/25 ≤ |(w i : ℝ) * t₀ - m| := by
    intro i hi m
    have h := hmargin i ((hTmem i).mpr hi) m
    rw [if_neg hi] at h
    have hcard : ((T.card : ℝ) + 1) ≤ 12 := by
      have : (T.card : ℝ) ≤ 11 := by exact_mod_cast hT11
      linarith
    have h12 : (1:ℝ)/12 ≤ 1 / ((T.card : ℝ) + 1) := by
      apply div_le_div_of_nonneg_left (by norm_num) ?_ hcard
      positivity
    calc (2:ℝ)/25 ≤ 1/12 := by norm_num
      _ ≤ 1 / ((T.card : ℝ) + 1) := h12
      _ ≤ |(w i : ℝ) * t₀ - m| := h
  -- the lifted, at fixed t₀, are distinct-frequency combs with shifts a i · t₀
  obtain ⟨θ₀, hθ₀⟩ := floor L r (fun i => (a i : ℝ) * t₀) rfl hr hdistinct
  refine ⟨t₀, θ₀, hbase, ?_⟩
  intro i hiL m
  exact hθ₀ i hiL m

/-- **UNCONDITIONAL for `|L| ≤ 6`**: a proper coupled 2-torus system with at
most 6 distinct lifted frequencies has a `2/25`-clear point. -/
theorem torus_A_window_empty_le6 (cite : LRCUpTo13)
    (L : Finset (Fin 12)) (hL : L.Nonempty) (hLne : L ≠ Finset.univ) (hL6 : L.card ≤ 6)
    (w r a : Fin 12 → ℤ)
    (hw : ∀ i, i ∉ L → w i ≠ 0) (hr : ∀ i, i ∈ L → 0 < r i)
    (hdistinct : ∀ i ∈ L, ∀ j ∈ L, r i = r j → i = j) :
    ∃ t θ : ℝ, (∀ i, i ∉ L → ∀ m : ℤ, 2/25 ≤ |(w i : ℝ) * t - m|) ∧
               (∀ i, i ∈ L → ∀ m : ℤ, 2/25 ≤ |(r i : ℝ) * θ + (a i : ℝ) * t - m|) :=
  torus_A_window_empty cite L hL hLne w r a hw hr hdistinct
    (circleClearFloor_of_le6 L.card hL6)

#print axioms circle_clear_of_density
#print axioms circleClearFloor_of_le6
#print axioms torus_A_window_empty
#print axioms torus_A_window_empty_le6

end CircleCover
end LonelyRunner
