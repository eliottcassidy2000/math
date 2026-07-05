/-
  TournamentH7.LRCGridAttainment — THM-592 GRID ATTAINMENT IN LEAN
  (kind-pasteur-2026-07-05-S4, HYP-4108).

  THE THEOREM: the loneliness margin `L(t) = min_i dist(vᵢt, ℤ)` attains its
  maximum at a rational point `m / (|vᵢ| + |vⱼ|)` — the equioscillation /
  merge-grid structure of THM-592 (mac-mini S93), machine-checked.

  THE PROOF (no derivatives, no piecewise-linear theory): take a global
  maximizer `t*` (mac-mini's `exists_max_margin` + periodicity).  Write each
  runner's offset `eᵢ = vᵢt* − round(vᵢt*)`, `|eᵢ| = distZ(vᵢt*)`.
  * If `M = 1/2`: every runner sits exactly at a half-integer, so `t*` is
    already on the grid with `i = j` (`2|vᵢ|·t* ∈ ℤ`).
  * If `0 < M < 1/2`: if every binding runner ascends in direction `d`
    (`eᵢ·vᵢd > 0`), the explicit push `t* + d·ε` with
    `ε = ½·minᵢ(binding ? (1−2M)/(2|vᵢ|) : (distZᵢ−M)/(2|vᵢ|))`
    STRICTLY increases the margin (binding runners climb by `|vᵢ|ε` and stay
    under the `1/2` roof; non-binding runners cannot fall to `M` by Lipschitz)
    — contradicting maximality (`no_uniform_push`).  So both directions have a
    binding runner: one descending (`eᵢvᵢ < 0`), one ascending (`eⱼvⱼ > 0`).
    Then `|vᵢ|t* = nᵢ − M` and `|vⱼ|t* = nⱼ + M` for integers `nᵢ, nⱼ`, and
    summing CANCELS `M`:  `(|vᵢ| + |vⱼ|)·t* = nᵢ + nⱼ ∈ ℤ`.

  CONSEQUENCES:
  * `grid_margin_domination`: margin `β` anywhere ⟹ margin `β` at a grid point
    `m/(|vᵢ|+|vⱼ|)` — every sweep's exact-M search space is now formally the
    merge grid (and the S3 non-reduced-fraction trap is a theorem, not lore).
  * `loose_branch_cert_2B`: the dichotomy's loose branch yields an integer
    certificate with modulus `s ≤ 2B` and ZERO slack — superseding S3's
    `s ≤ B/(2(β'−β)) + 1 ≈ 176B`.  The loose side of `TightLooseDichotomy`
    is a decidable search over `s ≤ 2B`, kernel-checkable via HYP-4102's atom.

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCWitnessAttainment

namespace LonelyRunner
namespace GridAttainment

open TournamentH7.LRCWitness

/-! ## distZ toolbox (on top of mac-mini's `distZ`) -/

/-- `distZ` is below the distance to every integer. -/
lemma distZ_le_abs_sub (y : ℝ) (m : ℤ) : distZ y ≤ |y - m| :=
  (le_distZ_iff (distZ y) y).1 le_rfl m

/-- The round point is the nearest integer. -/
lemma abs_sub_round_le_abs_sub (y : ℝ) (m : ℤ) : |y - round y| ≤ |y - m| := by
  rcases eq_or_ne m (round y) with rfl | hne
  · exact le_refl _
  · have h1 : |y - round y| ≤ 1 / 2 := abs_sub_round y
    have h2 : (1 : ℝ) ≤ |((round y : ℤ) : ℝ) - (m : ℝ)| := by
      rw [← Int.cast_sub, ← Int.cast_abs]
      exact_mod_cast Int.one_le_abs (sub_ne_zero.mpr (Ne.symm hne))
    have htri : |((round y : ℤ) : ℝ) - (m : ℝ)| ≤ |((round y : ℤ) : ℝ) - y| + |y - m| := by
      calc |((round y : ℤ) : ℝ) - (m : ℝ)|
          = |(((round y : ℤ) : ℝ) - y) + (y - m)| := by
            rw [show ((round y : ℤ) : ℝ) - (m : ℝ)
              = (((round y : ℤ) : ℝ) - y) + (y - m) by ring]
        _ ≤ |((round y : ℤ) : ℝ) - y| + |y - m| := abs_add_le _ _
    have hcomm : |((round y : ℤ) : ℝ) - y| = |y - round y| := abs_sub_comm _ _
    linarith

/-- `distZ` is exactly the distance to the round point. -/
lemma distZ_eq_round (y : ℝ) : distZ y = |y - round y| := by
  apply le_antisymm (distZ_le_abs_sub y (round y))
  exact (le_distZ_iff _ y).2 (abs_sub_round_le_abs_sub y)

/-- `distZ ≤ 1/2`. -/
lemma distZ_le_half (y : ℝ) : distZ y ≤ 1 / 2 := by
  rw [distZ_eq_round]; exact abs_sub_round y

/-- Round uniqueness: strictly inside a half-cell the round point is forced. -/
lemma round_eq_of_abs_sub_lt_half (y : ℝ) (m : ℤ) (h : |y - m| < 1 / 2) :
    round y = m := by
  have h1 := abs_lt.mp h
  rw [round_eq]
  apply Int.floor_eq_iff.mpr
  constructor
  · linarith [h1.1]
  · linarith [h1.2]

/-- **The anchor bound**: distance to ℤ is at least `min(|y−n|, 1−|y−n|)` for ANY
integer anchor `n`. -/
lemma anchor_le_distZ (y : ℝ) (n : ℤ) :
    min (|y - n|) (1 - |y - n|) ≤ distZ y := by
  apply (le_distZ_iff _ y).2
  intro m
  rcases eq_or_ne m n with rfl | hne
  · exact min_le_left _ _
  · refine le_trans (min_le_right _ _) ?_
    have h2 : (1 : ℝ) ≤ |((n : ℤ) : ℝ) - (m : ℝ)| := by
      rw [← Int.cast_sub, ← Int.cast_abs]
      exact_mod_cast Int.one_le_abs (sub_ne_zero.mpr (by exact_mod_cast (Ne.symm hne)))
    have htri : |((n : ℤ) : ℝ) - (m : ℝ)| ≤ |((n : ℤ) : ℝ) - y| + |y - m| := by
      calc |((n : ℤ) : ℝ) - (m : ℝ)| = |(((n : ℤ) : ℝ) - y) + (y - m)| := by
            rw [show ((n : ℤ) : ℝ) - (m : ℝ) = (((n : ℤ) : ℝ) - y) + (y - m) by ring]
        _ ≤ |((n : ℤ) : ℝ) - y| + |y - m| := abs_add_le _ _
    have hcomm : |((n : ℤ) : ℝ) - y| = |y - n| := abs_sub_comm _ _
    linarith

/-- One-sided Lipschitz: `distZ x − |x − y| ≤ distZ y`. -/
lemma distZ_lipschitz_lower (x y : ℝ) : distZ x - |x - y| ≤ distZ y := by
  apply (le_distZ_iff _ y).2
  intro m
  have h1 : distZ x ≤ |x - m| := distZ_le_abs_sub x m
  have htri : |x - m| ≤ |x - y| + |y - m| := by
    calc |x - m|
        = |(x - y) + (y - m)| := by
          rw [show x - (m : ℝ) = (x - y) + (y - m) by ring]
      _ ≤ |x - y| + |y - m| := abs_add_le _ _
  linarith

variable {k : ℕ} [Nonempty (Fin k)]

/-- The margin is below each runner's distance. -/
lemma margin_le_distZ (v : Fin k → ℤ) (t : ℝ) (i : Fin k) :
    margin v t ≤ distZ ((v i : ℝ) * t) :=
  Finset.inf'_le _ (Finset.mem_univ i)

/-- Integer-shift periodicity of the margin. -/
lemma margin_add_int (v : Fin k → ℤ) (t : ℝ) (j : ℤ) :
    margin v (t + (j : ℝ)) = margin v t := by
  unfold margin
  have hfun : (fun i => distZ ((v i : ℝ) * (t + (j : ℝ))))
      = (fun i => distZ ((v i : ℝ) * t)) := by
    funext i
    rw [show (v i : ℝ) * (t + (j : ℝ)) = (v i : ℝ) * t + ((v i * j : ℤ) : ℝ) by
      push_cast; ring, distZ_add_int]
  rw [hfun]

/-- The margin attains a GLOBAL maximum (Icc max + periodicity). -/
theorem exists_global_max_margin (v : Fin k → ℤ) :
    ∃ tstar : ℝ, ∀ t : ℝ, margin v t ≤ margin v tstar := by
  obtain ⟨t0, _, hmax⟩ := exists_max_margin v
  refine ⟨t0, fun t => ?_⟩
  have hmem : t - (⌊t⌋ : ℝ) ∈ Set.Icc (0 : ℝ) 1 := by
    constructor
    · linarith [Int.floor_le t]
    · linarith [Int.lt_floor_add_one t]
  have h1 : margin v (t - (⌊t⌋ : ℝ)) ≤ margin v t0 := hmax _ hmem
  have h2 : margin v (t - (⌊t⌋ : ℝ)) = margin v t := by
    rw [show t - (⌊t⌋ : ℝ) = t + ((-⌊t⌋ : ℤ) : ℝ) by push_cast; ring, margin_add_int]
  linarith

/-! ## The push lemma (the ε core) -/

/-- **No uniform push at a maximizer**: at a global maximizer with `0 < M < 1/2`,
the binding runners cannot ALL be ascending in a fixed direction `d = ±1` —
otherwise the explicit push `t* + d·ε` strictly increases the margin. -/
theorem no_uniform_push (v : Fin k → ℤ) (hv : ∀ i, v i ≠ 0) (tstar : ℝ)
    (hmax : ∀ t, margin v t ≤ margin v tstar)
    (hM0 : 0 < margin v tstar) (hM2 : margin v tstar < 1 / 2)
    (d : ℝ) (hd : d = 1 ∨ d = -1)
    (hpush : ∀ i, distZ ((v i : ℝ) * tstar) = margin v tstar →
      0 < ((v i : ℝ) * tstar - round ((v i : ℝ) * tstar)) * ((v i : ℝ) * d)) :
    False := by
  classical
  set M : ℝ := margin v tstar with hMdef
  have habsd : |d| = 1 := by rcases hd with rfl | rfl <;> norm_num
  -- speed magnitudes ≥ 1
  have hv1 : ∀ i, (1 : ℝ) ≤ |(v i : ℝ)| := by
    intro i
    rw [← Int.cast_abs]
    exact_mod_cast Int.one_le_abs (hv i)
  have hv0 : ∀ i, (0 : ℝ) < |(v i : ℝ)| := fun i => lt_of_lt_of_le one_pos (hv1 i)
  -- the per-runner budget
  set c : Fin k → ℝ := fun i =>
    if M < distZ ((v i : ℝ) * tstar)
    then (distZ ((v i : ℝ) * tstar) - M) / (2 * |(v i : ℝ)|)
    else (1 - 2 * M) / (2 * |(v i : ℝ)|) with hcdef
  set ε : ℝ := Finset.univ.inf' Finset.univ_nonempty c with hεdef
  have hε0 : 0 < ε := by
    rw [hεdef, Finset.lt_inf'_iff]
    intro i _
    rw [hcdef]
    dsimp only
    split_ifs with hib
    · exact div_pos (by linarith) (by linarith [hv0 i])
    · exact div_pos (by linarith) (by linarith [hv0 i])
  -- the pushed margin strictly exceeds M
  have hpushed : M < margin v (tstar + d * ε) := by
    unfold margin
    rw [Finset.lt_inf'_iff]
    intro i _
    have hεi : ε ≤ c i := by
      rw [hεdef]
      exact Finset.inf'_le _ (Finset.mem_univ i)
    set y : ℝ := (v i : ℝ) * tstar with hydef
    set δ : ℝ := (v i : ℝ) * d * ε with hδdef
    have hyδ : (v i : ℝ) * (tstar + d * ε) = y + δ := by rw [hydef, hδdef]; ring
    rw [hyδ]
    have habsδ : |δ| = |(v i : ℝ)| * ε := by
      rw [hδdef, abs_mul, abs_mul, habsd, mul_one, abs_of_pos hε0]
    have hMle : M ≤ distZ y := margin_le_distZ v tstar i
    by_cases hib : M < distZ y
    · -- non-binding: Lipschitz keeps it above M
      have hci : c i = (distZ y - M) / (2 * |(v i : ℝ)|) := by
        rw [hcdef]; exact if_pos hib
      have hδle : |δ| ≤ (distZ y - M) / 2 := by
        rw [habsδ]
        calc |(v i : ℝ)| * ε ≤ |(v i : ℝ)| * ((distZ y - M) / (2 * |(v i : ℝ)|)) := by
              apply mul_le_mul_of_nonneg_left _ (le_of_lt (hv0 i))
              rw [← hci]; exact hεi
          _ = (distZ y - M) / 2 := by
              have hvne : |(v i : ℝ)| ≠ 0 := ne_of_gt (hv0 i)
              field_simp
      have hlip : distZ y - |y - (y + δ)| ≤ distZ (y + δ) := distZ_lipschitz_lower y (y + δ)
      have habs2 : |y - (y + δ)| = |δ| := by
        rw [show y - (y + δ) = -δ by ring, abs_neg]
      rw [habs2] at hlip
      linarith
    · -- binding: climb along the sign, stay under the roof
      have hb : distZ y = M := le_antisymm (not_lt.mp hib) hMle
      set e : ℝ := y - round y with hedef
      have habse : |e| = M := by rw [hedef, ← distZ_eq_round]; exact hb
      have hsign : 0 < e * ((v i : ℝ) * d) := hpush i hb
      have hci : c i = (1 - 2 * M) / (2 * |(v i : ℝ)|) := by
        rw [hcdef]; exact if_neg hib
      have hδle : |δ| ≤ (1 - 2 * M) / 2 := by
        rw [habsδ]
        calc |(v i : ℝ)| * ε ≤ |(v i : ℝ)| * ((1 - 2 * M) / (2 * |(v i : ℝ)|)) := by
              apply mul_le_mul_of_nonneg_left _ (le_of_lt (hv0 i))
              rw [← hci]; exact hεi
          _ = (1 - 2 * M) / 2 := by
              have hvne : |(v i : ℝ)| ≠ 0 := ne_of_gt (hv0 i)
              field_simp
      have hδ0 : 0 < |δ| := by
        rw [habsδ]
        exact mul_pos (hv0 i) hε0
      -- same-sign addition: |e + δ| = M + |δ|
      have habs_sum : |e + δ| = M + |δ| := by
        rcases lt_trichotomy ((v i : ℝ) * d) 0 with hw | hw | hw
        · have he : e < 0 := by nlinarith
          have hδneg : δ < 0 := by
            rw [hδdef]
            exact mul_neg_of_neg_of_pos hw hε0
          rw [abs_of_neg (by linarith : e + δ < 0), abs_of_neg hδneg, ← habse,
            abs_of_neg he]
          ring
        · rw [hw, mul_zero] at hsign
          exact absurd hsign (lt_irrefl 0)
        · have he : 0 < e := by nlinarith
          have hδpos : 0 < δ := by
            rw [hδdef]
            exact mul_pos hw hε0
          rw [abs_of_pos (by linarith : 0 < e + δ), abs_of_pos hδpos, ← habse,
            abs_of_pos he]
      -- anchor at the old round point
      have hanchor := anchor_le_distZ (y + δ) (round y)
      have heq : y + δ - (round y : ℝ) = e + δ := by rw [hedef]; ring
      rw [heq, habs_sum] at hanchor
      have hmin : M < min (M + |δ|) (1 - (M + |δ|)) := by
        apply lt_min
        · linarith
        · linarith
      exact lt_of_lt_of_le hmin hanchor
  have := hmax (tstar + d * ε)
  linarith

/-! ## The main theorem -/

/-- **THM-592 GRID ATTAINMENT.**  A global maximizer of the margin with positive
value sits on the merge grid: `(|vᵢ| + |vⱼ|) · t* ∈ ℤ` for some runners `i, j`
(possibly equal). -/
theorem maximizer_on_grid (v : Fin k → ℤ) (hv : ∀ i, v i ≠ 0) (tstar : ℝ)
    (hmax : ∀ t, margin v t ≤ margin v tstar) (hM0 : 0 < margin v tstar) :
    ∃ (i j : Fin k) (m : ℤ), ((|v i| + |v j| : ℤ) : ℝ) * tstar = (m : ℝ) := by
  set M : ℝ := margin v tstar with hMdef
  rcases lt_or_ge M (1 / 2) with hM2 | hM2
  · -- 0 < M < 1/2: two binding runners with opposite drifts
    have hdesc : ∃ i, distZ ((v i : ℝ) * tstar) = M ∧
        ((v i : ℝ) * tstar - round ((v i : ℝ) * tstar)) * (v i : ℝ) < 0 := by
      by_contra hcon
      push_neg at hcon
      apply no_uniform_push v hv tstar hmax hM0 hM2 1 (Or.inl rfl)
      intro i hbi
      have h1 := hcon i hbi
      have he0 : ((v i : ℝ) * tstar - round ((v i : ℝ) * tstar)) ≠ 0 := by
        intro h0
        have : distZ ((v i : ℝ) * tstar) = 0 := by
          rw [distZ_eq_round, h0, abs_zero]
        rw [hbi] at this
        linarith
      have hv0 : ((v i : ℤ) : ℝ) ≠ 0 := by exact_mod_cast hv i
      have hne : ((v i : ℝ) * tstar - round ((v i : ℝ) * tstar)) * (v i : ℝ) ≠ 0 :=
        mul_ne_zero he0 hv0
      rw [mul_one]
      exact lt_of_le_of_ne h1 (Ne.symm hne)
    have hasc : ∃ j, distZ ((v j : ℝ) * tstar) = M ∧
        0 < ((v j : ℝ) * tstar - round ((v j : ℝ) * tstar)) * (v j : ℝ) := by
      by_contra hcon
      push_neg at hcon
      apply no_uniform_push v hv tstar hmax hM0 hM2 (-1) (Or.inr rfl)
      intro j hbj
      have h1 := hcon j hbj
      have he0 : ((v j : ℝ) * tstar - round ((v j : ℝ) * tstar)) ≠ 0 := by
        intro h0
        have : distZ ((v j : ℝ) * tstar) = 0 := by
          rw [distZ_eq_round, h0, abs_zero]
        rw [hbj] at this
        linarith
      have hv0 : ((v j : ℤ) : ℝ) ≠ 0 := by exact_mod_cast hv j
      have hne : ((v j : ℝ) * tstar - round ((v j : ℝ) * tstar)) * (v j : ℝ) ≠ 0 :=
        mul_ne_zero he0 hv0
      have hlt : ((v j : ℝ) * tstar - round ((v j : ℝ) * tstar)) * (v j : ℝ) < 0 :=
        lt_of_le_of_ne h1 hne
      nlinarith [hlt]
    obtain ⟨i, hbi, hdi⟩ := hdesc
    obtain ⟨j, hbj, hdj⟩ := hasc
    -- descending runner: |vᵢ|·t* = nᵢ − M
    have hi : ∃ n : ℤ, ((|v i| : ℤ) : ℝ) * tstar = (n : ℝ) - M := by
      set e : ℝ := (v i : ℝ) * tstar - round ((v i : ℝ) * tstar) with hedef
      have habse : |e| = M := by rw [hedef, ← distZ_eq_round]; exact hbi
      rcases lt_trichotomy (v i) 0 with hvi | hvi | hvi
      · -- v i < 0: e > 0, e = M, |vᵢ|t* = −vᵢt* = −(round + M) = (−round) − M
        have hviR : ((v i : ℤ) : ℝ) < 0 := by exact_mod_cast hvi
        have he : 0 < e := by nlinarith
        have heM : e = M := by rw [← habse, abs_of_pos he]
        refine ⟨-round ((v i : ℝ) * tstar), ?_⟩
        have habs : ((|v i| : ℤ) : ℝ) = -((v i : ℤ) : ℝ) := by
          rw [Int.cast_abs]
          exact abs_of_neg hviR
        rw [habs]
        push_cast
        nlinarith [heM, hedef]
      · exact absurd hvi (hv i)
      · -- v i > 0: e < 0, e = −M, |vᵢ|t* = round − M
        have hviR : (0 : ℝ) < ((v i : ℤ) : ℝ) := by exact_mod_cast hvi
        have he : e < 0 := by nlinarith
        have heM : e = -M := by
          have := abs_of_neg he
          rw [habse] at this
          linarith
        refine ⟨round ((v i : ℝ) * tstar), ?_⟩
        have habs : ((|v i| : ℤ) : ℝ) = ((v i : ℤ) : ℝ) := by
          rw [Int.cast_abs]
          exact abs_of_pos hviR
        rw [habs]
        nlinarith [heM, hedef]
    -- ascending runner: |vⱼ|·t* = nⱼ + M
    have hj : ∃ n : ℤ, ((|v j| : ℤ) : ℝ) * tstar = (n : ℝ) + M := by
      set e : ℝ := (v j : ℝ) * tstar - round ((v j : ℝ) * tstar) with hedef
      have habse : |e| = M := by rw [hedef, ← distZ_eq_round]; exact hbj
      rcases lt_trichotomy (v j) 0 with hvj | hvj | hvj
      · -- v j < 0: e < 0, e = −M, |vⱼ|t* = −vⱼt* = −(round − M) = (−round) + M
        have hvjR : ((v j : ℤ) : ℝ) < 0 := by exact_mod_cast hvj
        have he : e < 0 := by nlinarith
        have heM : e = -M := by
          have := abs_of_neg he
          rw [habse] at this
          linarith
        refine ⟨-round ((v j : ℝ) * tstar), ?_⟩
        have habs : ((|v j| : ℤ) : ℝ) = -((v j : ℤ) : ℝ) := by
          rw [Int.cast_abs]
          exact abs_of_neg hvjR
        rw [habs]
        push_cast
        nlinarith [heM, hedef]
      · exact absurd hvj (hv j)
      · -- v j > 0: e > 0, e = M, |vⱼ|t* = round + M
        have hvjR : (0 : ℝ) < ((v j : ℤ) : ℝ) := by exact_mod_cast hvj
        have he : 0 < e := by nlinarith
        have heM : e = M := by rw [← habse, abs_of_pos he]
        refine ⟨round ((v j : ℝ) * tstar), ?_⟩
        have habs : ((|v j| : ℤ) : ℝ) = ((v j : ℤ) : ℝ) := by
          rw [Int.cast_abs]
          exact abs_of_pos hvjR
        rw [habs]
        nlinarith [heM, hedef]
    obtain ⟨ni, hni⟩ := hi
    obtain ⟨nj, hnj⟩ := hj
    refine ⟨i, j, ni + nj, ?_⟩
    push_cast
    push_cast at hni hnj
    linarith
  · -- M = 1/2: every runner at a half-integer; i = j grid form
    have hM : M = 1 / 2 := by
      have hle : M ≤ 1 / 2 := by
        have := margin_le_distZ v tstar (Classical.arbitrary (Fin k))
        have h2 := distZ_le_half ((v (Classical.arbitrary (Fin k)) : ℝ) * tstar)
        linarith
      linarith
    set i0 : Fin k := Classical.arbitrary (Fin k) with hi0
    have hbd : distZ ((v i0 : ℝ) * tstar) = 1 / 2 := by
      have h1 : M ≤ distZ ((v i0 : ℝ) * tstar) := margin_le_distZ v tstar i0
      have h2 := distZ_le_half ((v i0 : ℝ) * tstar)
      linarith
    have habs : |(v i0 : ℝ) * tstar - round ((v i0 : ℝ) * tstar)| = 1 / 2 := by
      rw [← distZ_eq_round]; exact hbd
    rcases abs_eq (by norm_num : (0:ℝ) ≤ 1/2) |>.mp habs with hpm | hpm
    · -- vᵢ₀t* = round + 1/2 ⟹ 2vᵢ₀t* = 2·round + 1
      rcases lt_trichotomy (v i0) 0 with hvi | hvi | hvi
      · refine ⟨i0, i0, -(2 * round ((v i0 : ℝ) * tstar) + 1), ?_⟩
        have habsc : ((|v i0| : ℤ) : ℝ) = -((v i0 : ℤ) : ℝ) := by
          rw [Int.cast_abs]
          exact abs_of_neg (by exact_mod_cast hvi)
        push_cast [habsc]
        nlinarith [hpm]
      · exact absurd hvi (hv i0)
      · refine ⟨i0, i0, 2 * round ((v i0 : ℝ) * tstar) + 1, ?_⟩
        have habsc : ((|v i0| : ℤ) : ℝ) = ((v i0 : ℤ) : ℝ) := by
          rw [Int.cast_abs]
          exact abs_of_pos (by exact_mod_cast hvi)
        push_cast [habsc]
        nlinarith [hpm]
    · -- vᵢ₀t* = round − 1/2 ⟹ 2vᵢ₀t* = 2·round − 1
      rcases lt_trichotomy (v i0) 0 with hvi | hvi | hvi
      · refine ⟨i0, i0, -(2 * round ((v i0 : ℝ) * tstar) - 1), ?_⟩
        have habsc : ((|v i0| : ℤ) : ℝ) = -((v i0 : ℤ) : ℝ) := by
          rw [Int.cast_abs]
          exact abs_of_neg (by exact_mod_cast hvi)
        push_cast [habsc]
        nlinarith [hpm]
      · exact absurd hvi (hv i0)
      · refine ⟨i0, i0, 2 * round ((v i0 : ℝ) * tstar) - 1, ?_⟩
        have habsc : ((|v i0| : ℤ) : ℝ) = ((v i0 : ℤ) : ℝ) := by
          rw [Int.cast_abs]
          exact abs_of_pos (by exact_mod_cast hvi)
        push_cast [habsc]
        nlinarith [hpm]

/-! ## Domination and the sharp certificate -/

/-- **Grid margin domination**: margin `β > 0` anywhere ⟹ margin `β` at a grid
point `m/(|vᵢ|+|vⱼ|)`.  The exact-M search space IS the merge grid. -/
theorem grid_margin_domination (v : Fin k → ℤ) (hv : ∀ i, v i ≠ 0) (β : ℝ)
    (hβ : 0 < β) (t : ℝ) (hm : ∀ i, ∀ m : ℤ, β ≤ |(v i : ℝ) * t - m|) :
    ∃ (i j : Fin k) (m : ℤ), ∀ i', ∀ m' : ℤ,
      β ≤ |(v i' : ℝ) * ((m : ℝ) / ((|v i| + |v j| : ℤ) : ℝ)) - m'| := by
  obtain ⟨tstar, hmax⟩ := exists_global_max_margin v
  have hβm : β ≤ margin v tstar :=
    le_trans ((le_margin_iff v β t).2 hm) (hmax t)
  have hM0 : 0 < margin v tstar := lt_of_lt_of_le hβ hβm
  obtain ⟨i, j, m, hgrid⟩ := maximizer_on_grid v hv tstar hmax hM0
  refine ⟨i, j, m, ?_⟩
  have hs0 : (0 : ℤ) < |v i| + |v j| := by
    have h1 : (0 : ℤ) < |v i| := abs_pos.mpr (hv i)
    have h2 : (0 : ℤ) < |v j| := abs_pos.mpr (hv j)
    omega
  have hs0R : ((|v i| + |v j| : ℤ) : ℝ) ≠ 0 := by
    have : (0 : ℝ) < ((|v i| + |v j| : ℤ) : ℝ) := by exact_mod_cast hs0
    exact ne_of_gt this
  have ht : (m : ℝ) / ((|v i| + |v j| : ℤ) : ℝ) = tstar := by
    rw [← hgrid]
    field_simp
  rw [ht]
  exact (le_margin_iff v β tstar).1 hβm

/-- Margin at a rational point converts to the atom's mod conditions
(the S3 conversion, factored). -/
theorem mod_of_margin_at_rat {ι : Type*} (v : ι → ℤ) (β : ℝ) (kk s : ℤ)
    (hs : 0 < s)
    (h : ∀ i, ∀ m' : ℤ, β ≤ |(v i : ℝ) * ((kk : ℝ) / (s : ℝ)) - m'|) :
    ∀ i, ⌈β * s⌉ ≤ (v i * kk) % s ∧ (v i * kk) % s ≤ s - ⌈β * s⌉ := by
  intro i
  have hsR : (0 : ℝ) < (s : ℝ) := by exact_mod_cast hs
  have hint : ∀ m : ℤ, β * s ≤ ((|v i * kk - m * s| : ℤ) : ℝ) := by
    intro m
    have h1 := h i m
    have heq : (v i : ℝ) * ((kk : ℝ) / s) - m = ((v i * kk - m * s : ℤ) : ℝ) / s := by
      push_cast
      field_simp
    rw [heq, abs_div, abs_of_pos hsR, le_div_iff₀ hsR, ← Int.cast_abs] at h1
    exact h1
  set r : ℤ := (v i * kk) % s with hrdef
  have hr0 : 0 ≤ r := Int.emod_nonneg _ (ne_of_gt hs)
  have hrs : r < s := Int.emod_lt_of_pos _ hs
  have hX0 : v i * kk - ((v i * kk) / s) * s = r := by
    rw [hrdef, Int.emod_def]; ring
  have hX1 : v i * kk - ((v i * kk) / s + 1) * s = r - s := by
    rw [hrdef, Int.emod_def]; ring
  constructor
  · have h0 := hint ((v i * kk) / s)
    rw [hX0, abs_of_nonneg hr0] at h0
    exact Int.ceil_le.mpr h0
  · have h1 := hint ((v i * kk) / s + 1)
    rw [hX1, abs_of_nonpos (by omega), neg_sub] at h1
    have : ⌈β * s⌉ ≤ s - r := Int.ceil_le.mpr (by exact_mod_cast h1)
    omega

/-- **THE SHARP CERTIFICATE (modulus ≤ 2B, zero slack).**  The dichotomy's loose
branch at margin `β` yields an integer certificate with modulus a SUM OF TWO BASE
SPEEDS — in particular `s ≤ 2B`.  Supersedes the S3 bound `B/(2(β'−β)) + 1`:
no auxiliary margin `β'` needed.  The loose side of `TightLooseDichotomy` is a
decidable search over `s ≤ 2B`. -/
theorem loose_branch_cert_2B (v : Fin 13 → ℤ) (istar : Fin 13) (β : ℝ)
    (hβ : 0 < β) (B : ℤ)
    (hv : ∀ i, i ≠ istar → v i ≠ 0)
    (hBb : ∀ i, i ≠ istar → |v i| ≤ B)
    (hloose : ∃ t : ℝ, ∀ i, i ≠ istar → ∀ m : ℤ, β ≤ |(v i : ℝ) * t - m|) :
    ∃ kk s : ℤ, 0 < s ∧ s ≤ 2 * B ∧ ∀ i, i ≠ istar →
      ⌈β * s⌉ ≤ (v i * kk) % s ∧ (v i * kk) % s ≤ s - ⌈β * s⌉ := by
  classical
  obtain ⟨t, hm⟩ := hloose
  -- a base witness index (the base is nonempty)
  have hne : ∃ i : Fin 13, i ≠ istar := by
    by_contra hcon
    push_neg at hcon
    have h0 := hcon 0
    have h1 := hcon 1
    rw [← h1] at h0
    exact absurd h0 (by decide)
  obtain ⟨iw, hiw⟩ := hne
  -- THE REPLACEMENT TRICK: overwrite the killer slot with a copy of a base
  -- runner; the modified family ranges over base values only, so the Fin-13
  -- machinery applies with all bounds intact.
  set u : Fin 13 → ℤ := fun i => if i = istar then v iw else v i with hudef
  have hu_base : ∀ i, i ≠ istar → u i = v i := by
    intro i hi
    rw [hudef]
    exact if_neg hi
  have hu_star : u istar = v iw := by rw [hudef]; exact if_pos rfl
  have huv : ∀ i, u i ≠ 0 := by
    intro i
    by_cases hi : i = istar
    · rw [hi, hu_star]; exact hv iw hiw
    · rw [hu_base i hi]; exact hv i hi
  have huB : ∀ i, |u i| ≤ B := by
    intro i
    by_cases hi : i = istar
    · rw [hi, hu_star]; exact hBb iw hiw
    · rw [hu_base i hi]; exact hBb i hi
  have hum : ∀ i, ∀ m : ℤ, β ≤ |(u i : ℝ) * t - m| := by
    intro i m
    by_cases hi : i = istar
    · rw [hi, hu_star]; exact hm iw hiw m
    · rw [hu_base i hi]; exact hm i hi m
  obtain ⟨i, j, m, hdom⟩ := grid_margin_domination u huv β hβ t hum
  set s : ℤ := |u i| + |u j| with hsdef
  have hs0 : 0 < s := by
    have h1 : (0 : ℤ) < |u i| := abs_pos.mpr (huv i)
    have h2 : (0 : ℤ) < |u j| := abs_pos.mpr (huv j)
    omega
  refine ⟨m, s, hs0, ?_, ?_⟩
  · have h1 := huB i
    have h2 := huB j
    omega
  · intro i' hi'
    have hmods := mod_of_margin_at_rat u β m s hs0 (fun i'' m' => hdom i'' m') i'
    rw [hu_base i' hi'] at hmods
    exact hmods

#print axioms maximizer_on_grid
#print axioms grid_margin_domination
#print axioms loose_branch_cert_2B

end GridAttainment
end LonelyRunner
