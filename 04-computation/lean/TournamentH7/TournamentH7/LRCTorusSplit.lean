/-
  TournamentH7.LRCTorusSplit — THE COUPLED-TORUS SPLIT RUNG, ρ-parametric
  (kind-pasteur-2026-07-06-S20, HYP-4247).

  THE OBJECT (mac-mini HYP-4262's extracted reduction): (G) — gap emptiness in
  (1/13, 2/25) — decomposes as (A) coupled proper 2-subtorus window emptiness
  + (C) a finite 1-dim census.  The dangerous 2-tori are COUPLED lift-family
  limits: base speeds `w i` with teeth `‖w i · t‖`, lifted slopes `r i > 0`
  with couplings `a i` and teeth `‖r i · θ + a i · t‖` on the (t, θ)-torus.

  THE RUNG (`torus_split_rung`): if such a system is ρ-covered at EVERY
  (t, θ) with ρ ≤ 1/12, then `2ρ·(#lifted) ≥ 1` — the fee-mean/density wall
  as a THEOREM on the torus.  Proof, measure-free:
    · citation on the base (`cite_margin_gen`, ≤ 11 base runners, margin
      1/12 ≥ ρ at some t₀) forces the LIFTED combs alone to ρ-cover the
      θ-circle at t₀;
    · on a `D`-point θ-grid each lifted comb hosts at most `2ρ·D + 3·r i`
      points — the S19 sharp visit count made ρ-PARAMETRIC
      (`tooth_visit_count_rho`, additive-offset form);
    · `D ≤ Σ counts` ⟹ `D(1 − 2ρl) ≤ 3Σr` ⟹ contradiction for `D` large.

  COROLLARIES:
    · `torus_clear_point`: `2ρl < 1` ⟹ some (t, θ) is ρ-clear.
    · `torus_split_gap` (ρ = 2/25): a covered system needs ≥ 7 lifted;
      `torus_clear_gap`: EVERY coupled 2-torus system with ≤ 6 lifted has a
      2/25-clear point — the (A) window is EMPTY of ≤ 6-lifted coupled
      values (all `l ≤ 6` lift-floor LIMIT hypotheses, qualitatively, in one
      theorem; product tori = the `a = 0` special case).
    · `torus_split_cell38` (ρ = 3/38): ≥ 7 lifted — mac-mini-S3's density
      wall (`2ρc ≥ 1`) in the torus dialect.

  The `l ≥ 7` coupled residual is the genuine pole stratum (25/4 again),
  attacked by opus's two-band transport lane; numerics
  (lrc_torus_split_rung_kps_S20.py): even l = 7 systems leave ≥ 14%
  θ-uncovered at small frequencies.

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCClusterGcdSharp

namespace LonelyRunner
namespace TorusSplit

open TeethR ClusterGcd ClusterGcdSharp
open scoped Classical

/-! ## the ρ-parametric wrapped-arc count -/

/-- Grid points of `range p` in the wrapped arc `∃ m, |X + k/p − m| < ρ`:
at most `2ρp + 3`. -/
lemma arc_grid_count_rho (p : ℕ) (hp : 0 < p) (X : ℝ) (ρ : ℝ)
    (hρ0 : 0 < ρ) (hρh : ρ ≤ 1/2) :
    (((Finset.range p).filter
        (fun k : ℕ => ∃ m : ℤ, |X + (k : ℝ) / p - m| < ρ)).card : ℝ)
      ≤ 2 * ρ * p + 3 := by
  have hpR : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  set M₀ : ℤ := ⌈X - ρ⌉ with hM₀
  -- every witness integer is M₀ or M₀ + 1
  have hwitness : ∀ k : ℕ, k ∈ Finset.range p →
      (∃ m : ℤ, |X + (k : ℝ) / p - m| < ρ) →
      (∃ m : ℤ, (m = M₀ ∨ m = M₀ + 1) ∧ |X + (k : ℝ) / p - m| < ρ) := by
    intro k hk ⟨m, hm⟩
    refine ⟨m, ?_, hm⟩
    have hk0 : (0 : ℝ) ≤ (k : ℝ) / p := by positivity
    have hk1 : (k : ℝ) / p < 1 := by
      rw [div_lt_one hpR]
      exact_mod_cast Finset.mem_range.mp hk
    have habs := abs_lt.mp hm
    have hlo : X - ρ < (m : ℝ) := by linarith
    have hhi : (m : ℝ) < X + 1 + ρ := by linarith
    have hM₀le : M₀ ≤ m := by
      rw [hM₀]
      exact Int.ceil_le.mpr (le_of_lt hlo)
    have hM₀ge : (X - ρ : ℝ) ≤ (M₀ : ℝ) := by
      rw [hM₀]
      exact Int.le_ceil _
    have hmM₀ : (m : ℝ) - (M₀ : ℝ) < 2 := by linarith
    have : m - M₀ < 2 := by exact_mod_cast hmM₀
    omega
  -- split by the two witnesses
  have hsub : (Finset.range p).filter
      (fun k : ℕ => ∃ m : ℤ, |X + (k : ℝ) / p - m| < ρ)
      ⊆ ((Finset.range p).filter
          (fun k : ℕ => |X + (k : ℝ) / p - M₀| < ρ))
        ∪ ((Finset.range p).filter
          (fun k : ℕ => |X + (k : ℝ) / p - (M₀ + 1 : ℤ)| < ρ)) := by
    intro k hk
    rw [Finset.mem_filter] at hk
    obtain ⟨hkr, hex⟩ := hk
    obtain ⟨m, hm01, hm⟩ := hwitness k hkr hex
    rw [Finset.mem_union, Finset.mem_filter, Finset.mem_filter]
    rcases hm01 with rfl | rfl
    · exact Or.inl ⟨hkr, hm⟩
    · exact Or.inr ⟨hkr, by exact_mod_cast hm⟩
  -- the two pieces: k-intervals (A, B) and (A+p, B+p), clipped to [0, p)
  set A : ℝ := ((M₀ : ℝ) - X - ρ) * p with hA
  set B : ℝ := ((M₀ : ℝ) - X + ρ) * p with hB
  have hAlow : -(2 * ρ) * p ≤ A := by
    rw [hA]
    have h1 : (X - ρ : ℝ) ≤ (M₀ : ℝ) := by rw [hM₀]; exact Int.le_ceil _
    nlinarith [hpR]
  have hBp : B < p := by
    rw [hB]
    have h2 : (M₀ : ℝ) < X - ρ + 1 := by
      rw [hM₀]
      exact Int.ceil_lt_add_one _
    nlinarith [hpR]
  have hBA : B - A = 2 * ρ * p := by rw [hA, hB]; ring
  have h2rp : (0 : ℝ) ≤ 2 * ρ * p := by positivity
  -- piece 1: interval (max A (−1), B)
  have hc1 : (((Finset.range p).filter
      (fun k : ℕ => |X + (k : ℝ) / p - M₀| < ρ)).card : ℝ)
      ≤ B - max A (-1) + 1 := by
    apply range_filter_interval_card p (max A (-1)) B
      (by
        rcases le_or_gt A (-1) with h | h
        · rw [max_eq_right h]; linarith [hAlow, hBA]
        · rw [max_eq_left (le_of_lt h)]; linarith [hBA])
    intro k hk hP
    have habs := abs_lt.mp hP
    have hup : (k : ℝ) < B := by
      rw [hB]
      have h1 : (k : ℝ) / p < (M₀ : ℝ) - X + ρ := by linarith
      have h2 := (div_lt_iff₀ hpR).mp h1
      linarith
    have hdown : A < (k : ℝ) := by
      rw [hA]
      have h1 : (M₀ : ℝ) - X - ρ < (k : ℝ) / p := by linarith
      have h2 := (lt_div_iff₀ hpR).mp h1
      linarith
    have hneg1 : (-1 : ℝ) < (k : ℝ) := by
      have h0 : (0 : ℝ) ≤ (k : ℝ) := Nat.cast_nonneg k
      linarith
    exact ⟨max_lt_iff.mpr ⟨hdown, hneg1⟩, hup⟩
  -- piece 2: the p-translate (A+p, p), empty when 0 ≤ A
  have hc2 : (((Finset.range p).filter
      (fun k : ℕ => |X + (k : ℝ) / p - (M₀ + 1 : ℤ)| < ρ)).card : ℝ)
      ≤ max (1 - A) 0 := by
    rcases le_or_gt 0 A with hApos | hAneg
    · have hempty : (Finset.range p).filter
          (fun k : ℕ => |X + (k : ℝ) / p - (M₀ + 1 : ℤ)| < ρ) = ∅ := by
        rw [Finset.filter_eq_empty_iff]
        intro k hk hP
        have habs := abs_lt.mp hP
        have hk1 : (k : ℝ) / p < 1 := by
          rw [div_lt_one hpR]
          exact_mod_cast Finset.mem_range.mp hk
        have hcoef : 0 ≤ (M₀ : ℝ) - X - ρ := by
          have h := div_nonneg hApos hpR.le
          rw [hA] at h
          rwa [mul_div_cancel_right₀ _ hpR.ne'] at h
        have hlow : ((M₀ : ℝ) + 1) - X - ρ < (k : ℝ) / p := by
          have h1 := habs.1
          push_cast at h1
          linarith
        linarith
      rw [hempty]
      simp only [Finset.card_empty, Nat.cast_zero]
      exact le_max_right _ _
    · have hbound : (((Finset.range p).filter
          (fun k : ℕ => |X + (k : ℝ) / p - (M₀ + 1 : ℤ)| < ρ)).card : ℝ)
          ≤ (p : ℝ) - (A + p) + 1 := by
        apply range_filter_interval_card p (A + p) (p : ℝ) (by linarith)
        intro k hk hP
        have habs := abs_lt.mp hP
        constructor
        · have h1 : ((M₀ : ℝ) + 1) - X - ρ < (k : ℝ) / p := by
            have h1' := habs.1
            push_cast at h1'
            linarith
          have h2 := (lt_div_iff₀ hpR).mp h1
          have h3 : (((M₀ : ℝ) + 1) - X - ρ) * p = A + p := by rw [hA]; ring
          linarith
        · exact_mod_cast Finset.mem_range.mp hk
      have hlen : (p : ℝ) - (A + p) + 1 = 1 - A := by ring
      rw [hlen] at hbound
      exact hbound.trans (le_max_left _ _)
  -- assemble and clip
  have hle := Finset.card_le_card hsub
  have hu := Finset.card_union_le
    ((Finset.range p).filter (fun k : ℕ => |X + (k : ℝ) / p - M₀| < ρ))
    ((Finset.range p).filter (fun k : ℕ => |X + (k : ℝ) / p - (M₀ + 1 : ℤ)| < ρ))
  have h1 : (((Finset.range p).filter
      (fun k : ℕ => ∃ m : ℤ, |X + (k : ℝ) / p - m| < ρ)).card : ℝ)
      ≤ (((Finset.range p).filter
          (fun k : ℕ => |X + (k : ℝ) / p - M₀| < ρ)).card : ℝ)
        + (((Finset.range p).filter
          (fun k : ℕ => |X + (k : ℝ) / p - (M₀ + 1 : ℤ)| < ρ)).card : ℝ) := by
    exact_mod_cast le_trans hle hu
  have hclip : (B - max A (-1) + 1) + max (1 - A) 0 ≤ 2 * ρ * p + 3 := by
    rcases le_or_gt 0 A with h0 | h0
    · rw [max_eq_left (by linarith : (-1:ℝ) ≤ A)]
      have h1' : max (1 - A) 0 ≤ 1 := max_le (by linarith) (by norm_num)
      linarith [hBA]
    · have h1' : -1 ≤ max A (-1) := le_max_right A (-1)
      have h2' : max (1 - A) 0 = 1 - A := max_eq_left (by linarith)
      rw [h2']
      linarith [hBA]
  linarith [h1, hc1, hc2, hclip]

/-! ## the ρ-parametric sharp visit count (additive-offset form) -/

/-- Among `D` grid points, the inhomogeneous comb `∃ m, |c + w·(j/D) − m| < ρ`
hosts at most `2ρD + 3w`. -/
lemma tooth_visit_count_rho (w : ℤ) (hw : 0 < w) (D : ℕ) (hD : 0 < D)
    (c : ℝ) (ρ : ℝ) (hρ0 : 0 < ρ) (hρh : ρ ≤ 1/2) :
    (((Finset.range D).filter
        (fun j : ℕ => ∃ m : ℤ, |c + (w : ℝ) * ((j : ℝ) / D) - m| < ρ)).card : ℝ)
      ≤ 2 * ρ * D + 3 * w := by
  set W : ℕ := w.toNat with hW
  have hW0 : 0 < W := by rw [hW]; omega
  have hWw : (W : ℤ) = w := by rw [hW]; omega
  have hWR : ((W : ℕ) : ℝ) = ((w : ℤ) : ℝ) := by exact_mod_cast hWw
  set g : ℕ := Nat.gcd W D with hg
  have hg0 : 0 < g := Nat.gcd_pos_of_pos_left _ hW0
  set p : ℕ := D / g with hp
  set w' : ℕ := W / g with hw'
  have hgdvdD : g ∣ D := by rw [hg]; exact Nat.gcd_dvd_right W D
  have hgdvdW : g ∣ W := by rw [hg]; exact Nat.gcd_dvd_left W D
  have hDgp : D = g * p := by rw [hp]; exact (Nat.mul_div_cancel' hgdvdD).symm
  have hWgw : W = g * w' := by rw [hw']; exact (Nat.mul_div_cancel' hgdvdW).symm
  have hp0 : 0 < p := by
    rcases Nat.eq_zero_or_pos p with h0 | h
    · rw [h0, Nat.mul_zero] at hDgp
      omega
    · exact h
  have hco : Nat.Coprime w' p := by
    rw [hw', hp, hg]
    exact Nat.coprime_div_gcd_div_gcd (by rw [← hg]; exact hg0)
  have hpRl : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp0
  have hDRl : (0 : ℝ) < (D : ℝ) := by exact_mod_cast hD
  have hWpD : (W : ℝ) * p = (w' : ℝ) * D := by
    have h1 : W * p = w' * D := by
      calc W * p = (g * w') * p := by rw [← hWgw]
        _ = w' * (g * p) := by ring
        _ = w' * D := by rw [← hDgp]
    exact_mod_cast h1
  -- the value identity: w·(j/D) = (w′j)/p
  have hval : ∀ j : ℕ, (w : ℝ) * ((j : ℝ) / D) = ((w' : ℝ) * j) / p := by
    intro j
    rw [← hWR]
    rw [eq_div_iff hpRl.ne']
    calc (W : ℝ) * ((j : ℝ) / D) * p = ((j : ℝ) / D) * ((W : ℝ) * p) := by ring
      _ = ((j : ℝ) / D) * ((w' : ℝ) * D) := by rw [hWpD]
      _ = ((w' : ℝ) * j) * ((D : ℝ) / D) := by ring
      _ = ((w' : ℝ) * j) * 1 := by rw [div_self hDRl.ne']
      _ = (w' : ℝ) * j := by ring
  -- split the fraction at the mod
  have hmoddiv : ∀ j : ℕ, ((w' : ℝ) * j) / p
      = (((w' * j) % p : ℕ) : ℝ) / p + (((w' * j) / p : ℕ) : ℝ) := by
    intro j
    have hde : (w' * j) % p + p * ((w' * j) / p) = w' * j := Nat.mod_add_div (w' * j) p
    have hdeR : (((w' * j) % p : ℕ) : ℝ) + (p : ℝ) * (((w' * j) / p : ℕ) : ℝ)
        = (w' : ℝ) * j := by exact_mod_cast hde
    rw [← hdeR, add_div, mul_div_cancel_left₀ _ hpRl.ne']
  -- the comb condition depends on j only through (w′j) mod p
  have hcond : ∀ j : ℕ, (∃ m : ℤ, |c + (w : ℝ) * ((j : ℝ) / D) - m| < ρ)
      ↔ (∃ m : ℤ, |c + (((w' * j) % p : ℕ) : ℝ) / p - m| < ρ) := by
    intro j
    set q : ℕ := (w' * j) / p with hq
    constructor
    · rintro ⟨m, hm⟩
      refine ⟨m - q, ?_⟩
      rw [hval j, hmoddiv j] at hm
      have hswap : c + (((w' * j) % p : ℕ) : ℝ) / p - ((m - q : ℤ) : ℝ)
          = c + ((((w' * j) % p : ℕ) : ℝ) / p + (q : ℝ)) - (m : ℝ) := by
        push_cast
        ring
      rw [hswap]
      exact hm
    · rintro ⟨m, hm⟩
      refine ⟨m + q, ?_⟩
      rw [hval j, hmoddiv j]
      have hswap : c + ((((w' * j) % p : ℕ) : ℝ) / p + (q : ℝ)) - ((m + q : ℤ) : ℝ)
          = c + (((w' * j) % p : ℕ) : ℝ) / p - (m : ℝ) := by
        push_cast
        ring
      rw [hswap]
      exact hm
  have hfe : (Finset.range D).filter
      (fun j : ℕ => ∃ m : ℤ, |c + (w : ℝ) * ((j : ℝ) / D) - m| < ρ)
      = (Finset.range D).filter
        (fun j : ℕ => ∃ m : ℤ, |c + (((w' * j) % p : ℕ) : ℝ) / p - m| < ρ) :=
    Finset.filter_congr (fun j _ => hcond j)
  rw [hfe]
  -- block split + coprime permutation
  have hper : ∀ j : ℕ,
      (∃ m : ℤ, |c + (((w' * (j + p)) % p : ℕ) : ℝ) / p - m| < ρ)
      ↔ (∃ m : ℤ, |c + (((w' * j) % p : ℕ) : ℝ) / p - m| < ρ) := by
    intro j
    have hmm : (w' * (j + p)) % p = (w' * j) % p := by
      have hexp : w' * (j + p) = w' * j + w' * p := by ring
      rw [hexp, Nat.add_mul_mod_self_right]
    rw [hmm]
  conv_lhs => rw [hDgp]
  rw [filter_range_mul_periodic
    (fun j : ℕ => ∃ m : ℤ, |c + (((w' * j) % p : ℕ) : ℝ) / p - m| < ρ)
    p hper g]
  rw [filter_range_coprime_mul
    (fun k : ℕ => ∃ m : ℤ, |c + (k : ℝ) / p - m| < ρ)
    p w' hp0 hco]
  -- the arc count, times g blocks
  have harc := arc_grid_count_rho p hp0 c ρ hρ0 hρh
  have hgw : (g : ℝ) ≤ ((w : ℤ) : ℝ) := by
    rw [← hWR]
    exact_mod_cast Nat.le_of_dvd hW0 hgdvdW
  have hgD : (g : ℝ) * p = (D : ℝ) := by exact_mod_cast hDgp.symm
  have hg0R : (0 : ℝ) < (g : ℝ) := by exact_mod_cast hg0
  push_cast
  calc (g : ℝ) * (((Finset.range p).filter
      (fun k : ℕ => ∃ m : ℤ, |c + (k : ℝ) / p - m| < ρ)).card : ℝ)
      ≤ (g : ℝ) * (2 * ρ * p + 3) := by
        exact mul_le_mul_of_nonneg_left harc hg0R.le
    _ = 2 * ρ * ((g : ℝ) * p) + 3 * g := by ring
    _ = 2 * ρ * D + 3 * g := by rw [hgD]
    _ ≤ 2 * ρ * D + 3 * w := by linarith [hgw]

/-! ## the split rung -/

/-- **THE COUPLED-TORUS SPLIT RUNG**: a coupled 2-torus system (base teeth
`|w i · t − m|`, lifted teeth `|r i · θ + a i · t − m|`, slopes `r i > 0`)
that is ρ-covered at every (t, θ) with `ρ ≤ 1/12` has `2ρ·(#lifted) ≥ 1`. -/
theorem torus_split_rung (cite : LRCUpTo13) (ρ : ℝ) (hρ0 : 0 < ρ) (hρ12 : ρ ≤ 1/12)
    (L : Finset (Fin 12)) (hL : L.Nonempty)
    (w r a : Fin 12 → ℤ)
    (hw : ∀ i, i ∉ L → w i ≠ 0) (hr : ∀ i, i ∈ L → 0 < r i)
    (hcover : ∀ t θ : ℝ,
      (∃ i, i ∉ L ∧ ∃ m : ℤ, |(w i : ℝ) * t - m| < ρ) ∨
      (∃ i, i ∈ L ∧ ∃ m : ℤ, |(r i : ℝ) * θ + (a i : ℝ) * t - m| < ρ)) :
    (1 : ℝ) ≤ 2 * ρ * L.card := by
  classical
  -- the base, cited (padded to a 12-family)
  set T : Finset (Fin 12) := Finset.univ \ L with hT
  have hTmem : ∀ i, i ∈ T ↔ i ∉ L := by
    intro i
    rw [hT, Finset.mem_sdiff]
    simp only [Finset.mem_univ, true_and]
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
  -- the base margin at t₀ beats ρ
  have hmargin' : ∀ i, i ∉ L → ∀ m : ℤ, ρ ≤ |(w i : ℝ) * t₀ - m| := by
    intro i hi m
    have h := hmargin i ((hTmem i).mpr hi) m
    rw [if_neg hi] at h
    have hcard : ((T.card : ℝ) + 1) ≤ 12 := by
      have : (T.card : ℝ) ≤ 11 := by exact_mod_cast hT11
      linarith
    have h12 : (1:ℝ)/12 ≤ 1 / ((T.card : ℝ) + 1) := by
      apply div_le_div_of_nonneg_left (by norm_num) ?_ hcard
      positivity
    calc ρ ≤ 1/12 := hρ12
      _ ≤ 1 / ((T.card : ℝ) + 1) := h12
      _ ≤ |(w i : ℝ) * t₀ - m| := h
  -- at t₀ the lifted combs alone ρ-cover the θ-circle
  have hθcover : ∀ θ : ℝ, ∃ i, i ∈ L ∧
      ∃ m : ℤ, |(r i : ℝ) * θ + (a i : ℝ) * t₀ - m| < ρ := by
    intro θ
    rcases hcover t₀ θ with ⟨i, hiL, m, hm⟩ | h
    · exfalso
      have := hmargin' i hiL m
      linarith
    · exact h
  -- the count contradiction at a large grid
  set R : ℤ := ∑ i ∈ L, r i with hR
  have hR0 : 0 < R := Finset.sum_pos (fun i hi => hr i hi) hL
  have hRR : (0 : ℝ) < (R : ℝ) := by exact_mod_cast hR0
  by_contra hcon
  rw [not_le] at hcon
  have hgap : 0 < 1 - 2 * ρ * (L.card : ℝ) := by linarith
  obtain ⟨D, hD⟩ := exists_nat_gt (3 * (R : ℝ) / (1 - 2 * ρ * (L.card : ℝ)))
  have hDpos0 : (0 : ℝ) < 3 * (R : ℝ) / (1 - 2 * ρ * (L.card : ℝ)) :=
    div_pos (by linarith) hgap
  have hD0 : 0 < D := by
    have : (0 : ℝ) < (D : ℝ) := lt_trans hDpos0 hD
    exact_mod_cast this
  -- classify each grid point to a lifted comb
  have hchoice : ∀ j : ℕ, ∃ i : Fin 12, j ∈ Finset.range D →
      (i ∈ L ∧ ∃ m : ℤ, |(r i : ℝ) * ((j : ℝ) / D) + (a i : ℝ) * t₀ - m| < ρ) := by
    intro j
    by_cases hj : j ∈ Finset.range D
    · obtain ⟨i, hiL, hm⟩ := hθcover ((j : ℝ) / D)
      exact ⟨i, fun _ => ⟨hiL, hm⟩⟩
    · exact ⟨0, fun hc => absurd hc hj⟩
  choose F hFspec using hchoice
  have hFmem : ∀ j ∈ Finset.range D, F j ∈ L := fun j hj => (hFspec j hj).1
  have hsplit : (Finset.range D).card
      = ∑ i ∈ L, ((Finset.range D).filter (fun j => F j = i)).card :=
    Finset.card_eq_sum_card_fiberwise hFmem
  -- per-fiber: inside the ρ-parametric visit filter of r i
  have hfiber : ∀ i ∈ L, (((Finset.range D).filter (fun j => F j = i)).card : ℝ)
      ≤ 2 * ρ * D + 3 * ((r i : ℤ) : ℝ) := by
    intro i hi
    have hsub : (Finset.range D).filter (fun j => F j = i)
        ⊆ (Finset.range D).filter
          (fun j : ℕ => ∃ m : ℤ, |((a i : ℝ) * t₀) + (r i : ℝ) * ((j : ℝ) / D) - m| < ρ) := by
      intro j hj
      rw [Finset.mem_filter] at hj ⊢
      obtain ⟨hjr, hjF⟩ := hj
      have hm := (hFspec j hjr).2
      rw [hjF] at hm
      obtain ⟨m, hm⟩ := hm
      refine ⟨hjr, m, ?_⟩
      have harg : (a i : ℝ) * t₀ + (r i : ℝ) * ((j : ℝ) / D) - m
          = (r i : ℝ) * ((j : ℝ) / D) + (a i : ℝ) * t₀ - m := by ring
      rw [harg]
      exact hm
    have hcount := tooth_visit_count_rho (r i) (hr i hi) D hD0
      ((a i : ℝ) * t₀) ρ hρ0 (by linarith)
    calc (((Finset.range D).filter (fun j => F j = i)).card : ℝ)
        ≤ (((Finset.range D).filter
          (fun j : ℕ => ∃ m : ℤ, |((a i : ℝ) * t₀) + (r i : ℝ) * ((j : ℝ) / D) - m| < ρ)).card : ℝ) := by
          exact_mod_cast Finset.card_le_card hsub
      _ ≤ 2 * ρ * D + 3 * ((r i : ℤ) : ℝ) := hcount
  -- assemble over L
  have hsum : ((D : ℕ) : ℝ)
      ≤ 2 * ρ * D * L.card + 3 * ((R : ℤ) : ℝ) := by
    have hcards : ((D : ℕ) : ℝ)
        = ∑ i ∈ L, (((Finset.range D).filter (fun j => F j = i)).card : ℝ) := by
      have h0 : ((Finset.range D).card : ℝ)
          = ∑ i ∈ L, (((Finset.range D).filter (fun j => F j = i)).card : ℝ) := by
        exact_mod_cast congrArg (Nat.cast : ℕ → ℝ) hsplit
      rw [← h0, Finset.card_range]
    conv_lhs => rw [hcards]
    have hRcast : ((R : ℤ) : ℝ) = ∑ i ∈ L, ((r i : ℤ) : ℝ) := by
      rw [hR]
      push_cast
      ring
    calc (∑ i ∈ L, (((Finset.range D).filter (fun j => F j = i)).card : ℝ))
        ≤ ∑ i ∈ L, (2 * ρ * (D : ℝ) + 3 * ((r i : ℤ) : ℝ)) :=
          Finset.sum_le_sum hfiber
      _ = 2 * ρ * D * L.card + 3 * ((R : ℤ) : ℝ) := by
          rw [Finset.sum_add_distrib, Finset.sum_const, ← Finset.mul_sum, hRcast]
          ring
  -- D(1 − 2ρl) ≤ 3R against D > 3R/(1 − 2ρl)
  have hDbig : 3 * (R : ℝ) < (1 - 2 * ρ * (L.card : ℝ)) * D := by
    rw [div_lt_iff₀ hgap] at hD
    linarith [hD]
  nlinarith [hsum, hDbig]

/-! ## corollaries -/

/-- **THE CLEAR-POINT COROLLARY**: a coupled 2-torus system with
`2ρ·(#lifted) < 1` has a ρ-clear point. -/
theorem torus_clear_point (cite : LRCUpTo13) (ρ : ℝ) (hρ0 : 0 < ρ) (hρ12 : ρ ≤ 1/12)
    (L : Finset (Fin 12)) (hL : L.Nonempty)
    (w r a : Fin 12 → ℤ)
    (hw : ∀ i, i ∉ L → w i ≠ 0) (hr : ∀ i, i ∈ L → 0 < r i)
    (hsmall : 2 * ρ * (L.card : ℝ) < 1) :
    ∃ t θ : ℝ, (∀ i, i ∉ L → ∀ m : ℤ, ρ ≤ |(w i : ℝ) * t - m|) ∧
               (∀ i, i ∈ L → ∀ m : ℤ, ρ ≤ |(r i : ℝ) * θ + (a i : ℝ) * t - m|) := by
  by_contra hcon
  have hcover : ∀ t θ : ℝ,
      (∃ i, i ∉ L ∧ ∃ m : ℤ, |(w i : ℝ) * t - m| < ρ) ∨
      (∃ i, i ∈ L ∧ ∃ m : ℤ, |(r i : ℝ) * θ + (a i : ℝ) * t - m| < ρ) := by
    intro t θ
    by_contra hpoint
    apply hcon
    refine ⟨t, θ, ?_, ?_⟩
    · intro i hi m
      by_contra hlt
      rw [not_le] at hlt
      exact hpoint (Or.inl ⟨i, hi, m, hlt⟩)
    · intro i hi m
      by_contra hlt
      rw [not_le] at hlt
      exact hpoint (Or.inr ⟨i, hi, m, hlt⟩)
  have := torus_split_rung cite ρ hρ0 hρ12 L hL w r a hw hr hcover
  linarith

/-- At the gap ceiling ρ = 2/25: a covered coupled system needs ≥ 7 lifted. -/
theorem torus_split_gap (cite : LRCUpTo13)
    (L : Finset (Fin 12)) (hL : L.Nonempty)
    (w r a : Fin 12 → ℤ)
    (hw : ∀ i, i ∉ L → w i ≠ 0) (hr : ∀ i, i ∈ L → 0 < r i)
    (hcover : ∀ t θ : ℝ,
      (∃ i, i ∉ L ∧ ∃ m : ℤ, |(w i : ℝ) * t - m| < 2/25) ∨
      (∃ i, i ∈ L ∧ ∃ m : ℤ, |(r i : ℝ) * θ + (a i : ℝ) * t - m| < 2/25)) :
    7 ≤ L.card := by
  have h := torus_split_rung cite (2/25) (by norm_num) (by norm_num)
    L hL w r a hw hr hcover
  have h1 : (25 : ℝ) ≤ 4 * (L.card : ℝ) := by linarith
  have h2 : (25 : ℕ) ≤ 4 * L.card := by exact_mod_cast h1
  omega

/-- **THE (A) WINDOW l ≤ 6 KILL**: every coupled 2-torus system with at most
6 lifted runners has a 2/25-clear point — its value is ≥ 2/25, outside the
open gap (1/13, 2/25).  Product tori are the `a = 0` special case. -/
theorem torus_clear_gap (cite : LRCUpTo13)
    (L : Finset (Fin 12)) (hL : L.Nonempty) (hL6 : L.card ≤ 6)
    (w r a : Fin 12 → ℤ)
    (hw : ∀ i, i ∉ L → w i ≠ 0) (hr : ∀ i, i ∈ L → 0 < r i) :
    ∃ t θ : ℝ, (∀ i, i ∉ L → ∀ m : ℤ, 2/25 ≤ |(w i : ℝ) * t - m|) ∧
               (∀ i, i ∈ L → ∀ m : ℤ, 2/25 ≤ |(r i : ℝ) * θ + (a i : ℝ) * t - m|) := by
  apply torus_clear_point cite (2/25) (by norm_num) (by norm_num) L hL w r a hw hr
  have h6 : (L.card : ℝ) ≤ 6 := by exact_mod_cast hL6
  linarith

/-- At the 3/38 cell (mac-mini-S3's density wall in the torus dialect):
a covered coupled system needs ≥ 7 lifted. -/
theorem torus_split_cell38 (cite : LRCUpTo13)
    (L : Finset (Fin 12)) (hL : L.Nonempty)
    (w r a : Fin 12 → ℤ)
    (hw : ∀ i, i ∉ L → w i ≠ 0) (hr : ∀ i, i ∈ L → 0 < r i)
    (hcover : ∀ t θ : ℝ,
      (∃ i, i ∉ L ∧ ∃ m : ℤ, |(w i : ℝ) * t - m| < 3/38) ∨
      (∃ i, i ∈ L ∧ ∃ m : ℤ, |(r i : ℝ) * θ + (a i : ℝ) * t - m| < 3/38)) :
    7 ≤ L.card := by
  have h := torus_split_rung cite (3/38) (by norm_num) (by norm_num)
    L hL w r a hw hr hcover
  have h1 : (19 : ℝ) ≤ 3 * (L.card : ℝ) := by linarith
  have h2 : (19 : ℕ) ≤ 3 * L.card := by exact_mod_cast h1
  omega

/-- **THE FORCED RECTANGLE** (the l ≥ 7 residual's stage): in a 2/25-covered
coupled system, the lifted combs must θ-cover an ENTIRE explicit t-interval
around the base citation point — width `|t − t₀|·V ≤ 1/300` where `V` bounds
the base heights.  (`1/12 − 1/300 = 2/25` exactly: the S6 compression seam.)
The residual is thus a rectangle-covering rigidity question for ≥ 7 bands of
distinct slopes. -/
theorem torus_forced_rectangle (cite : LRCUpTo13)
    (L : Finset (Fin 12)) (hL : L.Nonempty)
    (w r a : Fin 12 → ℤ) (V : ℤ)
    (hw : ∀ i, i ∉ L → w i ≠ 0) (hwV : ∀ i, i ∉ L → |w i| ≤ V)
    (hcover : ∀ t θ : ℝ,
      (∃ i, i ∉ L ∧ ∃ m : ℤ, |(w i : ℝ) * t - m| < 2/25) ∨
      (∃ i, i ∈ L ∧ ∃ m : ℤ, |(r i : ℝ) * θ + (a i : ℝ) * t - m| < 2/25)) :
    ∃ t₀ : ℝ, ∀ t : ℝ, |t - t₀| * (V : ℝ) ≤ 1/300 → ∀ θ : ℝ,
      ∃ i, i ∈ L ∧ ∃ m : ℤ, |(r i : ℝ) * θ + (a i : ℝ) * t - m| < 2/25 := by
  classical
  set T : Finset (Fin 12) := Finset.univ \ L with hT
  have hTmem : ∀ i, i ∈ T ↔ i ∉ L := by
    intro i
    rw [hT, Finset.mem_sdiff]
    simp only [Finset.mem_univ, true_and]
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
  have hmargin' : ∀ i, i ∉ L → ∀ m : ℤ, (1:ℝ)/12 ≤ |(w i : ℝ) * t₀ - m| := by
    intro i hi m
    have h := hmargin i ((hTmem i).mpr hi) m
    rw [if_neg hi] at h
    have hcard : ((T.card : ℝ) + 1) ≤ 12 := by
      have : (T.card : ℝ) ≤ 11 := by exact_mod_cast hT11
      linarith
    have h12 : (1:ℝ)/12 ≤ 1 / ((T.card : ℝ) + 1) := by
      apply div_le_div_of_nonneg_left (by norm_num) ?_ hcard
      positivity
    calc (1:ℝ)/12 ≤ 1 / ((T.card : ℝ) + 1) := h12
      _ ≤ |(w i : ℝ) * t₀ - m| := h
  refine ⟨t₀, ?_⟩
  intro t ht θ
  rcases hcover t θ with ⟨i, hiL, m, hm⟩ | h
  · exfalso
    have h1 := hmargin' i hiL m
    have hAB : ((w i : ℝ) * t₀ - m) - ((w i : ℝ) * t - m) = (w i : ℝ) * (t₀ - t) := by
      ring
    have htri := abs_sub_abs_le_abs_sub ((w i : ℝ) * t₀ - m) ((w i : ℝ) * t - m)
    rw [hAB, abs_mul] at htri
    have h4 : |(w i : ℝ)| ≤ (V : ℝ) := by
      have := hwV i hiL
      have hcast : ((|w i| : ℤ) : ℝ) ≤ ((V : ℤ) : ℝ) := by exact_mod_cast this
      rwa [show ((|w i| : ℤ) : ℝ) = |(w i : ℝ)| by push_cast; ring] at hcast
    have h5 : |t₀ - t| = |t - t₀| := abs_sub_comm _ _
    have h6 : |(w i : ℝ)| * |t₀ - t| ≤ 1/300 := by
      rw [h5]
      calc |(w i : ℝ)| * |t - t₀| ≤ (V : ℝ) * |t - t₀| :=
            mul_le_mul_of_nonneg_right h4 (abs_nonneg _)
        _ = |t - t₀| * (V : ℝ) := by ring
        _ ≤ 1/300 := ht
    -- |w t − m| ≥ 1/12 − 1/300 = 2/25 exactly
    have h7 : (2:ℝ)/25 ≤ |(w i : ℝ) * t - m| := by linarith [h1, h6, htri]
    linarith [hm]
  · exact h

/-- **PRODUCT TORI ARE DEAD** (double citation): a PROPER coupled system
(nonempty base) with all couplings zero cannot be 2/25-covered — the lifted
combs would have to θ-cover statically, against the citation on the slopes.
(Properness matters: with an empty base and zero couplings the system
degenerates to a plain 12-family, where `{1..12}` IS 2/25-covered.)  So the
(A) residual after `torus_clear_gap` needs `l ≥ 7` AND a nonzero coupling. -/
theorem torus_product_dead (cite : LRCUpTo13)
    (L : Finset (Fin 12)) (hL : L.Nonempty) (hLne : L ≠ Finset.univ)
    (w r a : Fin 12 → ℤ)
    (hw : ∀ i, i ∉ L → w i ≠ 0) (hrne : ∀ i, i ∈ L → r i ≠ 0)
    (ha0 : ∀ i, i ∈ L → a i = 0)
    (hcover : ∀ t θ : ℝ,
      (∃ i, i ∉ L ∧ ∃ m : ℤ, |(w i : ℝ) * t - m| < 2/25) ∨
      (∃ i, i ∈ L ∧ ∃ m : ℤ, |(r i : ℝ) * θ + (a i : ℝ) * t - m| < 2/25)) :
    False := by
  classical
  set T : Finset (Fin 12) := Finset.univ \ L with hT
  have hTmem : ∀ i, i ∈ T ↔ i ∉ L := by
    intro i
    rw [hT, Finset.mem_sdiff]
    simp only [Finset.mem_univ, true_and]
  have hTcard : T.card ≤ 12 := by
    calc T.card ≤ Finset.univ.card := Finset.card_le_card (Finset.subset_univ T)
      _ = 12 := by rw [Finset.card_univ, Fintype.card_fin]
  have hT11 : T.card ≤ 11 := by
    have hL1 : 1 ≤ L.card := Finset.card_pos.mpr hL
    have : T.card = 12 - L.card := by
      rw [hT, Finset.card_univ_diff, Fintype.card_fin]
    omega
  -- first citation: the base clears at t₀
  obtain ⟨t₀, hmargin⟩ := cite_margin_gen cite
    (fun i => if i ∈ L then 1 else w i) T hTcard
    (fun i hi => by
      show (if i ∈ L then 1 else w i) ≠ 0
      rw [if_neg ((hTmem i).mp hi)]
      exact hw i ((hTmem i).mp hi))
  have hmargin' : ∀ i, i ∉ L → ∀ m : ℤ, (2:ℝ)/25 < |(w i : ℝ) * t₀ - m| := by
    intro i hi m
    have h := hmargin i ((hTmem i).mpr hi) m
    rw [if_neg hi] at h
    have hcard : ((T.card : ℝ) + 1) ≤ 12 := by
      have : (T.card : ℝ) ≤ 11 := by exact_mod_cast hT11
      linarith
    have h12 : (1:ℝ)/12 ≤ 1 / ((T.card : ℝ) + 1) := by
      apply div_le_div_of_nonneg_left (by norm_num) ?_ hcard
      positivity
    calc (2:ℝ)/25 < 1/12 := by norm_num
      _ ≤ 1 / ((T.card : ℝ) + 1) := h12
      _ ≤ |(w i : ℝ) * t₀ - m| := h
  -- second citation: the slopes clear at θ*
  have hL11 : L.card ≤ 11 := by
    have hlt : L.card < Finset.univ.card :=
      Finset.card_lt_card (Finset.ssubset_univ_iff.mpr hLne)
    rw [Finset.card_univ, Fintype.card_fin] at hlt
    omega
  obtain ⟨θs, hmargin2⟩ := cite_margin_gen cite
    (fun i => if i ∈ L then r i else 1) L
    (by
      calc L.card ≤ Finset.univ.card := Finset.card_le_card (Finset.subset_univ L)
        _ = 12 := by rw [Finset.card_univ, Fintype.card_fin])
    (fun i hi => by
      show (if i ∈ L then r i else 1) ≠ 0
      rw [if_pos hi]
      exact hrne i hi)
  have hmargin2' : ∀ i, i ∈ L → ∀ m : ℤ, (2:ℝ)/25 < |(r i : ℝ) * θs - m| := by
    intro i hi m
    have h := hmargin2 i hi m
    rw [if_pos hi] at h
    have hcard : ((L.card : ℝ) + 1) ≤ 12 := by
      have : (L.card : ℝ) ≤ 11 := by exact_mod_cast hL11
      linarith
    have h12 : (1:ℝ)/12 ≤ 1 / ((L.card : ℝ) + 1) := by
      apply div_le_div_of_nonneg_left (by norm_num) ?_ hcard
      positivity
    calc (2:ℝ)/25 < 1/12 := by norm_num
      _ ≤ 1 / ((L.card : ℝ) + 1) := h12
      _ ≤ |(r i : ℝ) * θs - m| := h
  -- the cover at (t₀, θ*) has nowhere to go
  rcases hcover t₀ θs with ⟨i, hiL, m, hm⟩ | ⟨i, hiL, m, hm⟩
  · exact absurd hm (not_lt.mpr (le_of_lt (hmargin' i hiL m)))
  · rw [ha0 i hiL] at hm
    push_cast at hm
    rw [zero_mul, add_zero] at hm
    exact absurd hm (not_lt.mpr (le_of_lt (hmargin2' i hiL m)))

#print axioms torus_split_rung
#print axioms torus_clear_gap
#print axioms torus_split_cell38
#print axioms torus_forced_rectangle
#print axioms torus_product_dead

end TorusSplit
end LonelyRunner
