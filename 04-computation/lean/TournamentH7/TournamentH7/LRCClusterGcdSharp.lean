/-
  TournamentH7.LRCClusterGcdSharp — THE SHARP VISIT COUNT AND THE FULL
  CLUSTER-GCD LADDER, |S| ≤ 6
  (kind-pasteur-2026-07-06-S19, HYP-4237).

  CRITICAL-PATH: mac-mini-S1's k-reduction for the 3/38 Farey cell (HYP-4232)
  composes the cluster-gcd ladder at |S| up to 6; the S18 file formalizes only
  `|S| ≤ 3` (lossy count).  This file delivers the sharp count and the ladder:

  * `tooth_visit_count_sharp`: among `D` equally-spaced copies, a comb `w`
    hosts at most `(4/25)·D + 3·w` — via
      (a) the visit condition is `p`-periodic in `j`, `p = D/gcd(w,D)`
          (`w·p/D = w/gcd ∈ ℤ`): block-split `D = gcd·p`
          (`filter_range_mul_periodic`);
      (b) inside one block the dependence is through `(w′·j) mod p`,
          `w′ = w/gcd` COPRIME to `p`: `k ↦ w′k mod p` permutes `range p`
          (ZMod-free: injective self-map + image extensionality,
          `filter_range_coprime_mul`);
      (c) the permuted condition is membership in ONE wrapped arc of length
          `(4/25)p`: at most two witness integers `M₀, M₀+1`, whose `ρ`-intervals
          are `p`-translates; the clipped two-piece count is `≤ (4/25)p + 3`
          (`arc_grid_count`).

  * `gap_gcd_rung_sharp`:  `(25 − 4·|S|)·d ≤ 75·Σ_{i∈S}|vᵢ|`  for every
    `1 ≤ |S| ≤ 6` and every positive common divisor `d` of the complement —
    the POLE at `|S| = 25/4` now formal.  (The draft's constant 50 reflects a
    tighter clip; 75 is what the uniform two-piece clip yields — the pole is
    what matters.)

  With this, mac-mini's `k ≥ 2` strata are height-bounded FORMALLY at every
  `|S| ≤ 6`; the `|S| ≥ 7` residual stays census-shaped as mapped.

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCClusterGcd

namespace LonelyRunner
namespace ClusterGcdSharp

open TeethR TripleWalk ClusterGcd
open scoped Classical

/-! ## (a) the block split -/

/-- Iterating a period: `P (n·p + j) ↔ P j`. -/
lemma periodic_iterate (P : ℕ → Prop) (p : ℕ)
    (hper : ∀ j, P (j + p) ↔ P j) : ∀ n j, P (n * p + j) ↔ P j := by
  intro n
  induction n with
  | zero => intro j; simp
  | succ m ihm =>
    intro j
    rw [show (m + 1) * p + j = (m * p + j) + p by ring, hper (m * p + j), ihm j]

/-- A `p`-periodic predicate's count over `range (g·p)` is `g` blocks' worth. -/
lemma filter_range_mul_periodic (P : ℕ → Prop) [DecidablePred P] (p : ℕ)
    (hper : ∀ j, P (j + p) ↔ P j) (g : ℕ) :
    ((Finset.range (g * p)).filter P).card
      = g * ((Finset.range p).filter P).card := by
  induction g with
  | zero => simp
  | succ n ih =>
    have hiter : ∀ j, P (n * p + j) ↔ P j := periodic_iterate P p hper n
    rw [show (n + 1) * p = n * p + p by ring, Finset.range_add, Finset.filter_union]
    rw [Finset.card_union_of_disjoint]
    · have hmap : ((Finset.range p).map (addLeftEmbedding (n * p))).filter P
          = ((Finset.range p).filter (fun j => P (n * p + j))).map
              (addLeftEmbedding (n * p)) := by
        rw [Finset.filter_map]
        rfl
      rw [ih, hmap, Finset.card_map]
      have hcongr : (Finset.range p).filter (fun j => P (n * p + j))
          = (Finset.range p).filter P :=
        Finset.filter_congr (fun j _ => hiter j)
      rw [hcongr]
      ring
    · apply Finset.disjoint_filter_filter
      rw [Finset.disjoint_left]
      intro a ha hb
      rw [Finset.mem_range] at ha
      rw [Finset.mem_map] at hb
      obtain ⟨b, hbr, hab⟩ := hb
      simp only [addLeftEmbedding_apply] at hab
      omega

/-! ## (b) the coprime permutation -/

/-- Counting through the unit multiplication `k ↦ (w′·k) mod p`. -/
lemma filter_range_coprime_mul (Q : ℕ → Prop) [DecidablePred Q] (p w' : ℕ)
    (hp : 0 < p) (hco : Nat.Coprime w' p) :
    ((Finset.range p).filter (fun k => Q ((w' * k) % p))).card
      = ((Finset.range p).filter Q).card := by
  have hinj : Set.InjOn (fun k => (w' * k) % p) ↑(Finset.range p) := by
    intro a ha b hb hab
    simp only [Finset.coe_range, Set.mem_Iio] at ha hb
    have hmod : w' * a ≡ w' * b [MOD p] := hab
    have hcan : a ≡ b [MOD p] := Nat.ModEq.cancel_left_of_coprime hco.symm hmod
    have hcan' : a % p = b % p := hcan
    rwa [Nat.mod_eq_of_lt ha, Nat.mod_eq_of_lt hb] at hcan'
  have himg : (Finset.range p).image (fun k => (w' * k) % p) = Finset.range p := by
    apply Finset.eq_of_subset_of_card_le
    · intro x hx
      rw [Finset.mem_image] at hx
      obtain ⟨k, hk, rfl⟩ := hx
      exact Finset.mem_range.mpr (Nat.mod_lt _ hp)
    · rw [Finset.card_image_of_injOn hinj]
  have hkey : ((Finset.range p).filter (fun k => Q ((w' * k) % p))).image
      (fun k => (w' * k) % p) = (Finset.range p).filter Q := by
    ext y
    simp only [Finset.mem_image, Finset.mem_filter, Finset.mem_range]
    constructor
    · rintro ⟨k, ⟨hk, hQk⟩, rfl⟩
      exact ⟨Nat.mod_lt _ hp, hQk⟩
    · rintro ⟨hy, hQy⟩
      have hy' : y ∈ (Finset.range p).image (fun k => (w' * k) % p) := by
        rw [himg]
        exact Finset.mem_range.mpr hy
      rw [Finset.mem_image] at hy'
      obtain ⟨k, hk, hky⟩ := hy'
      refine ⟨k, ⟨Finset.mem_range.mp hk, ?_⟩, hky⟩
      rw [hky]
      exact hQy
  have hfilter_inj : Set.InjOn (fun k => (w' * k) % p)
      ↑((Finset.range p).filter (fun k => Q ((w' * k) % p))) := by
    apply hinj.mono
    intro x hx
    simp only [Finset.coe_filter, Set.mem_setOf_eq] at hx
    simp only [Finset.coe_range, Set.mem_Iio]
    exact Finset.mem_range.mp hx.1
  calc ((Finset.range p).filter (fun k => Q ((w' * k) % p))).card
      = (((Finset.range p).filter (fun k => Q ((w' * k) % p))).image
          (fun k => (w' * k) % p)).card :=
        (Finset.card_image_of_injOn hfilter_inj).symm
    _ = ((Finset.range p).filter Q).card := by rw [hkey]

/-! ## (c) the wrapped-arc count -/

/-- Grid points of `range p` in the wrapped arc `∃ m, |X + ρ/p − m| < 2/25`:
at most `(4/25)p + 3`. -/
lemma arc_grid_count (p : ℕ) (hp : 0 < p) (X : ℝ) :
    (((Finset.range p).filter
        (fun ρ : ℕ => ∃ m : ℤ, |X + (ρ : ℝ) / p - m| < 2/25)).card : ℝ)
      ≤ (4/25) * p + 3 := by
  have hpR : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  set M₀ : ℤ := ⌈X - 2/25⌉ with hM₀
  -- every witness integer is M₀ or M₀ + 1
  have hwitness : ∀ ρ : ℕ, ρ ∈ Finset.range p →
      (∃ m : ℤ, |X + (ρ : ℝ) / p - m| < 2/25) →
      (∃ m : ℤ, (m = M₀ ∨ m = M₀ + 1) ∧ |X + (ρ : ℝ) / p - m| < 2/25) := by
    intro ρ hρ ⟨m, hm⟩
    refine ⟨m, ?_, hm⟩
    have hρ0 : (0 : ℝ) ≤ (ρ : ℝ) / p := by positivity
    have hρ1 : (ρ : ℝ) / p < 1 := by
      rw [div_lt_one hpR]
      exact_mod_cast Finset.mem_range.mp hρ
    have habs := abs_lt.mp hm
    -- X − 2/25 < m < X + 1 + 2/25
    have hlo : X - 2/25 < (m : ℝ) := by linarith
    have hhi : (m : ℝ) < X + 1 + 2/25 := by linarith
    have hM₀le : M₀ ≤ m := by
      rw [hM₀]
      exact Int.ceil_le.mpr (le_of_lt hlo)
    have hM₀ge : (X - 2/25 : ℝ) ≤ (M₀ : ℝ) := by
      rw [hM₀]
      exact Int.le_ceil _
    have hmM₀ : (m : ℝ) - (M₀ : ℝ) < 2 := by linarith
    have : m - M₀ < 2 := by exact_mod_cast hmM₀
    omega
  -- split by the two witnesses
  have hsub : (Finset.range p).filter
      (fun ρ : ℕ => ∃ m : ℤ, |X + (ρ : ℝ) / p - m| < 2/25)
      ⊆ ((Finset.range p).filter
          (fun ρ : ℕ => |X + (ρ : ℝ) / p - M₀| < 2/25))
        ∪ ((Finset.range p).filter
          (fun ρ : ℕ => |X + (ρ : ℝ) / p - (M₀ + 1 : ℤ)| < 2/25)) := by
    intro ρ hρ
    rw [Finset.mem_filter] at hρ
    obtain ⟨hρr, hex⟩ := hρ
    obtain ⟨m, hm01, hm⟩ := hwitness ρ hρr hex
    rw [Finset.mem_union, Finset.mem_filter, Finset.mem_filter]
    rcases hm01 with rfl | rfl
    · exact Or.inl ⟨hρr, hm⟩
    · exact Or.inr ⟨hρr, by exact_mod_cast hm⟩
  -- the two pieces: ρ-intervals (A, B) and (A+p, B+p), clipped to [0, p)
  set A : ℝ := ((M₀ : ℝ) - X - 2/25) * p with hA
  set B : ℝ := ((M₀ : ℝ) - X + 2/25) * p with hB
  -- structural bounds on A, B from the ceiling
  have hAlow : -(4/25) * p ≤ A := by
    rw [hA]
    have h1 : (X - 2/25 : ℝ) ≤ (M₀ : ℝ) := by rw [hM₀]; exact Int.le_ceil _
    nlinarith [hpR]
  have hBp : B < p := by
    rw [hB]
    have h2 : (M₀ : ℝ) < X - 2/25 + 1 := by
      rw [hM₀]
      exact Int.ceil_lt_add_one _
    nlinarith [hpR]
  have hBA : B - A = (4/25) * p := by rw [hA, hB]; ring
  -- piece 1: interval (max A (−1), B)
  have hc1 : (((Finset.range p).filter
      (fun ρ : ℕ => |X + (ρ : ℝ) / p - M₀| < 2/25)).card : ℝ)
      ≤ B - max A (-1) + 1 := by
    apply range_filter_interval_card p (max A (-1)) B
      (by
        rcases le_or_gt A (-1) with h | h
        · rw [max_eq_right h]; linarith [hAlow, hBA]
        · rw [max_eq_left (le_of_lt h)]; linarith [hBA, hpR])
    intro ρ hρ hP
    have habs := abs_lt.mp hP
    have hup : (ρ : ℝ) < B := by
      rw [hB]
      have h1 : (ρ : ℝ) / p < (M₀ : ℝ) - X + 2/25 := by linarith
      have h2 := (div_lt_iff₀ hpR).mp h1
      linarith
    have hdown : A < (ρ : ℝ) := by
      rw [hA]
      have h1 : (M₀ : ℝ) - X - 2/25 < (ρ : ℝ) / p := by linarith
      have h2 := (lt_div_iff₀ hpR).mp h1
      linarith
    have hneg1 : (-1 : ℝ) < (ρ : ℝ) := by
      have h0 : (0 : ℝ) ≤ (ρ : ℝ) := Nat.cast_nonneg ρ
      linarith
    exact ⟨max_lt_iff.mpr ⟨hdown, hneg1⟩, hup⟩
  -- piece 2: the (M₀+1)-interval is the p-translate (A+p, B+p); inside [0, p)
  -- it is (A+p, p): empty when 0 ≤ A, and of length −A when A < 0.
  have hc2 : (((Finset.range p).filter
      (fun ρ : ℕ => |X + (ρ : ℝ) / p - (M₀ + 1 : ℤ)| < 2/25)).card : ℝ)
      ≤ max (1 - A) 0 := by
    rcases le_or_gt 0 A with hApos | hAneg
    · -- the translate misses [0, p) entirely
      have hempty : (Finset.range p).filter
          (fun ρ : ℕ => |X + (ρ : ℝ) / p - (M₀ + 1 : ℤ)| < 2/25) = ∅ := by
        rw [Finset.filter_eq_empty_iff]
        intro ρ hρ hP
        have habs := abs_lt.mp hP
        have hρ1 : (ρ : ℝ) / p < 1 := by
          rw [div_lt_one hpR]
          exact_mod_cast Finset.mem_range.mp hρ
        have hcoef : 0 ≤ (M₀ : ℝ) - X - 2/25 := by
          have h := div_nonneg hApos hpR.le
          rw [hA] at h
          rwa [mul_div_cancel_right₀ _ hpR.ne'] at h
        have hlow : ((M₀ : ℝ) + 1) - X - 2/25 < (ρ : ℝ) / p := by
          have h1 := habs.1
          push_cast at h1
          linarith
        linarith
      rw [hempty]
      simp only [Finset.card_empty, Nat.cast_zero]
      exact le_max_right _ _
    · -- interval (A + p, p), length −A; count ≤ 1 − A
      have hbound : (((Finset.range p).filter
          (fun ρ : ℕ => |X + (ρ : ℝ) / p - (M₀ + 1 : ℤ)| < 2/25)).card : ℝ)
          ≤ (p : ℝ) - (A + p) + 1 := by
        apply range_filter_interval_card p (A + p) (p : ℝ) (by linarith)
        intro ρ hρ hP
        have habs := abs_lt.mp hP
        constructor
        · have h1 : ((M₀ : ℝ) + 1) - X - 2/25 < (ρ : ℝ) / p := by
            have h1' := habs.1
            push_cast at h1'
            linarith
          have h2 := (lt_div_iff₀ hpR).mp h1
          have h3 : (((M₀ : ℝ) + 1) - X - 2/25) * p = A + p := by rw [hA]; ring
          linarith
        · exact_mod_cast Finset.mem_range.mp hρ
      have hlen : (p : ℝ) - (A + p) + 1 = 1 - A := by ring
      rw [hlen] at hbound
      exact hbound.trans (le_max_left _ _)
  -- assemble the two pieces, then clip
  have hle := Finset.card_le_card hsub
  have hu := Finset.card_union_le
    ((Finset.range p).filter (fun ρ : ℕ => |X + (ρ : ℝ) / p - M₀| < 2/25))
    ((Finset.range p).filter (fun ρ : ℕ => |X + (ρ : ℝ) / p - (M₀ + 1 : ℤ)| < 2/25))
  have h1 : (((Finset.range p).filter
      (fun ρ : ℕ => ∃ m : ℤ, |X + (ρ : ℝ) / p - m| < 2/25)).card : ℝ)
      ≤ (((Finset.range p).filter
          (fun ρ : ℕ => |X + (ρ : ℝ) / p - M₀| < 2/25)).card : ℝ)
        + (((Finset.range p).filter
          (fun ρ : ℕ => |X + (ρ : ℝ) / p - (M₀ + 1 : ℤ)| < 2/25)).card : ℝ) := by
    exact_mod_cast le_trans hle hu
  have hclip : (B - max A (-1) + 1) + max (1 - A) 0 ≤ (4/25) * p + 3 := by
    rcases le_or_gt 0 A with h0 | h0
    · rw [max_eq_left (by linarith : (-1:ℝ) ≤ A)]
      have h1' : max (1 - A) 0 ≤ 1 := max_le (by linarith) (by norm_num)
      linarith [hBA]
    · have h1' : -1 ≤ max A (-1) := le_max_right A (-1)
      have h2' : max (1 - A) 0 = 1 - A := max_eq_left (by linarith)
      rw [h2']
      linarith [hBA]
  linarith [h1, hc1, hc2, hclip]

/-! ## the sharp visit count -/

/-- **THE SHARP VISIT COUNT**: among `D` equally-spaced copies, a comb `w`
hosts at most `(4/25)·D + 3·w`. -/
lemma tooth_visit_count_sharp (w : ℤ) (hw : 0 < w) (D : ℕ) (hD : 0 < D) (c : ℝ) :
    (((Finset.range D).filter
        (fun j : ℕ => ∃ m : ℤ, |(w : ℝ) * (c + (j : ℝ) / D) - m| < 2/25)).card : ℝ)
      ≤ (4/25) * D + 3 * w := by
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
  -- the value identity: w(c + j/D) = Wc + (w′j)/p
  have hval : ∀ j : ℕ, (w : ℝ) * (c + (j : ℝ) / D)
      = (W : ℝ) * c + ((w' : ℝ) * j) / p := by
    intro j
    rw [← hWR]
    have hkey : ((W : ℕ) : ℝ) * ((j : ℝ) / D) = ((w' : ℝ) * j) / p := by
      rw [eq_div_iff hpRl.ne']
      calc (W : ℝ) * ((j : ℝ) / D) * p = ((j : ℝ) / D) * ((W : ℝ) * p) := by ring
        _ = ((j : ℝ) / D) * ((w' : ℝ) * D) := by rw [hWpD]
        _ = ((w' : ℝ) * j) * ((D : ℝ) / D) := by ring
        _ = ((w' : ℝ) * j) * 1 := by rw [div_self hDRl.ne']
        _ = (w' : ℝ) * j := by ring
    calc ((W : ℕ) : ℝ) * (c + (j : ℝ) / D)
        = (W : ℝ) * c + (W : ℝ) * ((j : ℝ) / D) := by ring
      _ = (W : ℝ) * c + ((w' : ℝ) * j) / p := by rw [hkey]
  -- split the fraction at the mod: (w′j)/p = ((w′j) mod p)/p + (w′j) div p
  have hmoddiv : ∀ j : ℕ, ((w' : ℝ) * j) / p
      = (((w' * j) % p : ℕ) : ℝ) / p + (((w' * j) / p : ℕ) : ℝ) := by
    intro j
    have hde : (w' * j) % p + p * ((w' * j) / p) = w' * j := Nat.mod_add_div (w' * j) p
    have hdeR : (((w' * j) % p : ℕ) : ℝ) + (p : ℝ) * (((w' * j) / p : ℕ) : ℝ)
        = (w' : ℝ) * j := by exact_mod_cast hde
    rw [← hdeR, add_div, mul_div_cancel_left₀ _ hpRl.ne']
  -- the visit condition depends on j only through (w′j) mod p
  have hcond : ∀ j : ℕ, (∃ m : ℤ, |(w : ℝ) * (c + (j : ℝ) / D) - m| < 2/25)
      ↔ (∃ m : ℤ, |(W : ℝ) * c + (((w' * j) % p : ℕ) : ℝ) / p - m| < 2/25) := by
    intro j
    set q : ℕ := (w' * j) / p with hq
    constructor
    · rintro ⟨m, hm⟩
      refine ⟨m - q, ?_⟩
      rw [hval j, hmoddiv j] at hm
      have hswap : (W : ℝ) * c + (((w' * j) % p : ℕ) : ℝ) / p - ((m - q : ℤ) : ℝ)
          = (W : ℝ) * c + ((((w' * j) % p : ℕ) : ℝ) / p + (q : ℝ)) - (m : ℝ) := by
        push_cast
        ring
      rw [hswap]
      exact hm
    · rintro ⟨m, hm⟩
      refine ⟨m + q, ?_⟩
      rw [hval j, hmoddiv j]
      have hswap : (W : ℝ) * c + ((((w' * j) % p : ℕ) : ℝ) / p + (q : ℝ))
            - ((m + q : ℤ) : ℝ)
          = (W : ℝ) * c + (((w' * j) % p : ℕ) : ℝ) / p - (m : ℝ) := by
        push_cast
        ring
      rw [hswap]
      exact hm
  have hfe : (Finset.range D).filter
      (fun j : ℕ => ∃ m : ℤ, |(w : ℝ) * (c + (j : ℝ) / D) - m| < 2/25)
      = (Finset.range D).filter
        (fun j : ℕ => ∃ m : ℤ, |(W : ℝ) * c + (((w' * j) % p : ℕ) : ℝ) / p - m| < 2/25) :=
    Finset.filter_congr (fun j _ => hcond j)
  rw [hfe]
  -- block split + coprime permutation
  have hper : ∀ j : ℕ,
      (∃ m : ℤ, |(W : ℝ) * c + (((w' * (j + p)) % p : ℕ) : ℝ) / p - m| < 2/25)
      ↔ (∃ m : ℤ, |(W : ℝ) * c + (((w' * j) % p : ℕ) : ℝ) / p - m| < 2/25) := by
    intro j
    have hmm : (w' * (j + p)) % p = (w' * j) % p := by
      have hexp : w' * (j + p) = w' * j + w' * p := by ring
      rw [hexp, Nat.add_mul_mod_self_right]
    rw [hmm]
  conv_lhs => rw [hDgp]
  rw [filter_range_mul_periodic
    (fun j : ℕ => ∃ m : ℤ, |(W : ℝ) * c + (((w' * j) % p : ℕ) : ℝ) / p - m| < 2/25)
    p hper g]
  rw [filter_range_coprime_mul
    (fun ρ : ℕ => ∃ m : ℤ, |(W : ℝ) * c + (ρ : ℝ) / p - m| < 2/25)
    p w' hp0 hco]
  -- the arc count, times g blocks
  have harc := arc_grid_count p hp0 ((W : ℝ) * c)
  have hgw : (g : ℝ) ≤ ((w : ℤ) : ℝ) := by
    rw [← hWR]
    exact_mod_cast Nat.le_of_dvd hW0 hgdvdW
  have hgD : (g : ℝ) * p = (D : ℝ) := by exact_mod_cast hDgp.symm
  have hg0R : (0 : ℝ) < (g : ℝ) := by exact_mod_cast hg0
  push_cast
  calc (g : ℝ) * (((Finset.range p).filter
      (fun ρ : ℕ => ∃ m : ℤ, |(W : ℝ) * c + (ρ : ℝ) / p - m| < 2/25)).card : ℝ)
      ≤ (g : ℝ) * ((4/25) * p + 3) := by
        exact mul_le_mul_of_nonneg_left harc hg0R.le
    _ = (4/25) * ((g : ℝ) * p) + 3 * g := by ring
    _ = (4/25) * D + 3 * g := by rw [hgD]
    _ ≤ (4/25) * D + 3 * w := by linarith [hgw]

/-! ## the full ladder -/

/-- **THE SHARP CLUSTER-GCD RUNG (|S| ≤ 6)**: a 12-family with no `2/25`-margin
point satisfies `(25 − 4|S|)·d ≤ 75·Σ_S |vᵢ|` for every `S` with `1 ≤ |S| ≤ 6`
and every positive common divisor `d` of the complement. -/
theorem gap_gcd_rung_sharp (cite : LRCUpTo13) (v : Fin 12 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25)
    (S : Finset (Fin 12)) (hS1 : S.Nonempty) (hS6 : S.card ≤ 6)
    (d : ℤ) (hd0 : 0 < d) (hdvd : ∀ i, i ∉ S → d ∣ v i) :
    (25 - 4 * (S.card : ℤ)) * d ≤ 75 * (∑ i ∈ S, |v i|) := by
  classical
  -- the big side T, cited
  set T : Finset (Fin 12) := Finset.univ \ S with hT
  have hTmem : ∀ i, i ∈ T ↔ i ∉ S := by
    intro i
    rw [hT, Finset.mem_sdiff]
    simp only [Finset.mem_univ, true_and]
  have hTcard : T.card ≤ 12 := by
    calc T.card ≤ Finset.univ.card := Finset.card_le_card (Finset.subset_univ T)
      _ = 12 := by rw [Finset.card_univ, Fintype.card_fin]
  obtain ⟨t₀, hmargin⟩ := cite_margin_gen cite v T hTcard (fun i _ => hv i)
  -- the citation margin exceeds 2/25 (T.card ≤ 11 since S nonempty)
  have hT11 : T.card ≤ 11 := by
    have hScard1 : 1 ≤ S.card := Finset.card_pos.mpr hS1
    have : T.card = 12 - S.card := by
      rw [hT, Finset.card_univ_diff, Fintype.card_fin]
    omega
  have hmargin' : ∀ i ∈ T, ∀ m : ℤ, (2:ℝ)/25 < |(v i : ℝ) * t₀ - m| := by
    intro i hi m
    have h := hmargin i hi m
    have hcard : ((T.card : ℝ) + 1) ≤ 12 := by
      have : (T.card : ℝ) ≤ 11 := by exact_mod_cast hT11
      linarith
    have h12 : (1:ℝ)/12 ≤ 1 / ((T.card : ℝ) + 1) := by
      apply div_le_div_of_nonneg_left (by norm_num) ?_ hcard
      positivity
    calc (2:ℝ)/25 < 1/12 := by norm_num
      _ ≤ 1 / ((T.card : ℝ) + 1) := h12
      _ ≤ |(v i : ℝ) * t₀ - m| := h
  -- every shifted copy is T-good; the violator at each copy is in S
  set D : ℕ := d.toNat with hD
  have hDpos : 0 < D := by rw [hD]; omega
  have hDd : (D : ℤ) = d := by rw [hD]; omega
  have hDdR : ((D : ℕ) : ℝ) = ((d : ℤ) : ℝ) := by exact_mod_cast hDd
  have hpick : ∀ j : ℕ, j ∈ Finset.range D →
      ∃ i ∈ S, ∃ m : ℤ, |(v i : ℝ) * (t₀ + (j : ℝ) / D) - m| < 2/25 := by
    intro j _
    obtain ⟨i, m, hm⟩ := hnl (t₀ + (j : ℝ) / D)
    by_cases hiS : i ∈ S
    · exact ⟨i, hiS, m, hm⟩
    · exfalso
      have hiT : i ∈ T := (hTmem i).mpr hiS
      have hper := periodic_margin (v i) d hd0 (hdvd i hiS) t₀ (j : ℤ)
        ((2:ℝ)/25) (fun m' => le_of_lt (hmargin' i hiT m'))
      have hh := hper m
      rw [show ((j : ℤ) : ℝ) = (j : ℝ) by push_cast; ring, ← hDdR] at hh
      linarith
  -- the classifier and the fiberwise count
  have hchoice : ∀ j : ℕ, ∃ i : Fin 12, j ∈ Finset.range D →
      (i ∈ S ∧ ∃ m : ℤ, |(v i : ℝ) * (t₀ + (j : ℝ) / D) - m| < 2/25) := by
    intro j
    by_cases hj : j ∈ Finset.range D
    · obtain ⟨i, hiS, hm⟩ := hpick j hj
      exact ⟨i, fun _ => ⟨hiS, hm⟩⟩
    · exact ⟨0, fun hc => absurd hc hj⟩
  choose F hFspec using hchoice
  have hFmem : ∀ j ∈ Finset.range D, F j ∈ S := fun j hj => (hFspec j hj).1
  have hsplit : (Finset.range D).card
      = ∑ i ∈ S, ((Finset.range D).filter (fun j => F j = i)).card :=
    Finset.card_eq_sum_card_fiberwise hFmem
  -- per-fiber: inside the sharp visit filter of |v i|
  have hfiber : ∀ i ∈ S, (((Finset.range D).filter (fun j => F j = i)).card : ℝ)
      ≤ (4/25) * D + 3 * ((|v i| : ℤ) : ℝ) := by
    intro i _
    have hsub : (Finset.range D).filter (fun j => F j = i)
        ⊆ (Finset.range D).filter
          (fun j : ℕ => ∃ m : ℤ, |((|v i| : ℤ) : ℝ) * (t₀ + (j : ℝ) / D) - m| < 2/25) := by
      intro j hj
      rw [Finset.mem_filter] at hj ⊢
      obtain ⟨hjr, hjF⟩ := hj
      have hm := (hFspec j hjr).2
      rw [hjF] at hm
      exact ⟨hjr, (cov_abs_iff (v i) (t₀ + (j : ℝ) / D)).mp hm⟩
    have hcount := tooth_visit_count_sharp (|v i|) (abs_pos.mpr (hv i)) D hDpos t₀
    calc (((Finset.range D).filter (fun j => F j = i)).card : ℝ)
        ≤ (((Finset.range D).filter
          (fun j : ℕ => ∃ m : ℤ, |((|v i| : ℤ) : ℝ) * (t₀ + (j : ℝ) / D) - m| < 2/25)).card : ℝ) := by
          exact_mod_cast Finset.card_le_card hsub
      _ ≤ (4/25) * D + 3 * ((|v i| : ℤ) : ℝ) := hcount
  -- assemble over S
  have hsum : ((D : ℕ) : ℝ)
      ≤ (4/25) * D * S.card + 3 * (∑ i ∈ S, ((|v i| : ℤ) : ℝ)) := by
    have hcards : ((D : ℕ) : ℝ)
        = ∑ i ∈ S, (((Finset.range D).filter (fun j => F j = i)).card : ℝ) := by
      have h0 : ((Finset.range D).card : ℝ)
          = ∑ i ∈ S, (((Finset.range D).filter (fun j => F j = i)).card : ℝ) := by
        exact_mod_cast congrArg (Nat.cast : ℕ → ℝ) hsplit
      rw [← h0, Finset.card_range]
    conv_lhs => rw [hcards]
    calc (∑ i ∈ S, (((Finset.range D).filter (fun j => F j = i)).card : ℝ))
        ≤ ∑ i ∈ S, ((4/25) * (D : ℝ) + 3 * ((|v i| : ℤ) : ℝ)) :=
          Finset.sum_le_sum hfiber
      _ = (4/25) * D * S.card + 3 * (∑ i ∈ S, ((|v i| : ℤ) : ℝ)) := by
          rw [Finset.sum_add_distrib, Finset.sum_const, ← Finset.mul_sum]
          push_cast
          ring
  -- cast down to the integer conclusion
  have hcast : (((25 - 4 * (S.card : ℤ)) * d : ℤ) : ℝ)
      ≤ ((75 * (∑ i ∈ S, |v i|) : ℤ) : ℝ) := by
    have hdD : ((d : ℤ) : ℝ) = ((D : ℕ) : ℝ) := hDdR.symm
    push_cast
    push_cast at hsum
    rw [hdD]
    nlinarith [hsum]
  exact_mod_cast hcast

#print axioms tooth_visit_count_sharp
#print axioms gap_gcd_rung_sharp

end ClusterGcdSharp
end LonelyRunner
