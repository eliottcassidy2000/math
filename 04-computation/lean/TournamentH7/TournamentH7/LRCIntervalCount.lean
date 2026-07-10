/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-10-S211)
-/
import Mathlib

/-!
# Counting integers in a real interval (opus-S211)

The foundational counting lemma for THM-678's `d = 2` orbit count: a real open interval of width `W`
contains at most `⌊W⌋ + 1` integers. This is the de-circled crux — the two-detuned "bad branch" count
`|badⱼ| ≤ dⱼ·(⌊qⱼ/7⌋+1)` reduces to it via the injection `c ↦ δⱼc − g·round(δⱼ(u+c)/g)`.

Kernel-pure: no `sorry`, no `native_decide`.
-/

namespace LRCIntervalCount

open Finset
open scoped Classical

/-- **Integers in a real interval.** Any finset of integers all lying strictly inside a real interval
`(x, x + W)` (with `W ≥ 0`) has at most `⌊W⌋ + 1` elements. -/
theorem card_le_of_mem_Ioo (S : Finset ℤ) (x W : ℝ) (hW : 0 ≤ W)
    (hS : ∀ m ∈ S, x < (m : ℝ) ∧ (m : ℝ) < x + W) :
    S.card ≤ ⌊W⌋.toNat + 1 := by
  -- every `m ∈ S` sits in the integer window `Finset.Ioo ⌊x⌋ ⌈x + W⌉`
  have hsub : S ⊆ Finset.Ioo ⌊x⌋ ⌈x + W⌉ := by
    intro m hm
    obtain ⟨hlo, hhi⟩ := hS m hm
    rw [Finset.mem_Ioo]
    constructor
    · -- ⌊x⌋ < m : since x < m and m is an integer, ⌊x⌋ ≤ x < m, and ⌊x⌋ < m as integers
      have : (⌊x⌋ : ℝ) ≤ x := Int.floor_le x
      exact_mod_cast lt_of_le_of_lt this hlo
    · -- m < ⌈x + W⌉ : since m < x + W ≤ ⌈x+W⌉
      have : (m : ℝ) < x + W := hhi
      have hc : x + W ≤ (⌈x + W⌉ : ℝ) := Int.le_ceil _
      exact_mod_cast lt_of_lt_of_le this hc
  refine le_trans (Finset.card_le_card hsub) ?_
  rw [Int.card_Ioo]
  -- `(⌈x+W⌉ - ⌊x⌋ - 1).toNat ≤ ⌊W⌋.toNat + 1`
  have hceil : (⌈x + W⌉ : ℤ) ≤ ⌈x⌉ + ⌈(W : ℝ)⌉ := Int.ceil_add_le x W
  have hcx : (⌈x⌉ : ℤ) ≤ ⌊x⌋ + 1 := Int.ceil_le_floor_add_one x
  have hcw : (⌈(W : ℝ)⌉ : ℤ) ≤ ⌊W⌋ + 1 := Int.ceil_le_floor_add_one W
  have hkey : (⌈x + W⌉ : ℤ) - ⌊x⌋ - 1 ≤ ⌊W⌋ + 1 := by omega
  have hWnn : (0 : ℤ) ≤ ⌊W⌋ := Int.le_floor.mpr (by exact_mod_cast hW)
  omega

/-- **The two-detuned bad-branch count.** For any `δ, g` (`g ≥ 1`) and real `u`, the number of branches
`c ∈ [0, g)` at which `δ(u+c)/g` fails to clear `1/14` is at most `d·(⌊q/7⌋ + 1)`, where `d = gcd(δ,g)`,
`q = g/d`. Proof: the injection `c ↦ p·c − q·round(δ(u+c)/g)` (`p = δ/d`) sends each bad branch into a real
interval of width `q/7`, and is `≤ d`-to-1 (`gcd(p,q) = 1 ⟹ q ∣ c − c'`). -/
theorem bad_count_le (δ g : ℤ) (hg : 1 ≤ g) (u : ℝ) :
    ((Finset.Ico (0 : ℤ) g).filter
        (fun (c : ℤ) => ∃ n : ℤ, |(δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n| < 1 / 14)).card
      ≤ (Int.gcd δ g) * ((g / (Int.gcd δ g : ℤ)).toNat / 7 + 1) := by
  have hg0 : (0 : ℤ) < g := hg
  have hgR : (0 : ℝ) < (g : ℝ) := by exact_mod_cast hg0
  set d : ℤ := (Int.gcd δ g : ℤ) with hd
  have hdpos : 0 < d := by
    rw [hd]
    have : Int.gcd δ g ≠ 0 := by
      intro h; rw [Int.gcd_eq_zero_iff] at h; omega
    exact_mod_cast Nat.pos_of_ne_zero this
  have hd0 : d ≠ 0 := ne_of_gt hdpos
  have hddδ : d ∣ δ := Int.gcd_dvd_left δ g
  have hddg : d ∣ g := Int.gcd_dvd_right δ g
  set p : ℤ := δ / d with hp
  set q : ℤ := g / d with hq
  have hδdp : δ = d * p := (Int.mul_ediv_cancel' hddδ).symm
  have hgdq : g = d * q := (Int.mul_ediv_cancel' hddg).symm
  have hqpos : 0 < q := by
    rcases lt_or_ge 0 q with h | h
    · exact h
    · exfalso
      have hle : d * q ≤ 0 := mul_nonpos_iff.mpr (Or.inl ⟨le_of_lt hdpos, h⟩)
      rw [← hgdq] at hle; omega
  -- coprimality of `p, q` via Bézout `d = δ·A + g·B`
  have hcop : IsCoprime p q := by
    set A := Int.gcdA δ g with hA
    set B := Int.gcdB δ g with hB
    have hbez : (d : ℤ) = δ * A + g * B := by rw [hd, hA, hB]; exact Int.gcd_eq_gcd_ab δ g
    refine ⟨A, B, ?_⟩
    have hfac : d * (A * p + B * q) = d * 1 := by
      have h1 : d * (A * p + B * q) = (d * p) * A + (d * q) * B := by ring
      rw [h1, ← hδdp, ← hgdq, ← hbez, mul_one]
    exact mul_left_cancel₀ hd0 hfac
  -- the bad set and the injection ψ
  set S := (Finset.Ico (0 : ℤ) g).filter
      (fun (c : ℤ) => ∃ n : ℤ, |(δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n| < 1 / 14) with hS
  set ψ : ℤ → ℤ := fun c => p * c - q * ⌊(δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) + 1 / 2⌋ with hψ
  -- for each bad `c`, the rounding hits the witness `n`, so `|δ(u+c)/g − round| < 1/14`
  have hround : ∀ c ∈ S, ∃ n : ℤ, ⌊(δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) + 1 / 2⌋ = n ∧
      |(δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n| < 1 / 14 := by
    intro c hc
    rw [hS, Finset.mem_filter] at hc
    obtain ⟨n, hn⟩ := hc.2
    refine ⟨n, ?_, hn⟩
    rw [abs_lt] at hn
    rw [Int.floor_eq_iff]
    push_cast; constructor <;> linarith [hn.1, hn.2]
  -- `ψ c` lands in a width-`q/7` real interval
  have himg : ∀ c ∈ S, (ψ c : ℝ) ∈ Set.Ioo (-(p : ℝ) * u - (q : ℝ) / 14)
      ((-(p : ℝ) * u - (q : ℝ) / 14) + (q : ℝ) / 7) := by
    intro c hc
    obtain ⟨n, hfl, hlt⟩ := hround c hc
    have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hqpos
    have hval : (ψ c : ℝ) = (q : ℝ) * ((δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n) - (p : ℝ) * u := by
      rw [hψ]; simp only [hfl]
      have hgR' : (g : ℝ) ≠ 0 := ne_of_gt hgR
      have hgdqR : (g : ℝ) = (d : ℝ) * (q : ℝ) := by exact_mod_cast hgdq
      have hδdpR : (δ : ℝ) = (d : ℝ) * (p : ℝ) := by exact_mod_cast hδdp
      have hdR : (0 : ℝ) < (d : ℝ) := by exact_mod_cast hdpos
      push_cast
      rw [hgdqR, hδdpR]; field_simp; ring
    rw [abs_lt] at hlt
    rw [Set.mem_Ioo, hval]
    constructor <;> nlinarith [hlt.1, hlt.2, hqR]
  -- `ψ` is `≤ d`-to-1 on `S`:  the maps-to and the injectivity of `c ↦ (c/q).toNat`
  have hmaps : ∀ a : ℤ, ∀ c ∈ S.filter (fun c => ψ c = a),
      (fun c => (c / q).toNat) c ∈ Finset.range d.toNat := by
    intro a c hc
    rw [Finset.mem_filter, hS, Finset.mem_filter] at hc
    obtain ⟨⟨hcIco, _⟩, _⟩ := hc
    rw [Finset.mem_Ico] at hcIco
    rw [Finset.mem_range]
    have h1 : c / q < d := by
      rw [Int.ediv_lt_iff_lt_mul hqpos]
      have h := hcIco.2; rw [hgdq] at h; exact h
    have h2 : 0 ≤ c / q := Int.ediv_nonneg hcIco.1 (le_of_lt hqpos)
    show (c / q).toNat < d.toNat
    omega
  have hinj : ∀ a : ℤ, Set.InjOn (fun c => (c / q).toNat) ↑(S.filter (fun c => ψ c = a)) := by
    intro a c hc c' hc' heq
    simp only [Finset.mem_coe, Finset.mem_filter] at hc hc'
    have hqdvd : q ∣ (c - c') := by
      have hpc : p * c - q * ⌊(δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) + 1 / 2⌋
          = p * c' - q * ⌊(δ : ℝ) * ((u + (c' : ℝ)) / (g : ℝ)) + 1 / 2⌋ := hc.2.trans hc'.2.symm
      have hqp : q ∣ p * (c - c') :=
        ⟨⌊(δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) + 1 / 2⌋
          - ⌊(δ : ℝ) * ((u + (c' : ℝ)) / (g : ℝ)) + 1 / 2⌋, by linear_combination hpc⟩
      exact hcop.symm.dvd_of_dvd_mul_left hqp
    obtain ⟨hcS, _⟩ := hc; obtain ⟨hc'S, _⟩ := hc'
    rw [hS, Finset.mem_filter, Finset.mem_Ico] at hcS hc'S
    have hcq0 : 0 ≤ c / q := Int.ediv_nonneg hcS.1.1 (le_of_lt hqpos)
    have hc'q0 : 0 ≤ c' / q := Int.ediv_nonneg hc'S.1.1 (le_of_lt hqpos)
    have hdiv : c / q = c' / q := by
      simp only [] at heq; omega
    have hmod : c % q = c' % q := Int.modEq_iff_dvd.mpr (dvd_sub_comm.mp hqdvd)
    have e1 := Int.emod_add_ediv c q
    have e2 := Int.emod_add_ediv c' q
    rw [hdiv, hmod] at e1
    exact e1.symm.trans e2
  have hfib : ∀ a ∈ S.image ψ, (S.filter (fun c => ψ c = a)).card ≤ d.toNat := by
    intro a _
    rw [← Finset.card_range d.toNat]
    exact Finset.card_le_card_of_injOn _ (hmaps a) (hinj a)
  have h1 : S.card ≤ d.toNat * (S.image ψ).card := by
    apply Finset.card_le_mul_card_image; exact hfib
  -- the image lands in the width-`q/7` interval, so `|image| ≤ ⌊q/7⌋ + 1 = q.toNat/7 + 1`
  have hfloorR : ⌊(q : ℝ) / 7⌋ = q / 7 := by
    have h7 : (0 : ℝ) < 7 := by norm_num
    have e := Int.ediv_add_emod q 7
    have m1 := Int.emod_nonneg q (by norm_num : (7 : ℤ) ≠ 0)
    have m2 := Int.emod_lt_of_pos q (by norm_num : (0 : ℤ) < 7)
    rw [Int.floor_eq_iff]
    refine ⟨?_, ?_⟩
    · rw [le_div_iff₀ h7]
      have hle : ((q / 7) * 7 : ℤ) ≤ q := by omega
      calc ((q / 7 : ℤ) : ℝ) * 7 = (((q / 7) * 7 : ℤ) : ℝ) := by push_cast; ring
        _ ≤ (q : ℝ) := by exact_mod_cast hle
    · rw [div_lt_iff₀ h7]
      have hlt : (q : ℤ) < ((q / 7) + 1) * 7 := by omega
      have : (q : ℝ) < (((q / 7) + 1) * 7 : ℤ) := by exact_mod_cast hlt
      push_cast at this ⊢; linarith
  have hconv : (⌊(q : ℝ) / 7⌋).toNat = q.toNat / 7 := by rw [hfloorR]; omega
  have h2 : (S.image ψ).card ≤ q.toNat / 7 + 1 := by
    rw [← hconv]
    apply card_le_of_mem_Ioo (S.image ψ) (-(p : ℝ) * u - (q : ℝ) / 14) ((q : ℝ) / 7) (by positivity)
    intro m hm
    rw [Finset.mem_image] at hm
    obtain ⟨c, hc, rfl⟩ := hm
    have hmem := himg c hc
    rw [Set.mem_Ioo] at hmem
    exact hmem
  -- assemble
  have hdtoNat : d.toNat = Int.gcd δ g := by rw [hd]; exact Int.toNat_natCast _
  calc S.card ≤ d.toNat * (S.image ψ).card := h1
    _ ≤ d.toNat * (q.toNat / 7 + 1) := Nat.mul_le_mul_left _ h2
    _ = (Int.gcd δ g) * ((g / (Int.gcd δ g : ℤ)).toNat / 7 + 1) := by
        rw [hdtoNat, hq, hd]

end LRCIntervalCount
