/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-09-S196)
-/
import TournamentH7.LRCFreimanAP
import TournamentH7.LRCFreimanBurden

/-!
# The Finset bridge for the Freiman AP step (opus-S196)

`LRCFreimanAP.ap_of_min_burden` (opus-S195) is stated for an indexed `a : ℕ → ℤ`. This file wraps it for
a `Finset ℤ` so it is directly citable by THM-675 / THM-681: a set `s` of `n = |s| ≥ 5` integers with the
MINIMAL restricted sumset `|A +̂ A| = 2n − 3` is *literally* an arithmetic progression `{c, c+d, …,
c+(n−1)d}`.

The bridge is the sorted enumeration `enum s : ℕ → ℤ` — `s.orderEmbOfFin` on `[0, n)`, extended
`StrictMono`ly to all of `ℕ` by stepping `+1` above `∑|x|` (which dominates every element of `s`). We show
`enum s` is `StrictMono`, `Rset (enum s) n = restrictedSum s`, and then transport `ap_of_min_burden`.

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LRCFreimanAP

open Finset

/-- The sorted enumeration of `s`, extended `StrictMono`ly to all of `ℕ`: `s.orderEmbOfFin` on `[0,|s|)`,
and `(∑|x|) + 1 + i` above (a strictly increasing tail dominating every element of `s`). -/
noncomputable def enum (s : Finset ℤ) (i : ℕ) : ℤ :=
  if h : i < s.card then s.orderEmbOfFin rfl ⟨i, h⟩ else (∑ x ∈ s, |x|) + 1 + (i : ℤ)

theorem enum_lt {s : Finset ℤ} {i : ℕ} (h : i < s.card) :
    enum s i = s.orderEmbOfFin rfl ⟨i, h⟩ := dif_pos h

theorem enum_mem {s : Finset ℤ} {i : ℕ} (h : i < s.card) : enum s i ∈ s := by
  rw [enum_lt h]; exact s.orderEmbOfFin_mem rfl ⟨i, h⟩

theorem enum_le_sum {s : Finset ℤ} {i : ℕ} (h : i < s.card) : enum s i ≤ ∑ x ∈ s, |x| :=
  le_trans (le_abs_self _) (Finset.single_le_sum (fun y _ => abs_nonneg y) (enum_mem h))

/-- The extended sorted enumeration is strictly monotone on all of `ℕ`. -/
theorem enum_strictMono (s : Finset ℤ) : StrictMono (enum s) := by
  intro i j hij
  by_cases hj : j < s.card
  · have hi : i < s.card := by omega
    rw [enum_lt hi, enum_lt hj]
    exact (s.orderEmbOfFin rfl).strictMono (Fin.mk_lt_mk.mpr hij)
  · have hje : enum s j = (∑ x ∈ s, |x|) + 1 + (j : ℤ) := dif_neg hj
    by_cases hi : i < s.card
    · have h1 := enum_le_sum hi
      have h2 : (0 : ℤ) ≤ (j : ℤ) := Int.natCast_nonneg j
      rw [hje]; linarith
    · have hie : enum s i = (∑ x ∈ s, |x|) + 1 + (i : ℤ) := dif_neg hi
      have hlt : (i : ℤ) < (j : ℤ) := by exact_mod_cast hij
      rw [hie, hje]; linarith

/-- Every element of `s` is enumerated: surjectivity of `enum` onto `s` over `[0, |s|)`. -/
theorem enum_surj {s : Finset ℤ} {x : ℤ} (hx : x ∈ s) : ∃ i, i < s.card ∧ enum s i = x := by
  have hr : x ∈ Set.range (s.orderEmbOfFin rfl) := by
    rw [Finset.range_orderEmbOfFin]; exact hx
  obtain ⟨k, hk⟩ := hr
  exact ⟨k.val, k.isLt, by rw [enum_lt k.isLt]; exact hk⟩

/-- **The sumset bridge:** the indexed restricted sumset of the enumeration equals the finset restricted
sumset. -/
theorem Rset_enum_eq (s : Finset ℤ) : Rset (enum s) s.card = LRCFreiman.restrictedSum s := by
  ext z
  rw [mem_Rset, LRCFreiman.mem_restrictedSum]
  constructor
  · rintro ⟨i, j, hi, hj, hlt, rfl⟩
    exact ⟨enum s i, enum_mem hi, enum s j, enum_mem hj, enum_strictMono s hlt, rfl⟩
  · rintro ⟨x, hx, y, hy, hxy, rfl⟩
    obtain ⟨i, hi, rfl⟩ := enum_surj hx
    obtain ⟨j, hj, rfl⟩ := enum_surj hy
    exact ⟨i, j, hi, hj, (enum_strictMono s).lt_iff_lt.mp hxy, rfl⟩

/-- **The Finset Freiman AP step (`n ≥ 5`).** A set of `n = |s| ≥ 5` integers whose restricted sumset has
the MINIMAL size `2n − 3` is an arithmetic progression: `s = {c + k·d : 0 ≤ k < n}` for some `c` and
`d > 0`. (False for `n ≤ 4`; MISTAKE-133.) This is the near-AP characterization THM-675 / THM-681 apply to
majority-parity 7-classes, now directly citable at the `Finset` level. -/
theorem finset_min_burden_isAP (s : Finset ℤ) (hn : 5 ≤ s.card)
    (hcard : (LRCFreiman.restrictedSum s).card = 2 * s.card - 3) :
    ∃ c d : ℤ, 0 < d ∧ s = (range s.card).image (fun k : ℕ => c + (k : ℤ) * d) := by
  have hcard' : (Rset (enum s) s.card).card = 2 * s.card - 3 := by
    rw [Rset_enum_eq]; exact hcard
  have hap := ap_of_min_burden (enum_strictMono s) hn hcard'
  refine ⟨enum s 0, enum s 1 - enum s 0, by have := enum_strictMono s (show (0:ℕ) < 1 by norm_num); omega,
    ?_⟩
  have hval : ∀ i, i < s.card → enum s i = enum s 0 + (i : ℤ) * (enum s 1 - enum s 0) := by
    intro i
    induction i with
    | zero => intro _; simp
    | succ k IH =>
      intro hi
      have hk : k < s.card := by omega
      have hstep := hap k (by omega)
      have hIHk := IH hk
      have hexp : enum s (k + 1) = enum s k + (enum s 1 - enum s 0) := by linarith
      rw [hexp, hIHk]; push_cast; ring
  ext x
  simp only [mem_image, mem_range]
  constructor
  · intro hx
    obtain ⟨i, hi, hix⟩ := enum_surj hx
    exact ⟨i, hi, by rw [← hval i hi]; exact hix⟩
  · rintro ⟨k, hk, rfl⟩
    rw [← hval k hk]; exact enum_mem hk

end LRCFreimanAP
