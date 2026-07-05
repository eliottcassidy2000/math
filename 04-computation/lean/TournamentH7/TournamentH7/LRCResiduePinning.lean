/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-05-S75)
-/
import Mathlib.Tactic.NormNum.GCD
import Mathlib.Tactic.NormNum.DivMod
import TournamentH7.LonelyRunnerMathlib

/-!
# Residue pinning at the prime 13 (hdich, Half 1 — HYP-4097)

The n = 12 rigidity leaf (`hdich`) decomposes as residue pinning + lift rigidity
(opus-S74).  This file formalizes the pinning half: a 13-tight-from-above 12-family with no
multiple of 13 has residue image EXACTLY the full nonzero system `{1, …, 12}` — at the PRIME
modulus every nonzero residue is a unit, so the unit-residue lemma (`exists_residue_one_of_tight`,
kernel-pure since S34) forces all twelve classes, and a cardinality pigeonhole finishes.

Downstream (klein HYP-4096 assembly): `hdich` then only needs lift rigidity — the
sieve / citation-window / finite-CRT trichotomy of the S74 write-up.
-/

namespace LonelyRunner
namespace ResiduePinning

open Finset

/-- Transport of the unit-residue lemma to a NAMED residue: if `u * a ≡ 1 (mod 13)` then the
tight family contains a runner in class `u`. -/
theorem exists_mem_class_of_tight (S : Finset ℕ)
    (hnd : ∀ v ∈ S, ¬ (13 : ℕ) ∣ v)
    (htight : ∀ ε : ℝ, 0 < ε → ¬ ∃ t, IsLonelyAt S (1 / (13 : ℕ) + ε) t)
    {u a : ℕ} (hu13 : u < 13) (hua : u * a % 13 = 1) (ha : 0 < a)
    (hco : Nat.Coprime a 13) :
    ∃ v ∈ S, v % 13 = u := by
  obtain ⟨v, hv, hva⟩ := exists_residue_one_of_tight S (by norm_num) hnd htight ha hco
  refine ⟨v, hv, ?_⟩
  -- v ≡ v·(u·a) = u·(v·a) ≡ u (mod 13)
  have h1 : v * (u * a) ≡ v * 1 [MOD 13] := Nat.ModEq.mul_left v (by
    show u * a ≡ 1 [MOD 13]
    unfold Nat.ModEq
    simp [hua])
  have h2 : u * (v * a) ≡ u * 1 [MOD 13] := Nat.ModEq.mul_left u (by
    show v * a ≡ 1 [MOD 13]
    unfold Nat.ModEq
    simp [hva])
  have h3 : v * (u * a) = u * (v * a) := by ring
  have h4 : v ≡ u [MOD 13] := by
    calc v = v * 1 := by ring
      _ ≡ v * (u * a) [MOD 13] := h1.symm
      _ = u * (v * a) := h3
      _ ≡ u * 1 [MOD 13] := h2
      _ = u := by ring
  have := h4.symm
  unfold Nat.ModEq at this
  rw [← this, Nat.mod_eq_of_lt hu13]

/-- **Residue pinning (hdich Half 1).**  A 13-tight-from-above 12-family with no multiple of
13 has residue image exactly `{1, …, 12}`: at the prime modulus, tightness forces every
nonzero class. -/
theorem residue_pinning_13 (S : Finset ℕ) (hcard : S.card = 12)
    (hnd : ∀ v ∈ S, ¬ (13 : ℕ) ∣ v)
    (htight : ∀ ε : ℝ, 0 < ε → ¬ ∃ t, IsLonelyAt S (1 / (13 : ℕ) + ε) t) :
    S.image (· % 13) = Finset.Icc 1 12 := by
  have hsub : S.image (· % 13) ⊆ Finset.Icc 1 12 := by
    intro r hr
    obtain ⟨v, hv, rfl⟩ := Finset.mem_image.mp hr
    have hlt : v % 13 < 13 := Nat.mod_lt _ (by norm_num)
    have hne : v % 13 ≠ 0 := fun h => hnd v hv (Nat.dvd_of_mod_eq_zero h)
    exact Finset.mem_Icc.mpr ⟨by omega, by omega⟩
  have hsup : Finset.Icc 1 12 ⊆ S.image (· % 13) := by
    intro u hu
    obtain ⟨hu1, hu2⟩ := Finset.mem_Icc.mp hu
    -- the twelve inverse pairs (u, a) with u·a ≡ 1 (mod 13)
    have key : ∃ v ∈ S, v % 13 = u := by
      interval_cases u
      · exact exists_mem_class_of_tight S hnd htight (by norm_num) (by norm_num)
          (by norm_num) (by norm_num) (a := 1)
      · exact exists_mem_class_of_tight S hnd htight (by norm_num) (by norm_num)
          (by norm_num) (by norm_num) (a := 7)
      · exact exists_mem_class_of_tight S hnd htight (by norm_num) (by norm_num)
          (by norm_num) (by norm_num) (a := 9)
      · exact exists_mem_class_of_tight S hnd htight (by norm_num) (by norm_num)
          (by norm_num) (by norm_num) (a := 10)
      · exact exists_mem_class_of_tight S hnd htight (by norm_num) (by norm_num)
          (by norm_num) (by norm_num) (a := 8)
      · exact exists_mem_class_of_tight S hnd htight (by norm_num) (by norm_num)
          (by norm_num) (by norm_num) (a := 11)
      · exact exists_mem_class_of_tight S hnd htight (by norm_num) (by norm_num)
          (by norm_num) (by norm_num) (a := 2)
      · exact exists_mem_class_of_tight S hnd htight (by norm_num) (by norm_num)
          (by norm_num) (by norm_num) (a := 5)
      · exact exists_mem_class_of_tight S hnd htight (by norm_num) (by norm_num)
          (by norm_num) (by norm_num) (a := 3)
      · exact exists_mem_class_of_tight S hnd htight (by norm_num) (by norm_num)
          (by norm_num) (by norm_num) (a := 4)
      · exact exists_mem_class_of_tight S hnd htight (by norm_num) (by norm_num)
          (by norm_num) (by norm_num) (a := 6)
      · exact exists_mem_class_of_tight S hnd htight (by norm_num) (by norm_num)
          (by norm_num) (by norm_num) (a := 12)
    obtain ⟨v, hv, hvu⟩ := key
    exact Finset.mem_image.mpr ⟨v, hv, hvu⟩
  exact Finset.Subset.antisymm hsub hsup

/-- The pinning is exact: the residue map is injective on the family (12 runners, 12 classes). -/
theorem residue_injOn_13 (S : Finset ℕ) (hcard : S.card = 12)
    (hnd : ∀ v ∈ S, ¬ (13 : ℕ) ∣ v)
    (htight : ∀ ε : ℝ, 0 < ε → ¬ ∃ t, IsLonelyAt S (1 / (13 : ℕ) + ε) t) :
    Set.InjOn (· % 13) (S : Set ℕ) := by
  have himg := residue_pinning_13 S hcard hnd htight
  have hcards : (S.image (· % 13)).card = S.card := by
    rw [himg, hcard]
    simp
  exact Finset.injOn_of_card_image_eq hcards

end ResiduePinning
end LonelyRunner
