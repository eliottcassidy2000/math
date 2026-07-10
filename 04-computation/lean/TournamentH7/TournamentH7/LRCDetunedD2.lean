/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-10-S210)
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRC13Citation

/-!
# THM-678, the `d = 2` case (generic): the two-detuned dispatch (opus-S210)

A 13-family `v = g·H ∪ {δ₁, δ₂}` with `g ≥ 2` dividing all but two coordinates is lonely, in the generic
regime `(q₁, q₂) ≠ (2, 2)` (`qⱼ = g / gcd(δⱼ, g)`). This discharges the generic part of `MultiDetunedDispatch`
(opus-S209), which peels the near-dilate families the residual floor cannot bound (opus-S208).

The mechanism generalises the `d = 1` dispatch (`DetunedDispatch.lonely14_of_detuned`, THM-668):

* the 11 harmonic speeds are `g·w`; the LRC(11) citation gives `u` with clearance `≥ 1/12`, and at EVERY
  branch `(u + c)/g` the harmonic clearance is unchanged (`(g·w)(u+c)/g = w·u + w·c`, integer shift);
* the two detuned speeds must be cleared by a SINGLE branch `c ∈ [0, g)`. The counting dispatch: the set of
  `c` failing coordinate `j` has size `≤ (g/qⱼ)·Nⱼ` with `Nⱼ = ⌊qⱼ/7⌋ + 1` (orbit points of `δⱼ` in the
  `1/7` danger arc); `Σ Nⱼ/qⱼ < 1` for `(q₁,q₂) ≠ (2,2)`, so a good branch exists.

This file builds the **construction core** (`lonely14_of_two_detuned_good`): given the harmonic clearances
and a good branch, the family is lonely. The LRC(11) reduction and the counting existence are the remaining
pieces (see the docstring TODO). The residual `(2,2)` congruent-half-harmonic pair needs the mod-`2g` lift
(THM-678's residual), separate.

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace DetunedD2

open LonelyRunner

/-- **The `d = 2` construction core.** Given `g ≥ 1` dividing every coordinate except `i₁, i₂`, harmonic
clearances `‖(vⱼ/g)·u‖ ≥ 1/14` (from LRC(11), `1/12 ≥ 1/14`), and a branch `c` at which BOTH detuned
speeds clear `≥ 1/14`, the family is lonely at `(u + c)/g`. -/
theorem lonely14_of_two_detuned_good (v : Fin 13 → ℤ) (g : ℤ) (hg0 : 0 < g)
    (i₁ i₂ : Fin 13) (u : ℝ) (c : ℤ)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → g ∣ v j)
    (hharm : ∀ j, j ≠ i₁ → j ≠ i₂ → ∀ n : ℤ, (1 : ℝ) / 14 ≤ |((v j / g : ℤ) : ℝ) * u - n|)
    (hd1 : ∀ n : ℤ, (1 : ℝ) / 14 ≤ |(v i₁ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|)
    (hd2 : ∀ n : ℤ, (1 : ℝ) / 14 ≤ |(v i₂ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hgR : (0 : ℝ) < (g : ℝ) := by exact_mod_cast hg0
  have hgR' : (g : ℝ) ≠ 0 := ne_of_gt hgR
  refine ⟨(u + (c : ℝ)) / (g : ℝ), fun i n => ?_⟩
  rcases eq_or_ne i i₁ with h1 | h1
  · rw [h1]; exact hd1 n
  rcases eq_or_ne i i₂ with h2 | h2
  · rw [h2]; exact hd2 n
  · -- harmonic coordinate: the branch shift is an integer, so the clearance equals `‖(vᵢ/g)·u‖`
    have hd := hdvd i h1 h2
    have hvi : (v i : ℝ) = (g : ℝ) * ((v i / g : ℤ) : ℝ) := by
      have : v i = g * (v i / g) := (Int.mul_ediv_cancel' hd).symm
      exact_mod_cast this
    have hval : (v i : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n
        = ((v i / g : ℤ) : ℝ) * u - (((n - (v i / g) * c : ℤ)) : ℝ) := by
      rw [hvi]; push_cast; field_simp; ring
    rw [hval]
    exact hharm i h1 h2 (n - (v i / g) * c)

/-- **The two-detuned clearing obligation (the counting core of THM-678, `d = 2`).** For distinct nonzero
`δ₁, δ₂` neither divisible by `g ≥ 1`, in the generic regime `(q₁, q₂) ≠ (2, 2)` (`qⱼ = g/gcd(δⱼ, g)`), a
branch `c` clears BOTH detuned phases at `(u + c)/g`, for every real `u`. This is exactly the
`Σ Nⱼ/qⱼ < 1 ⟹ good branch` counting dispatch — the sole remaining analytic piece of the `d = 2` case. -/
def TwoDetunedClearing : Prop :=
  ∀ (δ₁ δ₂ g : ℤ), 0 < g → ¬ (g ∣ δ₁) → ¬ (g ∣ δ₂) →
    ¬ (g / (Int.gcd δ₁ g : ℤ) = 2 ∧ g / (Int.gcd δ₂ g : ℤ) = 2) →
    ∀ u : ℝ, ∃ c : ℤ,
      (∀ n : ℤ, (1 : ℝ) / 14 ≤ |(δ₁ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|) ∧
      (∀ n : ℤ, (1 : ℝ) / 14 ≤ |(δ₂ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|)

/-- **THM-678, `d = 2` generic — reduced to the counting.** A family `v = g·H ∪ {δ₁, δ₂}` (`g ≥ 2`
dividing all but `i₁, i₂`, both `q ≠ (2,2)`) is lonely, GIVEN the LRC(≤13) citation and the two-detuned
clearing obligation. The LRC(11) reduction supplies the harmonic clearances (`1/12 ≥ 1/14` at every
branch); `TwoDetunedClearing` supplies the good branch; `lonely14_of_two_detuned_good` assembles them. -/
theorem lonely14_of_two_detuned (cite : LRCUpTo13) (hclear : TwoDetunedClearing)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₁ i₂ : Fin 13) (hne : i₁ ≠ i₂)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → g ∣ v j) (hδ1 : ¬ g ∣ v i₁) (hδ2 : ¬ g ∣ v i₂)
    (hq : ¬ (g / (Int.gcd (v i₁) g : ℤ) = 2 ∧ g / (Int.gcd (v i₂) g : ℤ) = 2)) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hg0 : (0 : ℤ) < g := by omega
  -- the 11 harmonic coordinates, reindexed by an order embedding of the complement
  have hcard : (Finset.univ \ ({i₁, i₂} : Finset (Fin 13))).card = 11 := by
    rw [Finset.card_sdiff, Finset.inter_univ, Finset.card_pair hne, Finset.card_univ, Fintype.card_fin]
  set emb : Fin 11 ↪o Fin 13 :=
    (Finset.univ \ ({i₁, i₂} : Finset (Fin 13))).orderEmbOfFin hcard with hemb
  have hemb_mem : ∀ k, emb k ∈ (Finset.univ \ ({i₁, i₂} : Finset (Fin 13))) :=
    fun k => (Finset.univ \ ({i₁, i₂} : Finset (Fin 13))).orderEmbOfFin_mem hcard k
  have hemb_ne : ∀ k, emb k ≠ i₁ ∧ emb k ≠ i₂ := by
    intro k
    have := hemb_mem k
    rw [Finset.mem_sdiff, Finset.mem_insert, Finset.mem_singleton] at this
    push_neg at this
    exact this.2
  set w : Fin 11 → ℤ := fun k => v (emb k) / g with hw
  have hwnz : ∀ k, w k ≠ 0 := by
    intro k hk
    exact hv (emb k) (by
      have hd := hdvd (emb k) (hemb_ne k).1 (hemb_ne k).2
      have : v (emb k) = g * w k := (Int.mul_ediv_cancel' hd).symm
      rw [this, hk, mul_zero])
  obtain ⟨u, hu⟩ := cite 11 (by norm_num) w hwnz
  -- harmonic clearances: every `j ≠ i₁, i₂` is `emb k`, with clearance `1/12 ≥ 1/14`
  have hharm : ∀ j, j ≠ i₁ → j ≠ i₂ → ∀ n : ℤ, (1 : ℝ) / 14 ≤ |((v j / g : ℤ) : ℝ) * u - n| := by
    intro j hj1 hj2 n
    have hjmem : j ∈ (Finset.univ \ ({i₁, i₂} : Finset (Fin 13))) := by
      rw [Finset.mem_sdiff, Finset.mem_insert, Finset.mem_singleton]
      exact ⟨Finset.mem_univ j, by push_neg; exact ⟨hj1, hj2⟩⟩
    obtain ⟨k, hk⟩ : ∃ k, emb k = j := by
      have hr : j ∈ Set.range emb := by
        rw [hemb, Finset.range_orderEmbOfFin]; exact hjmem
      exact hr
    have hwk : (w k : ℤ) = v j / g := by show v (emb k) / g = v j / g; rw [hk]
    calc (1 : ℝ) / 14 ≤ 1 / 12 := by norm_num
      _ ≤ |(w k : ℝ) * u - n| := hu k n
      _ = |((v j / g : ℤ) : ℝ) * u - n| := by rw [hwk]
  -- the good branch from the clearing obligation, then the construction core
  obtain ⟨c, hc1, hc2⟩ := hclear (v i₁) (v i₂) g hg0 hδ1 hδ2 hq u
  exact lonely14_of_two_detuned_good v g hg0 i₁ i₂ u c hdvd hharm hc1 hc2

end DetunedD2
end LonelyRunner
