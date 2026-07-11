/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: kind-pasteur (LRC multi-agent project, 2026-07-10-S127)
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRC13Citation
import TournamentH7.LRCIntervalCount

/-!
# THM-678, the `d = 3` case (generic): the three-detuned dispatch (kps-S127)

The `d = 3` generalization of opus's two-detuned dispatch (`LRCDetunedD2` S210 + `LRCTwoDetunedClearing`
S211).  A 13-family `v = g·H ∪ {δ₁, δ₂, δ₃}` with `g ≥ 2` dividing all but THREE coordinates is lonely,
in the generic regime where the three bad-branch counts sum to `< g` — i.e. `Σⱼ Nⱼ/qⱼ < 1` with
`Nⱼ = ⌊qⱼ/7⌋ + 1`, `qⱼ = g/gcd(δⱼ,g)`.

The mechanism is identical to `d = 1, 2`, one coordinate wider:

* the 10 harmonic speeds are `g·w`; the LRC(10) citation gives `u` with clearance `≥ 1/11 ≥ 1/14`, and at
  EVERY branch `(u + c)/g` the harmonic clearance is unchanged (integer shift);
* the three detuned speeds must be cleared by a SINGLE branch `c ∈ [0, g)`. **opus's per-coordinate count
  `LRCIntervalCount.bad_count_le` is reused verbatim for the third coordinate**: the set of `c` failing
  coordinate `j` has size `≤ dⱼ·(⌊qⱼ/7⌋+1)`, and if the three sizes sum to `< g` the three bad sets cannot
  cover `[0, g)`, so a good branch exists (a three-set union bound).

Contents:
* `lonely14_of_three_detuned_good` — the construction core (three cleared detuned + harmonic ⟹ lonely);
* `ThreeDetunedClearing` — the counting obligation, stated with the explicit `Σ dⱼNⱼ < g` hypothesis
  (the exact shape `bad_count_le` produces);
* `threeDetunedClearing` — **PROVED** by the three-set union bound (kernel-pure);
* `lonely14_of_three_detuned` / `lonely14_of_three_detuned'` — THM-678 `d = 3` generic, reduced to /
  unconditional from the LRC(≤13) citation.

The EXCEPTIONAL triples (where `Σ dⱼNⱼ ≥ g`, so this generic bound fails) are enumerated in
`lrc14_three_detuned_exceptional_kps_S127` — the `(q₁,q₂,q₃)` needing the mod-lift residual, the `d = 3`
analogue of opus's `(2,2)` pair. Those are NOT closed here (separate, like THM-678's `d = 2` residual).

Kernel-pure: no `sorry`, no `native_decide`.  Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace DetunedD3

open LonelyRunner
open scoped Classical

/-- **The `d = 3` construction core.** Given `g ≥ 1` dividing every coordinate except `i₁, i₂, i₃`, harmonic
clearances `‖(vⱼ/g)·u‖ ≥ 1/14` (from LRC(10), `1/11 ≥ 1/14`), and a branch `c` at which ALL THREE detuned
speeds clear `≥ 1/14`, the family is lonely at `(u + c)/g`. -/
theorem lonely14_of_three_detuned_good (v : Fin 13 → ℤ) (g : ℤ) (hg0 : 0 < g)
    (i₁ i₂ i₃ : Fin 13) (u : ℝ) (c : ℤ)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (hharm : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → ∀ n : ℤ, (1 : ℝ) / 14 ≤ |((v j / g : ℤ) : ℝ) * u - n|)
    (hd1 : ∀ n : ℤ, (1 : ℝ) / 14 ≤ |(v i₁ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|)
    (hd2 : ∀ n : ℤ, (1 : ℝ) / 14 ≤ |(v i₂ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|)
    (hd3 : ∀ n : ℤ, (1 : ℝ) / 14 ≤ |(v i₃ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hgR : (0 : ℝ) < (g : ℝ) := by exact_mod_cast hg0
  refine ⟨(u + (c : ℝ)) / (g : ℝ), fun i n => ?_⟩
  rcases eq_or_ne i i₁ with h1 | h1
  · rw [h1]; exact hd1 n
  rcases eq_or_ne i i₂ with h2 | h2
  · rw [h2]; exact hd2 n
  rcases eq_or_ne i i₃ with h3 | h3
  · rw [h3]; exact hd3 n
  · -- harmonic coordinate: the branch shift is an integer, so the clearance equals `‖(vᵢ/g)·u‖`
    have hd := hdvd i h1 h2 h3
    have hvi : (v i : ℝ) = (g : ℝ) * ((v i / g : ℤ) : ℝ) := by
      have : v i = g * (v i / g) := (Int.mul_ediv_cancel' hd).symm
      exact_mod_cast this
    have hval : (v i : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n
        = ((v i / g : ℤ) : ℝ) * u - (((n - (v i / g) * c : ℤ)) : ℝ) := by
      rw [hvi]; push_cast; field_simp; ring
    rw [hval]
    exact hharm i h1 h2 h3 (n - (v i / g) * c)

/-- The THM-678 bad-branch count for one detuned speed: `d·(⌊q/7⌋+1)` with `d = gcd(δ,g)`, `q = g/d`
(exactly the bound `LRCIntervalCount.bad_count_le` proves). -/
def badCount (δ g : ℤ) : ℕ := (Int.gcd δ g) * ((g / (Int.gcd δ g : ℤ)).toNat / 7 + 1)

/-- **The three-detuned clearing obligation (the counting core of THM-678, `d = 3`).** For distinct nonzero
`δ₁, δ₂, δ₃` none divisible by `g ≥ 1`, in the generic regime where the three bad-branch counts sum to
`< g` (`Σⱼ badCount δⱼ g < g`, i.e. `Σⱼ Nⱼ/qⱼ < 1`), a single branch `c` clears ALL THREE detuned phases at
`(u + c)/g`, for every real `u`. -/
def ThreeDetunedClearing : Prop :=
  ∀ (δ₁ δ₂ δ₃ g : ℤ), 0 < g → ¬ (g ∣ δ₁) → ¬ (g ∣ δ₂) → ¬ (g ∣ δ₃) →
    badCount δ₁ g + badCount δ₂ g + badCount δ₃ g < g.toNat →
    ∀ u : ℝ, ∃ c : ℤ,
      (∀ n : ℤ, (1 : ℝ) / 14 ≤ |(δ₁ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|) ∧
      (∀ n : ℤ, (1 : ℝ) / 14 ≤ |(δ₂ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|) ∧
      (∀ n : ℤ, (1 : ℝ) / 14 ≤ |(δ₃ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n|)

/-- **THM-678 `d = 3`, the counting — PROVED.** The three-detuned clearing holds whenever the three
per-coordinate counts sum to `< g`.  Three-set union bound over `LRCIntervalCount.bad_count_le`. -/
theorem threeDetunedClearing : ThreeDetunedClearing := by
  intro δ₁ δ₂ δ₃ g hg _ _ _ hsum u
  have hg1 : (1 : ℤ) ≤ g := hg
  set bad₁ := (Finset.Ico (0 : ℤ) g).filter
      (fun (c : ℤ) => ∃ n : ℤ, |(δ₁ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n| < 1 / 14) with hbad₁
  set bad₂ := (Finset.Ico (0 : ℤ) g).filter
      (fun (c : ℤ) => ∃ n : ℤ, |(δ₂ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n| < 1 / 14) with hbad₂
  set bad₃ := (Finset.Ico (0 : ℤ) g).filter
      (fun (c : ℤ) => ∃ n : ℤ, |(δ₃ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n| < 1 / 14) with hbad₃
  have hb₁ : bad₁.card ≤ badCount δ₁ g := LRCIntervalCount.bad_count_le δ₁ g hg1 u
  have hb₂ : bad₂.card ≤ badCount δ₂ g := LRCIntervalCount.bad_count_le δ₂ g hg1 u
  have hb₃ : bad₃.card ≤ badCount δ₃ g := LRCIntervalCount.bad_count_le δ₃ g hg1 u
  -- the three bad sets cannot cover `[0, g)`
  have hIco : (Finset.Ico (0 : ℤ) g).card = g.toNat := by
    rw [Int.card_Ico]; congr 1; omega
  have hcard_lt : (bad₁ ∪ bad₂ ∪ bad₃).card < (Finset.Ico (0 : ℤ) g).card := by
    calc (bad₁ ∪ bad₂ ∪ bad₃).card ≤ (bad₁ ∪ bad₂).card + bad₃.card := Finset.card_union_le _ _
      _ ≤ (bad₁.card + bad₂.card) + bad₃.card := Nat.add_le_add_right (Finset.card_union_le _ _) _
      _ ≤ (badCount δ₁ g + badCount δ₂ g) + badCount δ₃ g :=
          Nat.add_le_add (Nat.add_le_add hb₁ hb₂) hb₃
      _ = badCount δ₁ g + badCount δ₂ g + badCount δ₃ g := by ring
      _ < g.toNat := hsum
      _ = (Finset.Ico (0 : ℤ) g).card := hIco.symm
  have hnotsub : ¬ (Finset.Ico (0 : ℤ) g ⊆ bad₁ ∪ bad₂ ∪ bad₃) := by
    intro h; exact absurd (Finset.card_le_card h) (by omega)
  rw [Finset.not_subset] at hnotsub
  obtain ⟨c, hcIco, hcnot⟩ := hnotsub
  rw [Finset.mem_union, not_or, Finset.mem_union, not_or] at hcnot
  obtain ⟨⟨hc1, hc2⟩, hc3⟩ := hcnot
  -- at `c` all three detuned phases clear `1/14`
  refine ⟨c, ?_, ?_, ?_⟩
  · intro n; refine not_lt.mp (fun h => hc1 ?_)
    rw [hbad₁, Finset.mem_filter]; exact ⟨hcIco, n, h⟩
  · intro n; refine not_lt.mp (fun h => hc2 ?_)
    rw [hbad₂, Finset.mem_filter]; exact ⟨hcIco, n, h⟩
  · intro n; refine not_lt.mp (fun h => hc3 ?_)
    rw [hbad₃, Finset.mem_filter]; exact ⟨hcIco, n, h⟩

/-- **THM-678, `d = 3` generic — reduced to the counting.** A family `v = g·H ∪ {δ₁, δ₂, δ₃}` (`g ≥ 2`
dividing all but `i₁, i₂, i₃`, generic counting) is lonely, GIVEN the LRC(≤13) citation and the
three-detuned clearing obligation.  The LRC(10) reduction supplies the harmonic clearances
(`1/11 ≥ 1/14` at every branch); `ThreeDetunedClearing` supplies the good branch;
`lonely14_of_three_detuned_good` assembles them. -/
theorem lonely14_of_three_detuned (cite : LRCUpTo13) (hclear : ThreeDetunedClearing)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₁ i₂ i₃ : Fin 13) (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (hδ1 : ¬ g ∣ v i₁) (hδ2 : ¬ g ∣ v i₂) (hδ3 : ¬ g ∣ v i₃)
    (hcount : badCount (v i₁) g + badCount (v i₂) g + badCount (v i₃) g < g.toNat) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hg0 : (0 : ℤ) < g := by omega
  -- the 10 harmonic coordinates, reindexed by an order embedding of the complement
  have hc3 : ({i₁, i₂, i₃} : Finset (Fin 13)).card = 3 :=
    Finset.card_eq_three.mpr ⟨i₁, i₂, i₃, h12, h13, h23, rfl⟩
  have hcard : (Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))).card = 10 := by
    rw [Finset.card_sdiff, Finset.inter_univ, hc3, Finset.card_univ, Fintype.card_fin]
  set emb : Fin 10 ↪o Fin 13 :=
    (Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))).orderEmbOfFin hcard with hemb
  have hemb_mem : ∀ k, emb k ∈ (Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))) :=
    fun k => (Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))).orderEmbOfFin_mem hcard k
  have hemb_ne : ∀ k, emb k ≠ i₁ ∧ emb k ≠ i₂ ∧ emb k ≠ i₃ := by
    intro k
    have := hemb_mem k
    rw [Finset.mem_sdiff, Finset.mem_insert, Finset.mem_insert, Finset.mem_singleton] at this
    push_neg at this
    exact this.2
  set w : Fin 10 → ℤ := fun k => v (emb k) / g with hw
  have hwnz : ∀ k, w k ≠ 0 := by
    intro k hk
    exact hv (emb k) (by
      have hd := hdvd (emb k) (hemb_ne k).1 (hemb_ne k).2.1 (hemb_ne k).2.2
      have : v (emb k) = g * w k := (Int.mul_ediv_cancel' hd).symm
      rw [this, hk, mul_zero])
  obtain ⟨u, hu⟩ := cite 10 (by norm_num) w hwnz
  -- harmonic clearances: every `j ∉ {i₁,i₂,i₃}` is `emb k`, with clearance `1/11 ≥ 1/14`
  have hharm : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ →
      ∀ n : ℤ, (1 : ℝ) / 14 ≤ |((v j / g : ℤ) : ℝ) * u - n| := by
    intro j hj1 hj2 hj3 n
    have hjmem : j ∈ (Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))) := by
      rw [Finset.mem_sdiff, Finset.mem_insert, Finset.mem_insert, Finset.mem_singleton]
      exact ⟨Finset.mem_univ j, by push_neg; exact ⟨hj1, hj2, hj3⟩⟩
    obtain ⟨k, hk⟩ : ∃ k, emb k = j := by
      have hr : j ∈ Set.range emb := by
        rw [hemb, Finset.range_orderEmbOfFin]; exact hjmem
      exact hr
    have hwk : (w k : ℤ) = v j / g := by show v (emb k) / g = v j / g; rw [hk]
    calc (1 : ℝ) / 14 ≤ 1 / 11 := by norm_num
      _ ≤ |(w k : ℝ) * u - n| := hu k n
      _ = |((v j / g : ℤ) : ℝ) * u - n| := by rw [hwk]
  -- the good branch from the clearing obligation, then the construction core
  obtain ⟨c, hc1, hc2, hc3'⟩ := hclear (v i₁) (v i₂) (v i₃) g hg0 hδ1 hδ2 hδ3 hcount u
  exact lonely14_of_three_detuned_good v g hg0 i₁ i₂ i₃ u c hdvd hharm hc1 hc2 hc3'

/-- **THM-678 `d = 3` generic, unconditional.** A family `v = g·H ∪ {δ₁, δ₂, δ₃}` (`g ≥ 2` dividing all but
`i₁, i₂, i₃`, generic counting `Σ badCount < g`) is lonely, from the LRC(≤13) citation alone — the counting
is discharged by `threeDetunedClearing`. -/
theorem lonely14_of_three_detuned' (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₁ i₂ i₃ : Fin 13) (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (hδ1 : ¬ g ∣ v i₁) (hδ2 : ¬ g ∣ v i₂) (hδ3 : ¬ g ∣ v i₃)
    (hcount : badCount (v i₁) g + badCount (v i₂) g + badCount (v i₃) g < g.toNat) :
    ∃ t : ℝ, Lonely 14 v t :=
  lonely14_of_three_detuned cite threeDetunedClearing v hv g hg i₁ i₂ i₃ h12 h13 h23
    hdvd hδ1 hδ2 hδ3 hcount

end DetunedD3
end LonelyRunner
