/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-09-S195)
-/
import Mathlib

/-!
# The Freiman AP step: minimal restricted sumset ⟹ AP, for `n ≥ 5` (opus-S195)

At the MINIMAL descent burden `|A +̂ A| = 2n − 3`, a set of `n ≥ 5` reals is an arithmetic progression.
This is the near-AP characterization THM-676 needs for majority-parity 7-classes (`n = 7 ≥ 5`).

**The `n ≥ 5` hypothesis is ESSENTIAL** (MISTAKE-133): for `n ≤ 4` the equality is achieved by NON-AP
"bi-arithmetic" sets (e.g. `{0,1,3,4}`: `A +̂ A = {1,3,4,5,7}`, card `5 = 2·4−3`, differences `1,2,1`).
The census `lrc14_ap_step_verify_opus_S195.out` confirms: non-AP minimal sets number 44, 38, 0, 0, 0, 0
for `n = 3..8`. The step holds **iff `n ≥ 5`**.

**Proof (the interleaved-chain route, cleaner than the S189 row-bijection+reflection blueprint).** The
`2n−3` sums `sᵢ = aᵢ + aᵢ₊₁` (`0 ≤ i ≤ n−2`) and `tᵢ = aᵢ + aᵢ₊₂` (`0 ≤ i ≤ n−3`) are distinct
(they strictly interleave `s₀ < t₀ < s₁ < t₁ < ⋯`), so `Rset = IC := {sᵢ} ∪ {tᵢ}` at minimal burden.
Then:

* **Step 1** (`d_{i+2} = d_i`): `aᵢ + aᵢ₊₃ ∈ Rset = IC` lies strictly in `(tᵢ, tᵢ₊₁)`, whose only chain
  element is `sᵢ₊₁`, forcing `aᵢ + aᵢ₊₃ = aᵢ₊₁ + aᵢ₊₂`, i.e. `d_{i+2} = d_i`.
* **Step 2** (`d₀ = d₁`, needs `n ≥ 5`): `a₀ + a₄ ∈ Rset = IC`; the strict order rules out every `sₖ` and
  every `tₖ` except `t₁ = a₁ + a₃`, forcing `a₀ + a₄ = a₁ + a₃`, i.e. `d₃ = d₀`; with `d₃ = d₁` (Step 1)
  this gives `d₀ = d₁`, breaking the even/odd ambiguity.

Together with `d_{i+2} = d_i` this makes all consecutive differences equal — an AP.

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LRCFreimanAP

open Finset

/-- The **restricted sumset** over the first `n` indices of a sequence `a`: `{aᵢ + aⱼ : i < j < n}`. -/
def Rset (a : ℕ → ℤ) (n : ℕ) : Finset ℤ :=
  ((range n ×ˢ range n).filter (fun p => p.1 < p.2)).image (fun p => a p.1 + a p.2)

theorem mem_Rset {a : ℕ → ℤ} {n : ℕ} {z : ℤ} :
    z ∈ Rset a n ↔ ∃ i j, i < n ∧ j < n ∧ i < j ∧ a i + a j = z := by
  simp only [Rset, mem_image, mem_filter, mem_product, mem_range, Prod.exists]
  constructor
  · rintro ⟨i, j, ⟨⟨hi, hj⟩, hlt⟩, heq⟩; exact ⟨i, j, hi, hj, hlt, heq⟩
  · rintro ⟨i, j, hi, hj, hlt, heq⟩; exact ⟨i, j, ⟨⟨hi, hj⟩, hlt⟩, heq⟩

/-- The consecutive-pair sums `sᵢ = aᵢ + aᵢ₊₁`, `0 ≤ i < n−1`. -/
def sImg (a : ℕ → ℤ) (n : ℕ) : Finset ℤ := (range (n - 1)).image (fun i => a i + a (i + 1))
/-- The skip-pair sums `tᵢ = aᵢ + aᵢ₊₂`, `0 ≤ i < n−2`. -/
def tImg (a : ℕ → ℤ) (n : ℕ) : Finset ℤ := (range (n - 2)).image (fun i => a i + a (i + 2))
/-- The **interleaved chain** = consecutive ∪ skip sums. At minimal burden this IS the restricted sumset. -/
def IC (a : ℕ → ℤ) (n : ℕ) : Finset ℤ := sImg a n ∪ tImg a n

theorem sImg_strictMono {a : ℕ → ℤ} (ha : StrictMono a) :
    StrictMono (fun i => a i + a (i + 1)) := fun i j h => add_lt_add (ha h) (ha (by omega))

theorem tImg_strictMono {a : ℕ → ℤ} (ha : StrictMono a) :
    StrictMono (fun i => a i + a (i + 2)) := fun i j h => add_lt_add (ha h) (ha (by omega))

theorem sImg_card {a : ℕ → ℤ} (ha : StrictMono a) (n : ℕ) : (sImg a n).card = n - 1 := by
  rw [sImg, card_image_of_injective _ (sImg_strictMono ha).injective, card_range]

theorem tImg_card {a : ℕ → ℤ} (ha : StrictMono a) (n : ℕ) : (tImg a n).card = n - 2 := by
  rw [tImg, card_image_of_injective _ (tImg_strictMono ha).injective, card_range]

/-- Consecutive and skip sums are disjoint (the interleaving `sᵢ ≠ tₖ`). -/
theorem sImg_tImg_disjoint {a : ℕ → ℤ} (ha : StrictMono a) (n : ℕ) :
    Disjoint (sImg a n) (tImg a n) := by
  rw [disjoint_left]
  intro z hz1 hz2
  rw [sImg, mem_image] at hz1; obtain ⟨i, _, rfl⟩ := hz1
  rw [tImg, mem_image] at hz2; obtain ⟨k, _, hk2⟩ := hz2
  rcases Nat.lt_or_ge i (k + 1) with h | h
  · have h1 : a i ≤ a k := ha.monotone (by omega)
    have h2 : a (i + 1) < a (k + 2) := ha (by omega)
    omega
  · have h1 : a (k + 1) ≤ a i := ha.monotone (by omega)
    have h2 : a k < a (k + 1) := ha (by omega)
    have h3 : a (k + 2) ≤ a (i + 1) := ha.monotone (by omega)
    omega

theorem IC_card {a : ℕ → ℤ} (ha : StrictMono a) {n : ℕ} (hn : 2 ≤ n) :
    (IC a n).card = 2 * n - 3 := by
  rw [IC, card_union_of_disjoint (sImg_tImg_disjoint ha n), sImg_card ha, tImg_card ha]
  omega

theorem IC_subset_Rset {a : ℕ → ℤ} {n : ℕ} : IC a n ⊆ Rset a n := by
  rw [IC]
  apply union_subset
  · intro z hz; rw [sImg, mem_image] at hz; obtain ⟨i, hi, rfl⟩ := hz; rw [mem_range] at hi
    exact mem_Rset.mpr ⟨i, i + 1, by omega, by omega, by omega, rfl⟩
  · intro z hz; rw [tImg, mem_image] at hz; obtain ⟨i, hi, rfl⟩ := hz; rw [mem_range] at hi
    exact mem_Rset.mpr ⟨i, i + 2, by omega, by omega, by omega, rfl⟩

/-- **At minimal burden the restricted sumset IS the interleaved chain.** -/
theorem Rset_eq_IC {a : ℕ → ℤ} (ha : StrictMono a) {n : ℕ} (hn : 2 ≤ n)
    (hcard : (Rset a n).card = 2 * n - 3) : Rset a n = IC a n := by
  refine (eq_of_subset_of_card_le IC_subset_Rset ?_).symm
  rw [IC_card ha hn]; exact le_of_eq hcard

/-- **Step 1 (`d_{i+2} = d_i`).** At minimal burden, `aᵢ + aᵢ₊₃` (a valid restricted sum) lies strictly
between the skip sums `tᵢ` and `tᵢ₊₁`, whose only interleaved-chain element is `sᵢ₊₁ = aᵢ₊₁ + aᵢ₊₂` —
forcing `aᵢ + aᵢ₊₃ = aᵢ₊₁ + aᵢ₊₂`. Equivalently the consecutive differences two apart are equal. -/
theorem diff_two {a : ℕ → ℤ} (ha : StrictMono a) {n : ℕ} (hn : 2 ≤ n)
    (hcard : (Rset a n).card = 2 * n - 3) (i : ℕ) (hi : i + 3 < n) :
    a i + a (i + 3) = a (i + 1) + a (i + 2) := by
  have hmem : a i + a (i + 3) ∈ Rset a n :=
    mem_Rset.mpr ⟨i, i + 3, by omega, by omega, by omega, rfl⟩
  rw [Rset_eq_IC ha hn hcard, IC, mem_union] at hmem
  rcases hmem with hs | ht
  · rw [sImg, mem_image] at hs
    obtain ⟨k, hk, hks⟩ := hs
    rw [mem_range] at hk
    change a k + a (k + 1) = a i + a (i + 3) at hks
    rcases lt_trichotomy k (i + 1) with h | h | h
    · exfalso
      have h1 : a k ≤ a i := ha.monotone (by omega)
      have h2 : a (k + 1) < a (i + 3) := ha (by omega)
      omega
    · rw [h] at hks
      have e : i + 1 + 1 = i + 2 := by omega
      rw [e] at hks
      omega
    · exfalso
      have h1 : a (i + 2) ≤ a k := ha.monotone (by omega)
      have h2 : a (i + 3) ≤ a (k + 1) := ha.monotone (by omega)
      have h3 : a i < a (i + 2) := ha (by omega)
      omega
  · exfalso
    rw [tImg, mem_image] at ht
    obtain ⟨k, hk, hkt⟩ := ht
    rw [mem_range] at hk
    change a k + a (k + 2) = a i + a (i + 3) at hkt
    rcases lt_trichotomy k i with h | h | h
    · have h1 : a k < a i := ha (by omega)
      have h2 : a (k + 2) < a (i + 3) := ha (by omega)
      omega
    · rw [h] at hkt
      have h3 : a (i + 2) < a (i + 3) := ha (by omega)
      omega
    · have h1 : a (i + 1) ≤ a k := ha.monotone (by omega)
      have h2 : a (i + 3) ≤ a (k + 2) := ha.monotone (by omega)
      have h3 : a i < a (i + 1) := ha (by omega)
      omega

/-- **Step 2 (`d₀ = d₁`, needs `n ≥ 5`).** The sum `a₀ + a₄` is a restricted sum, hence in the
interleaved chain; the strict order eliminates every consecutive sum `sₖ` (via Step 1 at `i = 0, 1`)
and every skip sum `tₖ` except `t₁ = a₁ + a₃` — so `a₀ + a₄ = a₁ + a₃`. -/
theorem sum04 {a : ℕ → ℤ} (ha : StrictMono a) {n : ℕ} (hn : 5 ≤ n)
    (hcard : (Rset a n).card = 2 * n - 3) : a 0 + a 4 = a 1 + a 3 := by
  have hmem : a 0 + a 4 ∈ Rset a n :=
    mem_Rset.mpr ⟨0, 4, by omega, by omega, by omega, rfl⟩
  rw [Rset_eq_IC ha (by omega) hcard, IC, mem_union] at hmem
  rcases hmem with hs | ht
  · exfalso
    rw [sImg, mem_image] at hs; obtain ⟨k, hk, hks⟩ := hs; rw [mem_range] at hk
    change a k + a (k + 1) = a 0 + a 4 at hks
    rcases Nat.lt_or_ge k 3 with hlt | hge
    · interval_cases k
      · have h := ha (show (1 : ℕ) < 4 by omega); simp only [Nat.reduceAdd] at hks; omega
      · have hd0 := diff_two ha (by omega) hcard 0 (by omega)
        have h := ha (show (3 : ℕ) < 4 by omega)
        simp only [Nat.reduceAdd] at hks hd0; omega
      · have hd1 := diff_two ha (by omega) hcard 1 (by omega)
        have h := ha (show (0 : ℕ) < 1 by omega)
        simp only [Nat.reduceAdd] at hks hd1; omega
    · have h1 : a 3 ≤ a k := ha.monotone (by omega)
      have h2 : a 4 ≤ a (k + 1) := ha.monotone (by omega)
      have h3 : a 0 < a 3 := ha (by omega)
      omega
  · rw [tImg, mem_image] at ht; obtain ⟨k, hk, hkt⟩ := ht; rw [mem_range] at hk
    change a k + a (k + 2) = a 0 + a 4 at hkt
    rcases Nat.lt_or_ge k 2 with hlt | hge
    · interval_cases k
      · exfalso; have h := ha (show (2 : ℕ) < 4 by omega)
        simp only [Nat.reduceAdd] at hkt; omega
      · simp only [Nat.reduceAdd] at hkt; omega
    · exfalso
      have h1 : a 2 ≤ a k := ha.monotone (by omega)
      have h2 : a 4 ≤ a (k + 2) := ha.monotone (by omega)
      have h3 : a 0 < a 2 := ha (by omega)
      omega

/-- **The Freiman AP step (`n ≥ 5`).** A strictly increasing sequence whose first `n` terms have the
MINIMAL restricted-sumset size `2n − 3` is an arithmetic progression: every consecutive difference
equals `a₁ − a₀`. (False for `n ≤ 4`; MISTAKE-133.) This is the near-AP characterization THM-676 uses
for majority-parity 7-classes. -/
theorem ap_of_min_burden {a : ℕ → ℤ} (ha : StrictMono a) {n : ℕ} (hn : 5 ≤ n)
    (hcard : (Rset a n).card = 2 * n - 3) :
    ∀ i, i + 1 < n → a (i + 1) - a i = a 1 - a 0 := by
  have hd01 : a 2 - a 1 = a 1 - a 0 := by
    have h04 := sum04 ha hn hcard
    have h14 := diff_two ha (by omega) hcard 1 (by omega)
    simp only [Nat.reduceAdd] at h14
    omega
  intro i
  induction i using Nat.strong_induction_on with
  | _ i IH =>
    intro hi
    rcases Nat.lt_or_ge i 2 with h2 | h2
    · interval_cases i
      · rfl
      · exact hd01
    · obtain ⟨j, rfl⟩ : ∃ j, i = j + 2 := ⟨i - 2, by omega⟩
      have hdj := diff_two ha (by omega) hcard j (by omega)
      have hIH := IH j (by omega) (by omega)
      show a (j + 3) - a (j + 2) = a 1 - a 0
      omega

end LRCFreimanAP
