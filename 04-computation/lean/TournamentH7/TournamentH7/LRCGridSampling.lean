/-
  TournamentH7.LRCGridSampling — THE WEIGHT-1 SAMPLING LEMMA (boxeph-2026-07-17-S78;
  the general form of kind-pasteur THM-984(I)'s hand-off).

  A modulus-q grid meets an interval (a, b) ⊆ [0, 1] in more than q(b−a) − 1
  points, and a SEPARATED family of n such intervals in at least
  q·(total length) − n points.  This is the abstract endpoint-counting engine
  of the live-floor bridge.  Applying it to LRC still requires a finite
  separated-interval decomposition of the strict safe interior, its component
  and length bounds, exact ceiling arithmetic, and the identification of its
  sampled grid points with `liveCount`.  Once that last positivity is known,
  `LRCLiveCountLonely` is the direct consumer; this file is deliberately
  interface-free (imports Mathlib only).

  Kernel-pure: no `native_decide`, no `sorry`.
-/
import Mathlib

namespace LonelyRunner.LRC14.Sampling

open Finset

/-- **The per-interval grid count**: the modulus-`q` grid points strictly inside
`(a, b) ⊆ [0, 1]` number more than `q(b−a) − 1`. -/
theorem card_grid_Ioo (q : ℕ) (hq : 0 < q) (a b : ℝ) (ha : 0 ≤ a) (hb : b ≤ 1) :
    (q : ℝ) * (b - a) - 1
      ≤ (((Finset.Ioo 0 q)).filter fun p : ℕ => a < (p : ℝ) / q ∧ (p : ℝ) / q < b).card := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  rcases le_or_gt ((q : ℝ) * (b - a)) 1 with htriv | hmain
  · calc (q : ℝ) * (b - a) - 1 ≤ 0 := by linarith
      _ ≤ _ := Nat.cast_nonneg _
  · -- main case: at least ⌈qb⌉ − ⌊qa⌋ − 1 grid points
    set n₁ : ℕ := (⌊(q : ℝ) * a⌋).toNat + 1 with hn₁
    set n₂ : ℕ := (⌈(q : ℝ) * b⌉).toNat with hn₂
    have hfloor_nn : 0 ≤ ⌊(q : ℝ) * a⌋ := Int.floor_nonneg.mpr (by positivity)
    have hceil_nn : 0 ≤ ⌈(q : ℝ) * b⌉ := by
      apply Int.ceil_nonneg
      nlinarith
    have hfl : (((⌊(q : ℝ) * a⌋).toNat : ℕ) : ℝ) = ((⌊(q : ℝ) * a⌋ : ℤ) : ℝ) := by
      exact_mod_cast Int.toNat_of_nonneg hfloor_nn
    have hcl : (((⌈(q : ℝ) * b⌉).toNat : ℕ) : ℝ) = ((⌈(q : ℝ) * b⌉ : ℤ) : ℝ) := by
      exact_mod_cast Int.toNat_of_nonneg hceil_nn
    have hsub : Finset.Ico n₁ n₂
        ⊆ ((Finset.Ioo 0 q)).filter fun p : ℕ => a < (p : ℝ) / q ∧ (p : ℝ) / q < b := by
      intro p hp
      rw [Finset.mem_Ico] at hp
      have hp1 : ((⌊(q : ℝ) * a⌋).toNat + 1 : ℕ) ≤ p := hp.1
      have hp2 : p < (⌈(q : ℝ) * b⌉).toNat := hp.2
      have hpa : (q : ℝ) * a < (p : ℝ) := by
        have h1 : (q : ℝ) * a < (⌊(q : ℝ) * a⌋ : ℝ) + 1 := Int.lt_floor_add_one _
        have h2 : ((⌊(q : ℝ) * a⌋ : ℤ) : ℝ) + 1 ≤ (p : ℝ) := by
          have hc : (((⌊(q : ℝ) * a⌋).toNat + 1 : ℕ) : ℝ) ≤ (p : ℝ) := by
            exact_mod_cast hp1
          rw [Nat.cast_add, Nat.cast_one, hfl] at hc
          linarith
        linarith
      have hpb : (p : ℝ) < (q : ℝ) * b := by
        have h1 : (p : ℝ) + 1 ≤ ((⌈(q : ℝ) * b⌉ : ℤ) : ℝ) := by
          have h : p + 1 ≤ (⌈(q : ℝ) * b⌉).toNat := hp2
          have hc : ((p + 1 : ℕ) : ℝ) ≤ (((⌈(q : ℝ) * b⌉).toNat : ℕ) : ℝ) := by
            exact_mod_cast h
          rw [Nat.cast_add, Nat.cast_one, hcl] at hc
          linarith
        have h2 : ((⌈(q : ℝ) * b⌉ : ℤ) : ℝ) < (q : ℝ) * b + 1 := Int.ceil_lt_add_one _
        linarith
      rw [Finset.mem_filter, Finset.mem_Ioo]
      have hppos : 0 < p := by
        by_contra h
        push_neg at h
        interval_cases p
        simp only [Nat.cast_zero] at hpa
        nlinarith
      have hpq : p < q := by
        have : (p : ℝ) < (q : ℝ) := by nlinarith
        exact_mod_cast this
      exact ⟨⟨hppos, hpq⟩, by rw [lt_div_iff₀ hqR]; linarith,
        by rw [div_lt_iff₀ hqR]; linarith⟩
    have hcard := Finset.card_le_card hsub
    have hcardIco : (Finset.Ico n₁ n₂).card = n₂ - n₁ := Nat.card_Ico n₁ n₂
    have hn2R : (q : ℝ) * b ≤ (n₂ : ℝ) := by
      have h := Int.le_ceil ((q : ℝ) * b)
      rw [hn₂]
      rw [hcl]
      linarith
    have hn1R : (n₁ : ℝ) ≤ (q : ℝ) * a + 1 := by
      have h := Int.floor_le ((q : ℝ) * a)
      rw [hn₁]
      rw [Nat.cast_add, Nat.cast_one, hfl]
      linarith
    have hn12 : n₁ ≤ n₂ := by
      by_contra h
      push_neg at h
      have : (n₂ : ℝ) < (n₁ : ℝ) := by exact_mod_cast h
      nlinarith
    calc (q : ℝ) * (b - a) - 1
        ≤ (n₂ : ℝ) - (n₁ : ℝ) := by nlinarith
      _ = ((n₂ - n₁ : ℕ) : ℝ) := by
          rw [Nat.cast_sub hn12]
      _ = ((Finset.Ico n₁ n₂).card : ℝ) := by rw [hcardIco]
      _ ≤ _ := by exact_mod_cast hcard

/-- **The separated-family sampling lemma** (the weight-1 bridge's general form):
`n` separated intervals inside `[0, 1]` with total length `L` meet the modulus-`q`
grid in more than `qL − n` points. -/
theorem card_grid_family (q : ℕ) (hq : 0 < q) (n : ℕ) (a b : Fin n → ℝ)
    (ha : ∀ i, 0 ≤ a i) (hb : ∀ i, b i ≤ 1)
    (hsep : ∀ i j, i ≠ j → b i ≤ a j ∨ b j ≤ a i) :
    (q : ℝ) * (∑ i, (b i - a i)) - n
      ≤ (((Finset.Ioo 0 q)).filter fun p : ℕ =>
          ∃ i, a i < (p : ℝ) / q ∧ (p : ℝ) / q < b i).card := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  set T : Fin n → Finset ℕ := fun i =>
    ((Finset.Ioo 0 q)).filter fun p : ℕ => a i < (p : ℝ) / q ∧ (p : ℝ) / q < b i with hT
  have hdisj : ∀ i ∈ Finset.univ, ∀ j ∈ Finset.univ, i ≠ j →
      Disjoint (T i) (T j) := by
    intro i _ j _ hij
    rw [Finset.disjoint_left]
    intro p hpi hpj
    simp only [hT, Finset.mem_filter] at hpi hpj
    rcases hsep i j hij with h | h
    · linarith [hpi.2.2, hpj.2.1]
    · linarith [hpj.2.2, hpi.2.1]
  have hsubU : Finset.univ.biUnion T
      ⊆ ((Finset.Ioo 0 q)).filter fun p : ℕ =>
          ∃ i, a i < (p : ℝ) / q ∧ (p : ℝ) / q < b i := by
    intro p hp
    rw [Finset.mem_biUnion] at hp
    obtain ⟨i, _, hpi⟩ := hp
    simp only [hT, Finset.mem_filter] at hpi
    rw [Finset.mem_filter]
    exact ⟨hpi.1, ⟨i, hpi.2⟩⟩
  calc (q : ℝ) * (∑ i, (b i - a i)) - n
      = ∑ i : Fin n, ((q : ℝ) * (b i - a i) - 1) := by
        rw [Finset.mul_sum, Finset.sum_sub_distrib, Finset.sum_const,
          Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one]
    _ ≤ ∑ i : Fin n, ((T i).card : ℝ) := by
        apply Finset.sum_le_sum
        intro i _
        exact card_grid_Ioo q hq (a i) (b i) (ha i) (hb i)
    _ = ((∑ i : Fin n, (T i).card : ℕ) : ℝ) := by
        rw [Nat.cast_sum]
    _ = (((Finset.univ.biUnion T).card : ℕ) : ℝ) := by
        rw [Finset.card_biUnion hdisj]
    _ ≤ _ := by
        exact_mod_cast Finset.card_le_card hsubU

#print axioms card_grid_Ioo
#print axioms card_grid_family

end LonelyRunner.LRC14.Sampling
