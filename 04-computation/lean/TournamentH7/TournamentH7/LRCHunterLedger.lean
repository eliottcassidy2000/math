/-
  TournamentH7.LRCHunterLedger — THE PATH-HUNTER (BONFERRONI) MEASURE INEQUALITY
  (klein-2026-07-02-S116, HYP-4021) — the combinatorial heart of the 7-wall crossing
  designed in kps-S22 (LRCFatBlockChain, HYP-3979).

  The union bound `μ(⋃ Dᵢ) ≤ Σ μ(Dᵢ)` dies at j = 7 (each danger arc has measure 1/7, so
  seven tile).  Hunter's inequality (the second-order Bonferroni along a spanning PATH)
  recovers a pair credit:

      μ(⋃_{i<n} Aᵢ) + Σ_{i=1}^{n-1} μ(Aᵢ ∩ Aᵢ₋₁)  ≤  Σ_{i<n} μ(Aᵢ).

  Proof: pure disjointification — `μ(S ∪ T) + μ(S ∩ T) = μ S + μ T`, and `Aᵢ₋₁ ⊆ ⋃_{j<i}`,
  so the running-intersection term dominates the path term.  NO ENNReal subtraction, NO
  measurability beyond each `Aᵢ`.  This is exactly the input kps's path-Bonferroni ledger
  needs; combined with the pair-floor `μ(I ∩ Dᵢ ∩ Dᵢ₊₁) ≥ |I|/49 − err` (mac-mini's
  JointRateCore per-cell obligation) it gives the ledger

      good ≥ |I|·(1 − c/7 + (c−1)/49) − fees = |I|·(48 − 6c)/49 − fees,

  POSITIVE for c ≤ 7 (the `ledger_coeff` identity below), crossing the 7-wall.
-/
import Mathlib

open MeasureTheory Finset

namespace LonelyRunner.LRC14.Hunter

variable {α : Type*} [MeasurableSpace α]

/-- **The path-Hunter / second-order path-Bonferroni inequality** (additive form, no
subtraction). -/
theorem path_hunter_add_le (μ : Measure α) (A : ℕ → Set α)
    (hA : ∀ i, MeasurableSet (A i)) (n : ℕ) :
    μ (⋃ i ∈ Finset.range n, A i) + ∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (i - 1))
      ≤ ∑ i ∈ Finset.range n, μ (A i) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rcases Nat.eq_zero_or_pos n with hn | hn
    · subst hn; simp
    · -- split the top terms off range (n+1) and Ico 1 (n+1)
      have hUunion : (⋃ i ∈ Finset.range (n + 1), A i)
          = (⋃ i ∈ Finset.range n, A i) ∪ A n := by
        ext x
        simp only [Set.mem_iUnion, Finset.mem_range, Set.mem_union]
        constructor
        · rintro ⟨i, hi, hx⟩
          rcases Nat.lt_succ_iff_lt_or_eq.mp hi with h | h
          · exact Or.inl ⟨i, h, hx⟩
          · exact Or.inr (h ▸ hx)
        · rintro (⟨i, hi, hx⟩ | hx)
          · exact ⟨i, Nat.lt_succ_of_lt hi, hx⟩
          · exact ⟨n, Nat.lt_succ_self n, hx⟩
      set S : Set α := ⋃ i ∈ Finset.range n, A i with hS
      have hSmeas : MeasurableSet S := by
        rw [hS]; exact Finset.measurableSet_biUnion _ (fun i _ => hA i)
      -- key exact identity: μ(S ∪ A n) + μ(S ∩ A n) = μ S + μ (A n)
      have hkey : μ (S ∪ A n) + μ (S ∩ A n) = μ S + μ (A n) :=
        measure_union_add_inter S (hA n)
      -- the path term at the top is dominated by the running-intersection term
      have hsub : A (n - 1) ∩ A n ⊆ S ∩ A n := by
        apply Set.inter_subset_inter_left
        rw [hS]
        exact Set.subset_biUnion_of_mem (Finset.mem_range.mpr (Nat.sub_lt hn Nat.one_pos))
      have hmono : μ (A (n - 1) ∩ A n) ≤ μ (S ∩ A n) := measure_mono hsub
      have htop : μ (A n ∩ A (n - 1)) ≤ μ (S ∩ A n) := by
        rw [Set.inter_comm (A n) (A (n - 1))]; exact hmono
      -- split the sums
      have hIco : ∑ i ∈ Finset.Ico 1 (n + 1), μ (A i ∩ A (i - 1))
          = (∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (i - 1))) + μ (A n ∩ A (n - 1)) := by
        rw [Finset.sum_Ico_succ_top hn]
      have hRange : ∑ i ∈ Finset.range (n + 1), μ (A i)
          = (∑ i ∈ Finset.range n, μ (A i)) + μ (A n) := by
        rw [Finset.sum_range_succ]
      -- assemble (all ≤ and + in ENNReal)
      rw [hUunion, hIco, hRange]
      calc μ (S ∪ A n) + ((∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (i - 1))) + μ (A n ∩ A (n - 1)))
          ≤ μ (S ∪ A n) + ((∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (i - 1))) + μ (S ∩ A n)) :=
            add_le_add (le_refl _) (add_le_add (le_refl _) htop)
        _ = (μ (S ∪ A n) + μ (S ∩ A n)) + ∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (i - 1)) := by
            ring
        _ = (μ S + μ (A n)) + ∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (i - 1)) := by rw [hkey]
        _ = (μ S + ∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (i - 1))) + μ (A n) := by ring
        _ ≤ (∑ i ∈ Finset.range n, μ (A i)) + μ (A n) := add_le_add ih (le_refl _)

/-- **The ledger coefficient**: the path-Bonferroni credit `1 − c/7 + (c−1)/49` equals
`(48 − 6c)/49`, strictly positive for every block size `c ≤ 7` — the 7-wall crossing. -/
theorem ledger_coeff (c : ℝ) : 1 - c / 7 + (c - 1) / 49 = (48 - 6 * c) / 49 := by ring

theorem ledger_coeff_pos {c : ℕ} (hc : c ≤ 7) : 0 < (48 - 6 * (c : ℝ)) / 49 := by
  have : (c : ℝ) ≤ 7 := by exact_mod_cast hc
  apply div_pos (by linarith) (by norm_num)

#print axioms path_hunter_add_le

end LonelyRunner.LRC14.Hunter
