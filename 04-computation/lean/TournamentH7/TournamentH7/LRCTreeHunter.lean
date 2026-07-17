/-
  TournamentH7.LRCTreeHunter — THE TREE-HUNTER INEQUALITY AND THE
  CONSECUTIVE-OVERLAP CLOSED FORM (boxeph-2026-07-17-S72; LEM-042/043/044).

  (1) `tree_hunter_add_le` — Hunter's second-order Bonferroni over an ARBITRARY
  rooted spanning tree, encoded by parent pointers `p : ℕ → ℕ` with
  `p i < i` for `i ≥ 1` (every finite rooted tree admits such a labelling by
  peeling leaves; in this encoding the top index is always a leaf):

      μ(⋃_{i<n} Aᵢ) + Σ_{i=1}^{n-1} μ(Aᵢ ∩ A_{p i}) ≤ Σ_{i<n} μ(Aᵢ).

  Subsumes klein's `path_hunter_add_le` (`p i = i-1`) and `star_hunter_add_le`
  (`p i = 0`) — same disjointification, with `A_{p n} ⊆ ⋃_{j<n}` the only new
  step.  This is the formal engine of the c ≤ 12 pair-credit route (LEM-044:
  max-weight spanning trees; the edge cap 1/14 makes c = 13 unreachable).

  (2) `muNum` — the LEM-042 pair-overlap integer sum: `muNum a b` equals
  `μ(D_a ∩ D_b) · 14ab/g` in units of `g/(14ab)` (danger sets at margin 1/14).
  `consecutive_form_upto_63` verifies IN-KERNEL the LEM-044 closed form
      49·muNum k (k+1) = 14·k·(k+1) + 14·(k mod 7)·(6 − k mod 7)
  for two full residue periods (k = 1..63); the general k is the named
  induction.

  (3) `window7_unique_zero` — the c = 8 pigeonhole core: every window of seven
  consecutive integers contains exactly one multiple of 7 (so exactly one
  vanishing credit excess — the strict-positivity heart of LEM-044(B)).

  Kernel-pure: no `native_decide`, no `sorry`.
-/
import Mathlib

open MeasureTheory Finset

namespace LonelyRunner.LRC14.Hunter

variable {α : Type*} [MeasurableSpace α]

/-- **The tree-Hunter inequality** (additive form): second-order Bonferroni over an
arbitrary parent-pointer tree (`p i < i` for `i ≥ 1`). -/
theorem tree_hunter_add_le (μ : Measure α) (A : ℕ → Set α) (p : ℕ → ℕ)
    (hp : ∀ i, 1 ≤ i → p i < i)
    (hA : ∀ i, MeasurableSet (A i)) (n : ℕ) :
    μ (⋃ i ∈ Finset.range n, A i) + ∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (p i))
      ≤ ∑ i ∈ Finset.range n, μ (A i) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rcases Nat.eq_zero_or_pos n with hn | hn
    · subst hn; simp
    · have hUunion : (⋃ i ∈ Finset.range (n + 1), A i)
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
      have hkey : μ (S ∪ A n) + μ (S ∩ A n) = μ S + μ (A n) :=
        measure_union_add_inter S (hA n)
      have hsub : A (p n) ∩ A n ⊆ S ∩ A n := by
        apply Set.inter_subset_inter_left
        rw [hS]
        exact Set.subset_biUnion_of_mem (Finset.mem_range.mpr (hp n hn))
      have htop : μ (A n ∩ A (p n)) ≤ μ (S ∩ A n) := by
        rw [Set.inter_comm (A n) (A (p n))]
        exact measure_mono hsub
      have hIco : ∑ i ∈ Finset.Ico 1 (n + 1), μ (A i ∩ A (p i))
          = (∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (p i))) + μ (A n ∩ A (p n)) := by
        rw [Finset.sum_Ico_succ_top hn]
      have hRange : ∑ i ∈ Finset.range (n + 1), μ (A i)
          = (∑ i ∈ Finset.range n, μ (A i)) + μ (A n) := by
        rw [Finset.sum_range_succ]
      rw [hUunion, hIco, hRange]
      calc μ (S ∪ A n) + ((∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (p i))) + μ (A n ∩ A (p n)))
          ≤ μ (S ∪ A n) + ((∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (p i))) + μ (S ∩ A n)) :=
            add_le_add (le_refl _) (add_le_add (le_refl _) htop)
        _ = (μ (S ∪ A n) + μ (S ∩ A n)) + ∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (p i)) := by
            ring
        _ = (μ S + μ (A n)) + ∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (p i)) := by rw [hkey]
        _ = (μ S + ∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (p i))) + μ (A n) := by ring
        _ ≤ (∑ i ∈ Finset.range n, μ (A i)) + μ (A n) := add_le_add ih (le_refl _)

/-- The path case is the tree with `p i = i - 1`. -/
theorem path_of_tree (μ : Measure α) (A : ℕ → Set α)
    (hA : ∀ i, MeasurableSet (A i)) (n : ℕ) :
    μ (⋃ i ∈ Finset.range n, A i) + ∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A (i - 1))
      ≤ ∑ i ∈ Finset.range n, μ (A i) :=
  tree_hunter_add_le μ A (fun i => i - 1) (fun _ hi => Nat.sub_lt hi Nat.one_pos) hA n

/-- The star case is the tree with `p i = 0`. -/
theorem star_of_tree (μ : Measure α) (A : ℕ → Set α)
    (hA : ∀ i, MeasurableSet (A i)) (n : ℕ) :
    μ (⋃ i ∈ Finset.range n, A i) + ∑ i ∈ Finset.Ico 1 n, μ (A i ∩ A 0)
      ≤ ∑ i ∈ Finset.range n, μ (A i) :=
  tree_hunter_add_le μ A (fun _ => 0) (fun _ hi => hi) hA n

/-! ## The LEM-042 pair-overlap integer sum and the LEM-044 closed form -/

/-- `muNum a b` = the pair-overlap sum in units of `gcd(a,b)/(14ab)`:
`μ(D_a ∩ D_b)·14ab/gcd = Σ_j min((a+b − 14·gcd·j)₊, 2·min a b)` with the `j = 0`
term counted once and `j ≥ 1` twice.  Here specialized to coprime pairs
(`gcd = 1`), which is all the closed form needs. -/
def muNum (a b : ℕ) : ℕ :=
  ∑ j ∈ Finset.range ((a + b) / 14 + 1),
    (if j = 0 then min (a + b) (2 * min a b)
     else 2 * min (a + b - 14 * j) (2 * min a b))

/-- **The consecutive closed form, in-kernel for two full residue periods**
(LEM-044(A)): `49·μ·14k(k+1)-units = 14k(k+1) + 14r(6−r)`, `r = k mod 7`. -/
theorem consecutive_form_upto_63 :
    ∀ k ∈ Finset.Icc 1 63,
      49 * muNum k (k + 1) = 14 * k * (k + 1) + 14 * (k % 7) * (6 - k % 7) := by
  decide

/-- Spot values: `μ(D_6 ∩ D_7) = 1/49` and `μ(D_1 ∩ D_2) = 1/14` in their units
(`muNum 6 7 = 12` means `12/(14·42) = 1/49`; `muNum 1 2 = 2` means `2/28 = 1/14`). -/
theorem muNum_6_7 : muNum 6 7 = 12 := by decide
theorem muNum_1_2 : muNum 1 2 = 2 := by decide

/-! ## The c = 8 pigeonhole core (LEM-044(B)) -/

/-- Every window of seven consecutive integers contains exactly one multiple of 7 —
the unique vanishing credit excess in the c = 8 consecutive theorem. -/
theorem window7_unique_zero (v : ℕ) :
    ∃! i, i < 7 ∧ (v + i) % 7 = 0 := by
  refine ⟨(7 - v % 7) % 7, ⟨?_, ?_⟩, ?_⟩
  · omega
  · omega
  · rintro j ⟨hj, hz⟩
    omega

#print axioms tree_hunter_add_le
#print axioms consecutive_form_upto_63
#print axioms window7_unique_zero

end LonelyRunner.LRC14.Hunter
