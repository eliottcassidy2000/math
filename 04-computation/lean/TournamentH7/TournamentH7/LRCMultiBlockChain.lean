/-
  TournamentH7.LRCMultiBlockChain — MULTI-BLOCK CHAINS through the citation transport
  (death-star-2026-07-17-S37, HYP-7181; generalizes THM-941's single-block engine).

  THM-941 composed ONE kps-S22 fat block with a singles tail under the citation.  This
  module chains ARBITRARILY MANY blocks: each level `(ws, ε)` pays the fat-mass fee at
  the current window and hands the `ε`-window to the next level; the S19 singles tail
  (cheap, ratio-3 fees) runs inside the LAST window.

  * `MultiBlockOK` — the recursive ledger; `finalWidth` — the last window's width.
  * `multiblock_window` — the window induction over `block_window_step`: a compliant
    level list digs a nested final window on which EVERY block runner is good.
  * `lonely_of_multiblock_split` — the cited composition: cite `k ≤ 12` runners ≤ B,
    run the level chain in the transported `2δ`-window, finish with a `ChainOK`
    singles tail at `finalWidth`.  THM-941's `lonely_of_block_split` is the
    one-level instance; the S36 corner closures are the one-level, empty-tail,
    `ε = 0` instances.
  * `lonely_of_two_block_split` — the two-level instance, fees exposed: the first
    formal statement in the corpus in which TWO separated dense clusters are dodged
    in one certificate.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCBlockSplitLift

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner

/-- The recursive multi-block ledger: each level `(ws, ε)` pays the kps-S22 fat-mass
fee at the current window `L` and recurses at `ε`. -/
def MultiBlockOK : List (List ℤ × ℝ) → ℝ → Prop
  | [], _ => True
  | (ws, ε) :: rest, L =>
      (∀ w ∈ ws, 0 < w) ∧ 0 ≤ ε ∧
      ((ws.map fun (w : ℤ) => L / 7 + 3 / (7 * (w : ℝ)) + ε * ((w : ℝ) * L + 3))).sum
        < L - ε ∧
      MultiBlockOK rest ε

/-- The width of the last window a level list hands to its tail. -/
def finalWidth : List (List ℤ × ℝ) → ℝ → ℝ
  | [], L => L
  | (_, ε) :: rest, _ => finalWidth rest ε

theorem finalWidth_nonneg : ∀ (levels : List (List ℤ × ℝ)) (L : ℝ), 0 ≤ L →
    MultiBlockOK levels L → 0 ≤ finalWidth levels L := by
  intro levels
  induction levels with
  | nil => intro L hL _; exact hL
  | cons lv rest ih =>
    intro L _ hOK
    obtain ⟨_, hε, _, hrest⟩ := hOK
    exact ih lv.2 hε hrest

/-- **The window induction**: a compliant level list digs a nested final window of
width `finalWidth` inside `[a, a+L]` on which every block runner is good. -/
theorem multiblock_window : ∀ (levels : List (List ℤ × ℝ)) (a L : ℝ), 0 ≤ L →
    MultiBlockOK levels L →
    ∃ x : ℝ, a ≤ x ∧ x + finalWidth levels L ≤ a + L ∧
      ∀ t : ℝ, x ≤ t → t ≤ x + finalWidth levels L →
        ∀ pair ∈ levels, ∀ w ∈ pair.1, ∀ m : ℤ, (1 : ℝ) / 14 ≤ |(w : ℝ) * t - m| := by
  intro levels
  induction levels with
  | nil =>
    intro a L hL _
    exact ⟨a, le_refl _, by simp [finalWidth], fun t _ _ pair hpair =>
      absurd hpair (List.not_mem_nil)⟩
  | cons lv rest ih =>
    intro a L hL hOK
    obtain ⟨hpos, hε, hmass, hrest⟩ := hOK
    obtain ⟨x₁, hx₁a, hx₁b, hgood⟩ := LRC14.block_window_step lv.1 hpos lv.2
      a (a + L) hε (by linarith)
      (by simpa [add_sub_cancel_left] using hmass)
    obtain ⟨x, hxa, hxb, hall⟩ := ih x₁ lv.2 hε hrest
    have hfw : finalWidth (lv :: rest) L = finalWidth rest lv.2 := by
      cases lv
      rfl
    refine ⟨x, by linarith, ?_, ?_⟩
    · rw [hfw]
      linarith
    · intro t ht1 ht2 pair hpair w hw m
      rw [hfw] at ht2
      rcases List.mem_cons.mp hpair with rfl | hmem
      · exact hgood t (by linarith) (by linarith) w hw m
      · exact hall t ht1 ht2 pair hmem w hw m

/-- **The cited multi-block composition**: cite `k ≤ 12` runners ≤ B; run the level
chain in the transported `2δ`-window; finish with a `ChainOK` singles tail at
`finalWidth`.  Every runner is cited, in some level, or in the tail. -/
theorem lonely_of_multiblock_split (cite : LRCUpTo13) (v : Fin 13 → ℤ)
    (hv : ∀ i, v i ≠ 0) (k : ℕ) (hk : k ≤ 12) (B : ℤ) (hB : 0 < B)
    (hcited : ∀ i : Fin 13, (i : ℕ) < k → |v i| ≤ B)
    (levels : List (List ℤ × ℝ))
    (tail : List ℤ) (htpos : ∀ w ∈ tail, 0 < w)
    (hOK : MultiBlockOK levels
      (2 * (((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ)))))
    (htail : LRC14.ChainOK tail
      (finalWidth levels (2 * (((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ))))))
    (hsplit : ∀ i : Fin 13, (i : ℕ) < k ∨
      (∃ pair ∈ levels, |v i| ∈ pair.1) ∨ |v i| ∈ tail) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hk13 : k ≤ 13 := by omega
  have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB
  have hkR : (k : ℝ) ≤ 12 := by exact_mod_cast hk
  set δ : ℝ := ((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ)) with hδ
  have hδpos : 0 < δ := by
    rw [hδ]
    apply div_pos (by linarith) (by positivity)
  set wc : Fin k → ℤ := fun j => v (Fin.castLE hk13 j) with hwc
  have hwcne : ∀ j, wc j ≠ 0 := fun j => hv _
  obtain ⟨t₀, hcite⟩ := cite k hk wc hwcne
  obtain ⟨x, hxa, hxb, hblocks⟩ := multiblock_window levels (t₀ - δ) (2 * δ)
    (by linarith) hOK
  have hfw0 : 0 ≤ finalWidth levels (2 * δ) :=
    finalWidth_nonneg levels (2 * δ) (by linarith) hOK
  obtain ⟨τ, hτ1, hτ2, hτtail⟩ := LRC14.good_chain tail x
    (finalWidth levels (2 * δ)) hfw0 htpos htail
  refine ⟨τ, ?_⟩
  intro i m
  rcases hsplit i with hlt | ⟨pair, hpair, hmem⟩ | hmem
  · -- cited: Lipschitz transport
    have hidx : Fin.castLE hk13 ⟨(i : ℕ), hlt⟩ = i := by
      apply Fin.ext
      rfl
    have h0 : (1 : ℝ) / (k + 1 : ℕ) ≤ |(v i : ℝ) * t₀ - m| := by
      simpa [hwc, hidx] using hcite ⟨(i : ℕ), hlt⟩ m
    have hvB : |(v i : ℝ)| ≤ (B : ℝ) := by
      rw [← Int.cast_abs]
      exact_mod_cast hcited i hlt
    have htri : |(v i : ℝ) * t₀ - m| ≤ |(v i : ℝ) * τ - m| + |(v i : ℝ)| * |t₀ - τ| := by
      calc |(v i : ℝ) * t₀ - m|
          = |((v i : ℝ) * τ - m) + (v i : ℝ) * (t₀ - τ)| := by congr 1; ring
        _ ≤ |(v i : ℝ) * τ - m| + |(v i : ℝ) * (t₀ - τ)| := abs_add_le _ _
        _ = |(v i : ℝ) * τ - m| + |(v i : ℝ)| * |t₀ - τ| := by rw [abs_mul]
    have hwin : |t₀ - τ| ≤ δ := by
      rw [abs_le]
      constructor <;> linarith
    have hlip : |(v i : ℝ)| * |t₀ - τ| ≤ (B : ℝ) * δ := by
      apply mul_le_mul hvB hwin (abs_nonneg _) hBR.le
    have hBδ : (B : ℝ) * δ = ((13 : ℝ) - k) / (14 * ((k : ℝ) + 1)) := by
      rw [hδ]
      field_simp
    have hmargin : (1 : ℝ) / (k + 1 : ℕ) - ((13 : ℝ) - k) / (14 * ((k : ℝ) + 1))
        = 1 / 14 := by
      have hcast : ((k + 1 : ℕ) : ℝ) = (k : ℝ) + 1 := by push_cast; ring
      rw [hcast]
      field_simp
      ring
    calc (1 : ℝ) / (14 : ℕ) = 1 / 14 := by norm_num
      _ = 1 / (k + 1 : ℕ) - ((13 : ℝ) - k) / (14 * ((k : ℝ) + 1)) := hmargin.symm
      _ ≤ |(v i : ℝ) * τ - m| := by
          rw [hBδ] at hlip
          linarith
  · -- in a block level
    have hgood := hblocks τ hτ1 hτ2 pair hpair |v i| hmem
    rcases abs_cases ((v i : ℝ)) with ⟨habs, _⟩ | ⟨habs, _⟩
    · have hcast : ((|v i| : ℤ) : ℝ) = (v i : ℝ) := by
        rw [Int.cast_abs, habs]
      have hg := hgood m
      rw [hcast] at hg
      exact (by norm_num : (1 : ℝ) / (14 : ℕ) = 1 / 14) ▸ hg
    · have hcast : ((|v i| : ℤ) : ℝ) = -(v i : ℝ) := by
        rw [Int.cast_abs, habs]
      have hg := hgood (-m)
      rw [hcast] at hg
      have heq : |-(v i : ℝ) * τ - (-m : ℤ)| = |(v i : ℝ) * τ - m| := by
        rw [show -(v i : ℝ) * τ - ((-m : ℤ) : ℝ) = -((v i : ℝ) * τ - m) by push_cast; ring,
          abs_neg]
      rw [heq] at hg
      exact (by norm_num : (1 : ℝ) / (14 : ℕ) = 1 / 14) ▸ hg
  · -- in the tail
    have hgood := hτtail |v i| hmem
    rcases abs_cases ((v i : ℝ)) with ⟨habs, _⟩ | ⟨habs, _⟩
    · have hcast : ((|v i| : ℤ) : ℝ) = (v i : ℝ) := by
        rw [Int.cast_abs, habs]
      have hg := hgood m
      rw [hcast] at hg
      exact (by norm_num : (1 : ℝ) / (14 : ℕ) = 1 / 14) ▸ hg
    · have hcast : ((|v i| : ℤ) : ℝ) = -(v i : ℝ) := by
        rw [Int.cast_abs, habs]
      have hg := hgood (-m)
      rw [hcast] at hg
      have heq : |-(v i : ℝ) * τ - (-m : ℤ)| = |(v i : ℝ) * τ - m| := by
        rw [show -(v i : ℝ) * τ - ((-m : ℤ) : ℝ) = -((v i : ℝ) * τ - m) by push_cast; ring,
          abs_neg]
      rw [heq] at hg
      exact (by norm_num : (1 : ℝ) / (14 : ℕ) = 1 / 14) ▸ hg

/-- **The two-block instance**, fees exposed: TWO separated dense clusters dodged in
one certificate — the first such statement in the corpus. -/
theorem lonely_of_two_block_split (cite : LRCUpTo13) (v : Fin 13 → ℤ)
    (hv : ∀ i, v i ≠ 0) (k : ℕ) (hk : k ≤ 12) (B : ℤ) (hB : 0 < B)
    (hcited : ∀ i : Fin 13, (i : ℕ) < k → |v i| ≤ B)
    (ws₁ ws₂ : List ℤ) (hw₁ : ∀ w ∈ ws₁, 0 < w) (hw₂ : ∀ w ∈ ws₂, 0 < w)
    (ε₁ ε₂ : ℝ) (hε₁ : 0 ≤ ε₁) (hε₂ : 0 ≤ ε₂)
    (tail : List ℤ) (htpos : ∀ w ∈ tail, 0 < w)
    (hmass₁ : ((ws₁.map fun (w : ℤ) =>
        (2 * (((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ)))) / 7
          + 3 / (7 * (w : ℝ))
          + ε₁ * ((w : ℝ) * (2 * (((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ)))) + 3))).sum
      < 2 * (((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ))) - ε₁)
    (hmass₂ : ((ws₂.map fun (w : ℤ) =>
        ε₁ / 7 + 3 / (7 * (w : ℝ)) + ε₂ * ((w : ℝ) * ε₁ + 3))).sum < ε₁ - ε₂)
    (htail : LRC14.ChainOK tail ε₂)
    (hsplit : ∀ i : Fin 13, (i : ℕ) < k ∨ |v i| ∈ ws₁ ∨ |v i| ∈ ws₂ ∨ |v i| ∈ tail) :
    ∃ t : ℝ, Lonely 14 v t := by
  apply lonely_of_multiblock_split cite v hv k hk B hB hcited
    [(ws₁, ε₁), (ws₂, ε₂)] tail htpos
    ⟨hw₁, hε₁, hmass₁, hw₂, hε₂, hmass₂, trivial⟩
    (by
      show LRC14.ChainOK tail (finalWidth [(ws₂, ε₂)] ε₁)
      show LRC14.ChainOK tail (finalWidth [] ε₂)
      exact htail)
  intro i
  rcases hsplit i with h | h | h | h
  · exact Or.inl h
  · exact Or.inr (Or.inl ⟨(ws₁, ε₁), List.mem_cons_self .., h⟩)
  · exact Or.inr (Or.inl ⟨(ws₂, ε₂), List.mem_cons_of_mem _ (List.mem_singleton.mpr rfl), h⟩)
  · exact Or.inr (Or.inr h)

/-! ## Axiom audit -/
#print axioms multiblock_window
#print axioms lonely_of_multiblock_split
#print axioms lonely_of_two_block_split

end LRC14Grand
end LonelyRunner
