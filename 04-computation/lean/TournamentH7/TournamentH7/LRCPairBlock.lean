/-
  TournamentH7.LRCPairBlock — THE PAIR DODGE and the BLOCK CHAIN
  (kind-pasteur-2026-07-02-S20, HYP-3977).  Near-equal-block attack, stage 1.

  After S19, the open core is families with a near-equal block (consecutive ratio < 3)
  somewhere above every split.  This file breaks the SMALLEST blocks: PAIRS.

  THE PAIR DODGE (`pair_window_step`): inside any interval `[a, a+L]` with
  `w·L ≥ 13/7`, for a pair `w ≤ w' < 3w`, one of the three QUARTER CANDIDATES
  `τ_k = (j + k/4)/w` (k = 1, 2, 3) has `‖w'·τ_k‖ ≥ 1/7`.  Reason: the candidate
  positions for `w'` are spaced `Δ = w'/(4w) ∈ [1/4, 3/4)` apart; if two candidates
  were both within `1/7` of ℤ their spacing would be within `2/7` of an integer —
  and `Δ ∈ [1/4, 2/7) ∪ (5/7, 3/4)` forces `‖2Δ‖ ≥ 3/7 > 2/7`, so all three bad is
  impossible.  The surviving candidate's sub-window of length `1/(14·w')` keeps BOTH
  margins ≥ 1/14 (`w`: 1/4 − 1/14 ≥ 1/14; `w'`: 1/7 − 1/14 = 1/14, exactly tight).

  THE BLOCK CHAIN (`BLevel`/`BChainOK`/`block_chain`): chain levels are now singles
  (S19 window step, entry `w·L ≥ 3/2`, output `1/(2w)`) or pairs (the dodge, entry
  `w·L ≥ 13/7`, output `1/(14w')`), composed by the same nesting induction; and
  `cite_blockchain_lonely` re-runs the S19 citation composition over block levels.

  Remaining after this file: blocks of THREE-plus near-equal runners above every
  split — the residual core, now visibly thinner.
-/
import Mathlib
import TournamentH7.LRCChainPeel

namespace LonelyRunner
namespace LRC14

/-! ## The pair dodge -/

/-- Two integers within `1/7` of two multiples of a real transfer their spacing:
helper extracting the near-integrality of `Δ` from two bad candidates. -/
theorem near_int_of_two_bad {y₁ y₂ : ℝ}
    (h1 : ∃ m : ℤ, |y₁ - m| < 1 / 7) (h2 : ∃ m : ℤ, |y₂ - m| < 1 / 7) :
    ∃ M : ℤ, |(y₂ - y₁) - M| < 2 / 7 := by
  obtain ⟨m₁, hm₁⟩ := h1
  obtain ⟨m₂, hm₂⟩ := h2
  refine ⟨m₂ - m₁, ?_⟩
  have : (y₂ - y₁) - ((m₂ : ℝ) - m₁) = (y₂ - m₂) - (y₁ - m₁) := by ring
  rw [show ((m₂ - m₁ : ℤ) : ℝ) = (m₂ : ℝ) - m₁ by push_cast; ring, this]
  calc |(y₂ - m₂) - (y₁ - m₁)| ≤ |y₂ - m₂| + |y₁ - m₁| := abs_sub _ _
    _ < 2 / 7 := by linarith

/-- **THE PAIR DODGE.**  Inside `[a, a+L]` with `w·L ≥ 13/7`, a near-equal pair
`w ≤ w' < 3w` has a common sub-window of length `1/(14w')` on which both runners
keep the `1/14` band. -/
theorem pair_window_step (w w' : ℤ) (hw : 0 < w) (hww' : w ≤ w') (hw3 : w' < 3 * w)
    (a L : ℝ) (hL : 13 / 7 ≤ (w : ℝ) * L) :
    ∃ x : ℝ, a ≤ x ∧ x + 1 / (14 * (w' : ℝ)) ≤ a + L ∧
      ∀ t : ℝ, x ≤ t → t ≤ x + 1 / (14 * (w' : ℝ)) →
        (∀ m : ℤ, (1 : ℝ) / 14 ≤ |(w : ℝ) * t - m|) ∧
        (∀ m : ℤ, (1 : ℝ) / 14 ≤ |(w' : ℝ) * t - m|) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  have hw'R : (0 : ℝ) < (w' : ℝ) := by
    have : 0 < w' := lt_of_lt_of_le hw hww'
    exact_mod_cast this
  have hwwR : (w : ℝ) ≤ (w' : ℝ) := by exact_mod_cast hww'
  have hw3R : (w' : ℝ) < 3 * (w : ℝ) := by exact_mod_cast hw3
  set j : ℤ := ⌈(w : ℝ) * a⌉ with hj
  have hjge : (w : ℝ) * a ≤ (j : ℝ) := Int.le_ceil _
  have hjlt : (j : ℝ) < (w : ℝ) * a + 1 := Int.ceil_lt_add_one _
  -- the three quarter candidates
  set τ : ℕ → ℝ := fun k => ((j : ℝ) + k / 4) / w with hτ
  -- the dodge: some k ∈ {1,2,3} has ‖w'·τ k‖ ≥ 1/7
  have hdodge : ∃ k : ℕ, 1 ≤ k ∧ k ≤ 3 ∧ ∀ m : ℤ, 1 / 7 ≤ |(w' : ℝ) * τ k - m| := by
    by_contra hcon
    push Not at hcon
    -- all three candidates are bad
    have hbad : ∀ k : ℕ, 1 ≤ k → k ≤ 3 → ∃ m : ℤ, |(w' : ℝ) * τ k - m| < 1 / 7 := by
      intro k h1 h3
      obtain ⟨m, hm⟩ := hcon k h1 h3
      exact ⟨m, by linarith⟩
    set Δ : ℝ := (w' : ℝ) / (4 * w) with hΔ
    have hΔval : ∀ k : ℕ, (w' : ℝ) * τ (k + 1) - (w' : ℝ) * τ k = Δ := by
      intro k
      rw [hτ, hΔ]
      push_cast
      field_simp
      ring
    have hΔlo : (1 : ℝ) / 4 ≤ Δ := by
      rw [hΔ, le_div_iff₀ (by positivity)]
      linarith
    have hΔhi : Δ < 3 / 4 := by
      rw [hΔ, div_lt_iff₀ (by positivity)]
      linarith
    -- spacing transfers
    obtain ⟨M1, hM1⟩ := near_int_of_two_bad (hbad 1 (by norm_num) (by norm_num))
      (hbad 2 (by norm_num) (by norm_num))
    obtain ⟨M2, hM2⟩ := near_int_of_two_bad (hbad 1 (by norm_num) (by norm_num))
      (hbad 3 (by norm_num) (by norm_num))
    have e1 : (w' : ℝ) * τ 2 - (w' : ℝ) * τ 1 = Δ := hΔval 1
    have e2 : (w' : ℝ) * τ 3 - (w' : ℝ) * τ 1 = 2 * Δ := by
      have h12 := hΔval 1
      have h23 := hΔval 2
      linarith
    rw [e1] at hM1
    rw [e2] at hM2
    -- M1 ∈ {0, 1}
    have hM1cases : M1 = 0 ∨ M1 = 1 := by
      have hlo : -(2 : ℝ) / 7 < (M1 : ℝ) - Δ := by
        have := abs_lt.mp hM1
        linarith [this.1]
      have hhi : (M1 : ℝ) - Δ < 2 / 7 := by
        have := abs_lt.mp hM1
        linarith [this.2]
      have h1 : (-1 : ℝ) < (M1 : ℝ) := by linarith
      have h2 : (M1 : ℝ) < 2 := by linarith
      have h1' : (-1 : ℤ) < M1 := by exact_mod_cast h1
      have h2' : M1 < 2 := by exact_mod_cast h2
      omega
    have habs1 := abs_lt.mp hM1
    have habs2 := abs_lt.mp hM2
    rcases hM1cases with rfl | rfl
    · -- Δ < 2/7, so 2Δ ∈ [1/2, 4/7): |2Δ − M2| < 2/7 forces M2 ∈ ∅
      simp only [Int.cast_zero] at habs1
      have hΔsmall : Δ < 2 / 7 := by linarith [habs1.2]
      have h2Δlo : (1 : ℝ) / 2 ≤ 2 * Δ := by linarith
      have h2Δhi : 2 * Δ < 4 / 7 := by linarith
      have hM2lo : (2 : ℝ) / 7 + -(2/7) < (M2 : ℝ) := by linarith [habs2.1]
      have hM2hi : (M2 : ℝ) < 4 / 7 + 2 / 7 := by linarith [habs2.2]
      have h0 : (0 : ℝ) < (M2 : ℝ) := by linarith
      have h1 : (M2 : ℝ) < 1 := by linarith
      have h0' : (0 : ℤ) < M2 := by exact_mod_cast h0
      have h1' : M2 < 1 := by exact_mod_cast h1
      omega
    · -- Δ > 5/7, so 2Δ ∈ (10/7, 3/2): |2Δ − M2| < 2/7 forces M2 ∈ ∅
      simp only [Int.cast_one] at habs1
      have hΔbig : 5 / 7 < Δ := by linarith [habs1.1]
      have h2Δlo : (10 : ℝ) / 7 < 2 * Δ := by linarith
      have h2Δhi : 2 * Δ < 3 / 2 := by linarith
      have hM2lo : (10 : ℝ) / 7 - 2 / 7 < (M2 : ℝ) := by linarith [habs2.1]
      have hM2hi : (M2 : ℝ) < 3 / 2 + 2 / 7 := by linarith [habs2.2]
      have h1 : (1 : ℝ) < (M2 : ℝ) := by linarith
      have h2 : (M2 : ℝ) < 2 := by linarith
      have h1' : (1 : ℤ) < M2 := by exact_mod_cast h1
      have h2' : M2 < 2 := by exact_mod_cast h2
      omega
  obtain ⟨k, hk1, hk3, hkgood⟩ := hdodge
  refine ⟨τ k, ?_, ?_, ?_⟩
  · -- a ≤ τ k
    rw [hτ]
    simp only
    rw [le_div_iff₀ hwR]
    have hk0 : (0 : ℝ) ≤ (k : ℝ) / 4 := by positivity
    have hc : a * (w : ℝ) = (w : ℝ) * a := mul_comm _ _
    linarith
  · -- τ k + 1/(14 w') ≤ a + L
    rw [hτ]
    simp only
    have hstep : ((j : ℝ) + (k : ℝ) / 4) / w + 1 / (14 * (w' : ℝ))
        ≤ ((j : ℝ) + (k : ℝ) / 4) / w + 1 / (14 * (w : ℝ)) := by
      have h1 : (0 : ℝ) < 14 * (w : ℝ) := by positivity
      have h2 : (0 : ℝ) < 14 * (w' : ℝ) := by positivity
      have : 14 * (w : ℝ) ≤ 14 * (w' : ℝ) := by linarith
      have := one_div_le_one_div_of_le h1 this
      linarith
    have hkR : (k : ℝ) ≤ 3 := by exact_mod_cast hk3
    have hmain : ((j : ℝ) + (k : ℝ) / 4) / w + 1 / (14 * (w : ℝ)) ≤ a + L := by
      have hsum : ((j : ℝ) + (k : ℝ) / 4) / w + 1 / (14 * (w : ℝ))
          = ((j : ℝ) + (k : ℝ) / 4 + 1 / 14) / w := by
        field_simp
      rw [hsum, div_le_iff₀ hwR]
      have hc : (a + L) * (w : ℝ) = (w : ℝ) * a + (w : ℝ) * L := by ring
      nlinarith [hjlt, hL, hkR]
    linarith
  · intro t ht1 ht2
    have hτk : (w : ℝ) * τ k = (j : ℝ) + (k : ℝ) / 4 := by
      rw [hτ]
      simp only
      field_simp
    have hwt1 : (j : ℝ) + (k : ℝ) / 4 ≤ (w : ℝ) * t := by
      have := mul_le_mul_of_nonneg_left ht1 hwR.le
      rwa [hτk] at this
    have hwt2 : (w : ℝ) * t ≤ (j : ℝ) + (k : ℝ) / 4 + 1 / 14 := by
      have := mul_le_mul_of_nonneg_left ht2 hwR.le
      rw [mul_add, hτk] at this
      have hww14 : (w : ℝ) * (1 / (14 * (w' : ℝ))) ≤ 1 / 14 := by
        rw [mul_one_div, div_le_div_iff₀ (by positivity) (by norm_num)]
        linarith
      linarith
    constructor
    · -- the w margin: w·t ∈ [j + k/4, j + k/4 + 1/14], k/4 ∈ [1/4, 3/4]
      intro m
      have hk1R : (1 : ℝ) ≤ (k : ℝ) := by exact_mod_cast hk1
      have hk3R : (k : ℝ) ≤ 3 := by exact_mod_cast hk3
      rcases le_or_gt (m : ℝ) (j : ℝ) with hm | hm
      · rw [abs_of_nonneg (by linarith)]
        linarith
      · have hjm : j + 1 ≤ m := by
          have : j < m := by exact_mod_cast hm
          omega
        have hjmR : (j : ℝ) + 1 ≤ (m : ℝ) := by exact_mod_cast hjm
        rw [abs_of_nonpos (by linarith)]
        linarith
    · -- the w' margin: |w'·t − w'·τk| ≤ w'/(14w') = 1/14, center margin ≥ 1/7
      intro m
      have hcenter := hkgood m
      have hdist : |(w' : ℝ) * t - (w' : ℝ) * τ k| ≤ 1 / 14 := by
        rw [show (w' : ℝ) * t - (w' : ℝ) * τ k = (w' : ℝ) * (t - τ k) by ring, abs_mul,
          abs_of_pos hw'R]
        have habs : |t - τ k| ≤ 1 / (14 * (w' : ℝ)) := by
          rw [abs_le]
          constructor <;> linarith
        calc (w' : ℝ) * |t - τ k| ≤ (w' : ℝ) * (1 / (14 * (w' : ℝ))) :=
              mul_le_mul_of_nonneg_left habs hw'R.le
          _ = 1 / 14 := by field_simp
      calc (1 : ℝ) / 14 = 1 / 7 - 1 / 14 := by norm_num
        _ ≤ |(w' : ℝ) * τ k - m| - |(w' : ℝ) * t - (w' : ℝ) * τ k| := by linarith
        _ ≤ |(w' : ℝ) * t - m| := by
            have := abs_sub_abs_le_abs_sub ((w' : ℝ) * τ k - m) ((w' : ℝ) * τ k - (w' : ℝ) * t)
            have heq : (w' : ℝ) * τ k - m - ((w' : ℝ) * τ k - (w' : ℝ) * t)
                = (w' : ℝ) * t - m := by ring
            rw [heq] at this
            have hcomm : |(w' : ℝ) * τ k - (w' : ℝ) * t| = |(w' : ℝ) * t - (w' : ℝ) * τ k| :=
              abs_sub_comm _ _
            linarith

/-! ## The block chain -/

/-- A chain level: a single runner or a near-equal pair. -/
inductive BLevel
  | single (w : ℤ)
  | pair (w w' : ℤ)

/-- The speeds a level certifies. -/
def BLevel.speeds : BLevel → List ℤ
  | .single w => [w]
  | .pair w w' => [w, w']

/-- Positivity and pair-shape side conditions. -/
def BLevel.OK : BLevel → Prop
  | .single w => 0 < w
  | .pair w w' => 0 < w ∧ w ≤ w' ∧ w' < 3 * w

/-- The output window length of a level. -/
noncomputable def BLevel.out : BLevel → ℝ
  | .single w => 1 / (2 * (w : ℝ))
  | .pair _ w' => 1 / (14 * (w' : ℝ))

/-- The entry condition of a level at window length `L`. -/
def BLevel.entry : BLevel → ℝ → Prop
  | .single w, L => 3 / 2 ≤ (w : ℝ) * L
  | .pair w _, L => 13 / 7 ≤ (w : ℝ) * L

/-- The block-chain ledger. -/
def BChainOK : List BLevel → ℝ → Prop
  | [], _ => True
  | ℓ :: ℓs, L => ℓ.OK ∧ ℓ.entry L ∧ BChainOK ℓs ℓ.out

/-- **The block-chain engine**: singles and pairs nest alike. -/
theorem block_chain : ∀ (ℓs : List BLevel) (a L : ℝ), 0 ≤ L → BChainOK ℓs L →
    ∃ t : ℝ, a ≤ t ∧ t ≤ a + L ∧
      ∀ ℓ ∈ ℓs, ∀ w ∈ ℓ.speeds, ∀ m : ℤ, (1 : ℝ) / 14 ≤ |(w : ℝ) * t - m| := by
  intro ℓs
  induction ℓs with
  | nil =>
      intro a L hL _
      exact ⟨a, le_refl _, by linarith, fun ℓ hℓ => absurd hℓ (List.not_mem_nil)⟩
  | cons ℓ ℓs ih =>
      intro a L hL hchain
      obtain ⟨hOK, hentry, hrest⟩ := hchain
      cases ℓ with
      | single w =>
          have hw : 0 < w := hOK
          have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
          obtain ⟨x, hx1, hx2, hgood⟩ := good_window_step w hw a L hentry
          obtain ⟨t, ht1, ht2, hℓs⟩ := ih x (BLevel.out (.single w)) (by
            unfold BLevel.out
            positivity) hrest
          refine ⟨t, by linarith, ?_, ?_⟩
          · have : t ≤ x + 1 / (2 * (w : ℝ)) := ht2
            linarith
          · intro ℓ' hℓ' w' hw' m
            rcases List.mem_cons.mp hℓ' with rfl | hmem
            · rcases List.mem_singleton.mp hw' with rfl
              exact hgood t ht1 ht2 m
            · exact hℓs ℓ' hmem w' hw' m
      | pair w w' =>
          obtain ⟨hw, hww', hw3⟩ := hOK
          have hw'R : (0 : ℝ) < (w' : ℝ) := by
            have : 0 < w' := lt_of_lt_of_le hw hww'
            exact_mod_cast this
          obtain ⟨x, hx1, hx2, hgood⟩ := pair_window_step w w' hw hww' hw3 a L hentry
          obtain ⟨t, ht1, ht2, hℓs⟩ := ih x (BLevel.out (.pair w w')) (by
            unfold BLevel.out
            positivity) hrest
          refine ⟨t, by linarith, ?_, ?_⟩
          · have : t ≤ x + 1 / (14 * (w' : ℝ)) := ht2
            linarith
          · intro ℓ' hℓ' u hu m
            rcases List.mem_cons.mp hℓ' with rfl | hmem
            · have hg := hgood t ht1 ht2
              rcases List.mem_cons.mp hu with rfl | hu'
              · exact hg.1 m
              · rcases List.mem_singleton.mp hu' with rfl
                exact hg.2 m
            · exact hℓs ℓ' hmem u hu m

/-- **THE CITE–BLOCK-CHAIN THEOREM**: the S19 composition with pair levels.  Cite any
`k ≤ 12` runners bounded by `B` at gap `1/(k+1)`; every other runner's absolute value
appears in some level of a `BChainOK` block chain for the transported window. -/
theorem cite_blockchain_lonely (cite : LRCUpTo13) (v : Fin 13 → ℤ)
    (hv : ∀ i, v i ≠ 0) (k : ℕ) (hk : k ≤ 12) (B : ℤ) (hB : 0 < B)
    (hcited : ∀ i : Fin 13, (i : ℕ) < k → |v i| ≤ B)
    (ℓs : List BLevel)
    (hsplit : ∀ i : Fin 13, (i : ℕ) < k ∨ ∃ ℓ ∈ ℓs, |v i| ∈ ℓ.speeds)
    (hchain : BChainOK ℓs
      (2 * (((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ))))) :
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
  obtain ⟨τ, hτ1, hτ2, hτgood⟩ := block_chain ℓs (t₀ - δ) (2 * δ) (by linarith) hchain
  refine ⟨τ, ?_⟩
  intro i m
  rcases hsplit i with hlt | ⟨ℓ, hℓ, hmem⟩
  · -- cited leg (verbatim S19)
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
  · -- chained leg via the level containing |v i|
    have hgood := hτgood ℓ hℓ |v i| hmem
    rcases abs_cases ((v i : ℝ)) with ⟨habs, _⟩ | ⟨habs, _⟩
    · have hcast : ((|v i| : ℤ) : ℝ) = (v i : ℝ) := by
        rw [Int.cast_abs, habs]
      have := hgood m
      rw [hcast] at this
      exact (by norm_num : (1:ℝ)/(14:ℕ) = 1/14) ▸ this
    · have hcast : ((|v i| : ℤ) : ℝ) = -(v i : ℝ) := by
        rw [Int.cast_abs, habs]
      have := hgood (-m)
      rw [hcast] at this
      have heq : |-(v i : ℝ) * τ - (-m : ℤ)| = |(v i : ℝ) * τ - m| := by
        rw [show -(v i : ℝ) * τ - ((-m : ℤ) : ℝ) = -((v i : ℝ) * τ - m) by push_cast; ring,
          abs_neg]
      rw [heq] at this
      exact (by norm_num : (1:ℝ)/(14:ℕ) = 1/14) ▸ this

end LRC14
end LonelyRunner
