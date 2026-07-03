/-
  TournamentH7.LRCChainPeel — THE CITE-CHAIN THEOREM (kind-pasteur-2026-07-02-S19,
  HYP-3976).  The compressed core's structural attack.

  Three kernel-pure layers over the citation node, generalizing the S18 91-peel from
  "one dominant runner" to "any ratio-3 chain over any cited base":

  * `spread7_lonely` — the DILATION WINDOW: max ≤ 7·min ⟹ lonely at `t = 1/(14W)`
    outright (every |vᵢ|·t ∈ [1/14, 1/2]).  Zero citations.

  * `good_window_step` / `good_chain` — the NESTED-WINDOW ENGINE: inside any interval
    of length L, a runner w with w·L ≥ 3/2 owns the full sub-window
    [(j + 1/14)/w, (j + 1/14 + 1/2)/w], keeping the 1/14 band with room (≥ 1/14 on
    the left of the cell, ≥ 3/7 on the right); the next runner nests whenever
    w' ≥ 3·w.  `ChainOK` is the recursive window ledger.

  * `cite_chain_lonely` — the COMPOSITION: cite any k ≤ 12 of the runners at gap
    1/(k+1); the margin transports along cited speeds ≤ B to [t₀ ± δ] with
    δ = (13−k)/(14(k+1)B) — consumed EXACTLY at radius 1/14 for every k:
    1/(k+1) − B·δ = 1/14; the remaining runners, listed as a `ChainOK` chain for the
    base length 2δ, nest inside.  k = 12 recovers a 91-peel variant; smaller k
    trades a longer chain for a fatter margin; k = 0 is the pure lacunary case.

  Remaining after this file: families where below every candidate split some
  consecutive ratio drops under 3 — the near-equal-block regime (difference-lattice
  land; the arc/measure route's domain).
-/
import Mathlib
import TournamentH7.LRC13Citation
import TournamentH7.LRC14ConcreteSurface

namespace LonelyRunner
namespace LRC14

/-! ## A. The dilation window -/

/-- **The spread-7 dilation window**: all absolute speeds in `[W, 7W]` ⟹ lonely at
`t = 1/(14W)`.  No citations, no census. -/
theorem spread7_lonely (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (W : ℤ) (hW : 0 < W)
    (hlo : ∀ i, W ≤ |v i|) (hhi : ∀ i, |v i| ≤ 7 * W) : ∃ t : ℝ, Lonely 14 v t := by
  have hWR : (0 : ℝ) < (W : ℝ) := by exact_mod_cast hW
  refine ⟨1 / (14 * W), ?_⟩
  intro i m
  have h1 : (W : ℝ) ≤ |(v i : ℝ)| := by
    rw [← Int.cast_abs]; exact_mod_cast hlo i
  have h2 : |(v i : ℝ)| ≤ 7 * (W : ℝ) := by
    rw [← Int.cast_abs]; exact_mod_cast hhi i
  have htpos : (0 : ℝ) < 1 / (14 * W) := by positivity
  have habs : |(v i : ℝ) * (1 / (14 * W))| = |(v i : ℝ)| * (1 / (14 * W)) := by
    rw [abs_mul, abs_of_pos htpos]
  have hx14 : (1 : ℝ) / 14 ≤ |(v i : ℝ)| * (1 / (14 * W)) := by
    rw [mul_one_div, le_div_iff₀ (by positivity : (0:ℝ) < 14 * W)]
    nlinarith
  have hx12 : |(v i : ℝ)| * (1 / (14 * W)) ≤ 1 / 2 := by
    rw [mul_one_div, div_le_iff₀ (by positivity : (0:ℝ) < 14 * W)]
    nlinarith
  rcases eq_or_ne m 0 with rfl | hm
  · simp only [Int.cast_zero, sub_zero]
    calc (1 : ℝ) / (14 : ℕ) = 1 / 14 := by norm_num
      _ ≤ |(v i : ℝ)| * (1 / (14 * W)) := hx14
      _ = |(v i : ℝ) * (1 / (14 * W))| := habs.symm
  · have hm1 : (1 : ℝ) ≤ |(m : ℝ)| := by
      rw [← Int.cast_abs]
      exact_mod_cast Int.one_le_abs hm
    have htri : |(m : ℝ)| - |(v i : ℝ) * (1 / (14 * W))| ≤
        |(v i : ℝ) * (1 / (14 * W)) - m| := by
      calc |(m : ℝ)| - |(v i : ℝ) * (1 / (14 * W))|
          ≤ |(m : ℝ) - (v i : ℝ) * (1 / (14 * W))| := abs_sub_abs_le_abs_sub _ _
        _ = |(v i : ℝ) * (1 / (14 * W)) - m| := abs_sub_comm _ _
    rw [habs] at htri
    calc (1 : ℝ) / (14 : ℕ) ≤ 1 / 2 := by norm_num
      _ ≤ |(v i : ℝ) * (1 / (14 * W)) - m| := by linarith

/-! ## B. The nested-window engine -/

/-- **The window step**: inside `[a, a+L]` with `w·L ≥ 3/2`, the runner `w` owns the
full sub-window `[(j + 1/14)/w, (j + 1/14 + 1/2)/w]` — every point of it keeps the
`1/14` band for `w`. -/
theorem good_window_step (w : ℤ) (hw : 0 < w) (a L : ℝ) (hL : 3 / 2 ≤ (w : ℝ) * L) :
    ∃ x : ℝ, a ≤ x ∧ x + 1 / (2 * (w : ℝ)) ≤ a + L ∧
      ∀ t : ℝ, x ≤ t → t ≤ x + 1 / (2 * (w : ℝ)) →
        ∀ m : ℤ, (1 : ℝ) / 14 ≤ |(w : ℝ) * t - m| := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  set j : ℤ := ⌈(w : ℝ) * a - 1 / 14⌉ with hj
  have hjge : (w : ℝ) * a - 1 / 14 ≤ (j : ℝ) := Int.le_ceil _
  have hjlt : (j : ℝ) < (w : ℝ) * a - 1 / 14 + 1 := Int.ceil_lt_add_one _
  refine ⟨((j : ℝ) + 1 / 14) / w, ?_, ?_, ?_⟩
  · rw [le_div_iff₀ hwR]
    nlinarith
  · have hsum : ((j : ℝ) + 1 / 14) / w + 1 / (2 * (w : ℝ))
        = ((j : ℝ) + 1 / 14 + 1 / 2) / w := by
      field_simp
    rw [hsum, div_le_iff₀ hwR]
    nlinarith
  · intro t ht1 ht2 m
    have hwt1 : (j : ℝ) + 1 / 14 ≤ (w : ℝ) * t := by
      rw [div_le_iff₀ hwR] at ht1
      have hc : t * (w : ℝ) = (w : ℝ) * t := mul_comm _ _
      linarith
    have hwt2 : (w : ℝ) * t ≤ (j : ℝ) + 1 / 14 + 1 / 2 := by
      have hxw : (w : ℝ) * (((j : ℝ) + 1 / 14) / w) = (j : ℝ) + 1 / 14 := by
        field_simp
      have h2w : (w : ℝ) * (1 / (2 * (w : ℝ))) = 1 / 2 := by
        field_simp
      have hmul := mul_le_mul_of_nonneg_left ht2 hwR.le
      rw [mul_add, hxw, h2w] at hmul
      exact hmul
    rcases le_or_gt (m : ℝ) (j : ℝ) with hm | hm
    · rw [abs_of_nonneg (by linarith)]
      calc (1 : ℝ) / 14 ≤ (w : ℝ) * t - j := by linarith
        _ ≤ (w : ℝ) * t - m := by linarith
    · have hjm : j + 1 ≤ m := by
        have : j < m := by exact_mod_cast hm
        omega
      have hjmR : (j : ℝ) + 1 ≤ (m : ℝ) := by exact_mod_cast hjm
      rw [abs_of_nonpos (by linarith)]
      linarith

/-- The recursive window ledger: the head half-covers the current window; each
successor works inside the head's `1/(2w)`-window. -/
def ChainOK : List ℤ → ℝ → Prop
  | [], _ => True
  | w :: ws, L => 3 / 2 ≤ (w : ℝ) * L ∧ ChainOK ws (1 / (2 * (w : ℝ)))

/-- **The chain engine**: a `ChainOK` list digs a common good point out of any base
interval. -/
theorem good_chain : ∀ (ws : List ℤ) (a L : ℝ), 0 ≤ L → (∀ w ∈ ws, 0 < w) →
    ChainOK ws L →
    ∃ t : ℝ, a ≤ t ∧ t ≤ a + L ∧
      ∀ w ∈ ws, ∀ m : ℤ, (1 : ℝ) / 14 ≤ |(w : ℝ) * t - m| := by
  intro ws
  induction ws with
  | nil =>
      intro a L hL _ _
      exact ⟨a, le_refl _, by linarith, fun w hw => absurd hw (List.not_mem_nil)⟩
  | cons w ws ih =>
      intro a L hL hpos hchain
      obtain ⟨hentry, hrest⟩ := hchain
      have hwpos : 0 < w := hpos w (List.mem_cons_self ..)
      have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hwpos
      obtain ⟨x, hx1, hx2, hgood⟩ := good_window_step w hwpos a L hentry
      obtain ⟨t, ht1, ht2, hws⟩ := ih x (1 / (2 * (w : ℝ))) (by positivity)
        (fun w' hw' => hpos w' (List.mem_cons_of_mem _ hw')) hrest
      refine ⟨t, by linarith, by linarith, ?_⟩
      intro w' hw' m
      rcases List.mem_cons.mp hw' with rfl | hw''
      · exact hgood t ht1 ht2 m
      · exact hws w' hw'' m

/-! ## C. The composition: cite k, chain the rest -/

/-- **THE CITE-CHAIN THEOREM.**  Cite any `k ≤ 12` of the runners (those with index
`< k`, bounded by `B`) at gap `1/(k+1)`; every other runner's absolute value sits in
a `ChainOK` chain for the transported window of length `2δ`, `δ = (13−k)/(14(k+1)B)`.
Then the family is lonely.  The margin is consumed exactly:
`1/(k+1) − B·δ = 1/14`. -/
theorem cite_chain_lonely (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (k : ℕ) (hk : k ≤ 12) (B : ℤ) (hB : 0 < B)
    (hcited : ∀ i : Fin 13, (i : ℕ) < k → |v i| ≤ B)
    (ws : List ℤ) (hwpos : ∀ w ∈ ws, 0 < w)
    (hsplit : ∀ i : Fin 13, (i : ℕ) < k ∨ |v i| ∈ ws)
    (hchain : ChainOK ws
      (2 * (((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ))))) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hk13 : k ≤ 13 := by omega
  have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB
  have hkR : (k : ℝ) ≤ 12 := by exact_mod_cast hk
  set δ : ℝ := ((13 : ℝ) - k) / (14 * ((k : ℝ) + 1) * (B : ℝ)) with hδ
  have hδpos : 0 < δ := by
    rw [hδ]
    apply div_pos (by linarith) (by positivity)
  -- cite the first k
  set wc : Fin k → ℤ := fun j => v (Fin.castLE hk13 j) with hwc
  have hwcne : ∀ j, wc j ≠ 0 := fun j => hv _
  obtain ⟨t₀, hcite⟩ := cite k hk wc hwcne
  -- chain the rest inside [t₀ − δ, t₀ + δ]
  obtain ⟨τ, hτ1, hτ2, hτgood⟩ := good_chain ws (t₀ - δ) (2 * δ) (by linarith)
    hwpos hchain
  refine ⟨τ, ?_⟩
  intro i m
  rcases hsplit i with hlt | hmem
  · -- cited: Lipschitz transport of the 1/(k+1) margin
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
  · -- chained: the good point serves |v i|, transfer through the sign
    have hgood := hτgood |v i| hmem
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
