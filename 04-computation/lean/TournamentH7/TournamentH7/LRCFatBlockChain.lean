/-
  TournamentH7.LRCFatBlockChain — THE WINDOWED BLOCK STEP and THE UNIFIED LEVEL CHAIN
  (kind-pasteur-2026-07-02-S22, HYP-3979).

  S21's `block_point_step` yields a POINT, so blocks could only close a chain.  This
  file upgrades it to a WINDOW and unifies the S19–S21 level zoo:

  * `teeth_count_le` — the tooth count over ANY window of length `L` is ≤ `wL + 3`
    (position-free).

  * `fatTeeth` / `fatTeeth_mass` / `window_avoids_of_avoid_fat` — fatten every tooth
    to the LEFT by `ε`; a point avoiding fat teeth owns the whole `[t, t+ε]` window;
    the fat mass costs `ε·(wL + 3)` extra, position-free.

  * `block_window_step` — ≤ 6 runners (arbitrary internal ratios) own a common good
    SUB-WINDOW of any width `ε` whose fat mass fits.  Blocks now sit ANYWHERE in a
    chain.

  * `GLevel`/`GChainOK`/`glevel_chain`/`cite_glevel_lonely` — ONE level constructor
    `(ws, out)`: singles, pairs, and ≤6-blocks unified; the ledger is a pure
    position-free arithmetic check per level.

  RESIDUAL (the 7-wall, precisely): levels of 7+ runners fail every fat-mass check
  (union density ≥ 1).  The designed crossing — the interface pinned for mac-mini's
  `JointRateCore` — is the PATH-BONFERRONI (Hunter) ledger: sequential subtraction
  with pair credits gives `good ≥ |I|·(1 − c/7 + (c−1)/49) − fees`, POSITIVE AT
  `c = 7` (49 − 7·7 + 6 = 6 > 0), boundary at `c = 8`; the pair-overlap floors
  `|I ∩ Dᵢ ∩ Dᵢ₊₁| ≥ |I|/49 − err` are exactly the joint-rate per-cell obligation.
-/
import Mathlib
import TournamentH7.LRCBlockSix

namespace LonelyRunner
namespace LRC14

/-! ## Position-free tooth counting -/

/-- The tooth count of runner `w` over `[a, b]` is at most `w·(b−a) + 3`,
independently of the window's position. -/
theorem teeth_count_le {w : ℤ} (hw : 0 < w) (a b : ℝ) (hab : a ≤ b) :
    (((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℕ) : ℝ)
      ≤ (w : ℝ) * (b - a) + 3 := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  have h1 : (⌊(w : ℝ) * b⌋ : ℝ) ≤ (w : ℝ) * b := Int.floor_le _
  have h2 : (w : ℝ) * a ≤ (⌈(w : ℝ) * a⌉ : ℝ) := Int.le_ceil _
  have hmax : ((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℤ)
      = max (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1) 0 := Int.toNat_eq_max _
  have hcast : (((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℕ) : ℝ)
      = ((max (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1) 0 : ℤ) : ℝ) := by
    rw [← hmax]
    push_cast
    ring
  rw [hcast]
  rcases max_cases (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1) (0 : ℤ)
    with ⟨he, _⟩ | ⟨he, _⟩
  · rw [he]
    push_cast
    linarith
  · rw [he]
    push_cast
    nlinarith

/-! ## Fat teeth: window-valued avoidance -/

/-- The teeth of runner `w`, each fattened to the left by `ε`. -/
noncomputable def fatTeeth (w : ℤ) (ε a b : ℝ) : List (ℝ × ℝ) :=
  (teeth w a b).map fun p => (p.1 - ε, p.2)

/-- Fat mass, position-free: `(b−a)/7 + 3/(7w) + ε·(w(b−a) + 3)`. -/
theorem fatTeeth_mass {w : ℤ} (hw : 0 < w) (ε a b : ℝ) (hε : 0 ≤ ε) (hab : a ≤ b) :
    ((fatTeeth w ε a b).map fun p => clipLen p a b).sum
      ≤ (b - a) / 7 + 3 / (7 * (w : ℝ)) + ε * ((w : ℝ) * (b - a) + 3) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  have hpoint : ∀ p : ℝ × ℝ, clipLen (p.1 - ε, p.2) a b ≤ clipLen p a b + ε := by
    intro p
    unfold clipLen
    simp only
    rcases le_total (min p.2 b - max (p.1 - ε) a) 0 with h1 | h1
    · rw [max_eq_left h1]
      have h2 : (0 : ℝ) ≤ max 0 (min p.2 b - max p.1 a) := le_max_left _ _
      linarith
    · rw [max_eq_right h1]
      have h2 : max p.1 a - ε ≤ max (p.1 - ε) a := by
        rcases max_cases p.1 a with ⟨he, _⟩ | ⟨he, _⟩
        · rw [he]
          exact le_max_left _ _
        · rw [he]
          calc a - ε ≤ a := by linarith
            _ ≤ max (p.1 - ε) a := le_max_right _ _
      rcases le_total (min p.2 b - max p.1 a) 0 with h3 | h3
      · rw [max_eq_left h3]
        linarith
      · rw [max_eq_right h3]
        linarith
  have hsum : ((fatTeeth w ε a b).map fun p => clipLen p a b).sum
      ≤ ((teeth w a b).map fun p => clipLen p a b).sum + (teeth w a b).length * ε := by
    unfold fatTeeth
    rw [List.map_map]
    have hadd : ((teeth w a b).map fun p => clipLen ((fun p : ℝ × ℝ => (p.1 - ε, p.2)) p) a b).sum
        ≤ ((teeth w a b).map fun p => clipLen p a b + ε).sum := by
      apply List.sum_le_sum
      intro p _
      exact hpoint p
    have hsplit : ((teeth w a b).map fun p => clipLen p a b + ε).sum
        = ((teeth w a b).map fun p => clipLen p a b).sum + (teeth w a b).length * ε := by
      induction teeth w a b with
      | nil => simp
      | cons p ps ihp =>
          simp only [List.map_cons, List.sum_cons, List.length_cons] at ihp ⊢
          rw [ihp]
          push_cast
          ring
    calc ((teeth w a b).map fun p =>
            clipLen ((fun p : ℝ × ℝ => (p.1 - ε, p.2)) p) a b).sum
        ≤ ((teeth w a b).map fun p => clipLen p a b + ε).sum := hadd
      _ = ((teeth w a b).map fun p => clipLen p a b).sum
          + (teeth w a b).length * ε := hsplit
  have hlen : ((teeth w a b).length : ℝ)
      = (((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℕ) : ℝ) := by
    unfold teeth
    rw [List.length_map, List.length_range]
  have hcount := teeth_count_le hw a b hab
  have hεcount : (teeth w a b).length * ε ≤ ε * ((w : ℝ) * (b - a) + 3) := by
    rw [hlen]
    calc (((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℕ) : ℝ) * ε
        ≤ ((w : ℝ) * (b - a) + 3) * ε := by
          apply mul_le_mul_of_nonneg_right hcount hε
      _ = ε * ((w : ℝ) * (b - a) + 3) := mul_comm _ _
  have hteeth := teeth_mass hw a b hab
  linarith

/-- Avoiding a fat tooth puts the whole `[t, t+ε]` window past the true tooth. -/
theorem window_avoids_of_avoid_fat {w : ℤ} {ε a b t : ℝ}
    (havoid : ∀ p ∈ fatTeeth w ε a b, t ≤ p.1 ∨ p.2 ≤ t) :
    ∀ t' : ℝ, t ≤ t' → t' ≤ t + ε → ∀ p ∈ teeth w a b, t' ≤ p.1 ∨ p.2 ≤ t' := by
  intro t' ht1 ht2 p hp
  have hfat : (p.1 - ε, p.2) ∈ fatTeeth w ε a b := by
    unfold fatTeeth
    rw [List.mem_map]
    exact ⟨p, hp, rfl⟩
  rcases havoid _ hfat with h | h
  · left
    simp only at h
    linarith
  · right
    simp only at h
    linarith

/-- FlatMap clip-sum splitter for fat teeth. -/
theorem fatFlatMap_clip_sum (ws : List ℤ) (ε a b : ℝ) :
    ((ws.flatMap fun w => fatTeeth w ε a b).map fun p => clipLen p a b).sum
      = (ws.map fun (w : ℤ) =>
          ((fatTeeth w ε a b).map fun p => clipLen p a b).sum).sum := by
  induction ws with
  | nil => simp
  | cons w ws ihw =>
      simp only [List.flatMap_cons, List.map_cons, List.map_append,
        List.sum_append, List.sum_cons]
      rw [ihw]

/-! ## The windowed block step -/

/-- **THE WINDOWED BLOCK STEP**: runners whose fat mass fits share a good sub-window
of width `ε`. -/
theorem block_window_step (ws : List ℤ) (hpos : ∀ w ∈ ws, 0 < w) (ε a b : ℝ)
    (hε : 0 ≤ ε) (hab : a ≤ b)
    (hmass : ((ws.map fun (w : ℤ) => (b - a) / 7 + 3 / (7 * (w : ℝ))
        + ε * ((w : ℝ) * (b - a) + 3))).sum < b - a - ε) :
    ∃ x : ℝ, a ≤ x ∧ x + ε ≤ b ∧
      ∀ t : ℝ, x ≤ t → t ≤ x + ε →
        ∀ w ∈ ws, ∀ m : ℤ, (1 : ℝ) / 14 ≤ |(w : ℝ) * t - m| := by
  set bads : List (ℝ × ℝ) := ws.flatMap fun w => fatTeeth w ε a b with hbads
  have hbadsum : (bads.map fun p => clipLen p a b).sum < b - a - ε := by
    rw [hbads, fatFlatMap_clip_sum ws ε a b]
    calc (ws.map fun (w : ℤ) =>
          ((fatTeeth w ε a b).map fun p => clipLen p a b).sum).sum
        ≤ ((ws.map fun (w : ℤ) => (b - a) / 7 + 3 / (7 * (w : ℝ))
            + ε * ((w : ℝ) * (b - a) + 3))).sum := by
          apply List.sum_le_sum
          intro w hw
          exact fatTeeth_mass (hpos w hw) ε a b hε hab
      _ < b - a - ε := hmass
  have hab' : a ≤ b - ε := by
    have h0 : (0 : ℝ) ≤ (bads.map fun p => clipLen p a b).sum := by
      apply List.sum_nonneg
      intro x hx
      rw [List.mem_map] at hx
      obtain ⟨q, _, rfl⟩ := hx
      exact clipLen_nonneg q a b
    linarith
  have hclip_shrink : (bads.map fun p => clipLen p a (b - ε)).sum
      ≤ (bads.map fun p => clipLen p a b).sum := by
    apply List.sum_le_sum
    intro p _
    unfold clipLen
    rcases le_total (min p.2 (b - ε) - max p.1 a) 0 with h | h
    · rw [max_eq_left h]
      exact le_max_left _ _
    · rw [max_eq_right h]
      have : min p.2 (b - ε) ≤ min p.2 b := min_le_min (le_refl _) (by linarith)
      calc min p.2 (b - ε) - max p.1 a ≤ min p.2 b - max p.1 a := by linarith
        _ ≤ max 0 (min p.2 b - max p.1 a) := le_max_right _ _
  obtain ⟨t, ht1, ht2, havoid⟩ := gap_exists bads.length bads a (b - ε) (le_refl _)
    hab' (by
      calc (bads.map fun p => clipLen p a (b - ε)).sum
          ≤ (bads.map fun p => clipLen p a b).sum := hclip_shrink
        _ < b - a - ε := hbadsum
        _ = (b - ε) - a := by ring)
  refine ⟨t, ht1, by linarith, ?_⟩
  intro t' ht'1 ht'2 w hw m
  have hfat : ∀ p ∈ fatTeeth w ε a b, t ≤ p.1 ∨ p.2 ≤ t := by
    intro p hp
    exact havoid p (List.mem_flatMap.mpr ⟨w, hw, hp⟩)
  apply good_of_avoid_teeth (w := w) (a := a) (b := b) (t := t') (hpos w hw)
    (by linarith) (by linarith)
  exact window_avoids_of_avoid_fat hfat t' ht'1 ht'2

/-! ## The unified level chain -/

/-- One chain level: its runners and its output window width. -/
structure GLevel where
  ws  : List ℤ
  out : ℝ

/-- The unified position-free ledger. -/
def GChainOK : List GLevel → ℝ → Prop
  | [], _ => True
  | ℓ :: ℓs, L =>
      (∀ w ∈ ℓ.ws, 0 < w) ∧ 0 ≤ ℓ.out ∧
      ((ℓ.ws.map fun (w : ℤ) => L / 7 + 3 / (7 * (w : ℝ))
        + ℓ.out * ((w : ℝ) * L + 3))).sum < L - ℓ.out ∧
      GChainOK ℓs ℓ.out

/-- **The unified chain engine**: levels of arbitrary mass-fitting blocks nest. -/
theorem glevel_chain : ∀ (ℓs : List GLevel) (a L : ℝ), 0 ≤ L → GChainOK ℓs L →
    ∃ t : ℝ, a ≤ t ∧ t ≤ a + L ∧
      ∀ ℓ ∈ ℓs, ∀ w ∈ ℓ.ws, ∀ m : ℤ, (1 : ℝ) / 14 ≤ |(w : ℝ) * t - m| := by
  intro ℓs
  induction ℓs with
  | nil =>
      intro a L hL _
      exact ⟨a, le_refl _, by linarith, fun ℓ hℓ => absurd hℓ (List.not_mem_nil)⟩
  | cons ℓ ℓs ih =>
      intro a L hL hchain
      obtain ⟨hpos, hout, hmass, hrest⟩ := hchain
      obtain ⟨x, hx1, hx2, hgood⟩ := block_window_step ℓ.ws hpos ℓ.out a (a + L)
        hout (by linarith) (by simpa [add_sub_cancel_left] using hmass)
      obtain ⟨t, ht1, ht2, hℓs⟩ := ih x ℓ.out hout hrest
      refine ⟨t, by linarith, by linarith, ?_⟩
      intro ℓ' hℓ' w hw m
      rcases List.mem_cons.mp hℓ' with rfl | hmem
      · exact hgood t ht1 ht2 w hw m
      · exact hℓs ℓ' hmem w hw m

/-- **THE CITE–UNIFIED-CHAIN THEOREM**: cite `k ≤ 12` bounded runners; every other
runner's absolute value sits in some level of a `GChainOK` chain for the transported
window `2δ`. -/
theorem cite_glevel_lonely (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (k : ℕ) (hk : k ≤ 12) (B : ℤ) (hB : 0 < B)
    (hcited : ∀ i : Fin 13, (i : ℕ) < k → |v i| ≤ B)
    (ℓs : List GLevel)
    (hsplit : ∀ i : Fin 13, (i : ℕ) < k ∨ ∃ ℓ ∈ ℓs, |v i| ∈ ℓ.ws)
    (hchain : GChainOK ℓs
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
  obtain ⟨τ, hτ1, hτ2, hτgood⟩ := glevel_chain ℓs (t₀ - δ) (2 * δ) (by linarith) hchain
  refine ⟨τ, ?_⟩
  intro i m
  rcases hsplit i with hlt | ⟨ℓ, hℓ, hmem⟩
  · have hidx : Fin.castLE hk13 ⟨(i : ℕ), hlt⟩ = i := by
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
  · have hgood := hτgood ℓ hℓ |v i| hmem
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
