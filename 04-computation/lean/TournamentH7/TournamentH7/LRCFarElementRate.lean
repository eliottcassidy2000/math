/-
  TournamentH7.LRCFarElementRate — THE FAR-ELEMENT RATE LEMMA (klein HYP-4001, Lean face).
  kind-pasteur-2026-07-02-S11 (HYP-3968).  Single-writer satellite.

  THE LEMMA (klein S94, wrap-counting; the peel leg of the census+peel surface — the ONE
  statement the fleet ledgers listed as unformalized): adding a far element `w` to a shape
  changes each sector-miss measure by the damped factor up to `O(1/w)`:

      | meas{x ∈ I : frac(wx) ∈ sector} − h·|I| |  ≤  2h/w        (per interval I),

  "each sector contributes at most two partial wraps at I's ends".

  THE FORMALIZATION (counting-free).  In tooth-aligned coordinates (teeth `[k, k+ρ)`,
  `ρ` = sector width, interval `[A, B)`), group the line into unit cells `[m, m+1)`.
  Per cell, the tooth term `t_m` and the coverage term `c_m` satisfy the EXACT identity
  `t_m = ρ·c_m` — except possibly at the TWO indices `m = ⌊A⌋` and `m = ⌊B⌋`, where
  `|t_m − ρ·c_m| ≤ ρ`.  No integer counting anywhere: the exceptional indices are NAMED
  by floors, and the coverage total is the `clip_chain_sum` telescope (LRCSevenTranslate).

  Main results:
   * `rate_core`  : `|Σ_{k<w} [min B (k+ρ) − max A k]⁺ − ρ(B−A)| ≤ 2ρ`  (−1 ≤ A ≤ B ≤ w)
   * `length_inter_comb_near` : the comb form — for one interval `[a,b) ⊆ [0,1)` and the
     `w`-comb at radius `r` (non-wrapping band `r ≤ φ`, `φ + r ≤ 1`):
     `|length (inter [(a,b)] (comb w r φ)) − 2r(b−a)| ≤ 4r/w`.
-/
import TournamentH7.LRCSevenTranslate

namespace LonelyRunner
namespace RatIntervals

/-! ## List helpers -/

private theorem sum_indicator_le {i : ℤ} {c : ℚ} (hc : 0 ≤ c) :
    ∀ {l : List ℕ}, l.Nodup →
      ((l.map fun k : ℕ => if (k : ℤ) = i then c else 0).sum ≤ c) := by
  intro l
  induction l with
  | nil => intro _; simpa using hc
  | cons a l ih =>
      intro hnd
      simp only [List.map_cons, List.sum_cons]
      by_cases ha : ((a : ℤ) = i)
      · rw [if_pos ha]
        have hzero : ((l.map fun k : ℕ => if (k : ℤ) = i then c else 0).sum = 0) := by
          apply List.sum_eq_zero
          intro x hx
          simp only [List.mem_map] at hx
          obtain ⟨b, hb, rfl⟩ := hx
          have hbne : ((b : ℤ) ≠ i) := by
            intro hbi
            have hab : a = b := by exact_mod_cast ha.trans hbi.symm
            exact (List.nodup_cons.mp hnd).1 (hab ▸ hb)
          rw [if_neg hbne]
        rw [hzero]
        linarith
      · rw [if_neg ha]
        have := ih (List.nodup_cons.mp hnd).2
        linarith

private theorem sum_indicator_eq_zero_of_neg {i : ℤ} (hi : i < 0) {c : ℚ} :
    ∀ l : List ℕ, ((l.map fun k : ℕ => if (k : ℤ) = i then c else 0).sum = 0) := by
  intro l
  apply List.sum_eq_zero
  intro x hx
  simp only [List.mem_map] at hx
  obtain ⟨b, hb, rfl⟩ := hx
  have hne : ((b : ℤ) ≠ i) := by
    have : (0 : ℤ) ≤ (b : ℤ) := Int.natCast_nonneg b
    omega
  rw [if_neg hne]

private theorem abs_sum_le_sum_abs' (f : ℕ → ℚ) :
    ∀ l : List ℕ, |(l.map f).sum| ≤ (l.map fun k => |f k|).sum := by
  intro l
  induction l with
  | nil => simp
  | cons a l ih =>
      simp only [List.map_cons, List.sum_cons]
      calc |f a + (l.map f).sum| ≤ |f a| + |(l.map f).sum| := abs_add_le _ _
        _ ≤ |f a| + (l.map fun k => |f k|).sum := by linarith

private theorem sum_le_sum' (f g : ℕ → ℚ) :
    ∀ l : List ℕ, (∀ k ∈ l, f k ≤ g k) → (l.map f).sum ≤ (l.map g).sum := by
  intro l
  induction l with
  | nil => intro _; simp
  | cons a l ih =>
      intro h
      simp only [List.map_cons, List.sum_cons]
      have h1 := h a (List.mem_cons_self ..)
      have h2 := ih fun k hk => h k (List.mem_cons_of_mem _ hk)
      linarith

private theorem sum_map_add' (f g : ℕ → ℚ) :
    ∀ l : List ℕ, (l.map fun k => f k + g k).sum = (l.map f).sum + (l.map g).sum := by
  intro l
  induction l with
  | nil => simp
  | cons a l ih =>
      simp only [List.map_cons, List.sum_cons, ih]
      ring

private theorem sum_map_sub_mul (f g : ℕ → ℚ) (q : ℚ) :
    ∀ l : List ℕ, (l.map fun k => f k - q * g k).sum
      = (l.map f).sum - q * (l.map g).sum := by
  intro l
  induction l with
  | nil => simp
  | cons a l ih =>
      simp only [List.map_cons, List.sum_cons, ih]
      ring

private theorem sum_map_mul_left' (q : ℚ) (f : ℕ → ℚ) :
    ∀ l : List ℕ, (l.map fun k => q * f k).sum = q * (l.map f).sum := by
  intro l
  induction l with
  | nil => simp
  | cons a l ih =>
      simp only [List.map_cons, List.sum_cons, ih]
      ring

/-! ## The core: tooth mass vs damped length, two boundary groups -/

/-- **The rate-lemma core (counting-free).**  Unit-spaced teeth `[k, k+ρ)`, `k < w`,
clipped against `[A, B)` with `−1 ≤ A ≤ B ≤ w`, `0 ≤ B`: the tooth mass equals
`ρ·(B − A)` up to `2ρ` — the two boundary cells `⌊A⌋` and `⌊B⌋`. -/
theorem rate_core {ρ : ℚ} (hρ0 : 0 ≤ ρ) (hρ1 : ρ ≤ 1) (w : ℕ)
    {A B : ℚ} (hA1 : -1 ≤ A) (hAB : A ≤ B) (hB0 : 0 ≤ B) (hBw : B ≤ w) :
    |((List.range w).map fun k : ℕ =>
        max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ))).sum - ρ * (B - A)| ≤ 2 * ρ := by
  -- the coverage total: the unit chain telescope
  have hcov : ((List.range w).map fun k : ℕ =>
      max 0 (min B ((k : ℚ) + 1) - max A (k : ℚ))).sum = B - max A 0 := by
    have h := clip_chain_sum (w := (1 : ℚ)) (by norm_num) w 0 A B hAB
      (by push_cast; linarith)
    have hml : (List.range w).map (fun K : ℕ =>
        max 0 (min B ((0 : ℚ) + ((K : ℚ) + 1) * 1) - max A ((0 : ℚ) + (K : ℚ) * 1)))
        = (List.range w).map (fun K : ℕ =>
        max 0 (min B ((K : ℚ) + 1) - max A (K : ℚ))) := by
      apply List.map_congr_left
      intro K _
      norm_num
    rw [hml] at h
    rw [h]
    have hmax : max A 0 ≤ B := max_le hAB hB0
    exact max_eq_right (by linarith)
  -- pointwise: the difference vanishes except at the two boundary cells
  have hpt : ∀ k ∈ List.range w,
      |max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ))
        - ρ * max 0 (min B ((k : ℚ) + 1) - max A (k : ℚ))|
      ≤ (if (k : ℤ) = ⌊A⌋ then ρ else 0) + (if (k : ℤ) = ⌊B⌋ then ρ else 0) := by
    intro k _
    have hind1 : (0 : ℚ) ≤ (if (k : ℤ) = ⌊A⌋ then ρ else 0) := by
      split <;> simp [hρ0]
    have hind2 : (0 : ℚ) ≤ (if (k : ℤ) = ⌊B⌋ then ρ else 0) := by
      split <;> simp [hρ0]
    by_cases hcase1 : (k : ℚ) + 1 ≤ A
    · -- cell entirely below the interval: both terms vanish
      have ht : max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ)) = 0 := by
        apply max_eq_left
        have h1 : min B ((k : ℚ) + ρ) ≤ (k : ℚ) + ρ := min_le_right _ _
        have h2 : A ≤ max A (k : ℚ) := le_max_left _ _
        linarith
      have hc : max 0 (min B ((k : ℚ) + 1) - max A (k : ℚ)) = 0 := by
        apply max_eq_left
        have h1 : min B ((k : ℚ) + 1) ≤ (k : ℚ) + 1 := min_le_right _ _
        have h2 : A ≤ max A (k : ℚ) := le_max_left _ _
        linarith
      rw [ht, hc]
      simp only [mul_zero, sub_zero, abs_zero]
      linarith
    · by_cases hcase2 : B ≤ (k : ℚ)
      · -- cell entirely above the interval: both terms vanish
        have ht : max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ)) = 0 := by
          apply max_eq_left
          have h1 : min B ((k : ℚ) + ρ) ≤ B := min_le_left _ _
          have h2 : (k : ℚ) ≤ max A (k : ℚ) := le_max_right _ _
          linarith
        have hc : max 0 (min B ((k : ℚ) + 1) - max A (k : ℚ)) = 0 := by
          apply max_eq_left
          have h1 : min B ((k : ℚ) + 1) ≤ B := min_le_left _ _
          have h2 : (k : ℚ) ≤ max A (k : ℚ) := le_max_right _ _
          linarith
        rw [ht, hc]
        simp only [mul_zero, sub_zero, abs_zero]
        linarith
      · by_cases hcase3 : A ≤ (k : ℚ) ∧ (k : ℚ) + 1 ≤ B
        · -- interior cell: exact identity t = ρ·c
          obtain ⟨hAk, hkB⟩ := hcase3
          have hmax : max A (k : ℚ) = (k : ℚ) := max_eq_right hAk
          have hmin1 : min B ((k : ℚ) + ρ) = (k : ℚ) + ρ :=
            min_eq_right (by linarith)
          have hmin2 : min B ((k : ℚ) + 1) = (k : ℚ) + 1 :=
            min_eq_right hkB
          rw [hmax, hmin1, hmin2]
          have h1 : max 0 ((k : ℚ) + ρ - (k : ℚ)) = ρ := by
            rw [show (k : ℚ) + ρ - (k : ℚ) = ρ by ring]
            exact max_eq_right hρ0
          have h2 : max 0 ((k : ℚ) + 1 - (k : ℚ)) = 1 := by
            rw [show (k : ℚ) + 1 - (k : ℚ) = (1 : ℚ) by ring]
            exact max_eq_right (by norm_num)
          rw [h1, h2, mul_one, sub_self, abs_zero]
          linarith
        · -- boundary cell: k = ⌊A⌋ or k = ⌊B⌋; crude bound |t − ρc| ≤ ρ
          have hA' : A < (k : ℚ) + 1 := lt_of_not_ge hcase1
          have hB' : (k : ℚ) < B := lt_of_not_ge hcase2
          have hbnd : ((k : ℤ) = ⌊A⌋) ∨ ((k : ℤ) = ⌊B⌋) := by
            rcases le_or_gt A (k : ℚ) with hAk | hAk
            · -- then ¬(k+1 ≤ B): k ≤ B < k+1, so k = ⌊B⌋
              right
              have hBk1 : B < (k : ℚ) + 1 := by
                by_contra hcon
                exact hcase3 ⟨hAk, le_of_not_gt hcon⟩
              have hfl : ⌊B⌋ = (k : ℤ) := by
                rw [Int.floor_eq_iff]
                constructor
                · push_cast; linarith
                · push_cast; linarith
              omega
            · -- k < A < k+1, so k = ⌊A⌋
              left
              have hfl : ⌊A⌋ = (k : ℤ) := by
                rw [Int.floor_eq_iff]
                constructor
                · push_cast; linarith
                · push_cast; linarith
              omega
          -- crude bounds: 0 ≤ t ≤ ρ, 0 ≤ c ≤ 1
          have ht0 : (0 : ℚ) ≤ max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ)) :=
            le_max_left _ _
          have htρ : max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ)) ≤ ρ := by
            apply max_le hρ0
            have h1 : min B ((k : ℚ) + ρ) ≤ (k : ℚ) + ρ := min_le_right _ _
            have h2 : (k : ℚ) ≤ max A (k : ℚ) := le_max_right _ _
            linarith
          have hc0 : (0 : ℚ) ≤ max 0 (min B ((k : ℚ) + 1) - max A (k : ℚ)) :=
            le_max_left _ _
          have hc1 : max 0 (min B ((k : ℚ) + 1) - max A (k : ℚ)) ≤ 1 := by
            apply max_le (by norm_num)
            have h1 : min B ((k : ℚ) + 1) ≤ (k : ℚ) + 1 := min_le_right _ _
            have h2 : (k : ℚ) ≤ max A (k : ℚ) := le_max_right _ _
            linarith
          have habs : |max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ))
              - ρ * max 0 (min B ((k : ℚ) + 1) - max A (k : ℚ))| ≤ ρ := by
            rw [abs_le]
            constructor
            · nlinarith
            · nlinarith
          rcases hbnd with h | h
          · rw [if_pos h]
            linarith
          · rw [if_pos h]
            linarith
  -- assemble
  have hd : ((List.range w).map fun k : ℕ =>
        (max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ))
          - ρ * max 0 (min B ((k : ℚ) + 1) - max A (k : ℚ)))).sum
      = ((List.range w).map fun k : ℕ =>
          max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ))).sum - ρ * (B - max A 0) := by
    rw [sum_map_sub_mul, hcov]
  have hDbound : |((List.range w).map fun k : ℕ =>
        max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ))).sum - ρ * (B - max A 0)|
      ≤ (((List.range w).map fun k : ℕ =>
          if (k : ℤ) = ⌊A⌋ then ρ else 0).sum)
        + (((List.range w).map fun k : ℕ =>
          if (k : ℤ) = ⌊B⌋ then ρ else 0).sum) := by
    rw [← hd]
    calc |((List.range w).map fun k : ℕ =>
            (max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ))
              - ρ * max 0 (min B ((k : ℚ) + 1) - max A (k : ℚ)))).sum|
        ≤ ((List.range w).map fun k : ℕ =>
            |max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ))
              - ρ * max 0 (min B ((k : ℚ) + 1) - max A (k : ℚ))|).sum :=
          abs_sum_le_sum_abs' _ _
      _ ≤ ((List.range w).map fun k : ℕ =>
            (if (k : ℤ) = ⌊A⌋ then ρ else 0) + (if (k : ℤ) = ⌊B⌋ then ρ else 0)).sum :=
          sum_le_sum' _ _ _ hpt
      _ = _ := sum_map_add' _ _ _
  have hindB : ((List.range w).map fun k : ℕ =>
      if (k : ℤ) = ⌊B⌋ then ρ else 0).sum ≤ ρ :=
    sum_indicator_le hρ0 List.nodup_range
  rcases le_or_gt 0 A with hA0 | hA0
  · -- A ≥ 0: max A 0 = A, both indicator sums ≤ ρ
    have hmaxA : max A 0 = A := max_eq_left hA0
    have hindA : ((List.range w).map fun k : ℕ =>
        if (k : ℤ) = ⌊A⌋ then ρ else 0).sum ≤ ρ :=
      sum_indicator_le hρ0 List.nodup_range
    rw [hmaxA] at hDbound
    have habs := abs_le.mp hDbound
    rw [abs_le]
    constructor <;> linarith
  · -- −1 ≤ A < 0: the ⌊A⌋ indicator sum is ZERO; pay ρ·|A| ≤ ρ instead
    have hfloorA : ⌊A⌋ < 0 := Int.floor_lt.mpr (by push_cast; linarith)
    have hindA : ((List.range w).map fun k : ℕ =>
        if (k : ℤ) = ⌊A⌋ then ρ else 0).sum = 0 :=
      sum_indicator_eq_zero_of_neg hfloorA _
    have hmaxA : max A 0 = 0 := max_eq_right (le_of_lt hA0)
    rw [hmaxA, hindA] at hDbound
    have hb1 := abs_le.mp hDbound
    have hcorr : |ρ * (B - 0) - ρ * (B - A)| ≤ ρ := by
      rw [show ρ * (B - 0) - ρ * (B - A) = ρ * A by ring, abs_mul, abs_of_nonneg hρ0]
      have hAabs : |A| ≤ 1 := by
        rw [abs_le]
        constructor <;> linarith
      nlinarith
    have hc1 := abs_le.mp hcorr
    rw [abs_le]
    constructor <;> linarith

/-! ## The comb form -/

/-- **The far-element rate lemma, comb form** (klein HYP-4001's per-interval statement):
for one interval `[a, b) ⊆ [0, 1)` and the `w`-comb at radius `r` with a non-wrapping
band (`r ≤ φ`, `φ + r ≤ 1`), the intersection length is `2r(b−a)` up to `4r/w` —
two partial wraps at the interval's ends. -/
theorem length_inter_comb_near {w : ℕ} (hw : 0 < w) {r φ a b : ℚ}
    (hr : 0 ≤ r) (hrφ : r ≤ φ) (hφ1 : φ + r ≤ 1)
    (ha : 0 ≤ a) (hab : a ≤ b) (hb : b ≤ 1) :
    |length (inter [(a, b)] (comb w r φ)) - 2 * r * (b - a)| ≤ 4 * r / w := by
  have hwQ : (0 : ℚ) < (w : ℚ) := by exact_mod_cast hw
  -- unfold the single-interval intersection into the per-tooth clip sum
  have hunf : length (inter [(a, b)] (comb w r φ))
      = ((List.range w).map fun k : ℕ =>
          max 0 (min b (((k : ℚ) + φ + r) / w) - max a (((k : ℚ) + φ - r) / w))).sum := by
    unfold inter comb length clip
    simp only [List.flatMap_cons, List.flatMap_nil, List.append_nil, List.map_map]
    congr 1
  set A : ℚ := (w : ℚ) * a - φ + r with hA
  set B : ℚ := (w : ℚ) * b - φ + r with hB
  set ρ : ℚ := 2 * r with hρ
  have hρ1 : ρ ≤ 1 := by rw [hρ]; linarith
  have hρ0 : 0 ≤ ρ := by rw [hρ]; linarith
  rcases le_or_gt 0 B with hB0 | hB0
  · -- main path: change variables and apply the core
    have hterm : ∀ k : ℕ,
        max 0 (min b (((k : ℚ) + φ + r) / w) - max a (((k : ℚ) + φ - r) / w))
        = (1 / w) * max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ)) := by
      intro k
      have hb_iff : b ≤ ((k : ℚ) + φ + r) / w ↔ B ≤ (k : ℚ) + ρ := by
        rw [le_div_iff₀ hwQ, hB, hρ]
        constructor <;> intro h <;> nlinarith [h]
      have hb_iff2 : ((k : ℚ) + φ + r) / w ≤ b ↔ (k : ℚ) + ρ ≤ B := by
        rw [div_le_iff₀ hwQ, hB, hρ]
        constructor <;> intro h <;> nlinarith [h]
      have ha_iff : a ≤ ((k : ℚ) + φ - r) / w ↔ A ≤ (k : ℚ) := by
        rw [le_div_iff₀ hwQ, hA]
        constructor <;> intro h <;> nlinarith [h]
      have ha_iff2 : ((k : ℚ) + φ - r) / w ≤ a ↔ (k : ℚ) ≤ A := by
        rw [div_le_iff₀ hwQ, hA]
        constructor <;> intro h <;> nlinarith [h]
      have hminc : min b (((k : ℚ) + φ + r) / w)
          = (1 / w) * min B ((k : ℚ) + ρ) + (φ - r) / w := by
        rcases le_total b (((k : ℚ) + φ + r) / w) with h | h
        · rw [min_eq_left h, min_eq_left (hb_iff.mp h), hB]
          field_simp
          try ring
        · rw [min_eq_right h, min_eq_right (hb_iff2.mp h), hρ]
          field_simp
          try ring
      have hmaxc : max a (((k : ℚ) + φ - r) / w)
          = (1 / w) * max A (k : ℚ) + (φ - r) / w := by
        rcases le_total a (((k : ℚ) + φ - r) / w) with h | h
        · rw [max_eq_right h, max_eq_right (ha_iff.mp h)]
          field_simp
          try ring
        · rw [max_eq_left h, max_eq_left (ha_iff2.mp h), hA]
          field_simp
          try ring
      rw [hminc, hmaxc,
        show (1 / (w : ℚ)) * min B ((k : ℚ) + ρ) + (φ - r) / w
          - ((1 / (w : ℚ)) * max A (k : ℚ) + (φ - r) / w)
          = (1 / (w : ℚ)) * (min B ((k : ℚ) + ρ) - max A (k : ℚ)) by ring]
      rcases le_total (min B ((k : ℚ) + ρ) - max A (k : ℚ)) 0 with h | h
      · rw [max_eq_left (mul_nonpos_of_nonneg_of_nonpos (by positivity) h),
          max_eq_left h, mul_zero]
      · rw [max_eq_right (mul_nonneg (by positivity) h), max_eq_right h]
    have hsum : length (inter [(a, b)] (comb w r φ))
        = (1 / w) * ((List.range w).map fun k : ℕ =>
            max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ))).sum := by
      rw [hunf]
      have hml : (List.range w).map (fun k : ℕ =>
          max 0 (min b (((k : ℚ) + φ + r) / w) - max a (((k : ℚ) + φ - r) / w)))
          = (List.range w).map (fun k : ℕ =>
          (1 / (w : ℚ)) * max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ))) := by
        apply List.map_congr_left
        intro k _
        exact hterm k
      rw [hml, sum_map_mul_left']
    -- the core's hypotheses
    have hwa : (0 : ℚ) ≤ (w : ℚ) * a := mul_nonneg hwQ.le ha
    have hwab : (w : ℚ) * a ≤ (w : ℚ) * b := mul_le_mul_of_nonneg_left hab hwQ.le
    have hwb1 : (w : ℚ) * b ≤ (w : ℚ) := by
      calc (w : ℚ) * b ≤ (w : ℚ) * 1 := mul_le_mul_of_nonneg_left hb hwQ.le
        _ = (w : ℚ) := mul_one _
    have hcore := rate_core (ρ := ρ) hρ0 hρ1 w (A := A) (B := B)
      (by rw [hA]; linarith) (by rw [hA, hB]; linarith) hB0 (by rw [hB]; linarith)
    have hfinal : length (inter [(a, b)] (comb w r φ)) - 2 * r * (b - a)
        = (1 / w) * (((List.range w).map fun k : ℕ =>
            max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ))).sum - ρ * (B - A)) := by
      rw [hsum]
      have hBA : B - A = (w : ℚ) * (b - a) := by rw [hA, hB]; ring
      rw [hBA, hρ]
      field_simp
      try ring
    rw [hfinal, abs_mul, abs_of_pos (by positivity : (0 : ℚ) < 1 / w)]
    calc (1 / (w : ℚ)) * |((List.range w).map fun k : ℕ =>
            max 0 (min B ((k : ℚ) + ρ) - max A (k : ℚ))).sum - ρ * (B - A)|
        ≤ (1 / (w : ℚ)) * (2 * ρ) := by
          apply mul_le_mul_of_nonneg_left hcore (by positivity)
      _ = 4 * r / w := by rw [hρ]; field_simp; ring
  · -- B < 0: the comb misses [a,b) entirely, and 2r(b−a) itself is < 4r/w
    have hwb : (w : ℚ) * b < φ - r := by rw [hB] at hB0; linarith
    have hzero : length (inter [(a, b)] (comb w r φ)) = 0 := by
      rw [hunf]
      apply List.sum_eq_zero
      intro x hx
      simp only [List.mem_map] at hx
      obtain ⟨k, _, rfl⟩ := hx
      apply max_eq_left
      have hk0 : (0 : ℚ) ≤ (k : ℚ) := by positivity
      have h2 : b < ((k : ℚ) + φ - r) / w := by
        rw [lt_div_iff₀ hwQ]
        nlinarith
      have h3 : min b (((k : ℚ) + φ + r) / w) ≤ b := min_le_left _ _
      have h4 : ((k : ℚ) + φ - r) / w ≤ max a (((k : ℚ) + φ - r) / w) := le_max_right _ _
      linarith
    rw [hzero]
    have hba : b - a ≤ b := by linarith
    have hblt : (w : ℚ) * (b - a) ≤ 1 := by nlinarith
    have h2rb : 2 * r * (b - a) ≤ 2 * r / w := by
      rw [div_eq_mul_inv]
      have hinv : (0 : ℚ) < ((w : ℚ))⁻¹ := by positivity
      have : 2 * r * (b - a) * (w : ℚ) ≤ 2 * r * 1 := by nlinarith
      calc 2 * r * (b - a) = 2 * r * (b - a) * (w : ℚ) * ((w : ℚ))⁻¹ := by
            field_simp
            try ring
        _ ≤ 2 * r * 1 * ((w : ℚ))⁻¹ := by
            apply mul_le_mul_of_nonneg_right this hinv.le
        _ = 2 * r * ((w : ℚ))⁻¹ := by ring
    have hnn : (0 : ℚ) ≤ 2 * r * (b - a) := by positivity
    rw [show (0 : ℚ) - 2 * r * (b - a) = -(2 * r * (b - a)) by ring, abs_neg,
      abs_of_nonneg hnn]
    have h24 : 2 * r / (w : ℚ) ≤ 4 * r / w := by
      apply div_le_div_of_nonneg_right ?_ hwQ.le
      linarith
    linarith

end RatIntervals
end LonelyRunner
