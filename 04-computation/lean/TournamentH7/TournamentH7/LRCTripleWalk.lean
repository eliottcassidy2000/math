/-
  TournamentH7.LRCTripleWalk — THE TRIPLE WALK: the l = 3 walk by recursion into
  the pair walk (kind-pasteur-2026-07-05-S9, HYP-4118).

  THE RECURSION: in a 3-comb cover of an interval, a comb REVISIT exposes a full
  consecutive inter-tooth gap of the revisited comb — length `(21/25)/p` — whose
  interior that comb cannot touch, so the OTHER TWO combs must cover it, and the
  S8 pair walk (`walk_core`) bounds any 2-comb-covered interval by
  `(4/25)(2/q₁ + 1/q₂)`.  Under pairwise balance `4·big ≤ 7·small` this is at
  most `3·(4/25)·(7/(4p)) = (21/25)/p` — EXACTLY the gap length, and the walk
  bound is STRICT: contradiction (`revisit_kill`).  Hence a balanced triple
  cover never revisits: at most THREE teeth, one per comb (`triple_walk`):

      `b − a < (4/25)(1/u + 1/v + 1/w)`.

  THE RUNG (`gap_triple_rung`): in a gap violator (no `2/25`-margin point),
  every PAIRWISE-BALANCED 3-subset (all ratios ≤ 7/4) has `min ≤ 12·B`, `B` a
  bound on the other nine — the 9-complement is cited at `1/10`, the excess
  `1/10 − 2/25 = 1/50` opens a window of length `1/(25B)`, and three balanced
  teeth cannot span it unless the smallest speed is `≤ 12B`.

  The walk ladder: l = 1 → 24B (S6), l = 2 → 22B (S8), l = 3 balanced → 12B.
  The 7/4 balance is exactly where the pair-walk recursion saturates;
  unbalanced triples fall to the pair rung (S8) or the gap descent (opus-S81).

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCPairWalk

namespace LonelyRunner
namespace TripleWalk

open TournamentH7.LRCWitness GridAttainment CertCompleteness PairWalk TeethR

/-- Points in the closed consecutive inter-tooth gap of comb `p` are at scaled
distance `≥ 2/25` from every integer. -/
lemma in_gap_far (p : ℤ) (hp : 0 < p) (m : ℤ) (t : ℝ)
    (h1 : ((m : ℝ) + 2/25) / p ≤ t) (h2 : t ≤ ((m : ℝ) + 1 - 2/25) / p) :
    ∀ n : ℤ, (2 : ℝ)/25 ≤ |(p : ℝ) * t - n| := by
  intro n
  have hpR : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  have hlo : (m : ℝ) + 2/25 ≤ (p : ℝ) * t := by
    rw [div_le_iff₀ hpR] at h1
    linarith
  have hhi : (p : ℝ) * t ≤ (m : ℝ) + 1 - 2/25 := by
    rw [le_div_iff₀ hpR] at h2
    linarith
  rcases le_or_gt (n : ℝ) (m : ℝ) with hn | hn
  · rw [abs_of_nonneg (by linarith)]
    linarith
  · have hn1 : (m : ℝ) + 1 ≤ (n : ℝ) := by
      have : m < n := by exact_mod_cast hn
      exact_mod_cast Int.add_one_le_iff.mpr this
    rw [abs_of_nonpos (by linarith), neg_sub]
    linarith

/-- **THE REVISIT KILL**: a full consecutive gap of comb `p` covered by two combs
`q₁, q₂` (left end on `q₁`) with balance `4p ≤ 7q₁`, `4p ≤ 7q₂` is impossible —
the pair walk bounds the cover by exactly the gap length, strictly. -/
lemma revisit_kill (p q₁ q₂ : ℤ) (hp : 0 < p) (hq₁ : 0 < q₁) (hq₂ : 0 < q₂)
    (hbal₁ : 4 * p ≤ 7 * q₁) (hbal₂ : 4 * p ≤ 7 * q₂) (m : ℤ)
    (hcovgap : ∀ t, ((m : ℝ) + 2/25) / p ≤ t → t ≤ ((m : ℝ) + 1 - 2/25) / p →
      (∃ n : ℤ, |(q₁ : ℝ) * t - n| < 2/25) ∨ (∃ n : ℤ, |(q₂ : ℝ) * t - n| < 2/25))
    (hstart : ∃ n : ℤ, |(q₁ : ℝ) * (((m : ℝ) + 2/25) / p) - n| < 2/25) : False := by
  have hpR : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  have hq₁R : (0 : ℝ) < (q₁ : ℝ) := by exact_mod_cast hq₁
  have hq₂R : (0 : ℝ) < (q₂ : ℝ) := by exact_mod_cast hq₂
  have hwalk := walk_core q₁ q₂ hq₁ hq₂ _ _ hcovgap hstart
  have hlen : ((m : ℝ) + 1 - 2/25) / p - ((m : ℝ) + 2/25) / p = (21/25) / p := by
    field_simp
    ring
  rw [hlen] at hwalk
  have hb₁ : (4 : ℝ)/25 / q₁ ≤ (7/25) / p := by
    rw [div_le_div_iff₀ hq₁R hpR]
    have h4 : (4 : ℝ) * p ≤ 7 * q₁ := by exact_mod_cast hbal₁
    linarith
  have hb₂ : (4 : ℝ)/25 / q₂ ≤ (7/25) / p := by
    rw [div_le_div_iff₀ hq₂R hpR]
    have h4 : (4 : ℝ) * p ≤ 7 * q₂ := by exact_mod_cast hbal₂
    linarith
  have hcollect : (7/25 : ℝ)/p + (7/25)/p + (7/25)/p = (21/25)/p := by ring
  linarith

/-- **The walk tail** (the continuation after the first `u`-tooth, second comb `x`,
third comb `y`): under the balance needed for `u`- and `x`-revisit kills, the
interval spans less than one tooth of each comb. -/
lemma walk_tail (u x y : ℤ) (hu : 0 < u) (hx : 0 < x) (hy : 0 < y)
    (hux : 4 * u ≤ 7 * x) (huy : 4 * u ≤ 7 * y)
    (hxy : 4 * x ≤ 7 * y) (hxu : 4 * x ≤ 7 * u)
    (a b : ℝ)
    (hcov : ∀ t, a ≤ t → t ≤ b →
      (∃ n : ℤ, |(u : ℝ) * t - n| < 2/25) ∨ ((∃ n : ℤ, |(x : ℝ) * t - n| < 2/25) ∨
      (∃ n : ℤ, |(y : ℝ) * t - n| < 2/25)))
    (m₀ : ℤ) (ha_e₀ : a < ((m₀ : ℝ) + 2/25) / u)
    (he₀a : ((m₀ : ℝ) + 2/25) / u - a < (4:ℝ)/25 / u)
    (he₀b : ((m₀ : ℝ) + 2/25) / u ≤ b)
    (hm₁ : ∃ n : ℤ, |(x : ℝ) * (((m₀ : ℝ) + 2/25) / u) - n| < 2/25) :
    b - a < (4:ℝ)/25 / u + (4:ℝ)/25 / x + (4:ℝ)/25 / y := by
  have huR : (0 : ℝ) < (u : ℝ) := by exact_mod_cast hu
  have hxR : (0 : ℝ) < (x : ℝ) := by exact_mod_cast hx
  have hyR : (0 : ℝ) < (y : ℝ) := by exact_mod_cast hy
  have hposu : (0:ℝ) < (4:ℝ)/25 / u := by positivity
  have hposx : (0:ℝ) < (4:ℝ)/25 / x := by positivity
  have hposy : (0:ℝ) < (4:ℝ)/25 / y := by positivity
  set e₀ : ℝ := ((m₀ : ℝ) + 2/25) / u with he₀def
  have hue₀ : (u : ℝ) * e₀ = (m₀ : ℝ) + 2/25 := by rw [he₀def]; field_simp
  -- the shared u-revisit killer
  have hkillU : ∀ tt, e₀ < tt → tt ≤ b →
      (∃ n : ℤ, |(u : ℝ) * tt - n| < 2/25) → False := by
    rintro tt htt httb ⟨n, hn⟩
    have habs := abs_lt.mp hn
    have hult : (m₀ : ℝ) + 2/25 < (u : ℝ) * tt := by
      rw [← hue₀]
      exact mul_lt_mul_of_pos_left htt huR
    have hmrev : m₀ + 1 ≤ n := by
      have : (m₀ : ℝ) < (n : ℝ) := by linarith
      have hmn : m₀ < n := by exact_mod_cast this
      exact Int.add_one_le_iff.mpr hmn
    have hmb : ((n : ℝ) - 2/25) / u ≤ b := by
      rw [div_le_iff₀ huR]
      have h1 : ((n : ℝ) - 2/25) < (u : ℝ) * tt := by linarith
      have h2 : (u : ℝ) * tt ≤ b * u := by
        rw [mul_comm]
        exact mul_le_mul_of_nonneg_right httb huR.le
      linarith
    apply revisit_kill u x y hu hx hy hux huy m₀ ?_ (by rw [← he₀def]; exact hm₁)
    intro t ht1 ht2
    have hta : a ≤ t := by
      have : e₀ ≤ t := by rw [he₀def]; exact ht1
      linarith
    have htb : t ≤ b := by
      have h1 : ((m₀ : ℝ) + 1 - 2/25) / u ≤ ((n : ℝ) - 2/25) / u := by
        apply div_le_div_of_nonneg_right ?_ huR.le
        have : (m₀ : ℝ) + 1 ≤ (n : ℝ) := by exact_mod_cast hmrev
        linarith
      linarith
    have hnotu : ¬ ∃ nn : ℤ, |(u : ℝ) * t - nn| < 2/25 := by
      rintro ⟨nn, hnn⟩
      exact absurd hnn (not_lt.mpr (in_gap_far u hu m₀ t ht1 ht2 nn))
    exact (hcov t hta htb).resolve_left hnotu
  -- step 2: comb x from e₀
  obtain ⟨n₁, hn₁⟩ := hm₁
  have hm₁' : |(x : ℝ) * e₀ - n₁| < 2/25 := by rw [he₀def]; exact hn₁
  obtain ⟨he₁gt, he₁lt⟩ := step_data x hx e₀ n₁ hm₁'
  set e₁ : ℝ := ((n₁ : ℝ) + 2/25) / x with he₁def
  have hxe₁ : (x : ℝ) * e₁ = (n₁ : ℝ) + 2/25 := by rw [he₁def]; field_simp
  by_cases hb₁ : b < e₁
  · linarith
  push_neg at hb₁
  have he₁a : a ≤ e₁ := by linarith
  have hnotx₁ : ¬ ∃ n : ℤ, |(x : ℝ) * e₁ - n| < 2/25 := by
    rintro ⟨n, hn⟩
    exact absurd hn (not_lt.mpr (boundary_uncovered x hx n₁ n))
  have hcov₁ := hcov e₁ he₁a hb₁
  rcases hcov₁ with hcu | hcx | hcy
  · exact (hkillU e₁ he₁gt hb₁ hcu).elim
  · exact absurd hcx hnotx₁
  · -- step 3: comb y from e₁
    obtain ⟨n₂, hn₂⟩ := hcy
    obtain ⟨he₂gt, he₂lt⟩ := step_data y hy e₁ n₂ hn₂
    set e₂ : ℝ := ((n₂ : ℝ) + 2/25) / y with he₂def
    by_cases hb₂ : b < e₂
    · linarith
    push_neg at hb₂
    have he₂a : a ≤ e₂ := by linarith
    have hnoty₂ : ¬ ∃ n : ℤ, |(y : ℝ) * e₂ - n| < 2/25 := by
      rintro ⟨n, hn⟩
      exact absurd hn (not_lt.mpr (boundary_uncovered y hy n₂ n))
    have hcov₂ := hcov e₂ he₂a hb₂
    rcases hcov₂ with hcu | hcx | hcy'
    · exact (hkillU e₂ (by linarith) hb₂ hcu).elim
    · -- x-revisit at e₂: kill via the x-gap covered by {y, u}
      exfalso
      obtain ⟨n₃, hn₃⟩ := hcx
      have habs := abs_lt.mp hn₃
      have hxlt : (n₁ : ℝ) + 2/25 < (x : ℝ) * e₂ := by
        rw [← hxe₁]
        exact mul_lt_mul_of_pos_left he₂gt hxR
      have hmrev : n₁ + 1 ≤ n₃ := by
        have : (n₁ : ℝ) < (n₃ : ℝ) := by linarith
        have hmn : n₁ < n₃ := by exact_mod_cast this
        exact Int.add_one_le_iff.mpr hmn
      have hmb : ((n₃ : ℝ) - 2/25) / x ≤ b := by
        rw [div_le_iff₀ hxR]
        have h1 : ((n₃ : ℝ) - 2/25) < (x : ℝ) * e₂ := by linarith
        have h2 : (x : ℝ) * e₂ ≤ b * x := by
          rw [mul_comm]
          exact mul_le_mul_of_nonneg_right hb₂ hxR.le
        linarith
      apply revisit_kill x y u hx hy hu hxy hxu n₁ ?_ ⟨n₂, by rw [← he₁def]; exact hn₂⟩
      intro t ht1 ht2
      have hta : a ≤ t := by
        have : e₁ ≤ t := by rw [he₁def]; exact ht1
        linarith
      have htb : t ≤ b := by
        have h1 : ((n₁ : ℝ) + 1 - 2/25) / x ≤ ((n₃ : ℝ) - 2/25) / x := by
          apply div_le_div_of_nonneg_right ?_ hxR.le
          have : (n₁ : ℝ) + 1 ≤ (n₃ : ℝ) := by exact_mod_cast hmrev
          linarith
        linarith
      have hnotx : ¬ ∃ nn : ℤ, |(x : ℝ) * t - nn| < 2/25 := by
        rintro ⟨nn, hnn⟩
        exact absurd hnn (not_lt.mpr (in_gap_far x hx n₁ t ht1 ht2 nn))
      -- cover of the x-gap: u ∨ y from hcov, reordered to y ∨ u
      rcases hcov t hta htb with hcu | hcx' | hcy''
      · exact Or.inr hcu
      · exact absurd hcx' hnotx
      · exact Or.inl hcy''
    · exact absurd hcy' hnoty₂

/-- **THE TRIPLE WALK**: an interval covered by three pairwise-balanced combs
(`4·big ≤ 7·small` for every ordered pair), with the left endpoint on comb `u`,
spans less than one tooth of each comb. -/
theorem triple_walk (u v w : ℤ) (hu : 0 < u) (hv : 0 < v) (hw : 0 < w) (a b : ℝ)
    (huv : 4 * u ≤ 7 * v) (hvu : 4 * v ≤ 7 * u)
    (huw : 4 * u ≤ 7 * w) (hwu : 4 * w ≤ 7 * u)
    (hvw : 4 * v ≤ 7 * w) (hwv : 4 * w ≤ 7 * v)
    (hcov : ∀ t, a ≤ t → t ≤ b →
      (∃ n : ℤ, |(u : ℝ) * t - n| < 2/25) ∨ ((∃ n : ℤ, |(v : ℝ) * t - n| < 2/25) ∨
      (∃ n : ℤ, |(w : ℝ) * t - n| < 2/25)))
    (hstart : ∃ n : ℤ, |(u : ℝ) * a - n| < 2/25) :
    b - a < (4:ℝ)/25 / u + (4:ℝ)/25 / v + (4:ℝ)/25 / w := by
  have huR : (0 : ℝ) < (u : ℝ) := by exact_mod_cast hu
  have hposv : (0:ℝ) < (4:ℝ)/25 / v := by positivity
  have hposw : (0:ℝ) < (4:ℝ)/25 / w := by positivity
  obtain ⟨m₀, hm₀⟩ := hstart
  obtain ⟨he₀gt, he₀lt⟩ := step_data u hu a m₀ hm₀
  by_cases hb₀ : b < ((m₀ : ℝ) + 2/25) / u
  · linarith
  push_neg at hb₀
  have hnotu₀ : ¬ ∃ n : ℤ, |(u : ℝ) * (((m₀ : ℝ) + 2/25) / u) - n| < 2/25 := by
    rintro ⟨n, hn⟩
    exact absurd hn (not_lt.mpr (boundary_uncovered u hu m₀ n))
  have hcov₀ := (hcov _ (le_of_lt he₀gt) hb₀).resolve_left hnotu₀
  rcases hcov₀ with hcx | hcy
  · -- second comb v: tail with (u, v, w)
    exact walk_tail u v w hu hv hw huv huw hvw hvu a b hcov m₀ he₀gt he₀lt hb₀ hcx
  · -- second comb w: tail with (u, w, v), then reorder the conclusion
    have hcov' : ∀ t, a ≤ t → t ≤ b →
        (∃ n : ℤ, |(u : ℝ) * t - n| < 2/25) ∨ ((∃ n : ℤ, |(w : ℝ) * t - n| < 2/25) ∨
        (∃ n : ℤ, |(v : ℝ) * t - n| < 2/25)) := by
      intro t h1 h2
      rcases hcov t h1 h2 with h | h | h
      · exact Or.inl h
      · exact Or.inr (Or.inr h)
      · exact Or.inr (Or.inl h)
    have := walk_tail u w v hu hw hv huw huv hwv hwu a b hcov' m₀ he₀gt he₀lt hb₀ hcy
    linarith

/-- Sign absorption: covering by a runner ⟺ covering by its absolute speed. -/
lemma cov_abs_iff (z : ℤ) (t : ℝ) :
    (∃ n : ℤ, |(z : ℝ) * t - n| < 2/25) ↔
    (∃ n : ℤ, |((|z| : ℤ) : ℝ) * t - n| < 2/25) := by
  constructor <;> rintro ⟨n, hn⟩
  · rcases le_or_gt 0 z with hpos | hneg
    · exact ⟨n, by rwa [abs_of_nonneg hpos]⟩
    · refine ⟨-n, ?_⟩
      rw [abs_of_neg hneg]
      rw [show ((-z : ℤ) : ℝ) * t - ((-n : ℤ) : ℝ) = -((z : ℝ) * t - n) by push_cast; ring,
        abs_neg]
      exact hn
  · rcases le_or_gt 0 z with hpos | hneg
    · exact ⟨n, by rwa [abs_of_nonneg hpos] at hn⟩
    · refine ⟨-n, ?_⟩
      rw [abs_of_neg hneg] at hn
      rw [show ((-z : ℤ) : ℝ) * t - (n : ℝ) = -((z : ℝ) * t - ((-n : ℤ) : ℝ)) by
        push_cast; ring, abs_neg] at hn
      exact hn

/-- **THE BALANCED-TRIPLE RUNG (12B)**: in a 12-family with no `2/25`-margin
point, every pairwise-balanced 3-subset (`4·|v_p| ≤ 7·|v_q|` for all ordered
pairs) has `min ≤ 12·B`, `B` a bound on the other nine. -/
theorem gap_triple_rung (cite : LRCUpTo13) (v : Fin 12 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25)
    (i j k : Fin 12) (hij : i ≠ j) (hik : i ≠ k) (hjk : j ≠ k)
    (bij : 4 * |v i| ≤ 7 * |v j|) (bji : 4 * |v j| ≤ 7 * |v i|)
    (bik : 4 * |v i| ≤ 7 * |v k|) (bki : 4 * |v k| ≤ 7 * |v i|)
    (bjk : 4 * |v j| ≤ 7 * |v k|) (bkj : 4 * |v k| ≤ 7 * |v j|)
    (B : ℤ) (hB0 : 0 < B) (hB : ∀ n, n ≠ i → n ≠ j → n ≠ k → |v n| ≤ B) :
    min |v i| (min |v j| |v k|) ≤ 12 * B := by
  classical
  by_contra hbig
  push_neg at hbig
  have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB0
  -- the 9-complement, cited at 1/10
  set T : Finset (Fin 12) := ((Finset.univ.erase i).erase j).erase k with hT
  have hTcard : T.card = 9 := by
    rw [hT, Finset.card_erase_of_mem, Finset.card_erase_of_mem,
      Finset.card_erase_of_mem (Finset.mem_univ i)]
    · rw [Finset.card_univ, Fintype.card_fin]
    · exact Finset.mem_erase.mpr ⟨Ne.symm hij, Finset.mem_univ j⟩
    · exact Finset.mem_erase.mpr ⟨Ne.symm hjk,
        Finset.mem_erase.mpr ⟨Ne.symm hik, Finset.mem_univ k⟩⟩
  obtain ⟨t', hmargin⟩ := cite_margin_gen cite v T (by rw [hTcard]; norm_num)
    (fun n _ => hv n)
  have hmargin10 : ∀ n, n ≠ i → n ≠ j → n ≠ k → ∀ m : ℤ,
      (1:ℝ)/10 ≤ |(v n : ℝ) * t' - m| := by
    intro n hni hnj hnk m
    have hnT : n ∈ T := by
      rw [hT, Finset.mem_erase, Finset.mem_erase, Finset.mem_erase]
      exact ⟨hnk, hnj, hni, Finset.mem_univ n⟩
    have := hmargin n hnT m
    rw [hTcard] at this
    norm_num at this
    exact this
  set δ : ℝ := 1 / (50 * (B : ℝ)) with hδdef
  have hδ0 : 0 < δ := by rw [hδdef]; positivity
  have hTsafe : ∀ t, |t - t'| ≤ δ → ∀ n, n ≠ i → n ≠ j → n ≠ k → ∀ m : ℤ,
      (2:ℝ)/25 ≤ |(v n : ℝ) * t - m| := by
    intro t ht n hni hnj hnk m
    have htrans := margin_transfer
      (fun z : {z : Fin 12 // z ≠ i ∧ z ≠ j ∧ z ≠ k} => v z)
      t' t (2/25) (1/10) ((B : ℤ) : ℝ)
      (fun z => by
        rw [← Int.cast_abs]
        exact_mod_cast hB z z.2.1 z.2.2.1 z.2.2.2)
      (fun z m' => hmargin10 z z.2.1 z.2.2.1 z.2.2.2 m')
      (by
        calc ((B : ℤ) : ℝ) * |t - t'| ≤ (B : ℝ) * δ :=
              mul_le_mul_of_nonneg_left ht hBR.le
          _ = 1/50 := by rw [hδdef]; field_simp
          _ ≤ 1/10 - 2/25 := by norm_num)
    exact htrans ⟨n, hni, hnj, hnk⟩ m
  set wi : ℤ := |v i| with hwi
  set wj : ℤ := |v j| with hwj
  set wk : ℤ := |v k| with hwk
  have hwi0 : 0 < wi := abs_pos.mpr (hv i)
  have hwj0 : 0 < wj := abs_pos.mpr (hv j)
  have hwk0 : 0 < wk := abs_pos.mpr (hv k)
  have hcov : ∀ t, t' - δ ≤ t → t ≤ t' + δ →
      (∃ n : ℤ, |(wi : ℝ) * t - n| < 2/25) ∨ ((∃ n : ℤ, |(wj : ℝ) * t - n| < 2/25) ∨
      (∃ n : ℤ, |(wk : ℝ) * t - n| < 2/25)) := by
    intro t ht1 ht2
    have htδ : |t - t'| ≤ δ := by rw [abs_le]; constructor <;> linarith
    obtain ⟨n, m, hnm⟩ := hnl t
    by_cases hni : n = i
    · exact Or.inl ((cov_abs_iff (v i) t).mp ⟨m, by rwa [← hni]⟩)
    by_cases hnj : n = j
    · exact Or.inr (Or.inl ((cov_abs_iff (v j) t).mp ⟨m, by rwa [← hnj]⟩))
    by_cases hnk : n = k
    · exact Or.inr (Or.inr ((cov_abs_iff (v k) t).mp ⟨m, by rwa [← hnk]⟩))
    · exact absurd hnm (not_lt.mpr (hTsafe t htδ n hni hnj hnk m))
  -- min bounds and the final arithmetic, shared by the three start cases
  have hmin0 : (0 : ℝ) < (min wi (min wj wk) : ℝ) := by
    have := lt_min hwi0 (lt_min hwj0 hwk0)
    exact_mod_cast this
  have hminwi : (min wi (min wj wk) : ℝ) ≤ (wi : ℝ) := by
    exact_mod_cast min_le_left wi (min wj wk)
  have hminwj : (min wi (min wj wk) : ℝ) ≤ (wj : ℝ) := by
    exact_mod_cast le_trans (min_le_right wi (min wj wk)) (min_le_left wj wk)
  have hminwk : (min wi (min wj wk) : ℝ) ≤ (wk : ℝ) := by
    exact_mod_cast le_trans (min_le_right wi (min wj wk)) (min_le_right wj wk)
  have hwiR : (0:ℝ) < (wi:ℝ) := by exact_mod_cast hwi0
  have hwjR : (0:ℝ) < (wj:ℝ) := by exact_mod_cast hwj0
  have hwkR : (0:ℝ) < (wk:ℝ) := by exact_mod_cast hwk0
  have hfinish : ∀ (p q r : ℤ), 0 < p → 0 < q → 0 < r →
      (min wi (min wj wk) : ℝ) ≤ (p : ℝ) → (min wi (min wj wk) : ℝ) ≤ (q : ℝ) →
      (min wi (min wj wk) : ℝ) ≤ (r : ℝ) →
      (t' + δ) - (t' - δ) < (4:ℝ)/25 / p + (4:ℝ)/25 / q + (4:ℝ)/25 / r → False := by
    intro p q r hp hq hr hmp hmq hmr hwalk
    have hpR : (0:ℝ) < (p:ℝ) := by exact_mod_cast hp
    have hqR : (0:ℝ) < (q:ℝ) := by exact_mod_cast hq
    have hrR : (0:ℝ) < (r:ℝ) := by exact_mod_cast hr
    have hbp : (4:ℝ)/25 / p ≤ (4:ℝ)/25 / (min wi (min wj wk) : ℝ) :=
      div_le_div_of_nonneg_left (by norm_num) hmin0 hmp
    have hbq : (4:ℝ)/25 / q ≤ (4:ℝ)/25 / (min wi (min wj wk) : ℝ) :=
      div_le_div_of_nonneg_left (by norm_num) hmin0 hmq
    have hbr : (4:ℝ)/25 / r ≤ (4:ℝ)/25 / (min wi (min wj wk) : ℝ) :=
      div_le_div_of_nonneg_left (by norm_num) hmin0 hmr
    have h2δ : (t' + δ) - (t' - δ) = 2 * δ := by ring
    have hsum : 2 * δ < 3 * ((4:ℝ)/25 / (min wi (min wj wk) : ℝ)) := by
      rw [← h2δ]
      linarith
    have hδval : 2 * δ = 1 / (25 * (B:ℝ)) := by rw [hδdef]; ring
    rw [hδval, show (3:ℝ) * ((4:ℝ)/25 / (min wi (min wj wk) : ℝ))
        = (12/25) / (min wi (min wj wk) : ℝ) by ring,
      div_lt_div_iff₀ (by positivity) hmin0] at hsum
    have hminbig : (12 : ℝ) * B < (min wi (min wj wk) : ℝ) := by
      have hmm : 12 * B < min wi (min wj wk) := by
        rw [hwi, hwj, hwk] at *
        exact hbig
      exact_mod_cast hmm
    linarith
  -- the window's left endpoint is covered by one of the three: start the walk there
  rcases hcov (t' - δ) (le_refl _) (by linarith) with hs | hs | hs
  · exact hfinish wi wj wk hwi0 hwj0 hwk0 hminwi hminwj hminwk
      (triple_walk wi wj wk hwi0 hwj0 hwk0 (t' - δ) (t' + δ)
        bij bji bik bki bjk bkj hcov hs)
  · have hcov' : ∀ t, t' - δ ≤ t → t ≤ t' + δ →
        (∃ n : ℤ, |(wj : ℝ) * t - n| < 2/25) ∨ ((∃ n : ℤ, |(wi : ℝ) * t - n| < 2/25) ∨
        (∃ n : ℤ, |(wk : ℝ) * t - n| < 2/25)) := by
      intro t h1 h2
      rcases hcov t h1 h2 with h | h | h
      · exact Or.inr (Or.inl h)
      · exact Or.inl h
      · exact Or.inr (Or.inr h)
    exact hfinish wj wi wk hwj0 hwi0 hwk0 hminwj hminwi hminwk
      (triple_walk wj wi wk hwj0 hwi0 hwk0 (t' - δ) (t' + δ)
        bji bij bjk bkj bik bki hcov' hs)
  · have hcov' : ∀ t, t' - δ ≤ t → t ≤ t' + δ →
        (∃ n : ℤ, |(wk : ℝ) * t - n| < 2/25) ∨ ((∃ n : ℤ, |(wi : ℝ) * t - n| < 2/25) ∨
        (∃ n : ℤ, |(wj : ℝ) * t - n| < 2/25)) := by
      intro t h1 h2
      rcases hcov t h1 h2 with h | h | h
      · exact Or.inr (Or.inl h)
      · exact Or.inr (Or.inr h)
      · exact Or.inl h
    exact hfinish wk wi wj hwk0 hwi0 hwj0 hminwk hminwi hminwj
      (triple_walk wk wi wj hwk0 hwi0 hwj0 (t' - δ) (t' + δ)
        bki bik bkj bjk bij bji hcov' hs)

#print axioms triple_walk
#print axioms gap_triple_rung

end TripleWalk
end LonelyRunner
