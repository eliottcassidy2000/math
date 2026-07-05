/-
  TournamentH7.LRCPairWalk — THE PAIR WALK: the sharp l = 2 rung
  (kind-pasteur-2026-07-05-S8, HYP-4117).

  THE WALK: an interval `[a, b]` covered by TWO radius-`2/25` combs satisfies
  `min(w₁, w₂)·(b − a) ≤ 12/25`.  Proof is a 3-step boundary walk, measure-free
  and subcover-free:

  * a tooth's right end `e = (m + r)/w` sits at distance EXACTLY `r` from ℤ
    (times `w`), so it is never covered by its own comb — the cover must
    alternate combs at each boundary;
  * two same-comb teeth in the walk are distinct, so the comb advanced a full
    `1/w` while the other comb bridged with one tooth: `(1−2r)·u < 2r·v`;
  * the FOURTH boundary forces the symmetric imbalance `(1−2r)·v < 2r·u`;
    ADDING the two: `1 − 2r < 2r`, absurd at `r = 2/25` — pure `linarith`.
    So the walk exits `[a, b]` within three teeth: `b − a < 2r/u + 2r/v + 2r/u`.

  THE RUNG (`gap_pair_rung`): in a gap violator (no `2/25`-margin point), EVERY
  pair of runners has `min(|vᵢ|, |vⱼ|) ≤ 22·B` where `B` bounds the other ten —
  the 10-complement is cited at `1/11` (klein's `cite_margin_gen`), the excess
  `1/11 − 2/25 = 3/275` opens a window of length `6/(275B)`, and the walk caps
  the pair: `min ≤ (12/25)·(275B/6) = 22B`.  Three times sharper than the S7
  density rung `C₂ = 1100/17 ≈ 64.7` — and sharper than the S6 single-runner
  constant 24.  Complementary to opus-S81's gap descent (which dodges SPREAD
  tops; the walk kills BALANCED pairs).

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCGridAttainment
import TournamentH7.LRCCertCompleteness
import TournamentH7.LRCTeethR

namespace LonelyRunner
namespace PairWalk

open TournamentH7.LRCWitness GridAttainment CertCompleteness TeethR

/-- A tooth's right end is never covered by its own comb: at `e = (m₀ + r)/w`
every integer is at scaled distance `≥ r`. -/
private lemma boundary_uncovered (w : ℤ) (hw : 0 < w) (m₀ : ℤ) (m : ℤ) :
    (2 : ℝ) / 25 ≤ |(w : ℝ) * (((m₀ : ℝ) + 2/25) / w) - m| := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  have heq : (w : ℝ) * (((m₀ : ℝ) + 2/25) / w) = (m₀ : ℝ) + 2/25 := by
    field_simp
  rw [heq]
  rcases eq_or_ne m m₀ with rfl | hne
  · rw [show (m : ℝ) + 2/25 - m = 2/25 by ring]
    rw [abs_of_pos (by norm_num)]
  · have h1 : (1 : ℝ) ≤ |(m₀ : ℝ) - m| := by
      rw [← Int.cast_sub, ← Int.cast_abs]
      exact_mod_cast Int.one_le_abs (sub_ne_zero.mpr (Ne.symm hne))
    have htri : |(m₀ : ℝ) - m| ≤ |(m₀ : ℝ) + 2/25 - m| + 2/25 := by
      calc |(m₀ : ℝ) - m| = |((m₀ : ℝ) + 2/25 - m) + (-(2/25))| := by ring_nf
        _ ≤ |(m₀ : ℝ) + 2/25 - m| + |(-(2/25) : ℝ)| := abs_add_le _ _
        _ = |(m₀ : ℝ) + 2/25 - m| + 2/25 := by norm_num
    linarith

/-- One walk step: a point covered by comb `w` lies strictly left of its tooth's
right end `e := (round(w·p) + r)/w`, with `e − p < 2r/w`; and if the point sits at
a previous same-comb boundary height `(n + r)/w'`… (kept inline below). -/
private lemma step_data (w : ℤ) (hw : 0 < w) (p : ℝ) (m : ℤ)
    (hm : |(w : ℝ) * p - m| < 2/25) :
    p < ((m : ℝ) + 2/25) / w ∧ ((m : ℝ) + 2/25) / w - p < (4/25) / (w : ℝ) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  have h1 := abs_lt.mp hm
  constructor
  · rw [lt_div_iff₀ hwR]
    have : (w : ℝ) * p < (m : ℝ) + 2/25 := by linarith [h1.2]
    linarith [this]
  · rw [lt_div_iff₀ hwR, sub_mul, div_mul_cancel₀ _ (ne_of_gt hwR)]
    have hwp : (m : ℝ) - 2/25 < (w : ℝ) * p := by linarith [h1.1]
    nlinarith [hwp, hwR]

/-- **THE PAIR WALK BOUND**: an interval `[a, b]` on which every point is within
`2/25` of an integer (scaled) for comb `u` OR comb `v` — with the left endpoint
covered by `u` — satisfies `b − a < 2·(4/25)/u + (4/25)/v`. -/
theorem walk_core (u v : ℤ) (hu : 0 < u) (hv : 0 < v) (a b : ℝ)
    (hcov : ∀ t, a ≤ t → t ≤ b →
      (∃ m : ℤ, |(u : ℝ) * t - m| < 2/25) ∨ (∃ m : ℤ, |(v : ℝ) * t - m| < 2/25))
    (hstart : ∃ m : ℤ, |(u : ℝ) * a - m| < 2/25) :
    b - a < (4:ℝ)/25 / u + (4:ℝ)/25 / u + (4:ℝ)/25 / v := by
  have huR : (0 : ℝ) < (u : ℝ) := by exact_mod_cast hu
  have hvR : (0 : ℝ) < (v : ℝ) := by exact_mod_cast hv
  have hposu : (0:ℝ) < (4:ℝ)/25 / u := by positivity
  have hposv : (0:ℝ) < (4:ℝ)/25 / v := by positivity
  obtain ⟨m₀, hm₀⟩ := hstart
  obtain ⟨he₀gt, he₀lt⟩ := step_data u hu a m₀ hm₀
  set e₀ : ℝ := ((m₀ : ℝ) + 2/25) / u with he₀
  by_cases hb₀ : b < e₀
  · linarith
  push_neg at hb₀
  have hcov₀ := hcov e₀ (le_of_lt he₀gt) hb₀
  have hnotu₀ : ¬ ∃ m : ℤ, |(u : ℝ) * e₀ - m| < 2/25 := by
    rintro ⟨m, hm⟩
    exact absurd hm (not_lt.mpr (boundary_uncovered u hu m₀ m))
  obtain ⟨m₁, hm₁⟩ := hcov₀.resolve_left hnotu₀
  obtain ⟨he₁gt, he₁lt⟩ := step_data v hv e₀ m₁ hm₁
  set e₁ : ℝ := ((m₁ : ℝ) + 2/25) / v with he₁
  by_cases hb₁ : b < e₁
  · linarith
  push_neg at hb₁
  have he₁a : a ≤ e₁ := by linarith
  have hcov₁ := hcov e₁ he₁a hb₁
  have hnotv₁ : ¬ ∃ m : ℤ, |(v : ℝ) * e₁ - m| < 2/25 := by
    rintro ⟨m, hm⟩
    exact absurd hm (not_lt.mpr (boundary_uncovered v hv m₁ m))
  obtain ⟨m₂, hm₂⟩ := hcov₁.resolve_right hnotv₁
  obtain ⟨he₂gt, he₂lt⟩ := step_data u hu e₁ m₂ hm₂
  set e₂ : ℝ := ((m₂ : ℝ) + 2/25) / u with he₂
  -- the u-comb advanced a full tooth
  have hue₀ : (u : ℝ) * e₀ = (m₀ : ℝ) + 2/25 := by
    rw [he₀]; field_simp
  have hue₁gt : (m₀ : ℝ) + 2/25 < (u : ℝ) * e₁ := by
    rw [← hue₀]
    exact mul_lt_mul_of_pos_left he₁gt huR
  have hm₂m₀ : m₀ + 1 ≤ m₂ := by
    have h2 := (abs_lt.mp hm₂).2
    have : (m₀ : ℝ) < (m₂ : ℝ) := by linarith
    exact_mod_cast Int.add_one_le_iff.mpr (by exact_mod_cast this)
  have himbB : (1 - 4/25 : ℝ) * v < (4/25) * u := by
    have h1 := (abs_lt.mp hm₂).1
    have hm₂R : (m₀ : ℝ) + 1 ≤ (m₂ : ℝ) := by exact_mod_cast hm₂m₀
    have hgap : (1 - 4/25 : ℝ) < (u : ℝ) * (e₁ - e₀) := by
      have hexp : (u : ℝ) * (e₁ - e₀) = (u : ℝ) * e₁ - ((m₀ : ℝ) + 2/25) := by
        rw [mul_sub, hue₀]
      rw [hexp]
      linarith
    have hstep := mul_lt_mul_of_pos_left he₁lt huR
    have h2 : (1 - 4/25 : ℝ) < (u : ℝ) * ((4:ℝ)/25 / v) := by linarith
    rw [show (u:ℝ) * ((4:ℝ)/25 / v) = ((4:ℝ)/25 * u) / v by ring,
      lt_div_iff₀ hvR] at h2
    linarith
  by_cases hb₂ : b < e₂
  · linarith
  push_neg at hb₂
  have he₂a : a ≤ e₂ := by linarith
  have hcov₂ := hcov e₂ he₂a hb₂
  have hnotu₂ : ¬ ∃ m : ℤ, |(u : ℝ) * e₂ - m| < 2/25 := by
    rintro ⟨m, hm⟩
    exact absurd hm (not_lt.mpr (boundary_uncovered u hu m₂ m))
  obtain ⟨m₃, hm₃⟩ := hcov₂.resolve_left hnotu₂
  -- the v-comb advanced a full tooth: the symmetric imbalance
  have hve₁ : (v : ℝ) * e₁ = (m₁ : ℝ) + 2/25 := by
    rw [he₁]; field_simp
  have hve₂gt : (m₁ : ℝ) + 2/25 < (v : ℝ) * e₂ := by
    rw [← hve₁]
    exact mul_lt_mul_of_pos_left he₂gt hvR
  have hm₃m₁ : m₁ + 1 ≤ m₃ := by
    have h2 := (abs_lt.mp hm₃).2
    have : (m₁ : ℝ) < (m₃ : ℝ) := by linarith
    exact_mod_cast Int.add_one_le_iff.mpr (by exact_mod_cast this)
  have himbA : (1 - 4/25 : ℝ) * u < (4/25) * v := by
    have h1 := (abs_lt.mp hm₃).1
    have hm₃R : (m₁ : ℝ) + 1 ≤ (m₃ : ℝ) := by exact_mod_cast hm₃m₁
    have hgap : (1 - 4/25 : ℝ) < (v : ℝ) * (e₂ - e₁) := by
      have hexp : (v : ℝ) * (e₂ - e₁) = (v : ℝ) * e₂ - ((m₁ : ℝ) + 2/25) := by
        rw [mul_sub, hve₁]
      rw [hexp]
      linarith
    have hstep := mul_lt_mul_of_pos_left he₂lt hvR
    have h2 : (1 - 4/25 : ℝ) < (v : ℝ) * ((4:ℝ)/25 / u) := by linarith
    rw [show (v:ℝ) * ((4:ℝ)/25 / u) = ((4:ℝ)/25 * v) / u by ring,
      lt_div_iff₀ huR] at h2
    linarith
  -- THE KILL: add the imbalances
  linarith [himbA, himbB, huR, hvR]

/-- **THE SHARP PAIR RUNG (22B)**: in a 12-family with no `2/25`-margin point,
EVERY pair of runners has `min(|vᵢ|, |vⱼ|) ≤ 22·B` where `B` bounds the other
ten.  (S7's density rung gave `≈ 64.7·B`; the walk is 3× sharper.) -/
theorem gap_pair_rung (cite : LRCUpTo13) (v : Fin 12 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25)
    (i j : Fin 12) (hij : i ≠ j) (B : ℤ) (hB0 : 0 < B)
    (hB : ∀ k, k ≠ i → k ≠ j → |v k| ≤ B) :
    min |v i| |v j| ≤ 22 * B := by
  classical
  by_contra hbig
  push_neg at hbig
  have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB0
  -- the 10-complement, cited at 1/11
  set T : Finset (Fin 12) := (Finset.univ.erase i).erase j with hT
  have hTcard : T.card = 10 := by
    rw [hT, Finset.card_erase_of_mem, Finset.card_erase_of_mem (Finset.mem_univ i)]
    · rw [Finset.card_univ, Fintype.card_fin]
    · exact Finset.mem_erase.mpr ⟨Ne.symm hij, Finset.mem_univ j⟩
  have hTmem : ∀ k, k ∈ T ↔ (k ≠ j ∧ k ≠ i) := by
    intro k
    rw [hT, Finset.mem_erase, Finset.mem_erase]
    simp only [Finset.mem_univ, and_true]
  obtain ⟨t', hmargin⟩ := cite_margin_gen cite v T (by rw [hTcard]; norm_num)
    (fun k _ => hv k)
  have hmargin11 : ∀ k, k ≠ i → k ≠ j → ∀ m : ℤ, (1:ℝ)/11 ≤ |(v k : ℝ) * t' - m| := by
    intro k hki hkj m
    have hkT : k ∈ T := (hTmem k).mpr ⟨hkj, hki⟩
    have := hmargin k hkT m
    rw [hTcard] at this
    norm_num at this
    exact this
  -- the window
  set δ : ℝ := 3 / (275 * (B : ℝ)) with hδdef
  have hδ0 : 0 < δ := by rw [hδdef]; positivity
  -- complement margin ≥ 2/25 on the window (Lipschitz transfer)
  have hTsafe : ∀ t, |t - t'| ≤ δ → ∀ k, k ≠ i → k ≠ j → ∀ m : ℤ,
      (2:ℝ)/25 ≤ |(v k : ℝ) * t - m| := by
    intro t ht k hki hkj m
    have htrans := margin_transfer (fun x : {x : Fin 12 // x ≠ i ∧ x ≠ j} => v x)
      t' t (2/25) (1/11) ((B : ℤ) : ℝ)
      (fun x => by
        rw [← Int.cast_abs]
        exact_mod_cast hB x x.2.1 x.2.2)
      (fun x m' => hmargin11 x x.2.1 x.2.2 m')
      (by
        calc ((B : ℤ) : ℝ) * |t - t'| ≤ (B : ℝ) * δ :=
              mul_le_mul_of_nonneg_left ht hBR.le
          _ = 3/275 := by rw [hδdef]; field_simp
          _ ≤ 1/11 - 2/25 := by norm_num)
    exact htrans ⟨k, hki, hkj⟩ m
  -- on the window, every point is covered by the pair {i, j}
  set wi : ℤ := |v i| with hwi
  set wj : ℤ := |v j| with hwj
  have hwi0 : 0 < wi := abs_pos.mpr (hv i)
  have hwj0 : 0 < wj := abs_pos.mpr (hv j)
  have habsi : ∀ t, (∃ m : ℤ, |(v i : ℝ) * t - m| < 2/25) ↔
      (∃ m : ℤ, |(wi : ℝ) * t - m| < 2/25) := by
    intro t
    constructor <;> rintro ⟨m, hm⟩
    · rcases le_or_gt 0 (v i) with hpos | hneg
      · exact ⟨m, by rwa [hwi, abs_of_nonneg hpos]⟩
      · refine ⟨-m, ?_⟩
        rw [hwi, abs_of_neg hneg]
        rw [show (↑(-v i) : ℝ) * t - ↑(-m) = -((v i : ℝ) * t - m) by push_cast; ring,
          abs_neg]
        exact hm
    · rcases le_or_gt 0 (v i) with hpos | hneg
      · exact ⟨m, by rwa [hwi, abs_of_nonneg hpos] at hm⟩
      · refine ⟨-m, ?_⟩
        rw [hwi, abs_of_neg hneg] at hm
        rw [show (↑(-v i) : ℝ) * t - ↑m = -((v i : ℝ) * t - (-m : ℤ)) by push_cast; ring,
          abs_neg] at hm
        exact hm
  have habsj : ∀ t, (∃ m : ℤ, |(v j : ℝ) * t - m| < 2/25) ↔
      (∃ m : ℤ, |(wj : ℝ) * t - m| < 2/25) := by
    intro t
    constructor <;> rintro ⟨m, hm⟩
    · rcases le_or_gt 0 (v j) with hpos | hneg
      · exact ⟨m, by rwa [hwj, abs_of_nonneg hpos]⟩
      · refine ⟨-m, ?_⟩
        rw [hwj, abs_of_neg hneg]
        rw [show (↑(-v j) : ℝ) * t - ↑(-m) = -((v j : ℝ) * t - m) by push_cast; ring,
          abs_neg]
        exact hm
    · rcases le_or_gt 0 (v j) with hpos | hneg
      · exact ⟨m, by rwa [hwj, abs_of_nonneg hpos] at hm⟩
      · refine ⟨-m, ?_⟩
        rw [hwj, abs_of_neg hneg] at hm
        rw [show (↑(-v j) : ℝ) * t - ↑m = -((v j : ℝ) * t - (-m : ℤ)) by push_cast; ring,
          abs_neg] at hm
        exact hm
  have hcov : ∀ t, t' - δ ≤ t → t ≤ t' + δ →
      (∃ m : ℤ, |(wi : ℝ) * t - m| < 2/25) ∨ (∃ m : ℤ, |(wj : ℝ) * t - m| < 2/25) := by
    intro t ht1 ht2
    have htδ : |t - t'| ≤ δ := by rw [abs_le]; constructor <;> linarith
    obtain ⟨k, m, hkm⟩ := hnl t
    by_cases hki : k = i
    · left; exact (habsi t).mp ⟨m, by rwa [← hki]⟩
    by_cases hkj : k = j
    · right; exact (habsj t).mp ⟨m, by rwa [← hkj]⟩
    · exact absurd hkm (not_lt.mpr (hTsafe t htδ k hki hkj m))
  -- the left endpoint is covered by one of the two: walk with that comb first
  have hwindow : (2:ℝ) * δ = 6 / (275 * (B:ℝ)) := by rw [hδdef]; ring
  have hlen : ∀ (hcase : (∃ m : ℤ, |(wi : ℝ) * (t' - δ) - m| < 2/25) ∨
      (∃ m : ℤ, |(wj : ℝ) * (t' - δ) - m| < 2/25)), False := by
    rintro (hstart | hstart)
    · have hwalk := walk_core wi wj hwi0 hwj0 (t' - δ) (t' + δ) hcov hstart
      -- b − a = 2δ; min ≤ wi, wj; 2δ < 2(4/25)/wi + (4/25)/wj ≤ 3(4/25)/min
      have hminwi : (min wi wj : ℝ) ≤ (wi : ℝ) := by exact_mod_cast min_le_left wi wj
      have hminwj : (min wi wj : ℝ) ≤ (wj : ℝ) := by exact_mod_cast min_le_right wi wj
      have hmin0 : (0 : ℝ) < (min wi wj : ℝ) := by
        exact_mod_cast lt_min hwi0 hwj0
      have hwiR : (0:ℝ) < (wi:ℝ) := by exact_mod_cast hwi0
      have hwjR : (0:ℝ) < (wj:ℝ) := by exact_mod_cast hwj0
      have hb1 : ((4:ℝ)/25) / wi ≤ ((4:ℝ)/25) / (min wi wj : ℝ) := by
        apply div_le_div_of_nonneg_left (by norm_num) hmin0 hminwi
      have hb2 : ((4:ℝ)/25) / wj ≤ ((4:ℝ)/25) / (min wi wj : ℝ) := by
        apply div_le_div_of_nonneg_left (by norm_num) hmin0 hminwj
      have hsum : (t' + δ) - (t' - δ) < 3 * (((4:ℝ)/25) / (min wi wj : ℝ)) := by
        calc (t' + δ) - (t' - δ)
            < (4:ℝ)/25 / wi + (4:ℝ)/25 / wi + (4:ℝ)/25 / wj := hwalk
          _ ≤ 3 * (((4:ℝ)/25) / (min wi wj : ℝ)) := by linarith
      -- 2δ < (12/25)/min ⟹ min < (12/25)/(2δ) = 22B — contradict min > 22B
      have h2δ : (t' + δ) - (t' - δ) = 2 * δ := by ring
      rw [h2δ, hwindow] at hsum
      have hminbig : (22 : ℝ) * B < (min wi wj : ℝ) := by
        have := hbig
        have hmm : 22 * B < min wi wj := by
          rw [hwi, hwj] at *
          exact hbig
        exact_mod_cast hmm
      -- 6/(275B) < 3(4/25)/min ⟹ min·6/(275B) < 12/25 ⟹ min < 22B: contradiction
      rw [show (3:ℝ) * ((4:ℝ)/25 / (min (wi:ℝ) (wj:ℝ))) = ((12:ℝ)/25) / (min (wi:ℝ) (wj:ℝ)) by
        ring, div_lt_div_iff₀ (by positivity) hmin0] at hsum
      linarith [hsum, hminbig, hBR]
    · have hcov' : ∀ t, t' - δ ≤ t → t ≤ t' + δ →
          (∃ m : ℤ, |(wj : ℝ) * t - m| < 2/25) ∨ (∃ m : ℤ, |(wi : ℝ) * t - m| < 2/25) := by
        intro t h1 h2
        exact (hcov t h1 h2).symm
      have hwalk := walk_core wj wi hwj0 hwi0 (t' - δ) (t' + δ) hcov' hstart
      have hminwi : (min wi wj : ℝ) ≤ (wi : ℝ) := by exact_mod_cast min_le_left wi wj
      have hminwj : (min wi wj : ℝ) ≤ (wj : ℝ) := by exact_mod_cast min_le_right wi wj
      have hmin0 : (0 : ℝ) < (min wi wj : ℝ) := by
        exact_mod_cast lt_min hwi0 hwj0
      have hwiR : (0:ℝ) < (wi:ℝ) := by exact_mod_cast hwi0
      have hwjR : (0:ℝ) < (wj:ℝ) := by exact_mod_cast hwj0
      have hb1 : ((4:ℝ)/25) / wi ≤ ((4:ℝ)/25) / (min wi wj : ℝ) := by
        apply div_le_div_of_nonneg_left (by norm_num) hmin0 hminwi
      have hb2 : ((4:ℝ)/25) / wj ≤ ((4:ℝ)/25) / (min wi wj : ℝ) := by
        apply div_le_div_of_nonneg_left (by norm_num) hmin0 hminwj
      have hsum : (t' + δ) - (t' - δ) < 3 * (((4:ℝ)/25) / (min wi wj : ℝ)) := by
        calc (t' + δ) - (t' - δ)
            < (4:ℝ)/25 / wj + (4:ℝ)/25 / wj + (4:ℝ)/25 / wi := hwalk
          _ ≤ 3 * (((4:ℝ)/25) / (min wi wj : ℝ)) := by linarith
      have h2δ : (t' + δ) - (t' - δ) = 2 * δ := by ring
      rw [h2δ, hwindow] at hsum
      have hminbig : (22 : ℝ) * B < (min wi wj : ℝ) := by
        have hmm : 22 * B < min wi wj := by
          rw [hwi, hwj] at *
          exact hbig
        exact_mod_cast hmm
      rw [show (3:ℝ) * ((4:ℝ)/25 / (min (wi:ℝ) (wj:ℝ))) = ((12:ℝ)/25) / (min (wi:ℝ) (wj:ℝ)) by
        ring, div_lt_div_iff₀ (by positivity) hmin0] at hsum
      linarith [hsum, hminbig, hBR]
  exact hlen (hcov (t' - δ) (le_refl _) (by linarith))

#print axioms walk_core
#print axioms gap_pair_rung

end PairWalk
end LonelyRunner
