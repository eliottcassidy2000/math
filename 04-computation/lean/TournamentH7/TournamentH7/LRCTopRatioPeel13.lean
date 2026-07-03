/-
  TournamentH7.LRCTopRatioPeel13 — THE 13-RATIO PEEL (klein-2026-07-02-S114, HYP-4019).

  SHARPENS kps-S18's `top_ratio_lonely` (LRCTopRatioPeel.lean) from ratio 91 to ratio 13.

  kps's mechanism sends the far runner to a HALF-INTEGER (distance 1/2 from ℤ), which needs its
  dial to span a FULL unit over the transport window, forcing |v i0| ≥ 91·B.  But `Lonely 14`
  only requires distance ≥ 1/14, so the far runner's dial need only span the danger width
  2·(1/14) = 1/7 to be guaranteed a safe phase — forcing only |v i0| ≥ 91/7 · B = 13·B.

  MECHANISM (same transport, sharper far step):
   * the LRC(≤13) citation on the other TWELVE runners gives t₀, 1/13-lonely for them;
   * the 1/13 − 1/182 = 1/14 margin transports over [t₀−δ, t₀+δ], δ = 1/(182·B);
   * the far runner's phase sweeps an interval of length |v i0|·2δ ≥ 13·B·2/(182B) = 1/7, and
     EVERY real interval of length ≥ 1/7 contains a point at distance ≥ 1/14 from ℤ
     (`safe_point_in_interval`) — so a τ in the window has the far runner 1/14-safe too.

  CONSEQUENCE: `lrc14_of_compressed_lonely_13` — the remaining compressed families have top
  ratio below 13 (not 91), a 7× smaller gap.  Sorry-free.
-/
import Mathlib
import TournamentH7.LRC13Citation
import TournamentH7.LRC14CertRoute

namespace LonelyRunner
namespace LRC14

/-- For any integer `d` and any fractional offset `f ∈ [1/14, 13/14]`, the point `d + f` is at
distance at least `1/14` from every integer (here from `0`; used with `d = ⌊·⌋ - m`). -/
theorem dist_ge_of_fract (d : ℤ) (f : ℝ) (h1 : (1:ℝ)/14 ≤ f) (h2 : f ≤ 13/14) :
    (1:ℝ)/14 ≤ |(d : ℝ) + f| := by
  rcases lt_trichotomy d 0 with hd | hd | hd
  · have hdR : (d : ℝ) ≤ -1 := by exact_mod_cast (by omega : d ≤ -1)
    rw [abs_of_nonpos (by linarith)]; linarith
  · subst hd; rw [Int.cast_zero, zero_add, abs_of_nonneg (by linarith)]; linarith
  · have hdR : (1 : ℝ) ≤ (d : ℝ) := by exact_mod_cast (by omega : (1:ℤ) ≤ d)
    rw [abs_of_nonneg (by linarith)]; linarith

/-- The boundary safe point `d + 1/14` is at distance `≥ 1/14` from every integer. -/
theorem frac14_dist (d : ℤ) : (1:ℝ)/14 ≤ |(d : ℝ) + 1/14| :=
  dist_ge_of_fract d (1/14) (le_refl _) (by norm_num)

/-- **Every interval of length ≥ 1/7 contains a point at distance ≥ 1/14 from ℤ.**  The danger
set `dist(·,ℤ) < 1/14` is a union of open intervals of width `1/7`, separated by safe gaps, so an
interval of length `≥ 1/7` cannot be swallowed by one danger interval. -/
theorem safe_point_in_interval (lo hi : ℝ) (hlen : (1:ℝ)/7 ≤ hi - lo) :
    ∃ φ : ℝ, lo ≤ φ ∧ φ ≤ hi ∧ ∀ m : ℤ, (1:ℝ)/14 ≤ |φ - m| := by
  have hfa : (⌊lo⌋ : ℝ) + Int.fract lo = lo := Int.floor_add_fract lo
  have hf0 : 0 ≤ Int.fract lo := Int.fract_nonneg lo
  have hf1 : Int.fract lo < 1 := Int.fract_lt_one lo
  set f : ℝ := Int.fract lo with hfdef
  -- hfa : ↑⌊lo⌋ + f = lo  (used in linarith as an atom identity; never rewritten into ⌊lo⌋)
  rcases lt_or_ge f (1/14) with hflo | hfge
  · -- f < 1/14 : safe point = ⌊lo⌋ + 1/14
    refine ⟨(⌊lo⌋ : ℝ) + 1/14, by linarith, by linarith, ?_⟩
    intro m
    have hrw : ((⌊lo⌋ : ℝ) + 1/14) - (m : ℝ) = ((⌊lo⌋ - m : ℤ) : ℝ) + 1/14 := by push_cast; ring
    rw [hrw]; exact frac14_dist _
  · rcases le_or_gt f (13/14) with hfmid | hfhi
    · -- f ∈ [1/14, 13/14] : lo itself is safe
      refine ⟨lo, le_refl _, by linarith, ?_⟩
      intro m
      have hrw : lo - (m : ℝ) = ((⌊lo⌋ - m : ℤ) : ℝ) + f := by push_cast; linarith [hfa]
      rw [hrw]; exact dist_ge_of_fract _ f hfge hfmid
    · -- f > 13/14 : safe point = ⌊lo⌋ + 1 + 1/14
      refine ⟨(⌊lo⌋ : ℝ) + 1 + 1/14, by linarith, by linarith, ?_⟩
      intro m
      have hrw : ((⌊lo⌋ : ℝ) + 1 + 1/14) - (m : ℝ) = ((⌊lo⌋ + 1 - m : ℤ) : ℝ) + 1/14 := by
        push_cast; ring
      rw [hrw]; exact frac14_dist _

/-- **THE 13-RATIO PEEL.**  If one runner is at least `13×` faster than every other, the family
is 14-lonely — from the LRC(≤13) citation node alone.  (Sharpens kps's 91.) -/
theorem top_ratio_lonely_13 (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (i0 : Fin 13) {B : ℤ} (hB : 0 < B) (hrest : ∀ i, i ≠ i0 → |v i| ≤ B)
    (hfar : 13 * B ≤ |v i0|) : ∃ t : ℝ, Lonely 14 v t := by
  set w : Fin 12 → ℤ := fun j => v (i0.succAbove j) with hw
  have hwne : ∀ j, w j ≠ 0 := fun j => hv _
  obtain ⟨t₀, h13⟩ := cite 12 (by omega) w hwne
  have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB
  set δ : ℝ := 1 / (182 * (B : ℝ)) with hδ
  have hδpos : 0 < δ := by positivity
  -- far runner speed as a real, and its absolute value
  have hvi0R : (v i0 : ℝ) ≠ 0 := by exact_mod_cast hv i0
  have hWfar : 13 * (B : ℝ) ≤ |(v i0 : ℝ)| := by rw [← Int.cast_abs]; exact_mod_cast hfar
  -- swept phase length ≥ 1/7
  have hswept : (1:ℝ)/7 ≤ |(v i0 : ℝ)| * (2 * δ) := by
    rw [hδ]
    have h1 : |(v i0 : ℝ)| * (2 * (1 / (182 * (B : ℝ)))) = |(v i0 : ℝ)| / (91 * (B : ℝ)) := by
      field_simp; ring
    rw [h1, le_div_iff₀ (by positivity)]
    have hc : (1:ℝ)/7 * (91 * (B : ℝ)) = 13 * (B : ℝ) := by ring
    linarith [hWfar]
  -- get a τ in the window whose far-runner phase is 1/14-safe
  obtain ⟨τ, hτlo, hτhi, hτsafe⟩ :
      ∃ τ : ℝ, t₀ - δ ≤ τ ∧ τ ≤ t₀ + δ ∧ ∀ m : ℤ, (1:ℝ)/14 ≤ |(v i0 : ℝ) * τ - m| := by
    rcases lt_or_gt_of_ne hvi0R with hneg | hpos
    · -- v i0 < 0 : phase interval [v i0 (t₀+δ), v i0 (t₀-δ)]
      have hlen : (1:ℝ)/7 ≤ (v i0 : ℝ) * (t₀ - δ) - (v i0 : ℝ) * (t₀ + δ) := by
        have : (v i0 : ℝ) * (t₀ - δ) - (v i0 : ℝ) * (t₀ + δ) = |(v i0 : ℝ)| * (2 * δ) := by
          rw [abs_of_neg hneg]; ring
        rw [this]; exact hswept
      obtain ⟨φ, hφlo, hφhi, hφsafe⟩ :=
        safe_point_in_interval ((v i0 : ℝ) * (t₀ + δ)) ((v i0 : ℝ) * (t₀ - δ)) hlen
      refine ⟨φ / (v i0 : ℝ), ?_, ?_, ?_⟩
      · rw [le_div_iff_of_neg hneg]; linarith
      · rw [div_le_iff_of_neg hneg]; linarith
      · intro m
        have hcancel : (v i0 : ℝ) * (φ / (v i0 : ℝ)) = φ := by field_simp
        rw [hcancel]; exact hφsafe m
    · -- v i0 > 0 : phase interval [v i0 (t₀-δ), v i0 (t₀+δ)]
      have hlen : (1:ℝ)/7 ≤ (v i0 : ℝ) * (t₀ + δ) - (v i0 : ℝ) * (t₀ - δ) := by
        have : (v i0 : ℝ) * (t₀ + δ) - (v i0 : ℝ) * (t₀ - δ) = |(v i0 : ℝ)| * (2 * δ) := by
          rw [abs_of_pos hpos]; ring
        rw [this]; exact hswept
      obtain ⟨φ, hφlo, hφhi, hφsafe⟩ :=
        safe_point_in_interval ((v i0 : ℝ) * (t₀ - δ)) ((v i0 : ℝ) * (t₀ + δ)) hlen
      refine ⟨φ / (v i0 : ℝ), ?_, ?_, ?_⟩
      · rw [le_div_iff₀ hpos]; linarith
      · rw [div_le_iff₀ hpos]; linarith
      · intro m
        have hcancel : (v i0 : ℝ) * (φ / (v i0 : ℝ)) = φ := by field_simp
        rw [hcancel]; exact hφsafe m
  -- assemble Lonely 14 v τ
  refine ⟨τ, ?_⟩
  intro i m
  rcases eq_or_ne i i0 with rfl | hne
  · exact hτsafe m
  · -- transport of the 1/13 margin over speeds ≤ B
    obtain ⟨j', hj'⟩ := Fin.exists_succAbove_eq hne
    have h0 : (1 : ℝ) / 13 ≤ |(v i : ℝ) * t₀ - m| := by
      have := h13 j' m
      rw [hw] at this
      simpa [hj'] using this
    have hvB : |(v i : ℝ)| ≤ (B : ℝ) := by
      rw [← Int.cast_abs]; exact_mod_cast hrest i hne
    have htri : |(v i : ℝ) * t₀ - m| ≤ |(v i : ℝ) * τ - m| + |(v i : ℝ)| * |t₀ - τ| := by
      calc |(v i : ℝ) * t₀ - m|
          = |((v i : ℝ) * τ - m) + (v i : ℝ) * (t₀ - τ)| := by congr 1; ring
        _ ≤ |(v i : ℝ) * τ - m| + |(v i : ℝ) * (t₀ - τ)| := abs_add_le _ _
        _ = |(v i : ℝ) * τ - m| + |(v i : ℝ)| * |t₀ - τ| := by rw [abs_mul]
    have hwin : |t₀ - τ| ≤ δ := by rw [abs_le]; constructor <;> linarith
    have hlip : |(v i : ℝ)| * |t₀ - τ| ≤ (B : ℝ) * δ :=
      mul_le_mul hvB hwin (abs_nonneg _) hBR.le
    have hBδ : (B : ℝ) * δ = 1 / 182 := by rw [hδ]; field_simp
    calc (1 : ℝ) / (14 : ℕ) = 1 / 13 - 1 / 182 := by norm_num
      _ ≤ |(v i : ℝ) * τ - m| := by rw [hBδ] at hlip; linarith

/-- Convenience form with the max-entry hypothesis (matches kps's `top_ratio_lonely'`). -/
theorem top_ratio_lonely_13' (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (i0 : Fin 13) (hpeel : ∀ i, i ≠ i0 → 13 * |v i| ≤ |v i0|) :
    ∃ t : ℝ, Lonely 14 v t := by
  -- B = max over i ≠ i0 of |v i|; use B := |v i0| / 13 is awkward — take B as the max.
  -- Simpler: set B := (Finset.univ.erase i0).sup' ... ; but reuse the general form with any B ≥ all.
  classical
  -- pick B = max of the other |v i| (nonempty since 13 ≥ 2)
  set others : Finset (Fin 13) := Finset.univ.erase i0 with hoth
  have hne : others.Nonempty := by
    rw [hoth]; exact ⟨i0.succAbove 0, Finset.mem_erase.mpr ⟨(i0.succAbove_ne 0), Finset.mem_univ _⟩⟩
  set B : ℤ := others.sup' hne (fun i => |v i|) with hBdef
  have hBpos : 0 < B := by
    obtain ⟨i, hi, hlt⟩ : ∃ i ∈ others, 0 < |v i| := by
      obtain ⟨i, hi⟩ := hne
      exact ⟨i, hi, abs_pos.mpr (hv i)⟩
    rw [hBdef]; exact lt_of_lt_of_le hlt (Finset.le_sup' (fun j => |v j|) hi)
  have hrest : ∀ i, i ≠ i0 → |v i| ≤ B := by
    intro i hi
    rw [hBdef]; exact Finset.le_sup' (fun j => |v j|) (Finset.mem_erase.mpr ⟨hi, Finset.mem_univ i⟩)
  have hfar : 13 * B ≤ |v i0| := by
    obtain ⟨istar, histar, hEq⟩ := Finset.exists_mem_eq_sup' hne (fun i => |v i|)
    rw [hBdef, hEq]
    exact hpeel istar (Finset.mem_erase.mp histar).1
  exact top_ratio_lonely_13 cite v hv i0 hBpos hrest hfar

/-- **LRC(14) from the citation + window census + the COMPRESSED (ratio < 13) families.**
Strictly leaner remaining hypothesis than kps's ratio-91 form (7× smaller compressed band). -/
theorem lrc14_of_compressed_lonely_13 (cite : LRCUpTo13)
    (hwindow : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → (∀ i, |v i| ≤ 22) →
      ∃ t : ℝ, Lonely 14 v t)
    (hcompressed : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      ¬ (∀ i, |v i| ≤ 22) →
      (∀ i0 : Fin 13, ∃ i, i ≠ i0 ∧ |v i0| < 13 * |v i|) →
      ∃ t : ℝ, Lonely 14 v t) :
    LRC14Statement := by
  apply lrc14_of_covering_lonely
  intro v hv hcov
  by_cases hbnd : ∀ i, |v i| ≤ 22
  · exact hwindow v hv hbnd
  · by_cases hpeel : ∃ i0 : Fin 13, ∀ i, i ≠ i0 → 13 * |v i| ≤ |v i0|
    · obtain ⟨i0, hp⟩ := hpeel
      exact top_ratio_lonely_13' cite v hv i0 hp
    · push Not at hpeel
      exact hcompressed v hv hcov hbnd (fun i0 => by
        obtain ⟨i, hi, hlt⟩ := hpeel i0
        exact ⟨i, hi, by omega⟩)

-- Axiom-footprint check: pure analysis + the LRC(≤13) citation (a hypothesis), NO sorry, NO native_decide.
#print axioms top_ratio_lonely_13
#print axioms lrc14_of_compressed_lonely_13

end LRC14
end LonelyRunner
