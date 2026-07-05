/-
  TournamentH7.LRCPeelCompression — THE 24B TOP-COMPRESSION OF GAP VIOLATORS
  (kind-pasteur-2026-07-05-S6, HYP-4112).

  Attacking the spectral-gap crux itself.  THE OBSERVATION: the LRC(≤13)
  citation floor of an 11-runner family is `1/12`, and `1/12 > 2/25` — the
  citation floor one peel level down EXCEEDS the gap ceiling, with excess
  exactly `1/12 − 2/25 = 1/300`.

  Consequences, for any 12-family with NO margin-`2/25` point (in particular
  any spectral-gap violator `M ∈ (1/13, 2/25)`):

  * `peel_height_bound`: EVERY runner `v i₀` must single-handedly cover the
    `2/25`-good window that the citation gives its complementary 11-family
    (window half-width `1/(300B)` around the citation witness, by Lipschitz
    transfer).  A radius-`2/25` comb of speed `w` cannot cover an interval `J`
    with `w·|J| > 4/25` (the real-dialect interval escape, proved here on the
    S4 `distZ` anchor bound).  Hence  `|v i₀| ≤ 24·B`  where `B` bounds the
    complement.

  * `gap_compressed_24`: gap violators are **24-COMPRESSED** — every runner is
    within `24×` of some other runner.  In particular `w_max ≤ 24·w_2nd`.

  This is the first unconditional height-structure constraint on the crux:
  citation-powered, no tower recursion, no 11-level rigidity.  Combined with
  the formal profile (HYP-4102 spread gate, HYP-4105 filters, HYP-4110 big
  pair), a gap violator is now: covering in every direction, pinned `0, ±1`
  mod every `q ≤ 25`, spread (`w_max > 11.5·w_min`), big-paired
  (`some vᵢ+vⱼ ≥ 38`), and top-compressed (`every runner ≤ 24× another`).

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCGridAttainment
import TournamentH7.LRCCertCompleteness
import TournamentH7.LRC13Citation

namespace LonelyRunner
namespace PeelCompression

open TournamentH7.LRCWitness GridAttainment CertCompleteness

/-! ## The real-dialect interval escape (on the S4 distZ toolbox) -/

/-- **Interval escape over ℝ**: an interval of length `> 2r` (`r ≤ 1/2`) contains
a point at `distZ ≥ r`.  (mac-mini's `IntervalEscape` is all-ℚ; this is the ℝ
twin, per-case from the anchor bound.) -/
theorem distZ_escape (r A B : ℝ) (hr : 0 < r) (hr2 : r ≤ 1 / 2)
    (hlen : 2 * r < B - A) :
    ∃ x : ℝ, A ≤ x ∧ x ≤ B ∧ r ≤ distZ x := by
  set m : ℤ := ⌊A + 1 / 2⌋ with hm
  have hm1 : (m : ℝ) ≤ A + 1 / 2 := Int.floor_le _
  have hm2 : A + 1 / 2 < (m : ℝ) + 1 := Int.lt_floor_add_one _
  by_cases hin : (m : ℝ) + 1 / 2 ≤ B
  · -- the half-integer m + 1/2 lies in [A, B]
    refine ⟨(m : ℝ) + 1 / 2, by linarith, hin, ?_⟩
    have h1 : |(m : ℝ) + 1 / 2 - m| = 1 / 2 := by
      rw [show (m : ℝ) + 1 / 2 - m = 1 / 2 by ring]
      norm_num
    have h2 := anchor_le_distZ ((m : ℝ) + 1 / 2) m
    rw [h1] at h2
    have h3 : min (1 / 2 : ℝ) (1 - 1 / 2) = 1 / 2 := by norm_num
    rw [h3] at h2
    linarith
  · push_neg at hin
    by_cases hA : A < (m : ℝ) - r
    · -- left poke: x = A, with r < m − A ≤ 1/2
      refine ⟨A, le_refl A, by linarith, ?_⟩
      have habs : |A - (m : ℝ)| = (m : ℝ) - A := by
        rw [abs_of_nonpos (by linarith)]; ring
      have h2 := anchor_le_distZ A m
      rw [habs] at h2
      have h4 : r ≤ min ((m : ℝ) - A) (1 - ((m : ℝ) - A)) :=
        le_min (by linarith) (by linarith)
      linarith
    · -- right poke: x = B, with r < B − m < 1/2
      push_neg at hA
      refine ⟨B, by linarith, le_refl B, ?_⟩
      have habs : |B - (m : ℝ)| = B - (m : ℝ) := by
        rw [abs_of_nonneg (by linarith)]
      have h2 := anchor_le_distZ B m
      rw [habs] at h2
      have h4 : r ≤ min (B - (m : ℝ)) (1 - (B - (m : ℝ))) :=
        le_min (by linarith) (by linarith)
      linarith

/-- **The comb escape**: a radius-`r` comb of positive integer speed `w` cannot
cover an interval `[t₁, t₂]` with `w·(t₂ − t₁) > 2r`: some `t` inside keeps
`|w·t − m| ≥ r` for EVERY integer `m`. -/
theorem comb_escape (w : ℤ) (hw : 0 < w) (r t₁ t₂ : ℝ) (hr : 0 < r)
    (hr2 : r ≤ 1 / 2) (hlen : 2 * r < (w : ℝ) * (t₂ - t₁)) :
    ∃ t : ℝ, t₁ ≤ t ∧ t ≤ t₂ ∧ ∀ m : ℤ, r ≤ |(w : ℝ) * t - m| := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  obtain ⟨x, hx1, hx2, hxd⟩ := distZ_escape r ((w : ℝ) * t₁) ((w : ℝ) * t₂) hr hr2
    (by nlinarith)
  refine ⟨x / w, ?_, ?_, ?_⟩
  · rw [le_div_iff₀ hwR]; linarith
  · rw [div_le_iff₀ hwR]; linarith
  · intro m
    have heq : (w : ℝ) * (x / w) = x := by field_simp
    rw [heq]
    exact le_trans hxd (distZ_le_abs_sub x m)

/-! ## The peel -/

/-- **THE PEEL HEIGHT BOUND (the 24B top-compression).**  In a 12-family with no
margin-`2/25` point, every runner is at most `24×` a bound on its complementary
11-family: the citation hands the complement a `1/12`-margin witness; Lipschitz
keeps the complement above `2/25` on a window of half-width `1/(300B)`; the
runner must cover that window alone, and a `2/25`-comb of speed `> 24B` cannot. -/
theorem peel_height_bound (cite : LonelyRunner.LRCUpTo13) (v : Fin 12 → ℤ)
    (hv : ∀ i, v i ≠ 0) (i₀ : Fin 12) (B : ℤ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₀ → |v j| ≤ B)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25) :
    |v i₀| ≤ 24 * B := by
  by_contra hbig
  push_neg at hbig
  -- the complementary 11-family
  set C : Fin 11 → ℤ := fun j => v (i₀.succAbove j) with hCdef
  have hC0 : ∀ j, C j ≠ 0 := fun j => hv _
  have hCB : ∀ j, |C j| ≤ B := fun j => hB _ (Fin.succAbove_ne i₀ j)
  -- the citation witness: margin 1/12 for the complement
  obtain ⟨t', hL⟩ := cite 11 (by norm_num) C hC0
  have hL' : ∀ j, ∀ m : ℤ, (1 : ℝ) / 12 ≤ |(C j : ℝ) * t' - m| := by
    intro j m
    have h := hL j m
    have hcast : (1 : ℝ) / ((11 + 1 : ℕ) : ℝ) = 1 / 12 := by norm_num
    rwa [hcast] at h
  -- real bounds
  have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB0
  set δ : ℝ := 1 / (300 * (B : ℝ)) with hδdef
  have hδ0 : 0 < δ := by rw [hδdef]; positivity
  -- the escaping time for the peeled runner
  set w₀ : ℤ := |v i₀| with hw₀def
  have hw₀ : 0 < w₀ := abs_pos.mpr (hv i₀)
  have hw₀R : (24 : ℝ) * B < (w₀ : ℝ) := by exact_mod_cast hbig
  have hlen : 2 * (2 / 25 : ℝ) < (w₀ : ℝ) * ((t' + δ) - (t' - δ)) := by
    have h1 : (t' + δ) - (t' - δ) = 2 * δ := by ring
    rw [h1, hδdef]
    have h2 : (w₀ : ℝ) * (2 * (1 / (300 * (B : ℝ)))) = (w₀ : ℝ) / (150 * B) := by
      field_simp
      ring
    rw [h2, lt_div_iff₀ (by positivity : (0 : ℝ) < 150 * (B : ℝ))]
    nlinarith
  obtain ⟨t, ht1, ht2, htesc⟩ := comb_escape w₀ hw₀ (2 / 25) (t' - δ) (t' + δ)
    (by norm_num) (by norm_num) hlen
  -- sign absorption: margin for the actual (possibly negative) runner
  have hesc' : ∀ m : ℤ, 2 / 25 ≤ |(v i₀ : ℝ) * t - m| := by
    intro m
    rcases le_or_gt 0 (v i₀) with hpos | hneg
    · have hcast : (w₀ : ℝ) = (v i₀ : ℝ) := by
        rw [hw₀def]
        exact_mod_cast abs_of_nonneg hpos
      have h := htesc m
      rwa [hcast] at h
    · have hcast : (w₀ : ℝ) = -(v i₀ : ℝ) := by
        rw [hw₀def]
        exact_mod_cast abs_of_neg hneg
      have h := htesc (-m)
      rw [hcast, show -(v i₀ : ℝ) * t - ((-m : ℤ) : ℝ) = -((v i₀ : ℝ) * t - m) by
        push_cast; ring, abs_neg] at h
      exact h
  -- Lipschitz transfer: the complement stays above 2/25 on the window
  have hclose : ((B : ℤ) : ℝ) * |t - t'| ≤ 1 / 12 - 2 / 25 := by
    have habs : |t - t'| ≤ δ := by
      rw [abs_le]
      constructor <;> linarith
    have hstep : ((B : ℤ) : ℝ) * |t - t'| ≤ (B : ℝ) * δ :=
      mul_le_mul_of_nonneg_left habs hBR.le
    have hval : (B : ℝ) * δ = 1 / 300 := by
      rw [hδdef]
      field_simp
    rw [hval] at hstep
    linarith [hstep]
  have hCmargin : ∀ j, ∀ m : ℤ, 2 / 25 ≤ |(C j : ℝ) * t - m| :=
    margin_transfer C t' t (2 / 25) (1 / 12) ((B : ℤ) : ℝ)
      (fun j => by rw [← Int.cast_abs]; exact_mod_cast hCB j) hL' hclose
  -- every runner clears 2/25 at t: contradiction with hnl
  have hall : ∀ i, ∀ m : ℤ, 2 / 25 ≤ |(v i : ℝ) * t - m| := by
    intro i m
    by_cases hi : i = i₀
    · rw [hi]
      exact hesc' m
    · obtain ⟨j, hj⟩ := Fin.exists_succAbove_eq hi
      rw [← hj]
      exact hCmargin j m
  obtain ⟨i, m, hlt⟩ := hnl t
  exact absurd (hall i m) (not_le.mpr hlt)

/-- **GAP VIOLATORS ARE 24-COMPRESSED**: with no margin-`2/25` point, every
runner is within `24×` of some OTHER runner.  In particular the top ratio is
bounded: `w_max ≤ 24·w_2nd`.  (Echo of the family-level "compressed" predicate
at `13×` — the gap violator's base is itself compressed, one constant up.) -/
theorem gap_compressed_24 (cite : LonelyRunner.LRCUpTo13) (v : Fin 12 → ℤ)
    (hv : ∀ i, v i ≠ 0)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25) :
    ∀ i₀, ∃ j, j ≠ i₀ ∧ |v i₀| ≤ 24 * |v j| := by
  intro i₀
  have hne : (Finset.univ.erase i₀).Nonempty := by
    rw [← Finset.card_pos, Finset.card_erase_of_mem (Finset.mem_univ i₀)]
    simp
  obtain ⟨jm, hjmem, hjmax⟩ :=
    Finset.exists_max_image (Finset.univ.erase i₀) (fun j => |v j|) hne
  have hjne : jm ≠ i₀ := Finset.ne_of_mem_erase hjmem
  refine ⟨jm, hjne, ?_⟩
  apply peel_height_bound cite v hv i₀ (|v jm|) (abs_pos.mpr (hv jm)) ?_ hnl
  intro j hj
  exact hjmax j (Finset.mem_erase.mpr ⟨hj, Finset.mem_univ j⟩)

#print axioms distZ_escape
#print axioms comb_escape
#print axioms peel_height_bound
#print axioms gap_compressed_24

end PeelCompression
end LonelyRunner
