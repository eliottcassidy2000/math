/-
  TournamentH7.LRCTeethR — THE RADIUS-PARAMETRIC WINDOW STACK (klein-2026-07-05-S136,
  HYP-4107; the Fin-12/ρ=1/13 transcription requested by mac-mini-S52 and kps-S3).

  Everything the window/ledger machinery needs, with the danger radius ρ a PARAMETER
  (0 < ρ < 1/2) instead of the hardcoded 1/14: teeth, mass bounds (density AND sparse),
  the Hunter block step, and the ι-generic multi-top window whose conclusion is a
  MARGIN statement — `Lonely 14` is the ρ = 1/14, ι = Fin 13 instance, and the n = 13
  rigidity level (ρ = 1/13, β = 1/12, Fin 12) is another: mac-mini's fee table
  `T_l = 156(13−l)/(13−2l)` lives here, with the `2l < 13` wall as the ρ-generic
  `(1 − 2lρ) > 0` remainder.

  Sources mirrored: LRCBlockSix (tooth/teeth/mass, 1/14 → ρ), LRCRealRegions
  (hunter_block_step; the region calculus and hunter_ledger are radius-free and
  consumed as-is), LRCMultiKillerWindow (the sparse bound and the fee composition).
-/
import TournamentH7.LRCRealRegions
import TournamentH7.LRCMultiKillerWindow

namespace LonelyRunner
namespace TeethR

open LRC14 RealRegion

/-- One radius-ρ tooth of runner `w` at integer `m`. -/
noncomputable def toothR (ρ : ℝ) (w : ℤ) (m : ℤ) : ℝ × ℝ :=
  (((m : ℝ) - ρ) / w, ((m : ℝ) + ρ) / w)

/-- The radius-ρ teeth of runner `w` meeting the window `[a, b]`. -/
noncomputable def teethR (ρ : ℝ) (w : ℤ) (a b : ℝ) : List (ℝ × ℝ) :=
  (List.range (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat).map fun (i : ℕ) =>
    toothR ρ w (⌈(w : ℝ) * a⌉ - 1 + (i : ℤ))

theorem teethR_live {ρ : ℝ} (hρ : 0 ≤ ρ) {w : ℤ} (hw : 0 < w) (a b : ℝ) :
    ∀ r ∈ teethR ρ w a b, r.1 ≤ r.2 := by
  intro r hr
  rw [teethR, List.mem_map] at hr
  obtain ⟨i, -, rfl⟩ := hr
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  unfold toothR
  simp only
  gcongr
  linarith

theorem teethR_sortedSep {ρ : ℝ} (hρ : ρ ≤ 1/2) {w : ℤ} (hw : 0 < w) (a b : ℝ) :
    SortedSep (teethR ρ w a b) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  unfold teethR
  apply sortedSep_map_range
  intro i j hij _
  unfold toothR
  simp only
  have hij' : (i : ℝ) + 1 ≤ (j : ℝ) := by exact_mod_cast hij
  gcongr
  push_cast
  linarith

/-- Each clipped radius-ρ tooth carries at most `2ρ/w`. -/
theorem clipLen_toothR_le {ρ : ℝ} (hρ : 0 ≤ ρ) {w : ℤ} (hw : 0 < w) (m : ℤ) (a b : ℝ) :
    clipLen (toothR ρ w m) a b ≤ 2 * ρ / (w : ℝ) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  unfold clipLen toothR
  simp only
  rcases le_total (min (((m : ℝ) + ρ) / w) b - max (((m : ℝ) - ρ) / w) a) 0 with h | h
  · rw [max_eq_left h]
    positivity
  · rw [max_eq_right h]
    have h1 : min (((m : ℝ) + ρ) / w) b ≤ ((m : ℝ) + ρ) / w := min_le_left _ _
    have h2 : ((m : ℝ) - ρ) / w ≤ max (((m : ℝ) - ρ) / w) a := le_max_left _ _
    have h3 : ((m : ℝ) + ρ) / w - ((m : ℝ) - ρ) / w = 2 * ρ / (w : ℝ) := by
      field_simp
      ring
    linarith

/-- **Density mass bound**: the clipped radius-ρ teeth of `w` carry at most
`2ρ(b−a) + 6ρ/w` over the window. -/
theorem teethR_mass {ρ : ℝ} (hρ : 0 ≤ ρ) {w : ℤ} (hw : 0 < w) (a b : ℝ) (hab : a ≤ b) :
    ((teethR ρ w a b).map fun p => clipLen p a b).sum
      ≤ 2 * ρ * (b - a) + 6 * ρ / (w : ℝ) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  have hper : ∀ x ∈ (teethR ρ w a b).map fun p => clipLen p a b,
      x ≤ 2 * ρ / (w : ℝ) := by
    intro x hx
    rw [List.mem_map] at hx
    obtain ⟨p, hp, rfl⟩ := hx
    rw [teethR, List.mem_map] at hp
    obtain ⟨i, -, rfl⟩ := hp
    exact clipLen_toothR_le hρ hw _ a b
  have hlen : ((teethR ρ w a b).map fun p => clipLen p a b).length
      = (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat := by
    rw [List.length_map, teethR, List.length_map, List.length_range]
  have hsum : ((teethR ρ w a b).map fun p => clipLen p a b).sum
      ≤ (((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℝ)) * (2 * ρ / (w : ℝ)) := by
    rw [← hlen]
    exact List.sum_le_card_nsmul _ _ hper |>.trans (by rw [nsmul_eq_mul])
  have hcount : (((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℝ))
      ≤ (w : ℝ) * (b - a) + 3 := by
    have h1 : ((⌊(w : ℝ) * b⌋ : ℤ) : ℝ) ≤ (w : ℝ) * b := Int.floor_le _
    have h2 : (w : ℝ) * a ≤ ((⌈(w : ℝ) * a⌉ : ℤ) : ℝ) := Int.le_ceil _
    rcases le_or_gt (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1) 0 with h | h
    · rw [Int.toNat_of_nonpos h]
      norm_num
      nlinarith [mul_nonneg hwR.le (sub_nonneg.mpr hab)]
    · have hchain : (((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat : ℕ) : ℝ)
          = ((⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1 : ℤ) : ℝ) := by
        rw [← Int.cast_natCast, Int.toNat_of_nonneg (le_of_lt h)]
      rw [hchain]
      push_cast
      have hexp : (w : ℝ) * (b - a) = (w : ℝ) * b - (w : ℝ) * a := by ring
      linarith
  calc ((teethR ρ w a b).map fun p => clipLen p a b).sum
      ≤ _ * (2 * ρ / (w : ℝ)) := hsum
    _ ≤ ((w : ℝ) * (b - a) + 3) * (2 * ρ / (w : ℝ)) := by
        apply mul_le_mul_of_nonneg_right hcount
        positivity
    _ = 2 * ρ * (b - a) + 6 * ρ / (w : ℝ) := by
        field_simp
        ring

end TeethR
end LonelyRunner
