/- LRCPairOverlapFloor.lean -- opus-2026-07-17-S350 (HYP-7330).
   THE PAIR-OVERLAP LOWER BOUND BY CONTAINMENT — the surviving half of the
   OverlapMeasureBridge, obtained WITHOUT the exact sawtooth identity.

   Two containments, both elementary:

   (1) `pair_overlap_contains` : for `0 < a ≤ b`, the single arc around 0,
       `Ioo (-(lam/b)) (lam/b)`, lies in BOTH combs — its `b`-radius is the
       smaller one, and `a ≤ b` makes it fit inside the `a`-arc too.  Hence
       (windowed, since the periodic intersection has infinite measure on ℝ)
       `2*lam/b ≤ volume (badArcs a lam ∩ badArcs b lam ∩ Ioo (-(1/2)) (1/2))`.
       POSITIVE — which is all the S349 existence capstone needs.

   (2) `pair_overlap_contains_gcd` : THE GCD STRENGTHENING.  `1/g` with
       `g = gcd a b` is a COMMON PERIOD of both combs (g ∣ a and g ∣ b), so
       the same argument runs around EVERY point `j/g`, not just `0`:
       `badArcs g (lam*g/b) ⊆ badArcs a lam ∩ badArcs b lam`.
       Over a period this is a factor-`g` improvement (`2*lam*g/b`), and it
       degrades gracefully to (1) at `g = 1`.

   The sharp value (`≈ 1/49 = (2λ)²`, the independence heuristic) still needs
   the sawtooth identity THM-965/856 — but only for the SHARP floor that
   nesting wants, never for existence.  Kernel-pure: no sorry, no
   native_decide. -/
import Mathlib
import TournamentH7.FragmentationLemma

open MeasureTheory

namespace LonelyRunner.LRC14.Hunter

/-- Membership in a comb via its `k = 0` arc. -/
theorem mem_badArcs_of_abs_lt {w : ℕ} {lam t : ℝ} (h : |t| < lam / w) :
    t ∈ LRC14.badArcs w lam := by
  rw [LRC14.badArcs, Set.mem_iUnion]
  refine ⟨0, ?_⟩
  rw [Set.mem_Ioo]
  rw [abs_lt] at h
  constructor <;> · simp only [Int.cast_zero, zero_div]; linarith [h.1, h.2]

/-- **(1) THE CONTAINMENT**: the `0`-arc of the FASTER comb sits inside both. -/
theorem pair_overlap_contains {a b : ℕ} {lam : ℝ} (ha : 0 < a) (hab : a ≤ b)
    (hlam : 0 < lam) :
    Set.Ioo (-(lam / b)) (lam / b)
      ⊆ LRC14.badArcs a lam ∩ LRC14.badArcs b lam := by
  have haR : (0 : ℝ) < a := by exact_mod_cast ha
  have hbR : (0 : ℝ) < b := lt_of_lt_of_le haR (by exact_mod_cast hab)
  have hle : lam / b ≤ lam / a := by
    apply div_le_div_of_nonneg_left (le_of_lt hlam) haR
    exact_mod_cast hab
  intro t ht
  rw [Set.mem_Ioo] at ht
  have habs : |t| < lam / b := abs_lt.mpr ⟨ht.1, ht.2⟩
  exact ⟨mem_badArcs_of_abs_lt (lt_of_lt_of_le habs hle),
    mem_badArcs_of_abs_lt habs⟩

/-- **THE WINDOWED VOLUME FLOOR** — positive, which is what the existence
capstone consumes.  (Unwindowed the intersection is periodic, hence of
infinite measure on `ℝ`, so the bound must be stated on a fundamental
domain.) -/
theorem pair_overlap_volume_ge {a b : ℕ} {lam : ℝ} (ha : 0 < a) (hab : a ≤ b)
    (hlam : 0 < lam) (hsmall : lam / b ≤ 1 / 2) :
    ENNReal.ofReal (2 * (lam / b))
      ≤ volume (LRC14.badArcs a lam ∩ LRC14.badArcs b lam
          ∩ Set.Ioo (-(1 / 2) : ℝ) (1 / 2)) := by
  have hsub : Set.Ioo (-(lam / b)) (lam / b)
      ⊆ LRC14.badArcs a lam ∩ LRC14.badArcs b lam
        ∩ Set.Ioo (-(1 / 2) : ℝ) (1 / 2) := by
    intro t ht
    refine ⟨pair_overlap_contains ha hab hlam ht, ?_⟩
    rw [Set.mem_Ioo] at ht ⊢
    exact ⟨by linarith [ht.1], by linarith [ht.2]⟩
  calc ENNReal.ofReal (2 * (lam / b))
      = volume (Set.Ioo (-(lam / b)) (lam / b)) := by
        rw [Real.volume_Ioo]; congr 1; ring
    _ ≤ _ := measure_mono hsub

/-- **(2) THE GCD STRENGTHENING**: `1/g` is a common period of both combs, so
the `0`-arc argument runs around every `j/g`.  Stated as the comb inclusion
`badArcs g (lam*g/b) ⊆ badArcs a lam ∩ badArcs b lam`. -/
theorem pair_overlap_contains_gcd {a b : ℕ} {lam : ℝ} (ha : 0 < a) (hab : a ≤ b)
    (hlam : 0 < lam) :
    LRC14.badArcs (Nat.gcd a b) (lam * Nat.gcd a b / b)
      ⊆ LRC14.badArcs a lam ∩ LRC14.badArcs b lam := by
  set g := Nat.gcd a b with hg
  have hgpos : 0 < g := Nat.gcd_pos_of_pos_left b ha
  have hgR : (0 : ℝ) < g := by exact_mod_cast hgpos
  have haR : (0 : ℝ) < a := by exact_mod_cast ha
  have hbR : (0 : ℝ) < b := lt_of_lt_of_le haR (by exact_mod_cast hab)
  have hga : g ∣ a := Nat.gcd_dvd_left a b
  have hgb : g ∣ b := Nat.gcd_dvd_right a b
  obtain ⟨a', ha'⟩ := hga
  obtain ⟨b', hb'⟩ := hgb
  have ha'pos : 0 < a' := by
    rcases Nat.eq_zero_or_pos a' with h | h
    · rw [h, Nat.mul_zero] at ha'; omega
    · exact h
  have hb'pos : 0 < b' := by
    rcases Nat.eq_zero_or_pos b' with h | h
    · rw [h, Nat.mul_zero] at hb'; omega
    · exact h
  have ha'R : (0 : ℝ) < a' := by exact_mod_cast ha'pos
  have hb'R : (0 : ℝ) < b' := by exact_mod_cast hb'pos
  intro t ht
  rw [LRC14.badArcs, Set.mem_iUnion] at ht
  obtain ⟨j, hj⟩ := ht
  rw [Set.mem_Ioo] at hj
  have hradius : lam * (g : ℝ) / (b : ℝ) / (g : ℝ) = lam / (b : ℝ) := by
    field_simp
  rw [hradius] at hj
  -- the common near-point is j/g; index it in each comb
  refine ⟨?_, ?_⟩
  · rw [LRC14.badArcs, Set.mem_iUnion]
    refine ⟨(a' : ℤ) * j, ?_⟩
    rw [Set.mem_Ioo]
    have hkey : (((a' : ℤ) * j : ℤ) : ℝ) / (a : ℝ) = (j : ℝ) / (g : ℝ) := by
      rw [ha']
      push_cast
      field_simp
    have hle : lam / (b : ℝ) ≤ lam / (a : ℝ) := by
      apply div_le_div_of_nonneg_left (le_of_lt hlam) haR
      exact_mod_cast hab
    rw [hkey]
    constructor <;> linarith [hj.1, hj.2]
  · rw [LRC14.badArcs, Set.mem_iUnion]
    refine ⟨(b' : ℤ) * j, ?_⟩
    rw [Set.mem_Ioo]
    have hkey : (((b' : ℤ) * j : ℤ) : ℝ) / (b : ℝ) = (j : ℝ) / (g : ℝ) := by
      rw [hb']
      push_cast
      field_simp
    rw [hkey]
    constructor <;> linarith [hj.1, hj.2]

/-! ## Axiom audit -/
#print axioms mem_badArcs_of_abs_lt
#print axioms pair_overlap_contains
#print axioms pair_overlap_volume_ge
#print axioms pair_overlap_contains_gcd

end LonelyRunner.LRC14.Hunter
