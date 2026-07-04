/-
  TournamentH7.LRCFolding — THE FOLDING IDENTITY'S ENGINE, FORMALIZED (THM-615 Lemma 1)
  (klein-2026-07-04-S123, HYP-4075).

  opus-S64's THM-615 reduces the m=2, f=2 CONFINEMENT (the covering-min rigidity's hard core) to
  `M(2U ∪ {w₁,w₂}) ≥ 1/12`, via the FOLDING IDENTITY
     M(2U ∪ {w₁,w₂}) = max_t min( g_E(t), Ψ(t) ),   Ψ = max(min a b, ½ − max a b),  a=‖w₁t‖,b=‖w₂t‖,
  which pairs `t` with `t+½` using: for ODD w, `‖w(t+½)‖ = ½ − ‖wt‖`; for EVEN 2u, `‖2u(t+½)‖ =
  ‖2ut‖`; and the lattice identity `max(min g X)(min g Y) = min g (max X Y)`.

  This file formalizes that ENGINE, sorry-free, with `reach x = dist(x, ℤ)`:
   * `reach_add_half`      : reach (x + ½) = ½ − reach x
   * `reach_odd_add_half`  : w odd ⟹ reach (w*(t+½)) = ½ − reach (w*t)        -- the tightener flip
   * `reach_even_add_half` : reach ((2u)*(t+½)) = reach ((2u)*t)              -- even-part (+½)-periodicity
   * `min_min_fold`        : max (min g X) (min g Y) = min g (max X Y)        -- the lattice fold
   * `folded_pair`         : the two odd tighteners' folded value = Ψ at the pair {t, t+½}
   * `psi_ge_of_not_extremity` : Ψ ≥ b unless "extremity" (one reach < b and the other > ½−b)

  SCOPE (HONEST): this is the folding identity's algebraic engine (elementary, THM-615 Lemma 1),
  machine-checked.  It does NOT prove `M(2U∪2odd) ≥ 1/12` (the confinement core = argmax arithmetic,
  opus's active open line; the general `U` case is unproved and LRC(14)-equivalent).  Infrastructure
  for the eventual formal confinement proof, not the closure.
-/
import Mathlib

namespace LonelyRunner.Folding

/-- `reach x = distance from x to the nearest integer = min(fract x, 1 − fract x) ∈ [0, ½]`. -/
noncomputable def reach (x : ℝ) : ℝ := min (Int.fract x) (1 - Int.fract x)

theorem reach_nonneg (x : ℝ) : 0 ≤ reach x := by
  unfold reach
  have h0 := Int.fract_nonneg x
  have h1 := Int.fract_lt_one x
  exact le_min h0 (by linarith)

theorem reach_le_half (x : ℝ) : reach x ≤ 1 / 2 := by
  unfold reach
  rcases le_total (Int.fract x) (1 / 2) with h | h
  · exact le_trans (min_le_left _ _) h
  · exact le_trans (min_le_right _ _) (by linarith)

/-- `reach` is `ℤ`-periodic. -/
theorem reach_add_int (x : ℝ) (n : ℤ) : reach (x + n) = reach x := by
  unfold reach
  have : Int.fract (x + (n : ℝ)) = Int.fract x := Int.fract_eq_fract.mpr ⟨n, by push_cast; ring⟩
  rw [this]

/-- **The half-shift.**  `reach (x + ½) = ½ − reach x`.  (Case split on `fract x ≶ ½`.) -/
theorem reach_add_half (x : ℝ) : reach (x + 1 / 2) = 1 / 2 - reach x := by
  have h0 := Int.fract_nonneg x
  have h1 := Int.fract_lt_one x
  have hfl : x - Int.fract x = (⌊x⌋ : ℝ) := by have := Int.floor_add_fract x; linarith
  unfold reach
  rcases lt_or_ge (Int.fract x) (1 / 2) with hlt | hge
  · have hf : Int.fract (x + 1 / 2) = Int.fract x + 1 / 2 := by
      have heq : Int.fract (x + 1 / 2) = Int.fract (Int.fract x + 1 / 2) :=
        Int.fract_eq_fract.mpr ⟨⌊x⌋, by linarith [hfl]⟩
      rw [heq, Int.fract_eq_self.mpr ⟨by linarith, by linarith⟩]
    rw [hf, min_eq_right (by linarith), min_eq_left (by linarith)]; ring
  · have hf : Int.fract (x + 1 / 2) = Int.fract x - 1 / 2 := by
      have heq : Int.fract (x + 1 / 2) = Int.fract (Int.fract x - 1 / 2) :=
        Int.fract_eq_fract.mpr ⟨⌊x⌋ + 1, by push_cast; linarith [hfl]⟩
      rw [heq, Int.fract_eq_self.mpr ⟨by linarith, by linarith⟩]
    rw [hf, min_eq_left (by linarith), min_eq_right (by linarith)]; ring

/-- **Tightener flip.**  For ODD `w`, `reach (w·(t+½)) = ½ − reach (w·t)`. -/
theorem reach_odd_add_half {w : ℤ} (hw : Odd w) (t : ℝ) :
    reach ((w : ℝ) * (t + 1 / 2)) = 1 / 2 - reach ((w : ℝ) * t) := by
  obtain ⟨k, hk⟩ := hw
  have hwR : (w : ℝ) = 2 * (k : ℝ) + 1 := by exact_mod_cast hk
  have hrw : (w : ℝ) * (t + 1 / 2) = ((w : ℝ) * t + 1 / 2) + (k : ℝ) := by rw [hwR]; ring
  rw [hrw, reach_add_int, reach_add_half]

/-- **Even-part (+½)-periodicity.**  `reach ((2u)·(t+½)) = reach ((2u)·t)`. -/
theorem reach_even_add_half (u : ℤ) (t : ℝ) :
    reach (((2 * u : ℤ) : ℝ) * (t + 1 / 2)) = reach (((2 * u : ℤ) : ℝ) * t) := by
  have hrw : ((2 * u : ℤ) : ℝ) * (t + 1 / 2) = ((2 * u : ℤ) : ℝ) * t + (u : ℝ) := by
    push_cast; ring
  rw [hrw, reach_add_int]

/-- **The lattice fold.**  `max (min g X) (min g Y) = min g (max X Y)`. -/
theorem min_min_fold (g X Y : ℝ) : max (min g X) (min g Y) = min g (max X Y) :=
  (min_max_distrib_left g X Y).symm

/-- **The single-tightener fold.**  Over the pair `{t, t+½}`, the best of `min(g, reach(w·))` for
  odd `w` (with `g` a `(+½)`-invariant value `gv`) is `min(gv, max(reach(w t), ½ − reach(w t)))`. -/
theorem folded_pair {w : ℤ} (hw : Odd w) (t gv : ℝ) :
    max (min gv (reach ((w : ℝ) * t))) (min gv (reach ((w : ℝ) * (t + 1 / 2))))
      = min gv (max (reach ((w : ℝ) * t)) (1 / 2 - reach ((w : ℝ) * t))) := by
  rw [reach_odd_add_half hw, min_min_fold]

/-- **`Ψ` and extremity.**  `Ψ(a,b) := max(min a b, ½ − max a b)`.  For `0 ≤ b ≤ 1/4`, `Ψ ≥ b`
  UNLESS "extremity": one of `a₁,a₂ < b` and the other `> ½ − b`.  (Here stated as the contrapositive
  positive form used in THM-615: not-extremity ⟹ `Ψ ≥ b`.) -/
theorem psi_ge_of_not_extremity (a₁ a₂ b : ℝ)
    (hnot : ¬ ((a₁ < b ∧ a₂ > 1 / 2 - b) ∨ (a₂ < b ∧ a₁ > 1 / 2 - b))) :
    b ≤ max (min a₁ a₂) (1 / 2 - max a₁ a₂) := by
  by_contra hlt
  push_neg at hlt
  have h1 : min a₁ a₂ < b := lt_of_le_of_lt (le_max_left _ _) hlt
  have h2 : max a₁ a₂ > 1 / 2 - b := by
    have := lt_of_le_of_lt (le_max_right _ _) hlt
    linarith
  apply hnot
  rcases le_total a₁ a₂ with h | h
  · rw [min_eq_left h] at h1; rw [max_eq_right h] at h2
    exact Or.inl ⟨h1, h2⟩
  · rw [min_eq_right h] at h1; rw [max_eq_left h] at h2
    exact Or.inr ⟨h1, h2⟩

#print axioms reach_add_half
#print axioms reach_odd_add_half
#print axioms reach_even_add_half
#print axioms min_min_fold
#print axioms folded_pair

end LonelyRunner.Folding
