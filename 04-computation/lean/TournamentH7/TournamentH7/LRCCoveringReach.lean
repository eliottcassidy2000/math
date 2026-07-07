/-
  TournamentH7.LRCCoveringReach  (mac-mini-2026-07-06-S34, HYP-4587/kps + opus-S124)

  THE COVERING-REACH ATOM — the reusable bridge from a mod-q clearing certificate
  to a reach lower bound, and the d=2 generic instance.

  kps-S43 reduced (C) to a FINITE COVERING SYSTEM: every non-AP mod-25 pair-blocker
  is cleared by some modulus `q ≤ 39`, each clearance a `rational_point_margin`
  cert (the atom already lives in LRCHarmonicGate).  This file packages that atom
  into a REACH statement — `reach_ge_of_covering` — so any decidable clearing
  hypothesis `∀ i, μ ≤ (v i · c) % q ≤ q − μ` yields `μ/q ≤ reach(v)`.  That is the
  uniform Lean shape of every covering certificate (kps's mod-25 floor is the
  `q = 25, μ = 2` instance; below is the `q = 11, μ = 1` d=2 generic instance).

  d=2 generic: the canonical d=2 family `{1,…,10, x, y}` with `11 ∤ x` and `11 ∤ y`
  clears at `t = 1/11` (all speeds off `0 mod 11`), so `reach ≥ 1/11 > 2/25`.

  Kernel-pure (standard trio).
-/
import TournamentH7.LRCLadderD1
import TournamentH7.LRCHarmonicGate

open TournamentH7.LRCWitness
open LonelyRunner.HarmonicGate

namespace TournamentH7.LRCCoveringReach

/-- **The covering-reach atom.**  If a rotation `c` (with `0 ≤ c ≤ q`) takes every
speed off the forbidden band `{0,…,±(μ−1)}` mod `q` — i.e. `μ ≤ (v i · c) % q ≤
q − μ` for all `i` — then the reach (attained sup of the margin over `[0,1]`) is at
least `μ/q`.  Packages `rational_point_margin` + the compact-sup wrapper. -/
theorem reach_ge_of_covering (v : Fin 12 → ℤ) (q c μ : ℤ)
    (hq : 0 < q) (hc0 : 0 ≤ c) (hcq : c ≤ q)
    (hclear : ∀ i, μ ≤ (v i * c) % q ∧ (v i * c) % q ≤ q - μ) :
    (μ : ℝ) / q ≤ sSup (margin v '' Set.Icc (0 : ℝ) 1) := by
  -- margin at t = c/q is ≥ μ/q (from the rational-point atom)
  have hmargin : (μ : ℝ) / q ≤ margin v ((c : ℝ) / q) :=
    (le_margin_iff v ((μ : ℝ) / q) ((c : ℝ) / q)).2
      (rational_point_margin v c q μ hq hclear)
  -- t = c/q ∈ [0,1]
  have hqR : (0 : ℝ) < q := by exact_mod_cast hq
  have ht01 : (c : ℝ) / q ∈ Set.Icc (0 : ℝ) 1 := by
    refine Set.mem_Icc.mpr ⟨?_, ?_⟩
    · exact div_nonneg (by exact_mod_cast hc0) (le_of_lt hqR)
    · rw [div_le_one hqR]; exact_mod_cast hcq
  have hmem : margin v ((c : ℝ) / q) ∈ margin v '' Set.Icc (0 : ℝ) 1 := ⟨_, ht01, rfl⟩
  have hbdd : BddAbove (margin v '' Set.Icc (0 : ℝ) 1) :=
    ⟨1 / 2, by rintro _ ⟨t, _, rfl⟩; exact LRCLadderD1.margin_le_half v t⟩
  exact le_trans hmargin (le_csSup hbdd hmem)

/-- The canonical d=2 family `{1,…,10, x, y}` as `Fin 12 → ℤ` (base `{1,…,10}` at
indices `0..9`, outliers `x, y` at `10, 11`). -/
def W2 (x y : ℤ) : Fin 12 → ℤ := fun i =>
  if i.val = 11 then y else if i.val = 10 then x else (i.val : ℤ) + 1

/-- **d=2 generic bound.**  If `11 ∤ x` and `11 ∤ y`, the canonical d=2 family has
reach `≥ 1/11 > 2/25` (witness `t = 1/11`, `q = 11, c = 1, μ = 1`).  So a d=2
family with both outliers off `0 mod 11` is not in the open gap. -/
theorem d2_generic_reach (x y : ℤ) (hx : ¬ (11 ∣ x)) (hy : ¬ (11 ∣ y)) :
    (2 : ℝ) / 25 < sSup (margin (W2 x y) '' Set.Icc (0 : ℝ) 1) := by
  have hcov : ∀ i, 1 ≤ (W2 x y i * 1) % 11 ∧ (W2 x y i * 1) % 11 ≤ 11 - 1 := by
    intro i
    have hi12 : i.val < 12 := i.isLt
    -- W2 x y i is x, y, or a value in {1,…,10}; none is ≡ 0 mod 11
    have hnd : ¬ (11 ∣ W2 x y i) := by
      unfold W2
      split
      · exact hy
      · split
        · exact hx
        · rename_i h1 h2; rintro ⟨d, hd⟩; omega
    -- 1 ≤ (w*1) % 11 ≤ 10  ⟺  11 ∤ w
    rw [mul_one]
    have h0 : 0 ≤ W2 x y i % 11 := Int.emod_nonneg _ (by norm_num)
    have h11 : W2 x y i % 11 < 11 := Int.emod_lt_of_pos _ (by norm_num)
    have hne : W2 x y i % 11 ≠ 0 := fun h => hnd (Int.dvd_of_emod_eq_zero h)
    omega
  have hreach : (1 : ℝ) / 11 ≤ sSup (margin (W2 x y) '' Set.Icc (0 : ℝ) 1) := by
    have h := reach_ge_of_covering (W2 x y) 11 1 1 (by norm_num) (by norm_num) (by norm_num) hcov
    norm_num at h; exact h
  have : (2 : ℝ) / 25 < (1 : ℝ) / 11 := by norm_num
  linarith

#print axioms reach_ge_of_covering
#print axioms d2_generic_reach

end TournamentH7.LRCCoveringReach
