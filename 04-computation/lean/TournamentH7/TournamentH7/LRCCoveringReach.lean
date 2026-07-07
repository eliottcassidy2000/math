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

/-! ### The translate-block bound (opus-S127 branch 3, uniform-k L-lift escape) -/

/-- The consecutive 12-block `{m, m+1, …, m+11}` as `Fin 12 → ℤ`. -/
def Blk (m : ℤ) : Fin 12 → ℤ := fun i => m + (i.val : ℤ)

/-- **Translate-block bound.**  The consecutive block `{m,…,m+11}` for `m ≥ 2` has
reach `≥ m/(2m+11) ≥ 2/25` — witness `t = 1/(2m+11)`, clearance `m` (a covering cert
at `q = 2m+11`, `c = 1`, `μ = m`).  Only `m = 1` (the AP `{1,…,12}`) is tight
(`1/13`).  So every non-AP translate of the AP — in particular every uniform L-lift
`{1+Lk,…,12+Lk}` with `k ≥ 1` (`m = 1+Lk ≥ 2`) — is loose.  This closes opus-S127's
branch-3 uniform-k escape of the finite-covering crux. -/
theorem translate_block_reach (m : ℤ) (hm : 2 ≤ m) :
    (2 : ℝ) / 25 ≤ sSup (margin (Blk m) '' Set.Icc (0 : ℝ) 1) := by
  have hq : (0 : ℤ) < 2 * m + 11 := by linarith
  have hcov : ∀ i, m ≤ (Blk m i * 1) % (2 * m + 11) ∧
      (Blk m i * 1) % (2 * m + 11) ≤ (2 * m + 11) - m := by
    intro i
    have hi12 : i.val < 12 := i.isLt
    have hi11 : (i.val : ℤ) ≤ 11 := by exact_mod_cast Nat.lt_succ_iff.mp hi12
    have hi0 : (0 : ℤ) ≤ (i.val : ℤ) := Int.natCast_nonneg _
    have hval : Blk m i * 1 = m + (i.val : ℤ) := by simp [Blk]
    rw [hval, Int.emod_eq_of_lt (by linarith) (by linarith)]
    exact ⟨by linarith, by linarith⟩
  have hreach := reach_ge_of_covering (Blk m) (2 * m + 11) 1 m hq (by norm_num)
    (by linarith) hcov
  have hmR : (2 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm
  have hqR : (0 : ℝ) < 2 * (m : ℝ) + 11 := by linarith
  have hge : (2 : ℝ) / 25 ≤ (m : ℝ) / (2 * (m : ℝ) + 11) := by
    rw [le_div_iff₀ hqR]; nlinarith
  push_cast at hreach
  linarith

/-! ### The hardest r=2 double-lift cert — the true covering modulus Q0 = 37

Challenging kps-S47's `Q0 = 25` (and kps-S44's `≤ 14`): the r=2 double-13-lift
shapes need moduli up to **37**, confirming klein-S144's `≤ 38`.  The global worst
family (over all 66 shapes and lifts `a,b ∈ [0,80)`) is shape `(10,12)` at `a=2,
b=26`: `{1,…,9, 11, 36, 350}` — a mod-25 blocker that clears at NO `q ≤ 36`, only at
`q = 37` (`c = 3`, `μ = 3`). This is a genuine covering cert at the true Q0. -/

/-- The hardest r=2 family: `{1,…,9, 11, 36, 350}` (shape (10,12), a=2, b=26). -/
def hardR2 : Fin 12 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 36, 350]

/-- **The hardest r=2 family is loose** at the true covering modulus `Q0 = 37`:
reach `≥ 3/37 > 2/25` (cert `q = 37, c = 3, μ = 3`; residues `3,6,9,12,14,15,18,21,
24,27,33,34` all in `[3,34]`).  A machine-checked witness that the covering-reach
atom scales to the real `Q0 = 37`, not kps-S47's `25`. -/
theorem hardR2_reach : (2 : ℝ) / 25 < sSup (margin hardR2 '' Set.Icc (0 : ℝ) 1) := by
  have hcov : ∀ i, (3 : ℤ) ≤ (hardR2 i * 3) % 37 ∧ (hardR2 i * 3) % 37 ≤ 37 - 3 := by
    decide
  have h := reach_ge_of_covering hardR2 37 3 3 (by norm_num) (by norm_num) (by norm_num) hcov
  norm_num at h
  have hnum : (2 : ℝ) / 25 < (3 : ℝ) / 37 := by norm_num
  linarith

#print axioms reach_ge_of_covering
#print axioms d2_generic_reach
#print axioms translate_block_reach
#print axioms hardR2_reach

end TournamentH7.LRCCoveringReach
