/-
  TournamentH7.LRCSingleKillerLadder  (kind-pasteur-2026-07-11-S127 cont.61)

  THE SINGLE-KILLER LADDER, machine-checked — the whole family, not just the extremal.

  cont.60 (LRCDeepWellReach) certified `reach({1,…,12, 182}) ≥ 14/183` (the deep well,
  `c = 1`). The single-killer balance (opus-S253, kps cont.59) says the FULL ladder

        S_c = {1, …, 12, 182·c}   (c ≥ 1,  killer = c·lcm(13,14))

  has `reach(S_c) ≥ 14c/(182c+1) = [1/13]·(182c)/(182c+1)`, increasing in `c`, TIGHT at
  `c = 1` (`14/183`) and `→ 1/13` as `c → ∞`. Every rung is a primitive divisor-complete
  single-killer covering family, so this closes the ENTIRE single-killer class in the
  reach API: `reach(S_c) ≥ 14/183` for all `c ≥ 1`, equality only at the deep well.

  The witness is the covering certificate at `q = 182c+1`, rotation `14c`, band `[14c,
  168c+1] = [μ, q−μ]` (`μ = 14c`): AP runner `j` lands at `14jc ∈ [14c, 168c]`, and the
  killer `182c` lands at `168c+1 = q−14c` (the upper edge), via `182c ≡ −1 (mod 182c+1)`
  — the parametrized twin of klein-S119's `n² ≡ n−1 (mod Φ₆)` at `n = 14`.

  Kernel-pure (standard trio: propext, Classical.choice, Quot.sound).
-/
import TournamentH7.LRCDeepWellReach

open TournamentH7.LRCWitness TournamentH7.LRCDeepWellReach

namespace TournamentH7.LRCSingleKillerLadder

/-- The single-killer ladder family `S_c = {1,…,12, 182·c}` as `Fin 13 → ℤ`. -/
def singleKiller (c : ℤ) : Fin 13 → ℤ := fun i => if i.val = 12 then 182 * c else (i.val : ℤ) + 1

/-- **The ladder clearing certificate.**  For `c ≥ 1`, rotation `14c` takes every runner of
`S_c` into the band `[14c, (182c+1) − 14c]` mod `q = 182c+1`: AP runners `j ↦ 14jc ∈
[14c,168c]`, killer `182c ↦ 168c+1` (upper edge, via `182c ≡ −1 mod q`). -/
theorem singleKiller_cover (c : ℤ) (hc : 1 ≤ c) :
    ∀ i, (14 * c) ≤ (singleKiller c i * (14 * c)) % (182 * c + 1) ∧
         (singleKiller c i * (14 * c)) % (182 * c + 1) ≤ (182 * c + 1) - 14 * c := by
  intro i
  have hcpos : (0 : ℤ) ≤ 14 * c := by linarith
  by_cases h : i.val = 12
  · -- KILLER runner `182c`: residue `168c+1 = q−14c` (via `182c ≡ −1 mod q`).
    have hval : singleKiller c i = 182 * c := by simp [singleKiller, h]
    rw [hval]
    have hcong : (182 * c * (14 * c)) % (182 * c + 1) = (168 * c + 1) % (182 * c + 1) := by
      have hdvd : (182 * c + 1) ∣ (168 * c + 1) - 182 * c * (14 * c) := ⟨1 - 14 * c, by ring⟩
      exact Int.modEq_iff_dvd.mpr hdvd
    rw [hcong, Int.emod_eq_of_lt (by linarith) (by nlinarith)]
    exact ⟨by linarith, by linarith⟩
  · -- AP runner `j = i.val+1`, `1 ≤ j ≤ 12`: residue `14jc ∈ [14c, 168c]`.
    have hval : singleKiller c i = (i.val : ℤ) + 1 := by simp [singleKiller, h]
    rw [hval]
    have hi0 : (0 : ℤ) ≤ (i.val : ℤ) := Int.natCast_nonneg _
    have hi11 : (i.val : ℤ) ≤ 11 := by
      have hn : i.val ≤ 11 := by have := i.isLt; omega
      exact_mod_cast hn
    have ha0 : (0 : ℤ) ≤ ((i.val : ℤ) + 1) * (14 * c) := mul_nonneg (by linarith) hcpos
    have hab : ((i.val : ℤ) + 1) * (14 * c) < 182 * c + 1 := by nlinarith [hi11, hc]
    rw [Int.emod_eq_of_lt ha0 hab]
    have hkc : (0 : ℤ) ≤ (i.val : ℤ) * (14 * c) := mul_nonneg hi0 hcpos
    have hup : (0 : ℤ) ≤ (11 - (i.val : ℤ)) * (14 * c) := mul_nonneg (by linarith) hcpos
    exact ⟨by nlinarith [hkc], by nlinarith [hup]⟩

/-- **The single-killer ladder reach bound.**  For `c ≥ 1`,
`reach(S_c) = sSup (margin S_c '' [0,1]) ≥ 14c/(182c+1)` — the single-killer balance value
`[1/13]·(182c)/(182c+1)`, via the covering-reach atom at `q = 182c+1`, rotation `14c`. -/
theorem singleKiller_reach (c : ℤ) (hc : 1 ≤ c) :
    (14 * (c : ℝ)) / (182 * (c : ℝ) + 1) ≤ sSup (margin (singleKiller c) '' Set.Icc (0 : ℝ) 1) := by
  have hq : (0 : ℤ) < 182 * c + 1 := by linarith
  have hc0 : (0 : ℤ) ≤ 14 * c := by linarith
  have hcq : (14 : ℤ) * c ≤ 182 * c + 1 := by linarith
  have h := reach_ge_of_covering13 (singleKiller c) (182 * c + 1) (14 * c) (14 * c)
    hq hc0 hcq (singleKiller_cover c hc)
  push_cast at h
  convert h using 2 <;> push_cast <;> ring

/-- **The whole single-killer ladder clears the covering-min floor:** for every `c ≥ 1`,
`reach(S_c) ≥ 14/183`, with equality only at the deep well `c = 1`.  (`14c/(182c+1) ≥
14/183 ⟺ 14 ≤ 14c ⟺ c ≥ 1`.)  So no single-killer covering family — extremal or not —
dips below `14/183`: the single-killer class is closed in the reach API. -/
theorem singleKiller_reach_ge_floor (c : ℤ) (hc : 1 ≤ c) :
    (14 : ℝ) / 183 ≤ sSup (margin (singleKiller c) '' Set.Icc (0 : ℝ) 1) := by
  have hcR : (1 : ℝ) ≤ (c : ℝ) := by exact_mod_cast hc
  have hden : (0 : ℝ) < 182 * (c : ℝ) + 1 := by linarith
  have hfloor : (14 : ℝ) / 183 ≤ 14 * (c : ℝ) / (182 * (c : ℝ) + 1) := by
    rw [div_le_iff₀ (show (0 : ℝ) < 183 by norm_num), div_mul_eq_mul_div, le_div_iff₀ hden]
    nlinarith [hcR]
  exact le_trans hfloor (singleKiller_reach c hc)

#print axioms singleKiller_cover
#print axioms singleKiller_reach
#print axioms singleKiller_reach_ge_floor

end TournamentH7.LRCSingleKillerLadder
