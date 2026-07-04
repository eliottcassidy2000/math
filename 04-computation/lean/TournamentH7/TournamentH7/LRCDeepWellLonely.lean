/-
  TournamentH7.LRCDeepWellLonely — WIRING THE DEEP-WELL WITNESS TO `Lonely 14`
  (klein-2026-07-03-S119b, HYP-4065; renumbered from HYP-4062 -- collided with kind-pasteur-S37).

  Instantiates the general-n covering-min extremizer witness (LRCDeepWellWitness) at the
  concrete LRC(14) predicate `LonelyRunner.Lonely 14 v t` (`v : Fin 13 → ℤ`, `t : ℝ`,
  `∀ i m, 1/14 ≤ |v i * t - m|`).  Certifies the covering-min extremizer inside the chain:

    `deepWell14 = {1,2,…,12, 182}`  is  `Lonely 14`  at  `t* = 14/183`.

  `182 = 14·13 = n(n-1)` is the pronic defect; `183 = Phi6(14)`.  Every runner is at real
  distance `≥ 14/183 > 1/14` from the integers at `t*`, from the integer witness bounds
  `ap_runner_lonely` / `defect_runner_lonely` (n = 14) via the elementary real bridge below.
  Sorry-free.
-/
import TournamentH7.LonelyRunner
import TournamentH7.LRCDeepWellWitness

namespace LonelyRunner.DeepWell

open LonelyRunner

/-- `Phi6(14) = 183`. -/
theorem phi6_14 : phi6 14 = 183 := by unfold phi6; norm_num

/-- The deep well at `n = 14`: speeds `{1,…,12}` and the pronic defect `182 = 14·13`. -/
def deepWell14 : Fin 13 → ℤ := fun i => if i.val = 12 then 182 else (i.val : ℤ) + 1

/-- **Real bridge.**  An integer distance bound `∀ m, 14 ≤ |v·14 − m·183|` upgrades to the
  real loneliness threshold `1/14 ≤ |v·(14/183) − m|` (using `183/14 = 13.07… ≤ 14`). -/
theorem bridge (v : ℤ) (hb : ∀ m : ℤ, (14 : ℤ) ≤ |v * 14 - m * 183|) (m : ℤ) :
    (1 : ℝ) / 14 ≤ |(v : ℝ) * ((14 : ℝ) / 183) - m| := by
  have key : (14 : ℝ) ≤ |(v : ℝ) * 14 - m * 183| := by
    have h : (14 : ℝ) ≤ ((|v * 14 - m * 183| : ℤ) : ℝ) := by exact_mod_cast hb m
    rw [Int.cast_abs] at h; push_cast at h; exact h
  have hfac : (v : ℝ) * ((14 : ℝ) / 183) - m = ((v : ℝ) * 14 - m * 183) / 183 := by
    rw [eq_div_iff (by norm_num : (183 : ℝ) ≠ 0)]; ring
  rw [hfac, abs_div, show |(183 : ℝ)| = 183 by norm_num, le_div_iff₀ (by norm_num : (0 : ℝ) < 183)]
  have hthr : (1 : ℝ) / 14 * 183 ≤ 14 := by norm_num
  linarith [key]

/-- **The covering-min extremizer is `Lonely 14` at `t* = 14/183`.**  Wires the general-n
  deep-well witness (`ap_runner_lonely` / `defect_runner_lonely`, `n = 14`) into the concrete
  LRC(14) predicate.  This certifies `{1,…,12,182}` — the unique covering-min extremal
  family, `M = 14/183 > 1/14` — inside the chain's own vocabulary. -/
theorem deepWell14_lonely : Lonely 14 deepWell14 ((14 : ℝ) / 183) := by
  intro i m
  simp only [Nat.cast_ofNat]
  by_cases h : i.val = 12
  · -- defect runner 182 = 14·13 = n(n-1)
    have hv : deepWell14 i = 182 := by simp [deepWell14, h]
    rw [hv]
    refine bridge 182 (fun m => ?_) m
    have hd := defect_runner_lonely (n := 14) (by norm_num) m
    rw [phi6_14] at hd
    -- hd : 14 ≤ |14*(14-1)*14 - m*183| ; goal needs |182*14 - m*183|
    have he : (14 : ℤ) * (14 - 1) * 14 = 182 * 14 := by norm_num
    rw [he] at hd; exact hd
  · -- AP runner j = i.val + 1 ∈ [1, 12]
    have hi : i.val ≤ 11 := by omega
    have hv : deepWell14 i = (i.val : ℤ) + 1 := by simp [deepWell14, h]
    rw [hv]
    refine bridge ((i.val : ℤ) + 1) (fun m => ?_) m
    have ha := ap_runner_lonely (n := 14) (j := (i.val : ℤ) + 1) (by norm_num)
      (by omega) (by omega) m
    rw [phi6_14] at ha
    exact ha

#print axioms deepWell14_lonely

end LonelyRunner.DeepWell
