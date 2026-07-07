/-
  TournamentH7.LRCEvenBranchWitness  (mac-mini-2026-07-06-S29, HYP-4572)

  THE EVEN-BRANCH CLEARANCE, formalized at N=12 (the LRC(14) case).

  opus-S118's canonical mediant family is  F(N) = {1,…,N}\{N−1} ∪ {3(N−1)}.
  HYP-4572 (mac-mini-S28) trichotomy: for N EVEN,
        M(F(N)) = 3/(3N−1)     (binder speed 2 at denominator Q = 3N−1),
  which is STRICTLY ABOVE the second-gap top 2/(2N+1) — so F(N) is NOT a gap
  member for even N.  The proof is a clean residue clearance: at t = b/(3N−1)
  with b = 3·2⁻¹, every speed v∈F(N) has  v·b ∉ {0,±1,±2} mod (3N−1)  because the
  only solution of 3v ∈ {0,±2,±4} mod (3N−1) inside {1,…,N}∪{3(N−1)} is v = N−1,
  which is exactly the removed element (using 3⁻¹ = N and 3(N−1) = Q−2).

  This file gives the N=12 instance (Q = 35, b = 19), fully:  the 12-speed family
  F12 = {1,…,10,12,33} clears the threshold 3/35 at t = 19/35, so its reach (the
  attained sup of the margin over the circle) is ≥ 3/35 > 2/25 = the gap top.
  Hence the canonical construction provably fails to be a gap member at N=12 —
  by PARITY (12 even), not by 38 = 2·19 compositeness.

  Kernel-pure (decide-free; `omega` + `norm_num`).  Axioms are the standard trio.
-/
import TournamentH7.LRCWitnessAttainment

open TournamentH7.LRCWitness
open scoped Topology

namespace TournamentH7.LRCEvenBranch

/-- The canonical mediant family at N=12:  `{1,…,10, 12, 33}` = `{1,…,12}\{11} ∪ {33}`.
Index `i = 0..9 ↦ i+1`, `i = 10 ↦ 12`, `i = 11 ↦ 33`. -/
def F12 : Fin 12 → ℤ := fun i =>
  if i.val = 11 then 33 else if i.val = 10 then 12 else (i.val : ℤ) + 1

/-- `distZ y ≤ 1/2` (the nearest integer is within `1/2`). -/
lemma distZ_le_half (y : ℝ) : distZ y ≤ 1 / 2 := by
  have hmem : ((round y : ℤ) : ℝ) ∈ Set.range ((↑) : ℤ → ℝ) := ⟨round y, rfl⟩
  calc distZ y ≤ dist y ((round y : ℤ) : ℝ) := Metric.infDist_le_dist_of_mem hmem
    _ = |y - (round y : ℝ)| := Real.dist_eq _ _
    _ ≤ 1 / 2 := abs_sub_round y

/-- The margin never exceeds `1/2` (bounds the reach for `BddAbove`). -/
lemma margin_le_half (v : Fin 12 → ℤ) (t : ℝ) : margin v t ≤ 1 / 2 := by
  unfold margin
  exact le_trans (Finset.inf'_le _ (Finset.mem_univ (0 : Fin 12))) (distZ_le_half _)

/-- **Integer → real bridge** at `Q = 35`, threshold `3`.  If `3 ≤ |v·19 − k·35|` for
all integers `k`, then `3/35 ≤ |v·(19/35) − m|`. -/
lemma bridge (v : ℤ) (hb : ∀ k : ℤ, (3 : ℤ) ≤ |v * 19 - k * 35|) (m : ℤ) :
    (3 : ℝ) / 35 ≤ |(v : ℝ) * (19 / 35) - (m : ℝ)| := by
  have key : (3 : ℝ) ≤ |(v : ℝ) * 19 - (m : ℝ) * 35| := by
    have h : (3 : ℝ) ≤ ((|v * 19 - m * 35| : ℤ) : ℝ) := by exact_mod_cast hb m
    rw [Int.cast_abs] at h; push_cast at h; exact h
  have hfac : (v : ℝ) * (19 / 35) - (m : ℝ) = ((v : ℝ) * 19 - (m : ℝ) * 35) / 35 := by
    rw [eq_div_iff (by norm_num : (35 : ℝ) ≠ 0)]; ring
  rw [hfac, abs_div, show |(35 : ℝ)| = 35 by norm_num,
    le_div_iff₀ (by norm_num : (0 : ℝ) < 35)]
  rw [show (3 : ℝ) / 35 * 35 = 3 by norm_num]; exact key

/-- **The even-branch witness at N=12.**  Every speed of `F12` is `≥ 3/35` from ℤ at
`t = 19/35`; equivalently `3/35 ≤ margin F12 (19/35)`. -/
theorem F12_margin : (3 : ℝ) / 35 ≤ margin F12 (19 / 35) := by
  rw [le_margin_iff]
  intro i m
  -- Reduce `F12 i` to its literal value (12 cases), keeping the product intact,
  -- then the integer clearance `3 ≤ |v·19 − k·35|` is one `omega` per value.
  have hval : ∀ k : ℤ, (3 : ℤ) ≤ |F12 i * 19 - k * 35| := by
    intro k
    have hcase : F12 i = 1 ∨ F12 i = 2 ∨ F12 i = 3 ∨ F12 i = 4 ∨ F12 i = 5 ∨
        F12 i = 6 ∨ F12 i = 7 ∨ F12 i = 8 ∨ F12 i = 9 ∨ F12 i = 10 ∨
        F12 i = 12 ∨ F12 i = 33 := by
      fin_cases i <;> simp only [F12] <;> norm_num
    rcases hcase with h | h | h | h | h | h | h | h | h | h | h | h <;> rw [h, le_abs] <;> omega
  exact bridge (F12 i) hval m

/-- **Reach lower bound.**  The attained sup of the margin over `[0,1]` exceeds the
second-gap top `2/25`.  So the canonical mediant family at N=12 is NOT a gap
member: its `M` lies above the gap `(1/13, 2/25)`. -/
theorem F12_reach_above_gap :
    (2 : ℝ) / 25 < sSup (margin F12 '' Set.Icc (0 : ℝ) 1) := by
  have hmem : margin F12 (19 / 35) ∈ margin F12 '' Set.Icc (0 : ℝ) 1 :=
    ⟨19 / 35, Set.mem_Icc.mpr ⟨by norm_num, by norm_num⟩, rfl⟩
  have hbdd : BddAbove (margin F12 '' Set.Icc (0 : ℝ) 1) := by
    refine ⟨1 / 2, ?_⟩
    rintro _ ⟨t, _, rfl⟩
    exact margin_le_half F12 t
  have hle : margin F12 (19 / 35) ≤ sSup (margin F12 '' Set.Icc (0 : ℝ) 1) :=
    le_csSup hbdd hmem
  have hw : (3 : ℝ) / 35 ≤ margin F12 (19 / 35) := F12_margin
  have hnum : (2 : ℝ) / 25 < (3 : ℝ) / 35 := by norm_num
  linarith

#print axioms F12_margin
#print axioms F12_reach_above_gap

end TournamentH7.LRCEvenBranch
