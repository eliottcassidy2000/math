/-
  TournamentH7.LRCOneSwapLadders — THE DEEP ONE-SWAP COVERING LADDERS, CLOSED BY FORMULA
  (kind-pasteur-2026-07-04-S5, HYP-4082/klein-S127).

  The one-swap covering stratum `F(j,X) = ({1,…,13} \ {j}) ∪ {X}` decomposes into 13 formula
  ladders, floored by the deep well `14/183` (klein-S127).  The FOUR deep ladders — where `j` is
  the unique base coverer of `q=j`, forcing `X` to be a multiple of `lcm(j,14)` — are:

    drop-13  X=182k  t*=14k/(182k+1)      M=14k/(182k+1)  (the deep well, covering-min floor)
    drop-12  X=84k   t*=(35k+2)/(84k+5)   M=7k/(84k+5)    (kps residue-liar; LRCResidueLiar.lean)
    drop-11  X=154k  t*=(56k+1)/(154k+3)  M=14k/(154k+3)
    drop-9   X=126k  t*=(28k+1)/(126k+5)  M=14k/(126k+5)

  Each is a rational-witness lonely certificate via a residue table (13 residues linear in k, each
  pinned into [14k, q−14k], binding at the lower runner and the coverer).  This file closes drop-9,
  drop-11, and drop-13 in Lean (drop-12 = LRCResidueLiar).  drop-13 is a bonus KERNEL-PURE deep-well
  certificate — the far-peel version (LRCFarPeelDeepWell) needed native_decide; this one does not.

  Kernel-pure (propext, Classical.choice, Quot.sound); no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRCSpread13
import TournamentH7.LRCResidueLiar

namespace LonelyRunner
namespace LRC14

/-- **THE RESIDUE-TABLE CERTIFICATE (single runner)**: at witness `p/q`, a runner `val` whose
residue `r = val·p − qq·q` is pinned into `[κ, q−κ]` with `q ≤ 14κ` clears the `q/14` bar. -/
lemma residue_key (p q kappa m val qq r : ℤ) (hq : 0 < q) (h14 : q ≤ 14 * kappa)
    (hEq : val * p = qq * q + r) (hlo : kappa ≤ r) (hhi : r ≤ q - kappa) :
    q ≤ 14 * |val * p - m * q| := by
  have hm : kappa ≤ |val * p - m * q| := lattice_dist_ge qq r hq hEq hlo hhi m
  linarith [hm, h14]

/-! ### drop-9 : `{1,…,8,10,11,12,13, 126k}`, lonely at `(28k+1)/(126k+5)`, `M = 14k/(126k+5)`. -/

def drop9 (k : ℤ) : Fin 13 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 126 * k]

theorem drop9_lonely (k : ℤ) (hk : 1 ≤ k) :
    Lonely 14 (drop9 k) (((28 * k + 1 : ℤ) : ℝ) / ((126 * k + 5 : ℤ) : ℝ)) := by
  apply lonely14_of_ratio (drop9 k) (28 * k + 1) (126 * k + 5) (by omega)
  intro i m
  fin_cases i <;> simp only [drop9]
  · exact residue_key _ _ (14*k) m 1  0 (28*k+1)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 2  0 (56*k+2)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 3  0 (84*k+3)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 4  0 (112*k+4) (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 5  1 (14*k)    (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 6  1 (42*k+1)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 7  1 (70*k+2)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 8  1 (98*k+3)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 10 2 (28*k)    (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 11 2 (56*k+1)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 12 2 (84*k+2)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 13 2 (112*k+3) (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m (126*k) (28*k-1) (112*k+5) (by omega) (by omega) (by ring) (by omega) (by omega)

/-! ### drop-11 : `{1,…,10,12,13, 154k}`, lonely at `(56k+1)/(154k+3)`, `M = 14k/(154k+3)`. -/

def drop11 (k : ℤ) : Fin 13 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 154 * k]

theorem drop11_lonely (k : ℤ) (hk : 1 ≤ k) :
    Lonely 14 (drop11 k) (((56 * k + 1 : ℤ) : ℝ) / ((154 * k + 3 : ℤ) : ℝ)) := by
  apply lonely14_of_ratio (drop11 k) (56 * k + 1) (154 * k + 3) (by omega)
  intro i m
  fin_cases i <;> simp only [drop11]
  · exact residue_key _ _ (14*k) m 1  0 (56*k+1)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 2  0 (112*k+2) (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 3  1 (14*k)    (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 4  1 (70*k+1)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 5  1 (126*k+2) (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 6  2 (28*k)    (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 7  2 (84*k+1)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 8  2 (140*k+2) (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 9  3 (42*k)    (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 10 3 (98*k+1)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 12 4 (56*k)    (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 13 4 (112*k+1) (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m (154*k) (56*k-1) (140*k+3) (by omega) (by omega) (by ring) (by omega) (by omega)

/-! ### drop-13 : `{1,…,12, 182k}` (the deep well), lonely at `14k/(182k+1)`, `M = 14k/(182k+1)`.
KERNEL-PURE — no native_decide (cf. LRCFarPeelDeepWell). -/

def drop13 (k : ℤ) : Fin 13 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 182 * k]

theorem drop13_lonely (k : ℤ) (hk : 1 ≤ k) :
    Lonely 14 (drop13 k) (((14 * k : ℤ) : ℝ) / ((182 * k + 1 : ℤ) : ℝ)) := by
  apply lonely14_of_ratio (drop13 k) (14 * k) (182 * k + 1) (by omega)
  intro i m
  fin_cases i <;> simp only [drop13]
  · exact residue_key _ _ (14*k) m 1  0 (14*k)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 2  0 (28*k)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 3  0 (42*k)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 4  0 (56*k)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 5  0 (70*k)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 6  0 (84*k)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 7  0 (98*k)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 8  0 (112*k) (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 9  0 (126*k) (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 10 0 (140*k) (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 11 0 (154*k) (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m 12 0 (168*k) (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (14*k) m (182*k) (14*k-1) (168*k+1) (by omega) (by omega) (by ring) (by omega) (by omega)

/-! ### drop-8 : `{1,…,7,9,…,13, 56k}`, lonely at `(7k+1)/(56k+7)`, `M = 7k/(56k+7)`. -/

def drop8 (k : ℤ) : Fin 13 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 56 * k]

theorem drop8_lonely (k : ℤ) (hk : 1 ≤ k) :
    Lonely 14 (drop8 k) (((7 * k + 1 : ℤ) : ℝ) / ((56 * k + 7 : ℤ) : ℝ)) := by
  apply lonely14_of_ratio (drop8 k) (7 * k + 1) (56 * k + 7) (by omega)
  intro i m
  fin_cases i <;> simp only [drop8]
  · exact residue_key _ _ (7*k) m 1  0 (7*k+1)   (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 2  0 (14*k+2)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 3  0 (21*k+3)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 4  0 (28*k+4)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 5  0 (35*k+5)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 6  0 (42*k+6)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 7  0 (49*k+7)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 9  1 (7*k+2)   (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 10 1 (14*k+3)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 11 1 (21*k+4)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 12 1 (28*k+5)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 13 1 (35*k+6)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m (56*k) (7*k) (7*k) (by omega) (by omega) (by ring) (by omega) (by omega)

/-! ### drop-10 : `{1,…,9,11,12,13, 70k}`, lonely at `(21k+2)/(70k+7)`, `M = 7k/(70k+7)`. -/

def drop10 (k : ℤ) : Fin 13 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 70 * k]

theorem drop10_lonely (k : ℤ) (hk : 1 ≤ k) :
    Lonely 14 (drop10 k) (((21 * k + 2 : ℤ) : ℝ) / ((70 * k + 7 : ℤ) : ℝ)) := by
  apply lonely14_of_ratio (drop10 k) (21 * k + 2) (70 * k + 7) (by omega)
  intro i m
  fin_cases i <;> simp only [drop10]
  · exact residue_key _ _ (7*k) m 1  0 (21*k+2)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 2  0 (42*k+4)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 3  0 (63*k+6)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 4  1 (14*k+1)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 5  1 (35*k+3)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 6  1 (56*k+5)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 7  2 (7*k)     (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 8  2 (28*k+2)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 9  2 (49*k+4)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 11 3 (21*k+1)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 12 3 (42*k+3)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m 13 3 (63*k+5)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (7*k) m (70*k) (21*k-1) (63*k+7) (by omega) (by omega) (by ring) (by omega) (by omega)

/-- Concrete deep well: `{1,…,12,182}` is lonely at `14/183` (kernel-pure). -/
theorem deepWell183_lonely : Lonely 14 (drop13 1) ((14 : ℝ) / 183) := by
  have h := drop13_lonely 1 (by norm_num); norm_num at h; exact h

#print axioms drop8_lonely
#print axioms drop9_lonely
#print axioms drop10_lonely
#print axioms drop11_lonely
#print axioms drop13_lonely
#print axioms deepWell183_lonely

end LRC14
end LonelyRunner
