/-
  TournamentH7.LRCTwoKillerLadder — THE TWO-KILLER (drop-two) COVERING-FAR LADDER
  (klein-2026-07-04-S130, HYP-4091).

  kps's LRCOneSwapLadders certifies the SINGLE-swap hexad: drop ONE element `j ∈ {8,…,13}`
  from `{1,…,13}` and add ONE killer `lcm(j,14)·k` (the unique large coverer of both `j`
  and `14`).  The NEXT stratum of `CoveringFarLonely` is the MULTI-killer families, which
  the single-swap hexad does not reach.  This file certifies the first two-killer class:

        twoKiller K := {1,…,11, 14, 156·K}                              (13 speeds)

  — drop BOTH `12` and `13`, add TWO killers: `14` (covers `q=14`) and `156K = 12·13·K`
  (covers `q=12` and `q=13`).  This is covering for every `K ≥ 1` and has the far entry
  `156K > 20`, so it is a genuine `CoveringFarLonely 20` instance beyond the hexad.

  The whole ladder is lonely: for every `K ≥ 1`,

        twoKiller K is lonely at t = 13K/(156K+1),  margin M = 13K/(156K+1) > 1/14.

  Same engine as the hexad: one rational witness, a residue table linear in `K`,
  `residue_key`/`lattice_dist_ge`.  Residue table at p = 13K, Q = 156K+1, floor κ = 13K:
    v = 1,…,11 : v·(13K) = 0·Q + 13K·v          (13K ≤ 13Kv ≤ 143K, all in [κ, Q−κ])
    v = 14     : 14·(13K) = 1·Q + (26K−1)        (26K−1 ∈ [13K, 143K+1])
    v = 156K   : 156K·(13K) = (13K−1)·Q + (143K+1)   (143K+1 = Q − κ, the binding end)
  and 14·κ = 182K ≥ 156K+1 = Q ⟺ 26K ≥ 1 ⟺ K ≥ 1.  Binding pair: runner 1 (bottom) and
  killer 156K (top) — the same 2-point equioscillation THM-618 identifies, here with the
  base `{1,…,11}` (optimum 1/12) and a split killer set.

  Kernel-pure (propext, Classical.choice, Quot.sound); no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRCSpread13
import TournamentH7.LRCResidueLiar
import TournamentH7.LRCOneSwapLadders

namespace LonelyRunner
namespace LRC14

/-- The two-killer / drop-`{12,13}` family `{1,…,11, 14, 156·K}` (13 speeds).
`156 = 12·13 = lcm(12,13)` covers `q=12,13`; `14` covers `q=14`; base `{1,…,11}`. -/
def twoKiller (K : ℤ) : Fin 13 → ℤ :=
  ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 156 * K]

/-- **THE TWO-KILLER LADDER IS LONELY** at `t = 13K/(156K+1)` for every `K ≥ 1`:
the first multi-killer (drop-two) `CoveringFarLonely` class, margin `13K/(156K+1) > 1/14`,
strictly increasing in `K`.  Kernel-pure, no `native_decide`. -/
theorem twoKiller_lonely (K : ℤ) (hK : 1 ≤ K) :
    Lonely 14 (twoKiller K) (((13 * K : ℤ) : ℝ) / ((156 * K + 1 : ℤ) : ℝ)) := by
  apply lonely14_of_ratio (twoKiller K) (13 * K) (156 * K + 1) (by omega)
  intro i m
  fin_cases i <;> simp only [twoKiller]
  · exact residue_key _ _ (13*K) m 1  0 (13*K)   (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (13*K) m 2  0 (26*K)   (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (13*K) m 3  0 (39*K)   (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (13*K) m 4  0 (52*K)   (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (13*K) m 5  0 (65*K)   (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (13*K) m 6  0 (78*K)   (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (13*K) m 7  0 (91*K)   (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (13*K) m 8  0 (104*K)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (13*K) m 9  0 (117*K)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (13*K) m 10 0 (130*K)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (13*K) m 11 0 (143*K)  (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (13*K) m 14 1 (26*K-1) (by omega) (by omega) (by ring) (by omega) (by omega)
  · exact residue_key _ _ (13*K) m (156*K) (13*K-1) (143*K+1) (by omega) (by omega) (by ring) (by omega) (by omega)

/-- Concrete base rung (`K=1`): the two-killer family `{1,…,11,14,156}` is lonely at
`13/157`, margin `13/157 > 1/14` — a `CoveringFarLonely` instance not on the single-swap
hexad. -/
theorem twoKiller156_lonely :
    Lonely 14 (twoKiller 1) ((13 : ℝ) / 157) := by
  have h := twoKiller_lonely 1 (by norm_num)
  norm_num at h
  exact h

-- Kernel-purity: expect only [propext, Classical.choice, Quot.sound].
#print axioms twoKiller_lonely
#print axioms twoKiller156_lonely

end LRC14
end LonelyRunner
