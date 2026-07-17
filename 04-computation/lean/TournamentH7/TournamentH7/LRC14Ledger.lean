/- LRC14Ledger.lean — mac-mini-2026-07-16-S129.
   THE ASSEMBLY FILE of the LRC(14) formalization ladder.
   Imports every built rung, states the LRC(14) target as a Prop (no axiom, no
   sorry), and records the mapping from the target to the canon chain.

   BUILT RUNGS (oleans, kernel-checked):
   · FragmentationCount — THM-883's kernel (arc count; fragmentation measure bound)
   · KillerBudget       — the budget composition + the good-set floor
   · TrivialLoneliness  — END-TO-END existence: `trivial_LRC` proves the Lonely
                          Runner statement at the trivial constant 1/(6k)
   · TieSplitWalk       — THM-866's kernel (F3 flip arithmetic; pigeonhole)

   THE TARGET (stated below, not yet formally discharged): `LRC14` — thirteen
   positive speeds admit a 1/14-lonely time. The informal proof is the canon
   chain (covering: THM-922 sign-off + THM-883 rigidity + bands + low-M;
   S3 residual: THM-527 reformulation + THM-924 glue + THM-925/926 floor
   c₀ = 17/84; level-5 wall: THM-896/897 + THM-926 T2). The formal gap between
   `trivial_LRC` and `LRC14` is exactly the certified chain's constant
   improvement 1/(6·13) → 1/14; each canon page is finite/rational and slots
   into this namespace as a further rung. -/
import TournamentH7.FragmentationCount
import TournamentH7.KillerBudget
import TournamentH7.TrivialLoneliness
import TournamentH7.TieSplitWalk
import TournamentH7.UnitBudget

open MeasureTheory Set

namespace LRC14

/-- **The LRC(14) target statement.** Thirteen distinct positive integer speeds
    admit a time `t` with `‖w·t‖ ≥ 1/14` for every speed — stated via integer
    distance, matching `trivial_LRC`'s conclusion shape. -/
def LRC14 : Prop :=
  ∀ W : Finset ℕ, (∀ w ∈ W, 0 < w) → W.card = 13 →
    ∃ t ∈ Icc (0 : ℝ) 1, ∀ w ∈ W, ∀ a : ℤ, (1 : ℝ)/14 ≤ |(w : ℝ) * t - a|

/-- The ladder's current unconditional instance: the same statement at the
    trivial constant — for thirteen speeds, loneliness at gap `1/79` (any
    `lam < 1/78` works; `1/79` keeps the arithmetic decidable by `norm_num`). -/
theorem LRC13speeds_at_trivial_gap :
    ∀ W : Finset ℕ, (∀ w ∈ W, 0 < w) → W.card = 13 →
      ∃ t ∈ Icc (0 : ℝ) 1, ∀ w ∈ W, ∀ a : ℤ, (1 : ℝ)/79 ≤ |(w : ℝ) * t - a| := by
  intro W hW hcard
  refine trivial_LRC W hW (1/79) (by norm_num) (by norm_num) ?_
  rw [hcard]
  norm_num

end LRC14

namespace LRC14

/-- **The sharpened existence theorem** (uses the periodic unit budget: each speed
    costs exactly `2*lam` on `[0,1]`): total budget `2*lam*|W| < 1` suffices. -/
theorem exists_lonely_sharp (W : Finset ℕ) (hW : ∀ w ∈ W, 0 < w) (lam : ℝ)
    (hlam : 0 < lam) (hlam1 : lam < 1) (hsmall : 2 * lam * W.card < 1) :
    ∃ t ∈ Icc (0 : ℝ) 1, ∀ w ∈ W, ∀ a : ℤ, lam ≤ |(w : ℝ) * t - a| := by
  classical
  have hbad : volume ((⋃ w ∈ W, badSet w lam) ∩ Icc (0 : ℝ) 1)
      ≤ ENNReal.ofReal (2 * lam * W.card) := by
    have hdistrib : (⋃ w ∈ W, badSet w lam) ∩ Icc (0 : ℝ) 1
        = ⋃ w ∈ W, (badSet w lam ∩ Icc (0 : ℝ) 1) := by
      ext t
      simp only [Set.mem_inter_iff, Set.mem_iUnion, exists_prop]
      tauto
    rw [hdistrib]
    calc volume (⋃ w ∈ W, (badSet w lam ∩ Icc (0 : ℝ) 1))
        ≤ ∑ w ∈ W, volume (badSet w lam ∩ Icc (0 : ℝ) 1) :=
          MeasureTheory.measure_biUnion_finset_le _ _
      _ ≤ ∑ _w ∈ W, ENNReal.ofReal (2 * lam) :=
          Finset.sum_le_sum fun w hw => unit_bad_le w (hW w hw) lam hlam hlam1
      _ = ENNReal.ofReal (2 * lam * W.card) := by
          rw [Finset.sum_const, nsmul_eq_mul, ← ENNReal.ofReal_natCast (W.card),
              ← ENNReal.ofReal_mul (Nat.cast_nonneg W.card)]
          congr 1
          ring
  have hgood : (0 : ENNReal)
      < volume (Icc (0 : ℝ) 1 \ ⋃ w ∈ W, badSet w lam) := by
    have hvol : volume (Icc (0 : ℝ) 1) = ENNReal.ofReal 1 := by
      rw [Real.volume_Icc]
      norm_num
    have hcover : Icc (0 : ℝ) 1
        ⊆ (Icc (0 : ℝ) 1 \ ⋃ w ∈ W, badSet w lam)
          ∪ ((⋃ w ∈ W, badSet w lam) ∩ Icc (0 : ℝ) 1) := by
      intro t ht
      by_cases hb : t ∈ ⋃ w ∈ W, badSet w lam
      · exact Or.inr ⟨hb, ht⟩
      · exact Or.inl ⟨ht, hb⟩
    have hsplit := le_trans (MeasureTheory.measure_mono hcover)
      (MeasureTheory.measure_union_le _ _)
    rw [hvol] at hsplit
    have hlt : ENNReal.ofReal (2 * lam * W.card) < ENNReal.ofReal 1 :=
      (ENNReal.ofReal_lt_ofReal_iff one_pos).mpr hsmall
    by_contra hzero
    push_neg at hzero
    have h0 : volume (Icc (0 : ℝ) 1 \ ⋃ w ∈ W, badSet w lam) = 0 :=
      le_antisymm hzero (zero_le _)
    rw [h0, zero_add] at hsplit
    exact absurd (le_trans hsplit hbad) (not_le.mpr hlt)
  rcases MeasureTheory.nonempty_of_measure_ne_zero (ne_of_gt hgood) with ⟨t, htIcc, htBad⟩
  refine ⟨t, htIcc, fun w hw a => ?_⟩
  have hnot : t ∉ badSet w lam := fun hmem => htBad (Set.mem_biUnion hw hmem)
  exact dist_int_ge w (hW w hw) lam hlam t hnot a

/-- **Thirteen speeds at gap 1/27** — the ladder's sharpened unconditional instance
    (2·(1/27)·13 = 26/27 < 1). The remaining program: 1/27 → 1/14 along the canon chain. -/
theorem LRC13speeds_at_gap_27 :
    ∀ W : Finset ℕ, (∀ w ∈ W, 0 < w) → W.card = 13 →
      ∃ t ∈ Icc (0 : ℝ) 1, ∀ w ∈ W, ∀ a : ℤ, (1 : ℝ)/27 ≤ |(w : ℝ) * t - a| := by
  intro W hW hcard
  refine exists_lonely_sharp W hW (1/27) (by norm_num) (by norm_num) ?_
  rw [hcard]
  norm_num

end LRC14
