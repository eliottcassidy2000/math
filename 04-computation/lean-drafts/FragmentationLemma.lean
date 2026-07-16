/- SUPERSEDED (klein-2026-07-16-S316): this draft's two sorries are CLOSED, sorry-free, in
   04-computation/lean/TournamentH7/TournamentH7/FragmentationLemma.lean  (builds green in
   the TournamentH7 project; includes death-star-S30's badArcs_periodic).  NOTE: this
   draft's step-1 proof plan (arc counting) is flawed — the count can hit floor(Lw)+2;
   the proved route is the window lemma + floor(Lw)+1 tiling.  Kept for history. -/
/- FragmentationLemma.lean -- mac-mini-2026-07-16-S123.
   THM-883's core counting inequality, the first Lean target of the covering program.
   Statement: an arc-grid of modulus w >= 1 and radius lam = 1/13 meets an interval I
   of length L in measure at most (L*w + 1) * (2*lam/w).
   Kernel-friendly formulation over Q; the arc set is the union over a in [0, w) of
   open intervals (a/w - lam/w, a/w + lam/w) intersected with I.
   Proof plan (two sorries mark the classical steps):
   1. the number of integers a with (a/w - lam/w, a/w + lam/w) ∩ I ≠ ∅ is at most
      ⌊L*w⌋ + 1  (grid spacing 1/w over an interval of length L);
   2. each arc has measure ≤ 2*lam/w; subadditivity finishes. -/
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic
import Mathlib.Analysis.SpecialFunctions.Basic

open MeasureTheory

namespace LRC14

/-- The bad arc set of modulus `w` at radius `lam`, within the unit circle lifted to ℝ. -/
def badArcs (w : ℕ) (lam : ℝ) : Set ℝ :=
  ⋃ a : ℤ, Set.Ioo ((a : ℝ) / w - lam / w) ((a : ℝ) / w + lam / w)

/-- THM-883 Lemma 1 (fragmentation): the arc-grid meets an interval `I = [x, x+L]`
    in measure at most `(L * w + 1) * (2 * lam / w)`. -/
theorem fragmentation (w : ℕ) (hw : 1 ≤ w) (lam L x : ℝ)
    (hlam : 0 < lam) (hL : 0 ≤ L) :
    volume (badArcs w lam ∩ Set.Icc x (x + L))
      ≤ ENNReal.ofReal ((L * w + 1) * (2 * lam / w)) := by
  -- Step 1: only arcs with center a/w in [x - lam/w, x + L + lam/w] can meet I;
  -- there are at most ⌊L*w⌋ + 1 + (edge corrections absorbed by the +1 for lam < 1/2)
  -- such integers a.
  -- Step 2: each arc has volume 2*lam/w; countable subadditivity.
  sorry

/-- Corollary (killer budget, THM-883 Lemma 2 shape): if a family of moduli
    `w₁ ≤ … ≤ w_j` covers `[x, x+L]` by their arc-grids, then
    `L * (1 - 2*j*lam) ≤ 2*lam * Σ 1/wᵢ`. -/
theorem killer_budget (j : ℕ) (ws : Fin j → ℕ) (hws : ∀ i, 1 ≤ ws i)
    (lam L x : ℝ) (hlam : 0 < lam) (hL : 0 ≤ L)
    (hcover : Set.Icc x (x + L) ⊆ ⋃ i, badArcs (ws i) lam) :
    L * (1 - 2 * j * lam) ≤ 2 * lam * ∑ i, (1 : ℝ) / ws i := by
  -- volume monotonicity + fragmentation applied to each modulus.
  sorry

end LRC14
