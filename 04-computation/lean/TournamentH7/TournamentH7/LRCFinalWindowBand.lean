/-
klein-2026-07-01-S89 (HYP-3844) -- the band arithmetic of the K=0 final-window lemma.

Concave (overtaking/engulfment) kinks of the lonely distribution function Λ_S at n = 14
occur only at radii d/(w-v) with v < w ∈ S, d ≥ 1 (THM-592). A kink radius lies in the
open final window (1/15, 1/14) iff 14d < w - v < 15d. This file formalizes the arithmetic
shell: the d = 1 band is EMPTY (no integer strictly between 14 and 15), hence any concave
kink in the window needs d ≥ 2 and w - v ≥ 29; in particular a set of diameter ≤ 28 has
NO concave kink radius in the window (K = 0 there — the geometric identification of kinks
with these radii is THM-592, cited, not formalized here).
-/
import Mathlib

namespace LRCFinalWindowBand

/-- The final-window band condition at n = 14: the kink radius `d/q` lies in the open
window `(1/15, 1/14)` iff `14d < q < 15d` (for positive `d`, `q`). -/
def inBand (d q : ℕ) : Prop := 14 * d < q ∧ q < 15 * d

/-- **d = 1 emptiness**: there is no difference `q` with a `d = 1` kink radius in the
open final window — no integer lies strictly between 14 and 15. -/
theorem d1_band_empty (q : ℕ) : ¬ inBand 1 q := by
  unfold inBand; omega

/-- Any band membership forces `d ≥ 2` and `q ≥ 29`. -/
theorem band_forces_29 {d q : ℕ} (hd : 0 < d) (h : inBand d q) : 2 ≤ d ∧ 29 ≤ q := by
  unfold inBand at h; omega

/-- **The diameter-28 corollary (arithmetic shell of the K = 0 lemma).** If every
pairwise difference of the speed set is at most 28, no concave kink radius `d/(w-v)`
lies in the open final window: the band condition fails for every pair and every `d`. -/
theorem no_band_kink_of_diam_le {q d : ℕ} (hq : q ≤ 28) (hd : 0 < d) :
    ¬ inBand d q := by
  unfold inBand; omega

/-- The bands are exactly `{29}, {43,44}, {57,58,59}, ...`: for `d ≥ 2` the band
`(14d, 15d)` contains precisely the `d - 1` integers `14d+1, ..., 15d-1`. -/
theorem band_card (d : ℕ) (hd : 2 ≤ d) :
    ((Finset.Ioo (14 * d) (15 * d)).card) = d - 1 := by
  rw [Nat.card_Ioo]; omega

end LRCFinalWindowBand
