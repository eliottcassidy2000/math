/-
  TournamentH7.LRCDeepCountExact — THE EXACT DEEP COUNT AND THE TIGHT FAMILY
  THROUGH THE FUNNEL (death-star-2026-07-17-S54, HYP-7295).

  Two capstones on THM-985's two-circle theorem:

  * `deep_count_exact` — the deep count of the canonical family in closed
    form, for EVERY q ≥ 2:
      `#deep(q) = 2·B + (B + 1 − (q + B) % 2)`,  `B = (q−1)/84`
    — circle I contributes the two arcs `[1, B] ∪ [q−B, q−1]` (2B), circle II
    the parity-corrected central arc `⌈(q−B)/2⌉ … ⌊(q+B)/2⌋`, and the two
    circles are PROVABLY DISJOINT for all q (3B < q).  Every constant is a
    literal, so the entire count is `omega`-native.  (Recon: exact on all
    q ∈ [2, 2000).)

  * `canonical_lonely` — **THE TIGHT FAMILY THROUGH THE FUNNEL**: at q = 14
    the formula gives ONE deep multiplier (p = 7 — the six even speeds at the
    half-integer, the mirror-fixed point) against SIX live multipliers
    (p ∈ {1,3,5,9,11,13} — the equality witnesses of the tight case), so
    THM-951's census pipeline fires by kernel `decide`:
    the canonical family (1,…,13) — the equality case of LRC(14) — is lonely,
    PROVED THROUGH THE B5 FUNNEL.  The loop from the S42 census machinery
    through the S52–53 circles closes.

  Kernel-pure: no `sorry`, no `native_decide` (`decide` on the q = 14 census
  gates only).
-/
import Mathlib
import TournamentH7.LRCTwoCircleII
import TournamentH7.LRCDeepCertificate

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- The pointwise ℤ→ℕ bridge for the two circles. -/
theorem circles_nat_iff (q p : ℕ) (hp : 0 < p) (hpq : p < q) :
    ((84 * (p : ℤ) < q ∨ 84 * ((q : ℤ) - p) < q) ∨ 84 * |2 * (p : ℤ) - q| < q)
      ↔ ((84 * p < q ∨ 84 * (q - p) < q) ∨
         (q ≤ 2 * p ∧ 84 * (2 * p - q) < q) ∨ (2 * p < q ∧ 84 * (q - 2 * p) < q)) := by
  rcases abs_cases (2 * (p : ℤ) - q) with ⟨he, hc⟩ | ⟨he, hc⟩ <;> rw [he] <;> omega

open Classical in
/-- **THE EXACT DEEP COUNT**: for every `q ≥ 2`,
`#deep(q) = 2·B + (B + 1 − (q + B) % 2)` with `B = (q−1)/84`. -/
theorem deep_count_exact (q : ℕ) (hq : 2 ≤ q) :
    ((Finset.Ioo 0 q).filter fun p => 6 ≤ bandCount canonical q p).card
      = 2 * ((q - 1) / 84) + ((q - 1) / 84 + 1 - (q + (q - 1) / 84) % 2) := by
  have hq0 : 0 < q := by omega
  -- rewrite the filter through the two-circle theorem and the ℕ-bridge
  have hfilter : ((Finset.Ioo 0 q).filter fun p => 6 ≤ bandCount canonical q p)
      = ((Finset.Ioo 0 q).filter fun p =>
          (84 * p < q ∨ 84 * (q - p) < q) ∨
          (q ≤ 2 * p ∧ 84 * (2 * p - q) < q) ∨ (2 * p < q ∧ 84 * (q - 2 * p) < q)) := by
    apply Finset.filter_congr
    intro p hp
    rw [Finset.mem_Ioo] at hp
    rw [deep_iff_circles q p hq0 hp.1 hp.2, circles_nat_iff q p hp.1 hp.2]
  rw [hfilter]
  -- split into the two disjoint circles
  have hsplit : ((Finset.Ioo 0 q).filter fun p =>
      (84 * p < q ∨ 84 * (q - p) < q) ∨
      (q ≤ 2 * p ∧ 84 * (2 * p - q) < q) ∨ (2 * p < q ∧ 84 * (q - 2 * p) < q))
      = ((Finset.Ioo 0 q).filter fun p => 84 * p < q ∨ 84 * (q - p) < q)
        ∪ ((Finset.Ioo 0 q).filter fun p =>
            (q ≤ 2 * p ∧ 84 * (2 * p - q) < q) ∨ (2 * p < q ∧ 84 * (q - 2 * p) < q)) := by
    rw [← Finset.filter_or]
  rw [hsplit]
  have hdisj : Disjoint
      ((Finset.Ioo 0 q).filter fun p => 84 * p < q ∨ 84 * (q - p) < q)
      ((Finset.Ioo 0 q).filter fun p =>
        (q ≤ 2 * p ∧ 84 * (2 * p - q) < q) ∨ (2 * p < q ∧ 84 * (q - 2 * p) < q)) := by
    rw [Finset.disjoint_left]
    intro p h1 h2
    rw [Finset.mem_filter, Finset.mem_Ioo] at h1 h2
    omega
  rw [Finset.card_union_of_disjoint hdisj]
  -- circle I: the two boundary arcs
  have hI : ((Finset.Ioo 0 q).filter fun p => 84 * p < q ∨ 84 * (q - p) < q)
      = Finset.Icc 1 ((q - 1) / 84) ∪ Finset.Icc (q - (q - 1) / 84) (q - 1) := by
    ext p
    simp only [Finset.mem_filter, Finset.mem_Ioo, Finset.mem_union, Finset.mem_Icc]
    omega
  have hIdisj : Disjoint (Finset.Icc 1 ((q - 1) / 84))
      (Finset.Icc (q - (q - 1) / 84) (q - 1)) := by
    rw [Finset.disjoint_left]
    intro p h1 h2
    rw [Finset.mem_Icc] at h1 h2
    omega
  -- circle II: the central arc
  have hII : ((Finset.Ioo 0 q).filter fun p =>
      (q ≤ 2 * p ∧ 84 * (2 * p - q) < q) ∨ (2 * p < q ∧ 84 * (q - 2 * p) < q))
      = Finset.Icc ((q - (q - 1) / 84 + 1) / 2) ((q + (q - 1) / 84) / 2) := by
    ext p
    simp only [Finset.mem_filter, Finset.mem_Ioo, Finset.mem_Icc]
    omega
  rw [hI, Finset.card_union_of_disjoint hIdisj, hII,
    Nat.card_Icc, Nat.card_Icc, Nat.card_Icc]
  omega

/-- **THE TIGHT FAMILY THROUGH THE FUNNEL**: the canonical family (1,…,13) —
the equality case of LRC(14) — is lonely, via the census pipeline at the
resonant modulus q = 14 (one deep multiplier p = 7 vs six live equality
witnesses). -/
theorem canonical_lonely : ∃ t : ℝ, Lonely 14 canonical t := by
  apply lonely_of_census canonical 14 (by norm_num)
    (fun i => by
      show (i : ℤ) + 1 ≠ 0
      have h := Int.natCast_nonneg (i : ℕ)
      omega)
  · -- CoverageCapped canonical 14 6
    decide
  · -- the race: #{bandCount = 6} < liveCount
    decide

/-! ## Axiom audit -/
#print axioms circles_nat_iff
#print axioms deep_count_exact
#print axioms canonical_lonely

end LRC14Concrete
end LonelyRunner
