/- LRCClusterGap.lean -- opus-2026-07-17-S338 (HYP-7210 / THM-955). DRAFT.
   THE CLUSTER GAP LEMMA: any k ≤ 6 positive speeds leave, in every window
   [a, b], a closed subinterval safe for all of them (arcSafe (1/14) w 0),
   of width ≥ [(1 − k/7)L − 2k/(7m)] / (1 + Σ(wL + 2)).

   Verified exactly 400/400 (cluster_gap_verify_opus_S336); paper THM-955.
   CONSUMER: LRCLacunaryNest.window_tail_glue / norm_glue — this lemma turns
   their base-certificate hypothesis into an analytic guarantee, so
   ≤6-block + gap + 7/3-tail families are lonely with no search.

   Three layers; (L1) is the one genuine formalization chunk (the sort):

   (L1) sorted_gap_pigeonhole: [a,b] minus N open rational intervals of
        total clipped length ≤ B contains a closed subinterval of width
        ≥ (L − B)/(N + 1).
        PLAN: mergeSort the 2N endpoints clipped to [a,b]; the consecutive
        segments of (a :: sorted ++ [b]) partition [a,b]; each segment is
        entirely inside some tooth or disjoint from all (no endpoint
        strictly interior); the tooth-contained segments have total length
        ≤ B; among the ≥ N+1 free segments... (count: 2N+2 points give
        2N+1 segments; ≤ N of them can be tooth-cores; pigeonhole on the
        rest). Finset.max' extraction.
   (L2) teeth enumeration: comb w in [a,b] = teeth centered j/w for
        j ∈ Finset.Icc ⌈w a⌉ ⌊w b⌋ (± the two edge teeth), each open of
        width 1/(7w) — the fragmentation counting, reused.
   (L3) assembly: the found free segment [a′, b′] avoids every tooth of
        every comb ⟹ arcSafe (1/14) w 0 a′ b′ per comb (the floor
        bookkeeping of nested_gap_step, in reverse). -/
import Mathlib

namespace LonelyRunner
namespace LRC14

/-- (L1) THE SORTED-GAP PIGEONHOLE (the chunk).  Removing `N` open rational
intervals of total (clipped) length ≤ `B` from `[a, b]` leaves a closed
subinterval of width ≥ `(b − a − B) / (N + 1)`. -/
theorem sorted_gap_pigeonhole (a b B : ℚ) (hab : a ≤ b)
    (teeth : List (ℚ × ℚ)) (hB : 0 ≤ B)
    (hmass : (teeth.map fun t => max 0 (min t.2 b - max t.1 a)).sum ≤ B) :
    ∃ a' b' : ℚ, a ≤ a' ∧ b' ≤ b ∧
      (b - a - B) / (teeth.length + 1) ≤ b' - a' ∧
      ∀ t ∈ teeth, b' ≤ t.1 ∨ t.2 ≤ a' := by
  sorry

/-- (L3 target shape) the cluster gap lemma, arcSafe form.  Stated for a
list of ≤ 6 positive moduli; the width constant matches THM-955. -/
theorem cluster_gap (Bs : List ℤ) (hk : Bs.length ≤ 6)
    (hpos : ∀ w ∈ Bs, 0 < w) (a b : ℚ) (hab : a ≤ b) :
    ∃ a' b' : ℚ, a ≤ a' ∧ b' ≤ b ∧
      ((1 - (Bs.length : ℚ)/7) * (b - a)
        - (Bs.map fun w => (2 : ℚ)/(7 * w)).sum)
        / (1 + (Bs.map fun w => (w : ℚ) * (b - a) + 2).sum)
      ≤ b' - a' ∧
      ∀ w ∈ Bs, WitnessCert.arcSafe (1/14) w 0 a' b' := by
  sorry

end LRC14
end LonelyRunner
