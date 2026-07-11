/-
  TournamentH7.LRCMomentCitation — the THM-661 citation node for `hMoment`
  (mac-mini-2026-07-09-S65 cont.24).

  THM-661 (canon, PROVED classically: mac-mini-S57/S58 + opus-S148 + the LEM-009 fleet):
  the unified covering-moment density floor.  On the covering reformulation
  `ν(E) = μ(goodSet E) = P(k arcs of length 1/7 fail to cover)`, the degree-≤4 moment
  bound clears the honest (A′) bars for EVERY admissible cluster, diameter-free:

  > for every set of k DISTINCT integer co-offsets `E`, `8 ≤ k ≤ 13`:
  >   `ν(E) ≥ bar_k = m_P + 1 − min_{|P| = 13−k} μ(safeSet P) = momentBar k`.

  Per project policy (like LRC(≤13)), this PROVED-in-canon analytic theorem enters Lean
  as a NAMED CITATION HYPOTHESIS (`THM661MomentFloor`), not a sorry.  Its companion
  (`SmallClusterFull`, the standard ≤7-arcs pigeonhole `ν = 1` — total arc length
  `k/7 ≤ 1` cannot cover the circle off a measure-zero set) covers dedup'd clusters
  that fall below 8 distinct teeth; lengths ≤ 2 need NO citation (brick (i),
  `goodSet_univ_of_length_le_two`, cont.19 — proved).

  TRANSCRIPTION VALIDATED (lrc14_nu_probe_macmini_S65cont24): exact `momentBar k`
  reproduces canon's six bars (0.6750/0.5622/0.4521/0.3312/0.1993/0.0565) and exact
  `ν(block_k)` reproduces canon's table (691/735, 247/294, 38/49, 1381/2205,
  13823/24255, 477/1078) to the rational.

  The bridge `hMoment_of_citations` dispatches an ARBITRARY reachable cluster list
  (duplicates allowed — `goodSet` only sees `toFinset`) through its dedup:
  ℓ ≤ 2 proved / 3 ≤ ℓ ≤ 7 pigeonhole / 8 ≤ ℓ ≤ 13 THM-661 + `momentBar` antitone.
  `lrc14_from_thm661_certs` then closes the moment route's Lean surface at exactly
  TWO citations + `hsmall` + `hpartA`.
-/
import Mathlib
import TournamentH7.LRCSafeCertDispatch
import TournamentH7.LRCGoodSetSmall

set_option maxHeartbeats 800000

namespace LonelyRunner
namespace LRC14
namespace MomentCitation

open DenseCovers MeasureTheory MomentFloor TournamentH7.GoodSet

/-- **THM-661 (citation).**  The unified covering-moment density floor: every cluster of
`k ∈ [8, 13]` DISTINCT integer co-offsets has good-set measure at least `momentBar k`.
PROVED in canon (THM-661 + its S58 uniform-floor addenda); enters Lean as a named
citation hypothesis per project policy. -/
def THM661MomentFloor : Prop :=
  ∀ E : List ℤ, E.Nodup → 8 ≤ E.length → E.length ≤ 13 →
    (momentBar E.length : ℝ) ≤ (slowμ (goodSet E)).toReal

/-- **The small-cluster pigeonhole (citation).**  With `3 ≤ k ≤ 7` distinct teeth the
`k` following-arcs of length `1/7` cover total length `≤ 1` and fail to cover the circle
off a measure-zero set of times: `ν = 1`.  (Lengths ≤ 2 are PROVED — brick (i).) -/
def SmallClusterFull : Prop :=
  ∀ E : List ℤ, E.Nodup → 3 ≤ E.length → E.length ≤ 7 →
    (slowμ (goodSet E)).toReal = 1

/-- `goodSet` sees only the `toFinset`, so dedup is invisible. -/
theorem goodSet_dedup (E : List ℤ) : goodSet E.dedup = goodSet E := by
  have h : E.dedup.toFinset = E.toFinset := by
    ext a
    simp [List.mem_toFinset, List.mem_dedup]
  unfold goodSet
  rw [h]

/-- The bars are genuine probabilities on the binding range (finite check). -/
theorem momentBar_le_one {k : ℕ} (h8 : 8 ≤ k) (h13 : k ≤ 13) :
    (momentBar k : ℝ) ≤ 1 := by
  interval_cases k <;> norm_num [momentBar, capRat, witnessMP]

/-- `momentBar` is ANTITONE on the binding range (finite check): fewer distinct teeth,
higher bar. -/
theorem momentBar_anti {l k : ℕ} (h8 : 8 ≤ l) (hlk : l ≤ k) (h13 : k ≤ 13) :
    (momentBar k : ℝ) ≤ (momentBar l : ℝ) := by
  have hl13 : l ≤ 13 := le_trans hlk h13
  interval_cases l <;> interval_cases k <;> norm_num [momentBar, capRat, witnessMP]

/-- `ν = 1` when the good set is everything. -/
theorem nu_of_univ {E : List ℤ} (h : goodSet E = Set.univ) :
    (slowμ (goodSet E)).toReal = 1 := by
  rw [h]
  simp

/-- **The bridge: `hMoment` from the two citations.**  An arbitrary reachable cluster
list dispatches through its dedup: `ℓ ≤ 2` is proved (brick (i)), `3 ≤ ℓ ≤ 7` is the
pigeonhole, `8 ≤ ℓ ≤ 13` is THM-661 plus `momentBar` antitonicity. -/
theorem hMoment_of_citations (h661 : THM661MomentFloor) (hsmall7 : SmallClusterFull) :
    ∀ s : Shape, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (momentBar (clusterSize s) : ℝ) ≤ nuShapeConcrete s := by
  intro s h8 h13
  have hν : nuShapeConcrete s = (slowμ (goodSet s.2.dedup)).toReal := by
    unfold nuShapeConcrete
    rw [goodSet_dedup]
  rw [hν]
  have hnd : s.2.dedup.Nodup := List.nodup_dedup s.2
  have hlen : s.2.dedup.length ≤ clusterSize s :=
    (List.dedup_sublist s.2).length_le
  have hne : s.2.dedup ≠ [] := by
    intro h
    have h0 : s.2 = [] := (List.dedup_eq_nil s.2).mp h
    have h8' : 8 ≤ s.2.length := h8
    rw [h0] at h8'
    simp at h8'
  rcases Nat.lt_or_ge s.2.dedup.length 3 with hsm | h3
  · -- ℓ ≤ 2: brick (i), proved
    have huniv := GoodSetSmall.goodSet_univ_of_length_le_two s.2.dedup hne (by omega)
    rw [nu_of_univ huniv]
    exact momentBar_le_one h8 h13
  rcases Nat.lt_or_ge s.2.dedup.length 8 with hmid | h8'
  · -- 3 ≤ ℓ ≤ 7: the pigeonhole citation
    rw [hsmall7 s.2.dedup hnd h3 (by omega)]
    exact momentBar_le_one h8 h13
  · -- 8 ≤ ℓ ≤ 13: THM-661 + antitonicity
    have hfloor := h661 s.2.dedup hnd h8' (le_trans hlen h13)
    have hanti := momentBar_anti h8' hlen h13
    linarith

/-- **The moment route closed at the citations.**  LRC(14) from: THM-661 (cited),
the small-cluster pigeonhole (cited), `hsmall` (`k ≤ 7` m_P floor), and `hpartA`
(the reach).  Lemma B is supplied by the cont.22/23 certificate table + dispatcher. -/
theorem lrc14_from_thm661_certs
    (h661 : THM661MomentFloor) (hsmall7 : SmallClusterFull)
    (hsmall : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 7 → (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  SafeCertDispatch.lrc14_from_momentfloor_certs
    (hMoment_of_citations h661 hsmall7) hsmall hpartA

end MomentCitation
end LRC14
end LonelyRunner
