/-
  TournamentH7.LRCSmallDischarge — hsmall DISCHARGED via the repaired legs
  (mac-mini-2026-07-09-S65 cont.25).

  The moment route's `hsmall` node (`k ≤ 7 ⟹ m_P ≤ witnessG2`) is UNSATISFIABLE at
  `k ≤ 2` (cont.16 finding, MISTAKE-136 genus).  This file discharges the SATISFIABLE
  content through the cont.17 repaired assembly instead:

  * `hk12` (k ∈ {1,2}, POSITIVITY): the cluster has ≤ 2 teeth so `goodSet = univ`
    (brick (i)); `witnessG2 = μ(safeSet P)`; the small part misses some `m ∈ [1,13]`
    (pigeonhole: ≤ 12 values), so it contains a twelve-element anchor family whose
    anchor interval is positive (`anchor_pos_of_missing`, 13 certificates).
  * `hsmall3` (3 ≤ k ≤ 7, m_P floor): `ν = 1` (brick (i) / the pigeonhole citation on
    the dedup) + Bonferroni + the m_P certificate layer: `μ(safeSet P) ≥ m_P` for every
    `|P| ≤ 10` — sizes ≤ 5 via the capRat table (`cap_le_of_canonical` at `k = 8`,
    `m_P ≤ capRat 8`), sizes 6..10 via the new completion trees (`mp_floor_sorted_len*`,
    5720 leaves; m_P is the EXACT min over `|S| ≤ 10`, attained at
    `{1,2,3,5,7,8,9,11,12,13}`).
  * `hlarge` (8 ≤ k ≤ 13): Bonferroni + `hMoment_of_citations` (cont.24) + `hB_certs`
    (cont.23) + the exact identity `momentBar k + capRat k − 1 = witnessMP`.

  **`lrc14_from_two_citations`: LRC(14) from THM-661 + the ≤7-arcs pigeonhole + the
  reach node `hpartA`.  Everything else is a theorem.**
-/
import Mathlib
import TournamentH7.LRCMPLeafTrees
import TournamentH7.LRCWitnessFloorRepair

set_option maxHeartbeats 1600000

namespace LonelyRunner
namespace LRC14
namespace SmallDischarge

open DenseCovers MeasureTheory MomentFloor MomentFloorRepair MomentCitation MPCert
  SafeCertDispatch TournamentH7.GoodSet

/-- `ν = 1` for any nonempty cluster list of length ≤ 7: dedup, then brick (i)
(≤ 2 teeth) or the pigeonhole citation (3..7 teeth). -/
theorem nu_eq_one_of_le7 (hsmall7 : SmallClusterFull) (E : List ℤ)
    (hne : E ≠ []) (hlen : E.length ≤ 7) :
    (slowμ (goodSet E)).toReal = 1 := by
  have hgd : goodSet E = goodSet E.dedup := (goodSet_dedup E).symm
  rw [hgd]
  have hnd := List.nodup_dedup E
  have hlend : E.dedup.length ≤ 7 := le_trans (List.dedup_sublist E).length_le hlen
  have hned : E.dedup ≠ [] := fun h => hne ((List.dedup_eq_nil E).mp h)
  rcases Nat.lt_or_ge E.dedup.length 3 with hsm | h3
  · exact nu_of_univ (GoodSetSmall.goodSet_univ_of_length_le_two E.dedup hned (by omega))
  · exact hsmall7 E.dedup hnd h3 hlend

/-- **The m_P canonical floor.**  A strictly sorted list of ≤ 10 speeds in `[1, 13]`
floors the safe measure at `m_P = 14249/252252` (sizes ≤ 5 via the capRat table at
`k = 8`; sizes 6..10 via the completion trees). -/
theorem mp_le_of_canonical (T : List ℤ) (hlt : T.Pairwise (· < ·))
    (hbounds : ∀ p ∈ T, 1 ≤ p ∧ p ≤ 13) (hlen : T.length ≤ 10) :
    ((14249 : ℝ)/252252) ≤ (slowμ (safeSet T)).toReal := by
  rcases Nat.lt_or_ge T.length 6 with h5 | h6
  · have h := cap_le_of_canonical T hlt hbounds (k := 8) (by norm_num) (by norm_num)
      (by omega)
    have hcap : ((14249 : ℝ)/252252) ≤ (capRat 8 : ℝ) := by norm_num [capRat]
    linarith
  · match T, hlt, hbounds, h6, hlen with
    | [p1, p2, p3, p4, p5, p6], hlt, hbounds, _, _ =>
        have hb1 := hbounds p1 (by simp)
        have hbt := hbounds p6 (by simp)
        obtain ⟨ha1, hlt2⟩ := List.pairwise_cons.mp hlt
        obtain ⟨ha2, hlt3⟩ := List.pairwise_cons.mp hlt2
        obtain ⟨ha3, hlt4⟩ := List.pairwise_cons.mp hlt3
        obtain ⟨ha4, hlt5⟩ := List.pairwise_cons.mp hlt4
        obtain ⟨ha5, _⟩ := List.pairwise_cons.mp hlt5
        exact mp_floor_sorted_len6 p1 p2 p3 p4 p5 p6 hb1.1
          (ha1 p2 (by simp)) (ha2 p3 (by simp)) (ha3 p4 (by simp))
          (ha4 p5 (by simp)) (ha5 p6 (by simp)) hbt.2
    | [p1, p2, p3, p4, p5, p6, p7], hlt, hbounds, _, _ =>
        have hb1 := hbounds p1 (by simp)
        have hbt := hbounds p7 (by simp)
        obtain ⟨ha1, hlt2⟩ := List.pairwise_cons.mp hlt
        obtain ⟨ha2, hlt3⟩ := List.pairwise_cons.mp hlt2
        obtain ⟨ha3, hlt4⟩ := List.pairwise_cons.mp hlt3
        obtain ⟨ha4, hlt5⟩ := List.pairwise_cons.mp hlt4
        obtain ⟨ha5, hlt6⟩ := List.pairwise_cons.mp hlt5
        obtain ⟨ha6, _⟩ := List.pairwise_cons.mp hlt6
        exact mp_floor_sorted_len7 p1 p2 p3 p4 p5 p6 p7 hb1.1
          (ha1 p2 (by simp)) (ha2 p3 (by simp)) (ha3 p4 (by simp))
          (ha4 p5 (by simp)) (ha5 p6 (by simp)) (ha6 p7 (by simp)) hbt.2
    | [p1, p2, p3, p4, p5, p6, p7, p8], hlt, hbounds, _, _ =>
        have hb1 := hbounds p1 (by simp)
        have hbt := hbounds p8 (by simp)
        obtain ⟨ha1, hlt2⟩ := List.pairwise_cons.mp hlt
        obtain ⟨ha2, hlt3⟩ := List.pairwise_cons.mp hlt2
        obtain ⟨ha3, hlt4⟩ := List.pairwise_cons.mp hlt3
        obtain ⟨ha4, hlt5⟩ := List.pairwise_cons.mp hlt4
        obtain ⟨ha5, hlt6⟩ := List.pairwise_cons.mp hlt5
        obtain ⟨ha6, hlt7⟩ := List.pairwise_cons.mp hlt6
        obtain ⟨ha7, _⟩ := List.pairwise_cons.mp hlt7
        exact mp_floor_sorted_len8 p1 p2 p3 p4 p5 p6 p7 p8 hb1.1
          (ha1 p2 (by simp)) (ha2 p3 (by simp)) (ha3 p4 (by simp))
          (ha4 p5 (by simp)) (ha5 p6 (by simp)) (ha6 p7 (by simp))
          (ha7 p8 (by simp)) hbt.2
    | [p1, p2, p3, p4, p5, p6, p7, p8, p9], hlt, hbounds, _, _ =>
        have hb1 := hbounds p1 (by simp)
        have hbt := hbounds p9 (by simp)
        obtain ⟨ha1, hlt2⟩ := List.pairwise_cons.mp hlt
        obtain ⟨ha2, hlt3⟩ := List.pairwise_cons.mp hlt2
        obtain ⟨ha3, hlt4⟩ := List.pairwise_cons.mp hlt3
        obtain ⟨ha4, hlt5⟩ := List.pairwise_cons.mp hlt4
        obtain ⟨ha5, hlt6⟩ := List.pairwise_cons.mp hlt5
        obtain ⟨ha6, hlt7⟩ := List.pairwise_cons.mp hlt6
        obtain ⟨ha7, hlt8⟩ := List.pairwise_cons.mp hlt7
        obtain ⟨ha8, _⟩ := List.pairwise_cons.mp hlt8
        exact mp_floor_sorted_len9 p1 p2 p3 p4 p5 p6 p7 p8 p9 hb1.1
          (ha1 p2 (by simp)) (ha2 p3 (by simp)) (ha3 p4 (by simp))
          (ha4 p5 (by simp)) (ha5 p6 (by simp)) (ha6 p7 (by simp))
          (ha7 p8 (by simp)) (ha8 p9 (by simp)) hbt.2
    | [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10], hlt, hbounds, _, _ =>
        have hb1 := hbounds p1 (by simp)
        have hbt := hbounds p10 (by simp)
        obtain ⟨ha1, hlt2⟩ := List.pairwise_cons.mp hlt
        obtain ⟨ha2, hlt3⟩ := List.pairwise_cons.mp hlt2
        obtain ⟨ha3, hlt4⟩ := List.pairwise_cons.mp hlt3
        obtain ⟨ha4, hlt5⟩ := List.pairwise_cons.mp hlt4
        obtain ⟨ha5, hlt6⟩ := List.pairwise_cons.mp hlt5
        obtain ⟨ha6, hlt7⟩ := List.pairwise_cons.mp hlt6
        obtain ⟨ha7, hlt8⟩ := List.pairwise_cons.mp hlt7
        obtain ⟨ha8, hlt9⟩ := List.pairwise_cons.mp hlt8
        obtain ⟨ha9, _⟩ := List.pairwise_cons.mp hlt9
        exact mp_floor_sorted_len10 p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 hb1.1
          (ha1 p2 (by simp)) (ha2 p3 (by simp)) (ha3 p4 (by simp))
          (ha4 p5 (by simp)) (ha5 p6 (by simp)) (ha6 p7 (by simp))
          (ha7 p8 (by simp)) (ha8 p9 (by simp)) (ha9 p10 (by simp)) hbt.2
    | p1 :: p2 :: p3 :: p4 :: p5 :: p6 :: p7 :: p8 :: p9 :: p10 :: p11 :: rest,
        _, _, _, hlen =>
        exfalso
        simp only [List.length_cons] at hlen
        omega

/-- **The m_P bounded floor**: ANY speed list with values in `[1,13]` and length ≤ 10
floors at `m_P` (dedup + insertionSort canonicalize, then the canonical floor). -/
theorem mp_le_of_bounded (P : List ℤ) (hbounds : ∀ p ∈ P, 1 ≤ p ∧ p ≤ 13)
    (hlen : P.length ≤ 10) :
    ((14249 : ℝ)/252252) ≤ (slowμ (safeSet P)).toReal := by
  have hperm := List.perm_insertionSort (· ≤ ·) (P.dedup)
  have hmemT : ∀ p, p ∈ List.insertionSort (· ≤ ·) (P.dedup) ↔ p ∈ P := by
    intro p
    rw [hperm.mem_iff, List.mem_dedup]
  have hcongr : safeSet P = safeSet (List.insertionSort (· ≤ ·) (P.dedup)) :=
    safeSet_congr (fun p => (hmemT p).symm)
  have hlenT : (List.insertionSort (· ≤ ·) (P.dedup)).length ≤ 10 := by
    have h1 : (List.insertionSort (· ≤ ·) (P.dedup)).length = P.dedup.length :=
      hperm.length_eq
    have h2 : P.dedup.length ≤ P.length := (List.dedup_sublist P).length_le
    omega
  have hnodupT : (List.insertionSort (· ≤ ·) (P.dedup)).Nodup :=
    hperm.nodup_iff.mpr (List.nodup_dedup P)
  have hsortT : (List.insertionSort (· ≤ ·) (P.dedup)).Pairwise (· ≤ ·) :=
    List.pairwise_insertionSort _ _
  have hltT : (List.insertionSort (· ≤ ·) (P.dedup)).Pairwise (· < ·) :=
    (hsortT.and hnodupT).imp (fun h => lt_of_le_of_ne h.1 h.2)
  have hboundsT : ∀ p ∈ List.insertionSort (· ≤ ·) (P.dedup), 1 ≤ p ∧ p ≤ 13 :=
    fun p hp => hbounds p ((hmemT p).mp hp)
  rw [hcongr]
  exact mp_le_of_canonical _ hltT hboundsT hlenT

/-- **`hk12` DISCHARGED**: positivity at `k ∈ {1, 2}`. -/
theorem hk12_discharged : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
    1 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 2 →
    0 < witnessG2 (shapeOf v) := by
  intro v hv h1 h2
  have hEne : (shapeOf v).2 ≠ [] := by
    apply List.ne_nil_of_length_pos
    show 0 < (shapeOf v).2.length
    exact h1
  have huniv := GoodSetSmall.goodSet_univ_of_length_le_two (shapeOf v).2 hEne h2
  have hw : witnessG2 (shapeOf v) = (slowμ (safeSet (shapeOf v).1)).toReal := by
    unfold witnessG2
    rw [huniv, Set.univ_inter]
  rw [hw]
  have hbounds := shapeOf_fst_mem_bounds v hv
  have hsplit := length_split_shapeOf v
  have hex : ∃ m ∈ Finset.Icc (1 : ℤ) 13, m ∉ (shapeOf v).1.toFinset := by
    by_contra hcon
    push_neg at hcon
    have hsub : Finset.Icc (1 : ℤ) 13 ⊆ (shapeOf v).1.toFinset :=
      fun m hm => hcon m hm
    have hle := Finset.card_le_card hsub
    have hicc : (Finset.Icc (1 : ℤ) 13).card = 13 := by decide
    have hcard := List.toFinset_card_le (shapeOf v).1
    omega
  obtain ⟨m, hmIcc, hmnot⟩ := hex
  rw [Finset.mem_Icc] at hmIcc
  exact anchor_pos_of_missing (shapeOf v).1 hbounds m hmIcc.1 hmIcc.2
    (fun h => hmnot (List.mem_toFinset.mpr h))

/-- **`hsmall3` DISCHARGED** (given the pigeonhole citation): the m_P floor at
`3 ≤ k ≤ 7`, via `ν = 1` + Bonferroni + the m_P certificate layer. -/
theorem hsmall3_discharged (hsmall7 : SmallClusterFull) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
    3 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 7 →
    (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) := by
  intro v hv h3 h7
  have hbonf := bonferroni_concrete (shapeOf v)
  have hEne : (shapeOf v).2 ≠ [] := by
    apply List.ne_nil_of_length_pos
    show 0 < (shapeOf v).2.length
    exact Nat.lt_of_lt_of_le (by norm_num) h3
  have hν : nuShapeConcrete (shapeOf v) = 1 := by
    show (slowμ (goodSet (shapeOf v).2)).toReal = 1
    exact nu_eq_one_of_le7 hsmall7 (shapeOf v).2 hEne h7
  have hP : ((14249 : ℝ)/252252) ≤ measGPConcrete (shapeOf v) := by
    show ((14249 : ℝ)/252252) ≤ (slowμ (safeSet (shapeOf v).1)).toReal
    apply mp_le_of_bounded
    · exact shapeOf_fst_mem_bounds v hv
    · have := length_split_shapeOf v
      omega
  have hwmp : (witnessMP : ℝ) = ((14249 : ℝ)/252252) := by
    norm_num [witnessMP]
  rw [hwmp]
  rw [hν] at hbonf
  linarith

/-- **`hlarge` DISCHARGED** (given the two citations): the m_P floor at `8 ≤ k ≤ 13`,
via Bonferroni + the moment citation + the capRat certificates + the exact identity. -/
theorem hlarge_discharged (h661 : THM661MomentFloor) (hsmall7 : SmallClusterFull) :
    ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
    8 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 13 →
    (witnessMP : ℝ) ≤ witnessG2 (shapeOf v) := by
  intro v hv h8 h13
  have hbonf := bonferroni_concrete (shapeOf v)
  have hM := hMoment_of_citations h661 hsmall7 (shapeOf v) h8 h13
  have hB := hB_certs v hv h8 h13
  have hid := momentBar_add_capRat (clusterSize (shapeOf v))
  have hcast : (momentBar (clusterSize (shapeOf v)) : ℝ)
      + (capRat (clusterSize (shapeOf v)) : ℝ) - 1 = (witnessMP : ℝ) := by
    have := congrArg (fun q : ℚ => (q : ℝ)) hid
    push_cast at this ⊢
    linarith
  linarith

/-- **THE MOMENT ROUTE, CLOSED.**  LRC(14) from exactly TWO NAMED CITATIONS — THM-661
(the unified covering-moment density floor) and the ≤7-arcs pigeonhole — plus the reach
node `hpartA`.  Every other node is a THEOREM: `hk12`/`hsmall3`/`hlarge` above (the
`k = 0` sieve is internal to the repaired assembly), `hB` = cont.23's dispatcher,
`hMoment` = cont.24's bridge. -/
theorem lrc14_from_two_citations (h661 : THM661MomentFloor) (hsmall7 : SmallClusterFull)
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  Repair.lrc14_from_repaired_nodes hk12_discharged (hsmall3_discharged hsmall7)
    (hlarge_discharged h661 hsmall7) hpartA

end SmallDischarge
end LRC14
end LonelyRunner
