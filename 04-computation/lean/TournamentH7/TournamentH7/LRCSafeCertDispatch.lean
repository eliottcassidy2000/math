/-
  TournamentH7.LRCSafeCertDispatch — the hB shapeOf dispatcher
  (mac-mini-2026-07-09-S65 cont.23).

  The spine from the 2380-certificate table to the `hB` node itself:

  * `safeSet_congr` — `safeSet` depends only on list MEMBERSHIP;
  * `safe_nil_toReal` — the empty small part has full measure (`k = 13` row);
  * `length_filter_split` — the small/large filter split of the 13 speeds;
  * `shapeOf_fst_mem_bounds` — reachable small speeds lie in `[1, 13]`;
  * `cap_le_of_canonical` — a NODUP SORTED list of small speeds with
    `length + k ≤ 13` floors at `capRat k` (match on length, cite the leaf
    trees, close the constant gap by `capRat_mono`);
  * `hB_certs` — the dispatcher: dedup + insertionSort canonicalize any
    reachable `(shapeOf v).1`, then `cap_le_of_canonical`;
  * `lrc14_from_momentfloor_certs` — the moment route with `hB` DISCHARGED:
    obligations are now exactly `hMoment` (THM-661), `hsmall`, `hpartA`.
-/
import Mathlib
import TournamentH7.LRCSafeCertLeafTrees

set_option maxHeartbeats 1600000

namespace LonelyRunner
namespace LRC14
namespace SafeCertDispatch

open DenseCovers MeasureTheory SafeCert MomentFloor

/-- `safeSet` depends only on membership. -/
theorem safeSet_congr {P Q : List ℤ} (h : ∀ p, p ∈ P ↔ p ∈ Q) :
    safeSet P = safeSet Q := by
  ext x
  constructor <;> intro hx p hp
  · exact hx p ((h p).mpr hp)
  · exact hx p ((h p).mp hp)

/-- The `k = 13` row: an empty small part has full safe measure. -/
theorem safe_nil_toReal : (slowμ (safeSet ([] : List ℤ))).toReal = 1 := by
  have huniv : safeSet ([] : List ℤ) = Set.univ := by
    apply Set.eq_univ_of_forall
    intro x p hp
    exact absurd hp (List.not_mem_nil (a := p))
  rw [huniv]
  simp

/-- The small/large filter split: the two predicate filters partition the speeds. -/
theorem length_filter_split (l : List ℤ) :
    (l.filter (fun a => a ≤ 13)).length + (l.filter (fun a => (13 : ℤ) < a)).length
      = l.length := by
  induction l with
  | nil => rfl
  | cons x t ih =>
      by_cases h : x ≤ (13 : ℤ)
      · have h2 : ¬ ((13 : ℤ) < x) := not_lt.mpr h
        simp [h, h2]
        omega
      · have h2 : (13 : ℤ) < x := not_le.mp h
        simp [h, h2]
        omega

/-- The small part and the cluster split the 13 speeds. -/
theorem length_split_shapeOf (v : Fin 13 → ℤ) :
    (shapeOf v).1.length + clusterSize (shapeOf v) = 13 := by
  simp only [shapeOf, clusterSize, List.length_map]
  rw [length_filter_split, List.length_ofFn]

/-- Reachable small speeds lie in `[1, 13]`. -/
theorem shapeOf_fst_mem_bounds (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) :
    ∀ p ∈ (shapeOf v).1, 1 ≤ p ∧ p ≤ 13 := by
  intro p hp
  have hp' : p ∈ (List.ofFn fun i => |v i|).filter (fun a => a ≤ 13) := hp
  obtain ⟨hmem, hle⟩ := List.mem_filter.mp hp'
  have h13 : p ≤ 13 := of_decide_eq_true hle
  obtain ⟨i, hi⟩ := List.mem_ofFn.mp hmem
  have hpos : 0 < |v i| := abs_pos.mpr (hv i)
  constructor
  · omega
  · exact h13

/-- `capRat` is monotone on the binding range (finite check). -/
theorem capRat_mono {k k' : ℕ} (h8 : 8 ≤ k) (hkk : k ≤ k') (h13 : k' ≤ 13) :
    (capRat k : ℝ) ≤ (capRat k' : ℝ) := by
  have hk13 : k ≤ 13 := le_trans hkk h13
  interval_cases k <;> interval_cases k' <;> norm_num [capRat]

/-- **The canonical-list floor.**  A strictly sorted list of speeds in `[1, 13]` whose
length leaves room for a cluster of size `k ∈ [8, 13]` floors at `capRat k`. -/
theorem cap_le_of_canonical (T : List ℤ) (hlt : T.Pairwise (· < ·))
    (hbounds : ∀ p ∈ T, 1 ≤ p ∧ p ≤ 13) {k : ℕ} (h8 : 8 ≤ k) (h13 : k ≤ 13)
    (hlen : T.length + k ≤ 13) :
    (capRat k : ℝ) ≤ (slowμ (safeSet T)).toReal := by
  match T, hlt, hbounds, hlen with
  | [], _, _, _ =>
      rw [safe_nil_toReal]
      interval_cases k <;> norm_num [capRat]
  | [p1], _, hbounds, hlen =>
      have hb1 := hbounds p1 (by simp)
      have hfloor := safe_floor_sorted_len1 p1 hb1.1 hb1.2
      have hk : k ≤ 12 := by
        simp only [List.length_cons, List.length_nil] at hlen; omega
      have hcap : (capRat k : ℝ) ≤ ((6 : ℝ)/7) := by
        have hc : (capRat 12 : ℝ) = ((6 : ℝ)/7) := by norm_num [capRat]
        rw [← hc]
        exact capRat_mono h8 hk (by norm_num)
      linarith
  | [p1, p2], hlt, hbounds, hlen =>
      have hb1 := hbounds p1 (by simp)
      have hb2 := hbounds p2 (by simp)
      obtain ⟨ha1, _⟩ := List.pairwise_cons.mp hlt
      have h12 : p1 < p2 := ha1 p2 (by simp)
      have hfloor := safe_floor_sorted_len2 p1 p2 hb1.1 h12 hb2.2
      have hk : k ≤ 11 := by
        simp only [List.length_cons, List.length_nil] at hlen; omega
      have hcap : (capRat k : ℝ) ≤ ((66 : ℝ)/91) := by
        have hc : (capRat 11 : ℝ) = ((66 : ℝ)/91) := by norm_num [capRat]
        rw [← hc]
        exact capRat_mono h8 hk (by norm_num)
      linarith
  | [p1, p2, p3], hlt, hbounds, hlen =>
      have hb1 := hbounds p1 (by simp)
      have hb3 := hbounds p3 (by simp)
      obtain ⟨ha1, hlt2⟩ := List.pairwise_cons.mp hlt
      obtain ⟨ha2, _⟩ := List.pairwise_cons.mp hlt2
      have h12 : p1 < p2 := ha1 p2 (by simp)
      have h23 : p2 < p3 := ha2 p3 (by simp)
      have hfloor := safe_floor_sorted_len3 p1 p2 p3 hb1.1 h12 h23 hb3.2
      have hk : k ≤ 10 := by
        simp only [List.length_cons, List.length_nil] at hlen; omega
      have hcap : (capRat k : ℝ) ≤ ((55 : ℝ)/91) := by
        have hc : (capRat 10 : ℝ) = ((55 : ℝ)/91) := by norm_num [capRat]
        rw [← hc]
        exact capRat_mono h8 hk (by norm_num)
      linarith
  | [p1, p2, p3, p4], hlt, hbounds, hlen =>
      have hb1 := hbounds p1 (by simp)
      have hb4 := hbounds p4 (by simp)
      obtain ⟨ha1, hlt2⟩ := List.pairwise_cons.mp hlt
      obtain ⟨ha2, hlt3⟩ := List.pairwise_cons.mp hlt2
      obtain ⟨ha3, _⟩ := List.pairwise_cons.mp hlt3
      have h12 : p1 < p2 := ha1 p2 (by simp)
      have h23 : p2 < p3 := ha2 p3 (by simp)
      have h34 : p3 < p4 := ha3 p4 (by simp)
      have hfloor := safe_floor_sorted_len4 p1 p2 p3 p4 hb1.1 h12 h23 h34 hb4.2
      have hk : k ≤ 9 := by
        simp only [List.length_cons, List.length_nil] at hlen; omega
      have hcap : (capRat k : ℝ) ≤ ((1979 : ℝ)/4004) := by
        have hc : (capRat 9 : ℝ) = ((1979 : ℝ)/4004) := by norm_num [capRat]
        rw [← hc]
        exact capRat_mono h8 hk (by norm_num)
      linarith
  | [p1, p2, p3, p4, p5], hlt, hbounds, hlen =>
      have hb1 := hbounds p1 (by simp)
      have hb5 := hbounds p5 (by simp)
      obtain ⟨ha1, hlt2⟩ := List.pairwise_cons.mp hlt
      obtain ⟨ha2, hlt3⟩ := List.pairwise_cons.mp hlt2
      obtain ⟨ha3, hlt4⟩ := List.pairwise_cons.mp hlt3
      obtain ⟨ha4, _⟩ := List.pairwise_cons.mp hlt4
      have h12 : p1 < p2 := ha1 p2 (by simp)
      have h23 : p2 < p3 := ha2 p3 (by simp)
      have h34 : p3 < p4 := ha3 p4 (by simp)
      have h45 : p4 < p5 := ha4 p5 (by simp)
      have hfloor := safe_floor_sorted_len5 p1 p2 p3 p4 p5 hb1.1 h12 h23 h34 h45 hb5.2
      have hk : k ≤ 8 := by
        simp only [List.length_cons, List.length_nil] at hlen; omega
      have hcap : (capRat k : ℝ) ≤ ((2243 : ℝ)/5880) := by
        have hc : (capRat 8 : ℝ) = ((2243 : ℝ)/5880) := by norm_num [capRat]
        rw [← hc]
        exact capRat_mono h8 hk (by norm_num)
      linarith
  | p1 :: p2 :: p3 :: p4 :: p5 :: p6 :: rest, _, _, hlen =>
      exfalso
      simp only [List.length_cons] at hlen
      omega

/-- **The hB dispatcher.**  Every reachable shape with `8 ≤ k ≤ 13` floors its small-part
safe measure at `capRat k`: dedup + insertionSort canonicalize the small speeds, then the
leaf trees dispatch to the certificate table. -/
theorem hB_certs (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (h8 : 8 ≤ clusterSize (shapeOf v)) (h13 : clusterSize (shapeOf v) ≤ 13) :
    (capRat (clusterSize (shapeOf v)) : ℝ) ≤ MomentFloor.measGPConcrete (shapeOf v) := by
  have hbounds := shapeOf_fst_mem_bounds v hv
  have hsplit := length_split_shapeOf v
  have hperm := List.perm_insertionSort (· ≤ ·) ((shapeOf v).1.dedup)
  have hmemT : ∀ p, p ∈ List.insertionSort (· ≤ ·) ((shapeOf v).1.dedup)
      ↔ p ∈ (shapeOf v).1 := by
    intro p
    rw [hperm.mem_iff, List.mem_dedup]
  have hcongr : safeSet (shapeOf v).1
      = safeSet (List.insertionSort (· ≤ ·) ((shapeOf v).1.dedup)) :=
    safeSet_congr (fun p => (hmemT p).symm)
  have hlenT : (List.insertionSort (· ≤ ·) ((shapeOf v).1.dedup)).length
      + clusterSize (shapeOf v) ≤ 13 := by
    have h1 : (List.insertionSort (· ≤ ·) ((shapeOf v).1.dedup)).length
        = (shapeOf v).1.dedup.length := hperm.length_eq
    have h2 : (shapeOf v).1.dedup.length ≤ (shapeOf v).1.length :=
      (List.dedup_sublist (shapeOf v).1).length_le
    omega
  have hnodupT : (List.insertionSort (· ≤ ·) ((shapeOf v).1.dedup)).Nodup :=
    hperm.nodup_iff.mpr (List.nodup_dedup (shapeOf v).1)
  have hsortT : (List.insertionSort (· ≤ ·) ((shapeOf v).1.dedup)).Pairwise (· ≤ ·) :=
    List.pairwise_insertionSort _ _
  have hltT : (List.insertionSort (· ≤ ·) ((shapeOf v).1.dedup)).Pairwise (· < ·) :=
    (hsortT.and hnodupT).imp (fun h => lt_of_le_of_ne h.1 h.2)
  have hboundsT : ∀ p ∈ List.insertionSort (· ≤ ·) ((shapeOf v).1.dedup),
      1 ≤ p ∧ p ≤ 13 :=
    fun p hp => hbounds p ((hmemT p).mp hp)
  have hkey := cap_le_of_canonical _ hltT hboundsT h8 h13 hlenT
  show (capRat (clusterSize (shapeOf v)) : ℝ) ≤ (slowμ (safeSet (shapeOf v).1)).toReal
  rw [hcongr]
  exact hkey

/-- **The moment route with `hB` DISCHARGED.**  LRC(14) from `hMoment` (THM-661),
`hsmall` (`k ≤ 7` pigeonhole), and `hpartA` (the reach) — the certificate table plus
this dispatcher supply Lemma B in full. -/
theorem lrc14_from_momentfloor_certs
    (hMoment : ∀ s : Shape, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (momentBar (clusterSize s) : ℝ) ≤ MomentFloor.nuShapeConcrete s)
    (hsmall : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 7 → (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  MomentFloorRepair.lrc14_from_momentfloor_concrete_shapes hMoment
    (fun v hv h8 h13 => hB_certs v hv h8 h13) hsmall hpartA

end SafeCertDispatch
end LRC14
end LonelyRunner
