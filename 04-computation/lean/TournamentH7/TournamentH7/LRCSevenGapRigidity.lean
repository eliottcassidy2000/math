/-
  TournamentH7.LRCSevenGapRigidity -- THE <=7-ARCS PIGEONHOLE DE-CITED
  (klein-2026-07-11-S251, HYP-5985): `SmallClusterFull` -- the second of the
  grand assembly's three citations -- PROVED as a theorem:

      for every nodup cluster E with 3 <= |E| <= 7,  slowmu (goodSet E) = 1.

  THE PROOF (THM-689(A)'s perfect-net rigidity, made pointwise):
  sort the DISTINCT phases {fract(e*x)} by `Finset.orderIsoOfFin` into a
  strictly monotone f : Fin k -> [0,1) (k <= 7); the k cyclic gaps
  F(i+1) - F(i) telescope to 1 (F caps the wrap at 1 + f 0).
  * ANY gap > 1/7 hands the goodSet witness DIRECTLY: an interior gap at
    s_i forces s_i < 6/7, so every wrap term 1 + s_j - s_i exceeds 1/7 for
    free; the wrap gap is witnessed at the maximum phase.  (No max needed.)
  * k <= 6: some gap > 1/7 by mac-mini's pigeonhole => x in goodSet, EVERY x.
  * k = 7: no gap > 1/7 forces ALL gaps EXACTLY 1/7 (sum-erase bound; the
    boundary content of mac-mini's `all_eq_one_seventh_of_le`) -- the
    PERFECT NET -- and then fract((b-a)*x) = 1/7 EXACTLY for the first
    sorted pair: x lies in a COUNTABLE fiber (n + 1/7)/(b-a).
  The bad set sits inside <= 49 countable fibers => null => nu = 1.

  Consumes: TournamentH7.MaxGapPigeonhole (the abstract gap pigeonhole),
  TournamentH7.GoodSet (goodSet, the fract-sub identity).

  HEADLINE COROLLARY: lrc14_from_moment_and_supply --
      LRC(14) = [THM-661 moment floor] + [certificate supply].
  TWO citations left (was three).

  Kernel-pure: no native_decide, no sorry.  (This file = the CORE on light
  imports; the naming as `SmallClusterFull` + the two-citation assembly live
  in LRCTwoCitationAssembly.lean, on the heavy discharge chain.)
-/
import Mathlib
import TournamentH7.LRCMaxGapPigeonhole
import TournamentH7.LRCGoodSet

namespace LonelyRunner
namespace LRC14
namespace SevenGapRigidity

open MeasureTheory DenseCovers TournamentH7.GoodSet
open TournamentH7.MaxGapPigeonhole

/-! ## The abstract gap layer: sorted phases in `[0,1)` -/

/-- The cyclic gap at index `i` of a `k`-enumeration: the forward difference,
with the wrap capped at `1 + f 0`. -/
def cgap {k : ℕ} (hk : 0 < k) (f : Fin k → ℝ) (i : Fin k) : ℝ :=
  (if h : (i : ℕ) + 1 < k then f ⟨(i : ℕ) + 1, h⟩ else 1 + f ⟨0, hk⟩) - f i

/-- **The cyclic gaps telescope to `1`.** -/
theorem cgap_sum {k : ℕ} (hk : 0 < k) (f : Fin k → ℝ) :
    ∑ i : Fin k, cgap hk f i = 1 := by
  set F : ℕ → ℝ := fun n => if h : n < k then f ⟨n, h⟩ else 1 + f ⟨0, hk⟩ with hF
  have hstep : ∀ i : Fin k, cgap hk f i = F ((i : ℕ) + 1) - F (i : ℕ) := by
    intro i
    have hFi : F (i : ℕ) = f i := by
      rw [hF]
      simp only [i.isLt, dif_pos]
    rw [hF, cgap, hFi]
  rw [Finset.sum_congr rfl (fun i _ => hstep i),
    Fin.sum_univ_eq_sum_range (fun n => F (n + 1) - F n) k,
    Finset.sum_range_sub F]
  have hFk : F k = 1 + f ⟨0, hk⟩ := by
    rw [hF]; simp only [lt_irrefl, dif_neg, not_false_iff]
  have hF0 : F 0 = f ⟨0, hk⟩ := by
    rw [hF]; simp only [hk, dif_pos]
  rw [hFk, hF0]; ring

/-- **The gap bridge**: ANY cyclic gap `> 1/7` makes its left endpoint a
goodSet-style witness -- every phase difference avoids `(0, 1/7]`.
Interior gap at `s_i` forces `s_i < 6/7`, so wrap terms `1 + s_j - s_i`
clear `1/7` for free; the wrap gap is witnessed at the maximum phase. -/
theorem cgap_witness {k : ℕ} (hk : 0 < k) (f : Fin k → ℝ)
    (hmono : StrictMono f) (hrange : ∀ i, 0 ≤ f i ∧ f i < 1)
    (i : Fin k) (hgap : 1 / 7 < cgap hk f i) :
    ∀ j, Int.fract (f j - f i) ∉ Set.Ioc (0 : ℝ) (1 / 7) := by
  intro j
  by_cases hlt : (i : ℕ) + 1 < k
  · -- interior gap: the successor phase sits above s_i + 1/7
    rw [cgap, dif_pos hlt] at hgap
    have hi67 : f i < 6 / 7 := by
      have h1 := (hrange ⟨(i : ℕ) + 1, hlt⟩).2
      linarith
    rcases lt_trichotomy j i with hji | rfl | hji
    · -- j < i: the wrap term 1 + s_j - s_i > 1/7 since s_i < 6/7
      have hfj : f j < f i := hmono hji
      have hd1 : 0 ≤ f j - f i + 1 := by
        have := (hrange j).1
        have := (hrange i).2
        linarith
      have hd2 : f j - f i + 1 < 1 := by linarith
      have hfr : Int.fract (f j - f i) = f j - f i + 1 := by
        rw [show f j - f i = (f j - f i + 1) + ((-1 : ℤ) : ℝ) by push_cast; ring,
          Int.fract_add_int, Int.fract_eq_self.mpr ⟨hd1, hd2⟩]
      rw [hfr]
      rintro ⟨-, h2⟩
      have := (hrange j).1
      linarith
    · -- j = i: fract 0 = 0 is not in the open-closed arc
      rw [sub_self, Int.fract_zero]
      rintro ⟨h1, -⟩
      exact lt_irrefl 0 h1
    · -- i < j: the forward term is at least the gap
      have hij1 : f ⟨(i : ℕ) + 1, hlt⟩ ≤ f j := by
        apply hmono.monotone
        rw [Fin.le_def]
        exact hji
      have hd0 : 0 ≤ f j - f i := by
        have := hmono.monotone (le_of_lt hji)
        linarith
      have hd1 : f j - f i < 1 := by
        have := (hrange j).2
        have := (hrange i).1
        linarith
      rw [Int.fract_eq_self.mpr ⟨hd0, hd1⟩]
      rintro ⟨-, h2⟩
      linarith
  · -- wrap gap: i is the maximum index, witness the top phase
    rw [cgap, dif_neg hlt] at hgap
    have hik : (i : ℕ) + 1 = k := by
      have := i.isLt
      omega
    rcases lt_trichotomy j i with hji | rfl | hji
    · have hfj : f j < f i := hmono hji
      have hf0 : f ⟨0, hk⟩ ≤ f j := by
        apply hmono.monotone
        rw [Fin.le_def]
        exact Nat.zero_le _
      have hd1 : 0 ≤ f j - f i + 1 := by
        have := (hrange j).1
        have := (hrange i).2
        linarith
      have hd2 : f j - f i + 1 < 1 := by linarith
      have hfr : Int.fract (f j - f i) = f j - f i + 1 := by
        rw [show f j - f i = (f j - f i + 1) + ((-1 : ℤ) : ℝ) by push_cast; ring,
          Int.fract_add_int, Int.fract_eq_self.mpr ⟨hd1, hd2⟩]
      rw [hfr]
      rintro ⟨-, h2⟩
      linarith
    · rw [sub_self, Int.fract_zero]
      rintro ⟨h1, -⟩
      exact lt_irrefl 0 h1
    · -- i < j is impossible: i is the top index
      exfalso
      have hj := j.isLt
      have hij : (i : ℕ) < (j : ℕ) := hji
      omega

/-- **The gap dichotomy** for at most seven phases: either some phase
witnesses (all differences avoid `(0, 1/7]`), or the phases form a PERFECT
NET and the first sorted pair sits at cyclic distance EXACTLY `1/7`. -/
theorem gap_dichotomy {k : ℕ} (hk : 0 < k) (hk7 : k ≤ 7) (f : Fin k → ℝ)
    (hmono : StrictMono f) (hrange : ∀ i, 0 ≤ f i ∧ f i < 1) :
    (∃ i, ∀ j, Int.fract (f j - f i) ∉ Set.Ioc (0 : ℝ) (1 / 7)) ∨
    (∃ i j, Int.fract (f j - f i) = 1 / 7) := by
  by_cases hbig : ∃ i, 1 / 7 < cgap hk f i
  · obtain ⟨i, hi⟩ := hbig
    exact Or.inl ⟨i, cgap_witness hk f hmono hrange i hi⟩
  push_neg at hbig
  -- no big gap: k <= 6 contradicts the pigeonhole, so k = 7
  have hk7' : k = 7 := by
    by_contra hne
    obtain ⟨i, hi⟩ := exists_gap_gt_one_seventh hk (by omega)
      (cgap hk f) (cgap_sum hk f)
    exact absurd hi (not_lt.mpr (hbig i))
  -- all seven gaps <= 1/7 summing to 1: the FIRST gap is exactly 1/7
  have h2k : 1 < k := by omega
  have h0k : 0 < k := hk
  have hg0 : cgap hk f ⟨0, h0k⟩ = 1 / 7 := by
    have hsplit : cgap hk f ⟨0, h0k⟩ +
        ∑ j ∈ Finset.univ.erase ⟨0, h0k⟩, cgap hk f j = 1 := by
      rw [Finset.add_sum_erase _ _ (Finset.mem_univ _)]
      exact cgap_sum hk f
    have herase : ∑ j ∈ Finset.univ.erase ⟨0, h0k⟩, cgap hk f j ≤
        ((k : ℝ) - 1) * (1 / 7) := by
      calc ∑ j ∈ Finset.univ.erase ⟨0, h0k⟩, cgap hk f j
          ≤ ∑ _j ∈ Finset.univ.erase ⟨0, h0k⟩, (1 / 7 : ℝ) :=
            Finset.sum_le_sum (fun j _ => hbig j)
        _ = ((Finset.univ.erase (⟨0, h0k⟩ : Fin k)).card : ℝ) * (1 / 7) := by
            rw [Finset.sum_const, nsmul_eq_mul]
        _ = ((k : ℝ) - 1) * (1 / 7) := by
            rw [Finset.card_erase_of_mem (Finset.mem_univ _), Finset.card_univ,
              Fintype.card_fin]
            have : ((k - 1 : ℕ) : ℝ) = (k : ℝ) - 1 := by
              have : 1 ≤ k := hk
              push_cast [Nat.cast_sub this]
              ring
            rw [this]
    have hcap : ((k : ℝ) - 1) * (1 / 7) ≤ 6 / 7 := by
      have : (k : ℝ) ≤ 7 := by exact_mod_cast hk7
      linarith
    have := hbig ⟨0, h0k⟩
    linarith
  -- extract the exact pair: gap 0 is interior (k = 7 > 1)
  have hlt01 : (0 : ℕ) + 1 < k := by omega
  have h01 : f ⟨1, by omega⟩ - f ⟨0, h0k⟩ = 1 / 7 := by
    rw [cgap, dif_pos hlt01] at hg0
    convert hg0 using 3
  refine Or.inr ⟨⟨0, h0k⟩, ⟨1, by omega⟩, ?_⟩
  rw [h01]
  exact Int.fract_eq_self.mpr ⟨by norm_num, by norm_num⟩

/-! ## The phase level: distinct phases of a cluster -/

/-- **Good-or-net**: at every time `x`, a cluster with at most seven distinct
teeth EITHER lies in `goodSet` OR realizes a phase difference of EXACTLY
`1/7` (the perfect-net event). -/
theorem good_or_pair (E : List ℤ) (hne : E.toFinset.Nonempty)
    (hcard : E.toFinset.card ≤ 7) (x : ℝ) :
    x ∈ goodSet E ∨
    ∃ a ∈ E.toFinset, ∃ b ∈ E.toFinset,
      Int.fract (((b - a : ℤ) : ℝ) * x) = 1 / 7 := by
  classical
  set S : Finset ℝ := E.toFinset.image (fun e : ℤ => Int.fract ((e : ℝ) * x))
    with hS
  have hSne : S.Nonempty := hne.image _
  have hk : 0 < S.card := Finset.card_pos.mpr hSne
  have hk7 : S.card ≤ 7 := le_trans Finset.card_image_le hcard
  set f : Fin S.card → ℝ := fun i => ((S.orderIsoOfFin rfl) i : ℝ) with hf
  have hmono : StrictMono f := fun i j hij =>
    Subtype.coe_lt_coe.mpr ((S.orderIsoOfFin rfl).lt_iff_lt.mpr hij)
  have hmem : ∀ i, f i ∈ S := fun i => ((S.orderIsoOfFin rfl) i).2
  have hrange : ∀ i, 0 ≤ f i ∧ f i < 1 := by
    intro i
    obtain ⟨e, -, hei⟩ := Finset.mem_image.mp (hmem i)
    rw [← hei]
    exact ⟨Int.fract_nonneg _, Int.fract_lt_one _⟩
  have hsurj : ∀ q ∈ S, ∃ j, f j = q := by
    intro q hq
    refine ⟨(S.orderIsoOfFin rfl).symm ⟨q, hq⟩, ?_⟩
    rw [hf]
    simp
  rcases gap_dichotomy hk hk7 f hmono hrange with ⟨i, hwit⟩ | ⟨i, j, hpair⟩
  · -- the witness phase pulls back to a goodSet witness tooth
    obtain ⟨a, haE, hai⟩ := Finset.mem_image.mp (hmem i)
    left
    unfold goodSet
    simp only [Set.mem_iUnion, Set.mem_iInter, Set.mem_setOf_eq]
    refine ⟨a, haE, ?_⟩
    intro b hbE
    have hq : Int.fract ((b : ℝ) * x) ∈ S :=
      Finset.mem_image_of_mem _ hbE
    obtain ⟨j, hj⟩ := hsurj _ hq
    have hid : Int.fract (((b - a : ℤ) : ℝ) * x) = Int.fract (f j - f i) := by
      have harg : (((b - a : ℤ) : ℝ) * x) = (b : ℝ) * x - (a : ℝ) * x := by
        push_cast; ring
      rw [harg, fract_sub_eq_fract_fract_sub_fract, hj, hai]
    rw [hid]
    exact hwit j
  · -- the perfect-net pair pulls back to a tooth pair
    obtain ⟨a, haE, hai⟩ := Finset.mem_image.mp (hmem i)
    obtain ⟨b, hbE, hbj⟩ := Finset.mem_image.mp (hmem j)
    refine Or.inr ⟨a, haE, b, hbE, ?_⟩
    have harg : (((b - a : ℤ) : ℝ) * x) = (b : ℝ) * x - (a : ℝ) * x := by
      push_cast; ring
    rw [harg, fract_sub_eq_fract_fract_sub_fract, hbj, hai]
    exact hpair

/-! ## The null set: perfect-net times are countable -/

/-- Each nonzero tooth difference pins the perfect-net times to a countable
fiber: `fract(d*x) = 1/7` forces `x = (n + 1/7)/d`. -/
theorem pair_fiber_countable {d : ℤ} (hd : d ≠ 0) :
    {x : ℝ | Int.fract ((d : ℝ) * x) = 1 / 7}.Countable := by
  apply Set.Countable.mono ?_ (Set.countable_range
    (fun n : ℤ => ((n : ℝ) + 1 / 7) / (d : ℝ)))
  intro x hx
  simp only [Set.mem_setOf_eq] at hx
  refine ⟨⌊(d : ℝ) * x⌋, ?_⟩
  have hdr : (d : ℝ) ≠ 0 := Int.cast_ne_zero.mpr hd
  have hfl := Int.floor_add_fract ((d : ℝ) * x)
  rw [hx] at hfl
  field_simp
  linarith

/-- **The perfect-net event is null**: the union over tooth pairs of the
exact-`1/7` fibers is countable, hence Lebesgue-null. -/
theorem pairSet_countable (E : List ℤ) :
    (⋃ a ∈ (E.toFinset : Set ℤ), ⋃ b ∈ (E.toFinset : Set ℤ),
      {x : ℝ | Int.fract (((b - a : ℤ) : ℝ) * x) = 1 / 7}).Countable := by
  refine Set.Countable.biUnion (Finset.countable_toSet _) (fun a _ => ?_)
  refine Set.Countable.biUnion (Finset.countable_toSet _) (fun b _ => ?_)
  by_cases hab : b - a = 0
  · have hempty : {x : ℝ | Int.fract (((b - a : ℤ) : ℝ) * x) = 1 / 7} = ∅ := by
      ext x
      simp only [Set.mem_setOf_eq, Set.mem_empty_iff_false, iff_false]
      rw [hab]
      push_cast
      rw [zero_mul, Int.fract_zero]
      norm_num
    rw [hempty]
    exact Set.countable_empty
  · exact pair_fiber_countable hab

/-! ## The de-citation -/

/-- **THE ≤7-ARCS PIGEONHOLE, PROVED (`SmallClusterFull` de-cited)**:
every nodup cluster with `3 ≤ |E| ≤ 7` has full good-set measure.  Off a
countable set of times the distinct phases have a cyclic gap `> 1/7`
(pigeonhole for `k ≤ 6`; for `k = 7` the boundary case is the perfect net,
which pins `x` to countably many fibers). -/
theorem goodSet_ae_full (E : List ℤ) (hnd : E.Nodup) (h3 : 3 ≤ E.length)
    (_h7 : E.length ≤ 7) : (slowμ (goodSet E)).toReal = 1 := by
  have hcardE : E.toFinset.card ≤ 7 := by
    rw [List.toFinset_card_of_nodup hnd]
    omega
  have hne : E.toFinset.Nonempty := by
    apply Finset.card_pos.mp
    rw [List.toFinset_card_of_nodup hnd]
    omega
  set N : Set ℝ := ⋃ a ∈ (E.toFinset : Set ℤ), ⋃ b ∈ (E.toFinset : Set ℤ),
    {x : ℝ | Int.fract (((b - a : ℤ) : ℝ) * x) = 1 / 7} with hN
  have hNnull : volume N = 0 := (pairSet_countable E).measure_zero _
  have hcover : (Set.univ : Set ℝ) ⊆ goodSet E ∪ N := by
    intro x _
    rcases good_or_pair E hne hcardE x with hg | ⟨a, ha, b, hb, hp⟩
    · exact Or.inl hg
    · refine Or.inr ?_
      rw [hN]
      exact Set.mem_iUnion₂.mpr ⟨a, by simpa using ha,
        Set.mem_iUnion₂.mpr ⟨b, by simpa using hb, hp⟩⟩
  have huniv : slowμ (Set.univ : Set ℝ) = 1 := by
    rw [slowμ, Measure.restrict_apply_univ, Real.volume_Ico]
    norm_num
  have hNnull' : slowμ N = 0 := by
    rw [slowμ, Measure.restrict_apply' measurableSet_Ico]
    exact measure_mono_null Set.inter_subset_left hNnull
  have hle1 : slowμ (goodSet E) ≤ 1 := by
    calc slowμ (goodSet E) ≤ slowμ Set.univ := measure_mono (Set.subset_univ _)
      _ = 1 := huniv
  have hge1 : 1 ≤ slowμ (goodSet E) := by
    calc (1 : ℝ≥0∞) = slowμ Set.univ := huniv.symm
      _ ≤ slowμ (goodSet E ∪ N) := measure_mono hcover
      _ ≤ slowμ (goodSet E) + slowμ N := measure_union_le _ _
      _ = slowμ (goodSet E) := by rw [hNnull', add_zero]
  rw [le_antisymm hle1 hge1]
  exact ENNReal.one_toReal

end SevenGapRigidity
end LRC14
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.LRC14.SevenGapRigidity.cgap_sum
#print axioms LonelyRunner.LRC14.SevenGapRigidity.cgap_witness
#print axioms LonelyRunner.LRC14.SevenGapRigidity.gap_dichotomy
#print axioms LonelyRunner.LRC14.SevenGapRigidity.good_or_pair
#print axioms LonelyRunner.LRC14.SevenGapRigidity.pair_fiber_countable
#print axioms LonelyRunner.LRC14.SevenGapRigidity.goodSet_ae_full
