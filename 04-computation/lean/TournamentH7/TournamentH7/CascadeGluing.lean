/- CascadeGluing.lean — klein-2026-07-16-S317, closing opus-S333's three sorries
   (draft: 04-computation/lean-drafts/CascadeGluing.lean).

   THM-928(A) THE CASCADE STEP and THM-932 THE LOCAL-DENSITY BLOCK-GLUING
   LEMMA, the next rungs above the proved FragmentationLemma (klein-S316).

   Rung structure:
     fragmentation       (PROVED, S316) : vol(badArcs w ∩ I) ≤ (L·w + 1)·(2λ/w)
     cascade_step        (PROVED here)  : vol(I \ badArcs w) ≥ L·(1 − 2λ) − 2λ/w
     window_floor_sample (PROVED here)  : the G1 sampling lemma, single interval
     union_floor_sample  (PROVED here)  : the block-gluing bound over k components
   Paper: THM-928 (cascade), THM-932 (gluing; exact-input loss terms κ(V_i)·l).
   Verification: block_gluing_opus_S333.py (G1/G2/G3 exact, 0 violations).
   opus's badArcs' duplicate is dropped — the FragmentationLemma `badArcs` is
   the same set, so the proved bound applies verbatim.  No sorries. -/
import TournamentH7.FragmentationLemma

open MeasureTheory

namespace LRC14

/-- `ofReal` is finitely subadditive over sums (no sign hypotheses). -/
theorem ofReal_sum_le {ι : Type*} (s : Finset ι) (f : ι → ℝ) :
    ENNReal.ofReal (∑ i ∈ s, f i) ≤ ∑ i ∈ s, ENNReal.ofReal (f i) := by
  induction s using Finset.cons_induction with
  | empty => simp
  | cons a s ha ih =>
      rw [Finset.sum_cons, Finset.sum_cons]
      calc ENNReal.ofReal (f a + ∑ i ∈ s, f i)
          ≤ ENNReal.ofReal (f a) + ENNReal.ofReal (∑ i ∈ s, f i) :=
            ENNReal.ofReal_add_le
        _ ≤ ENNReal.ofReal (f a) + ∑ i ∈ s, ENNReal.ofReal (f i) :=
            add_le_add le_rfl ih

/-- **THE CASCADE STEP (THM-928(A) L1, complement form)**: one speed `w` leaves
    at least `L·(1 − 2λ) − 2λ/w` of any interval of length `L` uncovered.
    Direct corollary of the proved `fragmentation` bound. -/
theorem cascade_step (w : ℕ) (hw : 1 ≤ w) (lam L x : ℝ)
    (hlam : 0 < lam) (hL : 0 ≤ L) :
    ENNReal.ofReal (L * (1 - 2 * lam) - 2 * lam / w)
      ≤ volume (Set.Icc x (x + L) \ badArcs w lam) := by
  set t := L * (1 - 2 * lam) - 2 * lam / w with ht
  rcases le_total t 0 with hneg | hpos
  · rw [ENNReal.ofReal_of_nonpos hneg]
    exact zero_le
  · have hw0 : (0:ℝ) < w := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hw
    have hwne : (w:ℝ) ≠ 0 := hw0.ne'
    have hfrag := fragmentation w hw lam L x hlam hL
    have hIcc : volume (Set.Icc x (x + L)) = ENNReal.ofReal L := by
      rw [Real.volume_Icc]
      congr 1
      ring
    have hfin : volume (badArcs w lam ∩ Set.Icc x (x + L)) ≠ ⊤ :=
      ne_top_of_le_ne_top ENNReal.ofReal_ne_top hfrag
    have hsplit : volume (Set.Icc x (x + L))
        ≤ volume (Set.Icc x (x + L) \ badArcs w lam)
          + volume (badArcs w lam ∩ Set.Icc x (x + L)) := by
      refine le_trans (measure_mono ?_) (measure_union_le _ _)
      intro y hy
      by_cases hb : y ∈ badArcs w lam
      · exact Or.inr ⟨hb, hy⟩
      · exact Or.inl ⟨hy, hb⟩
    have hnn : (0:ℝ) ≤ (L * w + 1) * (2 * lam / w) := by
      have h1 : (0:ℝ) ≤ L * w := mul_nonneg hL hw0.le
      have h2 : (0:ℝ) ≤ 2 * lam / w := div_nonneg (by linarith) hw0.le
      nlinarith
    have hbound_sum : t + (L * w + 1) * (2 * lam / w) = L := by
      rw [ht]
      field_simp
      ring
    have h2 : ENNReal.ofReal t + ENNReal.ofReal ((L * w + 1) * (2 * lam / w))
        = ENNReal.ofReal L := by
      rw [← ENNReal.ofReal_add hpos hnn, hbound_sum]
    have h3 : ENNReal.ofReal t + volume (badArcs w lam ∩ Set.Icc x (x + L))
        ≤ volume (Set.Icc x (x + L) \ badArcs w lam)
          + volume (badArcs w lam ∩ Set.Icc x (x + L)) :=
      calc ENNReal.ofReal t + volume (badArcs w lam ∩ Set.Icc x (x + L))
          ≤ ENNReal.ofReal t + ENNReal.ofReal ((L * w + 1) * (2 * lam / w)) :=
            add_le_add le_rfl hfrag
        _ = ENNReal.ofReal L := h2
        _ = volume (Set.Icc x (x + L)) := hIcc.symm
        _ ≤ _ := hsplit
    exact (ENNReal.add_le_add_iff_right hfin).mp h3

/-- **THE SAMPLING LEMMA (THM-932 G1, single-interval form)**: if every window
    of length `l` meets `W` in measure at least `η·l`, then an interval of
    length `L ≥ l` meets `W` in measure at least `η·(L − l)`.  Proof: tile from
    the left by `⌊L/l⌋` disjoint `Ico`-windows; superadditivity; `⌊L/l⌋·l ≥ L − l`. -/
theorem window_floor_sample (W : Set ℝ) (hW : MeasurableSet W)
    (l eta L x : ℝ) (hl : 0 < l) (heta : 0 ≤ eta) (hL : l ≤ L)
    (hfloor : ∀ a : ℝ, ENNReal.ofReal (eta * l)
                ≤ volume (W ∩ Set.Icc a (a + l))) :
    ENNReal.ofReal (eta * (L - l)) ≤ volume (W ∩ Set.Icc x (x + L)) := by
  set N : ℕ := Nat.floor (L / l) with hN
  have hL0 : (0:ℝ) ≤ L := le_trans hl.le hL
  have hLl : (0:ℝ) ≤ L / l := div_nonneg hL0 hl.le
  -- the Icc floor transfers to the half-open window (endpoint is null)
  have htile : ∀ a : ℝ, ENNReal.ofReal (eta * l) ≤ volume (W ∩ Set.Ico a (a + l)) := by
    intro a
    have hsplit : W ∩ Set.Icc a (a + l) ⊆ (W ∩ Set.Ico a (a + l)) ∪ {a + l} := by
      rintro y ⟨hyW, hy1, hy2⟩
      rcases lt_or_eq_of_le hy2 with h | h
      · exact Or.inl ⟨hyW, hy1, h⟩
      · exact Or.inr h
    calc ENNReal.ofReal (eta * l)
        ≤ volume (W ∩ Set.Icc a (a + l)) := hfloor a
      _ ≤ volume ((W ∩ Set.Ico a (a + l)) ∪ {a + l}) := measure_mono hsplit
      _ ≤ volume (W ∩ Set.Ico a (a + l)) + volume ({a + l} : Set ℝ) :=
          measure_union_le _ _
      _ = volume (W ∩ Set.Ico a (a + l)) := by
          rw [Real.volume_singleton, add_zero]
  -- the tiles are pairwise disjoint
  have hdisj : (↑(Finset.range N) : Set ℕ).PairwiseDisjoint
      (fun i : ℕ => W ∩ Set.Ico (x + i * l) (x + i * l + l)) := by
    intro i _ j _ hij
    have hIco : Disjoint (Set.Ico (x + (i:ℝ) * l) (x + i * l + l))
        (Set.Ico (x + (j:ℝ) * l) (x + j * l + l)) := by
      rcases lt_or_gt_of_ne hij with h | h
      · have hij1 : ((i:ℝ) + 1) ≤ (j:ℝ) := by exact_mod_cast h
        have := mul_le_mul_of_nonneg_right hij1 hl.le
        refine Set.Ico_disjoint_Ico.mpr ?_
        refine le_trans (min_le_left _ _) (le_trans ?_ (le_max_right _ _))
        linarith
      · have hij1 : ((j:ℝ) + 1) ≤ (i:ℝ) := by exact_mod_cast h
        have := mul_le_mul_of_nonneg_right hij1 hl.le
        refine Set.Ico_disjoint_Ico.mpr ?_
        refine le_trans (min_le_right _ _) (le_trans ?_ (le_max_left _ _))
        linarith
    exact hIco.mono Set.inter_subset_right Set.inter_subset_right
  have hmeas : ∀ i : ℕ, MeasurableSet (W ∩ Set.Ico (x + i * l) (x + i * l + l)) :=
    fun i => hW.inter measurableSet_Ico
  -- the tiled union sits inside the big interval
  have hNle : (N:ℝ) * l ≤ L := by
    have hfl : (N:ℝ) ≤ L / l := by
      have := Nat.floor_le hLl
      exact_mod_cast this
    have := mul_le_mul_of_nonneg_right hfl hl.le
    calc (N:ℝ) * l ≤ (L / l) * l := this
      _ = L := by field_simp
  have hsub : (⋃ i ∈ Finset.range N, Set.Ico (x + (i:ℝ) * l) (x + i * l + l))
      ⊆ Set.Icc x (x + L) := by
    intro y hy
    simp only [Set.mem_iUnion, Finset.mem_range, Set.mem_Ico] at hy
    obtain ⟨i, hiN, hy1, hy2⟩ := hy
    have hi0 : (0:ℝ) ≤ (i:ℝ) * l := by positivity
    have hi1 : ((i:ℝ) + 1) ≤ (N:ℝ) := by exact_mod_cast hiN
    have := mul_le_mul_of_nonneg_right hi1 hl.le
    constructor
    · linarith
    · linarith
  -- superadditivity over the tiles
  have hcount : ∑ i ∈ Finset.range N, volume (W ∩ Set.Ico (x + (i:ℝ) * l) (x + i * l + l))
      ≤ volume (W ∩ Set.Icc x (x + L)) := by
    rw [← measure_biUnion_finset hdisj (fun i _ => hmeas i)]
    refine measure_mono ?_
    intro y hy
    simp only [Set.mem_iUnion] at hy
    obtain ⟨i, hi, hyW, hyt⟩ := hy
    exact ⟨hyW, hsub (Set.mem_biUnion hi hyt)⟩
  -- the count: ⌊L/l⌋ tiles carry η·l each, and ⌊L/l⌋·l ≥ L − l
  have hNge : L - l ≤ (N:ℝ) * l := by
    have h2 : L / l < (N:ℝ) + 1 := by exact_mod_cast Nat.lt_floor_add_one (L / l)
    have := mul_lt_mul_of_pos_right h2 hl
    rw [div_mul_cancel₀ _ hl.ne'] at this
    nlinarith
  calc ENNReal.ofReal (eta * (L - l))
      ≤ ENNReal.ofReal (eta * ((N:ℝ) * l)) :=
        ENNReal.ofReal_le_ofReal (mul_le_mul_of_nonneg_left hNge heta)
    _ = (N : ENNReal) * ENNReal.ofReal (eta * l) := by
        rw [show eta * ((N:ℝ) * l) = (N:ℝ) * (eta * l) by ring,
            ENNReal.ofReal_mul (by positivity : (0:ℝ) ≤ (N:ℝ)),
            ENNReal.ofReal_natCast]
    _ = ∑ _i ∈ Finset.range N, ENNReal.ofReal (eta * l) := by
        rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
    _ ≤ ∑ i ∈ Finset.range N, volume (W ∩ Set.Ico (x + (i:ℝ) * l) (x + i * l + l)) :=
        Finset.sum_le_sum fun i _ => htile _
    _ ≤ volume (W ∩ Set.Icc x (x + L)) := hcount

/-- **THE BLOCK-GLUING BOUND (THM-932, k-component form)**: if `V` is a
    disjoint union of `k` closed intervals and `W` has local density floor
    `η` at scale `l`, then `vol(V ∩ W) ≥ η·(vol V − k·l)`.  Components shorter
    than `l` (or degenerate) are absorbed by `ofReal` subadditivity. -/
theorem union_floor_sample (W : Set ℝ) (hW : MeasurableSet W)
    (k : ℕ) (endpoints : Fin k → ℝ × ℝ)
    (hdisj : ∀ i j, i ≠ j → Set.Icc (endpoints i).1 (endpoints i).2 ∩
                            Set.Icc (endpoints j).1 (endpoints j).2 = ∅)
    (l eta : ℝ) (hl : 0 < l) (heta : 0 ≤ eta)
    (hfloor : ∀ a : ℝ, ENNReal.ofReal (eta * l)
                ≤ volume (W ∩ Set.Icc a (a + l))) :
    ENNReal.ofReal (eta * ((∑ i, ((endpoints i).2 - (endpoints i).1)) - k * l))
      ≤ volume (W ∩ ⋃ i, Set.Icc (endpoints i).1 (endpoints i).2) := by
  have hpair : Pairwise (Function.onFun Disjoint
      (fun i => W ∩ Set.Icc (endpoints i).1 (endpoints i).2)) := by
    intro i j hij
    have hd : Disjoint (Set.Icc (endpoints i).1 (endpoints i).2)
        (Set.Icc (endpoints j).1 (endpoints j).2) :=
      Set.disjoint_iff_inter_eq_empty.mpr (hdisj i j hij)
    exact hd.mono Set.inter_subset_right Set.inter_subset_right
  have hmeas : ∀ i : Fin k,
      MeasurableSet (W ∩ Set.Icc (endpoints i).1 (endpoints i).2) :=
    fun i => hW.inter measurableSet_Icc
  rw [Set.inter_iUnion, measure_iUnion hpair hmeas, tsum_fintype]
  -- per-component floor, short components absorbed
  have hstep : ∀ i : Fin k,
      ENNReal.ofReal (eta * (((endpoints i).2 - (endpoints i).1) - l))
        ≤ volume (W ∩ Set.Icc (endpoints i).1 (endpoints i).2) := by
    intro i
    rcases le_total l ((endpoints i).2 - (endpoints i).1) with hcase | hcase
    · have h := window_floor_sample W hW l eta
        ((endpoints i).2 - (endpoints i).1) (endpoints i).1 hl heta hcase hfloor
      rwa [show (endpoints i).1 + ((endpoints i).2 - (endpoints i).1)
            = (endpoints i).2 by ring] at h
    · have hnp : eta * (((endpoints i).2 - (endpoints i).1) - l) ≤ 0 := by
        have h1 : ((endpoints i).2 - (endpoints i).1) - l ≤ 0 := by linarith
        calc eta * (((endpoints i).2 - (endpoints i).1) - l)
            ≤ eta * 0 := mul_le_mul_of_nonneg_left h1 heta
          _ = 0 := mul_zero eta
      rw [ENNReal.ofReal_of_nonpos hnp]
      exact zero_le
  -- the real identity Σ η·(len_i − l) = η·(Σ len − k·l), then sum up
  have hid : eta * ((∑ i, ((endpoints i).2 - (endpoints i).1)) - k * l)
      = ∑ i, eta * (((endpoints i).2 - (endpoints i).1) - l) := by
    have h1 : ∀ i : Fin k, eta * (((endpoints i).2 - (endpoints i).1) - l)
        = eta * ((endpoints i).2 - (endpoints i).1) - eta * l := fun i => by ring
    simp only [h1]
    conv_rhs => rw [Finset.sum_sub_distrib]
    rw [← Finset.mul_sum, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
        nsmul_eq_mul]
    ring
  calc ENNReal.ofReal (eta * ((∑ i, ((endpoints i).2 - (endpoints i).1)) - k * l))
      = ENNReal.ofReal (∑ i, eta * (((endpoints i).2 - (endpoints i).1) - l)) := by
        rw [hid]
    _ ≤ ∑ i, ENNReal.ofReal (eta * (((endpoints i).2 - (endpoints i).1) - l)) :=
        ofReal_sum_le _ _
    _ ≤ ∑ i, volume (W ∩ Set.Icc (endpoints i).1 (endpoints i).2) :=
        Finset.sum_le_sum fun i _ => hstep i

end LRC14
