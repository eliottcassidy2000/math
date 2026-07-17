/- CascadeGluing.lean -- opus-2026-07-16-S334.
   The MEASURE-THEORETIC layer of THM-928(A)/THM-932/THM-933: connects real
   `volume` on ℝ to the algebraic gluing ledgers of
   LRCLocalDensityBlockGluing (codex S21, sorry-free).

   Three rungs above the proved `fragmentation` (FragmentationLemma.lean):

   * `cascade_step`      -- one speed leaves ≥ L(1−2λ) − 2λ/w of any length-L
                            interval uncovered (THM-928(A) L1, complement form);
                            a direct corollary of `fragmentation`.
   * `window_floor_sample` -- the G1 sampling lemma (THM-932): a local density
                            floor at scale l on all windows lifts to any
                            interval, price one window: η(L − l).
   * `union_floor_sample`  -- G1 summed over a disjoint interval family,
                            price one window per component: η(Σlen − k·l).

   Paper: THM-928/932; exact verification: block_gluing_opus_S333.py
   (0 violations, 40+40+12 configs).  Supersedes the statements drafted in
   04-computation/lean-drafts/CascadeGluing.lean. -/
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic
import TournamentH7.FragmentationLemma

open MeasureTheory

namespace LRC14

/-- `badArcs` is measurable (countable union of open intervals). -/
theorem badArcs_measurableSet (w : ℕ) (lam : ℝ) :
    MeasurableSet (badArcs w lam) :=
  MeasurableSet.iUnion fun _ => measurableSet_Ioo

/-- THE CASCADE STEP (THM-928(A) L1, complement form): one speed of modulus
`w` leaves at least `L·(1 − 2·lam) − 2·lam/w` of any interval of length `L`
uncovered.  Direct corollary of `fragmentation`. -/
theorem cascade_step (w : ℕ) (hw : 1 ≤ w) (lam L x : ℝ)
    (hlam : 0 < lam) (hL : 0 ≤ L) :
    ENNReal.ofReal (L * (1 - 2 * lam) - 2 * lam / w)
      ≤ volume (Set.Icc x (x + L) \ badArcs w lam) := by
  have hwpos : (0 : ℝ) < w := by exact_mod_cast hw
  have hIvol : volume (Set.Icc x (x + L)) = ENNReal.ofReal L := by
    rw [Real.volume_Icc]
    congr 1
    ring
  have hdiffeq : Set.Icc x (x + L) \ badArcs w lam
      = Set.Icc x (x + L) \ (badArcs w lam ∩ Set.Icc x (x + L)) := by
    ext t
    simp only [Set.mem_diff, Set.mem_inter_iff]
    tauto
  have hfin : volume (badArcs w lam ∩ Set.Icc x (x + L)) ≠ ⊤ := by
    refine ne_top_of_le_ne_top ?_ (measure_mono Set.inter_subset_right)
    rw [hIvol]
    exact ENNReal.ofReal_ne_top
  have hkey : volume (Set.Icc x (x + L) \ badArcs w lam)
      = volume (Set.Icc x (x + L))
        - volume (badArcs w lam ∩ Set.Icc x (x + L)) := by
    rw [hdiffeq]
    exact measure_diff Set.inter_subset_right
      ((badArcs_measurableSet w lam).inter measurableSet_Icc).nullMeasurableSet
      hfin
  have hfrag := fragmentation w hw lam L x hlam hL
  have hbnd_nonneg : (0 : ℝ) ≤ (L * w + 1) * (2 * lam / w) := by positivity
  calc ENNReal.ofReal (L * (1 - 2 * lam) - 2 * lam / w)
      ≤ ENNReal.ofReal (L - (L * w + 1) * (2 * lam / w)) := by
        apply ENNReal.ofReal_le_ofReal
        have hw0 : (w : ℝ) ≠ 0 := ne_of_gt hwpos
        have hexp : (L * w + 1) * (2 * lam / w) = 2 * lam * L + 2 * lam / w := by
          field_simp
          try ring
        rw [hexp]
        nlinarith [hlam.le, hL]
    _ = ENNReal.ofReal L - ENNReal.ofReal ((L * w + 1) * (2 * lam / w)) :=
        ENNReal.ofReal_sub _ hbnd_nonneg
    _ ≤ ENNReal.ofReal L - volume (badArcs w lam ∩ Set.Icc x (x + L)) :=
        tsub_le_tsub_left hfrag _
    _ = volume (Set.Icc x (x + L) \ badArcs w lam) := by rw [hkey, hIvol]

/-- `ofReal` is subadditive over finite sums. -/
theorem ofReal_sum_le {ι : Type*} (s : Finset ι) (f : ι → ℝ) :
    ENNReal.ofReal (∑ i ∈ s, f i) ≤ ∑ i ∈ s, ENNReal.ofReal (f i) := by
  classical
  refine Finset.induction_on s (by simp) ?_
  intro a t ha ih
  rw [Finset.sum_insert ha, Finset.sum_insert ha]
  exact le_trans ENNReal.ofReal_add_le (add_le_add le_rfl ih)

/-- Windows may be taken half-open at the same price: the floor hypothesis on
closed windows transfers to `Ico` windows (the endpoint is null). -/
theorem floor_on_Ico (W : Set ℝ) (l eta : ℝ)
    (hfloor : ∀ a : ℝ, ENNReal.ofReal (eta * l)
                ≤ volume (W ∩ Set.Icc a (a + l))) (a : ℝ) :
    ENNReal.ofReal (eta * l) ≤ volume (W ∩ Set.Ico a (a + l)) := by
  refine le_trans (hfloor a) ?_
  have hsplit : W ∩ Set.Icc a (a + l)
      ⊆ (W ∩ Set.Ico a (a + l)) ∪ {a + l} := by
    rintro t ⟨htW, hta, htb⟩
    rcases lt_or_eq_of_le htb with h | h
    · exact Or.inl ⟨htW, hta, h⟩
    · exact Or.inr (by simp [h])
  calc volume (W ∩ Set.Icc a (a + l))
      ≤ volume ((W ∩ Set.Ico a (a + l)) ∪ {a + l}) := measure_mono hsplit
    _ ≤ volume (W ∩ Set.Ico a (a + l)) + volume ({a + l} : Set ℝ) :=
        measure_union_le _ _
    _ = volume (W ∩ Set.Ico a (a + l)) := by
        rw [Real.volume_singleton, add_zero]

/-- THE SAMPLING LEMMA (THM-932 G1, single interval): a local density floor
`eta` at scale `l` (all windows) gives every interval of length `L ≥ l` at
least `eta · (L − l)` of `W`.  Proof: tile from the left by `⌊L/l⌋` disjoint
half-open windows. -/
theorem window_floor_sample (W : Set ℝ) (hW : MeasurableSet W)
    (l eta L x : ℝ) (hl : 0 < l) (heta : 0 ≤ eta) (hL : l ≤ L)
    (hfloor : ∀ a : ℝ, ENNReal.ofReal (eta * l)
                ≤ volume (W ∩ Set.Icc a (a + l))) :
    ENNReal.ofReal (eta * (L - l)) ≤ volume (W ∩ Set.Icc x (x + L)) := by
  set N : ℕ := ⌊L / l⌋₊ with hNdef
  have hLpos : 0 < L := lt_of_lt_of_le hl hL
  -- tile i occupies [x + i·l, x + i·l + l)
  set T : Fin N → Set ℝ := fun i => Set.Ico (x + i * l) (x + i * l + l)
    with hTdef
  -- (1) N·l ≤ L
  have hNl : (N : ℝ) * l ≤ L := by
    have hfl : (N : ℝ) ≤ L / l := Nat.floor_le (by positivity)
    calc (N : ℝ) * l ≤ (L / l) * l := by
          exact mul_le_mul_of_nonneg_right hfl (le_of_lt hl)
      _ = L := by field_simp
  -- (2) L − l ≤ N·l
  have hcount : L - l ≤ (N : ℝ) * l := by
    have h1 : L / l - 1 < (N : ℝ) := Nat.sub_one_lt_floor (L / l)
    have h2 : (L / l - 1) * l < (N : ℝ) * l :=
      mul_lt_mul_of_pos_right h1 hl
    have h3 : (L / l - 1) * l = L - l := by field_simp
    linarith
  -- (3) each tile sits inside the big interval
  have hsub : ∀ i : Fin N, W ∩ T i ⊆ W ∩ Set.Icc x (x + L) := by
    intro i t ⟨htW, ht1, ht2⟩
    refine ⟨htW, ?_, ?_⟩
    · have : (0 : ℝ) ≤ (i : ℝ) * l := by positivity
      linarith
    · have hi1 : ((i : ℕ) : ℝ) + 1 ≤ (N : ℝ) := by
        exact_mod_cast (i.isLt : (i : ℕ) + 1 ≤ N)
      have : (i : ℝ) * l + l ≤ (N : ℝ) * l := by nlinarith
      linarith
  -- (4) tiles are pairwise disjoint
  have hdisj : Pairwise (Function.onFun Disjoint fun i : Fin N => W ∩ T i) := by
    intro i j hij
    have key : ∀ a b : Fin N, (a : ℕ) < (b : ℕ) →
        Disjoint (W ∩ T a) (W ∩ T b) := by
      intro a b hab
      apply Set.disjoint_left.mpr
      rintro t ⟨-, -, hta2⟩ ⟨-, htb1, -⟩
      have hcast : ((a : ℕ) : ℝ) + 1 ≤ ((b : ℕ) : ℝ) := by exact_mod_cast hab
      have : x + (a : ℝ) * l + l ≤ x + (b : ℝ) * l := by nlinarith
      linarith
    have hne : (i : ℕ) ≠ (j : ℕ) := fun h => hij (Fin.val_injective h)
    rcases lt_or_gt_of_ne hne with h | h
    · exact key i j h
    · exact (key j i h).symm
  -- (5) assemble
  have hmeas : ∀ i : Fin N, MeasurableSet (W ∩ T i) := fun i =>
    hW.inter measurableSet_Ico
  calc ENNReal.ofReal (eta * (L - l))
      ≤ ENNReal.ofReal (eta * ((N : ℝ) * l)) := by
        apply ENNReal.ofReal_le_ofReal
        exact mul_le_mul_of_nonneg_left hcount heta
    _ = ENNReal.ofReal (∑ _i : Fin N, eta * l) := by
        congr 1
        rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin,
          nsmul_eq_mul]
        ring
    _ ≤ ∑ _i : Fin N, ENNReal.ofReal (eta * l) := ofReal_sum_le _ _
    _ ≤ ∑ i : Fin N, volume (W ∩ T i) :=
        Finset.sum_le_sum fun i _ => floor_on_Ico W l eta hfloor _
    _ = volume (⋃ i ∈ Finset.univ, (W ∩ T i)) :=
        (measure_biUnion_finset (hdisj.set_pairwise _)
          (fun i _ => hmeas i)).symm
    _ ≤ volume (W ∩ Set.Icc x (x + L)) := by
        apply measure_mono
        exact Set.iUnion₂_subset fun i _ => hsub i

/-- G1 over a disjoint family of intervals: price one window per component.
Components shorter than `l` (or degenerate) are handled by the sign of their
contribution — no ordering or length hypothesis is needed per component. -/
theorem union_floor_sample (W : Set ℝ) (hW : MeasurableSet W)
    (k : ℕ) (lo hi : Fin k → ℝ)
    (hdisj : Pairwise (Function.onFun Disjoint
      fun i => Set.Icc (lo i) (hi i)))
    (l eta : ℝ) (hl : 0 < l) (heta : 0 ≤ eta)
    (hfloor : ∀ a : ℝ, ENNReal.ofReal (eta * l)
                ≤ volume (W ∩ Set.Icc a (a + l))) :
    ENNReal.ofReal (eta * ((∑ i, (hi i - lo i)) - k * l))
      ≤ volume (W ∩ ⋃ i, Set.Icc (lo i) (hi i)) := by
  -- per-component bound, all cases
  have hcomp : ∀ i : Fin k, ENNReal.ofReal (eta * (hi i - lo i - l))
      ≤ volume (W ∩ Set.Icc (lo i) (hi i)) := by
    intro i
    by_cases hshort : hi i - lo i ≤ l
    · -- short component: LHS = 0
      have : eta * (hi i - lo i - l) ≤ 0 :=
        mul_nonpos_of_nonneg_of_nonpos heta (by linarith)
      rw [ENNReal.ofReal_eq_zero.mpr this]
      exact zero_le'
    · push_neg at hshort
      have := window_floor_sample W hW l eta (hi i - lo i) (lo i) hl heta
        (le_of_lt hshort) hfloor
      have hrw : lo i + (hi i - lo i) = hi i := by ring
      rwa [hrw] at this
  -- distribute W ∩ over the union, sum disjointly
  have hmeas : ∀ i : Fin k, MeasurableSet (W ∩ Set.Icc (lo i) (hi i)) :=
    fun i => hW.inter measurableSet_Icc
  have hdisj' : Pairwise (Function.onFun Disjoint
      fun i => W ∩ Set.Icc (lo i) (hi i)) := by
    intro i j hij
    exact Disjoint.mono Set.inter_subset_right Set.inter_subset_right
      (hdisj hij)
  have hunion : W ∩ ⋃ i, Set.Icc (lo i) (hi i)
      = ⋃ i, (W ∩ Set.Icc (lo i) (hi i)) := by
    rw [Set.inter_iUnion]
  have halg : (∑ i : Fin k, eta * (hi i - lo i - l))
      = eta * ((∑ i, (hi i - lo i)) - k * l) := by
    have h1 : (∑ i : Fin k, eta * (hi i - lo i - l))
        = ∑ i : Fin k, (eta * (hi i - lo i) - eta * l) :=
      Finset.sum_congr rfl fun i _ => by ring
    rw [h1, Finset.sum_sub_distrib, ← Finset.mul_sum, Finset.sum_const,
      Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    ring
  calc ENNReal.ofReal (eta * ((∑ i, (hi i - lo i)) - k * l))
      = ENNReal.ofReal (∑ i : Fin k, eta * (hi i - lo i - l)) := by
        rw [halg]
    _ ≤ ∑ i : Fin k, ENNReal.ofReal (eta * (hi i - lo i - l)) :=
        ofReal_sum_le _ _
    _ ≤ ∑ i : Fin k, volume (W ∩ Set.Icc (lo i) (hi i)) :=
        Finset.sum_le_sum fun i _ => hcomp i
    _ = volume (⋃ i ∈ Finset.univ, (W ∩ Set.Icc (lo i) (hi i))) :=
        (measure_biUnion_finset (hdisj'.set_pairwise _)
          (fun i _ => hmeas i)).symm
    _ = volume (W ∩ ⋃ i, Set.Icc (lo i) (hi i)) := by
        congr 1
        rw [hunion]
        ext t
        simp

/-! ## Axiom audit -/

#print axioms cascade_step
#print axioms floor_on_Ico
#print axioms window_floor_sample
#print axioms union_floor_sample

end LRC14
