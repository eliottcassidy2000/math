/-
  TournamentH7.LRCGoodDilation — DILATION INVARIANCE of the covering good-set
  (opus-2026-07-08-S162).

  The covering-side companion to `LRCDilationInvariance.iSup_margin_const_mul` (which proves
  dilation invariance of the loneliness minimax `M(c·v)=M(v)`).  Here the same gauge fact is
  proved for the GOOD SET `Good θ E` / its measure `muGood θ E` of `LRCTailDiameter`, and for
  the orbit-`Int.fract` form directly:
      Good θ (c·E) = (x ↦ c·x)⁻¹ (Good θ E),   and   Good θ E  is 1-periodic.

  Role.  This is the structural fact underpinning the k=11 covering-tail reduction to the
  *longest-AP* axis (opus-S155–S161): `muGood θ` (equivalently the degree-3 floor `D3`) depends
  only on the dilation class of the speed set, so every tail shape may be taken PRIMITIVE and
  the extremal analysis lives on the dilation-invariant longest AP.  The fixed-window "cluster
  size" of the refuted LEM-009 argument is NOT dilation-invariant; `muGood`/`D3` IS — the
  exact counterexample `(0,3,6,8,9,12,15,18,21,24,27)` (a *dilated* AP₁₀ + interior point,
  `D3 = 0.452986 < 0.4646`) is a dilate of the compact minimizer, invisible to window-cluster
  but not to `muGood` (MISTAKE-126 / CASE-tail-D3-min-is-not-block-outlier).

  Kernel-pure; mirrors the translation-invariance proofs in `LRCTailDiameter`.
-/
import TournamentH7.LRCTailDiameter

namespace LonelyRunner
namespace TailDiameter

open scoped ENNReal

/-! ### Dilation of the orbit (exact, elementary) -/

/-- **Orbit dilation.**  The `(c·e)`-phase at witness `a` and slow phase `x` equals the
`e`-phase at witness `a` and sped-up phase `c·x` — the same real number `(e)(c·x) − a`
inside `Int.fract`.  So an arc is empty of the `c·E`-orbit at `x` iff it is empty of the
`E`-orbit at `c·x`. -/
theorem emptyArc_dilate {θ : ℝ} (E : Finset ℤ) (c : ℤ) (x a : ℝ) :
    EmptyArc θ (E.image (fun e => c * e)) x a ↔ EmptyArc θ E ((c : ℝ) * x) a := by
  constructor
  · intro h e he
    have hmem : c * e ∈ E.image (fun e => c * e) := Finset.mem_image_of_mem _ he
    have := h (c * e) hmem
    have harg : ((c * e : ℤ) : ℝ) * x - a = (e : ℝ) * ((c : ℝ) * x) - a := by
      push_cast; ring
    rwa [harg] at this
  · intro h f hf
    rcases Finset.mem_image.mp hf with ⟨e, he, rfl⟩
    have := h e he
    have harg : ((c * e : ℤ) : ℝ) * x - a = (e : ℝ) * ((c : ℝ) * x) - a := by
      push_cast; ring
    rwa [harg]

/-- **`Good θ (c·E) = (·c)⁻¹ (Good θ E)`** : existence of an empty arc for the dilated set
at `x` is existence of one for the original set at `c·x`. -/
theorem good_dilate (θ : ℝ) (E : Finset ℤ) (c : ℤ) :
    Good θ (E.image (fun e => c * e)) = (fun x => (c : ℝ) * x) ⁻¹' (Good θ E) := by
  ext x
  constructor
  · rintro ⟨a, ha⟩
    exact ⟨a, (emptyArc_dilate E c x a).mp ha⟩
  · rintro ⟨a, ha⟩
    exact ⟨a, (emptyArc_dilate E c x a).mpr ha⟩

/-! ### 1-periodicity of the good set -/

/-- The `e`-phase is unchanged by shifting the slow phase by `1`: it adds the integer `e`
inside `Int.fract`, which `Int.fract` erases. -/
theorem emptyArc_add_one {θ : ℝ} (E : Finset ℤ) (x a : ℝ) :
    EmptyArc θ E (x + 1) a ↔ EmptyArc θ E x a := by
  constructor
  · intro h e he
    have := h e he
    have harg : (e : ℝ) * (x + 1) - a = ((e : ℝ) * x - a) + (e : ℤ) := by ring
    rw [harg, Int.fract_add_intCast] at this
    exact this
  · intro h e he
    have := h e he
    have harg : (e : ℝ) * (x + 1) - a = ((e : ℝ) * x - a) + (e : ℤ) := by ring
    rw [harg, Int.fract_add_intCast]
    exact this

/-- **`Good θ E` is 1-periodic.**  Together with `good_dilate` this is the structural core of
`muGood`-dilation invariance: `x ↦ c·x` is a measure-balanced `c`-fold cover of the period. -/
theorem good_add_one (θ : ℝ) (E : Finset ℤ) (x : ℝ) :
    x + 1 ∈ Good θ E ↔ x ∈ Good θ E := by
  constructor
  · rintro ⟨a, ha⟩; exact ⟨a, (emptyArc_add_one E x a).mp ha⟩
  · rintro ⟨a, ha⟩; exact ⟨a, (emptyArc_add_one E x a).mpr ha⟩

/-- `Good θ E` is invariant under natural-number translation (iterated `good_add_one`). -/
theorem good_add_natCast (θ : ℝ) (E : Finset ℤ) (n : ℕ) (x : ℝ) :
    x + (n : ℝ) ∈ Good θ E ↔ x ∈ Good θ E := by
  induction n with
  | zero => simp
  | succ k ih =>
    have hstep : x + ((k + 1 : ℕ) : ℝ) = (x + (k : ℝ)) + 1 := by push_cast; ring
    rw [hstep, good_add_one, ih]

/-! ### The measure fold and dilation invariance -/

open MeasureTheory Set

/-- **Periodicity fold.**  For the 1-periodic good set `S = Good θ E`,
`vol(S ∩ [0,n]) = n · vol(S ∩ [0,1])`.  Proved via `Measure.restrict` additivity (so no
measurability of `S` is needed) + Lebesgue translation invariance + the integer-translation
invariance of `Good`. -/
theorem muGood_fold (θ : ℝ) (E : Finset ℤ) (n : ℕ) :
    volume (Good θ E ∩ Icc 0 (n : ℝ)) = n • volume (Good θ E ∩ Icc 0 1) := by
  set S := Good θ E with hS
  induction n with
  | zero =>
    simp only [Nat.cast_zero, zero_smul, Set.Icc_self]
    exact measure_mono_null (Set.inter_subset_right) (by simp)
  | succ k ih =>
    have hcast : ((k + 1 : ℕ) : ℝ) = (k : ℝ) + 1 := by push_cast; ring
    -- disjoint measurable split  Icc 0 (k+1) = Icc 0 k ∪ Ioc k (k+1)
    have hdisj : Disjoint (Icc (0 : ℝ) (k : ℝ)) (Ioc (k : ℝ) ((k : ℝ) + 1)) := by
      rw [Set.disjoint_left]; rintro y ⟨_, hy2⟩ ⟨hy3, _⟩; exact absurd hy3 (not_lt.mpr hy2)
    have huni : Icc (0 : ℝ) ((k : ℝ) + 1) = Icc 0 (k : ℝ) ∪ Ioc (k : ℝ) ((k : ℝ) + 1) :=
      (Set.Icc_union_Ioc_eq_Icc (by positivity) (by linarith)).symm
    -- restrict-additivity: vol(S ∩ (A∪B)) = vol(S∩A) + vol(S∩B), A,B measurable
    have hrec : volume (S ∩ Icc 0 ((k : ℝ) + 1))
        = volume (S ∩ Icc 0 (k : ℝ)) + volume (S ∩ Ioc (k : ℝ) ((k : ℝ) + 1)) := by
      calc volume (S ∩ Icc 0 ((k : ℝ) + 1))
          = volume (Icc 0 ((k : ℝ) + 1) ∩ S) := by rw [Set.inter_comm]
        _ = (volume.restrict S) (Icc 0 ((k : ℝ) + 1)) :=
              (Measure.restrict_apply measurableSet_Icc).symm
        _ = (volume.restrict S) (Icc 0 (k : ℝ) ∪ Ioc (k : ℝ) ((k : ℝ) + 1)) := by rw [huni]
        _ = (volume.restrict S) (Icc 0 (k : ℝ)) + (volume.restrict S) (Ioc (k : ℝ) ((k : ℝ) + 1)) :=
              measure_union hdisj measurableSet_Ioc
        _ = volume (Icc 0 (k : ℝ) ∩ S) + volume (Ioc (k : ℝ) ((k : ℝ) + 1) ∩ S) := by
              rw [Measure.restrict_apply measurableSet_Icc, Measure.restrict_apply measurableSet_Ioc]
        _ = volume (S ∩ Icc 0 (k : ℝ)) + volume (S ∩ Ioc (k : ℝ) ((k : ℝ) + 1)) := by
              rw [Set.inter_comm S, Set.inter_comm S]
    -- translation: vol(S ∩ Ioc k (k+1)) = vol(S ∩ Ioc 0 1)
    have htr : volume (S ∩ Ioc (k : ℝ) ((k : ℝ) + 1)) = volume (S ∩ Ioc 0 1) := by
      have hpre : (fun h => (k : ℝ) + h) ⁻¹' (S ∩ Ioc (k : ℝ) ((k : ℝ) + 1)) = S ∩ Ioc 0 1 := by
        ext y
        simp only [Set.mem_preimage, Set.mem_inter_iff, Set.mem_Ioc]
        constructor
        · rintro ⟨hyS, h1, h2⟩
          refine ⟨?_, by linarith, by linarith⟩
          have : y + (k : ℝ) ∈ S := by rw [add_comm]; exact hyS
          exact (good_add_natCast θ E k y).mp this
        · rintro ⟨hyS, h1, h2⟩
          refine ⟨?_, by linarith, by linarith⟩
          rw [hS] at hyS ⊢
          have : y + (k : ℝ) ∈ Good θ E := (good_add_natCast θ E k y).mpr hyS
          rw [add_comm] at this; exact this
      rw [← hpre, measure_preimage_add]
    -- null: vol(S ∩ Ioc 0 1) = vol(S ∩ Icc 0 1)
    have hnull : volume (S ∩ Ioc 0 1) = volume (S ∩ Icc (0 : ℝ) 1) := by
      apply le_antisymm
      · exact measure_mono (Set.inter_subset_inter_right _ Set.Ioc_subset_Icc_self)
      · have hsub : S ∩ Icc (0 : ℝ) 1 ⊆ (S ∩ Ioc 0 1) ∪ {0} := by
          rintro x ⟨hxS, hx0, hx1⟩
          rcases eq_or_lt_of_le hx0 with h | h
          · exact Or.inr (by simp [h])
          · exact Or.inl ⟨hxS, h, hx1⟩
        calc volume (S ∩ Icc (0:ℝ) 1) ≤ volume ((S ∩ Ioc 0 1) ∪ {(0:ℝ)}) := measure_mono hsub
          _ ≤ volume (S ∩ Ioc 0 1) + volume ({(0:ℝ)}) := measure_union_le _ _
          _ = volume (S ∩ Ioc 0 1) := by simp
    rw [hcast, hrec, htr, hnull, ih, succ_nsmul]

/-- **DILATION INVARIANCE of the good-set measure.**  `muGood θ (c·E) = muGood θ E` for any
positive integer dilation `c`.  Combines `good_dilate` (set relation), `volume_preimage_mul_left`
(Lebesgue scaling by `c`), and `muGood_fold` (the `c`-fold cover of the period).  This is the
covering-side of the S155 dilation-invariance correction, now at the MEASURE level. -/
theorem muGood_dilate (θ : ℝ) (E : Finset ℤ) {c : ℤ} (hc : 0 < c) :
    muGood θ (E.image (fun e => c * e)) = muGood θ E := by
  have hcR : (0 : ℝ) < (c : ℝ) := by exact_mod_cast hc
  have hnat : (c : ℝ) = ((c.toNat : ℕ) : ℝ) := by exact_mod_cast (Int.toNat_of_nonneg hc.le).symm
  unfold muGood
  rw [good_dilate]
  -- (·c)⁻¹' (Good) ∩ [0,1] = (·c)⁻¹' (Good ∩ [0,c])
  have hL2 : (fun x => (c : ℝ) * x) ⁻¹' Good θ E ∩ Icc 0 1
           = (fun x => (c : ℝ) * x) ⁻¹' (Good θ E ∩ Icc 0 (c : ℝ)) := by
    ext x
    simp only [Set.mem_inter_iff, Set.mem_preimage, Set.mem_Icc]
    constructor
    · rintro ⟨hxS, hx0, hx1⟩
      exact ⟨hxS, mul_nonneg hcR.le hx0, by nlinarith⟩
    · rintro ⟨hxS, h0, h1⟩
      refine ⟨hxS, ?_, ?_⟩
      · exact le_of_mul_le_mul_left (by rwa [mul_zero]) hcR
      · exact le_of_mul_le_mul_left (by rwa [mul_one]) hcR
  rw [hL2, Real.volume_preimage_mul_left (ne_of_gt hcR) _]
  -- ofReal |c⁻¹| * volume(Good ∩ [0,c]) = volume(Good ∩ [0,1])
  have hfold : volume (Good θ E ∩ Icc 0 (c : ℝ))
             = (c.toNat : ℕ) • volume (Good θ E ∩ Icc 0 1) := by
    rw [hnat]; exact muGood_fold θ E c.toNat
  rw [hfold, nsmul_eq_mul]
  have hc0 : |(c : ℝ)⁻¹| = (c : ℝ)⁻¹ := abs_of_pos (by positivity)
  have hcast : ((c.toNat : ℕ) : ℝ≥0∞) = ENNReal.ofReal (c : ℝ) := by
    rw [hnat, ENNReal.ofReal_natCast]
  rw [hc0, hcast, ← mul_assoc, ← ENNReal.ofReal_mul (by positivity),
      inv_mul_cancel₀ (ne_of_gt hcR), ENNReal.ofReal_one, one_mul]

/-- **AFFINE (dilation ∘ translation) invariance — the WLOG-normalize primitive.**
`muGood θ (c·E + m) = muGood θ E` for any positive integer dilation `c` and shift `m`.  So `muGood`
(equivalently the degree-3 floor `D3`) depends ONLY on the affine-dilation class of the speed set;
every family reduces to its PRIMITIVE representative `(E − min E)/gcd`.  This is the reduction the
whole k=11 covering tail (opus-S155–S164) and the good-period capstone (dilation-invariant `j*`) rest
on.  Composes `muGood_dilate` (this file) with `muGood_translate` (`LRCTailDiameter`). -/
theorem muGood_affine (θ : ℝ) (E : Finset ℤ) {c : ℤ} (hc : 0 < c) (m : ℤ) :
    muGood θ (E.image (fun e => c * e + m)) = muGood θ E := by
  have himg : E.image (fun e => c * e + m)
            = (E.image (fun e => c * e)).image (fun e => e + m) := by
    ext y
    simp only [Finset.mem_image]
    constructor
    · rintro ⟨e, he, rfl⟩; exact ⟨c * e, ⟨e, he, rfl⟩, rfl⟩
    · rintro ⟨x, ⟨e, he, rfl⟩, rfl⟩; exact ⟨e, he, rfl⟩
  have htrans : (fun e : ℤ => e + m) = (fun e : ℤ => e - (-m)) := by funext e; ring
  rw [himg, htrans, muGood_translate θ (E.image (fun e => c * e)) (-m)]
  exact muGood_dilate θ E hc

-- kernel-purity audit (propext / Classical.choice / Quot.sound only)
#print axioms muGood_dilate
#print axioms muGood_affine

end TailDiameter
end LonelyRunner
