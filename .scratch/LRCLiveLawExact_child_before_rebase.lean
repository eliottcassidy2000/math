/-
  TournamentH7.LRCLiveLawExact -- exact resonant live law for the
  canonical family (HYP-7280, completion of THM-991).

  `LRCLiveLaw` proves that a live multiplier can occur only when `14 | q`
  and constructs the six unit multiples when `q = 14m`.  This file proves
  the missing converse.  At a resonant modulus, the fourteen residues

      k * p mod (14m),       0 <= k < 14,

  occupy distinct length-`m` blocks.  Hence their block map is a
  permutation.  Reading the unique residue in each successive block and
  applying `live_gap` makes the within-block offsets monotone.  The wrap
  gap closes the chain, so every offset is zero.  In particular `m | p`;
  the remaining residue-14 classification is the fixed six-unit fact.

  Kernel-pure: no `sorry`, no `native_decide`, and no new axioms.
-/
import TournamentH7.LRCLiveLaw

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- The block containing the residue `k*p mod (14m)`. -/
def liveBlockMap (m p : ℕ) (hm : 0 < m) (k : Fin 14) : Fin 14 :=
  ⟨((k : ℕ) * p % (14 * m)) / m, by
    have hr : ((k : ℕ) * p) % (14 * m) < 14 * m :=
      Nat.mod_lt _ (by omega)
    exact (Nat.div_lt_iff_lt_mul hm).2 (by simpa [Nat.mul_comm] using hr)⟩

/-- At a live resonant multiplier, distinct source indices occupy distinct
length-`m` blocks. -/
theorem liveBlockMap_injective (m p : ℕ) (hm : 0 < m)
    (hlive : bandCount canonical (14 * m) p = 0) :
    Function.Injective (liveBlockMap m p hm) := by
  intro k l hblock
  by_contra hne
  have hq : 0 < 14 * m := by omega
  have hk13 : (k : ℕ) ≤ 13 := by omega
  have hl13 : (l : ℕ) ≤ 13 := by omega
  have hdiv : (((k : ℕ) * p) % (14 * m)) / m =
      (((l : ℕ) * p) % (14 * m)) / m := by
    exact Fin.ext_iff.mp hblock
  have hkdecomp := Nat.div_add_mod (((k : ℕ) * p) % (14 * m)) m
  have hldecomp := Nat.div_add_mod (((l : ℕ) * p) % (14 * m)) m
  have hkmod : (((k : ℕ) * p) % (14 * m)) % m < m := Nat.mod_lt _ hm
  have hlmod : (((l : ℕ) * p) % (14 * m)) % m < m := Nat.mod_lt _ hm
  rcases lt_trichotomy (k : ℕ) (l : ℕ) with hkl | hkl | hlk
  · obtain ⟨hsep, -, -, -⟩ :=
      live_gap (14 * m) p hq hlive (k : ℕ) (l : ℕ) hkl hl13
    rcases hsep with hsep | hsep <;> omega
  · exact hne (Fin.ext hkl)
  · obtain ⟨hsep, -, -, -⟩ :=
      live_gap (14 * m) p hq hlive (l : ℕ) (k : ℕ) hlk hk13
    rcases hsep with hsep | hsep <;> omega

/-- The resonant block map is a permutation of the fourteen blocks. -/
noncomputable def liveBlockEquiv (m p : ℕ) (hm : 0 < m)
    (hlive : bandCount canonical (14 * m) p = 0) : Fin 14 ≃ Fin 14 :=
  Equiv.ofBijective (liveBlockMap m p hm)
    ((Finite.injective_iff_bijective).mp (liveBlockMap_injective m p hm hlive))

/-- The residue selected by the inverse block permutation. -/
noncomputable def liveResidueInBlock (m p : ℕ) (hm : 0 < m)
    (hlive : bandCount canonical (14 * m) p = 0) (b : Fin 14) : ℕ :=
  (((liveBlockEquiv m p hm hlive).symm b : ℕ) * p) % (14 * m)

/-- The selected residue really lies in its named block. -/
theorem liveResidueInBlock_div (m p : ℕ) (hm : 0 < m)
    (hlive : bandCount canonical (14 * m) p = 0) (b : Fin 14) :
    liveResidueInBlock m p hm hlive b / m = (b : ℕ) := by
  have h := (liveBlockEquiv m p hm hlive).apply_symm_apply b
  exact Fin.ext_iff.mp h

/-- The selected within-block offset. -/
noncomputable def liveBlockOffset (m p : ℕ) (hm : 0 < m)
    (hlive : bandCount canonical (14 * m) p = 0) (b : Fin 14) : ℕ :=
  liveResidueInBlock m p hm hlive b % m

/-- Quotient and offset reconstruct the selected residue. -/
theorem liveResidueInBlock_eq (m p : ℕ) (hm : 0 < m)
    (hlive : bandCount canonical (14 * m) p = 0) (b : Fin 14) :
    liveResidueInBlock m p hm hlive b =
      (b : ℕ) * m + liveBlockOffset m p hm hlive b := by
  have h := Nat.div_add_mod (liveResidueInBlock m p hm hlive b) m
  rw [liveResidueInBlock_div m p hm hlive b] at h
  simpa [liveBlockOffset, Nat.mul_comm] using h.symm

/-- If two live resonant residues are numerically ordered, both their gap and
their circular co-gap are at least one block.  This symmetric wrapper around
`live_gap` is convenient when the source-index order is unknown. -/
theorem live_residue_gap_of_lt (m p : ℕ) (hm : 0 < m)
    (hlive : bandCount canonical (14 * m) p = 0) (k l : Fin 14)
    (hres : ((k : ℕ) * p) % (14 * m) < ((l : ℕ) * p) % (14 * m)) :
    m ≤ ((l : ℕ) * p) % (14 * m) - ((k : ℕ) * p) % (14 * m) ∧
    m ≤ 14 * m -
      (((l : ℕ) * p) % (14 * m) - ((k : ℕ) * p) % (14 * m)) := by
  have hq : 0 < 14 * m := by omega
  have hc : (14 * m + 13) / 14 = m := by omega
  rcases lt_trichotomy (k : ℕ) (l : ℕ) with hkl | hkl | hlk
  · obtain ⟨hgap, -, hco, -⟩ :=
      live_gap (14 * m) p hq hlive (k : ℕ) (l : ℕ) hkl (by omega)
    rw [hc] at hgap hco
    rcases hgap with hgap | hgap
    · exact ⟨hgap, hco⟩
    · omega
  · subst l
    omega
  · obtain ⟨hgap, -, -, hco⟩ :=
      live_gap (14 * m) p hq hlive (l : ℕ) (k : ℕ) hlk (by omega)
    rw [hc] at hgap hco
    rcases hgap with hgap | hgap
    · omega
    · exact ⟨hgap, hco⟩

/-- Successive block offsets are monotone. -/
theorem liveBlockOffset_mono (m p : ℕ) (hm : 0 < m)
    (hlive : bandCount canonical (14 * m) p = 0) (b : Fin 13) :
    liveBlockOffset m p hm hlive ⟨(b : ℕ), by omega⟩ ≤
      liveBlockOffset m p hm hlive ⟨(b : ℕ) + 1, by omega⟩ := by
  let b₀ : Fin 14 := ⟨(b : ℕ), by omega⟩
  let b₁ : Fin 14 := ⟨(b : ℕ) + 1, by omega⟩
  have hδ0 : liveBlockOffset m p hm hlive b₀ < m := Nat.mod_lt _ hm
  have hδ1 : liveBlockOffset m p hm hlive b₁ < m := Nat.mod_lt _ hm
  have hr₀ := liveResidueInBlock_eq m p hm hlive b₀
  have hr₁ := liveResidueInBlock_eq m p hm hlive b₁
  have hres : liveResidueInBlock m p hm hlive b₀ <
      liveResidueInBlock m p hm hlive b₁ := by
    dsimp [b₀, b₁] at hr₀ hr₁ ⊢
    omega
  have hgap := live_residue_gap_of_lt m p hm hlive
    ((liveBlockEquiv m p hm hlive).symm b₀)
    ((liveBlockEquiv m p hm hlive).symm b₁) hres
  change m ≤ liveResidueInBlock m p hm hlive b₁ -
      liveResidueInBlock m p hm hlive b₀ ∧
    m ≤ 14 * m - (liveResidueInBlock m p hm hlive b₁ -
      liveResidueInBlock m p hm hlive b₀) at hgap
  dsimp [b₀, b₁] at hr₀ hr₁ ⊢
  omega

/-- The wrap gap closes the offset chain. -/
theorem liveBlockOffset_last_le_first (m p : ℕ) (hm : 0 < m)
    (hlive : bandCount canonical (14 * m) p = 0) :
    liveBlockOffset m p hm hlive (13 : Fin 14) ≤
      liveBlockOffset m p hm hlive (0 : Fin 14) := by
  have hδ0 : liveBlockOffset m p hm hlive (0 : Fin 14) < m := Nat.mod_lt _ hm
  have hδ13 : liveBlockOffset m p hm hlive (13 : Fin 14) < m := Nat.mod_lt _ hm
  have hr₀ := liveResidueInBlock_eq m p hm hlive (0 : Fin 14)
  have hr₁₃ := liveResidueInBlock_eq m p hm hlive (13 : Fin 14)
  have hres : liveResidueInBlock m p hm hlive (0 : Fin 14) <
      liveResidueInBlock m p hm hlive (13 : Fin 14) := by
    norm_num at hr₀ hr₁₃ ⊢
    omega
  have hgap := (live_residue_gap_of_lt m p hm hlive
    ((liveBlockEquiv m p hm hlive).symm (0 : Fin 14))
    ((liveBlockEquiv m p hm hlive).symm (13 : Fin 14)) hres).2
  change m ≤ 14 * m -
    (liveResidueInBlock m p hm hlive (13 : Fin 14) -
      liveResidueInBlock m p hm hlive (0 : Fin 14)) at hgap
  norm_num at hr₀ hr₁₃ ⊢
  omega

/-- Block zero contains the residue zero, so its offset vanishes. -/
theorem liveBlockOffset_zero (m p : ℕ) (hm : 0 < m)
    (hlive : bandCount canonical (14 * m) p = 0) :
    liveBlockOffset m p hm hlive (0 : Fin 14) = 0 := by
  let e := liveBlockEquiv m p hm hlive
  have he0 : e (0 : Fin 14) = 0 := by
    apply Fin.ext
    simp [e, liveBlockEquiv, liveBlockMap]
  have hes0 : e.symm (0 : Fin 14) = 0 := by
    apply e.injective
    simpa [he0]
  simp [liveBlockOffset, liveResidueInBlock, e, hes0]

/-- Every block offset is zero: the fourteen live residues form the exact
`m`-net, not merely one point in each block. -/
theorem liveBlockOffset_eq_zero (m p : ℕ) (hm : 0 < m)
    (hlive : bandCount canonical (14 * m) p = 0) (b : Fin 14) :
    liveBlockOffset m p hm hlive b = 0 := by
  have h0 := liveBlockOffset_zero m p hm hlive
  have h01 := liveBlockOffset_mono m p hm hlive (0 : Fin 13)
  have h12 := liveBlockOffset_mono m p hm hlive (1 : Fin 13)
  have h23 := liveBlockOffset_mono m p hm hlive (2 : Fin 13)
  have h34 := liveBlockOffset_mono m p hm hlive (3 : Fin 13)
  have h45 := liveBlockOffset_mono m p hm hlive (4 : Fin 13)
  have h56 := liveBlockOffset_mono m p hm hlive (5 : Fin 13)
  have h67 := liveBlockOffset_mono m p hm hlive (6 : Fin 13)
  have h78 := liveBlockOffset_mono m p hm hlive (7 : Fin 13)
  have h89 := liveBlockOffset_mono m p hm hlive (8 : Fin 13)
  have h910 := liveBlockOffset_mono m p hm hlive (9 : Fin 13)
  have h1011 := liveBlockOffset_mono m p hm hlive (10 : Fin 13)
  have h1112 := liveBlockOffset_mono m p hm hlive (11 : Fin 13)
  have h1213 := liveBlockOffset_mono m p hm hlive (12 : Fin 13)
  have h130 := liveBlockOffset_last_le_first m p hm hlive
  fin_cases b <;> omega

/-- Exact net extraction: the unique live residue in block `b` is `b*m`. -/
theorem liveResidueInBlock_exact (m p : ℕ) (hm : 0 < m)
    (hlive : bandCount canonical (14 * m) p = 0) (b : Fin 14) :
    liveResidueInBlock m p hm hlive b = (b : ℕ) * m := by
  rw [liveResidueInBlock_eq m p hm hlive b,
    liveBlockOffset_eq_zero m p hm hlive b, Nat.add_zero]

/-- Every source residue is exactly the left endpoint of its block. -/
theorem live_residue_eq_block_mul (m p : ℕ) (hm : 0 < m)
    (hlive : bandCount canonical (14 * m) p = 0) (k : Fin 14) :
    ((k : ℕ) * p) % (14 * m) =
      ((liveBlockMap m p hm k : Fin 14) : ℕ) * m := by
  let e := liveBlockEquiv m p hm hlive
  have hexact := liveResidueInBlock_exact m p hm hlive (e k)
  have hinv : e.symm (e k) = k := e.symm_apply_apply k
  simpa [liveResidueInBlock, e, liveBlockEquiv, hinv] using hexact

/-- The fixed six unit residues modulo fourteen. -/
def canonicalLiveUnits : Finset ℕ := {1, 3, 5, 9, 11, 13}

/-- Fixed residue-14 dichotomy.  A block label is either one of the six units,
or some canonical speed annihilates it modulo fourteen.  This is the only
bounded finite check in the proof. -/
theorem mem_canonicalLiveUnits_or_annihilator (o : Fin 14) :
    (o : ℕ) ∈ canonicalLiveUnits ∨
      ∃ v : Fin 13, ((((v : ℕ) + 1) * (o : ℕ)) % 14) = 0 := by
  fin_cases o <;> decide

/-- **Pointwise exact resonant classification.**  Every nonzero live
multiplier below `14m` is one of the six unit multiples of `m`. -/
theorem live_multiplier_eq_unit_mul (m p : ℕ) (hm : 0 < m)
    (hp : 0 < p) (hpq : p < 14 * m)
    (hlive : bandCount canonical (14 * m) p = 0) :
    ∃ o ∈ canonicalLiveUnits, p = o * m := by
  let o : Fin 14 := liveBlockMap m p hm (1 : Fin 14)
  have hres := live_residue_eq_block_mul m p hm hlive (1 : Fin 14)
  have hpmod : p % (14 * m) = p := Nat.mod_eq_of_lt hpq
  have hpform : p = (o : ℕ) * m := by
    simpa [o, hpmod] using hres
  rcases mem_canonicalLiveUnits_or_annihilator o with hunit | ⟨v, hv⟩
  · exact ⟨(o : ℕ), hunit, hpform⟩
  · have hsafe := (live_safe (14 * m) p hlive ((v : ℕ) + 1)
        (by omega) (by omega)).1
    have hscale : (((v : ℕ) + 1) * p) % (14 * m) =
        ((((v : ℕ) + 1) * (o : ℕ)) % 14) * m := by
      rw [hpform]
      rw [show ((v : ℕ) + 1) * ((o : ℕ) * m) =
          (((v : ℕ) + 1) * (o : ℕ)) * m by ring]
      exact Nat.mul_mod_mul_right _ _ _
    rw [hscale, hv] at hsafe
    simp at hsafe
    omega

/-- Exact equality of finite live-multiplier sets at a resonant modulus. -/
theorem canonical_liveMultipliers_eq (m : ℕ) (hm : 0 < m) :
    (Finset.Ioo 0 (14 * m)).filter
        (fun p => bandCount canonical (14 * m) p = 0) =
      canonicalLiveUnits.image (fun o => o * m) := by
  ext p
  constructor
  · intro hp
    rw [Finset.mem_filter, Finset.mem_Ioo] at hp
    obtain ⟨o, ho, rfl⟩ :=
      live_multiplier_eq_unit_mul m p hm hp.1.1 hp.1.2 hp.2
    exact Finset.mem_image.mpr ⟨o, ho, rfl⟩
  · intro hp
    obtain ⟨o, ho, rfl⟩ := Finset.mem_image.mp hp
    rw [Finset.mem_filter, Finset.mem_Ioo]
    have ho_cases : o = 1 ∨ o = 3 ∨ o = 5 ∨ o = 9 ∨ o = 11 ∨ o = 13 := by
      simpa [canonicalLiveUnits] using ho
    rcases ho_cases with rfl | rfl | rfl | rfl | rfl | rfl <;>
      exact ⟨⟨by omega, by omega⟩,
        resonant_live m _ hm (by norm_num) (by decide)⟩

/-- **Exact resonant live count.** -/
theorem liveCount_eq_six_of_dvd (m : ℕ) (hm : 0 < m) :
    liveCount canonical (14 * m) = 6 := by
  unfold liveCount
  rw [canonical_liveMultipliers_eq m hm]
  rw [Finset.card_image_of_injOn]
  · decide
  · intro a _ b _ hab
    exact Nat.eq_of_mul_eq_mul_right hm hab

/-- The positive-modulus all-`q` live law. -/
theorem canonical_liveCount (q : ℕ) (hq : 0 < q) :
    liveCount canonical q = if 14 ∣ q then 6 else 0 := by
  by_cases hdvd : 14 ∣ q
  · rw [if_pos hdvd]
    obtain ⟨m, rfl⟩ := hdvd
    exact liveCount_eq_six_of_dvd m (by omega)
  · rw [if_neg hdvd]
    exact liveCount_eq_zero_of_not_dvd q hq hdvd

/-- At the excluded zero modulus the live set is empty (whereas `14 | 0`,
so the unqualified conditional formula would deliberately be false). -/
theorem canonical_liveCount_zero : liveCount canonical 0 = 0 := by
  rfl

/-! ## Axiom audit -/
#print axioms liveBlockMap_injective
#print axioms live_residue_gap_of_lt
#print axioms liveBlockOffset_eq_zero
#print axioms live_residue_eq_block_mul
#print axioms live_multiplier_eq_unit_mul
#print axioms canonical_liveMultipliers_eq
#print axioms liveCount_eq_six_of_dvd
#print axioms canonical_liveCount

end LRC14Concrete
end LonelyRunner
