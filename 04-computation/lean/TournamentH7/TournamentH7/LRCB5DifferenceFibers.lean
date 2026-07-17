import Mathlib

/-!
Exact fiber count for the cyclic tooth-index difference map.  For speeds
`g*p` and `g*q` with coprime reduced parts, every class modulo `g*p*q` occurs
exactly `g` times among the two tooth-index sets.  This discharges the
gcd/multiplicity half of the remaining pair-overlap reindexer.
-/

namespace LonelyRunner
namespace LRCB5DifferenceFibers

noncomputable section

/-- Multiplication by `scale`, viewed as an additive map from `ZMod source`
into `ZMod target` when `source * scale = target`. -/
def zmodScaleHom (source scale target : ℕ) (h : source * scale = target) :
    ZMod source →+ ZMod target :=
  ZMod.lift source ⟨zmultiplesHom (ZMod target) (scale : ZMod target), by
    rw [zmultiplesHom_apply]
    rw [zsmul_eq_mul]
    have hcast := congrArg (fun value : ℕ => (value : ZMod target)) h
    simpa using hcast⟩

theorem zmodScaleHom_intCast
    (source scale target : ℕ) (h : source * scale = target) (value : ℤ) :
    zmodScaleHom source scale target h (value : ZMod source) =
      (value * scale : ℤ) := by
  unfold zmodScaleHom
  rw [ZMod.lift_coe]
  simp [zmultiplesHom_apply, zsmul_eq_mul]

/-- Difference of the two cyclic tooth indices. -/
def differenceHom (g p q : ℕ) :
    ZMod (g * p) × ZMod (g * q) →+ ZMod (g * p * q) :=
  (zmodScaleHom (g * p) q (g * p * q) (by ring)).comp
      (AddMonoidHom.fst _ _) -
    (zmodScaleHom (g * q) p (g * p * q) (by ring)).comp
      (AddMonoidHom.snd _ _)

theorem differenceHom_intCast
    (g p q : ℕ) (first second : ℤ) :
    differenceHom g p q
        ((first : ZMod (g * p)), (second : ZMod (g * q))) =
      (q * first - p * second : ℤ) := by
  unfold differenceHom
  simp [zmodScaleHom_intCast]
  ring

theorem differenceHom_surjective
    (g p q : ℕ) (hcoprime : Nat.Coprime p q) :
    Function.Surjective (differenceHom g p q) := by
  intro residue
  obtain ⟨value, rfl⟩ := ZMod.intCast_surjective residue
  let firstCoeff : ℤ := q.gcdA p
  let secondCoeff : ℤ := -q.gcdB p
  refine ⟨((firstCoeff * value : ℤ), (secondCoeff * value : ℤ)), ?_⟩
  rw [differenceHom_intCast]
  have hbezout : (q : ℤ) * firstCoeff - (p : ℤ) * secondCoeff = 1 := by
    dsimp [firstCoeff, secondCoeff]
    calc
      (q : ℤ) * q.gcdA p - (p : ℤ) * -q.gcdB p =
          (q : ℤ) * q.gcdA p + (p : ℤ) * q.gcdB p := by ring
      _ = (q.gcd p : ℕ) := (Nat.gcd_eq_gcd_ab q p).symm
      _ = 1 := by rw [hcoprime.symm.gcd_eq_one]; norm_num
  have hinteger :
      (q : ℤ) * (firstCoeff * value) -
          (p : ℤ) * (secondCoeff * value) = value := by
    calc
      (q : ℤ) * (firstCoeff * value) -
          (p : ℤ) * (secondCoeff * value) =
          ((q : ℤ) * firstCoeff - (p : ℤ) * secondCoeff) * value := by ring
      _ = value := by rw [hbezout, one_mul]
  simpa using
    congrArg (fun integer : ℤ => (integer : ZMod (g * p * q))) hinteger

theorem natCard_zmod (modulus : ℕ) [NeZero modulus] :
    Nat.card (ZMod modulus) = modulus := by
  rw [Nat.card_eq_fintype_card, ZMod.card]

theorem differenceHom_ker_card
    (g p q : ℕ) (hg : 0 < g) (hp : 0 < p) (hq : 0 < q)
    (hcoprime : Nat.Coprime p q) :
    Nat.card (differenceHom g p q).ker = g := by
  letI : NeZero (g * p) := ⟨(Nat.mul_pos hg hp).ne'⟩
  letI : NeZero (g * q) := ⟨(Nat.mul_pos hg hq).ne'⟩
  letI : NeZero (g * p * q) := ⟨(Nat.mul_pos (Nat.mul_pos hg hp) hq).ne'⟩
  let hom := differenceHom g p q
  have hsurjective : Function.Surjective hom :=
    differenceHom_surjective g p q hcoprime
  have hrange : hom.range = ⊤ := AddMonoidHom.range_eq_top.mpr hsurjective
  have hindex : hom.ker.index = g * p * q := by
    rw [AddSubgroup.index_ker, hrange]
    simp
  have hcard := hom.ker.card_mul_index
  rw [hindex, Nat.card_prod, natCard_zmod, natCard_zmod] at hcard
  have hproduct :
      Nat.card hom.ker * (g * p * q) = g * (g * p * q) := by
    calc
      Nat.card hom.ker * (g * p * q) = (g * p) * (g * q) := hcard
      _ = g * (g * p * q) := by ring
  exact Nat.mul_right_cancel (Nat.mul_pos (Nat.mul_pos hg hp) hq) hproduct

/-- Every cyclic difference class occurs exactly `g` times among the
`g*p` by `g*q` tooth-index pairs. -/
theorem sum_weight_differenceHom
    (g p q : ℕ) [NeZero g] [NeZero p] [NeZero q]
    (hg : 0 < g) (hp : 0 < p) (hq : 0 < q)
    (hcoprime : Nat.Coprime p q)
    (weight : ZMod (g * p * q) → ℝ) :
    (∑ pair : ZMod (g * p) × ZMod (g * q),
        weight (differenceHom g p q pair)) =
      (g : ℝ) * ∑ residue : ZMod (g * p * q), weight residue := by
  let hom := differenceHom g p q
  have hsurjective : Function.Surjective hom :=
    differenceHom_surjective g p q hcoprime
  have hkerCard : Nat.card hom.ker = g :=
    differenceHom_ker_card g p q hg hp hq hcoprime
  have hfiberCard (residue : ZMod (g * p * q)) :
      Fintype.card
          {pair : ZMod (g * p) × ZMod (g * q) // hom pair = residue} = g := by
    have hequiv := AddMonoidHom.fiberEquivKerOfSurjective hsurjective residue
    calc
      Fintype.card
          {pair : ZMod (g * p) × ZMod (g * q) // hom pair = residue} =
          Nat.card
            {pair : ZMod (g * p) × ZMod (g * q) // hom pair = residue} := by
        rw [Nat.card_eq_fintype_card]
      _ = Nat.card hom.ker := by
        apply Nat.card_congr
        exact hequiv
      _ = g := hkerCard
  rw [← Fintype.sum_fiberwise hom (fun pair => weight (hom pair))]
  rw [Finset.mul_sum]
  apply Fintype.sum_congr
  intro residue
  calc
    (∑ pair : {pair : ZMod (g * p) × ZMod (g * q) // hom pair = residue},
        weight (hom pair)) =
        ∑ _pair :
            {pair : ZMod (g * p) × ZMod (g * q) // hom pair = residue},
          weight residue := by
      apply Fintype.sum_congr
      intro pair
      rw [pair.property]
    _ = (g : ℝ) * weight residue := by
      rw [Finset.sum_const, Finset.card_univ, hfiberCard, nsmul_eq_mul]

#print axioms differenceHom_surjective
#print axioms differenceHom_ker_card
#print axioms sum_weight_differenceHom

end
end LRCB5DifferenceFibers
end LonelyRunner
