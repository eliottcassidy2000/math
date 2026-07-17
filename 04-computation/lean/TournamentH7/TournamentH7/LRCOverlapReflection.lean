/-
  TournamentH7.LRCOverlapReflection

  Reflection symmetry for the colored overlap movie.  On a bad runner the
  canonical witness at `q - p` is the reflected integer `v - n`; hence every
  overlap determinant changes sign while its support and absolute mass are
  preserved.  In particular, a genuinely colored event cannot be fixed by
  reflection, even at the half multiplier of an even modulus.  Colored
  fixed-triple and fixed-triangle activity therefore comes in exact pairs for
  every positive modulus.

  Tournament-analysis audit: the vertices are multiplier events, the pair
  observable is `overlapDet`, and reflection is the gauge involution.  It
  preserves the zero/nonzero support graph and absolute edge weights but
  reverses every signed color.  Quotienting to a static runner tournament
  destroys this orbit structure and identifies a sparse relation with its
  negative without recording that the two events are one projective relation.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.BucketBalance
import TournamentH7.LRCParityPairing
import TournamentH7.LRCOverlapColorFibers
import TournamentH7.LRCSevenOverlapActivity

namespace LonelyRunner
namespace LRC14Concrete

open scoped Classical

/-- Per-runner band membership is exactly reflection invariant on the open
multiplier window. -/
theorem inBand_reflect_iff
    (v : Fin 13 → ℤ) (q p : ℕ) (i : Fin 13)
    (hq : 0 < q) (hp0 : 0 < p) (hpq : p < q) :
    inBand v q (q - p) i ↔ inBand v q p i := by
  constructor
  · intro hreflected
    have hmirror := inBand_mirror v q (q - p) hq (by omega) i hreflected
    simpa [Nat.sub_sub_self hpq.le] using hmirror
  · exact inBand_mirror v q p hq hpq i

/-- Two integer multiples of a positive modulus cannot both lie strictly
inside the lonely-runner radius around the same integer. -/
private theorem close_witness_unique
    (x n₁ n₂ : ℤ) (q : ℕ) (hq : 0 < q)
    (h₁ : 14 * |x - n₁ * (q : ℤ)| < (q : ℤ))
    (h₂ : 14 * |x - n₂ * (q : ℤ)| < (q : ℤ)) :
    n₁ = n₂ := by
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  by_contra hne
  rcases lt_or_gt_of_ne hne with hlt | hgt
  · have hgap : (q : ℤ) ≤ (n₂ - n₁) * (q : ℤ) := by
      have : (1 : ℤ) ≤ n₂ - n₁ := by omega
      nlinarith
    have hleft : 14 * (x - n₁ * (q : ℤ)) < (q : ℤ) := by
      have habs := le_abs_self (x - n₁ * (q : ℤ))
      nlinarith
    have hright : -(q : ℤ) < 14 * (x - n₂ * (q : ℤ)) := by
      have habs := neg_abs_le (x - n₂ * (q : ℤ))
      nlinarith
    nlinarith
  · have hgap : (q : ℤ) ≤ (n₁ - n₂) * (q : ℤ) := by
      have : (1 : ℤ) ≤ n₁ - n₂ := by omega
      nlinarith
    have hleft : 14 * (x - n₂ * (q : ℤ)) < (q : ℤ) := by
      have habs := le_abs_self (x - n₂ * (q : ℤ))
      nlinarith
    have hright : -(q : ℤ) < 14 * (x - n₁ * (q : ℤ)) := by
      have habs := neg_abs_le (x - n₁ * (q : ℤ))
      nlinarith
    nlinarith

/-- On a bad runner the canonical nearest-integer witness reflects by
`n ↦ v - n`.  The badness hypothesis is essential: a safe-band residue need
not use this branch of `failWitness`. -/
theorem failWitness_reflect_of_bad
    (v : Fin 13 → ℤ) (q p : ℕ) (i : Fin 13)
    (hq : 0 < q) (hp0 : 0 < p) (hpq : p < q)
    (hbad : ¬ inBand v q p i) :
    failWitness v q (q - p) i = v i - failWitness v q p i := by
  have hbadReflected : ¬ inBand v q (q - p) i := by
    intro hreflected
    exact hbad ((inBand_reflect_iff v q p i hq hp0 hpq).mp hreflected)
  have hcanonical := bad_at_witness v q (q - p) i hq hbadReflected
  have horiginal := bad_at_witness v q p i hq hbad
  have hcast : ((q - p : ℕ) : ℤ) = (q : ℤ) - (p : ℤ) := by
    exact_mod_cast Nat.cast_sub hpq.le
  have hreflectedCandidate :
      14 * |v i * ((q - p : ℕ) : ℤ) -
          (v i - failWitness v q p i) * (q : ℤ)| < (q : ℤ) := by
    have hid :
        v i * ((q - p : ℕ) : ℤ) -
            (v i - failWitness v q p i) * (q : ℤ) =
          -(v i * (p : ℤ) - failWitness v q p i * (q : ℤ)) := by
      rw [hcast]
      ring
    rw [hid, abs_neg]
    exact horiginal
  exact close_witness_unique
    (v i * ((q - p : ℕ) : ℤ))
    (failWitness v q (q - p) i)
    (v i - failWitness v q p i) q hq hcanonical hreflectedCandidate

/-- Every signed overlap color reverses under multiplier reflection. -/
theorem overlapDet_reflect_of_bad
    (v : Fin 13 → ℤ) (q p : ℕ) (i j : Fin 13)
    (hq : 0 < q) (hp0 : 0 < p) (hpq : p < q)
    (hbadi : ¬ inBand v q p i) (hbadj : ¬ inBand v q p j) :
    overlapDet v q (q - p) i j = -overlapDet v q p i j := by
  unfold overlapDet
  rw [failWitness_reflect_of_bad v q p j hq hp0 hpq hbadj,
    failWitness_reflect_of_bad v q p i hq hp0 hpq hbadi]
  ring

/-- Reflection preserves absolute determinant activity. -/
theorem abs_overlapDet_reflect_of_bad
    (v : Fin 13 → ℤ) (q p : ℕ) (i j : Fin 13)
    (hq : 0 < q) (hp0 : 0 < p) (hpq : p < q)
    (hbadi : ¬ inBand v q p i) (hbadj : ¬ inBand v q p j) :
    |overlapDet v q (q - p) i j| = |overlapDet v q p i j| := by
  rw [overlapDet_reflect_of_bad v q p i j hq hp0 hpq hbadi hbadj, abs_neg]

/-- A sparse Plücker coefficient vector is also negated by reflection.  Thus
the two reflected events carry one projective integer relation, not two
independent relations. -/
theorem overlapTripleCoeff_reflect_of_bad
    (v : Fin 13 → ℤ) (q p : ℕ) (i j k index : Fin 13)
    (hq : 0 < q) (hp0 : 0 < p) (hpq : p < q)
    (hbadi : ¬ inBand v q p i)
    (hbadj : ¬ inBand v q p j)
    (hbadk : ¬ inBand v q p k) :
    overlapTripleCoeff v q (q - p) i j k index =
      -overlapTripleCoeff v q p i j k index := by
  have hjk := overlapDet_reflect_of_bad v q p j k hq hp0 hpq hbadj hbadk
  have hik := overlapDet_reflect_of_bad v q p i k hq hp0 hpq hbadi hbadk
  have hij := overlapDet_reflect_of_bad v q p i j hq hp0 hpq hbadi hbadj
  unfold overlapTripleCoeff
  rw [hjk, hik, hij]
  have hite_neg (P : Prop) [Decidable P] (a : ℤ) :
      (if P then -a else 0) = -(if P then a else 0) := by
    by_cases hP : P <;> simp [hP]
  rw [hite_neg (index = i) (overlapDet v q p j k),
    hite_neg (index = j) (-overlapDet v q p i k),
    hite_neg (index = k) (overlapDet v q p i j)]
  ring

/-- The only possible fixed multiplier of reflection carries no nonzero
overlap color. -/
theorem overlapDet_half_eq_zero_of_bad
    (v : Fin 13 → ℤ) (m : ℕ) (hm : 0 < m) (i j : Fin 13)
    (hbadi : ¬ inBand v (m + m) m i)
    (hbadj : ¬ inBand v (m + m) m j) :
    overlapDet v (m + m) m i j = 0 := by
  have hneg := overlapDet_reflect_of_bad v (m + m) m i j
    (by omega) hm (by omega) hbadi hbadj
  have hfixed : m + m - m = m := by omega
  rw [hfixed] at hneg
  omega

/-- Reflection sends the signed color fiber `d` into the signed color fiber
`-d`. -/
theorem reflect_mem_overlapColorFiber
    (v : Fin 13 → ℤ) (q : ℕ) (i j : Fin 13) (d : ℤ)
    (hq : 0 < q) (p : ℕ)
    (hp : p ∈ overlapColorFiber v q i j d) :
    q - p ∈ overlapColorFiber v q i j (-d) := by
  unfold overlapColorFiber at hp ⊢
  rw [Finset.mem_filter] at hp ⊢
  obtain ⟨hpWindow, hbadi, hbadj, hcolor⟩ := hp
  have hp0 : 0 < p := (Finset.mem_Ioo.mp hpWindow).1
  have hpq : p < q := (Finset.mem_Ioo.mp hpWindow).2
  have hreflectWindow : q - p ∈ Finset.Ioo 0 q := by
    rw [Finset.mem_Ioo]
    omega
  have hbadiReflected : ¬ inBand v q (q - p) i := by
    intro h
    exact hbadi ((inBand_reflect_iff v q p i hq hp0 hpq).mp h)
  have hbadjReflected : ¬ inBand v q (q - p) j := by
    intro h
    exact hbadj ((inBand_reflect_iff v q p j hq hp0 hpq).mp h)
  refine ⟨hreflectWindow, hbadiReflected, hbadjReflected, ?_⟩
  rw [overlapDet_reflect_of_bad v q p i j hq hp0 hpq hbadi hbadj,
    hcolor]

/-- Opposite signed color fibers have exactly the same activity. -/
theorem overlapColorFiber_card_neg
    (v : Fin 13 → ℤ) (q : ℕ) (i j : Fin 13) (d : ℤ)
    (hq : 0 < q) :
    (overlapColorFiber v q i j d).card =
      (overlapColorFiber v q i j (-d)).card := by
  apply Finset.card_bij (fun p _hp => q - p)
  · intro p hp
    exact reflect_mem_overlapColorFiber v q i j d hq p hp
  · intro p₁ hp₁ p₂ hp₂ heq
    unfold overlapColorFiber at hp₁ hp₂
    have hp₁q := (Finset.mem_Ioo.mp (Finset.mem_filter.mp hp₁).1).2
    have hp₂q := (Finset.mem_Ioo.mp (Finset.mem_filter.mp hp₂).1).2
    omega
  · intro p hp
    refine ⟨q - p, ?_, ?_⟩
    · have hmem := reflect_mem_overlapColorFiber v q i j (-d) hq p hp
      simpa using hmem
    · unfold overlapColorFiber at hp
      have hpq := (Finset.mem_Ioo.mp (Finset.mem_filter.mp hp).1).2
      omega

/-- Combining reflection with the coprime fixed-color bound: every nonzero
absolute color occurs either not at all or as exactly one reflected pair. -/
theorem overlapColorFiber_union_neg_card_eq_zero_or_two_of_coprime
    (v : Fin 13 → ℤ) (q : ℕ) (i j : Fin 13) (d : ℤ)
    (hq : 0 < q) (hvi : 0 < v i) (hvj : 0 < v j)
    (hwin : (q : ℤ) ≤ 7 * v i)
    (hcop : Int.gcd (v i) (v j) = 1) (hd : d ≠ 0) :
    let colors := overlapColorFiber v q i j d ∪
      overlapColorFiber v q i j (-d)
    colors.card = 0 ∨ colors.card = 2 := by
  dsimp only
  have hcard := overlapColorFiber_card_le_one_of_gcd_eq_one
    v q i j d hq hvi hvj hwin hcop
  have heq := overlapColorFiber_card_neg v q i j d hq
  have hdisjoint : Disjoint (overlapColorFiber v q i j d)
      (overlapColorFiber v q i j (-d)) := by
    rw [Finset.disjoint_left]
    intro p hp hpn
    have hcolor := (Finset.mem_filter.mp hp).2.2.2
    have hcolorNeg := (Finset.mem_filter.mp hpn).2.2.2
    have : d = -d := hcolor.symm.trans hcolorNeg
    have : d = 0 := by omega
    exact hd this
  rw [Finset.card_union_of_disjoint hdisjoint, ← heq]
  omega

/-- Every colored fixed-triple event lies in a two-element reflection orbit,
for odd and even moduli alike. -/
theorem coloredBadTripleMultipliers_card_even
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    (σ : Equiv.Perm (Fin 13)) (left right top : Fin 13) :
    Even (LRC14Grand.coloredBadTripleMultipliers
      v q σ left right top).card := by
  let s := LRC14Grand.coloredBadTripleMultipliers v q σ left right top
  let reflect : ℕ → ℕ := fun p => q - p
  apply Tournament.BucketBalance.even_card_of_fixedPointFree_involutiveOn
    s reflect
  · intro p hp
    unfold s LRC14Grand.coloredBadTripleMultipliers at hp ⊢
    rw [Finset.mem_filter] at hp ⊢
    obtain ⟨htriple, hcolor⟩ := hp
    unfold LRC14Grand.badTripleMultipliers at htriple ⊢
    rw [Finset.mem_filter] at htriple ⊢
    obtain ⟨hpWindow, hleft, hright, htop⟩ := htriple
    have hp0 : 0 < p := (Finset.mem_Ioo.mp hpWindow).1
    have hpq : p < q := (Finset.mem_Ioo.mp hpWindow).2
    have hreflectWindow : q - p ∈ Finset.Ioo 0 q := by
      rw [Finset.mem_Ioo]
      omega
    have hleftReflected : ¬ inBand v q (q - p) (σ left) := by
      intro h
      exact hleft ((inBand_reflect_iff v q p (σ left) hq hp0 hpq).mp h)
    have hrightReflected : ¬ inBand v q (q - p) (σ right) := by
      intro h
      exact hright ((inBand_reflect_iff v q p (σ right) hq hp0 hpq).mp h)
    have htopReflected : ¬ inBand v q (q - p) (σ top) := by
      intro h
      exact htop ((inBand_reflect_iff v q p (σ top) hq hp0 hpq).mp h)
    refine ⟨⟨hreflectWindow, hleftReflected, hrightReflected, htopReflected⟩, ?_⟩
    rw [overlapDet_reflect_of_bad v q p (σ left) (σ right)
      hq hp0 hpq hleft hright]
    exact neg_ne_zero.mpr hcolor
  · intro p hp
    unfold reflect
    unfold s LRC14Grand.coloredBadTripleMultipliers at hp
    have hpWindow := (Finset.mem_filter.mp hp).1
    unfold LRC14Grand.badTripleMultipliers at hpWindow
    have hpIoo := (Finset.mem_filter.mp hpWindow).1
    rw [Finset.mem_Ioo] at hpIoo
    omega
  · intro p hp hfixed
    unfold reflect at hfixed
    unfold s LRC14Grand.coloredBadTripleMultipliers at hp
    obtain ⟨htriple, hcolor⟩ := Finset.mem_filter.mp hp
    unfold LRC14Grand.badTripleMultipliers at htriple
    obtain ⟨hpWindow, hleft, hright, _htop⟩ := Finset.mem_filter.mp htriple
    have hp0 : 0 < p := (Finset.mem_Ioo.mp hpWindow).1
    have hpq : p < q := (Finset.mem_Ioo.mp hpWindow).2
    have hneg := overlapDet_reflect_of_bad v q p (σ left) (σ right)
      hq hp0 hpq hleft hright
    rw [hfixed] at hneg
    have : overlapDet v q p (σ left) (σ right) = 0 := by omega
    exact hcolor this

/-- The stronger nonzero-lower-triangle events also occur in reflection pairs
for every positive modulus. -/
theorem coloredBadTriangleMultipliers_card_even
    (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    (σ : Equiv.Perm (Fin 13)) (a b c top : Fin 13) :
    Even (LRC14Grand.coloredBadTriangleMultipliers
      v q σ a b c top).card := by
  let s := LRC14Grand.coloredBadTriangleMultipliers v q σ a b c top
  let reflect : ℕ → ℕ := fun p => q - p
  apply Tournament.BucketBalance.even_card_of_fixedPointFree_involutiveOn
    s reflect
  · intro p hp
    unfold s LRC14Grand.coloredBadTriangleMultipliers at hp ⊢
    rw [Finset.mem_filter] at hp ⊢
    obtain ⟨hpWindow, ha, hb, hc, htop, hab, hac, hbc⟩ := hp
    have hp0 : 0 < p := (Finset.mem_Ioo.mp hpWindow).1
    have hpq : p < q := (Finset.mem_Ioo.mp hpWindow).2
    have hreflectWindow : q - p ∈ Finset.Ioo 0 q := by
      rw [Finset.mem_Ioo]
      omega
    have haReflected : ¬ inBand v q (q - p) (σ a) := by
      intro h
      exact ha ((inBand_reflect_iff v q p (σ a) hq hp0 hpq).mp h)
    have hbReflected : ¬ inBand v q (q - p) (σ b) := by
      intro h
      exact hb ((inBand_reflect_iff v q p (σ b) hq hp0 hpq).mp h)
    have hcReflected : ¬ inBand v q (q - p) (σ c) := by
      intro h
      exact hc ((inBand_reflect_iff v q p (σ c) hq hp0 hpq).mp h)
    have htopReflected : ¬ inBand v q (q - p) (σ top) := by
      intro h
      exact htop ((inBand_reflect_iff v q p (σ top) hq hp0 hpq).mp h)
    refine ⟨hreflectWindow, haReflected, hbReflected, hcReflected,
      htopReflected, ?_, ?_, ?_⟩
    · rw [overlapDet_reflect_of_bad v q p (σ a) (σ b)
        hq hp0 hpq ha hb]
      exact neg_ne_zero.mpr hab
    · rw [overlapDet_reflect_of_bad v q p (σ a) (σ c)
        hq hp0 hpq ha hc]
      exact neg_ne_zero.mpr hac
    · rw [overlapDet_reflect_of_bad v q p (σ b) (σ c)
        hq hp0 hpq hb hc]
      exact neg_ne_zero.mpr hbc
  · intro p hp
    unfold reflect
    unfold s LRC14Grand.coloredBadTriangleMultipliers at hp
    have hpIoo := (Finset.mem_filter.mp hp).1
    rw [Finset.mem_Ioo] at hpIoo
    omega
  · intro p hp hfixed
    unfold reflect at hfixed
    unfold s LRC14Grand.coloredBadTriangleMultipliers at hp
    obtain ⟨hpWindow, ha, hb, _hc, _htop, hab, _hac, _hbc⟩ :=
      Finset.mem_filter.mp hp
    have hp0 : 0 < p := (Finset.mem_Ioo.mp hpWindow).1
    have hpq : p < q := (Finset.mem_Ioo.mp hpWindow).2
    have hneg := overlapDet_reflect_of_bad v q p (σ a) (σ b)
      hq hp0 hpq ha hb
    rw [hfixed] at hneg
    have : overlapDet v q p (σ a) (σ b) = 0 := by omega
    exact hab this

/-! ## Axiom audit -/

#print axioms inBand_reflect_iff
#print axioms failWitness_reflect_of_bad
#print axioms overlapDet_reflect_of_bad
#print axioms abs_overlapDet_reflect_of_bad
#print axioms overlapTripleCoeff_reflect_of_bad
#print axioms overlapDet_half_eq_zero_of_bad
#print axioms overlapColorFiber_card_neg
#print axioms overlapColorFiber_union_neg_card_eq_zero_or_two_of_coprime
#print axioms coloredBadTripleMultipliers_card_even
#print axioms coloredBadTriangleMultipliers_card_even

end LRC14Concrete
end LonelyRunner
