/-
  TournamentH7.LRCOverlapColorFibers

  Activity-weighted arithmetic for the colored overlap movie.  For one
  ordered pair of positive runners, equality of two `overlapDet` colors
  forces their first witnesses into one residue class modulo the reduced
  speed `a / gcd(a,b)`.  On the short window `q <= 7a`, multiplier-to-witness
  uniqueness then bounds every color fiber by `gcd(a,b)`.

  Tournament-analysis audit: the vertices are multiplier events rather than
  runners.  The pair observable is `overlapDet`, equality of colors is the
  switch, and the fiber cardinality is the preserved statistic.  Passing to
  the zero/nonzero support graph destroys both the gcd spacing and activity
  multiplicity; there is therefore no faithful static-runner tournament at
  this step.
-/

import TournamentH7.LRCDeepCount
import TournamentH7.LRCSevenOverlapRelations

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- Simultaneously bad multipliers on one ordered pair carrying one fixed
cross-product color. -/
def overlapColorFiber (v : Fin 13 → ℤ) (q : ℕ) (i j : Fin 13)
    (d : ℤ) : Finset ℕ :=
  (Finset.Ioo 0 q).filter fun p =>
    ¬ inBand v q p i ∧ ¬ inBand v q p j ∧ overlapDet v q p i j = d

/-- Every overlap color is a multiple of the gcd of the two speeds. -/
theorem gcd_dvd_overlapDet (v : Fin 13 → ℤ) (q p : ℕ) (i j : Fin 13) :
    (Int.gcd (v i) (v j) : ℤ) ∣ overlapDet v q p i j := by
  unfold overlapDet
  exact dvd_sub
    (dvd_mul_of_dvd_left (Int.gcd_dvd_left (v i) (v j)) _)
    (dvd_mul_of_dvd_left (Int.gcd_dvd_right (v i) (v j)) _)

/-- Consequently a nonzero color has magnitude at least the speed gcd. -/
theorem gcd_le_abs_overlapDet_of_ne (v : Fin 13 → ℤ) (q p : ℕ)
    (i j : Fin 13) (hne : overlapDet v q p i j ≠ 0) :
    (Int.gcd (v i) (v j) : ℤ) ≤ |overlapDet v q p i j| :=
  Int.le_abs_of_dvd hne (gcd_dvd_overlapDet v q p i j)

/-- The exact residue-class step behind the fiber bound.  Equal pair colors
force the first witnesses to agree modulo `a / gcd(a,b)`.  No badness or
window assumption is needed for this algebraic statement. -/
theorem failWitness_modEq_of_overlapDet_eq
    (v : Fin 13 → ℤ) (q p₁ p₂ : ℕ) (i j : Fin 13)
    (hvi : 0 < v i)
    (hcolor : overlapDet v q p₁ i j = overlapDet v q p₂ i j) :
    failWitness v q p₁ i ≡ failWitness v q p₂ i
      [ZMOD v i / (Int.gcd (v i) (v j) : ℤ)] := by
  have hmul :
      v j * failWitness v q p₁ i ≡
        v j * failWitness v q p₂ i [ZMOD v i] := by
    rw [Int.modEq_iff_dvd]
    have heq :
        v j * failWitness v q p₂ i - v j * failWitness v q p₁ i =
          v i * (failWitness v q p₂ j - failWitness v q p₁ j) := by
      unfold overlapDet at hcolor
      linear_combination hcolor
    rw [heq]
    exact dvd_mul_right _ _
  exact Int.ModEq.cancel_left_div_gcd hvi hmul

/-- Convert an integer congruence to the congruence of `toNat` values when all
three integers are nonnegative. -/
private theorem toNat_modEq_of_nonneg {m x y : ℤ}
    (hm : 0 ≤ m) (hx : 0 ≤ x) (hy : 0 ≤ y)
    (h : x ≡ y [ZMOD m]) :
    x.toNat ≡ y.toNat [MOD m.toNat] := by
  change x % m = y % m at h
  change x.toNat % m.toNat = y.toNat % m.toNat
  rw [← Int.natCast_inj, Int.natCast_mod, Int.natCast_mod]
  simpa [Int.toNat_of_nonneg hm, Int.toNat_of_nonneg hx,
    Int.toNat_of_nonneg hy] using h

/-- A finite subset of `[1,A*g]` contained in one residue class modulo the
positive spacing `A` has at most `g` elements. -/
private theorem card_le_of_pairwise_modEq
    (s : Finset ℕ) (A g : ℕ) (hA : 0 < A)
    (hmem : ∀ x ∈ s, 1 ≤ x ∧ x ≤ A * g)
    (hmod : ∀ x ∈ s, ∀ y ∈ s, x ≡ y [MOD A]) :
    s.card ≤ g := by
  classical
  by_cases hs : s = ∅
  · simp [hs]
  · obtain ⟨c, hc⟩ := Finset.nonempty_iff_ne_empty.mpr hs
    by_cases hc0 : c % A = 0
    · have hmap : ∀ x ∈ s, x / A ∈ Finset.Icc 1 g := by
        intro x hx
        have hxc := hmod x hx c hc
        change x % A = c % A at hxc
        have hx0 : x % A = 0 := hxc.trans hc0
        have hbounds := hmem x hx
        rw [Finset.mem_Icc]
        constructor
        · exact Nat.div_pos
            (Nat.le_of_dvd hbounds.1 (Nat.dvd_of_mod_eq_zero hx0)) hA
        · exact Nat.div_le_of_le_mul (by
            simpa [mul_comm] using hbounds.2)
      have hinj : ∀ x ∈ s, ∀ y ∈ s, x / A = y / A → x = y := by
        intro x hx y hy hdiv
        have hxy := hmod x hx y hy
        change x % A = y % A at hxy
        calc
          x = x % A + A * (x / A) := (Nat.mod_add_div x A).symm
          _ = y % A + A * (y / A) := by rw [hxy, hdiv]
          _ = y := Nat.mod_add_div y A
      calc
        s.card ≤ (Finset.Icc 1 g).card :=
          Finset.card_le_card_of_injOn (fun x => x / A) hmap hinj
        _ = g := by rw [Nat.card_Icc]; omega
    · have hmap : ∀ x ∈ s, x / A ∈ Finset.range g := by
        intro x hx
        have hxc := hmod x hx c hc
        change x % A = c % A at hxc
        have hx0 : x % A ≠ 0 := by
          intro hz
          exact hc0 (hxc.symm.trans hz)
        rw [Finset.mem_range]
        have hbounds := hmem x hx
        have hxne : x ≠ A * g := by
          intro hEq
          apply hx0
          rw [hEq]
          exact Nat.mod_eq_zero_of_dvd (dvd_mul_right A g)
        have hxlt : x < A * g := lt_of_le_of_ne hbounds.2 hxne
        exact (Nat.div_lt_iff_lt_mul hA).2 (by
          simpa [mul_comm] using hxlt)
      have hinj : ∀ x ∈ s, ∀ y ∈ s, x / A = y / A → x = y := by
        intro x hx y hy hdiv
        have hxy := hmod x hx y hy
        change x % A = y % A at hxy
        calc
          x = x % A + A * (x / A) := (Nat.mod_add_div x A).symm
          _ = y % A + A * (y / A) := by rw [hxy, hdiv]
          _ = y := Nat.mod_add_div y A
      simpa using
        (Finset.card_le_card_of_injOn (fun x => x / A) hmap hinj)

/-- **Fixed-color activity bound.**  For positive pair speeds `a = v i` and
`b = v j`, every fixed overlap-color fiber on the short window `q ≤ 7a`
contains at most `gcd(a,b)` simultaneously bad multipliers. -/
theorem overlapColorFiber_card_le_gcd
    (v : Fin 13 → ℤ) (q : ℕ) (i j : Fin 13) (d : ℤ)
    (hq : 0 < q) (hvi : 0 < v i) (hvj : 0 < v j)
    (hwin : (q : ℤ) ≤ 7 * v i) :
    (overlapColorFiber v q i j d).card ≤ Int.gcd (v i) (v j) := by
  classical
  let a : ℕ := (v i).toNat
  let b : ℕ := (v j).toNat
  let g : ℕ := Nat.gcd a b
  let A : ℕ := a / g
  have haZ : (a : ℤ) = v i := by
    simp [a, Int.toNat_of_nonneg hvi.le]
  have hbZ : (b : ℤ) = v j := by
    simp [b, Int.toNat_of_nonneg hvj.le]
  have hapos : 0 < a := by
    exact_mod_cast (show (0 : ℤ) < (a : ℤ) by rw [haZ]; exact hvi)
  have hnatAbsA : (v i).natAbs = a := by
    exact Int.ofNat.inj ((Int.natAbs_of_nonneg hvi.le).trans haZ.symm)
  have hnatAbsB : (v j).natAbs = b := by
    exact Int.ofNat.inj ((Int.natAbs_of_nonneg hvj.le).trans hbZ.symm)
  have hgcd : g = Int.gcd (v i) (v j) := by
    dsimp only [g]
    rw [Int.gcd_def, hnatAbsA, hnatAbsB]
  have hgcdZ : (g : ℤ) = (Int.gcd (v i) (v j) : ℤ) :=
    congrArg (fun n : ℕ => (n : ℤ)) hgcd
  have hgpos : 0 < g := by
    dsimp only [g]
    rw [Nat.gcd_pos_iff]
    exact Or.inl hapos
  have hgda : g ∣ a := Nat.gcd_dvd_left a b
  have hfactor : A * g = a := by
    dsimp only [A]
    exact Nat.div_mul_cancel hgda
  have hApos : 0 < A := by
    dsimp only [A]
    exact Nat.div_pos (Nat.le_of_dvd hapos hgda) hgpos
  have hmodulus : (A : ℤ) = v i / (Int.gcd (v i) (v j) : ℤ) := by
    calc
      (A : ℤ) = ((a / g : ℕ) : ℤ) := rfl
      _ = (a : ℤ) / (g : ℤ) := Int.natCast_div a g
      _ = v i / (Int.gcd (v i) (v j) : ℤ) := by rw [haZ, hgcdZ]
  let witnessSet : Finset ℕ :=
    (overlapColorFiber v q i j d).image fun p => (failWitness v q p i).toNat
  have hwmem : ∀ n ∈ witnessSet, 1 ≤ n ∧ n ≤ A * g := by
    intro n hn
    simp only [witnessSet, Finset.mem_image] at hn
    obtain ⟨p, hp, rfl⟩ := hn
    simp only [overlapColorFiber, Finset.mem_filter, Finset.mem_Ioo] at hp
    have hbad := hp.2.1
    have hwpos := failWitness_pos_of_window v q p i hq hp.1.1 hvi
      (by linarith) hbad
    have hwle := failWitness_le v q p i hq hp.1.2 hvi hbad
    rw [hfactor]
    have hw0 : 0 ≤ failWitness v q p i := by omega
    have hwcast : ((failWitness v q p i).toNat : ℤ) =
        failWitness v q p i := Int.toNat_of_nonneg hw0
    constructor
    · exact_mod_cast (show (1 : ℤ) ≤ ((failWitness v q p i).toNat : ℤ) by
        rw [hwcast]
        exact hwpos)
    · exact_mod_cast (show ((failWitness v q p i).toNat : ℤ) ≤ (a : ℤ) by
        rw [hwcast, haZ]
        exact hwle)
  have hwmod : ∀ x ∈ witnessSet, ∀ y ∈ witnessSet, x ≡ y [MOD A] := by
    intro x hx y hy
    simp only [witnessSet, Finset.mem_image] at hx hy
    obtain ⟨p₁, hp₁, rfl⟩ := hx
    obtain ⟨p₂, hp₂, rfl⟩ := hy
    simp only [overlapColorFiber, Finset.mem_filter, Finset.mem_Ioo] at hp₁ hp₂
    have hres := failWitness_modEq_of_overlapDet_eq v q p₁ p₂ i j hvi
      (hp₁.2.2.2.trans hp₂.2.2.2.symm)
    have hbad₁ := hp₁.2.1
    have hbad₂ := hp₂.2.1
    have hn₁ := failWitness_pos_of_window v q p₁ i hq hp₁.1.1 hvi
      (by linarith) hbad₁
    have hn₂ := failWitness_pos_of_window v q p₂ i hq hp₂.1.1 hvi
      (by linarith) hbad₂
    have hn₁0 : 0 ≤ failWitness v q p₁ i := by omega
    have hn₂0 : 0 ≤ failWitness v q p₂ i := by omega
    rw [← hmodulus] at hres
    simpa using toNat_modEq_of_nonneg (Int.natCast_nonneg A)
      hn₁0 hn₂0 hres
  have hfiberInj : Set.InjOn (fun p => (failWitness v q p i).toNat)
      (overlapColorFiber v q i j d) := by
    intro p₁ hp₁ p₂ hp₂ heq
    change p₁ ∈ overlapColorFiber v q i j d at hp₁
    change p₂ ∈ overlapColorFiber v q i j d at hp₂
    simp only [overlapColorFiber, Finset.mem_filter, Finset.mem_Ioo] at hp₁ hp₂
    have hbad₁ := hp₁.2.1
    have hbad₂ := hp₂.2.1
    have hn₁ := failWitness_pos_of_window v q p₁ i hq hp₁.1.1 hvi
      (by linarith) hbad₁
    have hn₂ := failWitness_pos_of_window v q p₂ i hq hp₂.1.1 hvi
      (by linarith) hbad₂
    have hn₁0 : 0 ≤ failWitness v q p₁ i := by omega
    have hn₂0 : 0 ≤ failWitness v q p₂ i := by omega
    have hcast₁ : ((failWitness v q p₁ i).toNat : ℤ) =
        failWitness v q p₁ i := Int.toNat_of_nonneg hn₁0
    have hcast₂ : ((failWitness v q p₂ i).toNat : ℤ) =
        failWitness v q p₂ i := Int.toNat_of_nonneg hn₂0
    have heqZ : failWitness v q p₁ i = failWitness v q p₂ i := by
      calc
        failWitness v q p₁ i = ((failWitness v q p₁ i).toNat : ℤ) := hcast₁.symm
        _ = ((failWitness v q p₂ i).toNat : ℤ) :=
          congrArg (fun n : ℕ => (n : ℤ)) heq
        _ = failWitness v q p₂ i := hcast₂
    exact window_unique_p v q i p₁ p₂ hq hvi hwin hbad₁ hbad₂ heqZ
  calc
    (overlapColorFiber v q i j d).card = witnessSet.card := by
      change (overlapColorFiber v q i j d).card =
        ((overlapColorFiber v q i j d).image fun p =>
          (failWitness v q p i).toNat).card
      exact (Finset.card_image_iff.mpr hfiberInj).symm
    _ ≤ g := card_le_of_pairwise_modEq witnessSet A g hApos hwmem hwmod
    _ = Int.gcd (v i) (v j) := hgcd

/-- Coprime positive speeds support at most one simultaneously bad event of
each fixed overlap color on the short window. -/
theorem overlapColorFiber_card_le_one_of_gcd_eq_one
    (v : Fin 13 → ℤ) (q : ℕ) (i j : Fin 13) (d : ℤ)
    (hq : 0 < q) (hvi : 0 < v i) (hvj : 0 < v j)
    (hwin : (q : ℤ) ≤ 7 * v i)
    (hcop : Int.gcd (v i) (v j) = 1) :
    (overlapColorFiber v q i j d).card ≤ 1 := by
  simpa [hcop] using overlapColorFiber_card_le_gcd
    v q i j d hq hvi hvj hwin

/-! ## Axiom audit -/

#print axioms gcd_dvd_overlapDet
#print axioms gcd_le_abs_overlapDet_of_ne
#print axioms failWitness_modEq_of_overlapDet_eq
#print axioms overlapColorFiber_card_le_gcd
#print axioms overlapColorFiber_card_le_one_of_gcd_eq_one

end LRC14Concrete
end LonelyRunner
