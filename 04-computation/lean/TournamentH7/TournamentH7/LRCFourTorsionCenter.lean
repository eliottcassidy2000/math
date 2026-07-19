import Mathlib

/-!
# Four-torsion centre gauges (THM-1206)

This file formalizes only the finite modular core of THM-1206.  For the
labelled centre residue vector `r = (1,2,3)` in `(ZMod 4)^3`, one common
gauge `a` satisfying

`-a * eᵢ = rᵢ`

exists exactly when `e = r` or `e = -r`.  The two possibilities have the
canonical gauges `a = 3` and `a = 1`, respectively.

The analytic bridge is deliberately external: this module does not formalize
fractional parts on the circle, the Bezout step which turns a real hitting
time into an integral common gauge, or the converse lift from a modular gauge
to a circle hitting time.  Thus no geometric or measure-theoretic conclusion
is hidden in the finite calculation below.
-/

namespace LRC14.FourTorsionCenter

/-- A labelled residue vector with three coordinates modulo four. -/
abbrev ResidueVector := Fin 3 → ZMod 4

/-- The labelled four-torsion centre numerator `(1,2,3)`. -/
def centre123 : ResidueVector := ![1, 2, 3]

/-- `a` is one common gauge taking every coordinate of `e` to the labelled
centre.  The quantifier over one `a` is the essential global datum. -/
def IsCommonGauge (a : ZMod 4) (e : ResidueVector) : Prop :=
  ∀ i, -a * e i = centre123 i

private lemma four_eq_zero : (4 : ZMod 4) = 0 :=
  ZMod.natCast_self 4

private lemma two_ne_zero : (2 : ZMod 4) ≠ 0 := by
  intro h
  have hval := congrArg ZMod.val h
  have htwo : ((2 : ZMod 4).val) = 2 :=
    ZMod.val_natCast_of_lt (by norm_num)
  rw [htwo] at hval
  norm_num at hval

private lemma neg_three_eq_one : -(3 : ZMod 4) = 1 := by
  calc
    -(3 : ZMod 4) = 1 - 4 := by ring
    _ = 1 := by rw [four_eq_zero]; ring

private lemma neg_one_eq_three : -(1 : ZMod 4) = 3 := by
  calc
    -(1 : ZMod 4) = 3 - 4 := by ring
    _ = 3 := by rw [four_eq_zero]; ring

private lemma seven_eq_three : (7 : ZMod 4) = 3 := by
  calc
    (7 : ZMod 4) = 3 + 4 := by ring
    _ = 3 := by rw [four_eq_zero]; ring

/-- Gauge `3 = -1` fixes the centre vector under `e ↦ -a e`. -/
theorem gauge_three_for_centre : IsCommonGauge 3 centre123 := by
  intro i
  rw [neg_three_eq_one, one_mul]

/-- Gauge `1` takes the negative centre vector back to the centre. -/
theorem gauge_one_for_neg_centre : IsCommonGauge 1 (-centre123) := by
  intro i
  change -(1 : ZMod 4) * (-centre123 i) = centre123 i
  ring

private theorem scalar_gauge_cases (a x : ZMod 4)
    (h : -a * x = 1) : a = 1 ∨ a = 3 := by
  have hle : a.val ≤ 3 := by
    have hlt := a.val_lt
    omega
  interval_cases hval : a.val
  · have ha : a = (0 : ZMod 4) := by
      calc
        a = (a.val : ZMod 4) := (ZMod.natCast_zmod_val a).symm
        _ = 0 := by rw [hval]; rfl
    subst a
    simp at h
    exfalso
    apply two_ne_zero
    calc
      (2 : ZMod 4) = 2 * 1 := by ring
      _ = 2 * 0 := by rw [← h]
      _ = 0 := by ring
  · left
    calc
      a = (a.val : ZMod 4) := (ZMod.natCast_zmod_val a).symm
      _ = 1 := by rw [hval]; rfl
  · have ha : a = (2 : ZMod 4) := by
      calc
        a = (a.val : ZMod 4) := (ZMod.natCast_zmod_val a).symm
        _ = 2 := by rw [hval]; rfl
    subst a
    exfalso
    apply two_ne_zero
    calc
      (2 : ZMod 4) = 2 * 1 := by ring
      _ = 2 * (-(2 : ZMod 4) * x) := by rw [h]
      _ = -(4 : ZMod 4) * x := by ring
      _ = 0 := by rw [four_eq_zero]; simp
  · right
    calc
      a = (a.val : ZMod 4) := (ZMod.natCast_zmod_val a).symm
      _ = 3 := by rw [hval]; rfl

private theorem gauge_three_action (x : ZMod 4) :
    -(3 : ZMod 4) * x = x := by
  rw [neg_three_eq_one, one_mul]

/-- The first odd centre coordinate already forces a common gauge to be one
of the two units modulo four. -/
lemma commonGauge_unit_cases {a : ZMod 4} {e : ResidueVector}
    (h : IsCommonGauge a e) : a = 1 ∨ a = 3 := by
  apply scalar_gauge_cases a (e 0)
  simpa [IsCommonGauge, centre123] using h (0 : Fin 3)

/-- The exact finite core: a single common gauge exists iff the residue
vector is the centre vector up to the two-unit action `e ↦ ±e`. -/
theorem exists_commonGauge_iff (e : ResidueVector) :
    (∃ a : ZMod 4, IsCommonGauge a e) ↔
      e = centre123 ∨ e = -centre123 := by
  constructor
  · rintro ⟨a, ha⟩
    rcases commonGauge_unit_cases ha with rfl | rfl
    · right
      funext i
      have hi : -e i = centre123 i := by
        simpa [IsCommonGauge] using ha i
      have hneg := congrArg Neg.neg hi
      simpa using hneg
    · left
      funext i
      simpa only [gauge_three_action] using ha i
  · rintro (rfl | rfl)
    · exact ⟨3, gauge_three_for_centre⟩
    · exact ⟨1, gauge_one_for_neg_centre⟩

/-- On the centre vector itself, the canonical gauge `3` is unique. -/
theorem commonGauge_centre_iff (a : ZMod 4) :
    IsCommonGauge a centre123 ↔ a = 3 := by
  constructor
  · intro h
    have hnega : -a = (1 : ZMod 4) := by
      simpa [IsCommonGauge, centre123] using h (0 : Fin 3)
    have haNeg : a = -(1 : ZMod 4) := by
      have hneg := congrArg Neg.neg hnega
      simpa using hneg
    exact haNeg.trans neg_one_eq_three
  · rintro rfl
    exact gauge_three_for_centre

/-- On the negative centre vector, the canonical gauge `1` is unique. -/
theorem commonGauge_neg_centre_iff (a : ZMod 4) :
    IsCommonGauge a (-centre123) ↔ a = 1 := by
  constructor
  · intro h
    simpa [IsCommonGauge, centre123] using h (0 : Fin 3)
  · rintro rfl
    exact gauge_one_for_neg_centre

/-- The primitive residue vector of the proportional direction `(1,2,3)`. -/
def residue123 : ResidueVector := ![1, 2, 3]

/-- The primitive residue vector of the first nonproportional witness
`(1,2,7)`, reduced in `ZMod 4`. -/
def residue127 : ResidueVector := ![1, 2, 7]

/-- The literal `(1,2,3)` residue vector has canonical common gauge `3`. -/
theorem residue123_gauge_three : IsCommonGauge 3 residue123 := by
  simpa [residue123] using gauge_three_for_centre

/-- The nonproportional integer vector `(1,2,7)` has the same residue vector
modulo four and therefore the same canonical common gauge `3`. -/
theorem residue127_gauge_three : IsCommonGauge 3 residue127 := by
  intro i
  fin_cases i
  · change -(3 : ZMod 4) * 1 = 1
    rw [neg_three_eq_one, one_mul]
  · change -(3 : ZMod 4) * 2 = 2
    rw [neg_three_eq_one, one_mul]
  · change -(3 : ZMod 4) * 7 = 3
    rw [neg_three_eq_one, one_mul, seven_eq_three]

/-- Explicit finite witness for the centre hit of `(1,2,7)` modulo four. -/
theorem residue127_has_commonGauge :
    ∃ a : ZMod 4, IsCommonGauge a residue127 :=
  ⟨3, residue127_gauge_three⟩

#print axioms exists_commonGauge_iff
#print axioms commonGauge_centre_iff
#print axioms commonGauge_neg_centre_iff
#print axioms residue127_has_commonGauge

end LRC14.FourTorsionCenter
