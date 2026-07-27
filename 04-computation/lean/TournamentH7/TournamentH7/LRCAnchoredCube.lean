import Mathlib.NumberTheory.FLT.Three

/-!
# The anchored cubic rectangle cannot vanish

This is the arithmetic kernel of THM-2517.  Fermat's theorem for exponent
three says that a positive rational anchor separates every nonnegative
rational point with a different owner coordinate after the cubic embedding.
-/

namespace LRCAnchoredCube

/-- A nonnegative rational cubic rectangle based at a positive anchor can
vanish only at the anchor itself. -/
theorem cube_rectangle_eq_zero_iff {a x y : ℚ} (ha : 0 < a)
    (hx : 0 ≤ x) (hy : 0 ≤ y) :
    x ^ 3 - y ^ 3 - a ^ 3 = 0 ↔ x = a ∧ y = 0 := by
  constructor
  · intro h
    have hEq : y ^ 3 + a ^ 3 = x ^ 3 := by linarith
    by_cases hy0 : y = 0
    · subst y
      have hpow : a ^ 3 = x ^ 3 := by simpa using hEq
      have hax : a = x := (show Odd 3 by decide).pow_injective hpow
      exact ⟨hax.symm, rfl⟩
    · have ha0 : a ≠ 0 := ne_of_gt ha
      have hx0 : x ≠ 0 := by
        intro hxzero
        subst x
        have hpos : 0 < y ^ 3 + a ^ 3 :=
          add_pos_of_nonneg_of_pos (pow_nonneg hy 3) (pow_pos ha 3)
        simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow] at hEq
        exact (ne_of_gt hpos) hEq
      have fltQ : FermatLastTheoremWith ℚ 3 :=
        fermatLastTheoremFor_iff_rat.mp fermatLastTheoremThree
      exact (fltQ y a x hy0 ha0 hx0 hEq).elim
  · rintro ⟨rfl, rfl⟩
    ring

/-- Local branch-free form: at any owner-nonflat target column, every source
leg has a nonzero cubic anchored rectangle. -/
theorem cube_rectangles_nonzero_at_nonflat_column
    (A : Fin 7 → Fin 13 → ℚ) (a : ℚ) (ha : 0 < a)
    (hA : ∀ ell s, 0 ≤ A ell s) (s : Fin 13) (hs : A 0 s ≠ a) :
    ∀ ell : Fin 7, ell ≠ 0 → A 0 s ^ 3 - A ell s ^ 3 - a ^ 3 ≠ 0 := by
  intro ell _hell hzero
  have hbase := (cube_rectangle_eq_zero_iff ha (hA 0 s) (hA ell s)).mp hzero
  exact hs hbase.1

/-- Existence form consumed by an anchored nonflat response table. -/
theorem exists_column_with_all_cube_rectangles
    (A : Fin 7 → Fin 13 → ℚ) (a : ℚ) (ha : 0 < a)
    (hA : ∀ ell s, 0 ≤ A ell s)
    (hnonflat : ∃ s : Fin 13, A 0 s ≠ a) :
    ∃ s : Fin 13, ∀ ell : Fin 7,
      ell ≠ 0 → A 0 s ^ 3 - A ell s ^ 3 - a ^ 3 ≠ 0 := by
  obtain ⟨s, hs⟩ := hnonflat
  exact ⟨s, cube_rectangles_nonzero_at_nonflat_column A a ha hA s hs⟩

/-- Global equality boundary for the simplified anchored cubic rectangles:
they all vanish exactly when the owner row is constantly the anchor and all
six nonowner rows vanish.  Under the delta-anchor hypothesis of THM-2517,
these expressions are the full four-cell anchored rectangles. -/
theorem all_cube_rectangles_eq_zero_iff_pure_baseline
    (A : Fin 7 → Fin 13 → ℚ) (a : ℚ) (ha : 0 < a)
    (hA : ∀ ell s, 0 ≤ A ell s) :
    (∀ s ell, ell ≠ 0 → A 0 s ^ 3 - A ell s ^ 3 - a ^ 3 = 0) ↔
      (∀ s, A 0 s = a) ∧ (∀ ell s, ell ≠ 0 → A ell s = 0) := by
  constructor
  · intro hzero
    constructor
    · intro s
      have hleg := hzero s (1 : Fin 7) (by decide)
      have hbase :=
        (cube_rectangle_eq_zero_iff ha (hA 0 s) (hA (1 : Fin 7) s)).mp hleg
      exact hbase.1
    · intro ell s hell
      have hbase :=
        (cube_rectangle_eq_zero_iff ha (hA 0 s) (hA ell s)).mp
          (hzero s ell hell)
      exact hbase.2
  · rintro ⟨howner, hnonowner⟩ s ell hell
    rw [howner s, hnonowner ell s hell]
    ring

end LRCAnchoredCube

#print axioms LRCAnchoredCube.cube_rectangle_eq_zero_iff
#print axioms LRCAnchoredCube.cube_rectangles_nonzero_at_nonflat_column
#print axioms LRCAnchoredCube.exists_column_with_all_cube_rectangles
#print axioms LRCAnchoredCube.all_cube_rectangles_eq_zero_iff_pure_baseline
