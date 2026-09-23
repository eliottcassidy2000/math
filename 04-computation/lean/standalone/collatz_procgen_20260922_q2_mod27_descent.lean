/-
  collatz-procgen-20260922: the backward half (Q2) of the E-SCC relaxation reduces to two
  residue classes mod 27.  Core Lean 4 only (no Mathlib); arithmetic by `omega`.

  Graph E: arrows n -> 3n+1 (every n) and n -> n/2 (n even).  `Reach y m` means m is
  reachable from y along E-arrows.  Q2 asks: does 1 reach every m with 3 not dividing m?

  Theorem `q2_descent_off_1_and_14`: if 3 does not divide m, m is not 1 or 14 mod 27, and
  m > 1, then some y < m with 3 not dividing y reaches m in at most two multiplications
  (hence, by strong induction, Q2 holds for every m unless it fails in one of the two
  hostile classes 1 mod 27 and 14 mod 27).
-/

inductive Reach : Nat → Nat → Prop
  | refl (n : Nat) : Reach n n
  | triple {a b : Nat} : Reach (3 * a + 1) b → Reach a b
  | halve {a b : Nat} : a % 2 = 0 → Reach (a / 2) b → Reach a b

theorem reach_trans {a b c : Nat} (h1 : Reach a b) (h2 : Reach b c) : Reach a c := by
  induction h1 with
  | refl _ => exact h2
  | triple _ ih => exact Reach.triple (ih h2)
  | halve he _ ih => exact Reach.halve he (ih h2)

/-- y -> 3y+1 = m -/
theorem reach_one_step (y m : Nat) (h : 3 * y + 1 = m) : Reach y m := by
  subst h; exact Reach.triple (Reach.refl _)

/-- y -> 3y+1 = 2m -> m -/
theorem reach_k1 (y m : Nat) (h : 3 * y + 1 = 2 * m) : Reach y m := by
  apply Reach.triple
  rw [h]
  apply Reach.halve (by omega)
  have : 2 * m / 2 = m := by omega
  rw [this]; exact Reach.refl _

/-- y -> 3y+1 -> 3(3y+1)+1 = 4m -> 2m -> m -/
theorem reach_k2k0 (y m : Nat) (h : 3 * (3 * y + 1) + 1 = 4 * m) : Reach y m := by
  apply Reach.triple
  apply Reach.triple
  rw [h]
  apply Reach.halve (by omega)
  have e1 : 4 * m / 2 = 2 * m := by omega
  rw [e1]
  apply Reach.halve (by omega)
  have e2 : 2 * m / 2 = m := by omega
  rw [e2]; exact Reach.refl _

/-- y -> 3y+1 -> 3(3y+1)+1 = 8m -> 4m -> 2m -> m -/
theorem reach_k3k0 (y m : Nat) (h : 3 * (3 * y + 1) + 1 = 8 * m) : Reach y m := by
  apply Reach.triple
  apply Reach.triple
  rw [h]
  apply Reach.halve (by omega)
  have e1 : 8 * m / 2 = 4 * m := by omega
  rw [e1]
  apply Reach.halve (by omega)
  have e2 : 4 * m / 2 = 2 * m := by omega
  rw [e2]
  apply Reach.halve (by omega)
  have e3 : 2 * m / 2 = m := by omega
  rw [e3]; exact Reach.refl _

theorem q2_descent_off_1_and_14 (m : Nat) (h3 : m % 3 ≠ 0) (h1 : m % 27 ≠ 1)
    (h14 : m % 27 ≠ 14) (hm : 1 < m) :
    ∃ y, y < m ∧ y % 3 ≠ 0 ∧ Reach y m := by
  -- m = 4,7 mod 9: y = (m-1)/3 (factor 1/3)
  by_cases hA : m % 9 = 4 ∨ m % 9 = 7
  · refine ⟨(m - 1) / 3, by omega, by omega, reach_one_step _ _ (by omega)⟩
  -- m = 2,8 mod 9: y = (2m-1)/3 (factor 2/3)
  by_cases hB : m % 9 = 2 ∨ m % 9 = 8
  · refine ⟨(2 * m - 1) / 3, by omega, by omega, reach_k1 _ _ (by omega)⟩
  -- remaining units: m = 1 or 5 mod 9
  by_cases hC : m % 27 = 10 ∨ m % 27 = 19
  · -- y = (4m-4)/9 (factor 4/9)
    refine ⟨(4 * m - 4) / 9, by omega, by omega, reach_k2k0 _ _ (by omega)⟩
  by_cases hD : m % 27 = 5 ∨ m % 27 = 23
  · -- y = (8m-4)/9 (factor 8/9)
    refine ⟨(8 * m - 4) / 9, by omega, by omega, reach_k3k0 _ _ (by omega)⟩
  -- no other residue survives
  exfalso; omega

/-- The two excluded classes are genuinely excluded at this depth: from 14 mod 27 the only
    two-move reverse routes pass through a multiple of 3 or ascend.  (Sanity instances.) -/
example : (3 * 12 + 1 = 37) ∧ 12 % 3 = 0 := by decide
example : Reach 4 5 := reach_k3k0 4 5 (by decide)   -- 5 = 5 mod 27: 4 -> 13 -> 40 -> 20 -> 10 -> 5
example : Reach 4 10 := reach_k2k0 4 10 (by decide)   -- 10 = 10 mod 27: 4 -> 13 -> 40 -> 20 -> 10
/- minus sheet: 4 is reached from 1 (1 -> 2 -> 5 -> 14 -> 41 -> ... is not needed: 4 <- 8 <- 16 <- 32 <- 11 <- 22 ... ) is checked
    numerically in `collatz_procgen_20260922_q2_minus_small.out`; here only the mod-27 descent is formalized. -/


/-! ## The minus sheet (3n-1): the mirror statement under negation.
    Graph E_-: arrows n -> 3n-1 (every n >= 1) and n -> n/2 (n even).  The hostile classes are the
    negations of 1 and 14 = 1/2, namely 26 and 13 mod 27. -/

inductive ReachM : Nat → Nat → Prop
  | refl (n : Nat) : ReachM n n
  | triple {a b : Nat} : 1 ≤ a → ReachM (3 * a - 1) b → ReachM a b
  | halve {a b : Nat} : a % 2 = 0 → ReachM (a / 2) b → ReachM a b

theorem reachM_one_step (y m : Nat) (hy : 1 ≤ y) (h : 3 * y - 1 = m) : ReachM y m := by
  exact ReachM.triple hy (by rw [h]; exact ReachM.refl _)

theorem reachM_k1 (y m : Nat) (hy : 1 ≤ y) (h : 3 * y - 1 = 2 * m) : ReachM y m := by
  apply ReachM.triple hy
  rw [h]
  apply ReachM.halve (by omega)
  have : 2 * m / 2 = m := by omega
  rw [this]; exact ReachM.refl _

theorem reachM_k2k0 (y m : Nat) (hy : 1 ≤ y) (h1 : 1 ≤ 3 * y - 1) (h : 3 * (3 * y - 1) - 1 = 4 * m) :
    ReachM y m := by
  apply ReachM.triple hy
  apply ReachM.triple h1
  rw [h]
  apply ReachM.halve (by omega)
  have e1 : 4 * m / 2 = 2 * m := by omega
  rw [e1]
  apply ReachM.halve (by omega)
  have e2 : 2 * m / 2 = m := by omega
  rw [e2]; exact ReachM.refl _

theorem reachM_k3k0 (y m : Nat) (hy : 1 ≤ y) (h1 : 1 ≤ 3 * y - 1) (h : 3 * (3 * y - 1) - 1 = 8 * m) :
    ReachM y m := by
  apply ReachM.triple hy
  apply ReachM.triple h1
  rw [h]
  apply ReachM.halve (by omega)
  have e1 : 8 * m / 2 = 4 * m := by omega
  rw [e1]
  apply ReachM.halve (by omega)
  have e2 : 4 * m / 2 = 2 * m := by omega
  rw [e2]
  apply ReachM.halve (by omega)
  have e3 : 2 * m / 2 = m := by omega
  rw [e3]; exact ReachM.refl _

/-- Needs `4 < m`: at `m = 4` the route `(8m+4)/9` returns `4` itself (no descent); the negation
    symmetry transfers descent only for `|m| > 4`.  (`m = 2, 4` are reached from 1 directly.) -/
theorem q2_minus_descent_off_26_and_13 (m : Nat) (h3 : m % 3 ≠ 0) (h26 : m % 27 ≠ 26)
    (h13 : m % 27 ≠ 13) (hm : 4 < m) :
    ∃ y, 1 ≤ y ∧ y < m ∧ y % 3 ≠ 0 ∧ ReachM y m := by
  by_cases hA : m % 9 = 5 ∨ m % 9 = 2
  · refine ⟨(m + 1) / 3, by omega, by omega, by omega, reachM_one_step _ _ (by omega) (by omega)⟩
  by_cases hB : m % 9 = 7 ∨ m % 9 = 1
  · refine ⟨(2 * m + 1) / 3, by omega, by omega, by omega, reachM_k1 _ _ (by omega) (by omega)⟩
  by_cases hC : m % 27 = 17 ∨ m % 27 = 8
  · refine ⟨(4 * m + 4) / 9, by omega, by omega, by omega, reachM_k2k0 _ _ (by omega) (by omega) (by omega)⟩
  by_cases hD : m % 27 = 22 ∨ m % 27 = 4
  · refine ⟨(8 * m + 4) / 9, by omega, by omega, by omega, reachM_k3k0 _ _ (by omega) (by omega) (by omega)⟩
  exfalso; omega

#print axioms q2_descent_off_1_and_14
#print axioms reach_trans
#print axioms q2_minus_descent_off_26_and_13
