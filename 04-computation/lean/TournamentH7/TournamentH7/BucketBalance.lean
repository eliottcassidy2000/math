/-
  TournamentH7.BucketBalance -- finite bucket half-line balances

  This module isolates the counting skeleton behind the merged tiling
  bucket constraints.  A finite population `alpha` is mapped to buckets
  by `q`; each move in a finite move set sends `x` to `step u x`.

  The oriented half-lines incident to a bucket split into internal
  half-lines and escaping half-lines, so their cardinalities add to
  `|fiber b| * |moves|`.  The unordered quotient balance
  `2 * internal_lines + escaping_lines = |fiber b| * |moves|` will be a
  later specialization once the fixed-point-free involution pairing for
  each move is formalized.
-/

import Mathlib.Tactic

namespace Tournament
namespace BucketBalance

universe u v w

variable {alpha : Type u} {beta : Type v} {move : Type w}

noncomputable section

/-! ### Fibers and oriented half-lines -/

/-- The finite fiber of bucket `b` under a bucket map `q`. -/
def fiber [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (b : beta) : Finset alpha := by
  classical
  exact Finset.univ.filter fun x => q x = b

/-- Oriented half-lines leaving points in bucket `b`, one for each move. -/
def incidentHalf [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (moves : Finset move) (b : beta) :
    Finset (alpha × move) := by
  classical
  exact (fiber q b).product moves

/-- Incident half-lines whose move remains inside bucket `b`. -/
def selfHalf [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) : Finset (alpha × move) := by
  classical
  exact (incidentHalf q moves b).filter fun xu => q (step xu.2 xu.1) = b

/-- Incident half-lines whose move leaves bucket `b`. -/
def crossHalf [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) : Finset (alpha × move) := by
  classical
  exact (incidentHalf q moves b).filter fun xu => q (step xu.2 xu.1) ≠ b

@[simp] theorem incidentHalf_card [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (moves : Finset move) (b : beta) :
    (incidentHalf q moves b).card = (fiber q b).card * moves.card := by
  classical
  simp [incidentHalf]

/-- The oriented half-line bucket balance.

This is the finite-set core of THM-346: every incident half-line from a
bucket either stays in the bucket or leaves it. -/
theorem halfLine_balance [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    (selfHalf q step moves b).card + (crossHalf q step moves b).card =
      (fiber q b).card * moves.card := by
  classical
  calc
    (selfHalf q step moves b).card + (crossHalf q step moves b).card
        = (incidentHalf q moves b).card := by
          simpa [selfHalf, crossHalf] using
            (Finset.card_filter_add_card_filter_not
              (s := incidentHalf q moves b)
              (p := fun xu : alpha × move => q (step xu.2 xu.1) = b))
    _ = (fiber q b).card * moves.card := incidentHalf_card q moves b

theorem selfHalf_card_le_total [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    (selfHalf q step moves b).card <= (fiber q b).card * moves.card := by
  classical
  have h := halfLine_balance q step moves b
  omega

theorem crossHalf_card_le_total [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    (crossHalf q step moves b).card <= (fiber q b).card * moves.card := by
  classical
  have h := halfLine_balance q step moves b
  omega

/-- No escaping half-lines iff every move from the bucket remains in it. -/
theorem crossHalf_card_eq_zero_iff [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    (crossHalf q step moves b).card = 0 <->
      forall x, x ∈ fiber q b -> forall u, u ∈ moves -> q (step u x) = b := by
  classical
  constructor
  · intro h x hx u hu
    by_contra hne
    have hmem : (x, u) ∈ crossHalf q step moves b := by
      simp [crossHalf, incidentHalf, hx, hu, hne]
    have hempty : crossHalf q step moves b = ∅ := Finset.card_eq_zero.mp h
    rw [hempty] at hmem
    simp at hmem
  · intro h
    apply Finset.card_eq_zero.mpr
    ext xu
    constructor
    · intro hmem
      rcases xu with ⟨x, u⟩
      have hfacts : (x ∈ fiber q b ∧ u ∈ moves) ∧ q (step u x) ≠ b := by
        simpa [crossHalf, incidentHalf] using hmem
      exact False.elim (hfacts.2 (h x hfacts.1.1 u hfacts.1.2))
    · intro hmem
      simp at hmem

/-- If every incident half-line stays inside the bucket, the balance
collapses to the full product count. -/
theorem selfHalf_card_eq_total_of_crossHalf_zero [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (h : (crossHalf q step moves b).card = 0) :
    (selfHalf q step moves b).card = (fiber q b).card * moves.card := by
  classical
  have hbal := halfLine_balance q step moves b
  omega

/-- If no incident half-line stays inside the bucket, all incident
half-lines escape. -/
theorem crossHalf_card_eq_total_of_selfHalf_zero [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (h : (selfHalf q step moves b).card = 0) :
    (crossHalf q step moves b).card = (fiber q b).card * moves.card := by
  classical
  have hbal := halfLine_balance q step moves b
  omega

/-! ### From oriented half-lines to unordered internal lines -/

/-- The natural partner of an oriented half-line `(x,u)` is `(step u x,u)`. -/
def pairHalf (step : move -> alpha -> alpha) (xu : alpha × move) : alpha × move :=
  (step xu.2 xu.1, xu.2)

@[simp] theorem pairHalf_fst (step : move -> alpha -> alpha) (xu : alpha × move) :
    (pairHalf step xu).1 = step xu.2 xu.1 := rfl

@[simp] theorem pairHalf_snd (step : move -> alpha -> alpha) (xu : alpha × move) :
    (pairHalf step xu).2 = xu.2 := rfl

/-- If every move is an involution, then `pairHalf` is an involution. -/
theorem pairHalf_pairHalf (step : move -> alpha -> alpha)
    (hstep : ∀ u, Function.Involutive (step u)) (xu : alpha × move) :
    pairHalf step (pairHalf step xu) = xu := by
  rcases xu with ⟨x, u⟩
  simp [pairHalf, hstep u x]

/-- Internal half-lines are closed under the partner map when the chosen move
is involutive. -/
theorem pairHalf_mem_selfHalf [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hstep : ∀ u, u ∈ moves -> Function.Involutive (step u))
    {xu : alpha × move} (hxu : xu ∈ selfHalf q step moves b) :
    pairHalf step xu ∈ selfHalf q step moves b := by
  classical
  rcases xu with ⟨x, u⟩
  have hfacts : (x ∈ fiber q b ∧ u ∈ moves) ∧ q (step u x) = b := by
    simpa [selfHalf, incidentHalf] using hxu
  have hxq : q x = b := by
    simpa [fiber] using hfacts.1.1
  have hback : q (step u (step u x)) = b := by
    rw [hstep u hfacts.1.2 x]
    exact hxq
  simp [selfHalf, incidentHalf, pairHalf, fiber, hfacts.1.2, hfacts.2, hback]

/-- If the selected move has no fixed points, an internal half-line and its
partner are distinct. -/
theorem pairHalf_ne_of_fixedPointFree [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hfixed : ∀ u, u ∈ moves -> ∀ x, step u x ≠ x)
    {xu : alpha × move} (hxu : xu ∈ selfHalf q step moves b) :
    pairHalf step xu ≠ xu := by
  classical
  rcases xu with ⟨x, u⟩
  have hfacts : (x ∈ fiber q b ∧ u ∈ moves) ∧ q (step u x) = b := by
    simpa [selfHalf, incidentHalf] using hxu
  intro hpair
  have hx : step u x = x := by
    simpa [pairHalf] using congrArg Prod.fst hpair
  exact hfixed u hfacts.1.2 x hx

/-- The internal unordered line count obtained by pairing internal half-lines. -/
def internalLineCount [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) : ℕ :=
  (selfHalf q step moves b).card / 2

theorem two_mul_internalLineCount_eq_selfHalf_card_of_even
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hself : Even (selfHalf q step moves b).card) :
    2 * internalLineCount q step moves b =
      (selfHalf q step moves b).card := by
  classical
  unfold internalLineCount
  rcases hself with ⟨k, hk⟩
  rw [hk]
  omega

/-- Unordered bucket balance once the internal oriented half-lines have been
paired into `internal` unordered lines. -/
theorem unordered_balance_of_selfHalf_card_eq_two_mul
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) (internal : ℕ)
    (hinternal : (selfHalf q step moves b).card = 2 * internal) :
    2 * internal + (crossHalf q step moves b).card =
      (fiber q b).card * moves.card := by
  classical
  have hbal := halfLine_balance q step moves b
  omega

/-- Unordered bucket balance in the common case where the internal oriented
half-lines have even cardinality. -/
theorem unordered_balance_of_even_selfHalf
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hself : Even (selfHalf q step moves b).card) :
    2 * internalLineCount q step moves b +
        (crossHalf q step moves b).card =
      (fiber q b).card * moves.card := by
  classical
  have htwo :=
    two_mul_internalLineCount_eq_selfHalf_card_of_even q step moves b hself
  have hbal := halfLine_balance q step moves b
  omega

end

end BucketBalance
end Tournament
