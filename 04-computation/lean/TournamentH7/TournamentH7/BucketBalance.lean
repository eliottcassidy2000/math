/-
  TournamentH7.BucketBalance -- finite bucket half-line balances

  This module isolates the counting skeleton behind the merged tiling
  bucket constraints.  A finite population `alpha` is mapped to buckets
  by `q`; each move in a finite move set sends `x` to `step u x`.

  The oriented half-lines incident to a bucket split into internal
  half-lines and escaping half-lines, so their cardinalities add to
  `|fiber b| * |moves|`.  The unordered quotient balance
  `2 * internal_lines + escaping_lines = |fiber b| * |moves|` follows once
  internal oriented half-lines are paired by fixed-point-free involutions.
-/

import Mathlib.Tactic

namespace Tournament
namespace BucketBalance

universe u v w z

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

/-! ### Transport rows between buckets -/

/-- Oriented half-lines starting in bucket `b` and landing in bucket `c`. -/
def transportHalf [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b c : beta) : Finset (alpha × move) := by
  classical
  exact (incidentHalf q moves b).filter fun xu => q (step xu.2 xu.1) = c

/-! ### Empty fibers and transport gaps -/

theorem fiber_eq_empty_iff [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (b : beta) :
    fiber q b = ∅ ↔ ∀ x, q x ≠ b := by
  classical
  constructor
  · intro h x hx
    have hxmem : x ∈ fiber q b := by
      simp [fiber, hx]
    rw [h] at hxmem
    simp at hxmem
  · intro h
    ext x
    constructor
    · intro hx
      have hq : q x = b := by
        simpa [fiber] using hx
      exact False.elim (h x hq)
    · intro hx
      simp at hx

theorem incidentHalf_eq_empty_of_fiber_eq_empty [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (moves : Finset move) (b : beta)
    (h : fiber q b = ∅) :
    incidentHalf q moves b = ∅ := by
  classical
  simp [incidentHalf, h]

theorem selfHalf_eq_empty_of_fiber_eq_empty [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) (h : fiber q b = ∅) :
    selfHalf q step moves b = ∅ := by
  classical
  simp [selfHalf, incidentHalf, h]

theorem crossHalf_eq_empty_of_fiber_eq_empty [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) (h : fiber q b = ∅) :
    crossHalf q step moves b = ∅ := by
  classical
  simp [crossHalf, incidentHalf, h]

theorem transportHalf_eq_empty_of_source_fiber_eq_empty
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b c : beta) (h : fiber q b = ∅) :
    transportHalf q step moves b c = ∅ := by
  classical
  simp [transportHalf, incidentHalf, h]

theorem transportHalf_eq_empty_of_target_fiber_eq_empty
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b c : beta) (h : fiber q c = ∅) :
    transportHalf q step moves b c = ∅ := by
  classical
  ext xu
  constructor
  · intro hxu
    have hfacts :
        xu ∈ incidentHalf q moves b ∧ q (step xu.2 xu.1) = c := by
      simpa [transportHalf] using hxu
    have htarget : step xu.2 xu.1 ∈ fiber q c := by
      simp [fiber, hfacts.2]
    rw [h] at htarget
    simp at htarget
  · intro hxu
    simp at hxu

theorem transportHalf_card_eq_zero_of_source_fiber_eq_empty
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b c : beta) (h : fiber q b = ∅) :
    (transportHalf q step moves b c).card = 0 := by
  rw [transportHalf_eq_empty_of_source_fiber_eq_empty q step moves b c h]
  simp

theorem transportHalf_card_eq_zero_of_target_fiber_eq_empty
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b c : beta) (h : fiber q c = ∅) :
    (transportHalf q step moves b c).card = 0 := by
  rw [transportHalf_eq_empty_of_target_fiber_eq_empty q step moves b c h]
  simp

theorem sum_transportHalf_card_eq_zero_of_source_fiber_eq_empty
    [Fintype alpha] [Fintype beta] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) (h : fiber q b = ∅) :
    (∑ c : beta, (transportHalf q step moves b c).card) = 0 := by
  classical
  have hzero : ∀ c : beta, (transportHalf q step moves b c).card = 0 := by
    intro c
    exact transportHalf_card_eq_zero_of_source_fiber_eq_empty q step moves b c h
  simp [hzero]

theorem sum_transportHalf_card_eq_zero_of_target_fiber_eq_empty
    [Fintype alpha] [Fintype beta] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (c : beta) (h : fiber q c = ∅) :
    (∑ b : beta, (transportHalf q step moves b c).card) = 0 := by
  classical
  have hzero : ∀ b : beta, (transportHalf q step moves b c).card = 0 := by
    intro b
    exact transportHalf_card_eq_zero_of_target_fiber_eq_empty q step moves b c h
  simp [hzero]

@[simp] theorem transportHalf_self [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    transportHalf q step moves b b = selfHalf q step moves b := by
  classical
  ext xu
  simp [transportHalf, selfHalf]

/-- The target-bucket transport row partitions all incident half-lines. -/
theorem sum_transportHalf_card_eq_incidentHalf_card
    [Fintype alpha] [Fintype beta] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    (∑ c : beta, (transportHalf q step moves b c).card) =
      (incidentHalf q moves b).card := by
  classical
  have h := Finset.card_eq_sum_card_fiberwise
    (s := incidentHalf q moves b) (t := (Finset.univ : Finset beta))
    (f := fun xu : alpha × move => q (step xu.2 xu.1))
    (by intro xu hxu; simp)
  simpa [transportHalf] using h.symm

/-- The off-diagonal target-bucket row is exactly the escaping half-line set. -/
theorem sum_transportHalf_card_offdiag_eq_crossHalf_card
    [Fintype alpha] [Fintype beta] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    (∑ c ∈ (Finset.univ.erase b),
        (transportHalf q step moves b c).card) =
      (crossHalf q step moves b).card := by
  classical
  have h := Finset.sum_card_fiberwise_eq_card_filter
    (s := incidentHalf q moves b) (t := (Finset.univ.erase b))
    (g := fun xu : alpha × move => q (step xu.2 xu.1))
  simpa [transportHalf, crossHalf] using h

/-! ### From oriented half-lines to unordered internal lines -/

/-- A finite set with a fixed-point-free involution has even cardinality.

This is the generic orbit-parity lemma behind the unordered bucket balance:
remove any element and its partner, then recurse on the remaining
partner-closed set. -/
theorem even_card_of_fixedPointFree_involutiveOn
    [DecidableEq alpha] (s : Finset alpha) (f : alpha -> alpha)
    (hmem : ∀ x, x ∈ s -> f x ∈ s)
    (hinv : ∀ x, x ∈ s -> f (f x) = x)
    (hfixed : ∀ x, x ∈ s -> f x ≠ x) :
    Even s.card := by
  classical
  let P : ℕ -> Prop := fun n =>
    ∀ s : Finset alpha, s.card = n ->
      (∀ x, x ∈ s -> f x ∈ s) ->
      (∀ x, x ∈ s -> f (f x) = x) ->
      (∀ x, x ∈ s -> f x ≠ x) ->
      Even s.card
  have hP : ∀ n, (∀ m, m < n -> P m) -> P n := by
    intro n ih s hs_card hmem hinv hfixed
    by_cases hs_empty : s = ∅
    · simp [hs_empty]
    · have hs_nonempty : s.Nonempty := Finset.nonempty_iff_ne_empty.mpr hs_empty
      rcases hs_nonempty with ⟨x, hx⟩
      let y := f x
      have hy : y ∈ s := hmem x hx
      have hyx : y ≠ x := hfixed x hx
      let t := (s.erase x).erase y
      have hy_erase : y ∈ s.erase x := by
        simp [hy, hyx]
      have ht_card_eq : t.card + 2 = s.card := by
        have h₁ : (s.erase x).card = s.card - 1 :=
          Finset.card_erase_of_mem hx
        have h₂ : t.card = (s.erase x).card - 1 :=
          Finset.card_erase_of_mem hy_erase
        have hpair : ({x, y} : Finset alpha).card = 2 := by
          simp [hyx.symm]
        have hsub : ({x, y} : Finset alpha) ⊆ s := by
          intro z hz
          simp at hz
          rcases hz with hz | hz
          · rwa [hz]
          · rwa [hz]
        have htwo : 2 ≤ s.card := by
          have := Finset.card_le_card hsub
          rwa [hpair] at this
        omega
      have ht_lt : t.card < n := by
        omega
      have ht_mem : ∀ z, z ∈ t -> f z ∈ t := by
        intro z hz
        have hz_ne_y : z ≠ y := (Finset.mem_erase.mp hz).1
        have hz_mem_erase_x : z ∈ s.erase x := (Finset.mem_erase.mp hz).2
        have hz_ne_x : z ≠ x := (Finset.mem_erase.mp hz_mem_erase_x).1
        have hz_s : z ∈ s := (Finset.mem_erase.mp hz_mem_erase_x).2
        have hfz_s : f z ∈ s := hmem z hz_s
        have hfz_ne_x : f z ≠ x := by
          intro hfz_x
          have hz_eq_y : z = y := by
            calc
              z = f (f z) := (hinv z hz_s).symm
              _ = f x := congrArg f hfz_x
              _ = y := rfl
          exact hz_ne_y hz_eq_y
        have hfz_ne_y : f z ≠ y := by
          intro hfz_y
          have hz_eq_x : z = x := by
            calc
              z = f (f z) := (hinv z hz_s).symm
              _ = f y := congrArg f hfz_y
              _ = x := by
                simpa [y] using hinv x hx
          exact hz_ne_x hz_eq_x
        exact Finset.mem_erase.mpr
          ⟨hfz_ne_y, Finset.mem_erase.mpr ⟨hfz_ne_x, hfz_s⟩⟩
      have ht_inv : ∀ z, z ∈ t -> f (f z) = z := by
        intro z hz
        have hz_s : z ∈ s := by
          exact (Finset.mem_erase.mp (Finset.mem_erase.mp hz).2).2
        exact hinv z hz_s
      have ht_fixed : ∀ z, z ∈ t -> f z ≠ z := by
        intro z hz
        have hz_s : z ∈ s := by
          exact (Finset.mem_erase.mp (Finset.mem_erase.mp hz).2).2
        exact hfixed z hz_s
      rcases ih t.card ht_lt t rfl ht_mem ht_inv ht_fixed with ⟨k, hk⟩
      refine ⟨k + 1, ?_⟩
      omega
  exact (Nat.strong_induction_on (p := P) s.card hP) s rfl hmem hinv hfixed

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

/-- Internal oriented half-lines have even cardinality when selected moves act
as fixed-point-free involutions. -/
theorem selfHalf_card_even_of_involutive_fixedPointFree
    [Fintype alpha] [DecidableEq alpha] [DecidableEq beta] [DecidableEq move]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hstep : ∀ u, u ∈ moves -> Function.Involutive (step u))
    (hfixed : ∀ u, u ∈ moves -> ∀ x, step u x ≠ x) :
    Even (selfHalf q step moves b).card := by
  classical
  refine even_card_of_fixedPointFree_involutiveOn
    (selfHalf q step moves b) (pairHalf step) ?_ ?_ ?_
  · intro xu hxu
    exact pairHalf_mem_selfHalf q step moves b hstep hxu
  · intro xu hxu
    rcases xu with ⟨x, u⟩
    have hfacts : (x ∈ fiber q b ∧ u ∈ moves) ∧ q (step u x) = b := by
      simpa [selfHalf, incidentHalf] using hxu
    simp [pairHalf, hstep u hfacts.1.2 x]
  · intro xu hxu
    exact pairHalf_ne_of_fixedPointFree q step moves b hfixed hxu

theorem two_mul_internalLineCount_eq_selfHalf_card_of_involutive_fixedPointFree
    [Fintype alpha] [DecidableEq alpha] [DecidableEq beta] [DecidableEq move]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hstep : ∀ u, u ∈ moves -> Function.Involutive (step u))
    (hfixed : ∀ u, u ∈ moves -> ∀ x, step u x ≠ x) :
    2 * internalLineCount q step moves b =
      (selfHalf q step moves b).card := by
  exact two_mul_internalLineCount_eq_selfHalf_card_of_even q step moves b
    (selfHalf_card_even_of_involutive_fixedPointFree q step moves b hstep hfixed)

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

/-- Unordered bucket balance for fixed-point-free involutive move systems. -/
theorem unordered_balance_of_involutive_fixedPointFree
    [Fintype alpha] [DecidableEq alpha] [DecidableEq beta] [DecidableEq move]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hstep : ∀ u, u ∈ moves -> Function.Involutive (step u))
    (hfixed : ∀ u, u ∈ moves -> ∀ x, step u x ≠ x) :
    2 * internalLineCount q step moves b +
        (crossHalf q step moves b).card =
      (fiber q b).card * moves.card := by
  exact unordered_balance_of_even_selfHalf q step moves b
    (selfHalf_card_even_of_involutive_fixedPointFree q step moves b hstep hfixed)

/-! ### Transport-row checksums -/

/-- Row-balance checksum for a target-bucket transport matrix once internal
oriented half-lines have even cardinality. -/
theorem transport_row_balance_of_even_selfHalf
    [Fintype alpha] [Fintype beta] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hself : Even (selfHalf q step moves b).card) :
    2 * internalLineCount q step moves b +
        (∑ c ∈ (Finset.univ.erase b),
          (transportHalf q step moves b c).card) =
      (fiber q b).card * moves.card := by
  classical
  rw [sum_transportHalf_card_offdiag_eq_crossHalf_card]
  exact unordered_balance_of_even_selfHalf q step moves b hself

/-- Row-balance checksum for fixed-point-free involutive move systems. -/
theorem transport_row_balance_of_involutive_fixedPointFree
    [Fintype alpha] [Fintype beta] [DecidableEq alpha] [DecidableEq beta]
    [DecidableEq move]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta)
    (hstep : ∀ u, u ∈ moves -> Function.Involutive (step u))
    (hfixed : ∀ u, u ∈ moves -> ∀ x, step u x ≠ x) :
    2 * internalLineCount q step moves b +
        (∑ c ∈ (Finset.univ.erase b),
          (transportHalf q step moves b c).card) =
      (fiber q b).card * moves.card := by
  classical
  exact transport_row_balance_of_even_selfHalf q step moves b
    (selfHalf_card_even_of_involutive_fixedPointFree q step moves b hstep hfixed)

/-! ### Boolean cube mask specialization -/

/-- A finite Boolean cube, used for tiling masks by taking
`index = StTile n`. -/
abbrev BoolCube (index : Type z) := index -> Bool

/-- Apply a Boolean mask to a cube point by coordinatewise xor. -/
def xorMask {index : Type z} (u x : BoolCube index) : BoolCube index :=
  fun i => Bool.xor (x i) (u i)

/-- A mask is nonzero when at least one coordinate is true. -/
def IsNonzeroMask {index : Type z} (u : BoolCube index) : Prop :=
  ∃ i, u i = true

@[simp] theorem xorMask_apply {index : Type z}
    (u x : BoolCube index) (i : index) :
    xorMask u x i = Bool.xor (x i) (u i) := rfl

/-- Xor by a fixed mask is an involution. -/
theorem xorMask_involutive {index : Type z} (u : BoolCube index) :
    Function.Involutive (xorMask u) := by
  intro x
  funext i
  cases hxi : x i <;> cases hui : u i <;> simp [xorMask, hxi, hui]

/-- Xor by a nonzero mask has no fixed point. -/
theorem xorMask_fixedPointFree_of_nonzero {index : Type z}
    {u : BoolCube index} (hu : IsNonzeroMask u) :
    ∀ x, xorMask u x ≠ x := by
  intro x h
  rcases hu with ⟨i, hi⟩
  have hcoord := congrFun h i
  rw [xorMask_apply, hi] at hcoord
  cases hx : x i <;> simp [hx] at hcoord

/-- The THM-346 Boolean-cube specialization: any bucket map out of a finite
Boolean cube satisfies the unordered balance for any finite family of nonzero
xor masks. -/
theorem unordered_balance_boolCube_masks
    {index : Type z} [Fintype index] [DecidableEq index]
    {beta : Type v} [DecidableEq beta]
    (q : BoolCube index -> beta) (moves : Finset (BoolCube index)) (b : beta)
    (hmoves : ∀ u, u ∈ moves -> IsNonzeroMask u) :
    2 * internalLineCount q xorMask moves b +
        (crossHalf q xorMask moves b).card =
      (fiber q b).card * moves.card := by
  classical
  refine unordered_balance_of_involutive_fixedPointFree q xorMask moves b ?_ ?_
  · intro u _hu
    exact xorMask_involutive u
  · intro u hu
    exact xorMask_fixedPointFree_of_nonzero (hmoves u hu)

/-- Target-bucket transport row checksum for finite Boolean cubes and nonzero
xor-mask families. -/
theorem transport_row_balance_boolCube_masks
    {index : Type z} [Fintype index] [DecidableEq index]
    {beta : Type v} [Fintype beta] [DecidableEq beta]
    (q : BoolCube index -> beta) (moves : Finset (BoolCube index)) (b : beta)
    (hmoves : ∀ u, u ∈ moves -> IsNonzeroMask u) :
    2 * internalLineCount q xorMask moves b +
        (∑ c ∈ (Finset.univ.erase b),
          (transportHalf q xorMask moves b c).card) =
      (fiber q b).card * moves.card := by
  classical
  refine transport_row_balance_of_involutive_fixedPointFree q xorMask moves b ?_ ?_
  · intro u _hu
    exact xorMask_involutive u
  · intro u hu
    exact xorMask_fixedPointFree_of_nonzero (hmoves u hu)

end

end BucketBalance
end Tournament
