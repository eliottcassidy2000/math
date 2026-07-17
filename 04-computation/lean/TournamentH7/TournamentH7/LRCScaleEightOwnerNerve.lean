/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: codex-c8-lean (LRC multi-agent project, 2026-07-17)
-/
import Mathlib

/-!
# The scale-eight owner-obligation obstruction

This file formalizes the terminal symbolic quotient in THM-963, equations
(17)--(23).  A signed doubling cycle is represented by a six-bit word `a`, with
odd total xor, and an order-eight unit word is represented by a six-digit
base-four word `x`.  The compact encodings `Fin 64` and `Fin 4096` avoid
materializing function spaces during the finite audit.  `ownerSat a x i` is
exactly the parity-complement plus four-symbol all-different predicate at owner
`i`.

The finite theorem `reducedDistanceTwoCore_false` shows that obligations two
steps apart on the six-cycle are disjoint after forgetting the two private
digits which occur in only one obligation each.  A separate six-vertex
combinatorial lemma then derives that no three distinct owner obligations, and
hence no six, can hold simultaneously.  This is precisely the terminal
obstruction supplied by the same-parity parts of the `K_{3,3}` nerve; proving
that the remaining pair intersections are nonempty is not needed for emptiness
and is not claimed here.

All finite tables use ordinary kernel `decide`.  This file does **not**
formalize the preceding literal-mask reduction from all scale-eight packets to
this quotient.  There is no `sorry`, no `native_decide`, and no added axiom.
-/

namespace LonelyRunner
namespace ScaleEightOwnerNerve

/-- Positions on the signed six-cycle. -/
abbrev Idx := Fin 6

/-- The base-four digit `x` in the odd unit `e = 2x + 1`. -/
abbrev Digit := Fin 4

/-- A six-bit sign word, with coordinate zero in the least significant bit. -/
abbrev SignWord := Fin 64

/-- A six-digit base-four unit word, with coordinate zero least significant. -/
abbrev UnitWord := Fin 4096

/-- Cyclic shift on the six owner positions. -/
def shift (i : Idx) (d : Nat) : Idx :=
  ⟨((i : Nat) + d) % 6, by omega⟩

/-- Read bit `i` of the compact sign-word encoding. -/
def signAt (a : SignWord) (i : Idx) : Bool :=
  (a.val / 2 ^ i.val) % 2 == 1

/-- Read base-four digit `i` of the compact unit-word encoding. -/
def digitAt (x : UnitWord) (i : Idx) : Digit :=
  ⟨(x.val / 4 ^ i.val) % 4, by omega⟩

/-- Read a sign bit at cyclic displacement `d` from `i`. -/
def signOffset (a : SignWord) (i : Idx) (d : Nat) : Bool :=
  signAt a (shift i d)

/-- Read a base-four digit at cyclic displacement `d` from `i`. -/
def digitOffset (x : UnitWord) (i : Idx) (d : Nat) : Digit :=
  digitAt x (shift i d)

/-- The xor `A_(i,d)` of the first `d` edge signs from owner `i`. -/
def prefixXor (a : SignWord) (i : Idx) : Nat → Bool
  | 0 => false
  | d + 1 => prefixXor a i d != signOffset a i d

/-- Equation (17): the signed doubling cycle has odd total sign. -/
def oddSign (a : SignWord) : Bool :=
  prefixXor a 0 6

/-- Reduction of a base-four digit modulo two, encoded as a bit. -/
def parityBit (x : Digit) : Bool :=
  x.val % 2 == 1

/-- Addition by a natural number in `Z/4Z`. -/
def add4 (x : Digit) (k : Nat) : Digit :=
  ⟨(x.val + k) % 4, by omega⟩

/-- The residue `k - x` in `Z/4Z`. -/
def sub4 (k : Nat) (x : Digit) : Digit :=
  ⟨(k + 4 - x.val) % 4, by omega⟩

/-- Four base-four symbols are all different.  Since there are exactly four
symbols, this is equivalent to their set being all of `Z/4Z`. -/
def allDifferent4 (q0 q1 q2 q3 : Digit) : Bool :=
  (q0 != q1) && (q0 != q2) && (q0 != q3) &&
  (q1 != q2) && (q1 != q3) && (q2 != q3)

/-- Equations (21)--(23), specialized at owner `i`.

The two Boolean symbols are the odd sheets.  The four `q` symbols are the
nonfixed even sheets after division by two. -/
def ownerSat (a : SignWord) (x : UnitWord) (i : Idx) : Bool :=
  let A1 := prefixXor a i 1
  let A2 := prefixXor a i 2
  let A3 := prefixXor a i 3
  let A4 := prefixXor a i 4
  let b2 := parityBit (digitOffset x i 2) != A2
  let b4 := parityBit (digitOffset x i 4) != A4
  let q0 := sub4 1 (digitOffset x i 0)
  let q1 := if A1 then add4 (digitOffset x i 1) 2 else sub4 1 (digitOffset x i 1)
  let q2 := if A2 then sub4 3 (digitOffset x i 2) else digitOffset x i 2
  let q3 := if A3 then digitOffset x i 3 else sub4 3 (digitOffset x i 3)
  (b2 != b4) && allDifferent4 q0 q1 q2 q3

/-- Three symbols are pairwise distinct. -/
def allDifferent3 (q0 q1 q2 : Digit) : Bool :=
  (q0 != q1) && (q0 != q2) && (q1 != q2)

/-- Xor of a six-bit function, written in the same left-associated form as
`prefixXor`. -/
def xorSix (b : Idx → Bool) : Bool :=
  ((((((false != b 0) != b 1) != b 2) != b 3) != b 4) != b 5)

/-- The reduced local truth table for owners `0` and `2`.

The first all-different constraint discards `q₁`, and the second discards
`q₃`.  Those symbols depend only on `x₁` and `x₅`, respectively.  Thus owner
satisfaction implies this four-digit predicate, while the discarded symbols
are irrelevant to the contradiction. -/
def reducedDistanceTwoCore
    (b : Idx → Bool) (x0 x2 x3 x4 : Digit) : Bool :=
  let A1 := false != b 0
  let A2 := A1 != b 1
  let A3 := A2 != b 2
  let A4 := A3 != b 3
  let b2 := parityBit x2 != A2
  let b4 := parityBit x4 != A4
  let q0 := sub4 1 x0
  let q2 := if A2 then sub4 3 x2 else x2
  let q3 := if A3 then x3 else sub4 3 x3
  let C1 := false != b 2
  let C2 := C1 != b 3
  let C3 := C2 != b 4
  let C4 := C3 != b 5
  let c2 := parityBit x4 != C2
  let c4 := parityBit x0 != C4
  let r0 := sub4 1 x2
  let r1 := if C1 then add4 x3 2 else sub4 1 x3
  let r2 := if C2 then sub4 3 x4 else x4
  (b2 != b4) && allDifferent3 q0 q2 q3 &&
    (c2 != c4) && allDifferent3 r0 r1 r2

/-- The terminal obstruction needs only `2^6 * 4^4 = 16,384` rows.  Ordinary
kernel `decide` closes this reduced table; no native-code axiom is involved. -/
theorem reducedDistanceTwoCore_false :
    ∀ (b : Idx → Bool) (x0 x2 x3 x4 : Digit),
      xorSix b = true → reducedDistanceTwoCore b x0 x2 x3 x4 = false := by
  decide

/-- Cyclic rotation does not change the xor of the six sign bits.  This tiny
`64 * 6` table is checked by ordinary kernel `decide`. -/
theorem xorSix_signOffset :
    ∀ (a : SignWord) (i : Idx),
      xorSix (fun d => signOffset a i d.val) = oddSign a := by
  decide

/-- Two full owner obligations imply the reduced four-digit core.  The proof
only forgets comparisons involving `q₁` at the first owner and `q₃` at the
second owner. -/
theorem ownerPair_implies_reducedCore
    (a : SignWord) (x : UnitWord) (i : Idx)
    (hi : ownerSat a x i = true)
    (hj : ownerSat a x (shift i 2) = true) :
    reducedDistanceTwoCore (fun d => signOffset a i d.val)
      (digitOffset x i 0) (digitOffset x i 2)
      (digitOffset x i 3) (digitOffset x i 4) = true := by
  fin_cases i <;>
    simp [ownerSat, reducedDistanceTwoCore, allDifferent4, allDifferent3,
      prefixXor, shift, signOffset, digitOffset] at hi hj ⊢ <;>
    aesop

/-- Owner obligations at cyclic distance two are disjoint.  This is the
finite terminal algebraic fact behind the same-parity side of the `K_{3,3}`
nerve in THM-963. -/
theorem distanceTwo_owner_disjoint
    (a : SignWord) (x : UnitWord) (i : Idx)
    (ha : oddSign a = true)
    (hi : ownerSat a x i = true) :
    ownerSat a x (shift i 2) = false := by
  cases hj : ownerSat a x (shift i 2) with
  | false => rfl
  | true =>
      have hcore := ownerPair_implies_reducedCore a x i hi hj
      have hxor : xorSix (fun d => signOffset a i d.val) = true := by
        rw [xorSix_signOffset, ha]
      have hfalse := reducedDistanceTwoCore_false
        (fun d => signOffset a i d.val)
        (digitOffset x i 0) (digitOffset x i 2)
        (digitOffset x i 3) (digitOffset x i 4) hxor
      simp [hfalse] at hcore

/-- Any three distinct positions on a six-cycle contain an ordered pair at
cyclic distance two.  This is the finite pigeonhole step: two positions have
the same parity. -/
theorem three_distinct_have_distanceTwo :
    ∀ i j k : Idx,
      i ≠ j → i ≠ k → j ≠ k →
      (j = shift i 2 ∨ i = shift j 2 ∨
       k = shift i 2 ∨ i = shift k 2 ∨
       k = shift j 2 ∨ j = shift k 2) := by
  decide

/-- No base-four unit word can satisfy three distinct owner obligations. -/
theorem no_three_distinct_owners
    (a : SignWord) (x : UnitWord)
    (ha : oddSign a = true)
    (i j k : Idx) (hij : i ≠ j) (hik : i ≠ k) (hjk : j ≠ k)
    (hi : ownerSat a x i = true)
    (hj : ownerSat a x j = true)
    (hk : ownerSat a x k = true) : False := by
  rcases three_distinct_have_distanceTwo i j k hij hik hjk with
    hji | hij' | hki | hik' | hkj | hjk'
  · subst j
    have hfalse := distanceTwo_owner_disjoint a x i ha hi
    simp [hfalse] at hj
  · subst i
    have hfalse := distanceTwo_owner_disjoint a x j ha hj
    simp [hfalse] at hi
  · subst k
    have hfalse := distanceTwo_owner_disjoint a x i ha hi
    simp [hfalse] at hk
  · subst i
    have hfalse := distanceTwo_owner_disjoint a x k ha hk
    simp [hfalse] at hi
  · subst k
    have hfalse := distanceTwo_owner_disjoint a x j ha hj
    simp [hfalse] at hk
  · subst j
    have hfalse := distanceTwo_owner_disjoint a x k ha hk
    simp [hfalse] at hj

/-- In particular, the intersection of all six owner obligations is empty. -/
theorem no_sixfold_owner_intersection
    (a : SignWord) (x : UnitWord)
    (ha : oddSign a = true)
    (hall : ∀ i : Idx, ownerSat a x i = true) : False := by
  exact no_three_distinct_owners a x ha 0 1 2
    (by decide) (by decide) (by decide) (hall 0) (hall 1) (hall 2)

#print axioms reducedDistanceTwoCore_false
#print axioms distanceTwo_owner_disjoint
#print axioms no_three_distinct_owners
#print axioms no_sixfold_owner_intersection

end ScaleEightOwnerNerve
end LonelyRunner
