set_option autoImplicit false

/-!
# A small polynomial normaliser for `Int` (core Lean, no Mathlib)

`pnorm` rewrites (the goal or the given hypotheses) into a sum of monomials: it expands
products, sorts the factors of every monomial by the commutativity/associativity lemmas
(simp's ordered rewriting), and then moves the numeral factors of every monomial to the
front, merging two numerals (pre-phase `lit_lit` and `Int.reduceMul`) before they could be
swapped. After `pnorm`, `omega` sees equal monomials as the same atom, so a polynomial
identity that is an integer linear combination of the (multiplied-out) hypotheses closes by
`omega`. Only `simp only` with equational lemmas and the numeral-evaluation simproc are
used, so no axioms beyond `propext` enter the proofs.
-/

namespace ProcgenSelfieEdim

theorem lit_mul_r (a : Int) (n : Nat) :
    a * (no_index (OfNat.ofNat n : Int)) = (OfNat.ofNat n : Int) * a := Int.mul_comm _ _

theorem lit_mul_m (a b : Int) (n : Nat) :
    a * ((no_index (OfNat.ofNat n : Int)) * b) = (OfNat.ofNat n : Int) * (a * b) :=
  Int.mul_left_comm _ _ _

theorem lit_lit (b : Int) (m n : Nat) :
    (no_index (OfNat.ofNat m : Int)) * ((no_index (OfNat.ofNat n : Int)) * b) =
      ((OfNat.ofNat m : Int) * (OfNat.ofNat n : Int)) * b := (Int.mul_assoc _ _ _).symm

/-- Casting a numeral. -/
theorem natCast_lit (n : Nat) :
    (((no_index (OfNat.ofNat n : Nat)) : Nat) : Int) = (OfNat.ofNat n : Int) := rfl

end ProcgenSelfieEdim

/-- Push `Nat → Int` casts through `+`, `*`, `^` and numerals. -/
syntax "pcast" (Lean.Parser.Tactic.location)? : tactic

macro_rules
  | `(tactic| pcast $[$loc]?) => `(tactic|
      simp only [Int.natCast_add, Int.natCast_mul, Int.natCast_pow, Int.natCast_one,
        Int.natCast_zero, ProcgenSelfieEdim.natCast_lit] $[$loc]?)

/-- Normalise integer polynomial expressions so that `omega` can finish (see the module doc). -/
syntax "pnorm" (Lean.Parser.Tactic.location)? : tactic

macro_rules
  | `(tactic| pnorm $[$loc]?) => `(tactic| (
      (try simp only [Int.mul_add, Int.add_mul, Int.mul_sub, Int.sub_mul, Int.mul_assoc,
        Int.mul_comm, Int.mul_left_comm, Int.mul_one, Int.one_mul, Int.mul_zero, Int.zero_mul,
        Int.mul_neg, Int.neg_mul, Int.neg_neg] $[$loc]?);
      (try simp only [↓ProcgenSelfieEdim.lit_lit, ↓Int.reduceMul, ProcgenSelfieEdim.lit_mul_r,
        ProcgenSelfieEdim.lit_mul_m] $[$loc]?)))
