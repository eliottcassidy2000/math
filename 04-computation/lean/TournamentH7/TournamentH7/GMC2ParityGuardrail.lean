import Mathlib

/-!
# A parity/Walsh guardrail for moment compression

The even- and odd-parity fibers of the Boolean 3-cube are the smallest
example in which all Walsh observables of degree at most two agree while the
top character changes sign.  Thus any GMC moment compression retaining only
mass, linear data, and pair data can lose a genuine signed cancellation
coordinate.

Everything here is an exact eight-point identity.  The proofs use ordinary
kernel reduction (`decide`), not `native_decide`.
-/

namespace GMC2.ParityGuardrail

/-- The Boolean 3-cube. -/
abbrev Cube3 := Fin 3 → Bool

/-- Addition in the Boolean cube. -/
def cubeXor (x y : Cube3) : Cube3 :=
  fun i => Bool.xor (x i) (y i)

/-- The parity label of a cube point: `false` is even and `true` is odd. -/
def parity (x : Cube3) : Bool :=
  Bool.xor (Bool.xor (x 0) (x 1)) (x 2)

/-- The counting weight of one parity fiber. -/
def parityWeight (ε : Bool) (x : Cube3) : ℤ :=
  if parity x = ε then 1 else 0

/-- The sign `(-1)^b` attached to one Boolean coordinate. -/
def bitSign (b : Bool) : ℤ :=
  if b then -1 else 1

/-- The Walsh character indexed by a coordinate set. -/
def walshCharacter (S : Finset (Fin 3)) (x : Cube3) : ℤ :=
  ∏ i ∈ S, bitSign (x i)

/-- The unnormalised Walsh transform of one parity-fiber counting measure. -/
def walshCoefficient (ε : Bool) (S : Finset (Fin 3)) : ℤ :=
  ∑ x : Cube3, parityWeight ε x * walshCharacter S x

/-- The squarefree Boolean monomial indexed by `S`. -/
def booleanMonomial (S : Finset (Fin 3)) (x : Cube3) : ℤ :=
  if ∀ i ∈ S, x i = true then 1 else 0

/-- The unnormalised squarefree moment of one parity fiber. -/
def booleanMoment (ε : Bool) (S : Finset (Fin 3)) : ℤ :=
  ∑ x : Cube3, parityWeight ε x * booleanMonomial S x

/-- Unnormalised integration against one parity-fiber counting measure. -/
def parityExpectation (ε : Bool) (f : Cube3 → ℤ) : ℤ :=
  ∑ x : Cube3, parityWeight ε x * f x

/-- The coordinate sets indexing Walsh characters of degree at most two. -/
def degreeTwoSets : Finset (Finset (Fin 3)) :=
  Finset.univ.filter (fun S => S.card ≤ 2)

/-- An arbitrary integer-valued linear combination of Walsh characters of
degree at most two. -/
def degreeTwoWalshObservable (c : Finset (Fin 3) → ℤ) (x : Cube3) : ℤ :=
  ∑ S ∈ degreeTwoSets, c S * walshCharacter S x

/-- XOR convolution of integer-valued functions on the Boolean cube. -/
def xorConvolution (f g : Cube3 → ℤ) (z : Cube3) : ℤ :=
  ∑ x : Cube3, f x * g (cubeXor x z)

/-- Both parity fibers have four points. -/
theorem parity_fiber_mass :
    walshCoefficient false ∅ = 4 ∧ walshCoefficient true ∅ = 4 := by
  decide

/-- The raw mass, singleton moments, and pair moments also agree.  This is the
coordinate form most directly consumed by moment/cut/Gram compressions. -/
theorem boolean_moments_agree_through_degree_two :
    ∀ S : Finset (Fin 3), S.card ≤ 2 →
      booleanMoment false S = booleanMoment true S := by
  decide

/-- The two fibers agree on every Walsh character of degree at most two. -/
theorem walsh_agree_through_degree_two :
    ∀ S : Finset (Fin 3), S.card ≤ 2 →
      walshCoefficient false S = walshCoefficient true S := by
  decide

/-- Consequently, the two parity fibers agree on every degree-at-most-two
Walsh polynomial, not only on the basis characters. -/
theorem degree_two_observables_agree (c : Finset (Fin 3) → ℤ) :
    parityExpectation false (degreeTwoWalshObservable c) =
      parityExpectation true (degreeTwoWalshObservable c) := by
  simp only [parityExpectation, degreeTwoWalshObservable, Finset.mul_sum]
  rw [Finset.sum_comm]
  conv_rhs => rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro S hS
  rw [← Finset.mul_sum, ← Finset.mul_sum]
  congr 1
  exact walsh_agree_through_degree_two S
    (Finset.mem_filter.mp hS).2

/-- In fact every nonconstant Walsh character below top degree vanishes on
both fibers. -/
theorem proper_nonempty_walsh_vanishes :
    ∀ S : Finset (Fin 3), S.Nonempty → S.card ≤ 2 →
      walshCoefficient false S = 0 ∧ walshCoefficient true S = 0 := by
  decide

/-- Complete Walsh spectrum of the two counting measures.  Only the constant
and cubic characters survive. -/
theorem complete_walsh_spectrum :
    ∀ ε : Bool, ∀ S : Finset (Fin 3),
      walshCoefficient ε S =
        if S = ∅ then 4
        else if S = Finset.univ then if ε then -4 else 4
        else 0 := by
  decide

/-- The cubic character is the missing coordinate: it is `+4` on the even
fiber and `-4` on the odd fiber. -/
theorem cubic_walsh_split :
    walshCoefficient false Finset.univ = 4 ∧
      walshCoefficient true Finset.univ = -4 := by
  decide

/-- In ordinary Boolean-moment coordinates, the missing cubic moment is `0`
on the even fiber and `1` on the odd fiber. -/
theorem cubic_boolean_moment_split :
    booleanMoment false Finset.univ = 0 ∧
      booleanMoment true Finset.univ = 1 := by
  decide

/-- The parity-fiber counting measures form the two-element group algebra up
to the fiber cardinality: `μ_ε * μ_η = 4 μ_(ε xor η)`. -/
theorem parity_weight_xor_convolution :
    ∀ ε η : Bool,
      xorConvolution (parityWeight ε) (parityWeight η) =
        fun z => 4 * parityWeight (Bool.xor ε η) z := by
  decide

end GMC2.ParityGuardrail

#print axioms GMC2.ParityGuardrail.parity_fiber_mass
#print axioms GMC2.ParityGuardrail.boolean_moments_agree_through_degree_two
#print axioms GMC2.ParityGuardrail.walsh_agree_through_degree_two
#print axioms GMC2.ParityGuardrail.degree_two_observables_agree
#print axioms GMC2.ParityGuardrail.proper_nonempty_walsh_vanishes
#print axioms GMC2.ParityGuardrail.complete_walsh_spectrum
#print axioms GMC2.ParityGuardrail.cubic_walsh_split
#print axioms GMC2.ParityGuardrail.cubic_boolean_moment_split
#print axioms GMC2.ParityGuardrail.parity_weight_xor_convolution
