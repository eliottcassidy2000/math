/-
  TournamentH7.LRCModularCuspLedger -- S247 q-Pochhammer / modular-cusp
  Lean-facing ledger for the LRC14 proof stack.

  This module is deliberately small.  It does not try to formalize the analytic
  theory of modular forms.  Instead it records the interface that a q-series
  sidecar must satisfy before it is allowed to act as an LRC14 proof packet:

    * a raw q-Pochhammer tail is only a product-tail datum;
    * a full modular-function exit must carry a finite Laurent principal part at
      the cusp infinity;
    * a Hurwitz-style limiting step must preserve the finite-principal-part
      certificate rather than silently changing zero/pole data;
    * the two sixth-power collision equations are related by a proved face map,
      so the two-term equation can be reused inside the three-term ledger.

  The analytic notions are exposed as named obligations.  Theorems below are
  sorry-free finite-tail and arithmetic glue.
-/

import Mathlib.Tactic

namespace LonelyRunner
namespace LRC14
namespace ModularCuspLedger

/-! ## Laurent principal parts at the cusp -/

/-- A q-expansion has only finitely many negative powers if all coefficients
below some pole order vanish. -/
def HasOnlyFiniteNegativePowers (coeff : Int -> Int) : Prop :=
  ∃ poleOrder : Nat, ∀ n : Int, n < -((poleOrder : Int)) -> coeff n = 0

/-- A q-expansion has no negative powers. -/
def HasNoNegativePowers (coeff : Int -> Int) : Prop :=
  ∀ n : Int, n < 0 -> coeff n = 0

/-- Holomorphic q-expansions are a special case of finite-principal-part
q-expansions. -/
theorem finite_negative_powers_of_no_negative_powers
    {coeff : Int -> Int} (h : HasNoNegativePowers coeff) :
    HasOnlyFiniteNegativePowers coeff := by
  refine ⟨0, ?_⟩
  intro n hn
  exact h n hn

/-- A proof-bearing finite-principal-part packet. -/
structure LaurentPrincipalPartPacket where
  coeff : Int -> Int
  poleOrder : Nat
  coeff_vanishes_below_pole :
    ∀ n : Int, n < -((poleOrder : Int)) -> coeff n = 0

/-- The packet form immediately gives the finite-negative-tail property. -/
theorem packet_has_only_finite_negative_powers
    (packet : LaurentPrincipalPartPacket) :
    HasOnlyFiniteNegativePowers packet.coeff :=
  ⟨packet.poleOrder, packet.coeff_vanishes_below_pole⟩

/-! ## Full modular-function and Hurwitz obligations -/

/-- Placeholder for the full `SL2(Z)` invariance obligation. -/
opaque IsFullModularGroupInvariant : (Int -> Int) -> Prop

/-- Placeholder for meromorphicity at the cusp `i infinity`. -/
opaque IsMeromorphicAtCuspInfinity : (Int -> Int) -> Prop

/-- The q-expansion theorem needed from the modular-function side: full modular
invariance plus meromorphicity at the cusp must produce a finite principal
part.  This is a named obligation, not asserted here. -/
def FullModularCuspExpansionObligation : Prop :=
  ∀ coeff : Int -> Int,
    IsFullModularGroupInvariant coeff ->
    IsMeromorphicAtCuspInfinity coeff ->
    HasOnlyFiniteNegativePowers coeff

/-- Once the modular-function theorem is supplied, its Lean-facing payload is
exactly the finite-negative-tail certificate. -/
theorem finite_negative_tail_from_full_modular_cusp
    (hmod : FullModularCuspExpansionObligation)
    {coeff : Int -> Int}
    (hinv : IsFullModularGroupInvariant coeff)
    (hmer : IsMeromorphicAtCuspInfinity coeff) :
    HasOnlyFiniteNegativePowers coeff :=
  hmod coeff hinv hmer

/-- A Hurwitz-style gate for q-expansion limits: the complex-analysis step is
allowed to identify a limit series with a certified series, and Lean then
transfers the finite-principal-part certificate. -/
structure HurwitzQExpansionGate (limit certified : Int -> Int) where
  coefficient_identity : ∀ n : Int, limit n = certified n
  certified_finite_tail : HasOnlyFiniteNegativePowers certified

/-- The finite-principal-part property is stable through a Hurwitz gate. -/
theorem finite_tail_of_hurwitz_gate
    {limit certified : Int -> Int}
    (gate : HurwitzQExpansionGate limit certified) :
    HasOnlyFiniteNegativePowers limit := by
  rcases gate.certified_finite_tail with ⟨poleOrder, htail⟩
  refine ⟨poleOrder, ?_⟩
  intro n hn
  rw [gate.coefficient_identity n]
  exact htail n hn

/-! ## S247 exact q-series readouts -/

/-- The S247 Euler-pentagonal truncation of `(q;q)_inf` through `q^12`.  The
initial `n < 0` branch makes the no-negative-powers proof tautological. -/
def qPochhammerEulerCoeffTo12 (n : Int) : Int :=
  if n < 0 then 0
  else if n = 0 then 1
  else if n = 1 then -1
  else if n = 2 then -1
  else if n = 5 then 1
  else if n = 7 then 1
  else if n = 12 then -1
  else 0

/-- The q-Pochhammer product tail has no negative q-powers. -/
theorem qPochhammerEulerCoeffTo12_no_negative_powers :
    HasNoNegativePowers qPochhammerEulerCoeffTo12 := by
  intro n hn
  simp [qPochhammerEulerCoeffTo12, hn]

/-- The S247 q-Pochhammer truncation is therefore legal as a holomorphic tail,
but not by itself a full modular-function certificate. -/
theorem qPochhammerEulerCoeffTo12_finite_negative_tail :
    HasOnlyFiniteNegativePowers qPochhammerEulerCoeffTo12 :=
  finite_negative_powers_of_no_negative_powers
    qPochhammerEulerCoeffTo12_no_negative_powers

/-- Principal part of the `j`-invariant readout in S247: only `q^-1`. -/
def jPrincipalCoeffS247 (n : Int) : Int :=
  if n < -1 then 0 else if n = -1 then 1 else 0

/-- The S247 `j` principal part has pole order one at the cusp. -/
theorem jPrincipalCoeffS247_finite_negative_tail :
    HasOnlyFiniteNegativePowers jPrincipalCoeffS247 := by
  refine ⟨1, ?_⟩
  intro n hn
  have hlt : n < (-1 : Int) := by omega
  simp [jPrincipalCoeffS247, hlt]

/-- Principal part of the `1/Delta` readout in S247: again only `q^-1`. -/
def invDeltaPrincipalCoeffS247 (n : Int) : Int :=
  if n < -1 then 0 else if n = -1 then 1 else 0

/-- The S247 `1/Delta` principal part has pole order one at the cusp. -/
theorem invDeltaPrincipalCoeffS247_finite_negative_tail :
    HasOnlyFiniteNegativePowers invDeltaPrincipalCoeffS247 := by
  refine ⟨1, ?_⟩
  intro n hn
  have hlt : n < (-1 : Int) := by omega
  simp [invDeltaPrincipalCoeffS247, hlt]

/-! ## Sixth-power collision ledger -/

/-- The two-term sixth-power collision equation
`a^6 + b^6 = d^6 + e^6`. -/
structure SixthPowerTwoCollision where
  a : Nat
  b : Nat
  d : Nat
  e : Nat
  equation : a ^ 6 + b ^ 6 = d ^ 6 + e ^ 6

/-- The three-term sixth-power collision equation
`a^6 + b^6 + c^6 = d^6 + e^6 + f^6`. -/
structure SixthPowerThreeCollision where
  a : Nat
  b : Nat
  c : Nat
  d : Nat
  e : Nat
  f : Nat
  equation : (a ^ 6 + b ^ 6) + c ^ 6 = (d ^ 6 + e ^ 6) + f ^ 6

/-- A two-term collision is a face of the three-term collision ledger by adding
the same sixth-power tail to both sides. -/
def twoCollisionToThreeWithTail
    (collision : SixthPowerTwoCollision) (tail : Nat) :
    SixthPowerThreeCollision where
  a := collision.a
  b := collision.b
  c := tail
  d := collision.d
  e := collision.e
  f := tail
  equation := by
    exact congrArg (fun x : Nat => x + tail ^ 6) collision.equation

/-- In particular, the zero-tail embedding is available without any extra
number-theoretic assumption. -/
def twoCollisionToThreeZeroTail
    (collision : SixthPowerTwoCollision) :
    SixthPowerThreeCollision :=
  twoCollisionToThreeWithTail collision 0

/-- Collision arities retained by the sidecar. -/
inductive SixthPowerCollisionArity where
  | twoTerm
  | threeTerm
deriving DecidableEq, Repr, Fintype

/-- Number of variables retained by each collision equation. -/
def retainedCollisionVariables : SixthPowerCollisionArity -> Nat
  | .twoTerm => 4
  | .threeTerm => 6

/-- The two-term equation is strictly smaller but embeds into the three-term
ledger. -/
theorem twoTerm_retained_variables_lt_threeTerm :
    retainedCollisionVariables .twoTerm <
      retainedCollisionVariables .threeTerm := by
  native_decide

/-! ## Tournament-analysis sidecar vertices -/

/-- S247 uses proof sidecars as tournament vertices, not runners or raw
q-coefficients. -/
inductive ModularSidecarVertex where
  | fullModularFunctionPacket
  | jRationalExit
  | finitePrincipalPartLedger
  | etaMultiplierBalance
  | hurwitzZeroPersistenceGate
  | qPochhammerTail
  | rawQCoefficients
deriving DecidableEq, Repr, Fintype

/-- The S247 sidecar tournament has seven proof-obligation vertices. -/
theorem modularSidecarVertex_count :
    Fintype.card ModularSidecarVertex = 7 := by
  native_decide

/-! ### Axiom audit -/

#print axioms packet_has_only_finite_negative_powers
#print axioms finite_negative_tail_from_full_modular_cusp
#print axioms finite_tail_of_hurwitz_gate
#print axioms qPochhammerEulerCoeffTo12_finite_negative_tail
#print axioms jPrincipalCoeffS247_finite_negative_tail
#print axioms invDeltaPrincipalCoeffS247_finite_negative_tail
#print axioms twoTerm_retained_variables_lt_threeTerm
#print axioms modularSidecarVertex_count

end ModularCuspLedger
end LRC14
end LonelyRunner
