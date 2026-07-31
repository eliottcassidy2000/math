import Mathlib.Algebra.Group.Fin.Basic
import Mathlib.Data.ZMod.Basic
import Mathlib.Logic.Equiv.Fin.Basic
import Mathlib.Tactic.Ring
import Lean.Elab.Tactic.Omega

open scoped Fin.NatCast

/-!
# Finite algebra behind THM-2878 and THM-2882

This scratch module isolates the theorem-pure part of the event-twisted
coefficient lift:

* the `u3 : D → S` event count is the base-thirteen carry;
* `Fin 169` has the nonsplit set chart `(a,q) ↦ 13a+q`;
* reduced transport is projective by one vertical carry, while the lifted
  transport composes honestly;
* the `q3 → q11 → q7` reduced seam has exponent-three clutch, and the
  inverse vertical edge closes it with exponent zero;
* the two-origin selector XOR coordinate is independent of the carry.

The phase is represented by its exponent in `ZMod 13`.  Therefore the
statement `seamDefect = 3` is exactly the algebraic content of the
cyclotomic statement `seamCurvature = ω^3`, without importing a particular
field or primitive root.
-/

namespace LonelyRunner
namespace C169Carry

/-- A base-thirteen residue. -/
abbrev Residue := Fin 13

/-- The base-thirteen chart coordinates `(ancestry, residue)`.  This is a
set chart; its coordinatewise addition is not the cyclic group law. -/
abbrev State := Residue × Residue

/-- Character exponents for a primitive thirteenth root. -/
abbrev Phase := ZMod 13

/-- The exact set equivalence `(a,q) ↦ 13a+q`. -/
def c169Chart : State ≃ Fin 169 :=
  finProdFinEquiv

/-- The base-thirteen carry `κ(q,h)=⌊(q+h)/13⌋`. -/
def kappa (q h : Residue) : ℕ :=
  (q.val + h.val) / 13

/-- Addition in the reduced residue chart. -/
def advance (q h : Residue) : Residue :=
  q + h

/-- The literal positive `u3 : D → S` event is based at residue `12`. -/
def u3Exit (q : Residue) : ℕ :=
  if q = 12 then 1 else 0

/-- Count positive `u3` exits along the `h` unit edges beginning at `q`. -/
def positiveExitCount (q h : Residue) : ℕ :=
  ∑ i ∈ Finset.range h.val, u3Exit (q + (i : Residue))

/-- Vertical ancestry translation in the base-thirteen chart. -/
def vertical (d : ℕ) (x : State) : State :=
  (x.1 + (d : Residue), x.2)

/-- The nonsplit lift of a reduced step. -/
def lift (h : Residue) (x : State) : State :=
  (x.1 + (kappa x.2 h : Residue), advance x.2 h)

/-- The event count is exactly the floor carry on every one of the 169
reduced paths. -/
theorem positiveExitCount_eq_kappa :
    ∀ q h : Residue, positiveExitCount q h = kappa q h := by
  decide

/-- Since both inputs are residues, the carry is Boolean-valued. -/
theorem kappa_eq_zero_or_one (q h : Residue) :
    kappa q h = 0 ∨ kappa q h = 1 := by
  unfold kappa
  have hq := q.isLt
  have hh := h.isLt
  omega

/-- Exact threshold characterization of a carry. -/
theorem kappa_eq_one_iff (q h : Residue) :
    kappa q h = 1 ↔ 13 ≤ q.val + h.val := by
  unfold kappa
  have hq := q.isLt
  have hh := h.isLt
  omega

/-- The chart really is the advertised base-thirteen numeral chart. -/
theorem c169Chart_val (a q : Residue) :
    (c169Chart (a, q)).val = 13 * a.val + q.val := by
  simp [c169Chart, finProdFinEquiv]
  omega

/-- Transport upstairs is literal addition by the reduced step in `Fin 169`. -/
theorem c169Chart_lift :
    ∀ x : State, ∀ h : Residue,
      c169Chart (lift h x) = c169Chart x + (h.val : Fin 169) := by
  rintro ⟨a, q⟩ h
  apply Fin.ext
  rw [c169Chart_val, Fin.val_add, c169Chart_val]
  simp only [lift, kappa, advance, Fin.val_add, Fin.val_natCast]
  omega

/-- Two lifted steps compose as honest addition in the cyclic 169-state
carrier.  The vertical term in `lift_comp` appears only after the total
step is reduced back to one base-thirteen digit. -/
theorem c169Chart_lift_two (x : State) (h k : Residue) :
    c169Chart (lift k (lift h x)) =
      (c169Chart x + (h.val : Fin 169)) + (k.val : Fin 169) := by
  rw [c169Chart_lift, c169Chart_lift]

/-- The carry witness at the end of the base-thirteen digit. -/
theorem lift_one_at_twelve :
    lift 1 ((0 : Residue), (12 : Residue)) =
      ((1 : Residue), (0 : Residue)) := by
  decide

/-- The base-thirteen set chart is not an additive product chart.  The
single witness is the same `12+1` carry as `lift_one_at_twelve`. -/
theorem c169Chart_product_add_failure :
    c169Chart
        (((0 : Residue), (12 : Residue)) +
          ((0 : Residue), (1 : Residue))) ≠
      c169Chart ((0 : Residue), (12 : Residue)) +
        c169Chart ((0 : Residue), (1 : Residue)) := by
  decide

/-- Consequently, the set equivalence cannot be misread as the
coordinatewise product-group identification. -/
theorem c169Chart_not_product_additive :
    ¬ ∀ x y : State, c169Chart (x + y) = c169Chart x + c169Chart y := by
  intro h
  exact c169Chart_product_add_failure
    (h ((0 : Residue), (12 : Residue))
      ((0 : Residue), (1 : Residue)))

/-- The exact base-thirteen carry cocycle.  The last term is the carry of
adding the two step digits themselves. -/
theorem kappa_cocycle :
    ∀ q h k : Residue,
      kappa q h + kappa (advance q h) k =
        kappa q (advance h k) + kappa h k := by
  intro q h k
  simp only [kappa, advance, Fin.val_add]
  omega

/-- Reduced step composition acquires exactly one vertical carry. -/
theorem lift_comp :
    ∀ x : State, ∀ h k : Residue,
      lift k (lift h x) =
        vertical (kappa h k) (lift (advance h k) x) := by
  rintro ⟨a, q⟩ h k
  apply Prod.ext
  · simp only [lift, vertical]
    rw [add_assoc, add_assoc]
    congr 1
    rw [← Nat.cast_add, ← Nat.cast_add, kappa_cocycle]
  · simp only [lift, vertical, advance]
    exact add_assoc q h k

/-! ## Exponent-valued coefficient transport -/

/-- The old flat frequency exponent. -/
def flatExp (h : Residue) : Phase :=
  3 * (h.val : Phase)

/-- The event exponent contributed by the `u3` carry. -/
def eventExp (q h : Residue) : Phase :=
  3 * (kappa q h : Phase)

/-- The event-twisted charged-channel exponent. -/
def twistedExp (q h : Residue) : Phase :=
  flatExp h + eventExp q h

/-- A vertical ancestry step has character-three exponent. -/
def verticalExp (d : ℕ) : Phase :=
  3 * (d : Phase)

/-- Reduction modulo thirteen in the residue chart agrees with reduction
in the exponent group. -/
theorem phase_advance (h k : Residue) :
    ((advance h k).val : Phase) =
      (h.val : Phase) + (k.val : Phase) := by
  simp only [advance, Fin.val_add, ZMod.natCast_mod, Nat.cast_add]

/-- Exponent-valued form of the carry cocycle. -/
theorem phase_kappa_cocycle (q h k : Residue) :
    (kappa q h : Phase) + (kappa (advance q h) k : Phase) =
      (kappa q (advance h k) : Phase) + (kappa h k : Phase) := by
  rw [← Nat.cast_add, ← Nat.cast_add, kappa_cocycle]

/-- The reduced coefficient law is projective by exactly the vertical
carry exponent. -/
theorem twistedExp_comp :
    ∀ q h k : Residue,
      twistedExp (advance q h) k + twistedExp q h =
        verticalExp (kappa h k) + twistedExp q (advance h k) := by
  intro q h k
  calc
    twistedExp (advance q h) k + twistedExp q h =
        3 * ((h.val : Phase) + (k.val : Phase)) +
          3 * ((kappa q h : Phase) +
            (kappa (advance q h) k : Phase)) := by
      simp only [twistedExp, flatExp, eventExp]
      ring
    _ = 3 * ((advance h k).val : Phase) +
          3 * ((kappa q (advance h k) : Phase) +
            (kappa h k : Phase)) := by
      rw [phase_kappa_cocycle, phase_advance]
    _ = verticalExp (kappa h k) + twistedExp q (advance h k) := by
      simp only [verticalExp, twistedExp, flatExp, eventExp]
      ring

/-- The seam route `q3 --8→ q11 --9→ q7` has event counts `(0,1)`;
the direct `q3 --4→ q7` route has event count `0`. -/
theorem seam_event_counts :
    kappa 3 8 = 0 ∧ kappa 11 9 = 1 ∧ kappa 3 4 = 0 := by
  decide

/-- The flat part has no reduced seam curvature. -/
theorem flat_seam_defect :
    flatExp 9 + flatExp 8 - flatExp 4 = 0 := by
  decide

/-- The event part carries the character-three seam defect. -/
theorem event_seam_defect :
    eventExp 11 9 + eventExp 3 8 - eventExp 3 4 = 3 := by
  decide

/-- The full reduced seam has exponent `3`, i.e. curvature `ω^3`. -/
theorem twisted_seam_defect :
    twistedExp 11 9 + twistedExp 3 8 - twistedExp 3 4 = 3 := by
  decide

/-- The via route and direct route have different lifted endpoints. -/
theorem seam_lifted_endpoints (a : Residue) :
    lift 9 (lift 8 (a, 3)) = (a + 1, 7) ∧
      lift 4 (a, 3) = (a, 7) := by
  decide +revert

/-- The endpoint discrepancy is precisely one positive vertical step. -/
theorem seam_endpoint_clutch (a : Residue) :
    lift 9 (lift 8 (a, 3)) = vertical 1 (lift 4 (a, 3)) := by
  decide +revert

/-- The inverse vertical step closes the two lifted seam routes at the
state level.  In `Fin 13`, adding `12` is the inverse of adding `1`. -/
theorem closed_lifted_endpoint (a : Residue) :
    vertical 12 (lift 9 (lift 8 (a, 3))) = lift 4 (a, 3) := by
  decide +revert

/-- Appending the inverse vertical clutch closes the exponent-valued loop.
This is the honest flatness statement upstairs. -/
theorem closed_lifted_loop_flat :
    (twistedExp 11 9 + twistedExp 3 8 - twistedExp 3 4) + (-3) = 0 := by
  decide

/-! ## Two-origin selector XOR calculus -/

/-- Indicator of the exceptional residue `q=3`. -/
def delta3 (q : Residue) : Bool :=
  q == 3

/-- E3 truth at the zero endpoint origin: `{0,3,11}`. -/
def truth0 (q : Residue) : Bool :=
  q == 0 || q == 3 || q == 11

/-- E3 truth at the stepped endpoint origin: `{0,11}`. -/
def truth1 (q : Residue) : Bool :=
  q == 0 || q == 11

/-- One origin contributes exactly when the selected block is the true block. -/
def contributes (selected truth : Bool) : Bool :=
  selected == truth

/-- A signed two-origin coefficient is nonzero exactly when one origin,
but not both, contributes. -/
def selectorNonzero (s0 s1 t0 t1 : Bool) : Bool :=
  Bool.xor (contributes s0 t0) (contributes s1 t1)

/-- Pure Boolean form of the selector law:
`nonzero ↔ s0 XOR s1 XOR t0 XOR t1`. -/
theorem selector_xor_law (s0 s1 t0 t1 : Bool) :
    selectorNonzero s0 s1 t0 t1 =
      Bool.xor (Bool.xor s0 s1) (Bool.xor t0 t1) := by
  decide +revert

/-- The two truth labels differ only at residue `3`. -/
theorem truth_origin_parity :
    ∀ q : Residue, Bool.xor (truth0 q) (truth1 q) = delta3 q := by
  decide

/-- The unique selector that keeps the zero-origin atom and rejects the
stepped-origin atom. -/
def positiveSelector (q : Residue) : Bool × Bool :=
  (truth0 q, Bool.not (truth1 q))

/-- This selector is indeed nonzero at every residue. -/
theorem positiveSelector_nonzero :
    ∀ q : Residue,
      selectorNonzero (positiveSelector q).1 (positiveSelector q).2
        (truth0 q) (truth1 q) = true := by
  decide

/-- Pointwise uniqueness of the zero-origin-positive selector. -/
theorem positiveSelector_unique :
    ∀ q : Residue, ∀ s0 s1 : Bool,
      (contributes s0 (truth0 q) = true ∧
          contributes s1 (truth1 q) = false) ↔
        (s0, s1) = positiveSelector q := by
  decide

/-- The origin-parity coordinate of the positive selector. -/
def selectorParity (q : Residue) : Bool :=
  Bool.xor (positiveSelector q).1 (positiveSelector q).2

/-- Exact parity law `p(q)=1 XOR δ(q,3)=not δ(q,3)`. -/
theorem selectorParity_eq_not_delta3 :
    ∀ q : Residue, selectorParity q = Bool.not (delta3 q) := by
  decide

/-- Selector parity change across a general reduced step. -/
def selectorToggle (q h : Residue) : Bool :=
  Bool.xor (selectorParity q) (selectorParity (advance q h))

/-- On unit edges, selector parity toggles exactly at `q2→q3` and
`q3→q4`. -/
theorem selector_unit_toggle_support :
    ∀ q : Residue,
      selectorToggle q 1 = (q == 2 || q == 3) := by
  decide

/-- The unique carry edge has constant selector parity. -/
theorem carry_edge_selector_constant :
    selectorParity 12 = selectorParity 0 ∧ kappa 12 1 = 1 := by
  decide

/-- Conversely, the `q0→q3` selector change has no carry. -/
theorem selector_change_without_carry :
    selectorParity 0 ≠ selectorParity 3 ∧ kappa 0 3 = 0 := by
  decide

/-- Even the complete selector edge signature does not determine carry:
the `q12→q0` and `q0→q1` signatures agree, but only the former carries. -/
theorem same_selector_signature_different_carry :
    (selectorParity 12, selectorParity (advance 12 1)) =
        (selectorParity 0, selectorParity (advance 0 1)) ∧
      kappa 12 1 ≠ kappa 0 1 := by
  decide

/-- Carry does not determine selector change: two carry-zero edges have
different selector toggles. -/
theorem same_carry_different_selector_toggle :
    kappa 0 3 = kappa 0 1 ∧
      selectorToggle 0 3 ≠ selectorToggle 0 1 := by
  decide

end C169Carry
end LonelyRunner

#print axioms LonelyRunner.C169Carry.positiveExitCount_eq_kappa
#print axioms LonelyRunner.C169Carry.c169Chart_not_product_additive
#print axioms LonelyRunner.C169Carry.kappa_cocycle
#print axioms LonelyRunner.C169Carry.lift_comp
#print axioms LonelyRunner.C169Carry.twistedExp_comp
#print axioms LonelyRunner.C169Carry.twisted_seam_defect
#print axioms LonelyRunner.C169Carry.closed_lifted_endpoint
#print axioms LonelyRunner.C169Carry.closed_lifted_loop_flat
#print axioms LonelyRunner.C169Carry.selector_xor_law
#print axioms LonelyRunner.C169Carry.positiveSelector_unique
#print axioms LonelyRunner.C169Carry.same_selector_signature_different_carry
#print axioms LonelyRunner.C169Carry.same_carry_different_selector_toggle
