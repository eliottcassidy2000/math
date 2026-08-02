# The modular free-product frame is a sidecar grammar, not yet an action

Date: 2026-08-02
Status: **STRUCTURAL SYNTHESIS / REFLECTION, CORRECTED 2026-08-02**, not a
theorem and not a proved dependency.  The correction below replaces the
earlier overstatement that cyclic substitution itself supplied a `C_3`
action.

## Inheritance pass

The proposed organizing picture is the free-product presentation
`PSL_2(Z) ~= C_2 * C_3`: the binary/parity tower and the ternary substitution
tower should be two operations on one state rather than two unrelated
analogies.  Current canon supplies exact pieces of both faces, but also proves
that their scalar shadows are insufficient.

On the quartic side, THM-2950 identifies the four half-system traces as a
`V_4` torsor with affine symmetry

```text
V_4 semidirect S_3 ~= S_4,
S_4/V_4 ~= S_3,                                         (1)
```

and proves equality of the quartic and cubic-resolvent discriminants.  Thus a
cubic resolvent retains the `S_3` quotient while forgetting the point of the
`V_4` lift torsor.  THM-2951 sharpens the loss: although the full fifth
compound reconstructs multiplication, there is no nonzero signed-pair
equivariant linear map from that compound to the balanced third-compound
sector.  A selected cofactor or norm therefore cannot manufacture the missing
quartic phase.

On the ternary side, THM-3121 gives the exact `C_3` path-cover substitution
kernel, and THM-3134 identifies its sufficient state as the complete endpoint
Taylor jet of the centered permutation polynomial.  The Hamiltonian-path
count is only the endpoint value.  Two tournaments with that same value but
different first jets already give different `C_3` substitutions.

The closest hostile examples are consequently parallel:

```text
quartic face: full norm / one cofactor, but different V_4 trace frames;
ternary face: same endpoint value, but different path-cover jets.          (2)
```

The least-used relevant sidecar is THM-2056's Farey fan: it gives a cheap
arithmetic test for whether proposed binary and ternary moves really act by a
common unimodular slope grammar rather than merely having the right arities.

## 1. What the modular-group proposal can mean exactly

The whole `V_4` torsor is **not** the order-two free factor.  It is a
three-direction kernel of pair flips.  A candidate order-two generator must
be a *pointed involution* on that torsor: for example, one chosen pair
reflection together with the orientation datum that tells which lift was
used.  Forgetting the point recovers only the cubic-resolvent quotient.

The order-three generator also needs a pointer.  Cyclic rotation acts on the
**labelled input triple** of full endpoint jets; for equal or canonically
identified factors it also becomes the genuine cyclic-wreath automorphism of
the substituted tournament.  The numerical THM-3121/3134 profile transform
is invariant under this rotation, so its output forgets the designated cyclic
slot origin even when the complete output jet is retained.  Thus the smallest
currently honest common state has the schematic form

```text
(pointed V_4 half-system, pair orientation;
 rooted triple of complete path-cover endpoint jets, quotient walk kernel). (3)
```

A genuine `C_2*C_3` action requires maps `S,R` on one carrier satisfying
`S^2=1` and `R^3=1`.  No commutation or braid relation should be inserted.
On the intrinsic four-point half-system torsor, however, after choosing an
origin, one nonzero pair flip, and an oriented three-cycle (none is canonical
on THM-2950's unbased torsor), the chosen generators already satisfy the extra relation
`(SR)^3=1`; they generate `V_4 semidirect C_3 ~= A_4 ~= PSL_2(F_3)`.
Therefore that torsor is the level-three congruence shadow, not a faithful
modular carrier.  More importantly, both maps must preserve the same physical
or algebraic carrier.  Current canon gives the two operations on different
objects; it does not yet give this common carrier map.

## 2. Source, target, preserved data, and destroyed data

For the ternary face the source is a **rooted labelled triple** of tournament
path-cover profiles, the target is the substituted profile, and the map is
THM-3134's Gregory--Newton transform.  It computes every ordered path-cover
count and, in the equal-factor case, the cyclic-wreath divisibility, but its
numerical readout is invariant under cyclic rotation of the input labels.
Projection to the endpoint value additionally
destroys the higher run counts.  The endpoint jet restores the run counts;
the rooted block coordinate separately restores the cyclic phase.

For the quartic face the source is multiplication on a reduced three-pair
real algebra, the target is the four half-system traces (or their cubic
resolvent quotient), and the map is the balanced exterior-power construction
of THM-2950/2951.  The full fifth compound preserves the operator, but scalar
norm/cofactor readouts destroy pair origin.  A pointed pair-respecting
bivector, or equivalently the full pair decomposition plus one lift point, is
the missing sidecar.

This gives a useful negative rule:

> Scalarize neither free-factor move before composing them.  The ternary move
> needs both a jet and a rooted cyclic block; the pointed binary move needs a
> torsor origin.

The two losses are orthogonal.  Restoring only one cannot recover the other.

## 3. Consequence for the degree-four Keller/resolvent route

The identity `S_4/V_4 ~= S_3` explains why grade-three anatomy can constrain
a hypothetical graph quartic's resolvent cubic.  It does **not** reconstruct
the quartic cover.  Even if the depressed cubic law, discriminant law,
Jelonek lead, and cuspidal geometry determine the three labelled resolvent
branches, the remaining fibre is a pointed `V_4` torsor.  THM-2951 rules out
recovering that point by a natural linear contraction of the fifth-compound
certificate.

Accordingly the next useful quartic test is not another scalar discriminant
identity.  It is a branchwise `V_4`-covariant observable whose value is
nonzero at a common physical atom and whose transport is compatible with the
resolved `S_3` monodromy.  Failure to supply that observable is exactly the
lift obstruction, not evidence against the cubic anatomy itself.

## 4. Cheapest decisive experiments

1. **Farey-flank test.**  Reduction
   `P^1(Q) -> P^1(F_3)` preserves the one-way implication “Farey neighbors
   remain distinct”, but not its converse and not THM-2056's fixed sufficient
   polar certificate.  This does not show that the transported slopes are
   physically unsafe.  A certificate-level lift must therefore transport the
   exact determinant, polar polygon, and quadratic metric, not merely the four
   residue points.
2. **Reduced-word fingerprint.**  The cheapest test is already decisive:
   `(SR)^3` is the identity on the four-point torsor, whereas in
   `PSL_2(Z)` it is `T^3` and is nontrivial.  For example `T^3(0,1)=(3,1)`;
   the two slopes become equal only after reduction modulo `3`.  The missing
   sidecar is consequently a `Gamma(3)`/Farey lift, not another search for a
   short torsor word.
3. **Two-sidecar hostile.**  Hold the endpoint jet fixed and vary the pointed
   `V_4` origin, then hold the origin fixed and vary a higher jet at constant
   endpoint value.  Any proposed common scalar invariant must distinguish
   both pairs before it can transport the combined state.
4. **Quartic physical atom.**  Search the existing Keller charts for one
   common atom carrying a nonzero `V_4`-covariant value.  The cubic resolvent
   supplies the quotient labels; this atom would supply the missing lift.

## 5. Strongest current survivor

The binary and ternary towers are not yet proved to be a faithful modular
action on one object.  What is proved is more precise and already useful:

```text
ternary composition is functorial after retaining the full endpoint jet,
but cyclic input phase additionally needs a rooted block;
quartic/cubic passage is functorial after retaining a pointed V_4 lift;
both scalar contractions have exact minimal hostile examples.             (4)
```

This reframes the modular proposal as a **sidecar-completion problem**.  The
intrinsic exact object is the congruence quotient
`PSL_2(Z) -> PSL_2(F_3) ~= A_4`; its first lost coordinate is the nontrivial
`Gamma(3)` action of `(SR)^3=T^3`.  The next theorem should construct one
carrier supporting that lift together with the rooted endpoint-jet move, or
prove that their physical predicates cannot coexist.  Until then,
`C_2*C_3` is a disciplined search grammar rather than a proved symmetry.
