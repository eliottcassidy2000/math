# Minkowski, Circuit, Ising, And De Moivre Carrier Atlas

codex-2026-06-27-S264 for HYP-3111.

## What Was Merged

This session merged four outside prompts into the current LRC14 proof frontier:

- Minkowski theorem as a symmetric convex-body threshold for a generated
  q-lattice of finite packet distributions.
- Circuit complexity as a size/depth/uniformity ledger for proof obligations.
- The Ising model as a finite partition-polynomial zero carrier tied to
  Lee-Yang whole-zero geometry.
- De Moivre's quintic as an exact algebraic fold:
  `x = u - a/u` turns the quintic into the auxiliary quadratic
  `y^2 + b*y - a^5 = 0` with `y=u^5`.

The important control choice is that the tournament vertices are not runners,
arcs, roots, or abstract graph names.  They are proof carriers:
finite-address packet, observer-gluing certificate, proof-circuit ledger,
zero-real ear map, Lee-Yang root curve, De Moivre fold, Minkowski q-body,
Ising partition zero, and raw scalar `p0`.

## Strongest Signals

The Minkowski sidecar made the Bravais/q-rank idea sharper.  The named packet
q-vectors span affine rank `6`, and the generated lattice has covolume proxy
`6.795578624e-12`.  The Euclidean ball threshold
`vol(B_R) = 2^rank*covolume` occurs at radius `0.020934`, just below the
shortest nonzero named q-difference `E_star_k12 <-> gw12` with norm
`0.021715`.  This is a real forcing threshold for a declared q-body, but not a
witness-time theorem.

The circuit sidecar is useful because it refuses to hide the route structure.
With ten input obligations, the scout circuit has gate-size `8`, depth `4`,
and five minimal monotone certificate families:

```text
q_witness_gate
level7_sieve AND dyadic_lift
demoivre_fold AND observer_gluing_certificate AND finite_address_packet
pgf_root_curve AND zero_real_ear_map AND minkowski_lattice_body
  AND observer_gluing_certificate AND finite_address_packet
pgf_root_curve AND zero_real_ear_map AND ising_partition_zero
  AND observer_gluing_certificate AND finite_address_packet
```

The De Moivre fold is the cleanest formalizable algebraic artifact.  The exact
Laurent identity

```text
(u - a/u)^5 + 5*a*(u - a/u)^3 + 5*a^2*(u - a/u)
  = u^5 - a^5/u^5
```

has no numeric debt.  Its weakness is not algebraic; it is that it preserves a
branch/fold normal form rather than an LRC witness predicate unless it is tied
to finite-address data.

The Ising packets confirmed the Lee-Yang moral in miniature: finite
ferromagnetic proof-carrier graphs have partition zeros numerically on the
unit circle.  The useful signal is again the whole zero set, not a single
moment or nearest-root scalar.

## What This Does Not Prove

Minkowski theorem does not prove LRC14 here.  It says that once a lattice and a
symmetric convex body are declared, a volume threshold forces a nonzero lattice
point.  The missing theorem is the legal declaration: the q-body must encode a
specific LRC predicate and its nonzero lattice point must map to a valid
witness, lift, or observer certificate.

Circuit complexity also does not prove lower bounds here.  The scout only uses
the stable complexity vocabulary: gate count, depth, essential inputs, and
minimal certificates.  The lower-bound use would be to show that any quotient
which removes finite-address or observer-gluing fields cannot compute the
residual predicate at the claimed depth.

The Ising model does not replace HYP-3109.  It reinforces HYP-3109's warning:
single-value scalarizations erase root geometry.  If Ising enters the proof,
it should enter through a finite partition-zero object and a declared mapping
from spins/interactions to LRC proof obligations.

De Moivre's quintic does not naturally encode the 14-runner row.  Its role is
as a reusable finite-depth cancellation template, especially for branch-orbit
or fifth-root fold bookkeeping.

## Current Frontier

The improved argument shape is:

```text
primitive residual row
  -> CRT/dyadic/level-7 finite address
  -> observer-gluing certificate
  -> retained sidecar:
       Lee-Yang root/zero-real ear geometry, or
       De Moivre-style exact fold plus finite branch address, or
       Minkowski q-body whose forced point is legally tied to an LRC predicate,
  -> LRC14Statement
```

The tournament readout puts `finite_address_packet` and
`observer_gluing_certificate` first.  That matches the older polynomial-method
bridge: the `14=2*7` descent still demands the level-7 and dyadic address
work.  The new carriers should be used as pressure gauges on that bridge, not
as replacements for it.

Post-rebase S262 added a root-lattice reachability supplement to HYP-3108, and
it changes how the Minkowski part should be read.  In the anchored `{0}+7`
bank, high `p0` is reciprocal-flat rather than Bragg-crystalline:
`corr(p0,Bravais_peak)=-0.430` and
`corr(p0,residue_entropy)=+0.541`, while
`corr(p0,nearest_root)=+0.899` and
`corr(p0,real_roots)=-0.483`.  The q-body should therefore be designed around
root stratum, entropy/flatness, segment clearance, and finite-address status.
A convex body around a large lattice peak would point the wrong way.

## S265 Bridge Update

After rebasing over HYP-3113, the clean way to use this lane is not "four
analogies that might prove LRC14."  It is a cut-payload ledger feeding the
two-map frontier:

```text
Minkowski q-body -> memory-lattice-ear map
Circuit DAG      -> packet-sheaf legal-exit regression test
Ising zeros      -> root-curve / Lee-Yang zero-arc portfolio
De Moivre fold   -> finite branch-address / first-obstruction fold test
```

The regenerated script adds a duodecimal audit: four carriers times three
legal cells, namely preserved predicate, destroyed coordinate, and handoff
payload.  All twelve cells are now named, but no carrier is proof-closed.

The resulting guardrail is sharper than the original atlas.  A future proof
shortcut cannot merely cite a Minkowski volume threshold, a small circuit
depth, a unit-circle Ising zero set, or the De Moivre quadratic fold.  It must
also provide the surviving packet field:

```text
q_body_inequality_word
proof_circuit_missing_input_vector
ising_zero_arc_signature
demoivre_branch_orbit_word
```

Those are the fields that can plausibly enter HYP-3113's packet-sheaf legal
exit.  Without them, the outside theorem is only a diagnostic pressure gauge.

## Next Programmatic Tests

1. Replace the toy Ising carrier graph with the actual HYP-3109 zero-real
   one-swap component and compare partition-zero angles before and after
   crossing the `#real=0 -> #real=2` wall.
2. Search for a q-body whose defining inequalities are LRC-native:
   component count, segment clearance, reciprocal-flat residue entropy,
   endpoint-owner debt, and finite-address status rather than arbitrary
   Euclidean radius or Bragg-peak concentration.
3. Turn the De Moivre Laurent identity into a small formal lemma, then test
   whether any HYP-3110 Jacobi/crystallographic branch packet needs exactly a
   fifth-root orbit sidecar.
4. Use the circuit ledger as a regression check: any proposed proof shortcut
   must either supply `q_witness_gate`, supply both `level7_sieve` and
   `dyadic_lift`, or route through `observer_gluing_certificate` and
   `finite_address_packet`.
