---
id: HYP-3243
title: LRC14 topology geometry and graph proof-route atlas
status: SYNTHESIS / executable proof-carrier tournament; not an LRC14 proof
source: codex-2026-06-28
tangent: T1340
technique: LTI-340
tournament_technique: LTT-240
script: 04-computation/lrc14_topology_geometry_graph_routes_codex_20260628.py
result: 05-knowledge/results/lrc14_topology_geometry_graph_routes_codex_20260628.out
reflection: 07-reflections/lrc14-topology-geometry-graph-proof-routes-codex-20260628.md
related:
  - HYP-3242
  - HYP-3241
  - HYP-3240
  - HYP-3238
  - HYP-3237
  - HYP-3236
  - HYP-3235
  - HYP-3234
  - HYP-3233
  - HYP-3232
  - HYP-3230
  - HYP-3228
  - HYP-3227
  - HYP-3225
  - HYP-3224
  - HYP-3223
  - HYP-3222
  - HYP-3220
  - HYP-3201
  - HYP-3128
  - HYP-3123
  - HYP-3108
  - THM-572
  - OPEN-Q-108
---

# HYP-3243: LRC14 Topology Geometry And Graph Proof-Route Atlas

## Claim

The useful visual route to LRC14 is not one picture.  It is a typed atlas of
compatible pictures, each preserving a different part of the LRC predicate:

```text
circle arrangement:
  danger arcs, endpoints, wall crossings, safe topes

oriented matroid:
  topes for strict lonely intervals, cocircuits for endpoint equality atoms

Cech / component topology:
  positive safe measure, closed-boundary atoms, barcode persistence

D7 / Borsuk-Ulam:
  antipodal Phi_14 witness pairs, index 3, parity/sign sidecar

cyclotomic witness strata:
  base Phi_14 core, dilation-promoted Phi_{14d} core

graph energy:
  Green conductance, Fiedler cuts, effective resistance, normal-fan faces
```

The current proof frontier after HYP-3240/HYP-3241, and now mac-mini's
HYP-3242/S78 Euler-characteristic nerve packet, is therefore sharper: the
witness construction side is mostly organized by `Phi_14` plus `Phi_{14d}`,
and the cover topology side has a canonical invariant
`cap = chi_meas(danger-cover nerve)`.  The remaining hard work is the
structural half:

```text
tight-locus finiteness + bulk equidistribution + legal gluing of the core/bulk wall.
```

## Why This Helps

The topology of the circle arrangement gives a visual criterion for what the
theorem needs:

- an open safe tope proves strict loneliness;
- a closed cocircuit at the boundary is an equality atom;
- AP and Goddyn-Wong share the same base `Phi_14` cocircuit witnesses;
- dilations move the same core idea to `Phi_{14d}`;
- a true bad row would have to kill every open tope and evade the named
  cocircuit/witness packets.

This is a proof-route statement rather than a metaphor.  Each picture has a
legal quotient test:

```text
carrier is proof-grade
iff it preserves the LRC predicate,
or every destroyed coordinate is reconstructible, dual-annihilated,
or named as sidecar debt.
```

This is HYP-3201's compression rule applied to topology and geometry.

## Route Atlas

**Visual topology route.**  Cut `R/Z` by all danger-arc endpoints
`(14k +/- 1)/(14v)`.  The cells are topes.  A positive-length all-safe tope is
a witness.  If no positive tope exists, inspect closed boundary cocircuits:
the known equality atoms route to the shared AP/Goddyn-Wong `Phi_14` witness
packet or to dilation-promoted `Phi_{14d}` packets.

**Cech route.**  Treat safe arcs and their overlaps as a finite good-cover
problem.  Positive `H0`/component mass is the bulk certificate; closed
boundary-only behavior must be separated from open components.  This route is
useful for visual proof pages because it can show where an attempted covering
really has a gap, but it must keep endpoint owners.

HYP-3242/S78 adds the clean formula for this route: the cap/lonely measure is
the measure-weighted Euler characteristic of the danger-cover nerve.  In this
atlas that identity is the Cech carrier's invariant, while the present HYP
adds the route ledger that says when the invariant must be supplemented by
owners, signs, witness strata, or state-lift sidecars.

**Oriented-matroid route.**  The endpoint arrangement has topes and
cocircuits.  A counterexample must define a no-open-tope/no-known-cocircuit
packet.  This suggests a finite chamber theorem: every primitive chamber is
open-safe, AP/GW equality, dilation equality, state-lift forbidden, or named
residual debt.

**Graph-energy route.**  HYP-3236's positive covariance conductance graph
visualizes the even/positive face: AP maximizes algebraic connectivity and
minimizes Green resistance in the bounded scout.  But the graph clips
negative covariance, so the usable proof packet is:

```text
Green slack + normal-fan/Toeplitz face + odd-negative sidecar.
```

This can discharge finite traps but cannot replace the core witness packet.

**Root-motion and ear route.**  Lee-Yang/PGF roots and ear decompositions
turn adding a runner into root movement plus retained payload.  They inspire a
structural proof in which root collisions identify residual chambers, while
directed/odd ears retain the parity data needed by HYP-3238.

**State-lift route.**  If the finite chamber atlas leaves a primitive bad
atom, the desired output is a typed state lift into the forbidden `H=7`
endpoint (THM-572).  This is the structural "last picture": a bad geometric
cell becomes a named tournament contradiction rather than a raw scalar
counterexample.

## Assumption Challenge

For this exploratory pass the tournament vertices are not runners.  The
candidate vertex sets considered were:

```text
proof carriers, topes, Cech components, endpoint owners, wall crossings,
Phi witness strata, conductance cuts, normal fan chambers, PGF roots,
ear payloads, and state-lift obligations.
```

The preserved predicate is the LRC14 packet:

```text
strict open witness OR known equality core OR finite forbidden residual.
```

The information destroyed by raw runner/scalar quotients is exactly the
payload the visual proof needs: endpoint owners, parity sign, cocircuit role,
bulk/core status, and witness-construction address.

## Tournament Analysis

Pairwise observable: for two visual/structural carriers, which one preserves
more of the current LRC proof packet while destroying less of the sidecar data
needed by the other.

Switch/gauge: orient `A -> B` when `A` keeps a higher weighted set of current
obligations, penalized by destroyed core/sign/boundary data.  Ties prefer
fewer destroyed coordinates and then the route that is closer to a finite
formal atlas.

The executable scout ranks `12` proof carriers.  Its Hamiltonian path is:

```text
oriented_matroid_topes_cocircuits
-> circle_endpoint_arrangement
-> cech_nerve_safe_components
-> d7_borsuk_ulam_index_packet
-> cyclotomic_witness_strata
-> normal_fan_toeplitz_fejer
-> lee_yang_pgf_root_motion
-> ear_payload_graph
-> state_lift_forbidden_H7
-> green_conductance_laplacian
-> bulk_spec_equidistribution
-> raw_scalar_shadow
```

This ranking should not be read as "topes win everything."  It says the
oriented cell/cocircuit language is currently the least lossy common visual
carrier because it can hold boundary equality, finiteness, and compression
legality at once.  The actual proof still needs the witness-core and
Cech/bulk routes to carry the complementary payload.

## Next Pull

Turn the atlas into a finite theorem schema:

```text
primitive row
-> open safe tope
   or AP/GW Phi_14 equality
   or dilation Phi_{14d} equality
   or finite chamber discharge by Toeplitz/Green/root-motion
   or state-lift H=7 contradiction
   or named residual debt.
```

The most promising visual proof page is a layered diagram: circle endpoints
on the bottom, topes/cocircuits as the cell complex, Green/Toeplitz faces as
energy labels on cells, and the `D7` antipodal index as the sign sidecar on
the core boundary.
