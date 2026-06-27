---
id: HYP-2995
title: LRC14 cocycle carrier atlas and packet-cocycle theorem target
status: PROOF-INTERFACE / cocycle unification atlas and theorem target; not a proof
source: codex-2026-06-24-S167
related:
  - HYP-2991
  - HYP-2990
  - HYP-2989
  - HYP-2988
  - HYP-2987
  - HYP-2986
  - HYP-2985
  - HYP-2984
  - HYP-2979
  - HYP-2978
  - HYP-2976
  - HYP-2974
  - HYP-2970
  - HYP-2969
  - HYP-2963
  - HYP-2230
  - HYP-2241
  - HYP-2618
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_cocycle_carrier_atlas_codex_s167.py
  - 05-knowledge/results/lrc14_cocycle_carrier_atlas_codex_s167.out
---

# LRC14 Cocycle Carrier Atlas

## Claim

A useful LRC14 cocycle is a signed local obstruction to forgetting a coordinate.
It lives on a labelled packet fiber, closes when local certificates glue, and
becomes harmless only by one of the allowed exits:

1. it is exact as a coboundary or potential;
2. it is annihilated by a dual or orthogonality certificate;
3. it is a known AP/Goddyn-Wong boundary atom;
4. it descends to a smaller labelled packet family; or
5. it is emitted as a named residual sector, currently F7/THM-572 state-lift
   debt.

This is a unification of HYP-2990's no-free-slider rule and HYP-2991's local
Haar zipper cocycle.  The atlas does not prove LRC14.  It names the cochain
that every quotient must either preserve or discharge.

## Carrier Dictionary

The S167 script records sixteen cocycle carriers:

| Carrier | Local cochain | LRC14 role |
|---|---|---|
| `labelled_packet_master_cocycle` | total lost-coordinate obstruction `omega_P` | theorem target for all packet quotients |
| `haar_zipper_zeta` | `zeta=T00-T01-T10+T11` | local mixed Haar/fixed-margin switch that margins erase |
| `endpoint_credit_winding` | endpoint credit transition `K(a,b)` | positive winding is the danger-cover obstruction; potentials certify safety |
| `carry_crt_cocycle` | `k mod 14` carry through the `Res_27` lift | links parity, apex divisibility, and pair sums |
| `owner_deletion_derivative` | owner obligations uncovered by deletion | repairs visible quotient leaks through private-owner derivatives |
| `ramanujan_exact_period_trace` | `c_q(a-b)` primitive exact-period trace | repairs squarefree/divisor quotient blindness |
| `fejer_toeplitz_moment` | Fejer quadratic form `v* M v` | dual certificate for positive-open packets |
| `boundary_moment_multichart` | chart-to-chart moment defect | blocks overreading one all-covered denominator chart |
| `product_rule_tiling` | product class: zero, atom, owner strip, handoff, descent, residual | keeps Haar/Walsh tiling products address-retained |
| `farey_k33_determinant` | excess `e=14p-q` and `K_{p,q}` incidence depth | routes unit-excess rows to AP/GW, C27, or K33/state lift |
| `c27_unital_transfer` | local q=3 pair-completion transfer defect | tracks petal/GW/K33 shell transfer |
| `root_packet_boundary` | root sum as endpoint boundary | clean cochain language for closed packets and OCF/tournament analogies |
| `metagraph_transfer_cocycle` | successor defect / transfer curvature | makes recursion paths telescope or expose SCC debt |
| `sequence_shadow_difference` | finite difference under `n+2` / `n*2` mode changes | separates useful recurrence shadow from analogy |
| `octahedral_hodge_current` | discrete curl/divergence defect | support-six current route and Hodge split |
| `ocf_activity_coimage` | hard-core activity residue | coimage/path-homology support repair |

## Tournament Analysis

Vertices are cocycle carriers, not runners.  The proof gauge orients toward the
carrier that better preserves the LRC predicate, packet fiber, local defect,
closedness/gluing law, legal exits, formal checkability, and LRC14 relevance.

Proof-gauge fingerprint:

```text
vertices=16
score_hist={0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 6: 2, 7: 1, 8: 2, 10: 1, 11: 1, 12: 1, 14: 3}
directed_3cycles=4
SCC_sizes=[1, 1, 1, 1, 1, 5, 1, 1, 1, 3]
Hamiltonian_path_count=27
```

Canonical proof-gauge path:

```text
labelled_packet_master_cocycle
> haar_zipper_zeta
> endpoint_credit_winding
> product_rule_tiling
> ramanujan_exact_period_trace
> carry_crt_cocycle
> fejer_toeplitz_moment
> owner_deletion_derivative
> boundary_moment_multichart
> farey_k33_determinant
> c27_unital_transfer
> root_packet_boundary
> metagraph_transfer_cocycle
> octahedral_hodge_current
> ocf_activity_coimage
> sequence_shadow_difference
```

Exploration gauge puts cross-domain reach higher and becomes transitive with
one Hamiltonian path.  The `12` edge flips between the gauges are useful: they
mark carriers that are fertile analogies but not yet theorem-safe.

## Packet-Cocycle Theorem Target

Let `P(S)` retain qdiv, exact `M=p/q` / Farey data, Haar/Baire topology,
endpoint owners, exact-period labels, Fejer certificate data, C27/K33/state
labels, and smoothing route.  For every quotient `Q:P(S)->Y`, construct a
fiber cochain `omega_Q` measuring the coordinate that `Q` forgets.

Then `Q` is theorem-safe only if the LRC predicate is constant on `Y`-fibers,
or `omega_Q` is reconstructed, exact as a coboundary/potential, annihilated by
a dual certificate, descended to a smaller packet family, identified as
AP/GW boundary equality, or emitted as a named F7/THM-572 residual.

The sharp missing piece is to make this statement executable on the HYP-2963
packet bank and then familywise: packet key, quotient map, coefficient ring,
fiber graph, cochain values, closedness law, and one legal exit per packet
fiber.

## F7 Interpretation

F7 should not mean "unknown row".  In the cocycle language it should mean:

```text
a non-exact, non-dual-annihilated, non-descended cocycle class left after the
known Haar, endpoint, Ramanujan, Fejer, Farey/C27/K33, carry-owner, and
boundary-moment teeth all fail.
```

To make this rigorous, the next pass must specify the coefficient ring, the
state-lift map into THM-572-style tournament obstruction data, and the packet
family on which the class is measured.

## Assumption Challenge

The session explicitly rejects the default assumption that tournament vertices
should be runners or arcs.  For cocycle work the better vertex sets are:
packet fibers, proof carriers, endpoint arrows, Haar squares, carry
derivatives, exact-period modes, transfer states, product-rule classes, and
residual proof obligations.  The preserved predicate is not a scalar score but
the ability to glue local certificates into a global LRC14 certificate.

## Next Work

1. Define an `omega_Q` record schema for every quotient used in HYP-2963.
2. Compute `zeta` and product-rule classes on actual named packet families.
3. Turn HYP-2970 endpoint-credit potentials into exact interval certificates.
4. Prove the carry-owner no-leak theorem: local/global carry cochains are
   reconstructed from paired owner derivatives or pay strict loneliness tax.
5. Define F7 as a cocycle class with a coefficient ring and a THM-572 state
   lift, not as an untyped failure bucket.
