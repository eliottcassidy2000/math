---
id: HYP-3061
title: LRC14 geometry-regime archive audit
status: SYNTHESIS / archive consolidation and guardrail; not a proof
source: codex-2026-06-26-S225
tangent: T1143
related:
  - HYP-3058
  - HYP-3057
  - HYP-3056
  - HYP-3055
  - HYP-3054
  - HYP-3047
  - HYP-3043
  - HYP-3039
  - HYP-3003
  - HYP-2963
  - HYP-2943
  - HYP-2934
  - HYP-2928
  - HYP-2900
  - HYP-2887
  - THM-572
  - OPEN-Q-108
---

# HYP-3061: LRC14 Geometry-Regime Archive Audit

## Claim

The old project motif that `5,6,7` marks a geometry transition is real, but
only after the axis is declared.

On the triangle-group and regular-tiling axes, the clean statement is

```text
(2,3,5)  spherical / Platonic
(2,3,6)  Euclidean / flat / planar wall
(2,3,7)  hyperbolic / first negative-curvature debt
```

and equivalently for `{3,q}`:

```text
q=5  spherical icosahedral triangle vertex
q=6  Euclidean triangular tiling wall
q=7  first hyperbolic triangular tiling
```

On the tournament-size axis, the nearby but different statement is

```text
n=5  last Platonic / icosahedral-looking boundary
n=6  first pivot: beta_3, rooted-perspective failure, DT-sufficiency failure
n=7  apex-prime / forbidden H=7 / seven-sector obstruction layer
```

The useful LRC14 import is a typed sidecar:

```text
geometry_regime_signature:
  axis
  input
  regime
  curvature_or_defect
  preserved_payload
  destroyed_payload
  lrc_handoff
  source_artifacts
```

This turns old geometry analogies into controlled-forgetting objects.  A
geometry slogan is admissible only after it says what LRC predicate it
preserves, what coordinate it forgets, and which packet route or certificate
will discharge the loss.

## Archive Axes

| Axis | `5` side | `6` side | `7` side | LRC14 reading |
|---|---|---|---|---|
| Schwarz triangle `(2,3,p)` | `(2,3,5)`, reciprocal sum `31/30`, spherical icosahedral/dodecahedral world | `(2,3,6)`, reciprocal sum `1`, Euclidean wall | `(2,3,7)`, reciprocal sum `41/42`, hyperbolic defect `1/42` | Use HYP-3058's reciprocal-curvature sidecar before flattening a three-lane packet. |
| Schlaefli `{3,q}` | `{3,5}` icosahedral positive defect | `{3,6}` flat triangular tiling | `{3,7}` first hyperbolic triangular tiling | `q=6` is the honest planar wall; `q=7` is first negative curvature. |
| Regular maps `{p,q}` | Platonic finite region `(p-2)(q-2)<4` | Euclidean tiling wall `(p-2)(q-2)=4` | Hyperbolic maps `(p-2)(q-2)>4` | HYP-2943: use regular-map curvature only with a concrete recursion carrier. |
| Tournament size `n` | `n=5`: five-cycles, Platonic boundary, `G_5` f-vector coincidence | `n=6`: first beta_3, first A000568/rooted-perspective defect, two 3-cycles can coexist | `n=7`: Paley/seven-sector/apex prime, H=7 obstruction, mod-7 barrier lifts | Do not confuse size `n` with triangle parameter `p`; attach value origin and observer-cut sidecars. |
| Fugacity / gap function | `g(phi)=-1`, spherical side | `g(tau)=0`, Euclidean/tribonacci wall | `g(2)=+1`, tournament/hyperbolic side | The axis is a model of pressure, not a proof unless packet payload is retained. |
| LRC Farey product lane | `p=1` unit-excess star / q-threshold parent | `p=2`, `2/27` planar two-block strip | `p>=3`, `3/41` first `K_{3,3}` product wall | S131: product incidence becomes serious at `p=3`; retain exact `M`, `p*q`, endpoint owners. |
| Support-six current lane | simplex/residue skeleton before the six-packet lift | six residues / octahedral current carrier | octahedron cycle rank `7` and apex-seven state-lift pressure | HYP-2887: useful only if divergence, face curl, and owner payload are explicit. |
| Annular prism lane | finite solid/Archimedean local word | prism/antiprism annulus | `n=14` gives `V=28`, curvature `1/14` | HYP-2943: possible AP/GW/LRC14 carrier, not yet a theorem. |
| Totient/recursion curvature | small prime content | exact-period floor | `14=2*7` and residual prime curvature at n=14 | HYP-2900: finite denominator bases are not enough; analytic equidistribution remains debt. |

## Source Threads

The archive contains several coherent but differently typed strands.

- `07-reflections/five-as-bridge.md` gives the core triangle-group spine:
  `(2,3,5)` spherical, `(2,3,6)` flat, `(2,3,7)` hyperbolic, and ties `n=5`
  to Platonic-looking tournament coincidences.
- `07-reflections/six-as-pivot.md` makes `n=6` the first high-dimensional
  obstruction: beta_3 appears, two disjoint 3-cycles can coexist, and small
  directed-triangle sufficiency fails.
- `07-reflections/the-five-platonic-tournaments.md` records the five Platonic
  analogies but already warns that the tournament graph is not literally the
  icosahedron.
- `07-reflections/the-tessellation.md`,
  `07-reflections/the-modular-tournament.md`, and
  `07-reflections/the-three-infinity-structure.md` supply the
  `{3,infinity}`, modular, and gap-function language.
- `05-knowledge/results/pivot_six_s90w.out`,
  `05-knowledge/results/platonic_primes_s90bm.out`,
  `05-knowledge/results/two_three_seven_s90az.out`, and
  `05-knowledge/results/coxeter_hyperbolic_s120.out` preserve older
  Coxeter/Schwarz/Platonic computations.
- HYP-2943 imports the regular-map and tiling wall into LRC14 recursion
  carriers, especially the `n=14` annular prism/antiprism possibility.
- HYP-2887 supplies the support-six/octahedral current carrier, where the
  octahedron has cycle rank `7` and no harmonic one-current once divergence
  and face curl are controlled.
- HYP-2934 and S131 identify the Farey product lane: `1/13` is the star
  parent, `2/27` is the planar two-block strip, and `3/41` is the first
  `K_{3,3}` product wall.
- HYP-2900 gives the totient-curvature side of `14=2*7`, with an
  even-half residual prime curvature signal but an analytic-equidistribution
  guardrail.
- HYP-3058 records the Fermat-Catalan reciprocal-sum threshold as the same
  hyperbolic sidecar.

## Guardrails

The archive also contains traps that should remain attached to the motif.

1. `5,6,7` is not one coordinate system.  The triangle axis is
   `(2,3,5)`, `(2,3,6)`, `(2,3,7)`, while the tournament axis is
   `n=5`, `n=6`, `n=7`.
2. The `G_5` f-vector matching the icosahedron is not an `A_5` or Galois
   icosahedron theorem.  The f-vector coincidence is evidence for a
   boundary motif, not an isomorphism claim.
3. `(2,3,7)`, `42`, and `1/42` are route pressure, not proof certificates.
   They must be coupled to exact `M`, endpoint-owner, topology, route, and
   certificate payload.
4. Paley `T_7` automorphism counts must distinguish tournament automorphisms
   from full symmetry including anti-automorphism/converse symmetry.
5. The seven-point hexagonal flower has only the center/petal graph edges; it
   is not the complete tournament on seven vertices.  It is a flat carrier
   metaphor unless cross-petal orientation is supplied.
6. HYP-2943's hexagonal norm `7` and centered-hex ring counts are different
   sequences.  Use the norm only when the corresponding state-lift packet has
   been built.
7. The old `1729` modular thread was later weakened by HYP-2306.  Do not use
   it as a modular law.

## Packet Sidecar

Recommended packet vocabulary:

```text
geometry_regime_signature:
  axis: schwarz_triangle | schlaefli_tiling | regular_map |
        tournament_size | fugacity_gap | farey_product |
        support_six_current | annular_prism | totient_curvature |
        arithmetic_motif
  input: string
  regime: spherical | euclidean | hyperbolic | pivot | forbidden | motif
  curvature_or_defect: rational/string
  preserved_payload: list
  destroyed_payload: list
  lrc_handoff: q-witness | AP/GW boundary | C27 petal |
               K33 state lift | Fejer/Toeplitz | F7 debt |
               archive motif only
  source_artifacts: list
  guardrail: string
```

Mandatory companion fields when this sidecar is used on LRC14 packets:

```text
exact_M
endpoint_owner_payload
safe_topology_payload
value_origin_type
observer_cut_payload
magnitude_spectrum_required
route_label
certificate_or_state_lift_payload
destroyed_coordinate
```

## Tournament Analysis

Vertices are geometry carriers and proof obligations, not runners:

```text
labelled_packet_sheaf
geometry_regime_signature
observer_cut_payload
exact_M_Farey_spectrum
annular_14_prism_payload
octahedral_current_carrier
K33_product_wall
hyperbolic_reciprocal_signature
totient_euler_curvature
raw_5_6_7_numerology
```

Pairwise observable: compare the retained LRC predicate tuple

```text
(exact M, binding scale, endpoint owner, topology, route handoff,
 magnitude spectrum, certificate/state-lift payload).
```

Gauge: orient `A -> B` when `A` preserves all fields preserved by `B` and at
least one additional field needed to separate a route/status-changing quotient
fiber.  Ties are broken by this conservative retention path:

```text
labelled_packet_sheaf >
geometry_regime_signature >
observer_cut_payload >
exact_M_Farey_spectrum >
annular_14_prism_payload >
octahedral_current_carrier >
K33_product_wall >
hyperbolic_reciprocal_signature >
totient_euler_curvature >
raw_5_6_7_numerology.
```

With this declared gauge the carrier tournament is intended as transitive.
Any edge flip against the path is evidence that a geometry analogy has hidden
packet content and should be split into multiple sidecars.

## Assumption Challenge

This audit explicitly considered runners, tournament sizes, primes, triangle
cone points, Schlaefli `q`, Platonic solids, tiling patches, Farey nodes,
support-six residues, octahedral face curls, packet proof obligations, and
observer cuts as possible vertices.

The selected vertices are proof carriers because raw `5`, `6`, `7` labels do
not preserve the LRC14 predicate.  They destroy exact scale, endpoint owner,
safe topology, magnitude spectrum, packet route, and certificate sidecars
unless those are recorded separately.

## Next Pull

1. Add `geometry_regime_signature` to selected HYP-2963 packet families and
   test whether the typed regime predicts route after exact `M` and owner
   sidecars are retained.
2. Build the HYP-2943 annular `14` prism/antiprism carrier and compare it
   with AP/GW and the known boundary atoms.
3. Compare S131's `2/27` planar two-block strip with the `3/41` `K_{3,3}`
   product wall inside the packet bank.
4. Connect the HYP-2887 octahedral face-curl/divergence carrier to
   support-six packet coimages.
5. Separate archive motifs from proof carriers in every future use of
   `5`, `6`, `7`, `14`, `28`, `42`, and `1/42`.
