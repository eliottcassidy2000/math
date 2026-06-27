# LRC14 Geometry-Regime Archive Audit

Session: `codex-2026-06-26-S225`

The old "5,6,7 geometry" idea is one of the better archive sparks, but it
needs a typed axis before it can be used in LRC14.  The clean triangle-group
version is not `5=planar, 6=spherical, 7=hyperbolic`; it is:

```text
(2,3,5)  spherical / Platonic
(2,3,6)  Euclidean / flat / planar
(2,3,7)  hyperbolic
```

For Schlaefli triangular tilings this becomes:

```text
{3,5}  spherical
{3,6}  Euclidean
{3,7}  hyperbolic
```

The nearby tournament-size story is different:

```text
n=5  Platonic-looking boundary
n=6  pivot and first obstruction layer
n=7  seven-sector / H=7 / apex-prime obstruction layer
```

## What The Archive Says

The strongest old sources:

- `07-reflections/five-as-bridge.md`: `(2,3,5)` vs `(2,3,7)`, with
  `(2,3,6)` as the flat wall.
- `07-reflections/six-as-pivot.md`: `n=6` is where beta_3 appears and two
  disjoint 3-cycles can coexist.
- `07-reflections/the-five-platonic-tournaments.md`: useful Platonic analogy,
  but not a literal icosahedron theorem.
- `07-reflections/the-tessellation.md`,
  `07-reflections/the-modular-tournament.md`, and
  `07-reflections/the-three-infinity-structure.md`: `{3,infinity}` and
  modular/hyperbolic tournament language.
- HYP-2943: regular maps, Euclidean tiling wall, and annular `n=14`
  prism/antiprism carrier.
- HYP-2887: support-six octahedral current carrier with cycle rank `7`.
- HYP-2934/S131: `1/13` star, `2/27` planar strip, `3/41` `K_{3,3}` wall.
- HYP-2900: totient curvature and `14=2*7` residual prime pressure.
- HYP-3058: reciprocal-sum/hyperbolic sidecar, including `(2,3,7)` and
  margin `1/42`.

## The Sidecar

The reusable object is:

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

Use it when an analogy says "this is spherical", "this is flat", "this is
hyperbolic", "this is the six pivot", or "this is the seven obstruction".
Do not use the analogy without also keeping:

```text
exact_M
endpoint_owner_payload
safe_topology_payload
value_origin_type
observer_cut_payload
magnitude_spectrum_required
route_label
certificate_or_state_lift_payload
```

## Guardrails

The f-vector coincidence around `G_5` is not an `A_5`/Galois icosahedron
theorem.  `(2,3,7)`, `42`, and `1/42` are route pressure, not proof.  Paley
`T_7` needs an automorphism-vs-converse-symmetry distinction.  The seven-point
hex flower is a flat carrier with missing cross-petal orientations, not a
complete seven-tournament proof.  HYP-2943's hexagonal norm `7` is distinct
from centered hex-ring counts.  The old `1729` modular spine was weakened by
later audit, so it should remain an archive motif unless a packet field
reconstructs it.

## LRC14 Transfer

The productive transfer is not a new scalar.  It is a routing discipline:

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

If a geometry analogy cannot survive that path, it should stay in the archive
as a motif.  If it can, it becomes a packet sidecar with a concrete handoff:
q-witness, AP/GW boundary, C27 petal, K33 state lift, Fejer/Toeplitz
certificate, or named F7 debt.

## Next Work

The next useful pass is finite and concrete: annotate a small HYP-2963 packet
sample with `geometry_regime_signature`, especially AP/GW, C27 petals, K33
rows, the `2/27` and `3/41` unit-excess nodes, and any candidate annular-14
carrier.  Then test whether the geometry sign predicts route after exact `M`,
endpoint-owner, topology, and certificate fields are retained.
