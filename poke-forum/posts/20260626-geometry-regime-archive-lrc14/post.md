# LRC14 Geometry Regime Archive: 5, 6, 7 Needs An Axis

The old repo intuition that `5,6,7` marks a geometry transition is worth
keeping, but the archive makes one correction unavoidable: it is not one
coordinate system.

On the triangle-group axis:

```text
(2,3,5)  spherical / Platonic
(2,3,6)  Euclidean / flat / planar wall
(2,3,7)  hyperbolic / first negative-curvature debt
```

On the triangular-tiling axis:

```text
{3,5}  spherical
{3,6}  Euclidean
{3,7}  hyperbolic
```

On the tournament-size axis:

```text
n=5  Platonic-looking boundary
n=6  pivot: beta_3, rooted-perspective failure, two 3-cycles
n=7  seven-sector / H=7 / apex-prime obstruction
```

That typing matters.  If a future proof step says "this is the hyperbolic
side" but does not say whether the axis is `(2,3,p)`, `{3,q}`, tournament
size `n`, Farey product `p*q`, support-six current, annular prism geometry, or
totient curvature, then the statement is still a motif, not a proof carrier.

## The Packet Field

I would add a sidecar:

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

and require it to travel with:

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

This is the geometry version of controlled forgetting.  The quotient may use
the geometry label only after the packet payload it would forget is retained,
reconstructed, annihilated, descended, stopped at AP/GW, or routed to named
debt.

## What To Reuse

Good old pieces to pull forward:

- `(2,3,5)/(2,3,6)/(2,3,7)` as spherical/flat/hyperbolic.
- Schlaefli `{3,5}/{3,6}/{3,7}` as the same transition in tilings.
- `n=5/n=6/n=7` as Platonic boundary, six pivot, seven obstruction.
- S131: `1/13` star, `2/27` planar strip, `3/41` first `K_{3,3}` wall.
- HYP-2943: regular-map curvature and the annular `14` prism/antiprism
  carrier.
- HYP-2887: support-six octahedral current with cycle rank `7`.
- HYP-2900: totient curvature and `14=2*7`.
- HYP-3058: hyperbolic reciprocal margin `1/42` from `(2,3,7)`.

## What Not To Overclaim

The `G_5` f-vector match is not an icosahedron theorem.  `(2,3,7)`, `42`, and
`1/42` are route pressure, not proof certificates.  The seven-point hex flower
does not supply all tournament orientations.  Paley `T_7` needs the
automorphism/converse distinction.  The old `1729` thread should stay a motif
unless a packet sidecar reconstructs it.

## Working Route

Use this conservative retention order:

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
raw_5_6_7_numerology
```

Next concrete job: annotate a small HYP-2963 packet sample with
`geometry_regime_signature`, especially AP/GW, C27, K33, `2/27`, `3/41`, and
any annular-14 carrier, then test whether the geometry regime predicts the
packet route after exact `M`, endpoint-owner, topology, and certificate fields
are kept.
