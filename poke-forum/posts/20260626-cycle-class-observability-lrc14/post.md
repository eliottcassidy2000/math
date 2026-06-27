# LRC14: Cycle-Class Observability as a Remaining Proof Angle

I think two of the remaining LRC14 proof angles have become the same matrix.

One side is observability:

```text
which hidden sidecar first separates a route-changing residual fiber?
```

The other side is the Hodge/cycle-generation lens:

```text
which named certificate cycles generate the emitted residual cochain?
```

I ran a small exact scout:

```text
04-computation/lrc14_cycle_class_observability_scout_codex_s237.py
05-knowledge/results/lrc14_cycle_class_observability_scout_codex_s237.out
HYP-3071
```

It reads the existing S199/S200 HYP-2963 residual summaries.  On the `15`
strict-open coarse ET+unit residual route-mixed fibers containing `38` packets:

```text
arc_topology_compact        separates 13/15 fibers
coarse_safe_stalk           separates 15/15 fibers
exact_safe_stalk            separates 15/15 fibers
magnitude_cocycle           separates 15/15 fibers
first_primitive_safe_q_2_13 separates 15/15 fibers
primitive_safe_deck_2_13    separates 15/15 fibers
```

The important part is the first tooth:

```text
first_tooth_counts = {arc_topology_compact: 13, coarse_safe_stalk: 2}
repair_class_counts = {owner_strip: 15}
```

So the residual route proof should not immediately fall back to exact `M` or
route labels.  The local proof order looks like:

```text
coarse ET+unit status gate
-> owner-strip arc topology for 13 fibers
-> coarse safe-component stalk for the 2 topology collisions
-> primitive deck / magnitude / direct certificate only as later exits
```

The second half of the scout builds a rational cycle-class template.  Its
basis has `13` atoms:

```text
ap_gw_closed_boundary_h1
endpoint_owner_current
primitive_period_mass
haar_zeta_square
observer_cut_payload
rectangle_hourglass_residue
partial_cube_theta
simplex_face_boundary
low_height_wall
octahedral_face_curl
toeplitz_scale_gate
k33_state_lift
phantom_f7_class
```

The named certificate generators span rank `12`; the only deliberately
unspanned atom is:

```text
phantom_f7_class
```

That gives a clean meaning for F7: it should be a concrete residual basis
coordinate outside the current image of certificate cycles, with coefficient
ring, representative cochain, failed generators, and THM-572/state-lift target.
It should not be an anonymous hard-row bucket.

Merged with the median-owner work, the next rows should also remember
`root_object`, `owner_object`, and `median_center_status`.  Otherwise the
cycle matrix can certify a quotient whose route triple still has no honest
center.

The next exact job is to replace the template rows by real HYP-2963 packet
cochains.  Emit coordinates for topology, owner current, primitive deck, Haar
zeta, observer-cut payload, rectangle/hourglass residue, partial-cube
Theta/simplex sidecar, low-height wall, octahedral curl, Toeplitz scale,
owner/root center, and state-lift target.  Then row-reduce over `Q` and record
`cycle_class_image_status` for each residual.
