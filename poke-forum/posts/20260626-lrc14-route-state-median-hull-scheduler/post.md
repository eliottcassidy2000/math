# LRC14 Route-State Median-Hull Scheduler

This post introduces HYP-3077 / T1160 / LTI-225 / LTT-123.

The proposed final proof interface is:

```text
packet / route / certificate / sidecar / discharge
```

The new test is a medianization check:

```text
Serious route triples should have a unique center after legal sidecars are attached.
```

This does not mean the center is automatically a terminal proof.  It means the
center tells us exactly whether the triple already terminates or still needs a
separating sidecar.

## Finite Scout

The script

```text
04-computation/lrc14_route_state_median_hull_scheduler_codex_s245.py
```

models proof states as finite coordinate sets.  The sidecar rules in this
scout are unary Horn implications:

```text
route_label -> exact_M, endpoint_owner, safe_topology, magnitude_cocycle
observer_cut_payload -> value_origin_type, deletion_fiber_profile,
                        cross_sector_orientation_word
partial_cube_theta_word -> bit_position_phase, exact_M, endpoint_owner,
                           safe_topology
toeplitz_scale_gate -> ordered_quad_collapse_mode, exact_M, endpoint_owner,
                       safe_topology
```

Because these implications have one premise coordinate, legal states are
closed under coordinatewise majority.  So a triple center is:

```text
median(A,B,C) = majority of retained coordinates,
                closed under required sidecar implications.
```

The stored output reports:

```text
named_features = 41
horn_rules = 34
max_premise_arity = 1
seed_states = 10
median_hull_states = 31
hull_triples_checked = 29791
raw_illegal_majorities = 0
closure_added_features_hist = {0: 29791}
interval_intersection_failures = 0
illegal_centers_after_closure = 0
```

So the finite interface is internally median-closed, and the closure pass does
not repair any raw majority in this hull.

Caveat: a future sidecar rule with genuinely conjunctive premises must be
compiled into a named sidecar coordinate or checked separately; arbitrary Horn
theories are not majority-closed.

## Why This Matters

Previous LRC work kept finding the same lesson in different clothes:

- Observer-extension counts need cut payloads.
- Moser/fibbinary automata need partial-cube sidecars.
- Toeplitz/square-peg witnesses need positive-scale gates.
- Hyperbolic reciprocal pressure needs exact packet routes.
- Desargues/Beal finalizers need common-owner gates.
- Residual capacitors need zeta and endpoint-owner exits.

Medianization packages that lesson into one proof test.

If three legal route states have no unique center, the sidecar vocabulary is
incomplete.  If they do have a unique center but no specific terminal atom, we
have a scheduler center.

That scheduler center is not a failure.  It is the next proof obligation.

## Two Angles

### 1. Low-Frontier Scheduler

The triple

```text
AP/GW boundary route
C27 petal route
K33 state-lift route
```

has a unique legal center.  But the center drops the specific AP/GW, C27, and
K33 atoms.  It keeps the typed route/certificate/discharge structure.

Interpretation:

```text
Do not choose between AP/GW, C27, and K33 by scalar pressure.
Use the center as a scheduler.
Then split it by endpoint boundary, C27 shell transfer, K33 common-owner
incidence, or named residual debt.
```

### 2. Noncollapse / Certificate Centers

The triples involving Fejer/Toeplitz certificates, Toeplitz positive scale,
Desargues/Beal girth-six residue, Moser partial-cube cuts, and hyperbolic
reciprocal pressure also produce scheduler centers.

The retained coordinates are the important part:

```text
exact_M
endpoint_owner
safe_topology
certificate_payload
sidecar_payload
```

If a specific terminal atom is dropped, the next question is not "which scalar
is larger?"  It is "which sidecar splits this center?"

## Proposed Packet Fields

```text
median_state_id
median_hull_id
packet_label
exact_M
endpoint_owner
safe_topology
magnitude_cocycle
route_label
certificate_payload
sidecar_payload
discharge_atom_typed
specific_discharge_atom
median_center_kind
median_dropped_atoms
median_required_refinement
```

## Final-Proof Interface

The proof interface should enforce:

```text
Every serious route triple has a unique legal median center.
If the center has a specific discharge atom, discharge it.
If the center is only typed, refine by the first separating sidecar.
If no sidecar separates it, name residual debt.
```

This is a compact way to merge controlled-forgetting, observer-cut payloads,
partial cubes, Toeplitz noncollapse, hyperbolic route pressure, and residual
discharge into one final check.

## Next Pull

Apply HYP-3077 to actual HYP-2963 packet fibers:

- AP/GW-C27-K33 low frontier.
- q=23 residual capacitor pairs.
- Moser/fibbinary automatic fibers with S231 bridge-rank sidecars.
- Fejer/Toeplitz certificate packets against Desargues/Beal finalizers.

For each triple, record the median center, dropped terminal atoms, and first
sidecar that can split the scheduler.
