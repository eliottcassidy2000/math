# LRC14 Coordinate-Resurrection Sheaf Reflection -- codex S267

This note records the abstract reframing behind HYP-3118.  It is not a proof.
It is a map for what to measure next.

## Core Reframe

The recent LRC14 frontier keeps repeating a pattern:

```text
good scalar or theorem-shaped signal
  -> quotient forgets proof coordinates
  -> missing-input ledger catches the failure
  -> sidecar proposal restores some but not all coordinates
```

The HYP-3118 move is to make the repair relation the object.  A shortcut is
not "bad" merely because it destroys information.  It is bad only when it does
not name which coordinate it destroyed, what proof operation next needs that
coordinate, and which minimal legal sidecar resurrects it.

In the scout, the base stalk is:

```text
finite_address + observer_gluing + endpoint_owner + uniformity
```

The live sections are:

```text
root_ear_payload, relation_lattice_shape, component_bound,
cocycle_exactness, state_lift, pde_weak_form
```

That split is the most important output.  The result is not "Lee-Yang beats
Minkowski" or "Savitch beats PDE".  The result is: all of those are sections,
but no section is terminal until the base stalk is present.

## Map 1: Quotient To Resurrection

For each shortcut, measure:

```text
kept coordinates
destroyed coordinates
first proof operation that needs a destroyed coordinate
minimal legal repair cover
terminal exit or named debt
```

The scout's strongest pattern is that scalar-like shortcuts with no live
section need four sidecars.  Shortcuts that already carry a live section need
three.  The proof-circuit compiler already carries uniformity, so it also
needs three, but for a different reason: it detects missing coordinates before
supplying endpoint owner, finite address, observer gluing, and a live payload.

This suggests a new audit format for every future LRC shortcut:

```text
shortcut_id
quotient_preserves
quotient_destroys
core_gap
live_section_present
minimal_repair_cover_rank
adjoint_section_status
exit_kind
```

## Map 2: External Pattern To New Signal

Savitch's theorem is not just "space is small".  It suggests a midpoint
certificate ladder: the proof route may be legal if it can recompute local
segments instead of storing a global path.

Bravais lattices are not just "Minkowski volume".  They say that relation
shape, basis choice, and wall ownership matter before volume pressure can be a
proof carrier.

Lee-Yang is not just "root count".  The single value `p0` is only one point of
the whole miss-count PGF.  The curve, zero trajectory, nearest-zero radius,
argument path, and discriminant walls are the carrier.

The `phi4` density is not just "energy".  It is a stress test for moment walls,
boundary regimes, and zero modes.  Without boundary/operator data it is a
pretty analogy and not a proof route.

Ear-decomposition theorems are not one analogy.  They split into three gluing
grammars:

- directed ears certify strong observer connectivity;
- odd ears certify parity or ownership debt;
- nested ears certify a low-complexity branch/series-parallel shadow.

This is where tournaments re-enter: the useful vertices may be proof
obligations, repair covers, or gluing grammars, not runners.

## Wild Hypotheses Worth Testing

1. The four core coordinates are coloops of a repair matroid.  Any proof route
   that loses one must explicitly re-add it before terminalizing.
2. The live sidecars form exchangeable sections over the same base stalk, so
   many broad analogies differ less by strength than by which coordinate they
   resurrect cheaply.
3. The PGF zero trajectory is the natural Lee-Yang section, while one-runner
   ear payloads are its transition maps.
4. Bravais shape and Savitch midpoint certificates may be dual repair
   languages: one normalizes space of relations, the other normalizes space of
   proof paths.
5. The obstruction cocycle is not merely a fallback sidecar.  It may be the
   boundary operator of the coordinate-resurrection sheaf.
6. State-lift appears missing everywhere in the shortcut table but should not
   be promoted to the base stalk yet.  It may instead be the first nontrivial
   higher section: absent from cheap quotients, decisive for hard residuals.
7. The right extremality problem may be "which quotient has the smallest
   legal adjoint section" rather than "which statistic is maximized".

Incoming S266 proof-carrier work adds a concrete compatibility test for this
language.  HYP-2112 `Phi`, HYP-2108 `P`, HYP-2109 `L/M/R`, HYP-3023 magnitude
cocycle, HYP-3077 Horn closure, and HYP-3082 protected branch status should be
treated as named repair sections for exact-gap, endpoint-activation,
wall-crossing, route-purity, legality, and bridge-safety coordinates.  If the
coordinate-resurrection story is right, these gates should plug into the same
ledger as PGF roots, Bravais walls, Savitch ladders, and ear grammars.

## Concrete Next Measurement

Build a small residual-row table, likely starting from HYP-2963/HYP-3098 and
the HYP-3112 root/ear rows, with these fields:

```text
destroyed_coordinate_vector
core_stalk_presence
live_section_type
coordinate_resurrection_cover
repair_cover_rank
adjoint_section_status
pgf_zero_trajectory_signature
bravais_shape_wall_signature
midpoint_certificate_depth_profile
observer_ear_certificate_type
odd_ear_payload_parity_debt
nested_ear_branchwidth_shadow
phi_gap_sum
P_max_activation
LMR_terminal_state
magnitude_cocycle_route_purity
horn_legality_status
protected_branch_status
quartic_moment_wall_profile
terminal_exit_or_named_debt
```

The falsifiable target is also simple: every live residual should either hit a
known proof minterm or admit a legal repair that strictly reduces the
missing-coordinate vector.  If not, HYP-3118 has found a real obstruction
rather than only a vocabulary upgrade.
