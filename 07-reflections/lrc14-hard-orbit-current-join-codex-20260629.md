# LRC14 Hard-Orbit Current-Join Reflection

HYP-3475 made the hard mirror-orbit family finite but still left it looking
like a separate seven-row problem.  HYP-3472 made the boundary-current
exception family finite but left nine separating-current exceptions.  Joining
the two ledgers is the useful move.

After rebasing over HYP-3476 and HYP-3477, this join should be read as the
Lean-backed ledger for their shared route split: HYP-3476 isolates the
singleton hard/currentless overlap, and HYP-3477 audits the hard-orbit
discharge family.

The exact readout is sharper than expected:

```text
hard_orbits_delta_ge_7=8
hard_orbits_with_separating_current=7/8
hard_rows_without_separating_current=[random_covering_031]
ap84_hard_rows=[]
```

So the hard orbit and current-cut debts are not the same broad obstruction.
They intersect only at `random_covering_031`.  That changes the next proof
shape:

```text
hard_orbit_discharge
  <= separating_current_transfer + random031_named_clause
```

The AP84 base rows are not hard-orbit rows.  They stay in the AP84 sidecar
packet, where HYP-3462/HYP-3470 already carry the relevant structure.  The six
random current exceptions besides `random_covering_031` are also not hard
mirror-orbit rows, so they should be treated as low-delta current exceptions
rather than as part of the HYP-3475 hard-family proof.

Assumption challenge: I did not use runners or raw row names as tournament
vertices.  The selected vertices were joined proof carriers: hard mirror orbit,
separating current, singleton intersection ledger, and terminal dispatch.  The
quotient preserves whether a hard gate debt has a current-transfer exit.  It
destroys exact internal interval geometry, so typed orbit and best-current
sidecars remain attached.

The next serious step is to make the separating-current transfer mathematical
instead of finite-bank telemetry.  If that transfer is proved, the hard-orbit
part of the colored-gate route is reduced to the already named HYP-3455 /
HYP-3460 `random_covering_031` clause.

I also added `TournamentH7.LRCHardOrbitCurrentJoin` as a small Lean ledger for
the dispatch arithmetic.  It does not prove the geometry, but it makes the
finite count partition and `7+1=8` hard-orbit dispatch sorry-free and visible
to future formalization work.
