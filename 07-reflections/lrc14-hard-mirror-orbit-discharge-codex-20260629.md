# LRC14 Hard Mirror-Orbit Discharge Reflection

HYP-3477 resolves the immediate HYP-3475 ambiguity: the hard mirror-orbit
family is not a new broad obstruction class.

The exact audit splits the eight delta-`>=7` mirror orbits into two proof
types.  Seven are ordinary projection-current cases: even when q=`14V` grid
witnesses hit the hard gate itself, the same row has a lower-delta E/branch
gate whose adjacent blocker labels already cut the HYP-3451 dead-cover
projection.  The unique delta-`8` orbit on `random_covering_022` is in this
ordinary class.

`random_covering_031` is the only exceptional hard orbit.  It has zero
compatible q=`14V` hits on the max-delta mirror pair, while its hard components
are still reached through lower-delta compatible gates.  The bounded
pair-current sidecar used here does not cut its dead projection, so the right
handoff is exactly HYP-3455 plus the active HYP-3476 pair-current diagnostic,
exception-frontier route sidecar, and later gluing/current packets, not another
scalar hard-delta search.

Proof-facing consequence: prove a lower-delta E/branch projection-current
discharge for non-random031 hard orbits, and isolate random031 as the finite
seven-owner gluing clause.  Do not state a universal phase-bypass theorem:
the phase-grid bypass is true for random031 only in this hard-family audit.

Rebase integration: HYP-3476 now has both the pair-current audit and the
row-level packet router, while HYP-3477 is the orbit-level hard-family
discharge.  Together they say AP pair-current debt is closed, random
edge-exception debt is zero-edge singleton debt, and hard mirror debt overlaps
failed boundary current only at `random_covering_031`.
