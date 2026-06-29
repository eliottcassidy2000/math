# LRC14 Colored Gate Mirror-Orbit Reflection

HYP-3471 made the color layer local: every dead row has a rank-`<=2`
E/branch survivor gate.  The next missing coordinate was not another color.  It
was the mirror orbit.

After rebasing, HYP-3472/HYP-3473/HYP-3474 occupy the boundary-current,
formalization, and quotient-legality lanes.  HYP-3475 should be read as the
mirror-orbit sibling: it supplies the involutive gate-extension sidecar those
lanes must either retain or discharge.

The exact bank result is useful:

```text
8702 survivor gates -> 4351 mirror orbits
fixed_orbits=0
unpaired_or_duplicate_gates=0
dead_rows_with_e_branch_low_rank_orbit=130/130
```

So the colored-extension quotient from HYP-3461 can be made concrete on the
current bank: work with two-gate mirror orbits carrying typed endpoint residues,
branch masks, adjacency, and cover deltas.

The strongest new proof obligation is smaller than the old random031 clause:

```text
discharge the 8 mirror orbits with cover-delta >= 7.
```

They lie on seven random rows.  `random_covering_022` has the unique delta-`8`
orbit and one additional delta-`7` orbit; `random_covering_031` is one of six
other delta-`7` rows.  This reframes HYP-3455: it is no longer a one-off
exception, but a representative member of the hard mirror-orbit family.

Assumption challenge: I considered runners, individual gates, endpoint colors,
typed mod-`14` words, component boundaries, and proof obligations.  The
mirror-orbit quotient preserves the gluing/debt predicate and exact low-rank
escape status, while destroying branch orientation of the individual member.
That orientation is retained as the ordered structural sidecar.
