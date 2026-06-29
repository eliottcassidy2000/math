# LRC14 Dead-Cover Boundary-Current Reflection

HYP-3472 pushed a different direction from the color-reservoir work: it treated
the HYP-3471 E/branch survivor gates as boundary-current objects against the
HYP-3451 dead-cover blocker projection.

The useful result is a three-level ladder, not a single theorem:

```text
dead rows with touching E/branch gate: 130/130
dead rows with projection edge-cut gate: 123/130
dead rows with separating-current gate: 121/130
```

The first line is exact on the current bank and is the cleanest successor to
HYP-3471:

```text
dead_components(row)>0
  => rank<=2 E/branch gate touching the dead-cover projection.
```

The stronger graph-current transfers have named exceptions.  Projection edge
cut fails on seven non-AP rows:

```text
random_covering_001
random_covering_031
random_covering_039
random_covering_062
random_covering_074
random_covering_086
random_covering_101
```

Separating current fails on those seven rows plus the two AP84 base rows:

```text
covering_AP_with_84
ap_omit_12_tail_84x01
```

This is a good sign rather than a dead end.  The AP84 base rows already belong
to HYP-3462/HYP-3470.  The row `random_covering_031` was already HYP-3455's
seven-owner gluing clause, and HYP-3472 independently puts it in the edge-cut
exception set.  The remaining six random rows are now small named debt rather
than part of a vague non-AP transfer problem.

Assumption challenge:

- Runners were not the right vertices for this audit: they hide whether a gate
  touches the dead-cover projection.
- Raw survivor gates were still too coarse: `130/130` dead rows have one, but
  only `123/130` get an edge cut and only `121/130` get a separating current.
- The selected quotient used dead-cover blocker labels plus E/branch gate
  boundary labels.  It preserves the predicate "dead obstruction can route to a
  boundary gate" and destroys full interval geometry, exact wall order, scalar
  color counts, and same-row gate multiplicity.
- Alternate vertices considered: runners, gaps, fixed sections, section
  boundaries, wall crossings, residues, cover arcs, Fourier modes, matroid
  circuits, survivor gates, dead components, blocker labels, current vectors,
  graph cuts, and proof obligations.

Tournament Analysis used proof carriers:

```text
B00_projection_cut_gate
-> B01_separating_boundary_current
-> B02_dead_positive_e_branch_implication
-> B03_closed_ap84_packet
-> B04_random031_seven_owner_clause
-> B05_typed_gate_word
-> B06_raw_gate_count
-> B07_raw_dead_fraction
```

There were no directed `3`-cycles and one Hamiltonian path by the declared
payload-retention gauge.  The scalar shadows are last because raw gate counts
and raw dead fraction cannot see the seven edge-cut exceptions.

Next proof work:

1. Formalize the universal touch lemma.
2. Prove projection-edge transfer outside the seven named random exceptions.
3. Prove separating-current transfer outside those seven plus the two AP84 base
   rows.
4. Route the exception rows through HYP-3455, AP84 sidecars, owner-current,
   two-adic descent, signed-SPEC, or state-lift debt.
