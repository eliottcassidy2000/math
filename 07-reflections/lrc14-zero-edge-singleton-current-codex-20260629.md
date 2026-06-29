# LRC14 Zero-Edge Singleton-Current Audit

HYP-3480 closes the gap left by the HYP-3478 geometry audit more tightly than
expected.  HYP-3478 split the six small-touch/no-hard rows into four
branch-unit minimum rows and two cover-delta minimum rows.  HYP-3480 checks the
component level instead of the row minimum and finds the same terminal carrier
on all six rows.

The exact result is:

```text
target_components_with_complete_branch_unit_touch=14/14
target_mirror_pairs_with_branch_unit_mirror_gate=7/7
control_components_with_complete_branch_unit_touch=0/4
control_mirror_pairs_with_branch_unit_mirror_gate=0/2
```

So the useful theorem target is no longer "four clean rows plus two sidecar
rows" unless the formal proof insists on absolute minimum gates.  Rows
`random_covering_039` and `random_covering_074` still have cover-delta
absolute minima, but each dead singleton component also has a complete
branch-unit E/branch touch, and each mirror component pair has a
mirror-compatible branch-unit gate pair.

This matters because HYP-3476 already showed projection-current cannot help
these rows: the dead-cover projection has no edges.  The right finite packet is
component-local:

```text
swapped singleton B0/B1 owner pair
+ mirror-compatible complete branch-unit E/branch gate pair
+ route sidecar R
+ owner residue/two-adic word
```

The control row explains why the sidecars must stay explicit.  `random031` also
has zero projection edges and singleton mirror pairs, but it has no complete
branch-unit touches and no mirror-unit component pair.  It remains the
HYP-3455/HYP-3460/HYP-3479 hard/currentless gluing clause, now refined by
HYP-3486/HYP-3485/HYP-3484/HYP-3483/HYP-3482/HYP-3481 into a
seam-complement/fiber/forbidden-seam/topology packet, not a member of the
small-touch packet.  HYP-3490 adds the upstream reason the seven-row family
cannot close by larger adjacent-label pair currents: all touched blocker labels
are private.

Formalization checkpoint: `TournamentH7.LRCSingletonCurrentLedger` now records
the seven-row audit, the `4+2+1=7` terminal split, and the exact
component/mirror-pair count arithmetic.  The stored build output has no axioms
for the target-row count, random031 control classification, audited-row
partition, complete component touch, complete mirror-pair gate, dispatch
completeness, and dispatch/count matching hooks.

Tournament Analysis used proof carriers rather than runners or raw row names:
mirror-unit singleton-current packet, complete component-touch certificate,
random031 control clause, route sidecar `R`, cover-delta shadow,
owner-residue/two-adic sidecar, raw zero-edge count, and raw row list.  The
guardrail is strict: a quotient that remembers only "zero projection edges" or
"singleton components" falsely merges the six discharged candidates with
random031.  It must retain component-local branch-unit touches and the route
label.

Next proof work: prove the mirror-unit singleton-current lemma and wire it into
the route partition after HYP-3476/HYP-3479.  Then the non-AP currentless
random frontier is reduced by HYP-3490 to one six-row singleton-current packet
plus the HYP-3486/HYP-3485/HYP-3484/HYP-3483/HYP-3482/HYP-3481 random031
hard-control packet.
