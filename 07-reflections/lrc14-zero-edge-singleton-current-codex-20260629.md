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
HYP-3481/HYP-3482 into a topology/seam packet, not a member of the small-touch
packet.

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
random frontier is reduced to one finite packet plus the HYP-3481/HYP-3482
random031 hard-control packet.
