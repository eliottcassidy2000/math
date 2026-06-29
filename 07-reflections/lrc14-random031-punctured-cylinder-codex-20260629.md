# LRC14 random031 punctured-cylinder reflection

This session reframed the unique hard/currentless overlap,
`random_covering_031`, as a punctured-cylinder/seam-complement packet rather
than a defective graph-current row.

The exact run found four isolated dead islands in two mirror pairs, owner
spans `(22,20,20,22)`, hard components `(43,54)`, a max-delta seam carrying
the full seven-owner layer `(23,45,93,113,147,169,173)`, and a connected rescue
graph on the six non-apex owners.  The HYP-3460 phase pullback supplies the
decisive asymmetry: `282` phase witnesses, `0` hard-seam gate hits, but `12`
same-component lower-delta hits on the hard components split as six per
mirror branch.

The useful mental picture is now:

```text
random031 = mirror-punctured cylinder
four dead islands = isolated punctures
hard pair = forbidden seam
phase witnesses = flow on the seam complement
lower-delta hits = mirror bypass on the same seam components
```

This connects cleanly to the older recursion dichotomy.  The owner rim is
additive, with pairs `(23,45)`, `(93,113)`, `(147,169)` plus apex `173`; the
phase coordinate is multiplicative/two-adic, `u=2t mod 1`.  The hard seam is
the place where these two recursions disagree.  That makes it a candidate
gluing obstruction, not a scalar wall.

Next experiment: build a seam-complement graph by removing the max-delta hard
pair, route the `282` phase witnesses through branch-compatible lower-delta
gates, and test whether each flow component reaches a low-rank escape before
it would need the forbidden seam.  Repeat for all HYP-3477 hard orbits to
separate genuine phase walls from zero-hit boundary seams.

Tournament Analysis used proof carriers as vertices:

```text
C00_forbidden_seam_complement_flow
-> C01_mirror_punctured_cylinder_model
-> C02_seven_owner_layered_seam_word
-> C03_phase_branch_bypass_worldlines
-> C04_additive_rim_vs_doubling_fold_lens
-> C05_dead_island_owner_pairs
-> C06_raw_counts_shadow
```

The quotient guardrail is sharp: raw counts forget the seam, punctures,
mirror pairing, owner layers, and bypass.  A legal proof packet must keep at
least those coordinates or reconstruct them.
