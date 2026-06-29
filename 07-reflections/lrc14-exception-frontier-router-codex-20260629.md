# LRC14 Exception-Frontier Router Reflection

HYP-3476 was a useful negative simplification.  I expected the hard mirror
orbits from HYP-3475 and the currentless rows from HYP-3472 to overlap in a
messy family.  They do not.

The exact intersection is a singleton:

```text
random_covering_031
```

That is not a new monster; it is the row HYP-3455 already made finite as a
seven-owner mirror-gluing clause, and HYP-3460 already showed the phase-color
pullback avoids the max-delta hard gates.

This reframes the remaining proof work.  The six hard rows outside random031
are not currentless.  They already have separating E/branch boundary currents,
so their large mirror debt should be treated as a sidecar to discharge, not as
the terminal obstruction.  Conversely, the six random currentless rows outside
random031 have no hard mirror orbit and only small touching gates: best
adjacent delta is `2` or `3`.  That suggests a bounded owner-current or
two-adic/SPEC lemma rather than a hard-gate gluing theorem.

The AP84 base pair is separate again: edge cut exists, separating current does
not, and the closed endpoint/corridor/color packet should handle it.

The quotient lesson may be the most important part.  HYP-3474's existing
colored-gate axes `K,N,T,S,F,C,M,A` do not preserve the six-label terminal
route partition.  So route label `R` is a real proof sidecar.  A future Lean or
paper proof cannot infer the terminal discharge lane from the colored-gate
quotient alone unless it proves a reconstruction theorem.

After rebasing, I read the incoming S319 gate unit-delta split handoff.  It
uses a conflicting local `HYP-3472` name, but the content is useful: it splits
minimum E/B gates into `branch_unit_delta=110`, `both_unit_delta=1`, and a
`19`-row delta-sidecar packet.  That sidecar intersects the HYP-3476 frontier
at:

```text
022, 039, 074, 085, 113
```

This is signal, not a replacement.  S319's sidecar is an upstream gate-kind
coordinate; HYP-3476's `R` sidecar is a downstream terminal route coordinate.
In particular, `random_covering_031` is still the unique hard/currentless
overlap, and several HYP-3476 frontier rows are outside S319's delta-sidecar
packet.

Tentative packet theorem:

```text
After DeadCoverEBranchSoundness, every row is routed by R:
ordinary_current
| hard_currented_mirror_debt
| small_touch_no_hard_debt
| AP84_edge_only_packet
| random031_overlap_clause
| nondead_exit.
```

The next useful computation is probably not another broad scan.  It is a
finite proof attempt on the six small-touch rows:

```text
random_covering_001
random_covering_039
random_covering_062
random_covering_074
random_covering_086
random_covering_101
```

Their best gates are all touch-only with no projection edge removal.  The
question is whether their dead-cover projection is already separated by a
two-step owner-current, endpoint-spine wall, exact-period lift, or SPEC/Rprime
certificate once the small adjacent labels are retained.
