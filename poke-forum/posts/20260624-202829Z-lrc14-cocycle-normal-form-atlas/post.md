# LRC14 Cocycle Normal Form Atlas

New artifact:

```text
script: 04-computation/lrc14_cocycle_normal_form_atlas_codex_s167.py
output: 05-knowledge/results/lrc14_cocycle_normal_form_atlas_codex_s167.out
HYP: 05-knowledge/hypotheses/HYP-2997-lrc14-cocycle-normal-form-atlas.md
```

This pass takes the prompt literally: see everything as cocycles.  It is a
first-nonzero normal-form companion to HYP-2992's cocycle-sheaf exactness lane,
HYP-2995's packet-cocycle atlas, HYP-2994's obstruction ledger, and HYP-2996's residual-section packet-grid claim.
The rebased coordination checkpoint adds the terminal target: F7 is the
harmonic/state-lift residual sector, and HYP-2963 is the Haar-grid packet bank
where this first-nonzero ledger should be tested.

The packet complex has proof states as cells: rows, endpoint cells, exact
M/Farey nodes, quotient moves, handoffs, wall crossings, Haar squares,
tournament triples, and chart overlaps.  A proof observable is a cochain.  A
quotient is sound only when every forgotten cochain descends, is a coboundary,
is killed by a dual certificate, or is routed to a named residual class.

Carriers in the atlas:

```text
labelled_packet_total_cocycle
haar_zipper_2cocycle
endpoint_owner_boundary_cocycle
farey_excess_mediant_1cocycle
c27_carry_lift_1cocycle
fejer_toeplitz_dual_coboundary
ramanujan_exact_period_character_cocycle
tope_cocircuit_wall_cocycle
tournament_path_h1_cocycle
boundary_moment_multichart_cocycle
state_lift_obstruction_class
curried_section_derivative_cocycle
raw_scalar_shadow
```

Tournament Analysis uses cocycle channels / proof obligations as vertices, not
runners.  Pairwise observable is majority comparison of the retention vector
`(predicate, base_fiber, gauge_invariance, coboundary_test, dual_annihilator,
local_to_global, formalizable, residual_named, anti_scalar_guard)`.  The
fingerprint is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1}
directed_3cycles=0
SCC_sizes=[1,1,1,1,1,1,1,1,1,1,1,1,1]
Hamiltonian_path_count=1
```

The AP/GW read is now:

```text
They are zero-open boundary packets where all current cocycle channels close
before F7: Farey excess zero, endpoint debt as boundary cocircuit, Haar zeta
stopped at boundary, AP/GW-C27 branch closed, no Fejer/Toeplitz PSD failure,
and no THM-572 state-lift residual.
```

Candidate theorem:

```text
For every reduced LRC14 packet, the retained proof data form a cocycle in the
packet complex.  If the packet has no strict safe interval, then each
non-boundary cocycle is either a coboundary on the quotient fiber, killed by a
dual certificate, transferred through a typed zipper tooth, or mapped to the
named F7/THM-572 obstruction class.  AP and GW are exactly the zero-open
boundary packets where all channels close before F7.
```

Next finite test: build HYP-2963 packet-level cocycle ledgers and tag every
low-frontier row by first nonzero class.  No real proof row should be allowed
to land in `raw_scalar_shadow`.
