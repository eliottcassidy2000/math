---
id: HYP-2991
title: LRC14 Haar zipper cocycle and discrepancy handoff theorem
status: PROOF-INTERFACE / local cocycle audit and zipper theorem target; not a proof
source: codex-2026-06-24-S166
related:
  - HYP-2990
  - HYP-2989
  - HYP-2988
  - HYP-2987
  - HYP-2986
  - HYP-2985
  - HYP-2984
  - HYP-2595
  - HYP-2594
  - HYP-2450
  - HYP-2736
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_haar_zipper_cocycle_codex_s166.py
  - 05-knowledge/results/lrc14_haar_zipper_cocycle_codex_s166.out
---

# HYP-2991: LRC14 Haar Zipper Cocycle

HYP-2991 is the local computable refinement of HYP-2990's abstract zipper
atlas.  HYP-2990 says a quotient may forget a coordinate only when the
coordinate is invariant, reconstructible, annihilated by a dual certificate, or
routed to a named residual.  The Haar-product square identifies the coordinate
that the fixed-margin LRC14 quotient is trying to forget:

```text
zeta(T) = T00 - T01 - T10 + T11.
```

This is the two-dimensional Haar product on dyadic children, the elementary
fixed-margin switch coordinate, and the local orientation of the tournament
tiling square.  Row and column margins are not enough.  The theorem-facing
object is the labelled packet fiber together with this mixed cocycle, unless
the cocycle has been cancelled, exposed, handed off, descended, or named as
state-lift debt.

## Computation

Script:

```text
04-computation/lrc14_haar_zipper_cocycle_codex_s166.py
```

Stored output:

```text
05-knowledge/results/lrc14_haar_zipper_cocycle_codex_s166.out
```

The local audit enumerates all nonnegative `2 x 2` tables of total at most
`10`.  It records `1001` tables, `506` row/column margin fibers, and `285`
nontrivial margin fibers.  The augmented key

```text
(row sums, column sums, total, zeta)
```

has `1001` keys and `0` duplicate augmented keys.  Inside nontrivial fixed
margin fibers, the zipper step gcd is `4`.  The largest audited span is the
balanced margin fiber `((5,5),(5,5),10)`, with zeta values

```text
-10, -6, -2, 2, 6, 10.
```

Conclusion: margins alone collide, while margins plus `zeta` are a complete
local coordinate on these `2 x 2` fixed-margin fibers.

## Dyadic Product Rule

The depth-4 dyadic Haar rectangle census gives `961` rectangles and `923521`
ordered products.  The interaction classes are:

```text
orthogonal_zero             871992
same_tile_boundary_atom        961
vertical_owner_strip          6076
horizontal_owner_strip        6076
cross_zipper_handoff         19208
nested_refinement            19208
```

The nonzero non-atom classes are sign-balanced before LRC packet labels break
the symmetry.  This is the discrepancy-theory form of the tournament tiling
model: the raw product algebra cancels, and only labelled asymmetry can produce
proof pressure.

## Zipper Teeth

HYP-2991 proposes the following local teeth for the LRC14 zipper theorem:

```text
Z0 orthogonal cancellation | disjoint Haar coordinate          | discard by discrepancy orthogonality
Z1 boundary stop           | same-tile indicator               | AP/GW cocircuit or named boundary atom
Z2 owner strip             | one coordinate fixed/refined      | endpoint owner, Fejer center, or Ramanujan label
Z3 cross handoff           | opposite coordinate nesting       | zipper arrow between proof clocks
Z4 nested descent          | same-direction nesting            | family compression or state-lift descent
Z5 residual cocycle        | unpaired zeta survives all exits  | F7 / THM-572 state-lift obligation
```

Candidate theorem:

```text
On every labelled LRC14 packet fiber, each local mixed Haar cocycle is either
paired by color-compatible discrepancy cancellation, stopped at a boundary
cocircuit, handed to an owner/period/certificate clock, descended to a smaller
packet family, or converted into a named state-lift residual.  If no tooth
applies, the packet is the F7 sector.
```

This makes HYP-2987's `F7` less mysterious.  It is not an anonymous leftover.
It is an unpaired local cocycle after every admissible zipper tooth has failed.

## Tournament Analysis

Tournament vertices are proof carriers plus the local zipper cocycle, not
runners.  The pairwise observable is retained LRC predicate, fixed fiber, mixed
sign, endpoint topology, certificate handoff, exposure test, state route, and
formal check.  The switch/gauge is majority comparison of the retention vector;
ties use the declared Hamiltonian path.

The computed carrier tournament is transitive:

```text
haar_zipper_cocycle > certificate_handoff_atlas > exposure_kernel_audit >
tope_cocircuit_wall > color_resonance_discrepancy >
admissible_smoothing_clock > fixed_margin_tiling_shadow >
raw_component_count_K
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles=0
SCC_sizes=[1,1,1,1,1,1,1,1]
Hamiltonian_path_count=1
```

The interpretation is strict: HYP-2594's `K` is a lossy continuous boundary
count, HYP-2595's color-compatible resonance is the global compatibility test
for mixed modes, and HYP-2991's `zeta` is the local mixed sign coordinate that
must be retained or discharged before either quotient can carry a theorem.

## Assumption Challenge

The pass explicitly considered runners, dyadic rectangles, endpoint walls,
row/column margin fibers, switch moves, color resonances, proof carriers, Fejer
atoms, state-lift obligations, and theorem arrows as possible tournament
vertices.  The chosen vertices are proof carriers plus the local cocycle.

Preserved data:

```text
strict-witness predicate
fixed labelled packet fiber
mixed Haar sign / zeta
endpoint and topology labels
named certificate exits
```

Destroyed data:

```text
raw runner identity
continuous component count K
row/column margin shadows without zeta
```

The challenged assumption is that the tournament tiling model should be built
from runner vertices.  Here the actionable tournament is a proof-carrier
tournament, with `zeta` as the local switchboard coordinate.

## Next Work

1. Lift the local `2 x 2` cocycle to the actual HYP-2963 labelled packet grid.
2. Count independent color-compatible `zeta` switches in the HYP-2595 banks and
   compare them to `k+c_GP`.
3. Attach each nonzero owner-strip/cross/nested coefficient to the HYP-2987
   zipper arrow it should feed: family compression, admissible smoothing,
   Fejer/Ramanujan handoff, or state lift.
4. Define `F7` as an unpaired mixed-cocycle residual sector, then test whether
   every audited row has a Z0-Z4 exit before reaching it.
5. Build a formal no-free-slider lemma specialized to fixed margins:
   row/column margins may be quotiented only after `zeta` is reconstructed,
   annihilated, or routed.
