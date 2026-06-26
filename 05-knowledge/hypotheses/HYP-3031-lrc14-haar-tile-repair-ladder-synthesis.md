---
id: HYP-3031
title: LRC14 Haar-tile repair ladder and tournament-tiling product dictionary
status: SYNTHESIS / proof-interface dictionary and theorem target; not a proof
source: codex-2026-06-26-S195
tangent: T1112
related:
  - HYP-3030
  - HYP-3029
  - HYP-3028
  - HYP-3027
  - HYP-3026
  - HYP-3025
  - HYP-3024
  - HYP-3023
  - HYP-3020
  - HYP-2992
  - HYP-2991
  - HYP-2989
  - HYP-2997
  - HYP-2995
  - HYP-2990
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3031: Haar-Tile Repair Ladder Synthesis

## Claim

The recent HYP-3023..HYP-3030 automatic-fiber/status/stalk/topology-gate work and the older
HYP-2989/HYP-2991/HYP-2992 Haar-product/tournament-tiling work are the same
controlled-forgetting pattern at different scales.

The local object is the `2 x 2` mixed cocycle

```text
zeta(T) = T00 - T01 - T10 + T11.
```

Equivalently:

```text
2D Haar product atom        h_I(x) h_J(y)
fixed-margin switch         [[+1,-1],[-1,+1]]
tournament staircase tile   local fixed-path tile flip / Walsh xor bit
quotient repair cochain     first lost coordinate on an automatic fiber
```

The recent automatic-word failures are therefore not separate from discrepancy
theory.  An unsafe quotient is a row/column shadow: it keeps one-dimensional
margins but forgets `zeta`.  A repair ladder is the global rule saying where
that forgotten mixed coefficient reappears: exact `M/q`, boundary topology,
packet labels, harmonic/dual certificates, or F7/THM-572 residual debt.

## Dictionary

The usable dictionary is:

```text
Haar orthogonal zero      -> dual/discrepancy annihilation
same-tile indicator      -> AP/Goddyn-Wong boundary atom
owner strip              -> endpoint owner, Fejer center, active normal fan
cross handoff            -> K33/THM-572 state lift or proof-clock transfer
nested refinement         -> family descent, covering moment, magnitude formula
unpaired zeta             -> F7 residual sector
```

The two-dimensional Haar product rule

```text
h_{I x J} h_{I' x J'} = (h_I h_{I'})(h_J h_{J'})
```

has the same proof discipline as fixed-path tournament tiling multiplication:

```text
chi_A chi_B = chi_{A xor B}.
```

Both are exact while the two-dimensional address is retained.  Both become
dangerous when collapsed to scalar counts, raw component counts, row/column
margins, automatic words, or tournament isomorphism classes without the lost
mixed coordinate.

## Recent-Agent Synthesis

HYP-2989 supplied the minimal square: diagonal and anti-diagonal packets have
the same row and column margins but opposite mixed Haar coefficient.  HYP-2991
named the lost coordinate `zeta` and made it a local zipper cocycle.  HYP-2992
expanded the square to dyadic rectangle interaction classes: most products are
orthogonal, and the surviving signed interactions are owner strips, cross
handoffs, and nested refinements.

HYP-3023 then found the same pattern inside the full HYP-2963 automatic-word
bank: automatic words and residue-terminal fibers are row/column shadows, while
the magnitude cocycle is the first tested non-route splitter with zero mixed
theorem-route fibers.

HYP-3024 added the discrepancy side: exact Erdos-Turan clocks split too finely,
but a coarse ET plus Henselian-unit gate has `0` mixed boundary/open fibers and
only `15` mixed route fibers on the full bank.  This says the coarse gate kills
the status-level mixed cocycle but does not classify every route.

HYP-3025 supplied the topology side: closed arc-Cech nerves and boundary
cocircuit facets say which same-tile atoms are AP/GW boundary structure and
which runner-quotient collapses forget topology.

HYP-3026 fused the carriers: barcode, safe-stick body, endpoint current, CRT
chart, magnitude, automatic sidecars, and ET/Hensel gates must travel as a
packet unless a quotient proves why one coordinate descends.

HYP-3027 supplied the global repair rule: exact `M` repairs open/boundary
status but not theorem route; `M+q` and boundary topology each leave one mixed
route pair; packet labels or the guarded non-route signature are route-pure in
the audited bank.  In Haar language, HYP-3027 is the first-nonzero `zeta`
locator after automatic row/column shadows are fixed.

HYP-3028 turned the remaining `15` coarse-gate mixed-route fibers into
certificate scheduling debt rather than counterexample pressure: the gate
already has `0` mixed boundary/open fibers, so route labels can be delayed
after status is protected.

HYP-3029 showed that the hard AP/GW automatic word also admits a local
safe-component stalk descent: exact largest-component stalks remove the target
fiber's mixed routes while retaining endpoint and peak-owner geometry.

HYP-3030 inserted the missing topology proof gate between those arithmetic and
local-stalk teeth: AP/GW are the residue-terminal boundary rows with the closed
arc-Cech full-cover cycle, while open cohabitants and coarse route-mixed
survivors have closed arc `beta1=0`.  For this synthesis, HYP-3030 says the
same-tile boundary case of `zeta` must be recognized before route scheduling.

## Theorem Target

Haar-tile repair ladder theorem:

```text
Fix a primitive LRC14 packet fiber after automatic/residue shadows.
For every local two-coordinate packet grid attached to that fiber,
the mixed Haar/tournament-tile cocycle is either
  (a) orthogonally annihilated,
  (b) a same-tile AP/GW boundary atom,
  (c) repaired by an owner strip / endpoint / Fejer / normal-fan coordinate,
  (d) handed across proof clocks by a cross-handoff,
  (e) descended by nested refinement to a family magnitude or covering formula,
  (f) or emitted as named F7/THM-572 residual debt.
```

Equivalently, a non-AP/GW zero-open residual cannot be invisible to all
owner-strip, cross-handoff, nested-refinement, magnitude, topology, and packet
repair coordinates.

## Tournament Analysis

Vertices are proof obligations / repair teeth, not runners:

```text
raw_automatic_shadow
row_column_margin_shadow
coarse_et_unit_status_gate
residual_status_gate
exact_M_q_repair
arc_cech_boundary_topology
haar_zeta_packet
safe_component_stalk
magnitude_cocycle
guarded_packet_signature
```

Pairwise observable:

```text
preserves_boundary_open,
preserves_route,
retains_zeta,
retains_topology,
retains_packet_labels,
low_proof_cost
```

Switch/gauge:

orient `A -> B` when `A` preserves more of the first five proof coordinates;
use lower proof cost only as the final tie-breaker.  The tie Hamiltonian path
is the displayed vertex order from raw shadow to guarded signature.

Resulting proof-carrier tournament is transitive:

```text
guarded_packet_signature
> safe_component_stalk
> magnitude_cocycle
> haar_zeta_packet
> residual_status_gate
> arc_cech_boundary_topology
> exact_M_q_repair
> coarse_et_unit_status_gate
> row_column_margin_shadow
> raw_automatic_shadow
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

This is a theorem-use ranking, not a creativity ranking.  The point is that
cheap status gates are still useful early teeth, but any quotient that forgets
the mixed cocycle must prove reconstruction, annihilation, descent, boundary
equality, or named residual debt.

## Assumption Challenge

Alternate vertex sets considered: runners, gaps, dyadic rectangles, fixed
circle sections, section boundaries, endpoint walls, wall-crossing events,
residues, Fourier modes, cover arcs, matroid circuits, fixed-path tournament
tiles, quotient fibers, and proof obligations.

The chosen vertices are repair teeth because the preserved LRC predicate is
not raw runner identity.  It is boundary/open status plus theorem route after
declared quotienting.  This choice destroys raw runner names, continuous
component counts, row/column margins without `zeta`, and automatic-word
language identity unless those labels are reattached by a later repair tooth.

The challenged assumption is that tournament tiling should enter LRC14 as a
new scalar or an isomorphism-class census.  Here it enters as a product-rule
guardrail: the fixed-path tile algebra tells us which two-dimensional
coordinate a quotient is trying to erase.

## Next Pull

1. Add a packet-level `zeta_repair_class` field with values
   `orthogonal_zero`, `same_tile_boundary`, `owner_strip`, `cross_handoff`,
   `nested_refinement`, and `residual`.
2. For the two HYP-3027 mixed pairs, build the actual two-coordinate packet
   grid and identify which Haar/tournament-tile class separates them.
3. Treat HYP-3030's topology gate and HYP-3024's coarse ET+unit gate as the
   cheap status teeth, then require HYP-3027's repair cochain to name the first
   route-level nonzero cocycle.
