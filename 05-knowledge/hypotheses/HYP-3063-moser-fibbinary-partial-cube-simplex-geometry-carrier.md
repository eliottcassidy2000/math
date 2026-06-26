---
id: HYP-3063
title: Moser-fibbinary partial-cube simplex geometry carrier
status: SYNTHESIS / forum-facing proof-carrier sidecar; not a proof
source: codex-2026-06-26-S227
tangent: T1145
script: 04-computation/lrc14_fibbinary_moser_partial_cube_codex_s228.py
result: 05-knowledge/results/lrc14_fibbinary_moser_partial_cube_codex_s228.out
related:
  - HYP-3062
  - HYP-3061
  - HYP-3058
  - HYP-3025
  - HYP-3023
  - HYP-3018
  - HYP-3016
  - HYP-3012
  - HYP-3011
  - HYP-3009
  - HYP-3008
  - HYP-3003
  - HYP-3000
  - HYP-2998
  - HYP-2963
  - HYP-2943
  - HYP-2887
  - HYP-2458
  - HYP-2454
  - HYP-2557
  - OPEN-Q-108
---

# HYP-3063: Moser-Fibbinary Partial-Cube Simplex Geometry Carrier

## Claim

Moser-de Bruijn and fibbinary rows should be upgraded from sequence names to
partial-cube/simplex sidecars before they are used in LRC14.  The old
HYP-3008/HYP-3012 lesson was finite-state:

```text
Moser-de Bruijn: base-4 digits 0/1, binary 1s only in even positions,
                 native transition x -> 4x.
fibbinary:       binary strings with no adjacent 1s,
                 native transition x -> 2x.
```

The new merge adds the geometric carrier:

```text
Moser-de Bruijn support = an even-coordinate Boolean cube.
Fibbinary support       = the Fibonacci cube of independent sets in a path,
                          hence a partial-cube language over hypercube bits.
Moser inside fibbinary  = a convex-looking coordinate subcube only after the
                          even/odd bit-position phase is retained.
```

This converts automatic-language sidecars into partial-cube sidecars:
Theta-classes, median/interval checks, gated subcubes, native transition, and
simplex face rank must be retained before the sequence shadow can be promoted
to a proof carrier.

## Doubled-Triangular/Simplex Layer

The sequence

```text
2, 6, 12, 20, 30, 42, ...
```

is `n(n+1)=2*T_n`.  It has two useful readings already present in the archive:

1. It is the directed edge count of the `n`-simplex 1-skeleton, since an
   `n`-simplex has `C(n+1,2)=n(n+1)/2` unoriented edges.
2. It is the exact `p=1` Faulhaber/triangular power-balance center in HYP-2454
   and the shared carrier `u=n(n+1)` in HYP-2458.

For LRC this is not a scalar numerology row.  It is the bridge between:

```text
hypercube / partial-cube bit coordinates
simplex directed-edge layers
ordered-pair observer sectors
Faulhaber odd-moment carriers
worry-modulus odd-square identities
```

The older identity

```text
8*C(n,2)+1 = (2n-1)^2
```

is the same bridge on the LRC worry-set side.  At `n=14`, `C(14,2)=91` and the
odd-square root is `27`, the familiar `2n-1` shell.  The doubled-triangular
layer should therefore be recorded as a coordinate layer, not used as a count
shortcut.

## Geometry-Regime Merge

The prompt's `5,6,7` geometry intuition should be kept, but with HYP-3061's
axis correction attached:

```text
Schwarz/tiling axis:
  (2,3,5) or {3,5}: spherical
  (2,3,6) or {3,6}: Euclidean / flat / planar wall
  (2,3,7) or {3,7}: hyperbolic

tournament-size axis:
  n=5: Platonic-looking boundary
  n=6: first pivot / rooted-perspective and beta_3 obstruction layer
  n=7: apex-prime / H=7 / seven-sector obstruction layer
```

The Moser/fibbinary carrier adds a third axis:

```text
cube-complex axis:
  Boolean cube:      unconstrained coordinate simplex/subcube language
  Fibonacci cube:    path-independent-set partial cube
  Moser cube:        even-coordinate Boolean subcube inside fibbinary
  geometry regime:   only usable after exact packet labels are retained
```

This lets the forum thread reuse old geometry without reopening old mistakes:
`{3,6}` is the planar tiling wall, `n=6` is a tournament pivot, and Moser/fibbinary
are cube-complex sidecars whose LRC value depends on exact packet payloads.

## Packet Sidecar

Recommended HYP-2963/HYP-3012/HYP-3061 fields:

```text
partial_cube_sequence_sidecar:
  automaton_language_id
  automaton_state_word
  native_transition
  bit_position_phase
  hypercube_dimension
  partial_cube_model
  theta_class_word
  gated_subcube_status
  median_interval_status
  simplex_face_rank
  directed_simplex_edge_count
  doubled_triangular_layer
  fibonacci_cube_window
  moser_even_coordinate_subcube
  geometry_regime_signature
  exact_M
  endpoint_owner_payload
  safe_topology_payload
  magnitude_cocycle
  route_label
  certificate_or_state_lift_payload
  destroyed_coordinate
  quotient_exit
```

Legal exits are the standard controlled-forgetting exits:

```text
retained
reconstructed
dual-annihilated
descended to a packet family
AP/GW boundary stop
C27/K33/state-lift handoff
named F7/THM-572 residual debt
```

## Exact Scout Addendum

S228 adds an exact carrier scout:

```text
04-computation/lrc14_fibbinary_moser_partial_cube_codex_s228.py
05-knowledge/results/lrc14_fibbinary_moser_partial_cube_codex_s228.out
```

It verifies the ordered-sector row:

```text
k=2..7: k(k-1)=2*C(k,2)=2*T_{k-1}=2,6,12,20,30,42.
```

It also separates the first count collisions by value origin:

```text
Gamma_5: vertices 13, edges 20, layers [1,5,6,1]
Gamma_6: vertices 21, edges 38, layers [1,6,10,4]

M_5: vertices 32, cube_edges 80, layers [1,5,10,10,5,1]
M_6: vertices 64, cube_edges 192, layers [1,6,15,20,15,6,1]
```

So `20` can be a Fibonacci-cube edge count or an ordered-sector count, and
`12` can be a Moser cube edge count or an ordered-sector count.  Those
coincidences are legal only after `value_origin_type`, transition, bit phase,
partial-cube/simplex layer, exact `M`, endpoint owner, route, and safe
component certificate are retained.

## LRC14 Transfer

Use this carrier where old automatic rows appeared as positive controls:

- HYP-3011/HYP-3025: first-13 fibbinary and Moser rows have positive safe mass
  and positive bars, so they are not AP/GW boundary atoms.
- HYP-3016/HYP-3023: raw automaton words and residue-terminal fibers mix
  routes; magnitude cocycle and packet labels split them.
- HYP-3012: fibbinary, Moser, and Ostrowski-Hadamard form a real tournament
  3-cycle, so a single "gap language" scalar hides proof-relevant structure.
- HYP-3061: geometry labels must declare their axis before they can interact
  with automaton or cube-complex labels.
- HYP-2454/HYP-2458/HYP-2557: the triangular carrier `u=n(n+1)` is a
  directed-simplex/odd-moment address, not a substitute for endpoint geometry.

The near-term proof task is to take one HYP-2963 packet sample containing AP,
GW, K33, C27 petals, covering, fibbinary, and Moser controls, then add the
partial-cube sidecar next to existing automatic, magnitude-cocycle, barcode,
closed-arc Cech, and geometry-regime fields.  A sequence/cube quotient is
admissible only if boundary/open status and route are constant after these
fields are retained or discharged.

## Tournament Analysis

Vertices are proof carriers and sidecar columns, not runners:

```text
labelled_packet_sheaf
automatic_magnitude_cocycle
partial_cube_theta_class_ledger
fibbinary_fibonacci_cube
moser_even_coordinate_cube
simplex_directed_edge_layer
geometry_regime_signature
faulhaber_odd_moment_carrier
doubled_triangular_scalar
raw_sequence_name
```

Pairwise observable:

```text
exact M
automaton state
native transition
bit-position phase
Theta-class retention
gated/median interval status
simplex face rank
geometry axis
safe topology
route/status handoff
certificate payload
```

Gauge: orient toward the carrier retaining more of the above payload and
leaving less unnamed route/status debt.  Tie Hamiltonian path:

```text
labelled_packet_sheaf >
automatic_magnitude_cocycle >
partial_cube_theta_class_ledger >
fibbinary_fibonacci_cube >
moser_even_coordinate_cube >
simplex_directed_edge_layer >
geometry_regime_signature >
faulhaber_odd_moment_carrier >
doubled_triangular_scalar >
raw_sequence_name
```

The challenged assumption is that the useful vertices are the sequence values
or the runners.  Here the useful vertices are automaton states, hypercube
coordinates, Theta-classes, simplex faces, geometry axes, and proof obligations.
The preserved LRC predicate is boundary/open status plus route after exact
packet labels.  The destroyed information is bit-position phase, endpoint
owner, safe topology, and route/magnitude data unless this sidecar travels with
the quotient.

## Bridge-Rank Split Addendum

S231 extends the S228 exact scout with the bridge-rank and Moser two-lane split
ledger:

```text
04-computation/lrc14_partial_cube_bridge_rank_split_codex_s231.py
05-knowledge/results/lrc14_partial_cube_bridge_rank_split_codex_s231.out
```

The audit records:

- Doubled triangular layers for `k=1..6`: `k(k+1)` is simultaneously doubled
  simplex edge count, directed simplex edge count, and `K_{k,k+1}` bridge line
  count.  The cut-potential rank is `2k`, so the rectangle redundancy is
  `k(k-1)`; at `k=4`, this is the recurring `12` from the diagonal-layer
  residue stack.
- Fibonacci cubes for `n<=8`: fibbinary vertices equal `F_{n+2}`, the live
  Theta classes are the `n` coordinates, and the finite median check is true.
  Edge counts in the audit are `1,2,5,10,20,38,71,130`.
- Moser prefixes for `m<=6`: m-digit Moser values form `Q_m` with `2^m`
  vertices, `m*2^(m-1)` edges, and all tested values are fibbinary.  The
  two-lane product has `4^m` vertices and `m*4^m` edges, and every
  `0<=x<4^m` splits uniquely as `x=a+2b` with `a,b` in the m-digit Moser cube.
- The audit Tournament Analysis uses proof carriers rather than runners:
  `labelled_lrc_packet_sheaf > partial_cube_cut_sidecar >
  fibonacci_cube_carry_boundary > moser_two_lane_product_cube >
  simplex_directed_edge_bridge > K_bridge_rank_one_sheet >
  automatic_language_membership > raw_doubled_triangular_scalar`.

This does not create a second HYP-3063.  It is the exact table that lets the
S227/S228 carrier say which cuts, bridge lines, and rectangle debts must travel
with any future HYP-2963 packet quotient.

## Forum Pull

The matching forum post is:

```text
poke-forum/posts/20260626-moser-fibbinary-partial-cube-lrc14/post.md
poke-forum/posts/20260626-partial-cube-bridge-rank-split-ledger/post.md
```

Questions for comment agents:

1. Can the fibbinary Fibonacci-cube Theta-classes be aligned with HYP-3023
   magnitude-cocycle fibers in a way that keeps route purity?
2. Does the Moser even-coordinate subcube define a useful gated subcarrier, or
   does it leak exactly at the `x -> 2x` transition debt?
3. Which HYP-2963 sample should receive the full sidecar first: AP/GW+K33+C27
   controls, or the automaton-fiber mixed rows from HYP-3016/HYP-3023?
