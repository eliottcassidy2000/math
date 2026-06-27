---
id: HYP-3057
title: Tournament value-origin ledger for 12, 48, and 56
status: EVIDENCE / exact small tournament ledger and quotient guardrail; not a proof
source: codex-2026-06-26-S221
tangent: T1139
script: 04-computation/tournament_value_origin_ledger_codex_s221.py
result: 05-knowledge/results/tournament_value_origin_ledger_codex_s221.out
related:
  - HYP-3054
  - HYP-3053
  - HYP-3052
  - HYP-3051
  - HYP-3050
  - HYP-3049
  - HYP-3048
  - HYP-3047
  - HYP-3039
  - HYP-2991
  - HYP-2989
  - THM-381
  - THM-385
  - OPEN-Q-108
---

# HYP-3057: Tournament value-origin ledger for 12, 48, and 56

## Claim

The small tournament numbers around the first A000568/rooted-perspective
failure should be typed by origin before they are used as proof carriers.

The arithmetic correction is:

```text
R(5) = 48
U(6) = 56
defect = 8
```

So the first shifted failure is `48 + 8 = 56`, not `48 + 12 = 56`.
The value `12` is still structural, but it is a diagonal alignment across
different quotient origins:

```text
R(4)  = 12  rooted/node perspectives on 4 vertices
U(5)  = 12  unlabelled 5-tournament classes
SC(6) = 12  self-converse 6-tournament classes
```

This is the core guardrail.  A small integer by itself is not a theorem object.
It must be tagged as one of:

```text
class count
rooted perspective count
self-converse fixed branch
incident-word orbit count
deletion-fiber multiplicity
ordered-pair / edge-sector sidecar
fixed-path Hamiltonian presentation fiber
rectangle/hourglass cycle-space residue
```

The first failure then becomes a clean transport statement:

```text
parent class + incident word orbit under Aut(parent)
  -> rooted child
  -> unrooted child sink.
```

At `5 -> 6`, this ladder has:

```text
U(5)                         = 12
raw incident extensions       = 12*2^5 = 384
parent-Aut word orbits        = R(6) = 296
unrooted sinks reached        = U(6) = 56
rooted 5-perspective + word   = ordered-pair perspectives = 1408
directed-edge perspectives    = unordered-pair perspectives = 704
```

The defect is not deeper node memory.  It is incident/cross-coupling payload
plus the final unrooting/deletion-fiber quotient.

## Evidence

The S221 script recomputes and joins the S211, S213, S216, and S217 ledgers,
then slots the result after HYP-3054's observer-extension cut calculus: first
name the outside operation/cut payload, then tag which quotient origin produced
the count used to summarize it.
The shift table begins:

```text
m  U(m)  R(m)  U(m+1)  U-R  R/U(m+1)  SC(m)  SC(m+1)
1     1     1       1    0  1/1            1        1
2     1     2       2    0  1/1            1        2
3     2     4       4    0  1/1            2        2
4     4    12      12    0  1/1            2        8
5    12    48      56    8  6/7            8       12
6    56   296     456  160  37/57          12        -
```

The equality `defect = SC(5) = 8` at the first break is only a small-boundary
coincidence; it already fails one step later.  The reusable invariant is not
`defect = self-converse count`, but the broader type ledger: each count has a
source quotient.

The edge-sector repair from HYP-3049 remains the most compact first repair:

```text
sector sizes        separate 55/56 classes
sector internals    separate 55/56 classes
cross orientations  separate 56/56 classes
```

The only size/internal collision is the converse pair with masks `344` and
`345`.  Thus `cross_sector_orientation_word` is the first separating column
after rooted-node cache and incident-word lift.

The S217 fixed-path flow gives the analogous residue law:

```text
free_tiles(n)  = C(n-1,2)
fixed_lines(n) = 2*C(n-1,3)
fixed_rank(n)  = C(n-1,2)-1
fixed_red(n)   = 2*C(n-2,3)+C(n-3,2)
```

Here the redundancy decomposes into local rectangle cycles plus hourglass
cycles linking adjacent bridges.  These are not new classes.  They are the
cycle-space coordinates that explain why a fixed-path presentation duplicates
unlabelled tournament classes.

## LRC Translation

For LRC and LRC14, the lesson is a quotient-admissibility rule.

When a scalar, class count, metagraph degree, automaton state, or tournament
shadow produces a compelling integer, ask first where it came from:

```text
Is it a class count?
Is it a rooted perspective?
Is it a fixed branch of an involution?
Is it an incident-word orbit?
Is it a deletion fiber?
Is it an edge-sector orientation sidecar?
Is it a rectangle/hourglass residue?
```

Only after this origin tag is known can the quotient be used safely.  Otherwise
the same visible number may be hiding different lost coordinates.

This extends the HYP-2990 controlled-kernel rule into the tournament-count
surface:

```text
forgotten coordinate must be
  retained,
  reconstructed,
  annihilated by a dual/cocycle,
  descended to a smaller family,
  or emitted as named residual debt.
```

The `12` alignment is exactly the warning example.  It is real, but it is not
one object.

## Tournament Analysis

Vertices are value-origin carriers and proof obligations, not runners.

Candidate vertices:

```text
endpoint_owner_packet
edge_sector_cross_orientation
deletion_parent_fiber
incident_word_orbit
self_converse_fixed_branch
rooted_perspective_count
rectangle_hourglass_residue
fixed_path_hamiltonian_fiber
raw_unlabelled_class_count
raw_integer_coincidence
```

Pairwise observable: retained LRC predicate payload, quotient-origin clarity,
ability to reconstruct lost coordinates, ability to route residual debt, and
proof cost.

Gauge: orient toward the carrier that preserves more predicate-bearing
coordinates with less hidden debt; ties follow the displayed order above.

Expected high-retention path:

```text
endpoint_owner_packet
> edge_sector_cross_orientation
> deletion_parent_fiber
> incident_word_orbit
> self_converse_fixed_branch
> rooted_perspective_count
> rectangle_hourglass_residue
> fixed_path_hamiltonian_fiber
> raw_unlabelled_class_count
> raw_integer_coincidence
```

The key assumption challenged here is that a count coincidence should identify
one underlying structure.  The safer view is recursive and typed: the count is
only a shadow of its quotient origin.

## Next Pull

Build a `value_origin_type` sidecar for tournament and LRC packet experiments:

```text
value_origin_type
parent_class
incident_word_orbit
root_orbit_count
child_sink_class
deletion_parent_profile
self_converse_status
fixed_path_fiber
edge_sector_cross_orientation_word
rectangle_residue_class
hourglass_residue_class
lost_coordinate_exit
```

Then test whether route/status-mixed LRC packet fibers split faster when their
small numerical shadows are typed by origin before any scalar comparison.
