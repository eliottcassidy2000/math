---
id: HYP-3055
title: Duodecimal observer-extension cut payload at the first A000568 defect
status: SYNTHESIS / exact count ledger and cross-carrier abstraction; not a proof
source: codex-2026-06-26-S219
tangent: T1137
script: 04-computation/duodecimal_observer_extension_payload_codex_s219.py
result: 05-knowledge/results/duodecimal_observer_extension_payload_codex_s219.out
related:
  - HYP-3054
  - HYP-3053
  - HYP-3052
  - HYP-3051
  - HYP-3050
  - HYP-3049
  - HYP-3048
  - HYP-3047
  - HYP-3043
  - HYP-3039
  - HYP-3031
  - HYP-2991
  - HYP-2989
  - HYP-2928
  - HYP-2120
  - HYP-2121
  - THM-381
  - THM-385
  - OPEN-Q-108
---

# HYP-3055: Duodecimal observer-extension cut payload at the first A000568 defect

## Claim

At the first shifted rooted-perspective failure,

```text
P(5)=48,   U(6)=56,   U(5)=12,   U(4)=4,
```

the recurring `12` is real but it is not the additive defect.  The arithmetic
guardrail is

```text
48 + 12 = 60,
```

so the exact first-failure identity is instead

```text
U(6) = P(5) + U(5) - U(4) = 48 + 12 - 4 = 56.
```

Equivalently,

```text
U(6)-P(5) = 8 = U(5)-U(4),
```

and in denominator-14 form,

```text
P(5)/U(6) + U(5)/U(6) - U(4)/U(6)
  = 12/14 + 3/14 - 1/14
  = 1.
```

The dozen is therefore a control/fold slice, while the additive defect is the
dozen after quotienting out the four-class overlap.

## Dozen bridge

The same `12` appears across the three neighboring sizes:

```text
P(4) = U(5) = SC(6) = 12.
```

Here `P(m)` is the rooted/node-perspective count on `m` vertices, `U(n)` is
A000568, and `SC(n)` is the self-converse class count.  At the first defect:

```text
P(5) = 5*U(5) - SC(6) = 60 - 12 = 48,
U(6) = 5*U(5) - U(4)  = 60 - 4  = 56,
U(6)-P(5) = SC(6)-U(4) = 12 - 4 = 8.
```

This should not be promoted to a recurrence.  The next row already rejects it:

```text
U(7)-P(6)=160,   while   U(6)-U(5)=44.
```

The status is local exact structure at the first controlled-forgetting
failure, not a new A000568 recurrence.

## Observer-extension interpretation

HYP-3047 through HYP-3053 say that the rooted node cache has become too small.
The missing coordinate is an observer-extension/cut payload:

```text
parent class
  + incident word modulo Aut(parent)
  + endpoint role / ordered-pair sector data
  + deletion-parent fiber profile
  + rectangle/hourglass consistency residues
  -> unrooted child sink.
```

The exact S213 edge-sector lift verifies that a rooted 5-perspective plus the
new observer's incident word gives the ordered-pair 6-perspective count
`1408=1408`; forgetting endpoint role gives `704` directed-edge/unordered-pair
perspectives.  Sector sizes and internal sector decks separate only `55/56`
six-classes; adding cross-sector orientation separates `56/56`.  Thus
`cross_sector_orientation_word` remains the first compact observability column
after the rooted cache.

The S216 deletion-fiber table gives the sink side.  At `5 -> 6`, parent
automorphism word orbits are `296`, matching rooted 6-perspectives, while the
unrooted sink has `56` classes.  The average rooted child presentations per
unrooted class is `296/56=37/7`, and unique deletion-parent class count
averages `215/56`.  The quotient is not source class plus a scalar word count;
the deletion fiber is part of the carrier.

The S217 diagonal-layer law supplies the geometric sidecar.  `K_{3,4}` is the
first dozen line block:

```text
12 lines = 6 independent cut-space rank + 6 rectangle residues.
```

Globally, full adjacent-layer redundancy splits as local rectangles plus
hourglass cycles:

```text
full redundancy = 2*C(n-1,3) + C(n-2,2).
```

Fixed-path redundancy has the analogous smaller law:

```text
fixed-path redundancy = 2*C(n-2,3) + C(n-3,2).
```

Zero residues descend to layer potentials.  Nonzero residues name the hidden
coordinate, just as cross-sector orientation names the first edge-sector
coordinate.

## Cross-project abstract

The same proof-interface pattern should be used outside A000568:

- LRC endpoint-owner packets: retain observer source, owner strip, route/status
  payload, and magnitude before quotienting an apex or residue tournament.
- Haar/fixed-margin rectangles: the 2-by-2 zeta residue is the same kind of
  rectangle defect as a `K_{k,k+1}` cycle-space residue.
- Pair-good decoys and residual capacitors: raw counts are scouts; blocker
  teeth, edge tail/tip sectors, and residual currents are the theorem-facing
  keys.
- Matrix sidecar observability: rows are coarse-fiber packet pairs and columns
  are hidden coordinates; a quotient is proof-safe only when every
  route/status-changing pair is separated, reconstructed, annihilated,
  descended, or routed to named residual debt.
- The LRC tournament spectrum: a single apex iso class forgets magnitude; the
  set of phase-swept winding classes plus binding scale is the observer-cut
  payload.

## Assumption challenge

This pass explicitly considered runners, arcs, gaps, fixed diagonal sections,
section boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
matroid circuits, and proof obligations.  The chosen tournament-analysis
vertices are proof carriers:

```text
endpoint_owner_packet_sheaf
observer_extension_cut_payload
cross_sector_orientation_word
deletion_parent_fiber_profile
rectangle_hourglass_cycle_residue
incident_word_orbit_under_aut
ordered_pair_edge_sector_deck
rooted_node_perspective_cache
fixed_path_half_tiling_shadow
raw_A000568_class_count
raw_labelled_word_count
```

The preserved predicate is observer/cut/cycle payload available for LRC-style
fiber safety.  The destroyed data is labelled identity, endpoint role, deletion
parent, rectangle/hourglass residue, and magnitude unless retained as named
sidecars.

## Evidence

Script:

```text
04-computation/duodecimal_observer_extension_payload_codex_s219.py
```

Stored output:

```text
05-knowledge/results/duodecimal_observer_extension_payload_codex_s219.out
```

The carrier tournament is transitive:

```text
endpoint_owner_packet_sheaf >
observer_extension_cut_payload >
cross_sector_orientation_word >
deletion_parent_fiber_profile >
rectangle_hourglass_cycle_residue >
incident_word_orbit_under_aut >
ordered_pair_edge_sector_deck >
rooted_node_perspective_cache >
fixed_path_half_tiling_shadow >
raw_A000568_class_count >
raw_labelled_word_count.
```

The point is not that `12` is a scalar answer.  The point is that the dozen
keeps reappearing exactly where controlled forgetting wants a fold/control
slice: rooted 4-perspectives, 5-classes, self-converse 6-classes, source/sink
extension slices, and the first dozen inter-layer line block.
