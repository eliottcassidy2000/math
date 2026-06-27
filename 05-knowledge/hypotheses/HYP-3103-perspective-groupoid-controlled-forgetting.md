---
id: HYP-3103
title: Perspective groupoid functors for A000568 and LRC14 controlled forgetting
status: RESERVED / synthesis and computation in progress; not a proof
source: codex-2026-06-27-S261
tangent: T1181
technique: LTI-242
tournament_technique: LTT-140
script: 04-computation/perspective_groupoid_forgetting_codex_s261.py
result: 05-knowledge/results/perspective_groupoid_forgetting_codex_s261.out
related:
  - HYP-3102
  - HYP-3101
  - HYP-3100
  - HYP-3057
  - HYP-3054
  - HYP-3050
  - HYP-3049
  - HYP-3048
  - HYP-3047
  - HYP-2121
  - HYP-2120
  - HYP-2087
  - HYP-2096
  - THM-381
  - THM-385
  - THM-573
  - OPEN-Q-108
---

# HYP-3103: Perspective Groupoid Functors For A000568 And LRC14 Controlled Forgetting

## Reserved Claim

This session reserves HYP-3103 for a formal extension of the A000568
perspective ladder into a small "perspective groupoid" calculus.

Known input from HYP-3047/HYP-3050:

```text
P(5)=48, U(6)=56, defect=8.
k-depth node memory reaches exact rooted node type by depth 3 at m=5.
```

Known input from HYP-3049/HYP-3057:

```text
rooted 5-perspective + incident word = ordered-pair perspective on U(6),
directed-edge sector decks separate 55/56 until cross-sector orientation is kept,
and small integers such as 12, 48, 56 must be typed by quotient origin.
```

The new layer is to treat node, directed-edge, cycle, clique, conflict, and
observer-cut perspectives as functors with:

```text
root object
acted-on automorphism group
depth rule
forgotten coordinate
sidecar needed by the next operation
dihedral/converse action, when present
LRC predicate preserved
```

The point is not to make a deeper node color.  The point is to specify which
quotient has legal transport through add/delete/observer/reflection operations.

## Pending Computation

`04-computation/perspective_groupoid_forgetting_codex_s261.py` will test this
on the first A000568 failure window and emit a proof-carrier tournament over
perspective functors.  The planned scout will record:

```text
node k-depth counts
directed-edge dual sector signatures
cycle and transitive-clique root carriers
conflict/Omega-style carriers once visible
dihedral/converse losses
controlled-forgetting sidecar obligations
```

## LRC14 Target

The intended LRC14 transfer is a typed observer-cut field for packet ledgers:

```text
perspective_functor
root_object
automorphism_action
dihedral_reflection_status
observer_endpoint_role
incident_word_orbit
cross_sector_orientation_word
cycle_chirality_payload
clique_insertion_cut
conflict_omega_payload
legal_discharge
```

This is the formal bridge between the prompt's k-depth / edge-perspective
intuition, HYP-3054's observer-extension cut calculus, and HYP-3100's
tournament-certificate grammar.

After rebasing over incoming HYP-3101 and HYP-3102, this lane is also a
bridge to the component-bound and first-obstruction proof routes: a perspective
functor is legal only when it can survive those two LRC stress tests or name
the sidecar that fails.

## Assumption Challenge

Vertices under consideration include runners, tournament nodes, directed
edges, ordered pairs, cycles, cliques, conflict pairs, wall-crossing events,
dihedral orbits, residues, endpoint-owner packets, matrix sidecar columns,
proof obligations, and formal certificate functors.

The selected object for this reserved pass is the perspective functor itself.
It preserves the relation between a quotient and the next operation that will
stress it.  It destroys raw labelled runner identity and full extension rows,
unless the sidecar named above carries them.
