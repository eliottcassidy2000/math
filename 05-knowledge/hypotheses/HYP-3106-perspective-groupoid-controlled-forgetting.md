---
id: HYP-3106
title: Perspective groupoid functors for A000568 and LRC14 controlled forgetting
status: EVIDENCE / exact scout and proof-carrier synthesis; not a proof
source: codex-2026-06-27-S261
tangent: T1183
technique: LTI-244
tournament_technique: LTT-142
script: 04-computation/perspective_groupoid_forgetting_codex_s261.py
result: 05-knowledge/results/perspective_groupoid_forgetting_codex_s261.out
related:
  - HYP-3103
  - HYP-3104
  - HYP-3105
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

# HYP-3106: Perspective Groupoid Functors For A000568 And LRC14 Controlled Forgetting

## Claim

The A000568 perspective ladder should be used as a family of quotient functors,
not as one scalar coincidence.

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

The S261 layer treats node, directed-edge, cycle, clique, conflict,
observer-cut, component-bound, obstruction-cocycle, obstruction-transfer,
dihedral, owner, miss-count-root, and maximizer-transfer perspectives as
functors with:

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
quotient has legal transport through add/delete/observer/reflection/component
bound/chart-gluing/moment-root operations.

## Computation

`04-computation/perspective_groupoid_forgetting_codex_s261.py` tests this on
the first two shifted A000568 failure windows and emits a proof-carrier
tournament over perspective functors.

```text
05-knowledge/results/perspective_groupoid_forgetting_codex_s261.out
```

Exact shifted counts:

```text
m  U(m)  P_node(m)  U(m+1)  defect  node_depths
1     1          1  1        0       [1, 1, 1, 1]
2     1          2  2        0       [2, 2, 2, 2]
3     2          4  4        0       [4, 4, 4, 4]
4     4         12  12       0       [10, 12, 12, 12]
5    12         48  56       8       [36, 47, 48, 48]
6    56        296  456      160     [196, 280, 294, 296]
```

Thus the first prompt failure is still `m=5`, but the next row makes the
lesson stronger: node-depth reaches exact rooted memory at `m=6` by depth `4`,
yet `P_node(6)=296` remains far below `U(7)=456`.

Exact non-node carrier totals:

```text
m  arc_orbits  triple_orbits  transitive_triples  cyclic_triples  conflict_pair_orbits
1           0              0                   0               0                    0
2           1              0                   0               0                    0
3           4              2                   1               1                    0
4          16             12                   8               4                    0
5          88             88                  64              24                    0
6         704            928                 688             240                   32
```

The conflict/Omega carrier first becomes visible at `m=6`, exactly where the
shifted rooted-node gap becomes too large to treat as a one-step curiosity.

The edge-sector and dihedral/converse audit reproduces HYP-3049:

```text
size/internal sector decks: 55/56
cross/full sector decks:   56/56
self-converse U(6):        12
only size/internal collision: masks 344 and 345, a converse pair
```

S261 also folds in the incoming S66 miss-count PGF signal:

```text
consec_8: #real=0/6, extreme_mass=0.3476, L_yK8=3.5823
```

The root-realness / zero-confinement field is another controlled-forgetting
sidecar: a moment scalar is unsafe if the PGF root signature is the coordinate
that separates the extremal sector-correlation regime.

## Perspective-Functor Tournament

S261 uses perspective functors as tournament vertices:

```text
raw_A000568_class
node_k_depth_cache
exact_rooted_node
ordered_pair_extension
directed_edge_sector
cycle_chirality
transitive_clique_insertion
conflict_omega
dihedral_reflection_quotient
normal_fan_component_packet
first_obstruction_cocycle
miss_count_pgf_roots
endpoint_owner_packet
```

Pairwise observable axes:

```text
root, extension, edge, cycle, topology, obstruction, analytic roots,
dihedral, owner.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,6:1,7:3,8:1,10:1,11:1,12:1}
directed_3_cycles=4
scc_sizes=[5,1,1,1,1,1,1,1,1]
hamiltonian_path_count=13
```

One priority path is:

```text
endpoint_owner_packet
-> first_obstruction_cocycle
-> directed_edge_sector
-> ordered_pair_extension
-> cycle_chirality
-> dihedral_reflection_quotient
-> normal_fan_component_packet
-> conflict_omega
-> transitive_clique_insertion
-> miss_count_pgf_roots
-> exact_rooted_node
-> node_k_depth_cache
-> raw_A000568_class
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
intuition, HYP-3054's observer-extension cut calculus, HYP-3100's
tournament-certificate grammar, incoming HYP-3101's component-bound route,
incoming HYP-3102's first-obstruction route, and S66's miss-count PGF root
ledger.  Incoming HYP-3104's maximizer signal atlas supplies the matching
transfer-risk sidecar: a perspective quotient must know whether it is erasing
the first live interaction order, an exchange trap, or an `H`-spectrum alarm.
Incoming HYP-3105's obstruction-transfer atlas adds the terminal-use test:
the quotient must state whether a forbidden spectrum has a faithful transfer
functor or only a destroyed-coordinate analogy.

A perspective functor is legal only when it can survive the next LRC stress
test or name the sidecar that fails.

## Assumption Challenge

Vertices under consideration include runners, tournament nodes, directed
edges, ordered pairs, cycles, cliques, conflict pairs, wall-crossing events,
dihedral orbits, residues, endpoint-owner packets, matrix sidecar columns,
proof obligations, and formal certificate functors.

The selected object for this scout pass is the perspective functor itself.
It preserves the relation between a quotient and the next operation that will
stress it.  It destroys raw labelled runner identity and full extension rows,
unless the sidecar named above carries them.
