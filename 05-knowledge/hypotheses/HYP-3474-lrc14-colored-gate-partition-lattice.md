---
id: HYP-3474
title: LRC14 colored gate partition-lattice
status: EVIDENCE / finite quotient-legality scout; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3471 colored gate-reservoir, complementary to HYP-3472 dead-cover boundary-current and HYP-3473 colored gate formalization, layered over HYP-3453 survivor-gate/component-cover join and the HYP-3462/HYP-3470 AP84 closure stack
tangent: T1434
technique: LTI-434
tournament_technique: LTT-334
script: 04-computation/lrc14_colored_gate_partition_lattice_codex_20260629.py
result: 05-knowledge/results/lrc14_colored_gate_partition_lattice_codex_20260629.out
reflection: 07-reflections/lrc14-colored-gate-partition-lattice-codex-20260629.md
related:
  - HYP-3473
  - HYP-3472
  - HYP-3471
  - HYP-3470
  - HYP-3462
  - HYP-3461
  - HYP-3460
  - HYP-3459
  - HYP-3458
  - HYP-3457
  - HYP-3456
  - HYP-3455
  - HYP-3454
  - HYP-3453
  - HYP-3451
  - HYP-3450
  - HYP-3439
  - HYP-3438
  - HYP-3436
  - HYP-3431
  - HYP-2595
  - HYP-2594
  - HYP-2593
  - THM-523
  - OPEN-Q-108
---

# HYP-3474: LRC14 Colored Gate Partition-Lattice

## Claim

The next useful abstraction after HYP-3471 is not another scalar color count;
it is a partition lattice of legal forgetful quotients.

Given a row-level quotient `Q`, the quotient is legal for a proof predicate
only when every fiber of `Q` is pure for that predicate.  HYP-3474 executes
this Myhill-Nerode-style guardrail on the current HYP-3471 colored-gate bank.

The finite theorem target becomes:

```text
dead-cover obstruction
  -> rank <= 2 E/branch survivor gate
  -> minimum E/branch structural gate sidecar
  -> route flags / AP-vs-random / same-branch / cross-branch dispatch
```

The quotient may forget a coordinate only after fiber purity is proved or after
a named sidecar reconstructs the lost coordinate.

## Exact Readout

Script:

```text
04-computation/lrc14_colored_gate_partition_lattice_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_colored_gate_partition_lattice_codex_20260629.out
```

Axis legend:

```text
K = row endpoint-kind set
N = row numeric mod-14 residue set
T = row typed mod-14 residue set
S = row structural sidecar set
F = row full colored-gate set
C = row low-rank gate count profile
M = minimum E/branch structural gate
A = canonical AP84 palette presence bit
```

Bank summary:

```text
rows=135
dead_rows=130
theorem_failures=0
ap84_packet_rows=68
same_branch_rows=102
cross_branch_rows=54
only_e_branch_rows=31
```

The audited theorem-failure bit is identically false:

```text
theorem_failure_bit = dead and not has_e_branch = False on 135/135 rows.
```

This means every quotient preserves the implication as a bank fact.  It does
not mean every quotient preserves the routing data needed for a proof.  For
that, route flags and minimum-gate families are the relevant targets.

Singleton quotient mixing:

```text
axis=K fibers=8   mixed_route_rows=132
axis=N fibers=123 mixed_route_rows=0   mixed_dead_min_struct_rows=12
axis=T fibers=123 mixed_route_rows=0   mixed_dead_min_struct_rows=12
axis=S fibers=107 mixed_route_rows=32  mixed_dead_min_struct_rows=32
axis=F fibers=123 mixed_route_rows=0   mixed_dead_min_struct_rows=12
axis=C fibers=122 mixed_route_rows=0   mixed_dead_min_struct_rows=15
axis=M fibers=17  mixed_route_rows=124 mixed_dead_min_struct_rows=85
axis=A fibers=2   mixed_route_rows=135 mixed_dead_min_struct_rows=135
```

The route flag distribution is:

```text
(dead, ap84, same, cross, only_e)
(True, True, True, True, False):   30
(True, False, True, False, False): 28
(True, False, True, True, False):  24
(True, True, True, False, False):  20
(True, True, False, False, True):  17
(True, False, False, False, True): 11
(False, False, False, False, True): 2
(False, False, False, False, False): 2
(False, True, False, False, True): 1
```

The coarsest route-pure quotients in the declared axis family are:

```text
route_flags_pure axes=('C',) fibers=122 max_fiber=10
route_flags_pure axes=('F',) fibers=123 max_fiber=10
route_flags_pure axes=('N',) fibers=123 max_fiber=10
route_flags_pure axes=('T',) fibers=123 max_fiber=10
```

This is a useful warning, not a proof shortcut.  The duplicate fibers explain
why the count-profile route purity is believable as finite telemetry: the
large duplicate classes are AP-tail copies, plus one nondead random pair.

```text
C duplicate size=10:
  covering_AP_with_84, ap_omit_12_tail_84x01, ..., ap_omit_12_tail_84x12
  route flags=(True, True, False, False, True)
  min_struct=(B1|E, branch1, left_bad_edge, 1, 1)

C duplicate size=2:
  random_covering_044, random_covering_053
  route flags=(False, False, False, False, False)
```

For the dead minimum gate family, the smallest pure two-axis packets are:

```text
dead_min_struct_pure axes=('C','M') fibers=125 max_fiber=7
dead_min_struct_pure axes=('F','M') fibers=125 max_fiber=7
dead_min_struct_pure axes=('N','M') fibers=125 max_fiber=7
dead_min_struct_pure axes=('T','M') fibers=125 max_fiber=7
```

The same two-axis packets are pure for the more detailed dead minimum full
family `(minimum typed word, minimum structural word)` on this bank.  That is
the strongest new lead: the minimum E/branch structural sidecar `M`, although
far too coarse by itself, becomes a legal obstruction-family classifier once
paired with a high-cardinality row shadow.

## Interpretation

HYP-3471 supplied the local colored gate.  HYP-3474 supplies the legal
forgetting test for using that gate.

Three levels now separate cleanly:

```text
bank fact:
  dead -> E/branch gate has no failures in the audited rows

route quotient:
  C, N, T, F are pure for route flags in this bank

gate-family quotient:
  M must be paired with C/N/T/F to be pure for the minimum-gate family
```

Endpoint kind `K` is excellent for proving same/cross/only-E branch bits as
coarse shape data, but it mixes AP-vs-random route flags across `132/135`
rows.  The AP84 bit `A` separates AP packet presence and nothing else.  The
minimum structural gate `M` is mathematically meaningful but not a legal route
quotient alone.

The tempting surprise is the count profile `C`.  On the current bank it is
route-pure and wins the tournament score tie.  The disciplined conclusion is:
`C` is a finite compression diagnostic that may hint at a hidden gate-count
recursion, but it is not yet a theorem quotient because it has no intrinsic
reconstruction map to interval geometry, owner currents, or component
adjacency.

## Proof Pull

The next rigorous finite lemma should be formulated as a labelled packet
theorem:

```text
For every primitive LRC14 cover-row in the target class:
  if a dead-cover component exists,
  then there is a rank <= 2 E/branch survivor gate G,
  and the quotient (row-shadow, min-struct(G)) is pure for the route class
  or emits one of the named debts:
    AP84 packet, same-branch gluing, cross-branch gluing,
    endpoint-owner current, two-adic floor descent, signed-SPEC/Rprime.
```

This is the finite analogue of the broader "guardrails for what a quotient is
allowed to forget" principle that has appeared across irreducibility,
unital-design, Faulhaber, Pollock, unit-distance, and tournament-tiling
threads.  Here the guardrail is executable.

## Tournament Analysis

Vertices are quotient carriers and sidecar selections, not runners or arcs.
The pairwise observable is:

```text
route-predicate purity
+ mixed-fiber penalty
+ compression bonus
+ max-fiber penalty
```

Switch: higher aggregate proof-facing score; ties use carrier name.

Fingerprint:

```text
score_hist={53:1,83:1,106:1,148:1,184:1,208:3,209:4}
directed_3cycles=0
hamiltonian_path=
  count_profile
  -> full_color_set
  -> numeric_mod14_set
  -> typed_mod14_set
  -> all_color_sidecars
  -> full_plus_min_struct
  -> typed_plus_min_struct
  -> structural_plus_ap_bit
  -> structural_set
  -> endpoint_kind_set
  -> ap84_presence_bit
  -> min_e_branch_struct
```

Assumption challenge: runners, arcs, raw colors, endpoint kinds, numeric
residues, typed residues, structural sidecars, row-level color sets, minimum
E/branch gates, cover components, and proof obligations were considered.  The
chosen carrier is row-level quotient fibers in the partition lattice generated
by colored-gate axes.  It preserves target-purity of proof-facing row
predicates.  It destroys exact interval geometry, owner-current routing,
component adjacency, and row identity unless a color/sidecar axis reconstructs
that information.
