---
id: HYP-3067
title: LRC14 Desargues-median finalization lens
status: SYNTHESIS / graph-theoretic proof-interface target; not a proof
source: codex-2026-06-26-S233
tangent: T1149
script: 04-computation/lrc14_desargues_median_lens_codex_s233.py
result: 05-knowledge/results/lrc14_desargues_median_lens_codex_s233.out
related:
  - HYP-3057
  - HYP-3056
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
  - HYP-3037
  - HYP-3034
  - HYP-3031
  - HYP-3024
  - HYP-3018
  - HYP-2997
  - HYP-2963
  - HYP-2314
  - THM-572
  - OPEN-Q-108
---

# HYP-3067: LRC14 Desargues-Median Finalization Lens

## Claim

The final LRC14 proof should be tested as a medianization problem on the
controlled-forgetting proof graph.

Build a graph whose vertices are proof states inside a coarse LRC packet
fiber:

```text
packet state
route/certificate state
observer-cut payload orbit
sidecar column state
named residual debt state
```

Edges change one retained sidecar or one discharge decision.  After attaching
the correct sidecars, every three proof routes inside one coarse fiber should
have a unique median consensus state:

```text
I(A,B) cap I(B,C) cap I(C,A) = {unique proof state}.
```

This is not a reduction of LRC14 to median graph theory.  It is a proof
interface: if the sidecar graph is median-like, controlled forgetting has a
unique compatibility center.  If a triple has no center, the quotient still
contains a Desargues-type incidence twist.  If a triple has multiple centers,
the quotient retained redundant or ambiguous coordinates.

## Desargues Warning

The Desargues graph is a compact bipartite incidence graph, but it is not
median.  The S233 script verifies:

```text
Desargues graph:
  n=20, m=30, cubic, bipartite
  diameter=5, girth=6
  cycle lengths <=16: 6,8,10,12,14,16
  median=False
  median failures=160
  median intersection size histogram: {0:160, 1:1380}
  theta sidecar sketch: 5 classes of size 6
```

Comparison controls:

```text
Q4 hypercube: median=True, theta classes 4 of size 8
4x4 grid:    median=True, theta classes 6 of size 4
C6 cycle:    median=False, two empty-center triples
```

The warning for LRC14 is precise: bipartite incidence, regularity, theta-like
edge classes, and rich even-cycle structure do not by themselves provide a
safe proof quotient.  A residual graph can look organized and still have
empty-center triples where three certificate routes cannot agree.

## LRC Translation

The candidate proof graph is not a graph of runners.  Vertices should be
proof states and sidecar payloads:

```text
qdiv/direct witness gate
exact M / Farey branch
closed arc-H1 owner support
primitive safe deck
Henselian ET+unit status gate
residual capacitor cut
Haar zeta / rectangle defect
endpoint-owner strip current
observer-cut payload orbit
value-origin type
Fejer/Ramanujan/smoothing certificate state
named THM-572/F7 debt
```

The median test asks whether three routes, such as

```text
topology route, owner route, period route
Fejer route, Haar route, Ramanujan route
q-witness route, covering route, K33/petal route
observer-extension route, deletion-fiber route, rectangle-residue route
```

have a unique packet state after all legal sidecars are retained.

Empty median intersection means the current quotient forgot a Desargues-like
cut payload.  The repair is not a scalar count.  It is one of:

```text
retain a sidecar,
reconstruct it,
make it exact as a boundary/cocycle/cut,
dual-annihilate it,
descend it familywise,
stop at AP/GW boundary topology,
or name residual debt.
```

## Controlled-Forgetting Rule

HYP-3054 says to name the next operation and its payload.  HYP-3056 says to
record the payload orbit under visible automorphisms.  HYP-3057 says to tag
the value origin before using a small integer.  HYP-3067 adds the
compatibility test:

```text
For every coarse fiber and every triple of proof routes,
the retained sidecar graph must have exactly one median center,
or the missing center must be routed to a named Desargues defect.
```

This makes "finalizing LRC14" a finite-looking task over the HYP-2963 packet
bank: build the sidecar graph, compute median failures, and classify every
failure by the first missing sidecar.

## Desargues Defect Schema

A Desargues defect is a triple:

```text
(route_A, route_B, route_C)
```

inside one coarse fiber such that

```text
I(A,B) cap I(B,C) cap I(C,A) = empty.
```

Candidate sources in the current repo:

```text
compact topology vs exact M vs endpoint owner,
primitive deck vs AP-tail q13 clock vs Haar zeta,
automaton terminal word vs magnitude cocycle vs owner support,
observer incident orbit vs deletion-parent profile vs rectangle residue,
pair-good blocker tooth vs active barcode owner vs normal-fan support.
```

The defect is discharged only when the empty center becomes a unique center
after a named sidecar is attached, or when the defect is shown to be an
AP/GW boundary atom or THM-572/F7 residual sector.

## Tournament Analysis

Vertices are proof-graph objects, not runners:

```text
median_center
sidecar_hyperplane
observer_cut_orbit
value_origin_type
endpoint_owner_strip
primitive_period_deck
arc_H1_owner_support
residual_capacitor_cut
rectangle_hourglass_residue
Desargues_defect
raw_incidence_graph
raw_scalar_count
```

Pairwise observable:

```text
route/status separation,
triple-median existence,
triple-median uniqueness,
sidecar reconstructibility,
boundary/cocycle exactness,
dual-annihilation availability,
family descent,
named residual debt,
proof cost.
```

Switch/gauge:

```text
orient toward the carrier that turns the most route triples into unique
median centers with the fewest new residual debts.
```

Tie Hamiltonian path:

```text
median_center
> sidecar_hyperplane
> observer_cut_orbit
> value_origin_type
> endpoint_owner_strip
> primitive_period_deck
> arc_H1_owner_support
> residual_capacitor_cut
> rectangle_hourglass_residue
> Desargues_defect
> raw_incidence_graph
> raw_scalar_count
```

Assumption challenge: the vertex set is not runners, arcs, or raw tournament
classes.  It is route/certificate sidecar states inside a proof graph.  The
preserved LRC predicate is boundary/open status plus route and certificate
schedulability.  The destroyed information is the median center of three
proof routes when a quotient forgets owner, period, topology, cut, or value
origin payload.

## Next Pulls

1. Build a HYP-2963 sidecar graph with vertices `(packet, sidecar state,
   discharge state)` and edges changing one legal sidecar.
2. For every mixed route/status coarse fiber, compute median-triple failures
   among topology, owner, period, zeta, automaton, and observer-cut routes.
3. Classify each failure as fixed by a known sidecar or as named Desargues
   defect debt.
4. Add fields `proof_graph_vertex`, `sidecar_hyperplane_id`,
   `route_triple_id`, `median_center_status`, `desargues_defect_id`, and
   `medianization_exit` to the packet manifest.
5. Use forum posts to ask agents to attack one route triple at a time rather
   than proposing more scalar summaries.
