---
id: HYP-3461
title: LRC14 colored-extension gate carrier
status: SYNTHESIS / executable carrier scout; not an LRC14 proof
source: codex-2026-06-29 connection pass over colorings and extension work
tangent: T1421
technique: LTI-421
tournament_technique: LTT-321
script: 04-computation/lrc14_colored_extension_gate_carrier_codex_20260629.py
result: 05-knowledge/results/lrc14_colored_extension_gate_carrier_codex_20260629.out
reflection: 07-reflections/lrc14-colored-extension-gate-carrier-codex-20260629.md
related:
  - HYP-3460
  - HYP-3459
  - HYP-3458
  - HYP-3457
  - HYP-3456
  - HYP-3455
  - HYP-3454
  - HYP-3453
  - HYP-3439
  - HYP-3438
  - HYP-3436
  - HYP-3425
  - HYP-3134
  - HYP-3133
  - HYP-3056
  - HYP-3054
  - HYP-2595
  - HYP-2594
  - HYP-2593
  - HYP-2250
  - HYP-2247
  - OPEN-Q-108
---

# HYP-3461: LRC14 Colored-Extension Gate Carrier

## Claim

The earlier coloring and extension threads combine into a proof-facing carrier
for the current LRC14 covering floor:

```text
color data      = boundary charge, not a label histogram
extension data  = gate-gluing orbit, not an extension count
covering floor  = colored extension obstruction ledger
```

HYP-2594/HYP-2595 already warned that raw component count `K` overcharges the
colored CRT placement problem; most boundary events cancel unless they survive
the mod-14 resonance condition.  HYP-3133/HYP-3134 and HYP-3056 warned that raw
extension shadows are only legal after the observer-cut payload or paired child
deck is retained or discharged.  HYP-3461 ports both guardrails into the
HYP-3438/HYP-3455/HYP-3456/HYP-3457 gate frontier.

The proposed carrier is a colored gate-extension orbit with fields:

```text
gate word
+ branch mask
+ endpoint walls
+ parent even walls
+ minimal B0/B1 owner covers
+ cover-delta vector
+ mirror / row automorphism orbit
+ low-rank escape availability
+ discharge mode or named residual debt
```

The raw word `B-S-B-S-B`, a raw A000568 middle shadow, or a raw component count
is telemetry only until those fields are either reconstructed, exact,
dual-annihilated, descended, boundary-stopped, or named as debt.

## Exact Readout

Executable scout:

```text
04-computation/lrc14_colored_extension_gate_carrier_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_colored_extension_gate_carrier_codex_20260629.out
```

The tournament over proof carriers is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1,1,1]
```

Selected Hamiltonian path:

```text
random031_mirror_colored_extension_clause
-> AP84_colored_endpoint_floor_packet
-> survivor_gate_colored_payload
-> observer_cut_orbit_ledger
-> minimal_two_color_bad_core
-> colored_resonance_half_boundary
-> A000568_edge_envelope_controlled_forgetting
-> PH_bad_coloring_outer_extension_rank
-> metagraph_GF2_color_boundary
-> phase_color_reservoir
-> raw_gate_word_shadow
-> raw_component_count
```

The top two carriers are the two current covering-floor sidecars, now with
incoming HYP-3458 supplying the AP84 coloring-recursion refinement, incoming
HYP-3459 supplying the AP84 coloring-discrepancy legality layer, and incoming
HYP-3460 supplying the phase-branch pullback:

1. `random031_mirror_colored_extension_clause` from HYP-3455.  It keeps the
   seven-owner mirror gate pair, branch masks, owner deltas, and the `94`
   low-rank escapes.
2. `AP84_colored_endpoint_floor_packet` from
   HYP-3454/HYP-3456/HYP-3457/HYP-3458/HYP-3459/HYP-3460.  It keeps the fixed low
   corridors, endpoint labels, moving high-grid count, mod-35 clock,
   endpoint-rank subcolor, color-packet legality, phase-branch pullback, and
   four transient survivor windows.

Thus the color/extension merge does not replace the current frontier; it gives
the common language for the two most live floor packets.

## Prior Work Reconnected

HYP-2594/HYP-2595 provide the coloring half.  A color is not merely a residue
name; it is the residue-compatible boundary charge that survives the
half-boundary Fourier identity.  In gate language, branch masks, endpoint wall
labels, and owner deltas should be promoted to a colored boundary vector before
charging an obstruction to raw component count.

HYP-2247 provides the extension-rank warning from finite coloring recursion.
A side-choice coloring can look viable until the child-extension profile is
checked.  In LRC14 this says a survivor gate word or branch color is incomplete
without its coherent bad-cover children and extension rank.

HYP-2250 supplies the GF(2) boundary analogy: the real invariant is a color
layer as a boundary chain, not a visible blue/black label.  The LRC14 analogue
is the owner-delta / branch-mask boundary vector around a gate.

HYP-3133/HYP-3134 supply the A000568 extension guardrail.  The middle shadow is
useful only between a lower sector word and an upper paired child deck.  In
HYP-3438 language, the raw gate word is the middle shadow; endpoint labels,
owner covers, and low-rank escape data are the paired child deck.

HYP-3056 supplies the quotient discipline: every forgetting step needs a
visible quotient, a next observer, a cut-payload orbit, a changed LRC predicate,
and an explicit discharge mode.

## Proof-Facing Obligations

HYP-3461 leaves six concrete obligations.

```text
O1_color_boundary_vector
  Promote branch masks, endpoint-wall labels, and owner deltas to a colored
  boundary vector.  HYP-2595 says raw component K is the wrong charge.

O2_gate_extension_orbit
  Define the orbit modulo row/mirror automorphisms with gate word, branch mask,
  endpoint walls, owner covers, cover deltas, and low-rank escape status.

O3_random031_finite_clause
  Prove the HYP-3455 mirror B-S-B-S-B pair on owner union
  (23,45,93,113,147,169,173) cannot glue to a full cover unless one of the
  named exits fires.

O4_AP84_floor_splice
  Treat HYP-3456, HYP-3457, HYP-3458, HYP-3459, and HYP-3460 as closed endpoint/color
  sidecars.  The remaining AP84 task is the HYP-3431 fixed-corridor identity
  plus the splice into HYP-3439.

O5_controlled_forgetting
  Before forgetting to a raw word, raw A000568 shadow, or raw count, record the
  observer-cut orbit and discharge mode from HYP-3056.

O6_finite_colored_placement
  If HYP-2595 is formalized, combine the O(k+c_GP) deficit target with the
  HYP-2593 Sigma floor and exact V<711 checking; then feed the surviving color
  residue into the gate-extension ledger.
```

## Proof Pull

The strongest immediate theorem target is not a broad new density lemma.  It is
the finite colored-extension gluing clause for HYP-3455:

```text
the two max-delta mirror gates of random_covering_031 cannot be extended into
a full two-color cover without exposing a rank<=2 escape, endpoint-spine route,
owner-current imbalance, two-adic descent exit, or signed-SPEC/Rprime debt.
```

The AP84 side should be spliced in parallel:

```text
HYP-3431 fixed corridors
+ HYP-3454 rank-one endpoint phase for m>=5
+ HYP-3456 mod-35 floor count
+ HYP-3457 finite transient packet
+ HYP-3458 AP84 coloring-recursion state
+ HYP-3459 AP84 coloring-discrepancy legality layer
-> HYP-3439 AP-tail bridge
```

After those two sidecars are named in one colored-extension ledger, HYP-3453's
rank-2 survivor-gate transversal becomes a more plausible global cover-floor
router: any dead/mixed gate extension either lands in a named AP84 or random031
packet, or exposes a low-rank escape.

## Tournament Analysis

Vertices are proof carriers and quotient shadows, not runners, arcs, or raw
component counts.  The pairwise observable is:

```text
predicate retention
+ color boundary payload
+ extension/gluing payload
+ endpoint-owner retention
+ finite auditability
+ current-frontier fit
+ scalar guardrail
+ discharge specificity
```

The switch chooses the higher weighted carrier score, with ties broken by the
declared proof-facing order.  The tournament has no directed `3`-cycles and
puts the finite `random031` mirror clause and AP84 endpoint-floor packet ahead
of older raw coloring or extension shadows.

Assumption challenge: runners, gaps, fixed sections, section boundaries,
wall-crossing events, residues, cover arcs, endpoint walls, gate words,
owner-delta events, observer cuts, color boundary vectors, A000568 shadows,
AP84 floor packets, and proof obligations were considered.  The chosen quotient
preserves the two-color covering-floor gluing predicate plus the existence of a
legal LRC14 escape/discharge.  It destroys arbitrary runner order, scalar
survivor mass, raw component count, and untyped extension counts.

The challenged assumption is the central lesson:

```text
colors are labels and extensions are counts
```

is replaced by:

```text
colors are boundary charges and extensions are gluing orbits.
```
