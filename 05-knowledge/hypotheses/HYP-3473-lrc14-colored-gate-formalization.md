---
id: HYP-3473
title: LRC14 colored gate formalization interface
status: FORMALIZATION / Lean interface and conditional assembly; not an LRC14 proof
source: codex-2026-06-29 formalization session turning HYP-3471's colored gate-reservoir theorem target into a sorry-free Lean interface
tangent: T1433
technique: LTI-433
tournament_technique: LTT-333
lean_module: 04-computation/lean/TournamentH7/TournamentH7/LRCColoredGateFormalization.lean
result: 05-knowledge/results/lrc14_colored_gate_formalization_lean_codex_20260629.out
reflection: 07-reflections/lrc14-colored-gate-formalization-codex-20260629.md
related:
  - HYP-3471
  - HYP-3462
  - HYP-3470
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
  - HYP-3438
  - HYP-3436
  - HYP-2595
  - HYP-2594
  - HYP-2593
  - THM-523
  - OPEN-Q-108
---

# HYP-3473: LRC14 Colored Gate Formalization Interface

## Claim

HYP-3473 is the Lean-facing formal interface for HYP-3471.  It does not prove
the finite geometric implication from the component-cover bank, and it does not
prove LRC14.  Its contribution is to make the current proof route exact:

```text
dead component
  -> low-rank E/branch survivor gate
  -> typed residue + branch + adjacency + delta payload retained
  -> one legal terminal exit supplies Mreach >= 1/14
  -> LRC14Statement by the existing Mreach skeleton
```

The open finite producer is named as:

```text
DeadCoverEBranchSoundness Row gates hasDeadComponent
```

and the global LRC-facing producer is named as:

```text
ColoredGateGlobalCoverage
```

This keeps the HYP-3471 result proof-facing without pretending that the sampled
bank is already a theorem.

## Lean Artifact

Lean module:

```text
04-computation/lean/TournamentH7/TournamentH7/LRCColoredGateFormalization.lean
```

Stored build output:

```text
05-knowledge/results/lrc14_colored_gate_formalization_lean_codex_20260629.out
```

Build command:

```text
cd 04-computation/lean/TournamentH7
lake build TournamentH7.LRCColoredGateFormalization
```

The module defines the endpoint/gate vocabulary:

```text
GateEndpointKind = E | B0 | B1
GateBranchMask = branch0 | branch1 | both
GateAdjacency = left_bad_edge | right_bad_edge | two_sided
TypedEndpointResidue
ColoredSurvivorGate
OneEOneBranch
LowRankGate
EBranchLowRankGate
```

It proves the basic logical separations:

```text
OneEOneBranch gate -> not SameBranchGate gate
OneEOneBranch gate -> not CrossBranchGate gate
EBranchLowRankGate gate -> gate.endpointRank <= 2
```

It also records the HYP-3471 count ledger exactly:

```text
rowsAudited=135
rowsWithDeadComponents=130
lowRankGates=8666
eBranchLowRankGates=7002
sameBranchLowRankGates=1482
crossBranchLowRankGates=182
deadRowsWithEBranchLowRankGate=130
deadRowsWithoutEBranchLowRankGate=0
deadRowsOnlyEBranchLowRankGates=28
deadRowsWithSameBranchLowRankGate=102
deadRowsWithCrossBranchLowRankGate=54
```

and proves the arithmetic ledger checks:

```text
deadRowsWithEBranchLowRankGate = rowsWithDeadComponents
deadRowsWithoutEBranchLowRankGate = 0
0 < eBranchLowRankGates
```

## Conditional Assembly

The terminal packet type is:

```text
ColoredGateTerminalPacket (v : Fin 13 -> Int)
```

with exits:

```text
ap84_corridor_splice
ap84_color_grid_placement
random031_gate_gluing
component_conductance_menger
owner_current_imbalance
two_adic_descent
signed_spec_rprime
named_residual_debt
```

Each terminal packet carries the final discharge:

```text
(1 : Real) / 14 <= Mreach v
```

so the top-level formal theorem is deliberately conditional:

```text
lrc14_from_colored_gate_global_coverage :
  ColoredGateGlobalCoverage -> LRC14Statement
```

This is the correct formal shape for the current proof frontier.  HYP-3462
should fill the AP84 corridor-splice exit, HYP-3470 fills exact AP84 CRT
placement when needed, HYP-3461/HYP-3460/HYP-3459/HYP-3458 carry color and
extension sidecars, HYP-3455 names the random031 gluing stress, HYP-3451 names
the Menger/Green-current route, and two-adic/SPEC exits remain legal named
debts rather than hidden assumptions.

## Axiom Audit

The stored build output includes `#print axioms` hooks.  The count and packet
interface lemmas are sorry-free.  In particular:

```text
eBranch_gate_of_dead_cover: no axioms
hyp3471_dead_rows_all_have_e_branch_gate_count: no axioms
fullColoredGatePayload_keeps_boundary_data: no axioms
```

The conditional top-level theorem uses the existing `LRCFourteenSkeleton`
machinery through `lonely_of_Mreach_ge`; its reported foundational axioms are
the inherited classical/propositional/quotient ones from that skeleton, not new
colored-gate assumptions.

## Tournament Analysis

Tournament vertices are proof carriers, not runners, arcs, or raw colors:

```text
dead_positive_e_branch_gate
full_colored_gate_word
structural_color_sidecar
typed_mod14_gate_word
ap84_four_color_packet
endpoint_kind_coloring
numeric_14_residue
raw_color_count
```

Lean records the HYP-3471 score ledger:

```text
66, 66, 61, 58, 54, 42, 31, 7
```

and proves in Lean that the carrier count is `8`, the E/branch theorem target
is maximal, raw color count is minimal, and the full colored gate word ties the
dead-positive target.

Pairwise observable: retained proof payload, ability to feed a terminal route,
compression legality, and scalar-loss penalty.

Switch/gauge: the carrier retaining the stronger theorem predicate wins; ties
use the declared colored-gate route order.

Tie Hamiltonian path:

```text
dead_positive_e_branch_gate
-> full_colored_gate_word
-> structural_color_sidecar
-> typed_mod14_gate_word
-> ap84_four_color_packet
-> endpoint_kind_coloring
-> numeric_14_residue
-> raw_color_count
```

## Assumption Challenge

Alternate vertex sets considered: runners, arcs, fixed sections, residues,
cover arcs, survivor gates, `E_safe` components, branch-wall events, Fourier
modes, matroid topes, and proof obligations.

Chosen vertex set: proof-carrier obligations and terminal packets.  This
preserves the LRC predicate because a terminal packet contains the actual
`Mreach >= 1/14` discharge and the global theorem consumes exactly that field.

Destroyed information: exact finite interval geometry, row enumeration,
owner-support multiplicities, and proof of universality for the bank.  HYP-3473
does not hide that destruction.  It exposes it as the producer obligations
`DeadCoverEBranchSoundness` and `ColoredGateGlobalCoverage`.

## Next Pull

1. Instantiate `DeadCoverEBranchSoundness` for the HYP-3438/HYP-3453 finite
   gate bank from exact interval data, not from informal row counts.
2. Prove a component-cover current/cut lemma: a dead island with no E/branch
   leak forces impossible branch-current divergence or a Menger cut with named
   terminal sidecars.
3. Fill terminal packet exits for AP84, random031, conductance, owner-current,
   two-adic descent, and signed-SPEC/Rprime routes.
4. Only after those producers exist, instantiate `ColoredGateGlobalCoverage`
   and apply `lrc14_from_colored_gate_global_coverage`.
