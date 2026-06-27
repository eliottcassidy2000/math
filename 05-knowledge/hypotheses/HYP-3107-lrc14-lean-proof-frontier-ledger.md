---
id: HYP-3107
title: LRC14 Lean proof-frontier ledger
status: FORMALIZATION / conditional frontier; not a proof
source: codex-2026-06-27-S259
tangent: T1184
technique: LTI-245
tournament_technique: LTT-143
related:
  - HYP-3083
  - HYP-3085
  - HYP-3088
  - HYP-3089
  - HYP-3090
  - HYP-3091
  - HYP-3092
  - HYP-3093
  - HYP-3094
  - HYP-3095
  - HYP-3096
  - HYP-3097
  - HYP-3098
  - HYP-3099
  - HYP-3100
  - HYP-3101
  - HYP-3102
  - HYP-3103
  - HYP-3104
  - HYP-3105
  - HYP-3106
  - THM-573
  - THM-575
  - THM-576
  - THM-577
  - OPEN-Q-108
lean:
  - 04-computation/lean/TournamentH7/TournamentH7/LRCProofFrontier.lean
---

# HYP-3107: LRC14 Lean Proof-Frontier Ledger

This pass formalizes the current LRC14 proof edge as a Lean-facing frontier
interface.  It does not prove LRC14.  It makes the bleeding edge explicit
enough that each open analytic or combinatorial step can later be replaced by
a theorem without changing the top-level wiring.

After the final rebase, HYP-3099 sharpens the tournament side of this frontier:
the cap-optimality exchange tournament is bounded but non-transitive, so the
right proof shape is a finite local-minima certificate rather than greedy
descent; and apex-7 versus forbidden-H=7 is a coincidence, so raw winding-H
must stay out of the terminal proof interface.

HYP-3100 then packages the broader tournament-contradiction grammar, and
HYP-3105 extends it into an obstruction-transfer atlas.  Together with the
`TournamentH7.LRCBleedingEdgeFrontier` wrapper, H/H-spectrum contradictions
are legal only after the encoding functor, preserved LRC predicate, destroyed
coordinate, sidecar discharge, and faithful transfer payload have all been
made explicit.

The later rebases over HYP-3101/HYP-3102/HYP-3103/HYP-3104/HYP-3106 sharpen
the open fields: `CoverageExtremality` should talk to the normal-fan/Cech
component packet, the miss-count PGF zero signal, and the maximizer signal
atlas; residual classifier glue should be audited by first-obstruction
cocycles; and HYP-3106 perspective quotients should retain the sidecar
required by the next observer/cut operation.

One more mainline fetch adds a naming guardrail and a negative theorem-shape
lesson.  The current tree still overloads `HYP-3101` between the normal-fan
component file and the S31ah tournament certificate toolkit in the index, but
the old HYP-3103 split has been repaired: `HYP-3103` now names the
miss-count PGF-zero signal, while `HYP-3106` names the perspective-groupoid
controlled-forgetting functors.  This frontier therefore cites the route name
when a namespace has multiple historical meanings.  The S31ah certificate
battery validates the H/Omega engine but also says its coarse LRC14 use is
vacuous:
apex-7 is the order-2 antipodal symmetry of `14 = 2*7`, not the forbidden
H=7/Omega-K3 event.  Any future tournament certificate must enter through a
fine-scale or packet-preserving observer, not through raw H numerology.

The new module is:

```text
TournamentH7.LRCProofFrontier
```

and it is imported by `TournamentH7.lean`.

## Current Split

The current state-of-proof split is:

```text
solved:
  q-witness gate
  level-7 lift sieve
  pair-Pascal cap RHS and THM-577 symbolic value for k=10,11
  terminal Mreach readout

open:
  bounded-core coverage extremality for k=8,9,10
  finite local-minima certificate for the non-transitive cap-exchange graph
  reflection-Perron / order-3/order-4 certificate
  Node-3 effective Erdos-Turan peel
  finite-ruler normalized arc glue
  fine-scale replacement for the degenerate coarse winding tournament

conditional:
  finite-address residual packets imply LRC14 once every residual row emits
  a proof-bearing terminal certificate.
```

In Lean this appears as:

```lean
FrontierNode
FrontierExperiment
PairPascalCapLedger
CoverageExtremality
Node3EffectivePeel
FineScaleWindingTransfer
ResidualToFiniteAddressPackets
ResidualToObserverGluingCertificates
BleedingEdgeFrontier
lrc14_from_bleeding_edge_frontier
lrc14_from_bleeding_edge_observer_gluing_frontier
lrc14_from_bleeding_edge_packet_wrapper_frontier
```

The formalization intentionally keeps the open nodes as `Prop` fields.  This
preserves the proof boundary: the module checks the wiring and exact rational
cap arithmetic, not the missing extremality theorem.

## Exact Cap Ledger

The module proves the cap-side arithmetic now considered solved.  After the
S258/S64 rebase, THM-577 upgrades the dense `k=10,11` cap values from
search-backed table facts to symbolic closed-form overlap consequences; the
Lean side records the exact rational ledger:

```text
capRat k = k*(k+1)/182            for k=10,11,12,13
capRat 9 = pairPascalMassRat 9 - 1/4004
capRat 8 = pairPascalMassRat 8 - 1081/76440
```

Thus the RHS cap is recorded as a pair-Pascal count shadow with explicit
higher-order debt at the two binding rows.  This imports HYP-3092/THM-576 and
the new THM-577 value closure into Lean without pretending that the LHS
coverage extremality or packet optimality is solved.

## Conditional Assembly

The new theorem:

```lean
lrc14_from_residual_packet_frontier
```

says that if every nonzero 13-speed row is either discharged by an early gate
or classified as a residual row, and every residual row emits a
`FiniteAddressBranchPacket`, then the existing
`lrc14_from_cutting_edge_branch_coverage` theorem proves `LRC14Statement`.

After rebasing over the observer-gluing Lean frontier, the module also imports
`TournamentH7.LRCObserverGluingLedger` and adds:

```lean
lrc14_from_residual_observer_gluing_frontier
lrc14_from_bleeding_edge_observer_gluing_frontier
```

These close the same top-level statement through `ObserverGluingCertificate`
coverage.  Thus HYP-3098/S258 is now a formal feeder interface, not just a
prose or Python-side scout.

The module also imports `TournamentH7.LRCBleedingEdgeFrontier` and exposes the
delegating theorem:

```lean
lrc14_from_bleeding_edge_packet_wrapper_frontier
```

which closes through `BleedingEdgeFrontierCoverage`, the conservative
finite-address packet wrapper added upstream.

The packaged theorem:

```lean
lrc14_from_bleeding_edge_frontier
```

wraps the same statement in a `BleedingEdgeFrontier` record that also carries
the current open scientific obligations:

- coverage extremality,
- Node-3 peel,
- fine-scale winding transfer,
- residual packet production.

The finite-address and observer-gluing producer fields are used in the final
Lean implications.  The other fields are retained so future proofs cannot
silently forget the analytic and tournament-side debt.

## Tournament Analysis

Tournament vertices are proof obligations and formal observer nodes:

```text
q_witness_gate
level7_lift_sieve
pair_pascal_cap_rhs
bounded_coverage_extremality
reflection_perron_certificate
node3_effective_peel
fine_modp_winding_transfer
finite_ruler_glue
finite_address_packet_glue
terminal_mreach_readout
```

Pairwise observable:

```text
does node A discharge, refine, or make formally checkable an obligation that
node B only names?
```

Switch/gauge:

```text
A -> B when A preserves the LRC predicate and either replaces an open node by
a proof-bearing packet, supplies a solved exact arithmetic ledger, or exposes
a false/degenerate quotient before it enters the proof.
```

Tie Hamiltonian path:

```text
terminal_mreach_readout
> finite_address_packet_glue
> level7_lift_sieve
> q_witness_gate
> pair_pascal_cap_rhs
> bounded_coverage_extremality
> reflection_perron_certificate
> node3_effective_peel
> fine_modp_winding_transfer
> finite_ruler_glue
```

This is deliberately not a runner tournament.  The coarse mod-14 winding
tournament is known to degenerate at exactly the binding scale by apex-7
antipodal ties, so it is demoted to an experiment target rather than used as a
proof vertex.

## Creative Tests Opened

The Lean file names eight experiment vertices:

```text
eberlein_hankel_dip
reflection_perron_dual
fine_modp_winding_H
observer_sheaf_overlap
finite_ruler_sampler
node3_erdos_turan_tail
pascal_scissors_nonequidecomp
equivalence_triad_audit
```

The best next tests are:

1. Instantiate `CoverageExtremality` with the exact `p0`/coverage observable
   and a consecutive-cluster constructor for `k=8,9,10`; use HYP-3099's
   bounded/non-transitive exchange verdict to enumerate and certify the finite
   local minima instead of assuming greedy descent.
2. Replace `FineScaleWindingTransfer` with a mod-`p` or sector-pair
   observable that avoids the mod-14 antipodal tie degeneracy.
3. Turn the `PairPascalCapLedger` debt at `k=8,9` into an Eberlein/Hankel
   degree-2 to degree-4 moment-gap certificate.
4. Build a residual classifier over HYP-2963/outside-bank rows whose output is
   either an `EarlyGateCertificate`, a concrete `ObserverGluingCertificate`, or
   a concrete `FiniteAddressBranchPacket`.
5. Use the observer-gluing ledger fields from HYP-3095/HYP-3097/HYP-3098 as
   the data required by `ResidualToObserverGluingCertificates`, then compress
   certified rows to `FiniteAddressBranchPacket` when that stronger route is
   available.
6. Add HYP-3100 tournament-certificate columns to
   `BleedingEdgeFrontierCoverage`: encoding functor, preserved predicate,
   destroyed coordinate, sidecar discharge, and terminal contradiction status.
7. Test the HYP-3103/HYP-3104 miss-count PGF zeros and maximizer signal atlas
   as fine-scale coverage invariants:
   does zero confinement separate the same rows as coverage extremality without
   forgetting the HYP-3101 component packet or HYP-3102 first obstruction?
8. Add the HYP-3093/HYP-3097 equivalence triad explicitly to frontier
   experiments: equinumerosity is the Pascal/base-count shadow,
   equidecomposability is the sector-pair/component/endpoint-owner scissors
   sidecar, and equidistribution is the legal Haar/Weyl limiting law only
   after resonance and observer debt are named.

## Assumption Challenge

Rejected assumption: the right formal vertex set is runners, residues, arcs,
or raw Hamiltonian-path counts.  HYP-3099 adds a sharper rejection: H=7 is not
the apex-7 proof carrier, and the cap-exchange tournament is useful as a
non-transitivity diagnostic rather than as a transitive proof engine.
The S31ah certificate-toolkit rebase strengthens this rejection: the H gap
engine is real and reusable, but it is a terminal contradiction only after the
encoding is a complete tournament shadow of the active LRC predicate.

Alternative vertex sets considered: proof obligations, solved gates,
coverage-mass functionals, moment-dual certificates, sector-pair matrices,
fine mod-`p` winding observers, residual classifiers, Node-3 tail rows,
finite-ruler sampling obligations, Lean packet constructors, and terminal
Mreach readouts.

Preserved LRC predicate: the ability to derive `LRC14Statement` from a
proof-bearing packet cover.

Destroyed information: all analytic content inside coverage extremality,
Node-3, fine-scale winding transfer, and finite-ruler sampling until those
fields are instantiated with real theorems.

The challenged hidden point is that formalization is not a final cleanup
stage.  Here it is a discovery instrument: it says exactly which creative
experiments can move the proof frontier because they would replace an open
`Prop` field by a theorem.
