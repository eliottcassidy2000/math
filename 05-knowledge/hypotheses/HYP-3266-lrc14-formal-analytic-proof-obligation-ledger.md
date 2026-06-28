---
id: HYP-3266
title: LRC14 formal and analytic proof-obligation ledger
status: LEDGER / proof-frontier audit; not an LRC14 proof
source: codex-2026-06-28
tangent: T1349
technique: LTI-349
tournament_technique: LTT-249
script: 04-computation/lrc14_formal_analytic_obligation_ledger_codex_20260628.py
result: 05-knowledge/results/lrc14_formal_analytic_obligation_ledger_codex_20260628.out
reflection: 07-reflections/lrc14-formal-analytic-proof-obligation-ledger-codex-20260628.md
related:
  - HYP-3267
  - HYP-3265
  - HYP-3259
  - HYP-3258
  - HYP-3257
  - HYP-3255
  - HYP-3254
  - HYP-3253
  - HYP-3252
  - HYP-3251
  - HYP-3250
  - HYP-3249
  - HYP-3248
  - HYP-3247
  - HYP-3136
  - HYP-3132
  - OPEN-Q-108
---

# HYP-3266: LRC14 Formal and Analytic Proof-Obligation Ledger

## Claim

The current LRC14 frontier is no longer a vague collection of proof ideas.  It
can be written as a finite obligation ledger with three layers:

```text
closed formal glue       = Lean already proves the theorem wiring
finite import boundary   = exact computations/certificates that need packaging
analytic/rigidity core   = the genuinely open mathematics
```

The executable ledger records `18` proof obligations.  The important correction
from HYP-3253/HYP-3252 is that the topological/index frame describes the AP
saddle, while the S-dependent proof still lives in the floor.  The important
correction from HYP-3253/S81 and HYP-3255/S82 is that the naive "large v removes
exactly 1/7" argument is false on resonant apex rows; the theorem-facing
analytic statement is off-grid bulk survivor positivity, not generic `6/7`
survival.  HYP-3257, HYP-3258, HYP-3259, and HYP-3265 further sharpen the
rigidity side: unit equioscillation is only a rank-3 coordinate, the census
splits into binding/covering layers, the tight locus has real-manifold flex
directions, and the six unit contacts form a contact graph that must survive
before scalarizing.

## Closed or Nearly Closed

Closed in Lean:

```text
O00 concrete Mreach compactness readout
O01 denominator-sieve saturation gates
O04 cap RHS / pair-Pascal ledger
O05 small-k cover-set pigeonhole and p0 monotonicity
O06 concrete witness floor from p0 wide bound
O07 gK8 per-shape Delsarte and arithmetic imports
```

Closed by exact finite computation but not yet packaged as theorem certificates:

```text
O02 unit-witness construction for AP/GW/dilations
O03 bounded AP/GW single-swap margin
O08 contact-holonomy quotient-curvature repair
```

Conditional formal glue already exists:

```text
O09 residual packet to finite-address / observer-gluing certificates
```

## Open Core

The remaining proof work is concentrated in five analytic/rigidity obligations:

```text
O10 wide cover bound hp0cap for binding k=8..12
O11 corrected witnessG2/rhoGlob floor
O12 Part A / off-grid bulk survivor positivity
O13 gK8 concentration extremality
O14 doublet R-tail uniform bound
O15 full tight-locus rigidity
```

The obsolete uniform `rhoStar` 2/7 route is retained only as a warning:

```text
O17 obsolete rhoStar 2/7 route
```

The new optional bridge is:

```text
O16 Qsqrt(-7) signed-floor reorganization
```

This is not yet a theorem.  It is a basis-change experiment suggested by the
failure of the `p mod 4` sign-cancellation bridge, by HYP-3254's Qsqrt(-7)
floor reorganization, and by HYP-3267's `zeta_7` contact holonomy.

## Formal Boundary

The Lean-facing names to keep green are:

```text
LRCMreachConcrete.lonely_of_Mreach_ge
LRCFourteenSkeleton.lrc14_from_witness_floor_cases_given_nodes
LRCWitnessFloorConcrete.goodSet_witness_margin_from_wide_bound
LRCCoverBound.slowmu_coverSet_lt_cap_of_decorrelation
LRCWitnessPartA.finite_witness_pos_from_goodSet_margin_uniform_arc_bound_shapes
LRCProofFrontier.lrc14_from_residual_packet_frontier
LRCFourteenSkeleton.gK8_concentration_extremality
LRCFourteenSkeleton.doublet_Rtail_uniform_bound
```

This pinpoints exactly where formalization can continue: replace opaque
objects such as `witnessG2`, `rhoGlob`, `LyVal`, `Rtail`, and `isResidual` by
concrete definitions and then prove the existing `Prop` obligations.

## Best Current Proof Route

The most direct route is now:

```text
prove O10 hp0cap decorrelation or O13 gK8 concentration
  -> use O06 concrete witness floor from p0 wide bound
  -> prove O11/O12 witnessG2 floor and off-grid bulk survivor positivity
  -> use O00 Mreach compactness readout
  -> LRC14
```

The finite-tight-locus route is complementary:

```text
prove O15 full rigidity from finite equioscillation plus blind sidecars
  -> AP/GW/dilations are the only tight rows
  -> HYP-3250 uniform margin handles the rest
```

## Assumption Challenge

Alternate vertices considered: runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes, matroid
cells, and proof obligations.

Chosen vertices: proof obligations.

Preserved predicate: a conditional route to `Mreach >= 1/14`, hence LRC14.

Destroyed information: row geometry, endpoint-owner cells, unit-blind
residue/height moves, and exact residual family labels.  Those must re-enter
through HYP-3267 contact holonomy, HYP-3265 contact-graph case split,
HYP-3257 blind ledgers, HYP-3258/HYP-3259 census and manifold splits,
observer-gluing packets, hp0cap decorrelation, Part A off-grid survivor
positivity, or tight-locus rigidity.

## Status

This is a ledger, not a proof.  It reduces the live frontier to named
obligations and identifies the places where formalization can proceed without
inventing new mathematics.
