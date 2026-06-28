---
id: HYP-3402
title: LRC14 endpoint-owner current and tropical height-wall proof angles
status: SYNTHESIS / executable theorem-target scout; not an LRC14 proof
source: codex-2026-06-28
tangent: T1363
technique: LTI-363
tournament_technique: LTT-263
script: 04-computation/lrc14_owner_current_tropical_wall_angles_codex_20260628.py
result: 05-knowledge/results/lrc14_owner_current_tropical_wall_angles_codex_20260628.out
reflection: 07-reflections/lrc14-owner-current-tropical-wall-angles-codex-20260628.md
related:
  - HYP-3311
  - HYP-3400
  - HYP-3310
  - HYP-3301
  - HYP-3300
  - HYP-3266
  - HYP-3265
  - HYP-3260
  - HYP-3253
  - HYP-3247
  - HYP-3243
  - HYP-3236
  - HYP-3225
  - HYP-3223
  - HYP-2969
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3402: LRC14 Endpoint-Owner Current And Tropical Height-Wall Angles

## Claim

After HYP-3311, the next non-repeating proof push should not reuse the
residue-word separator as if it were global.  The curated HYP-2969 bank says
the first actual-packet ambiguity is repaired by the nonunit residue word, but
HYP-3260/HYP-3310 already warn that same-residue height moves and owner loss
can escape residue-only data outside that bank.

HYP-3402 therefore selects two different theorem targets:

```text
1. endpoint_owner_boundary_current
2. tropical_height_discriminant_wall
```

They attack the two places HYP-3311 leaves explicitly open:

```text
endpoint-owner loss
same-residue or same-v2 height flex outside the curated bank
```

The composite `owner_valuation_bicurrent` is demoted to a bridge.  It should be
used only after one of the two separate programs finds the first leak, or if an
enlarged packet bank forces endpoint owner and valuation owner to travel
together.

## Exact Scout Readout

The executable scout ranks proof carriers by retained LRC obligations minus
destroyed scarce coordinates, with corpus hits from recent LRC notes only as a
capped triage bonus.  The two selected carriers tie:

```text
endpoint_owner_boundary_current    score=89
tropical_height_discriminant_wall  score=89
owner_valuation_bicurrent          score=85
cluster_mutation_wall_crossing     score=70
ramanujan_projector_backend        score=57
boolean_minterm_certificate_dag    score=43
raw_residue_word_reuse             score=-1
```

The model deliberately demotes `raw_residue_word_reuse`.  The residue word is a
good bank-local sidecar, but HYP-3311 already used it; HYP-3402 asks what proof
carrier survives once residue-only exactness first fails.

## Mixed-Fiber Sidecar Audit

The scout reuses the unique coarse mixed fiber from HYP-3311:

```text
rows=7
kernels={'positive-Haar-open': 3, 'unit-petal-named': 4}
```

Controlled-forgetting readout:

```text
coarse_base     fibers= 1 mixed_kernel_fibers=1 max_target_width=2
residue_word    fibers= 7 mixed_kernel_fibers=0 max_target_width=1
v2_word         fibers= 6 mixed_kernel_fibers=1 max_target_width=2
owner_current   fibers= 7 mixed_kernel_fibers=0 max_target_width=1
tropical_wall   fibers= 7 mixed_kernel_fibers=0 max_target_width=1
bicurrent       fibers= 7 mixed_kernel_fibers=0 max_target_width=1
```

Minimal single-sidecar repairs over the coarse base are:

```text
('residue_word',)
('owner_current',)
('tropical_wall',)
('bicurrent',)
```

This is not a proof that owner current or tropical walls globally separate the
LRC packet.  It is a concrete direction selector: they separate the same
HYP-3311 warning fiber while naming different theorem obligations than the
residue word.

## Angle 1: Endpoint-Owner Boundary Current

Read endpoint owners and theorem exits as a finite current ledger.

```text
unit-petal-named rows     = unit-boundary inflows
positive-Haar-open rows   = covering-floor outflows
AP/GW equality            = boundary H1 stop
K33/H7                    = forbidden-state sink
strict positive open cell = dissipative current exit
```

The proposed theorem is:

```text
Every residual packet fiber has conserved owner-current, a Farkas/Green dual
certificate, an AP/GW boundary-H1 stop, a forbidden H7 lift, or a named
owner-current debt.
```

Cross-disciplinary carriers:

```text
Hodge decomposition / exact-current split
Farkas certificates / endpoint credit duals
Green-Thomson currents / effective-resistance bottlenecks
oriented-matroid cocircuits / boundary owner facets
chip-firing / divisor-flow conservation
```

This is different from HYP-3300 observability because the row is not just
"which columns separate the fiber."  It asks for a conservation or duality law:
if endpoint owners are collapsed, where did the current go?

## Angle 2: Tropical Height/Discriminant Wall

Read covering flex as a valuation-bend problem.

HYP-3311 already exhibits the warning:

```text
drop6 fattening core add180
petal 10->20
```

share the same `v2` word but have different theorem exits.  The `tropical_wall`
sidecar adds the first bending owner / wall type and separates that collision
without using the full residue word.

The proposed theorem is:

```text
Every same-residue or same-v2 covering flex crosses a Newton/secondary-fan
wall with positive off-grid floor, lands on the AP/Goddyn-Wong 12->24 hinge,
or emits named height-discriminant debt.
```

Cross-disciplinary carriers:

```text
tropical discriminants and Newton polygons
secondary polytopes / regular subdivisions
valuated matroids and tropical Plucker circuits
cluster mutation / scattering wall consistency
Hensel-Krasner valuation stability
```

This is different from HYP-3310/HYP-3311 because it does not stop at residue or
`v2` words.  It asks for the first wall crossed by the height flex and what
LRC predicate that wall preserves.

## Tournament Analysis

Vertices are proof carriers and hidden sidecars, not runners/arcs.

```text
pairwise_observable =
  weighted retained proof obligations minus destroyed endpoint/height/off-grid
  coordinates

switch/gauge =
  higher proof-carrier score; ties use fewer destroyed scarce fields
```

Exact fingerprint:

```text
score_hist={-1: 1, 43: 1, 57: 1, 70: 1, 85: 1, 89: 2}
directed_3cycles=0
hamiltonian_path_count=1
priority_path=
  endpoint_owner_boundary_current
  -> tropical_height_discriminant_wall
  -> owner_valuation_bicurrent
  -> cluster_mutation_wall_crossing
  -> ramanujan_projector_backend
  -> boolean_minterm_certificate_dag
  -> raw_residue_word_reuse
```

## Assumption Challenge

Do not assume the vertices are runners or arcs.  This session considered:

```text
endpoint owners
owner-current boundaries
valuation owners
tropical walls
secondary-fan chambers
proof exits
hidden sidecars
Farkas/Green dual certificates
```

Preserved LRC predicates:

```text
endpoint-owner memory
height-flex legality
off-grid floor status
open/boundary theorem exit
quotient descent legality
state-lift exit
```

Destroyed information:

```text
owner-current alone can forget analytic zero control and height flex;
tropical walls alone can forget endpoint owners and odd sign sidecars;
raw residue-word reuse forgets owner, height, off-grid, analytic, and sign data.
```

Challenged assumption:

```text
Once nonunit residue words separate the curated HYP-2969 bank, the global
proof should keep using residue words as the next theorem currency.
```

The HYP-3402 answer is no: residue words remain useful, but the next global
failure should be priced as endpoint-current leakage or a tropical
height-discriminant wall.

## Next Pull

Extend the HYP-3311 actual-packet bank in two directions:

1. Add endpoint-owner/current fields to the broader HYP-2963 residual sample
   and find the first fiber where the owner-current conservation law fails.
2. Add a tropical wall word to same-residue and same-v2 covering flexes,
   especially height decoys such as the HYP-3260 `2->16` and hinge-adjacent
   `12->36`, `12->96`, and `12->24` families.

The desired output is a first-leak table:

```text
residue word exact
residue word fails but owner-current works
residue/v2 fail but tropical wall works
both fail -> named owner/height/off-grid debt
```

That table would turn HYP-3402 from a proof-angle scout into a theorem-facing
packet extension.
