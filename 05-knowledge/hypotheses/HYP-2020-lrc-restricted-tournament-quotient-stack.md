---
id: HYP-2020
status: OPEN
source: codex-2026-06-01-S535
related:
  - HYP-1981
  - HYP-1982
  - HYP-1986
  - HYP-1988
  - HYP-1999
  - HYP-2018
  - HYP-2019
  - HYP-2001
  - HYP-2008
  - HYP-2009
  - HYP-2015
  - THM-381
  - THM-384
  - THM-385
  - THM-387
---

# HYP-2020: LRC lives in a restricted tournament quotient stack

## Statement

The A000568 formulation should be sharpened from one quotient to a stack of
restricted tournament-class maps.  The proof question remains:

```text
Which tournament isomorphism classes can an arithmetic LRC clock exhibit?
```

but the class universe should not be all ordinary tournaments.  It should be a
fiber product of small images:

```text
raw circular body
  -> source-deleted target class
  -> observer-source marked class
  -> outside clasp/gap-threshold class
  -> kinetic gap-flow class
  -> endpoint blocker-deficit class
  -> apex-boundary/source-sink class
  -> labelled endpoint-pressure / channel-state class.
```

The LRC target is the intersection where the observer-source condition,
THM-384 long-long adjacent gaps, and zero blocker-deficit layer agree.  A proof
of LRC should show that every primitive speed clock reaches this restricted
target stack or its THM-383 compactified boundary.  A counterexample would have
to be a closed arithmetic walk in the stack that avoids every pure source
fiber while still carrying nonempty labelled endpoint pressure.

## Evidence

`lrc_restricted_tournament_mappings_s535.py` audits exact wall and open-cell
samples for primitive clocks:

```text
total n=4, max_speed=12, 196 speed sets
total n=5, max_speed=10, 205 speed sets
total n=6, max_speed=9,  126 speed sets
total n=7, max_speed=7,    7 speed sets
```

The raw phase quotient is still too coarse:

```text
n=4: phase_runner classes=2,  mixed=2,  certifies 0/196
n=5: phase_runner classes=4,  mixed=4,  certifies 0/205
n=6: phase_runner classes=11, mixed=11, certifies 0/126
n=7: phase_runner classes=17, mixed=7,  certifies 1/7
```

The source-deleted target menu is small and certifies every audited speed set:

```text
n=4: source_deleted_phase classes=2,  pure=2,  certifies 196/196
n=5: source_deleted_phase classes=4,  pure=4,  certifies 205/205
n=6: source_deleted_phase classes=11, pure=11, certifies 126/126
n=7: source_deleted_phase classes=9,  pure=9,  certifies 7/7
```

Threshold-decorated quotients trade a label tax for class-local correctness.
Every audited threshold-decorated map has zero mixed fibers and certifies every
sampled speed set.  Representative image sizes:

```text
observer_source_marked: 11,37,116,89 classes for n=4,5,6,7
gap_threshold_necklace: 20,42,77,68
gap_kinetic_flow:       20,42,77,50
blocker_deficit_shadow: 15,29,25,15
apex_boundary_runner:   22,53,129,90
```

The meta-tournament over quotient maps is transitive under the aggregate
profile used in S535:

```text
source_deleted_phase
> blocker_deficit_shadow
> gap_kinetic_flow
> observer_source_marked
> gap_threshold_necklace
> apex_boundary_runner
> phase_runner.
```

This is not a theorem about the best possible quotient.  It is evidence that
the useful direction is a suite of small, class-pure images rather than a
single scalar H-value or the raw A000568 base.

## Abstract Metrics

These metrics should be tracked for any proposed LRC tournamentization.

```text
image_size
  Number of tournament classes visited by the arithmetic clock.

target_image_size
  Number of classes visited at lonely/source states.

source_codimension_bits
  log2(A000568(m) / target_image_size), when the map lands in ordinary
  m-vertex tournament classes.

label_tax_bits
  Extra class bits paid by colors/marks/threshold labels; negative ambient
  bits are expected when correctness is moved into the fiber.

mixed_fiber_count
  Number of classes containing both good and bad states.

purity_curvature
  mixed_fiber_count times class-distribution entropy; high values indicate
  a quotient that bends good and bad chambers together.

certification_rate
  Fraction of speed sets with a sampled state in a good-only class.

compression_bits
  log2(sampled_states / image_size), measuring quotient strength.

kinetic_torsion
  Difference between the static gap-threshold image and the kinetic gap-flow
  image.  In S535 it first appears at n=7: 68 static gap classes versus 50
  kinetic classes.

apex_label_tax
  The additional image size caused by remembering which two runners bracket
  the observer-source gap and whether the two adjacent gaps are long.

blocker_viscosity
  The number and entropy of nonzero endpoint-deficit classes below the source
  layer.  This is a local proxy for HYP-1988 repair difficulty.

compactification_index
  Difference between open source menu size and wall-compactified source menu
  size, measuring how much THM-383 boundary data is needed.

pressure_survival_rank
  Size of the owner-compatible endpoint-pressure core after projecting to a
  chosen quotient.  A counterexample-shaped class must keep this positive.

monodromy_defect
  Number of quotient classes gained after completing one full arithmetic
  period compared with one local chamber pass; a measure of global clock
  arithmetic not visible in local geometry.
```

## Predictions

1. The source-deleted target image should be exactly the HYP-1999 Ferrers
   interval menu plus THM-383 boundary classes.
2. The kinetic gap-flow quotient should become strictly smaller than the static
   gap-threshold quotient as `n` grows, because THM-387 forbids two of the four
   two-gap transitions.
3. The blocker-deficit quotient should remain very small but too tautological
   alone; its value is as a repair layer over source-deleted and gap-flow
   classes.
4. In hard frontiers such as `n=14`, wall-only/AP behavior should appear as a
   zero-dimensional intersection of source target, apex-boundary target,
   diameter-channel state, and endpoint-pressure debt.
5. Non-AP clocks should either open the source-deleted target menu in positive
   measure or export debt to a labelled pressure/channel class; no unlabelled
   raw A000568 loop should be considered counterexample-shaped.

## Next Tests

1. Run the S535 quotient stack on selected `n=14` and `n=18` hard ladders.
2. Add the HYP-2015 diameter-channel state as an eighth quotient layer.
3. Compare open Ferrers menus against source-deleted sampled images after
   separating wall states from open cells.
4. Define kinetic_torsion for reset-race ledgers rather than instantaneous
   gap derivatives.
5. Search for a finite obstruction theorem: any closed walk avoiding the pure
   source stack has positive pressure_survival_rank and therefore triggers
   THM-380.

## Sources

- `04-computation/lrc_restricted_tournament_mappings_s535.py`
- `05-knowledge/results/lrc_restricted_tournament_mappings_s535.out`
- `07-reflections/lrc-restricted-tournament-quotient-stack-s535.md`
- HYP-1981
- HYP-1982
- HYP-2018
- HYP-2019
- HYP-1999
- HYP-2001
- HYP-2008
- HYP-2009
- HYP-2015
