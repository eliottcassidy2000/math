---
source: codex-2026-06-01-S502b
status: synthesis plus finite exact audit
tags:
  - lonely-runner
  - tournament-clock
  - tournament-analysis
  - n14
  - n18
  - wall-crossing
---

# LRC and the Tournament Clock: The Two-Clock Corridor

S24 made the runner tournament clock precise: include the observer as speed `0`,
orient each pair by the half-turn phase comparator, and let time trace a closed
walk through the circular-tournament menu.  LRC uses a different clock.  Its
walls are not pairwise half-turn events; they are anchored endpoint events

```text
||v t|| = 1/n.
```

The useful integration is therefore not "replace LRC by the tournament clock."
It is an overlay of two wall systems:

```text
half-turn clock:  t = m / (2 |s_i - s_j|)
LRC endpoint:     ||v t|| = 1/n
```

The new script `lrc_tournament_clock_overlay_s502.py` overlays these walls for
the initial and hard quotient-ladder rows at `n=14` and `n=18`.

## What Changes From S24

For the initial segments, the old picture is exact and beautiful.  The lonely
unit witnesses are half-turn clock walls:

```text
n14 initial: 6/6 lonely boundary witnesses are half-turn clock walls.
n18 initial: 6/6 lonely boundary witnesses are half-turn clock walls.
```

At `t=1/n`, the runners form a regular `n`-gon with antipodal half-turn ties
when `n` is even.  The LRC witness is a clock-wall top event in the circular
tournament menu.

The hard rows behave differently.  Their positive lonely intervals are not
single half-turn cells.  They are tiny corridors crossing a small number of
half-turn walls:

```text
n14 row-parent:  gap 5/12936 crosses 1 wall
n14 gate:        gap 5/25872 crosses 1 wall
n14 double-gate: gap 5/51744 crosses 2 walls

n18 row-parent:  gap 1/3168 crosses 2 walls
n18 gate:        gap 1/6336 crosses 1 wall
n18 double-gate: gap 1/12672 crosses 2 walls
```

So the selected LRC midpoint is a cell interior, but the witness interval itself
is a two- or three-cell corridor in the half-turn clock.

## The New Invariant

The alignment of LRC boundary witnesses with half-turn walls is stable along
the dyadic ladder:

```text
n14: 36/84 = 72/168 = 144/336 = 3/7
n18: 48/176 = 96/352 = 192/704 = 3/11
```

This ratio is now a candidate clock-corridor invariant.  It survives
row-parent/gate/double-gate scaling exactly in the audited rows, even though the
absolute endpoint debt doubles each step.

The half-turn tournament shape also distinguishes the rows.  `n14` hard rows
sit in interior circular cells with score width `3`.  The `n18` row-parent
midpoint still has score width `1` and `240` cyclic triples, matching the
near-regular clock-top profile of the `n18` initial row.  That supports the
current read of `n=18`: it is a square-core proof laboratory with more symmetry,
not a more plausible disproof target by scalar gap alone.

## Proof Implication

The LRC proof object should keep both clocks:

```text
coarse half-turn clock:
  spread, circular-menu cell, wall-crossing corridor, score/cycle profile

anchored endpoint clock:
  forbidden intervals, endpoint debt, private boundary witnesses, protection
```

A future feature vector should include

```text
phase_walls_inside_lonely_gap
boundary_phase_alignment_ratio
half_turn_score_width_sequence_across_corridor
endpoint_debt
pressure_peel_status
```

This reframes HYP-1942.  The bridge-fiber export problem is not a static
tournament snapshot.  It is a small corridor through adjacent circular
tournament cells, carrying endpoint debt from the anchored `1/n` clock.

## Next Tests

1. Record the actual left/middle/right half-turn tournament invariants across
   each hard lonely interval, not just the wall count.
2. Test whether the `3/7` and `3/11` alignment ratios can be derived from the
   odd core and skip fiber.
3. Add pressure-DAG data at each cell in the corridor; a counterexample-like
   row should create a nontrivial SCC somewhere in this corridor, not merely at
   the midpoint.
4. Extend the overlay to the one-slot nonmultiple repairs from S455/S493.

## Artifacts

```text
04-computation/lrc_tournament_clock_overlay_s502.py
05-knowledge/results/lrc_tournament_clock_overlay_s502.out
05-knowledge/hypotheses/HYP-1967-lrc-two-clock-corridor.md
```
