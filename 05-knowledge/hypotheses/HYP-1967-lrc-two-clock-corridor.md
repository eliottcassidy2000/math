---
id: HYP-1967
status: OPEN
source: codex-2026-06-01-S502b
related:
  - HYP-1942
  - HYP-1951
  - HYP-1952
  - HYP-1964
  - THM-373
  - THM-374
---

# HYP-1967: LRC hard rows are two-clock corridors through the tournament clock

## Statement

The half-turn tournament clock and the LRC endpoint clock should be kept as two
separate but overlaid wall systems.

For speeds with observer included as `0`, the half-turn tournament clock has
walls

```text
t = m / (2 |s_i - s_j|).
```

The LRC endpoint clock has walls

```text
||v t|| = 1/n.
```

The initial-segment LRC witnesses sit directly on half-turn clock walls at the
near-regular top of the circular tournament menu.  In contrast, the hard
first-even quotient ladders for `n=14` and `n=18` have positive lonely intervals
that are short corridors through a small number of adjacent half-turn clock
cells.  A proof should therefore track a clock-corridor vector, not only a
single tournament snapshot.

## Evidence

`04-computation/lrc_tournament_clock_overlay_s502.py` overlays the two exact
wall systems for the initial, row-parent, gate, and double-gate rows of `n=14`
and `n=18`.

Initial rows:

```text
n14 initial: 6/6 lonely boundary witnesses are half-turn clock walls.
n18 initial: 6/6 lonely boundary witnesses are half-turn clock walls.
```

Hard positive-gap rows:

```text
n14 row-parent:  gap 5/12936 crosses 1 half-turn wall
n14 gate:        gap 5/25872 crosses 1 half-turn wall
n14 double-gate: gap 5/51744 crosses 2 half-turn walls

n18 row-parent:  gap 1/3168 crosses 2 half-turn walls
n18 gate:        gap 1/6336 crosses 1 half-turn wall
n18 double-gate: gap 1/12672 crosses 2 half-turn walls
```

The half-turn/LRC boundary alignment ratio is invariant along each
row-parent/gate/double-gate ladder:

```text
n14: 36/84 = 72/168 = 144/336 = 3/7
n18: 48/176 = 96/352 = 192/704 = 3/11
```

At the selected LRC gap midpoints, the half-turn tournament remains strongly
connected.  `n18` row-parent is especially structured: it has score width `1`
and `240` cyclic triples, matching the near-regular clock-top profile of the
`n18` initial row even though it is a positive-gap quotient ladder.

## Predictions

1. Dyadic scaling of a fixed hard ladder preserves the boundary-alignment ratio
   with the half-turn clock.
2. A useful LRC Tournament Analysis feature vector should include:
   `phase_walls_inside_lonely_gap`, `boundary_phase_alignment_ratio`,
   half-turn score-width sequence across the corridor, and anchored endpoint
   debt.
3. The `n=18` square-core packet should look more regular than `n=14` in the
   half-turn clock, but not more dangerous, because the endpoint clock still
   shows exported debt.
4. A counterexample-like row should not merely have a small scalar gap; it
   should produce a clock corridor with no endpoint-private peel and a
   nontrivial pressure/SCC event across the adjacent half-turn cells.

## See Also

`04-computation/lrc_tournament_clock_overlay_s502.py`,
`05-knowledge/results/lrc_tournament_clock_overlay_s502.out`,
`07-reflections/lrc-tournament-clock-overlay-s502.md`,
HYP-1942, HYP-1951, HYP-1952, HYP-1964, THM-373, THM-374.
