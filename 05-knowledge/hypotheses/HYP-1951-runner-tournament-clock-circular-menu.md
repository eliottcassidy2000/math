---
id: HYP-1951
status: OPEN
source: oracle-2026-06-01-S24
related:
  - THM-373
  - THM-374
  - THM-369
  - HYP-1850
  - HYP-1931
  - HYP-1967
  - HYP-1971
  - HYP-1982
---

# HYP-1951: The runner tournament-clock realizes a fixed circular menu in G_n

## Statement

Lift integer-speed runners `s_0 < ... < s_{n-1}` to a tournament movie by the
half-turn phase comparator

```text
i -> j iff frac((s_j - s_i)t) in (0,1/2),
```

with ties resolved by the fixed base path.  As `t` ranges over `[0,1)`, this is
a closed walk through the tournament cube.

The open hypothesis is that, for each fixed `n`, the isomorphism classes
reachable by this clock are exactly the circular half-turn tournaments on `n`
points; the speed set selects a closed walk through this menu, but not a new
menu.

## Canonized Parts

S500 converted two pieces of the S24 clock into canon:

- THM-373 proves the finite wall decomposition.  The only arc-change times are
  `m/(2|s_i-s_j|)`, so the clock is a finite closed tournament walk.
- THM-374 proves the transitive-cell criterion.  A half-turn circular
  tournament is transitive exactly when all points lie in an open semicircle,
  equivalently when the circular gap is greater than `1/2`.

These results make the "H as loneliness meter" direction precise at the
transitive endpoint: `H=1` cells are exactly bunched/open-semicircle cells for
the half-turn clock.

S26b adds an important correction.  Beyond the transitive endpoint, `H` should
not be read as a scalar largest-gap meter: max-gap ranges overlap, rounded
max-gap buckets can have several `H` values, and higher `H` can occur with a
larger max gap than lower `H`.  The robust reading is that `H` is a
half-turn circular-tournament cell invariant, correlated with spread and finer
than score/3-cycle summaries, while LRC threshold data still needs the anchored
endpoint clock and two-neighbor gaps.

S506 refines the loneliness reading.  The transitive endpoint is an
unanchored bunching detector, while the LRC threshold `1/n` lives on the
marked, high-H side of the clock: near-regular circular cells can have many
vertices whose two adjacent gaps are both safe.  Thus `H` is useful, but the
LRC invariant must also remember the marked stationary vertex, the safe-gap
mask, and pressure/deletion data.  See HYP-1971.

S512 adds the A000568-class-fiber correction.  The open circular menu is a
real restriction, but equality walls with the fixed tie Hamiltonian path form a
larger boundary compactification.  For LRC, raw unmarked classes and
observer-marked classes are mixed in total `n=3,4`; threshold-decorated gap and
pair-cell fibers are the first class objects that separate good-only witness
states.  See HYP-1982.

## Evidence

The S24 scripts

```text
tournament_clock_s24.py
tournament_clock_identify_s24.py
tournament_clock_basketball_gap_s24.py
```

found small fixed menus in low dimensions.  For example, `n=5` visited
isomorphism classes with `H` in `{1,9,11,15}`, while the basketball pass-flux
lift on five starters reached `H=5` and `H=13`, outside the circular menu.
This supports the claim that the restriction comes from circle geometry, not
from the tournament-analysis switchboard mechanism alone.

S502 overlays this half-turn clock with the actual LRC endpoint clock for the
hard `n=14` and `n=18` rows.  Initial-segment lonely boundary witnesses are
half-turn clock walls (`6/6` for both `n=14` and `n=18`).  Hard quotient-ladder
positive gaps are instead short corridors crossing one or two half-turn walls,
with ladder-invariant boundary alignment ratios `3/7` for `n=14` and `3/11`
for `n=18`.  This refines the LRC bridge: half-turn spread data and anchored
`1/n` endpoint data are two clocks that must be overlaid, not identified.

The S506 audit found strong negative correlations between `H` and maximum
circular gap in sampled clocks, but also found bucket inversions at `n=6` and
`n=7`.  This makes the fixed-menu claim more important: `H` orders much of the
circle geometry, but it does not replace the circular order and marked LRC
data.

## Predictions

1. The circular menu for `n` is the full set of tournaments realizable by
   placing `n` points on a circle and orienting by clockwise half-turn.
2. Speed sets with the same difference-wall lattice should have the same clock
   signature up to time reparametrization and relabelling.
3. The alpha-comparator family `frac((s_j-s_i)t) in (0,1/k)` should recover
   LRC threshold cells at `1/k`, with open-arc containment replacing the
   semicircle criterion.
4. LRC hard-row proofs should record half-turn clock corridors across each
   positive lonely interval, including phase-wall count, boundary alignment
   ratio, and the tournament invariants on adjacent cells.

## Sources

- `07-reflections/tournament-clock-runner-walks-in-Gn-s24.md`
- `07-reflections/tournament-analysis-metric-lifts-s23.md`
- `04-computation/lrc_tournament_clock_overlay_s502.py`
- `07-reflections/lrc-tournament-clock-overlay-s502.md`
- THM-373
- THM-374
- THM-369
- HYP-1850
- HYP-1931
- HYP-1967
- HYP-1971
- HYP-1982
- `04-computation/h_loneliness_meter_s506.py`
- `05-knowledge/results/h_loneliness_meter_s506.out`
- `07-reflections/h-as-loneliness-meter-s506.md`
