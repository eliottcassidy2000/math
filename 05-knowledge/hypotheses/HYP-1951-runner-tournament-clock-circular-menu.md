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

## Predictions

1. The circular menu for `n` is the full set of tournaments realizable by
   placing `n` points on a circle and orienting by clockwise half-turn.
2. Speed sets with the same difference-wall lattice should have the same clock
   signature up to time reparametrization and relabelling.
3. The alpha-comparator family `frac((s_j-s_i)t) in (0,1/k)` should recover
   LRC threshold cells at `1/k`, with open-arc containment replacing the
   semicircle criterion.

## Sources

- `07-reflections/tournament-clock-runner-walks-in-Gn-s24.md`
- `07-reflections/tournament-analysis-metric-lifts-s23.md`
- THM-373
- THM-374
- THM-369
- HYP-1850
- HYP-1931
