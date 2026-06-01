---
id: THM-374
name: half-turn-transitive-semicircle-criterion
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S500
depends_on: []
---

# THM-374: Half-turn circular tournaments are transitive exactly in a semicircle

## Statement

Let `X` be a finite set of distinct points on the unit circle with no antipodal
pairs.  Define the half-turn circular tournament by

```text
x -> y  iff  y lies clockwise from x by distance in (0,1/2).
```

Then this tournament is transitive iff all points of `X` lie in some open
semicircle.

Equivalently, for a generic runner-clock cell, the transitive half-turn
tournament is exactly the cell in which the runner configuration has an empty
semicircle, or a circular gap of length greater than `1/2`.

## Proof

Suppose first that all points lie in an open arc of length `1/2`.  Choose the
clockwise-first point `x_0` in that arc and list the remaining points in
clockwise order:

```text
x_0, x_1, ..., x_{m}.
```

For `a<b`, the clockwise distance from `x_a` to `x_b` is less than `1/2`,
because both points lie in the same open semicircle.  Hence `x_a -> x_b`.
The displayed order is therefore a transitive ordering of the tournament.

Conversely, suppose the half-turn tournament is transitive.  A finite
transitive tournament has a source vertex `x`, beating every other vertex.
For every `y != x`, the arc `x -> y` means that `y` lies in the clockwise open
half-circle from `x` to `x+1/2`.

Because `X` is finite and no point is antipodal to `x`, the maximum clockwise
distance from `x` to any other point is strictly less than `1/2`.  Moving the
left endpoint slightly counterclockwise from `x` gives an open semicircle that
contains `x` and all other points.

The final equivalence follows because a finite point set lies in an open
semicircle iff the complementary circle has a gap of length greater than
`1/2`.

## Significance

This gives a rigorous bridge between tournament structure and circular
geometry.  In the runner phase clock, `H=1` transitive cells are not arbitrary
rank-collapse artifacts: they are exactly the bunched configurations with an
empty half-circle.  This is the cleanest proven piece of the "H as loneliness
meter" story.

## Related

- THM-373
- THM-370
- HYP-1951
- `04-computation/tournament_clock_basketball_gap_s24.py`
- `07-reflections/tournament-clock-runner-walks-in-Gn-s24.md`
