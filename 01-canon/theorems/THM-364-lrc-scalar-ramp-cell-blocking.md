---
id: THM-364
name: lrc-scalar-ramp-cell-blocking
status: PROVED
date: 2026-05-31
session: codex-2026-05-31-S371
depends_on:
  - THM-358
---

# THM-364: Scalar-Ramp Cell Blocking

## Statement

Let `n >= 2`, `k=n-1`, and let `C` be an open cell of the discontinuity
arrangement for

```text
alpha |-> (floor(n*{alpha}), floor(n*{2 alpha}), ..., floor(n*{(n-1) alpha})).
```

Write the constant floor pattern on `C` as

```text
b_i = floor(n*{i alpha}),  i=1,...,n-1.
```

For any multiplier `m mod n`, define the scalar ramp

```text
v_i = m i mod n.
```

Then for every shift `s mod n`, at least one coordinate is blocked:

```text
(s v_i + b_i) mod n in {0, n-1}
```

for some `1 <= i <= n-1`.

Equivalently, scalar ramps block every full micro-staircase candidate cell.
They cannot be required to have a generic tight-lift witness in a
micro-staircase lemma.

## Proof

Choose any `alpha` in the open cell `C`, and set

```text
beta = alpha + s m/n  in R/Z.
```

For each `i`, since `alpha` is not on a breakpoint, write

```text
n*{i alpha} = b_i + theta_i,  0 < theta_i < 1.
```

Then

```text
n*{i beta} = n*{i alpha + s m i/n}
```

has floor congruent modulo `n` to

```text
b_i + s m i.
```

Because `v_i == m i mod n`, this gives

```text
floor(n*{i beta}) == (s v_i + b_i) mod n.
```

Suppose no coordinate were blocked.  Then every value

```text
floor(n*{i beta})
```

would lie in `{1,2,...,n-2}`.  Hence

```text
||i beta|| >= 1/n
```

for every `1 <= i <= n-1`.

By THM-358, the only times with this closed initial-segment safety property
are the nonzero unit endpoints `a/n`.  Thus `beta=a/n` for some unit `a mod n`.
But then

```text
alpha = (a - s m)/n  in R/Z,
```

so `alpha` lies on an `i=1` breakpoint of the cell arrangement.  This
contradicts that `alpha` was chosen in an open cell.  Therefore some coordinate
is blocked.

## Verification Record

`04-computation/lonely_runner_scalar_excision_s371.py` verifies this identity
for the `n=14` cell system used in the fourteen-runner work.  It reconstructs
all `812` representative open cells, checks all `14` scalar multipliers and
all `11368` shifted candidates, and records zero midpoint identity failures.

The stored output is
`05-knowledge/results/lonely_runner_scalar_excision_s371.out`.

## Notes

This theorem proves the first item in the HYP-1818 test plan.  The scalar ramp
family is not a random computational obstruction; it is exactly the
Dirichlet-equality initial segment transported through the shift
`alpha -> alpha + s m/n`.

For the fourteen-runner case, S371 also shows that the best S364 non-scalar
near-blocker is a one-coordinate defect of the scalar ramp `m=8`, and its `56`
missed cells are precisely cells uniquely blocked by that scalar coordinate.
This suggests the next finite target: prove a fragility theorem saying that
every non-scalar deformation of a scalar ramp exposes a uniquely protected
micro-staircase cell.

## Related

- THM-358: initial-segment unit skeleton.
- THM-363: LRC scalar-gauge reindexing.
- THM-360: unit endpoint divisibility filter.
- HYP-1817: fourteen-runner micro-staircase lift.
- HYP-1818: scalar-ramp excision before composite micro-staircases.
- `04-computation/lonely_runner_scalar_excision_s371.py`.
