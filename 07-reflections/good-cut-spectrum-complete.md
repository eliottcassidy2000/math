# Good-Cut Spectrum Complete

**Session:** codex-2026-05-30
**Lean:** `TournamentH7.GoodCuts`, `TournamentH7.StaircaseConnectivity`, `TournamentH7.BucketBalance`

The good-cut coordinate has now split into two clean layers.

The first layer is pure interval geometry. A staircase tile `(hi, lo)` covers
the cut interval

```text
{lo+1, lo+2, ..., hi}.
```

Every tile interval has length at least 2. Therefore a tiling can open no cuts,
or at least two cuts, but never exactly one. The new Lean theorem completes
that story:

```lean
Tournament.StTiling.goodCutCount_spectrum
```

For `n >= 3`, a number `r` is realised as `goodCutCount` if and only if

```text
r = 0  or  2 <= r <= n-1.
```

So bucket 1 is not just absent from the examples, and not merely forbidden by
a lower bound. It is the only missing value in the abstract interval-union
model.

## What Changed Mathematically

The proof is constructive. Bucket 0 is realised by the all-down tiling. For
any `2 <= r <= n-1`, a single upward tile from vertex 0 to vertex r covers
exactly the cuts `{1, ..., r}`. Lean packages this as:

- `goodCuts_singleUp_eq_cutInterval`;
- `exists_goodCutCount_eq_of_allowed`;
- `goodCutCount_spectrum`.

This is a useful distinction. The gap theorem says "1 is impossible." The
spectrum theorem says "there are no other interval-geometric obstructions."
Any missing value in a later quotient, H-layer, or isomorphism-class statistic
must come from tournament structure, not from the cut-union coordinate itself.

## The Formal Bridge Is Now Closed

The spectrum theorem began at the tiling-coordinate level. The bridge is now
formalized:

```text
StTiling.goodCutCount = n-1
    <=> every cut has an upward crossing
    <=> the induced base-path tournament is strongly connected.
```

Lean now supplies the concrete constructor and translation lemmas:

- `StTiling.toTournament_hasBasePath`;
- `StTiling.isGoodCut_iff_crossesUpward_toTournament`;
- `StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected`.

The slogan is now a theorem:

> top good-cut bucket is exactly strong connectivity for the staircase model.

## Hypotheses Worth Testing

1. **Quotient purity.** Prior computations suggest `goodCutCount` may descend
   to merged tournament classes. If true, the spectrum theorem gives an
   honest height function on `G_n/Z_2` with exactly one missing level.

2. **Distribution by interval unions.** The spectrum is complete, but the
   distribution is not. Counting tilings by `goodCutCount = r` is a union-of-
   intervals problem on a path with interval lengths at least 2. The observed
   coefficients involving SC sub-tiling counts should be expressible by
   inclusion-exclusion over the complement gaps.

3. **Escape ratio monotonicity.** `BucketBalance.halfLine_balance` formalises
   the internal/escaping half-line split for any finite quotient. Apply it to
   the good-cut height. Does the escape ratio rise monotonically from bucket 0
   to the top bucket, or does the middle carry most of the perturbation
   pressure?

The pleasant part: the Lean theorem is small, but it narrows the research
question. From now on, failures of good-cut regularity are not failures of
the interval abstraction. They are where the tournament quotient begins to add
genuine structure.
