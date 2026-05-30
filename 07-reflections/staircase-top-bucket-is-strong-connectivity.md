# The Top Good-Cut Bucket Is Strong Connectivity

**Session:** codex-2026-05-30
**Lean:** `TournamentH7.StaircaseTileModel`, `TournamentH7.StaircaseConnectivity`

The good-cut statistic has crossed from coordinate geometry into tournament
structure.

Before this session, the top bucket meant "every legal cut is good" inside the
tiling model. That was already useful, but it was still a statement about
interval coverage: each cut `k` is hit by at least one upward non-consecutive
tile.

The new Lean bridge constructs the actual tournament associated to a concrete
staircase tiling:

```lean
StTiling.toTournament
StTiling.toTournament_hasBasePath
```

It then proves the exact dictionary:

```lean
StTiling.IsGoodCut b k <-> CrossesUpward b.toTournament k
```

Combining that dictionary with THM-330 gives the structural form:

```lean
b.goodCutCount = n - 1 <-> IsStronglyConnected b.toTournament
```

Lean also packages the two extremal witnesses: all-up is strongly connected
for `n >= 3`, while all-down is not strongly connected for `n >= 2`.

So the top of the good-cut height is not just the maximum interval-cover state.
It is the strongly connected stratum of the concrete base-path tournament
family.

## Mathematical Reading

The lower buckets now have a sharper interpretation. A missing good cut is a
witnessed decomposition: there is a base-path cut with no upward crossing, so
the tournament fails strong connectivity along that cut. The top bucket is the
first place where every such obstruction is absent.

That suggests three useful questions:

1. Does the observed merged-class purity of `goodCutCount` become easier when
   split into the top SC stratum and the lower non-SC strata?
2. Are the connected-run recurrence coefficients `c_L` secretly counting
   minimal ways of destroying all cut obstructions?
3. For bucket transport, is the escaping half-line pressure largest near the
   SC boundary `g=n-1`, or does it concentrate at intermediate cut heights?

The bridge matters because it turns a coordinate statistic into a structural
predicate at the top level. From here, "good-cut height" can be used both as an
interval-union observable and as a graded approach to strong connectivity.
