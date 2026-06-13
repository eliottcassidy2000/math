---
status: connector synthesis + exact quotient audit
tags: [LRC, n14, round-tournament, A000016, merged-metagraph, converse, tiling-isomorphism, HYP-2089]
source: oraclebox1-connector-2026-06-03-S575
---

# LRC's round necklace body already lives in the merged quotient

Three recent threads were saying the same thing in different languages:

1. HYP-1998/S574: open LRC runner tournaments are ROUND, hence counted by
   `A000016`.
2. HYP-2087: binary lonely time words are fixed by time reversal `t -> -t`.
3. The tiling directive: understand the map from structured tilings to merged
   isomorphism class nodes.

The missing connector is:

```text
T_runner(-t) = T_runner(t)^op
```

for the open half-turn circular comparator.  Reflection of time is converse of
the runner tournament.  Therefore the right open-body target is not raw
`A000016`; it is the converse-merged image of `A000016` inside `G_m/Z_2`.

The S575 audit makes this finite at the fourteen-runner frontier.  For
`n=14`, `m=13`:

```text
4096 labelled round d-vectors
  -> 316 round/A000016 classes
  -> 190 converse-merged round nodes
```

There are `64` self-converse round classes and `126` non-fixed converse pairs.

So the open part of LRC(14) is a `190`-node merged quotient, not the full
`A000568(13)` universe and not even the raw `316` round classes.  The boundary
seam is where the remaining proof data must live.

This also clarifies the role of HYP-2088.  The `D_q/U_a/N_j` clock-blocker
ledger is not another binary time-word statistic.  It is the label layer to
attach to the self-converse/boundary seam before quotienting forgets ownership:

```text
merged round or boundary node
  -> denominator blockers D_q
  -> unit-shell blockers U_a
  -> n-clock blockers N_j
  -> endpoint-owner/core labels
```

That is the concrete next fibre table.  It is small enough on the open side
(`190` nodes at `n=14`) and still compatible with the repo's primary merged
metagraph geometry.
