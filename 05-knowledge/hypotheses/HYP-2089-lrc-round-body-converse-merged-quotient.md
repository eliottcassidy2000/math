---
id: HYP-2089
status: SUPPORTED by exact S575 converse-merge audit through m=13; formula proof open
source: oraclebox1-connector-2026-06-03-S575
related:
  - HYP-1998
  - HYP-2086
  - HYP-2087
  - HYP-2088
  - THM-381
---

# HYP-2089: the LRC round body should enter the merged metagraph through the converse quotient

HYP-1998/S574 identify the open-time LRC runner sub-tournaments on `m=n-1`
vertices as exactly the ROUND tournaments, counted by `A000016`.  HYP-2087
shows that binary LRC time words are invariant under time reversal.

The connector observation is that these are the same quotient at the tournament
level:

```text
T_runner(-t) = T_runner(t)^op
```

for the half-turn circular comparator at open times.  Indeed, replacing every
phase by its negative reverses the cyclic orientation of every non-antipodal
pair.  Thus the natural open LRC class path is not the raw `A000016` round body
but its converse-merged image in `G_m/Z_2`.

## Exact audit

`04-computation/lrc_round_merged_converse_s575.py` reuses the already validated
S574 round generator and individualization-refinement canonical labeling.  For
each round isomorphism class it canonicalizes the converse tournament, counts
fixed classes, and forms converse-merged pairs.

Output:

```text
  m  valid-d   round  A000016  SC-round   2^floor  merged
  3        4       2        2         2         2       2
  4        8       2        2         2         2       2
  5       16       4        4         4         4       4
  6       32       6        6         4         4       5
  7       64      10       10         8         8       9
  8      128      16       16         8         8      12
  9      256      30       30        16        16      23
 10      512      52       52        16        16      34
 11     1024      94       94        32        32      63
 12     2048     172      172        32        32     102
 13     4096     316      316        64        64     190
```

The audit checks:

```text
round = A000016
SC-round = 2^floor((m-1)/2)     for 3 <= m <= 13
merged = (round + SC-round)/2
```

## n=14 consequence

For the LRC `n=14` runner body, `m=13`:

```text
4096 labelled round d-vectors -> 316 round classes -> 190 converse-merged nodes.
```

Of the `316` round classes, `64` are converse-fixed and `126` are non-fixed
converse pairs.

So the open LRC necklace body should be attached to the repo's merged metagraph
through a `190`-node quotient, not through all `316` round classes and not
through the full ambient `A000568(13)` set.

## Connection to HYP-2088

This gives a concrete target for the tiling-to-isomorphism directive:

```text
round labelled words -> converse-merged round nodes -> boundary/fixed seam
```

The HYP-2088 `D_q/U_a/N_j` blocker ledger should be attached at the boundary
seam, before the binary time-word quotient forgets denominator, shell, n-clock,
and owner labels.  In particular, the next useful computation is not another
binary Burnside count; it is a fibre table

```text
merged round/self-converse boundary node
  -> multiset of D/U/N blocker labels and endpoint-owner labels.
```

## Honest Scope

The time-reversal-to-converse statement is theorem-level for open half-turn
times.  The `m<=13` counts are exact computations using the previously pinned
S574 canonicalizer.  The closed formula `SC-round(m)=2^floor((m-1)/2)` is only
verified here through `m=13`; it should be proved directly from the binary
necklace/round-tournament model before being canonized.

## Files

- `04-computation/lrc_round_merged_converse_s575.py`
- `05-knowledge/results/lrc_round_merged_converse_s575.out`
