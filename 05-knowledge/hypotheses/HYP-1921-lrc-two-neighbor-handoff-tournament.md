---
id: HYP-1921
status: OPEN
source: codex-2026-05-31-S461
related:
  - THM-357
  - THM-365
  - HYP-1853
  - HYP-1900
  - HYP-1901
  - HYP-1902
  - HYP-1910
  - HYP-1920
  - HYP-1930
---

# HYP-1921: LRC counterexamples require a labelled two-neighbor handoff tournament

## Statement

The scalar nearest-runner function

```text
min_v ||v t||
```

forgets the proof-relevant boundary data.  A Lonely Runner counterexample
would have to produce an all-protected labelled handoff core in which every
endpoint owned by a threshold runner has at least one strict protector and no
two-sided boundary witness remains.

The minimal finite lift records, at every forbidden endpoint:

```text
left nearest runner,
right nearest runner,
nearest distance,
second distinct nearest distance,
threshold owner set,
strict protector set.
```

This data forms a labelled directed handoff graph:

```text
owner speed -> protector speed
```

with labels given by endpoint value, side, sign, slack, and denominator depth.
The proposed tournament import is that any all-protected LRC handoff core must
be peelable or have positive slack divergence, analogous to tournament
good-cut/SCC protection and endpoint-transfer rank.

## Evidence

`lrc_tournament_two_neighbor_s461.py` audits four exact endpoint systems.

```text
initial n=14:
  176 endpoints, 6 witnesses
  42 redundant handoffs, 128 single handoffs
  handoff SCC sizes (12,1)

n14 seven-ladder:
  1150 endpoints, 84 witnesses
  356 redundant handoffs, 710 single handoffs
  handoff SCC sizes (13)

n14 S380 gate ladder:
  2298 endpoints, 168 witnesses
  710 redundant handoffs, 1420 single handoffs
  handoff SCC sizes (13)

initial n=16:
  232 endpoints, 8 witnesses
  60 redundant handoffs, 164 single handoffs
  handoff SCC sizes (14,1)
```

In every audited row, the compressed pairwise protection-dominance tournament
on speeds is complete and transitive.  Thus speed-level tournament compression
is too coarse; it erases the witnesses.  The nontrivial object is the labelled
endpoint handoff graph, where dense SCCs can coexist with exposed boundary
points.

At representative witnesses, the stationary runner has one nearest neighbor
on each side exactly at threshold, e.g. in the initial `n=14` row:

```text
t=1/14:
  left = 13 at 1/14
  right = 1 at 1/14
```

The S380 row has many fragile single handoffs:

```text
owner at threshold + exactly one strict protector inside.
```

Those are the LRC analogue of individual backward arcs protecting cuts.

## Interpretation

Tournament structure should not be searched for in arbitrary tournaments on
the speed set.  It should be searched for in:

```text
round tournaments from circular runner order,
cut-protection handoff graphs at endpoints,
Omega-style compatibility graphs of labelled repair packets.
```

The two-neighbor lift is the finite object connecting these views.  It keeps
the left/right boundary strands that the scalar minimum discards.

Zeckendorf enters if the residual handoff core quotients to a path-like carry
graph.  Then adjacent repairs cannot coexist; they carry into the next layer,
giving a no-adjacent normal form for endpoint debt.

## Predictions

1. A true all-protected LRC endpoint core, if generated abstractly, will fail
   an arithmetic slack or two-sided nearest-neighbor condition before it
   becomes a realized speed set.
2. Single-handoff rows should supply private pivots in an endpoint-transfer
   matrix, while redundant handoffs should collapse into compatible repair
   packets.
3. Speed-level majority protection tournaments will often be transitive and
   therefore insufficient; labelled endpoint packets are necessary.
4. The `n=14` proof should combine HYP-1920's forced odd fan with a
   two-neighbor handoff peel: the even bridge fiber must eventually expose a
   left/right boundary witness or export Zeckendorf-height debt.
5. A useful LRC Omega graph should take labelled handoffs as vertices and put
   edges between handoffs that cannot coexist because of side, slack,
   denominator-depth, or speed-budget conflicts.

## Sources

- `04-computation/lrc_tournament_two_neighbor_s461.py`
- `05-knowledge/results/lrc_tournament_two_neighbor_s461.out`
- `04-computation/lrc_pairwise_tournament_s470.py`
- `05-knowledge/results/lrc_pairwise_tournament_s470.out`
- `07-reflections/lrc-tournament-two-neighbor-lift-s461.md`
- THM-357
- THM-365
- HYP-1853
- HYP-1900
- HYP-1902
- HYP-1920
