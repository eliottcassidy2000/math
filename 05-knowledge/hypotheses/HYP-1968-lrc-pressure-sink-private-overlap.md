---
id: HYP-1968
status: OPEN
source: codex-2026-06-01-S504
related:
  - HYP-1960
  - HYP-1961
  - HYP-1967
  - THM-379
  - THM-380
---

# HYP-1968: LRC endpoint-private debt is carried by pressure sinks

## Statement

In the hard first-even Lonely Runner rows, the first endpoint-peeling layer is
seen more naturally in the sink side of the strict deletion-relief pressure DAG
than in the source side.

With pressure oriented

```text
blocker -> blocked runner,
```

endpoint-private owners should be interpreted as debt-bearing blocked runners.
A proof by pressure peeling should charge source blockers through the pressure
DAG and peel sink/private endpoint layers, rather than identify endpoint-private
owners directly with pressure sources.

## Evidence

`04-computation/lrc_pressure_private_overlap_s504.py` compares the first S362
endpoint-peeling layer against S500 pressure DAG layers on five hard rows:

```text
n14-d7, n14-d14, n18-d3, n18-d9, n18-d18.
```

Across `4284` sampled exact time slices:

```text
pressure DAGs:          4284/4284
source/interval hit:    3913/4284
source/private hit:     2854/4284
sink/interval hit:      4244/4284
sink/private hit:       3464/4284
```

The transitive reduction retains almost all strict pressure edges:

```text
n14-d7:   99.6%
n14-d14:  99.8%
n18-d3:   97.0%
n18-d9:   98.3%
n18-d18:  99.1%
```

Thus the strict pressure shadow is sparse, acyclic, and mostly irredundant,
but unlabelled source/private matching is not universal.

## Predictions

1. The `40` sink/interval miss samples are boundary-label failures of the
   unlabelled pressure gauge, not genuine endpoint-core candidates.
2. Adding endpoint sign and center labels to transitive-reduction edges will
   turn most sink/interval misses into labelled hits.
3. A nonempty terminal endpoint core that is pressure-realized must appear as a
   labelled pressure SCC, in agreement with THM-379 and THM-380.
4. In two-clock corridor language, pressure-label misses should concentrate at
   half-turn wall crossings or adjacent cells where the endpoint sign changes.

## Proof Program

The sharpened induction target is:

```text
endpoint-private owner
  -> pressure sink or labelled down-closure sink
  -> charged by source blockers through transitive-reduction edges
  -> peel endpoint layer
  -> repeat unless a labelled pressure SCC appears.
```

This keeps both pieces of the current evidence: endpoint peeling supplies the
finite boundary debt layer, and Tournament Analysis supplies the acyclic
pairwise dependency order.

## See Also

`07-reflections/lrc-pressure-private-overlap-s504.md`,
`05-knowledge/results/lrc_pressure_private_overlap_s504.out`,
HYP-1960, HYP-1961, HYP-1967, THM-379, THM-380.
