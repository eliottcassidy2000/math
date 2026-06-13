---
source: codex-2026-06-01-S504
status: finite audit plus proof-strategy refinement
tags:
  - lonely-runner
  - endpoint-pressure
  - tournament-analysis
  - proof-search
  - n14
  - n18
---

# LRC Pressure/Private-Endpoint Overlap S504

This session tested the next proof move suggested by S500 and S503:

```text
pressure DAG layers should be matched against endpoint-peeling layers.
```

The goal was not another scalar gap scan.  It was to ask whether the current
strict deletion-relief pressure gauge is already the owner-realization needed
by THM-380, or whether it is only a useful unlabelled shadow.

## Tournament Analysis Setup

Pairwise observable:

```text
relief_i(j) = nearest_distance_i(after deleting j) - nearest_distance_i.
```

Switch/gauge:

```text
j -> i when relief_i(j) > relief_j(i).
```

This says that `j` is the more irreplaceable blocker of `i`.  Ties are left
unoriented.  The tie Hamiltonian path is not used to force missing arcs here;
the certificate is the strict pressure DAG itself.  When layer signatures are
printed, speeds are listed in the fixed numerical Hamiltonian path.

Endpoint side:

```text
build the exact S362 endpoint/interval system;
take the first peeling layer;
record interval-owner speeds and private endpoint-owner speeds.
```

Pressure side:

```text
sample the S500 pressure DAG times;
compute topological layers;
compute the transitive reduction;
compare pressure sources/sinks to endpoint owners.
```

## Computed Signal

The hard rows were the S500/S490 first-even ladders:

```text
n14-d7, n14-d14, n18-d3, n18-d9, n18-d18.
```

The first endpoint-peel owner layers were:

```text
n14-d7:   deadE=84,  deadI=102, interval owners {7,21,35,63,77,84,91}, private {77,84,91}
n14-d14:  deadE=168, deadI=205, interval owners {14,42,70,126,154,168,182}, private {154,168,182}
n18-d3:   deadE=56,  deadI=54,  interval owners {3,15,21,33,39,45,48,51}, private {33,39,45,48,51}
n18-d9:   deadE=176, deadI=168, interval owners {9,45,63,99,117,135,144,153}, private {99,117,135,144,153}
n18-d18:  deadE=352, deadI=337, interval owners {18,90,126,198,234,270,288,306}, private {198,234,270,288,306}
```

Across all `4284` pressure samples:

```text
pressure DAGs:          4284/4284
source/interval hit:    3913/4284
source/private hit:     2854/4284
sink/interval hit:      4244/4284
sink/private hit:       3464/4284
```

The sink side is stronger than the source side.  That is the key correction.
Since pressure arrows point from blocker to blocked runner, the endpoint-private
debt carriers are showing up more often as pressure sinks than as pressure
sources.

The transitive reductions are also nearly the full strict pressure relation:

```text
n14-d7:   mean TR retention 99.6%
n14-d14:  mean TR retention 99.8%
n18-d3:   mean TR retention 97.0%
n18-d9:   mean TR retention 98.3%
n18-d18:  mean TR retention 99.1%
```

So the pressure graph is already sparse and mostly irredundant.  The minimal
relation to label is not much smaller than the strict pressure shadow.

## Proof Interpretation

S500's slogan was:

```text
pressure sources are chargeable blockers.
```

S504 refines it:

```text
pressure sources are chargeable blockers,
but endpoint-private debt usually appears at pressure sinks.
```

That fits the direction of the pressure edge.  If `a -> b`, then `a` is the
blocker and `b` is the blocked/debt-bearing runner.  An endpoint-private row is
not necessarily the thing to delete first; it may be the thing whose endpoint
debt is being certified by incoming pressure.

This changes the induction shape.  A plausible layer proof is now:

```text
1. identify endpoint-private owners in the first endpoint peel layer;
2. show each private owner is either already a pressure sink or is one labelled
   edge away from a pressure sink in the transitive reduction;
3. charge incoming pressure from source blockers through the DAG;
4. peel the sink/private endpoint layer;
5. iterate, unless a labelled pressure SCC appears.
```

The current unlabelled pressure gauge is therefore not yet enough for the
THM-380 realization clause.  It gives a strong acyclic certificate, but not a
complete owner-protection realization certificate.

## What Failed Usefully

The source/private overlap is only `2854/4284`.  That rules out the overly
simple proof:

```text
endpoint-private row = pressure source.
```

The sink/interval overlap is `4244/4284`, almost universal but not exact.  The
remaining `40` samples are now the right local exceptions to inspect.  They
should reveal whether the pressure gauge needs:

```text
endpoint sign labels,
left/right boundary labels,
second-neighbor relief,
or two-clock corridor phase labels.
```

The result also says not to crush the pressure DAG to one scalar.  Its
transitive reduction keeps almost every strict edge, so the edge labels are the
data.

## Next Moves

1. For the `40` sink/interval misses, print the exact time, missing private
   owners, nearest-neighbor profile, and half-turn clock cell.
2. Reverse-audit the pressure orientation: sources are blockers, sinks are
   debt carriers.  The proof language should name both roles explicitly.
3. Add endpoint labels to transitive-reduction edges:

```text
protector speed, protected endpoint owner, endpoint value, sign, center.
```

4. Test whether every first-peel private owner is either a sink or lies in the
   down-closure of a sink in the pressure DAG.
5. Overlay this with the S502 two-clock corridor.  The `40` misses may occur
   exactly at half-turn boundary cells where the unlabelled pressure gauge loses
   an endpoint sign.

## Artifacts

```text
04-computation/lrc_pressure_private_overlap_s504.py
05-knowledge/results/lrc_pressure_private_overlap_s504.out
05-knowledge/hypotheses/HYP-1968-lrc-pressure-sink-private-overlap.md
```
