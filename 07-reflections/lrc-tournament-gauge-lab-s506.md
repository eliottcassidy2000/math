---
source: codex-2026-06-01-S506
status: finite gauge audit plus criterion menu
tags:
  - lonely-runner
  - tournament-analysis
  - hamiltonian-paths
  - arc-gauges
  - loneliness-metric
---

# LRC Tournament Gauge Lab S506

This session tried several ways to assign tournament arcs in LRC snapshots until
the resulting tournament shape looked useful as a loneliness-related metric.
The answer is not one universal arc rule.  It is a three-part gauge menu:

```text
ranker gauges    -> who is lonely, but H collapses;
phase gauges     -> global loneliness entropy, H is useful;
pressure gauges  -> blocker/debt dependency, score/SCC shape is useful.
```

The script audited the selected LRC times for the hard rows:

```text
n=14: initial, scale 7,14,28 with skip 6
n=18: initial, scale 9,18,36 with skip 8
```

## Criteria For A Useful LRC Tournament Gauge

A gauge is useful only after we say what it is supposed to measure.

### Criterion A: Coordinate Ranker

If the goal is to identify the loneliest runner, a scalar vertex score is fine:

```text
origin distance,
local adjacent clearance,
nearest-neighbor moat,
two-neighbor moat,
observer-blocker relief.
```

But these gauges produce transitive or nearly transitive tournaments.  Their
Hamiltonian-path count is therefore useless as an entropy metric.  They are
coordinate gauges, not shape gauges.

### Criterion B: H Entropy

If the goal is for `H(T)` to measure global loneliness/spread, the gauge must
keep cyclic phase information.  The best clean rule remains the half-turn phase
clock:

```text
i -> j iff j lies in the clockwise open half-circle from i.
```

This gives high H for spread/cyclic configurations and lower H for bunched or
hierarchical configurations.

### Criterion C: Endpoint Awareness

If the goal is to mix LRC threshold information into H, only gentle overrides
are safe.  The best candidate from this pass is:

```text
safe_phase_gate:
  if exactly one runner has adjacent clearance >= 1/n, the safe runner wins;
  otherwise use the half-turn phase rule.
```

This has strong spread correlation in the audit, but it can collapse H on some
hard rows.  So it should be treated as an endpoint-aware perturbation of
`phase_half`, not as a replacement.

### Criterion D: Pressure/Peeling

If the goal is to prove no counterexample core survives, H is not the main
statistic.  Use deletion-relief pressure:

```text
pressure_k1:
  blocker -> runner whose nearest-neighbor moat it improves more.

pressure_k2:
  blocker -> runner whose two-neighbor moat it improves more.

deficit_pressure:
  blocker -> runner whose threshold deficit it reduces more.
```

Here score width, SCC size, cyclic triples, and source/sink layers are more
meaningful than H.  `pressure_k2` was the best pressure-shape candidate in this
pass because it kept more nontrivial cycle/SCC range than `k1` or pure deficit
pressure.

## Exact Audit Highlights

The `phase_half` gauge:

```text
H_ratio range: 1.6752..2.3981
cyclic triples: 106..240
largest SCC: 14..18
spread correlation: 0.706

n14 scales 7->14->28: 2.0831 -> 1.6752 -> 1.6752
n18 scales 9->18->36: 2.3981 -> 2.0965 -> 2.0965
```

This is the clean H-loneliness metric.  It sees the first gate drop and the
post-gate H plateau from S505.

The `safe_phase_gate` hybrid:

```text
H_ratio range: 0.0204..2.3675
cyclic triples: 70..240
largest SCC: 12..18
spread correlation: 0.922

n14 scales 7->14->28: 0.0495 -> 0.2304 -> 0.2304
n18 scales 9->18->36: 0.2601 -> 0.0204 -> 0.0204
```

This is endpoint-aware, but sometimes too harsh.  It detects local safe/unsafe
threshold structure, but it can turn H from entropy into an alarm.

The scalar rankers:

```text
origin_lonely_rank,
local_clearance_rank,
nearest_rank,
two_neighbor_rank,
origin_blocker_rank
```

all had only one H value in the audit and `largest_scc=1`.  This confirms the
ranker/shape split: scalar loneliness is not tournament loneliness.

The pressure gauges:

```text
pressure_k1: H_ratio max 0.0003, cyclic triples 0..27, SCC 1..16
pressure_k2: H_ratio max 0.0172, cyclic triples 0..40, SCC 1..18
deficit_pressure: H collapsed, cyclic triples 0, SCC 1
```

Their H is small, but `pressure_k2` keeps nontrivial shape.  This is a
dependency metric, not a global loneliness entropy.

## Proposed Composite Metric

The useful LRC tournament metric should be a vector, not a scalar:

```text
LRC_TA_metric =
  (
    H_ratio(phase_half),          # global phase loneliness entropy
    H_ratio(safe_phase_gate),     # endpoint-aware alarm/perturbation
    score_width(pressure_k2),     # blocker-debt hierarchy
    largest_SCC(pressure_k2),     # pressure-core survival
    endpoint_debt,
    gap * endpoint_debt
  ).
```

Interpretation:

```text
high phase H + low endpoint alarm       = evenly spread but endpoint-safe;
phase H plateau + growing endpoint debt = pure denominator-depth recursion;
pressure SCC                            = possible counterexample-shaped core;
ranker extremes                         = individual lonely/debt carriers.
```

## Next Gauge Definitions To Try

1. Label `safe_phase_gate` arcs by which adjacent gap certified safety.
2. Replace binary safe/unsafe override by a three-zone threshold:

```text
clearance < 1/n,
1/n <= clearance < 2/n,
clearance >= 2/n.
```

3. Build a signed hybrid:

```text
phase_half edge, plus endpoint sign label, without flipping the arc.
```

This may keep H stable while adding endpoint data.

4. Try pressure only on the transitive reduction and leave half-turn phase as
the tie path, rather than numerical order.
5. Use the composite metric above as the standard row summary for future LRC
frontier scans.

## Artifacts

```text
04-computation/lrc_tournament_gauge_lab_s506.py
05-knowledge/results/lrc_tournament_gauge_lab_s506.out
05-knowledge/hypotheses/HYP-1973-lrc-gauge-trichotomy.md
```
