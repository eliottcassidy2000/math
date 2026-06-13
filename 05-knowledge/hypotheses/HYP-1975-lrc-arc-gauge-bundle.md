---
id: HYP-1975
status: OPEN
source: codex-2026-06-01-S506c
related:
  - HYP-1932
  - HYP-1940
  - HYP-1967
  - HYP-1968
  - HYP-1970
  - HYP-1971
  - HYP-1972
  - HYP-1973
  - HYP-1974
  - THM-372
  - THM-373
  - THM-374
---

# HYP-1975: LRC loneliness is a gauge bundle, not a single tournament

## Statement

In LRC Tournament Analysis, no single arc-assignment criterion should be
promoted to "the" loneliness metric.  The useful object is a bundle of
tournament gauges, each recording a different binary projection of the same
runner movie:

```text
global spread        phase_half H, open_arc_density scores/c3
local crowding       close_pair, danger_band, kinetic_close switches
pressure peeling     blocker, two_neighbor, threshold_deficit relief SCCs
marked observer      observer outdegree / two-sided guard data
finite boundary      labelled endpoint debt, not runner-level endpoint ranks
```

Scalar rank gauges still matter as readable diagnostics, but they usually
collapse to transitive tournaments and lose the edge-local obstruction.  Plain
runner-level endpoint tournaments also collapse in the hard rows; endpoint data
should stay labelled by owner, center, sign, and protector before projection.

## Evidence

`04-computation/lrc_arc_gauge_zoo_s506.py` tries many arc criteria on selected
`n=14` and `n=18` LRC rows: initial, row-parent, gate, and double-gate.  It
computes exact Hamiltonian-path counts for gauges on at most 14 vertices and
records score width, score histogram, directed 3-cycles, SCC size, source/sink
counts, and observer outdegree.

The n=14 selected rows separate the channels:

```text
phase_half:
  initial H=24104937, row-parent H=22168229, gate H=17826951
  correlation with observer margin: Hratio +0.393, c3 +0.416

open_arc_density:
  row-parent H=2952757, gate H=452385
  correlation with observer margin: Hratio +0.190, c3 +0.403

danger_band_switch:
  row-parent H=4231707, gate H=3826755, c3=71 in both

close_pair_switch:
  row-parent H=109371, gate H=207651
  observer-out correlation +0.506, but H/c3 correlations are negative

two_neighbor_relief:
  row-parent H=183523, gate H=55701
  correlation with observer margin: Hratio +0.184, observer-out +0.231

blocker_relief:
  row-parent H=3601 with largest SCC 12; gate H=161 with largest SCC 10

kinetic_close_switch:
  row-parent H=2703685, gate H=1249481
```

The scalar gauges `origin_danger_rank`, `local_moat_rank`, and
`observer_relief_rank` are readable but usually transitive (`H=1`).  The
velocity-only `kinetic_escape` gauge is also inert on these hard rows.  The
plain runner-level endpoint summaries `endpoint_cross_protect` and
`endpoint_debt_rank` are transitive on every selected row, so they are not the
right endpoint compression.

For n=18, exact `H` is skipped, but the same shape appears through score/cycle
data: open-arc density has `c3=184` on the row-parent and `123/124` on the gate
levels; danger-band has `c3=117,124,112`; two-neighbor relief keeps nontrivial
cycle/SCC information while scalar ranks collapse.

## Interpretation

The sign of a correlation is not by itself the metric.  Close-pair and
threshold-deficit gauges can become more active when the observer is less
lonely, because they are crowding alarms rather than spread meters.  Their
negative `H` or `c3` correlations are still useful: they tell the proof search
where local `1/n` crowding is concentrated.

This resolves a tension after HYP-1970.  Half-turn `H` is a genuine global
loneliness/spread meter, but LRC lives at threshold `1/n` and at a marked
observer.  The useful feature vector must therefore keep several tournament
movies in parallel instead of forcing H to do endpoint, pressure, and local
crowding work at once.

## Predictions

1. A counterexample-shaped LRC row should not merely lower scalar gap; it should
   create at least one nontrivial event in the bundle, such as pressure
   `largest_scc > 1` after endpoint-private leaves are peeled, or a labelled
   endpoint cycle not visible in runner-level endpoint ranks.
2. The most promising proof metric is a small vector, not a scalar:
   `phase_H`, `open_arc_density c3`, `danger_band c3`, `two_neighbor_relief SCC`,
   `observer_out`, and labelled endpoint debt.
3. Endpoint compression should be rebuilt as a labelled switchboard over
   endpoint-runner pair-cells.  The speed-level tournament is too coarse.
4. Open-arc density is a good bridge gauge between half-turn H and local
   crowding, because it distinguishes row-parent/gate density distortion while
   staying circular rather than anchored only at the observer.

## Sources

- `04-computation/lrc_arc_gauge_zoo_s506.py`
- `05-knowledge/results/lrc_arc_gauge_zoo_s506.out`
- `07-reflections/lrc-arc-gauge-zoo-s506.md`
- HYP-1932
- HYP-1940
- HYP-1967
- HYP-1968
- HYP-1970
