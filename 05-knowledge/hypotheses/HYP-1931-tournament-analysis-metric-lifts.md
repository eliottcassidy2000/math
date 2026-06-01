---
id: HYP-1931
status: OPEN
source: oracle-2026-06-01-S23
related:
  - THM-372
  - THM-373
  - THM-374
  - HYP-1895
  - HYP-1900
  - HYP-1910
  - HYP-1920
  - HYP-1932
  - HYP-1940
  - HYP-1941
  - THM-002
---

# HYP-1931: Tournament Analysis is metric lifting plus wall-crossing

## Statement

Tournament Analysis is the process

```text
raw pairwise or geometric data
  -> metric/sensor
  -> binary comparator plus tie-break path
  -> tournament
  -> tournament invariants and wall-crossing as variables move.
```

The critical distinction is between rank lifts and analyzer lifts.

S500 canonized the finite switchboard mechanism behind this process.  THM-372
proves that a pairwise switchboard plus a tie path is a well-defined
tournament.  THM-373 specializes the wall-crossing statement to integer runner
phase clocks, and THM-374 proves that the half-turn circular tournament is
transitive exactly in the open-semicircle case.

Rank lifts assign one scalar score `s_i` to each object and orient
`i -> j` when `s_i > s_j`, with a fixed Hamiltonian path as tie-break.  These
lifts are useful for order extraction but generically collapse to transitive
tournaments.

Analyzer lifts decide each edge by pair-specific data.  These include directed
flux comparisons, circular phase comparisons, symmetric metric switches,
label-shift lenses, and oriented area/volume signs.  These lifts preserve the
cyclic structure that the repo's Hamiltonian-path, OCF, good-cut, and LRC
machinery can see.

## Evidence

`tournament_analysis_metric_lifts_s23.py` implements several lifts.

The basketball pass example treats five starters as labelled vertices.  For
each pair, the edge points from the player who passed more to the other; exact
ties are completed by the fixed path `1 -> 2 -> 3 -> 4 -> 5`.  The synthetic
quarters move from a transitive guard-hub state (`H=1`) to cyclic states with
`H=11`, `H=13`, and `H=5`.

For continuous circle runners on eight vertices, the score/rank lifts all
collapse:

```text
anchor-arc-far, anchor-chord-far, row-sum-arc, row-sum-chord,
local-cell-radius, profile-entropy, receding-anchor
```

Every one of their `316/316` sampled continuous states was transitive.  In
contrast, edge-local and edge-switch lifts were transitive in only `1/290`
sampled states.  Their Hamiltonian path ranges were broad:

```text
circle phase-halfturn:     H=65..653
circle label-shift-arc:   H=107..653
circle arc-median-switch: H=51..653
cube label-shift-L2:      H=51..585
cube Linf-resonance:      H=15..659
```

The exact state-sequence clusters are also informative:

```text
phase-halfturn = pivot-area
anchor-arc-far = anchor-chord-far
```

Thus some geometric lenses are genuinely the same tournament path after
quotienting, while others define different paths through the tournament cube.

For active LRC witness rows, scalar anchor and local-cell lifts are transitive,
but phase/lens/switch lifts remain cycle-rich.  For example:

```text
n16 d=2 witness: phase 162 cycles, label-shift 158, switch 161
n16 d=8 witness: phase 168 cycles, label-shift 164, switch 156
```

This supports HYP-1895: the LRC two-neighbor bracket is a marked metric
constraint, but tournament information appears when that metric is lifted by
edge-local or edge-switch comparators.

S454 sharpens the middle layer as a switchboard: one bit-channel per unordered
pair plus a fixed Hamiltonian-path tie-break.  In the expanded samples, rank
shadows were transitive in `172/172` states, while analyzer shadows were
transitive in only `31/672` states.  The mean analyzer path had `22.95` live
edges out of `28` in the eight-vertex families.

S480 gives the same idea the language of gauges: a map
`(labels, pairwise observables, gauge, tie path) -> tournament movie`.  Its
runner chord threshold and annulus gauges stay cyclic and strongly connected
where vertex-score gauges collapse to total orders.

S471 packages the same process as a repo-level functor: pairwise observable plus
switch functional plus tie Hamiltonian path gives a tournament-valued
observable.  This aligns the applied metric language with the old tiling-cube
base-path coordinate.

## Predictions

1. A metric-to-tournament pipeline should always report whether its comparator
   is rank-based or edge-local.  Rank-only probes are likely to miss OCF
   structure.
2. Useful analogues of LRC should preserve at least one nontransitive analyzer
   lift: phase, switch, lens, flux, or volume.
3. The fixed Hamiltonian path used for tie-breaks is a structural layer, not
   bookkeeping.  It is the same labelled base path that appears in the tiling
   model and in lex completions of degenerate circular states.
4. Continuous problems should be studied as wall-crossing paths in the
   tournament cube `{0,1}^{binom(n,2)}`.  Comparator zeros are the walls.
5. Future feature extractors should record both a rank shadow and an analyzer
   shadow.  The difference between them measures how much cyclic information
   the metric carries beyond scalar ordering.
6. Switchboard features should be standard: live edges, wall count, Hamming
   diameter, Hamiltonian-path range, SCC range, and exact path clusters.

## Sources

- `04-computation/tournament_analysis_metric_lifts_s23.py`
- `05-knowledge/results/tournament_analysis_metric_lifts_s23.out`
- `07-reflections/tournament-analysis-metric-lifts-s23.md`
- `04-computation/tournament_analysis_switchboard_s454.py`
- `05-knowledge/results/tournament_analysis_switchboard_s454.out`
- `07-reflections/tournament-analysis-switchboard-s454.md`
- `04-computation/tournament_analysis_metric_lifts_s480.py`
- `05-knowledge/results/tournament_analysis_metric_lifts_s480.out`
- `07-reflections/tournament-analysis-metric-lifts-s480.md`
- `04-computation/tournament_analysis_framework_s471.py`
- `05-knowledge/results/tournament_analysis_framework_s471.out`
- `07-reflections/tournament-analysis-framework-s471.md`
- THM-372
- THM-373
- THM-374
- HYP-1895, HYP-1900
