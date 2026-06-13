---
id: HYP-1972
status: OPEN
source: codex-2026-06-01-S506b
related:
  - HYP-1951
  - HYP-1967
  - HYP-1971
  - HYP-1976
  - HYP-1950
  - HYP-1960
  - THM-370
  - THM-374
---

# HYP-1972: LRC loneliness needs an arc-criteria metric vector

## Statement

No single tournament arc rule turns the Lonely Runner Conjecture into a scalar
loneliness meter.  The useful object is a vector of marked Tournament Analysis
criteria:

```text
(phase_H, phase_score_width,
 origin_marked_score, safe_deficit_score_hist,
 origin_blocker_strict/tie_shape,
 pressure_largest_SCC, pressure_triangles,
 strict_tie_rate, edge_flip_rate)
```

The half-turn phase tournament supplies the cyclic spread coordinate.  LRC
threshold information comes from marked gauges: the stationary vertex score in
local-moat or deficit tournaments, the strict/tie rate of the close-sector
gauge, and SCC/triangle data in deletion-relief pressure gauges.

## Evidence

`lrc_arc_criteria_loneliness_s506.py` tests twelve arc-assignment criteria on
exact small runner clocks and selected `n=14`/`n=18` LRC rows.  Each criterion
declares its pairwise observable, switch/gauge, and tie Hamiltonian path
(increasing runner labels).

The exact small-clock scorecard ranks the strongest fingerprint/target pairs:

```text
lrc_close_sector       tie_rate -> safe_gap_count        mean |rho| = 0.881
local_moat_min         origin_score -> origin_clearance  mean |rho| = 0.821
phase_half             H -> max_gap                      mean |rho| = 0.790
safe_deficit           origin_score -> origin_clearance  mean |rho| = 0.758
local_moat_sum         origin_score -> origin_clearance  mean |rho| = 0.689
```

This separates three regimes.

1. `phase_half` is the best `H` gauge: `H` tracks unanchored half-turn spread
   and has mean Spearman `H~max_gap = -0.723`.
2. `local_moat_min` and `safe_deficit` are marked LRC gauges: their completed
   tournaments are often transitive or tie-path dominated, but the marked
   stationary score tracks origin clearance.
3. `k2_relief_pressure`, `threshold_deficit_pressure`, and
   `bracket_relief_pressure` are proof-shape gauges: they should be read by
   strict SCCs, directed triangles, and tie rates, not by completed `H`.

In the selected `n=14` and `n=18` hard rows, strengthened pressure gauges still
peel rather than closing into a strict pressure core: the reported strict SCCs
remain `1` and strict pressure triangle counts remain `0` in the sampled rows.
The phase-half snapshots, by contrast, stay high-H or near-regular at the
positive-gap midpoints, confirming that scalar phase spread is not the endpoint
proof object.

## Predictions

1. Any counterexample-shaped LRC search should require both endpoint debt and a
   nontrivial pressure shape event: strict SCC size greater than `1`, a strict
   pressure triangle, or a labelled failure of the current pressure gauge to
   realize endpoint-core protection.
2. The `lrc_close_sector` tie rate should be a cheap threshold-local proxy for
   safe-gap count across broader speed families.
3. `local_moat_min` and `safe_deficit` should be used as marked-score meters,
   not as `H` meters, because their completed tournaments mostly encode the tie
   path.
4. Future LRC scripts should emit the full metric vector above for selected
   times and corridors, not just `H` or scalar maximum gap.
5. The runner-level metric vector should be paired with the second-order
   pair-cell operation-grid vector from HYP-1976 when pairwise distance cells
   are the actual obstruction carriers.

## Sources

- `04-computation/lrc_arc_criteria_loneliness_s506.py`
- `05-knowledge/results/lrc_arc_criteria_loneliness_s506.out`
- `07-reflections/lrc-arc-criteria-loneliness-s506.md`
- HYP-1976
- HYP-1971
- HYP-1967
- HYP-1950
- THM-370
- THM-374
