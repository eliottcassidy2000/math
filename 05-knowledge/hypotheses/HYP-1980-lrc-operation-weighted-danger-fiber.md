---
id: HYP-1980
status: OPEN
source: codex-2026-06-01-S511
related:
  - HYP-1972
  - HYP-1975
  - HYP-1976
  - HYP-1977
  - HYP-1978
  - HYP-1979
  - HYP-1963
  - HYP-1965
  - HYP-1966
  - HYP-1967
  - THM-372
  - THM-373
  - THM-374
---

# HYP-1980: LRC operation-grid labels become metrics when weighted by current danger

## Statement

The LRC/A000568 analogy should be computed with an operation-weighted danger
fiber:

```text
A000568 base quotient
+ endpoint/gap labels
+ pair-cell operation-grid labels
+ current 1/N danger mass on those labelled pair-cells
```

Static operation-grid labels by themselves are not time-loneliness metrics.
They become useful loneliness coordinates only after they are multiplied by the
current pair-cell danger deficit

```text
max(0, 1/N - dist(i,j,t)).
```

Thus addition, multiplication, dyadic branch height, odd-core chains, and
product-sum proximity should be treated as fiber coordinates over the moving
tournament quotient, not as replacements for the LRC threshold geometry.

## Evidence

`lrc_operation_grid_arc_criteria_s511.py` compares twenty runner-level arc
criteria.  It imports the S506 phase/threshold gauges and adds runner scalars
induced from pair-cell data:

```text
incident_danger_sum
dyadic_branch_pressure
odd_core_branch_count
same_odd_chain_degree
additive_shadow_degree
multiplicative_shadow_degree
product_sum_interface
dyadic_danger_curvature
odd_chain_danger
additive_danger_interface
multiplicative_danger_interface
product_sum_danger
```

Every local criterion declares the same Tournament Analysis switch:

```text
objects: runners including observer 0
observable: runner scalar induced by pair-cells and operation labels
switch: larger scalar points to smaller scalar
tie path: 0 -> 1 -> ... -> N-1
```

The small exact clock scorecard separates three roles:

```text
s506_lrc_close_sector:
  tie_rate -> close_pair_count, mean |rho| = 1.000

additive_danger_interface:
  tie_rate -> safe_gap_count, mean |rho| = 0.882

product_sum_danger:
  tie_rate -> safe_gap_count, mean |rho| = 0.840

local_two_guard:
  origin_score_norm -> lonely_vertices, mean |rho| = 0.837

incident_danger_sum:
  origin_score_norm -> max_gap, mean |rho| = 0.815

s506_phase_half:
  H -> max_gap, mean |rho| = 0.798
```

The static operation labels have zero within-clock Spearman score:

```text
dyadic_branch_pressure
odd_core_branch_count
same_odd_chain_degree
additive_shadow_degree
multiplicative_shadow_degree
product_sum_interface
```

This is not a failure.  It confirms that static arithmetic labels are branch
coordinates.  The metric appears when those labels weight the moving danger
mass.  In the criteria-choice tournament, the top five profile winners are the
hybrid operation-danger gauges:

```text
additive_danger_interface
dyadic_danger_curvature
multiplicative_danger_interface
product_sum_danger
odd_chain_danger
```

The selected `n=14` and `n=18` rows show the same split.  Static
`product_sum_interface` and `dyadic_branch_pressure` stay fixed across unit,
half-unit, gap-midpoint, and boundary times, while the corresponding danger
weighted criteria change tie rates and observer scores with current close-pair
burden.

## Interpretation

Addition and multiplication should be read as two directions in the same LRC
fiber:

```text
addition       = horizontal transport among x+2 odd-core cells
multiplication = vertical refinement along x*2 doubling branches
product-sum    = interface where additive arity is traded for multiplicative
                 factor structure, e.g. (a-1)(b-1)=N-1
```

A000568 supplies the base quotient because it counts unlabelled complete
orientations and its Burnside formula only sees odd cycle partitions.  The
LRC proof object lives above that base: endpoint labels and operation-grid
labels decide whether a quotient cell has a safe section.

## Predictions

1. Future hard-row corridor scripts should emit the vector

   ```text
   phase_H,
   lrc_close_tie_rate,
   incident_danger_origin_score,
   dyadic_danger_score_hist,
   additive_danger_score_hist,
   multiplicative_danger_score_hist,
   product_sum_danger_score_hist,
   strict pressure SCCs,
   endpoint-labelled origin bracket.
   ```

2. Static operation labels should separate speed families and denominator
   types better than time slices within a fixed speed family.
3. Hybrid operation-danger gauges should be useful for prioritizing which
   close pairs are repair opportunities inside mixed A000568 fibers.
4. A counterexample-shaped row should require all operation-weighted danger
   fibers over its visited A000568 classes to remain endpoint-forbidden.
5. The next projection-defect audit should attach these hybrid danger vectors
   to each pointed class bucket, then test whether mixed safe/unsafe buckets
   are explained by operation-weighted repair channels.

## Sources

- `04-computation/lrc_operation_grid_arc_criteria_s511.py`
- `05-knowledge/results/lrc_operation_grid_arc_criteria_s511.out`
- `07-reflections/lrc-operation-grid-arc-criteria-s511.md`
- HYP-1976
- HYP-1977
- HYP-1972
