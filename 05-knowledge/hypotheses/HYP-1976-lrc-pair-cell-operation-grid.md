---
id: HYP-1976
status: OPEN
source: codex-2026-06-01-S509
related:
  - HYP-1982
  - HYP-1972
  - HYP-1975
  - HYP-1974
  - HYP-1973
  - HYP-1971
  - HYP-1967
  - HYP-1965
  - HYP-1966
  - HYP-1963
  - THM-370
  - THM-374
---

# HYP-1976: LRC needs pair-cell operation-grid Tournament Analysis

## Statement

The LRC Tournament Analysis metric vector should have two layers:

```text
runner-level arc criteria
+ pair-cell operation-grid criteria
```

The runner-level layer asks which runner wins a relation.  The pair-cell layer
uses unordered runner pairs `{i,j}` as tournament vertices and asks which
distance, wall, dyadic, odd-core, or product-sum cell is the dominant
obstruction.

In this layer, the best loneliness metric is not static arithmetic by itself.
The useful threshold signal is the dynamic pair-cell danger-deficit tournament:
its tie rate measures how many pair-cells are outside the `1/n` danger sector.
The static arithmetic criteria supply branch coordinates for interpreting hard
rows: dyadic row `v2(|v_i-v_j|)`, odd core of `|v_i-v_j|`, same odd doubling
chain, and product-sum proximity to `(a-1)(b-1)=n-1`.

## Evidence

`lrc_pair_cell_operation_grid_s509.py` tests seven pair-cell criteria:

```text
edge_danger_deficit
edge_phase_distance
edge_wall_frequency
edge_dyadic_row
edge_odd_core
edge_same_odd_chain
edge_product_sum_defect
```

Each criterion declares a pair-cell observable, a switch, and the lexicographic
tie Hamiltonian path on pair-cells.

The exact small-clock scorecard gives one dominant dynamic signal:

```text
edge_danger_deficit:
  tie_rate -> close_pair_count
  mean |rho| = 0.986
```

The next dynamic signal is weaker but still meaningful:

```text
edge_phase_distance:
  origin_pair_mean_norm -> max_gap
  mean |rho| = 0.536
```

Static criteria have zero within-clock Spearman score because they do not move
with time.  They are therefore not loneliness meters in the time direction.
Their role is coordinate labeling.  In selected `n=14` and `n=18` hard rows,
the dyadic/odd/product-sum criteria produce stable transitive or tie-path
shapes, while the danger-deficit tie rate changes across unit, half-unit,
gap-midpoint, and boundary times.

The A000568 audit in the same script verifies the repo analogy.  A000568 is
computed by Burnside using only odd cycle partitions:

```text
1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056
```

This is the unlabeled tournament-counting version of odd-core survival: even
cycle orbits cannot carry a coherent complete orientation.  The LRC operation
grid splits that odd-core rule into horizontal `x+2` movement along odd cores
and vertical `x*2` doubling branches.  Product-sum equations provide the
addition/multiplication critical-pair coordinate:

```text
k=14 best seed: 1^11 + 2 2 5
k=18 best seed: 1^15 + 2 3 4
k=21 best seed: 1^18 + 3 3 3
k=24 best seed: 1^22 + 2 24
```

## Predictions

1. Corridor-level LRC scripts should emit both the S506b runner metric vector
   and the S509 pair-cell vector.
2. The pair-cell danger-deficit tie rate should remain a cheap proxy for the
   close-pair burden across larger random and structured speed families.
3. Static dyadic/odd/product-sum pair-cell criteria should separate hard-row
   families better across denominators than across times within one fixed
   family.
4. A counterexample-shaped row should require not just small runner gaps but a
   nontrivial pair-cell core after quotienting the transitive tie-path
   completion.
5. HYP-1982 predicts that the pair-cell vector becomes proof-relevant when its
   classes are threshold-decorated by zero/nonzero danger colors; rank-only
   pair-deficit classes are still too coarse.

## Sources

- `04-computation/lrc_pair_cell_operation_grid_s509.py`
- `05-knowledge/results/lrc_pair_cell_operation_grid_s509.out`
- `07-reflections/lrc-pair-cell-operation-grid-s509.md`
- HYP-1972
- HYP-1982
- HYP-1963
- HYP-1965
- HYP-1966
