---
id: HYP-1835
status: EXPLORATORY
source: codex-2026-05-31-S378
related:
  - HYP-1820
  - HYP-1821
  - HYP-1823
  - HYP-1827
  - HYP-1829
  - HYP-1830
  - HYP-1834
  - HYP-1833
---

# HYP-1835: LRC recursion is governed by an operation-metagraph state vector

## Statement

The Lonely Runner problem should be tracked across `n` by a two-layer
metagraph, not by a single obstruction score.

For denominator `n=k+1`, the residue layer has:

```text
nodes  = scalar-gauge-normalized residue vectors in (Z/nZ)^k
edges  = one-coordinate residue changes
height = missed micro-staircase cell count
labels = torsion order, gcd layer, exposed-cell repair deficit
```

The endpoint layer has:

```text
nodes  = primitive k-speed sets
height = max-gap ratio, boundary count, unprotected endpoints, peel depth
labels = additive gates, multiplicative gates, divisor edges, unit skeleton
```

As `n` changes, the useful recursive state is therefore

```text
(missed_cells, repair_defect, endpoint_core, phi(n), tau(n),
 selected_torsion_order, operation_closure, product_sum_seed_type).
```

This is the LRC analogue of the repo's tournament metagraph program: the
geometry lives in the changing landscape, not in one scalar invariant.

## Evidence

S378 computes scalar-puncture moat metrics for `4 <= n <= 22`.  Composite
denominators select proper torsion layers:

```text
n=14: best delta order 2, gcd 7, missed 56
n=15: best delta order 3, gcd 5, missed 120
n=21: best delta order 3, gcd 7, missed 224
```

Prime denominators in the same scan select full-order delta layers:

```text
n=13,17,19: best delta order n, gcd 1
```

The one-step non-reverting repair ledger gives exposed-cell deficits:

```text
n=14: gain 56, new exposed 308, ratio 11/2
n=15: gain 60, new exposed 280, ratio 14/3
n=20: gain 180, new exposed 1220, ratio 61/9
n=22: gain 220, new exposed 1386, ratio 63/10
```

The endpoint sample shows initial segments are boundary-only equality spines,
while replacing `n-1` by `n` inserts the mandatory divisibility gate and
immediately moves to positive-gap, endpoint-core-empty examples.

The companion S378 natural-operation metric script records the matching
operation view:

- addition's one-shadow is the complete transitive order;
- multiplication's one-shadow is the divisor skeleton;
- product-sum minima change by factor-packing type as `n` crosses
  prime/composite seams;
- feature vectors for hard LRC families need both endpoint and operation
  closure coordinates.

## Interpretation

The product-sum defect law from THM-361 is a useful model.  There, `1`s repair
the gap between multiplicative product and additive sum.  In the LRC residue
system, local coordinate moves can repair old missed cells only by creating a
new exposed package.  The right invariant is a repair balance, not just a
coverage count.

The scalar ramps play the role of the additive complete shadow: they are exact
equality objects that must be quotiented.  The divisor/torsion layers play the
role of the multiplicative sparse shadow: they survive the quotient and govern
where the obstruction leaks.

## Predictions

1. Prime-to-composite and composite-to-prime transitions in LRC proof searches
   will be better predicted by `phi(n)`, `tau(n)`, and selected torsion order
   than by raw missed-cell counts.
2. Any `n=14` branch-and-bound certificate should improve sharply if nodes are
   ordered by repair deficit and torsion label, not only by Hamming distance
   from the scalar spine.
3. Endpoint-core searches should rank speed sets by the combined feature vector
   `(max_gap_ratio, unprotected, peel_depth, divisor_edges, operation_closure)`.
4. Product-sum factor-packing transitions should mark useful test denominators
   for LRC, because both systems are measuring dense additive completion
   meeting sparse multiplicative structure.

## Test Plan

1. Build the normalized residue metagraph for `n=14` on progressively larger
   support shells and record its repair-deficit edge weights.
2. Extend the S378 scalar-puncture atlas beyond `n=22` with cached cell systems
   and plot selected torsion order against `phi(n)` and `tau(n)`.
3. Add endpoint feature vectors to the S373/S374 disproof search ranking.
4. Compare product-sum minimal seed changes against LRC moat jumps to see
   whether factor-packing type predicts the strongest repair deficits.

## Sources

- `04-computation/lonely_runner_operation_metagraph_s378.py`
- `05-knowledge/results/lonely_runner_operation_metagraph_s378.out`
- `04-computation/natural_lrc_recursive_modes_s378.py`
- `05-knowledge/results/natural_lrc_recursive_modes_s378.out`
- `07-reflections/natural-lrc-recursive-metagraph-s378.md`
- HYP-1834
