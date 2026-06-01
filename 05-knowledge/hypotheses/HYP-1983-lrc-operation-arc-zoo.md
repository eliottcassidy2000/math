---
id: HYP-1983
status: OPEN
source: codex-2026-06-01-S511c
related:
  - HYP-1972
  - HYP-1975
  - HYP-1976
  - HYP-1977
  - HYP-1978
  - HYP-1979
  - HYP-1980
  - HYP-1981
  - HYP-1982
  - HYP-1963
  - HYP-1965
  - HYP-1966
---

# HYP-1983: LRC operation arc rules split into ledgers, phase perturbations, and marked bundle coordinates

## Statement

Runner-level operation arc criteria for LRC have three distinct roles:

```text
static arithmetic rankers      -> ledgers / branch coordinates
phase-overridden operation gates -> perturbations of the half-turn clock
majority operation bundles     -> new marked fiber coordinates
```

The useful A000568 analogy is therefore not "choose one tournament class and
read H."  It is a switchboard statement: a speed row gives pairwise data, and
different binary switches extract different complete relations over the same
runner set.  The unmarked tournament class counted by A000568 is only the base.
LRC-relevant information lives in the marked fiber: observer position, endpoint
thresholds, two-neighbor guard data, close-pair tie mass, and operation labels
from the `x+2`/`x*2` grid.

## Evidence

`04-computation/lrc_operation_arc_zoo_s511.py` tests twelve runner-level arc
criteria:

```text
phase_half
origin_clearance_rank
endpoint_deficit_rank
two_neighbor_guard
dyadic_height
odd_core_row
same_chain_dyadic_phase
adjacent_oddrow_phase
product_sum_gate
goldbach_lemoine_gate
danger_then_grid
operation_bundle_vote
```

The exact small-clock scorecard confirms the split.  The half-turn phase
tournament is still the global spread meter, with signed Spearman correlations

```text
phase_half: H -> max_gap  -0.733
phase_half: c3 -> max_gap -0.780
```

The sign is expected: a larger max gap means more semicircle bunching, hence
lower phase complexity.  Static arithmetic labels are fixed within a clock
movie, so their within-clock correlations vanish.  They are coordinates, not
time meters.

The hard rows show why the coordinates still matter.  On `n14 initial`,
operation gates perturb the phase tournament without discarding the phase
geometry:

```text
phase_half H              24,104,937
same_chain_dyadic_phase   18,528,105
product_sum_gate          19,613,219
goldbach_lemoine_gate     21,815,229
operation_bundle_vote         37,407
```

On the positive-gap hard rows, scalar endpoint and two-neighbor rankers collapse
to transitive tournaments (`H=1` for `N=14`) even though their observer scores
remain useful marked ledgers.  The phase-overridden arithmetic gates are often
identical to phase at the best witness time, while `danger_then_grid` and
`operation_bundle_vote` create different score histograms and strict SCC sizes:

```text
n14 row-parent best:
  phase_half H             17,826,951
  danger_then_grid H       17,022,581
  operation_bundle_vote H       4,621, strictSCC 8

n14 gate best:
  phase_half H             17,755,905
  danger_then_grid H       16,177,369
  operation_bundle_vote H      12,903, strictSCC 7
```

For `N=18`, H is not computed, but the same split is visible in directed
triangles, score width, observer outdegree, tie count, strict SCC size, and
score histograms.

The operation-grid report makes the user's prime/additive picture explicit.
Even denominators carry Goldbach gates `N=p+q`; odd denominators carry
Lemoine/Levy gates `N=p+2q`; multiplication supplies the branch coordinate
`N=2^h * odd_core` and the two-core product-sum equation

```text
(a-1)(b-1)=N-1.
```

Examples:

```text
N=14: Goldbach (3,11),(7,7); product-sum (2,14)
N=15: Lemoine (5,5),(11,2); product-sum (2,15),(3,8)
N=16: Goldbach (3,13),(5,11); product-sum (2,16),(4,6)
N=18: Goldbach (5,13),(7,11); product-sum (2,18)
```

## Interpretation

Addition is horizontal: fixed sums, adjacent odd rows, and `x+2` motion.  By
itself it tends to create rank/order data.  Multiplication is vertical:
divisibility, dyadic height, and branch refinement under `x*2`.  Product-sum
equations are the interface where additive arity is paid for by multiplicative
factor structure.

A000568 enters because every arc criterion produces a complete binary relation
before quotienting.  But the same speed row can project to many different
tournament shapes depending on the switch.  The LRC proof target is not an
unmarked class invariant; it is a marked switchboard fiber over the A000568
base.

## Predictions

1. Static operation labels should continue to separate denominator families
   better than time cells inside one fixed family.
2. Phase-overridden operation gates should be most useful at transition cells:
   half-turn walls, endpoint walls, and short hard-row lonely corridors.
3. Majority operation bundles should be read by score histogram, observer
   outdegree, strict SCC size, and directed triangles before H alone.
4. A counterexample-shaped row should survive all three layers: phase spread,
   endpoint/two-neighbor ledgers, and operation-gated marked fiber checks.

## Sources

- `04-computation/lrc_operation_arc_zoo_s511.py`
- `05-knowledge/results/lrc_operation_arc_zoo_s511.out`
- `07-reflections/lrc-operation-arc-zoo-s511.md`
- HYP-1980
- HYP-1981
- HYP-1976
- HYP-1979
