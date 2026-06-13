---
id: HYP-1904
status: OPEN
source: codex-2026-05-31-S452
related:
  - THM-370
  - HYP-1836
  - HYP-1890
  - HYP-1895
  - HYP-1900
  - HYP-1901
  - HYP-1902
  - HYP-1903
  - HYP-1910
---

# HYP-1904: LRC admits a runner-distance tournament lift

## Statement

The Lonely Runner problem should be lifted from a single stationary-runner
star to a three-level runner-distance structure:

```text
star     = one stationary frame
cycle    = adjacent circular gaps / two nearest neighbors
complete = all pairwise distances, indexed by K_n edges
```

For a fixed configuration of `n` runners on the circle, THM-370 says a runner
is lonely iff the two adjacent circular gaps around it are both at least
`1/n`.  Therefore a no-lonely configuration is exactly a configuration whose
safe gaps form an independent set in the cycle graph `C_n`, with the empty safe
mask removed by the average-gap constraint.

The complete pairwise-distance lift has `binom(n,2)=T_{n-1}` coordinates, the
same triangular count as tournament arcs.  A proof route should exploit
tournament-style data attached to these pairwise distances: score sequences,
directed triangle counts, cut imbalance, chamber swaps, and close-degree
profiles.

## Evidence

`lrc_runner_distance_tournament_s452.py` verifies and records the basic
counts:

```text
n=14: I(C_n,1)=843, feasible no-lonely masks=842, I(C_n,2)=16385
n=16: I(C_n,1)=2207, feasible no-lonely masks=2206, I(C_n,2)=65537
```

The `I(C_n,1)` count is the Lucas/cycle version of Zeckendorf independence.
The `I(C_n,2)` count is the same fugacity used in the tournament OCF thread.
Thus the two-neighbor LRC lift gives a direct independence-polynomial bridge
between Zeckendorf and tournament methods.

The same script compares exact snapshots:

```text
initial n=14 boundary: close-degree hist 0^14, score hist 6^14, 3-cycles 70
n14 seven-ladder:      close-degree hist 0^2 1^8 2^4, 3-cycles 110
n14 S380 gate ladder:  close-degree hist 0 1^12 2,   3-cycles 106
initial n=16 boundary: close-degree hist 0^16, score hist 7^16, 3-cycles 112
```

The near-counterexample rows are distinguishable by complete-lift tournament
features, not just by scalar maximum gap.

## Predictions

1. A chamber automaton whose states are circular orders and safe-gap masks
   should expose small forbidden patterns invisible to the one-star endpoint
   cover.
2. Known near-counterexamples should have safe-gap masks forming path/cycle
   quotients, making HYP-1902's Zeckendorf/Ostrowski normal forms concrete.
3. Pairwise distance tournaments should supply pruning features for LRC search:
   score-sequence imbalance, close-degree histograms, directed 3-cycle counts,
   and good-cut analogues.
4. A future `n=16` certificate may combine the dyadic endpoint-debt product
   law with the `C_16` safe-mask identity `I(C_16,2)=65537`.
5. In a counterexample, unsafe gaps must form an edge cover of `C_n` in every
   chamber.  Endpoint-protection rows should be reweighted as a Hall/Farkas
   obstruction to maintaining that edge cover without positive endpoint debt.

## Sources

- `04-computation/lrc_runner_distance_tournament_s452.py`
- `01-canon/theorems/THM-370-lrc-two-neighbor-cycle-independence.md`
- `05-knowledge/results/lrc_runner_distance_tournament_s452.out`
- `07-reflections/lrc-runner-distance-tournament-s452.md`
- `07-reflections/lrc-zeckendorf-bridge-s451.md`
- HYP-1903
- HYP-1902
- HYP-1901
- HYP-1900
- HYP-1895
