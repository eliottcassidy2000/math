# HYP-1853: LRC labelled endpoint cycles are circular tournament cut-protection cores

**Status:** EXPLORATORY; supported by S386 synthesis and exact tournament
cut-protection toy audits.

## Claim

The labelled arithmetic endpoint cycle from THM-365 is the circular, arithmetic
analogue of the tournament good-cut/SCC protection system from THM-354.

For a tournament with Hamiltonian path

```text
v_0 -> v_1 -> ... -> v_{n-1},
```

a backward arc

```text
v_j -> v_i,  i < j
```

protects every path cut

```text
i < k <= j.
```

Thus the good-cut theorem is a labelled interval-protection theorem on a line:
cuts are endpoints, backward arcs are protectors, arc labels are type-A roots
`e_j-e_i`, and strong components are exactly the surviving protection cores.

The Lonely Runner endpoint system replaces:

```text
path cuts       -> circular rational endpoints,
backward arcs   -> forbidden intervals strictly protecting endpoints,
root labels     -> integer speed labels (u,p,m,eps,a),
SCC condensation -> quotient-layer leak/endpoint-debt order.
```

## Evidence

`lrc_tournament_labelled_cycle_bridge_s386.py` computes the cut-protection
audit:

```text
transitive T4:      scc=4, goodCuts=0, badCuts=(1,2,3)
directed C3:        scc=1, goodCuts=2, badCuts=()
3-core plus tail:   scc=3, goodCuts=2, badCuts=(3,4)
circulant C5:       scc=1, goodCuts=4, badCuts=()
```

This is exactly THM-354 in protection language:

```text
goodCuts = n - scc_count.
```

A bad cut is an unprotected endpoint.  A strong component is a protection core.
The condensation order is the tournament version of endpoint-debt descent.

The analogy also matches the repo's endpoint-transfer work:

- THM-356 says private child witnesses and triangular minors prove rank, while
  support matching alone is insufficient.
- HYP-1792 says endpoint collision triples are incidence hyperedges, not
  ordinary parent-metagraph triangles.
- THM-365 says LRC counterexamples must contain labelled endpoint cycles, but
  S384 shows bare circular-arc cycles are easy abstract mirages.

The shared moral is:

```text
the proof lives in labelled incidence, not in the unlabeled support graph.
```

## Predictions

1. THM-354 can be rewritten as a small endpoint-protection theorem for
   Hamiltonian-path cuts, making the bridge formal rather than metaphorical.
2. LRC should have a condensation-like quotient-layer order.  Any attempted
   labelled endpoint cycle should be forced to move downward in this order,
   unless it leaks a positive gap or unprotected endpoint.
3. LRC cycle obstructions should be attacked with private endpoints,
   triangular incidence minors, or leaf-peelable labelled hypergraphs, not
   overlap-graph density alone.
4. The right LRC analogue of `Omega(T)` should have labelled protection arrows
   as vertices and conflicts from incompatible speed/endpoint labels; fugacity
   should measure repair cost or endpoint debt.
5. The `n=14` seven-ladder may fail for the same reason endpoint-collision
   triples fail to be parent-metagraph triangles: its projected graph looks
   cycle-like, but its labelled incidence layer has a leaf.

## See

- THM-354
- THM-356
- THM-365
- HYP-1792
- HYP-1793
- HYP-1841
- `04-computation/lrc_tournament_labelled_cycle_bridge_s386.py`
- `05-knowledge/results/lrc_tournament_labelled_cycle_bridge_s386.out`
- `07-reflections/lrc-labelled-cycles-and-tournament-protection-s386.md`
