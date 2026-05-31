---
id: HYP-1895
status: OPEN
source: oracle-2026-05-31-S22
related:
  - THM-357
  - THM-369
  - HYP-1866
  - HYP-1890
  - HYP-1894
  - HYP-1900
---

# HYP-1895: LRC is a two-neighbor cut problem in a circular tournament flow

## Statement

For an LRC speed set `V` with `n=|V|+1`, include the stationary runner as
vertex `0` and place all vertices at

```text
x_0(t)=0,        x_v(t)=v*t mod 1.
```

Let `L(t)` and `R(t)` be the circular predecessor and successor gaps adjacent
to `0`.  Then

```text
t is lonely  <=>  L(t) >= 1/n and R(t) >= 1/n.
```

Thus LRC can be read as a marked two-neighbor condition inside the circular
order flow of the runners.  The usual quantity `min_i ||v_i t||` is the
minimum of this two-sided bracket, but keeping both sides remembers which
labelled runners bracket the stationary vertex and how that bracket changes
under adjacent swaps, collisions, and half-turn ties.  This is the entry point
for tournament good-cut methods: every time state gives a circular tournament
completion, while a counterexample must keep at least one of the two marked
cut gaps below threshold in every order state.

## Evidence

`lrc_distance_tournament_lens_s22.py` audits boundary rows, quotient ladders,
and the S21 coarse-row completions for `n=14,15,16`.  All 15 sampled witness
rows satisfy the exact two-sided bracket condition.

The extra bracket data is not redundant with the older gap-debt ledger:

```text
n=15 d=3 ladder:   bracket margin = +7/360
n=15 d=5 ladder:   bracket margin = +7/360

n=16 d=2 ladder:   bracket margin = +1/176
n=16 d=4 ladder:   bracket margin = +1/176
n=16 d=8 ladder:   bracket margin = +1/176
n=16 d=16 ladder:  bracket margin = +1/176
```

For `n=16`, the time-domain Archimedean gap ratios run
`1/33,1/66,1/132,1/264`, while the stationary two-neighbor margin remains
`+1/176` across the dyadic ladder.  That says the dyadic export shrinks the
time interval but preserves a local circular-order clearance around the
stationary runner.

The global nearest-neighbor digraph is useful but cannot replace the marked
bracket.  In many positive-gap samples, no moving runner has the stationary
runner as its nearest neighbor, yet the stationary two-sided bracket certifies
loneliness.  The correct object is therefore not an unmarked nearest-neighbor
graph; it is a circular tournament/order plus a distinguished vertex and its
two incident gaps.

The lex-completed circular tournaments also show why score data alone is too
coarse.  For example, the `n=16 d=8` ladder has the same quasi-regular score
histogram as the `n=16` initial boundary row, but it belongs to the positive-gap
export branch.  Labelled bracket owners and metric gaps carry the missing
information.

## Predictions

1. A proof search should track the marked predecessor/successor pair of the
   stationary vertex as a state variable, not just the minimum distance.
2. Endpoint-protection rows should be recast as events where a bracket owner
   crosses one of the two threshold gates around `0`.
3. The S21 Hall/dual certificate should gain two marked cut rows for the
   stationary bracket; these rows may explain the repeated `n=16` margin
   `1/176` across dyadic depth.
4. Tournament invariants that ignore labels or metric gaps, such as score
   histograms alone, will miss LRC branch distinctions.  Useful tournament
   features should be labelled: bracket owners, outdegree of the stationary
   vertex, cyclic triples through it, and nearest/second-nearest incidence into
   it.
5. A counterexample can be attacked as a finite symbolic flow through circular
   order states: every state must force `L(t)<1/n` or `R(t)<1/n`, and every
   transition must be realizable by the arithmetic speed labels.

## Sources

- `04-computation/lrc_distance_tournament_lens_s22.py`
- `05-knowledge/results/lrc_distance_tournament_lens_s22.out`
- `07-reflections/lrc-distance-tournament-two-neighbor-s22.md`
- HYP-1894, HYP-1900
