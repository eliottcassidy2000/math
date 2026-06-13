---
source: oracle-2026-05-31-S22
status: exploratory synthesis
tags:
  - lonely-runner
  - tournaments
  - nearest-neighbor
  - circular-order
  - n14
  - n15
  - n16
---

# Two-neighbor distance tournaments for LRC

The useful shift this session is small but sharp: stop asking only for the
nearest runner to the stationary runner.  On the circle, loneliness is a
two-sided condition.

For a time `t`, put the stationary runner at `0` and the moving runners at
`v_i*t mod 1`.  Let `L(t)` be the clockwise gap from the nearest runner behind
`0`, and let `R(t)` be the clockwise gap to the nearest runner ahead of `0`.
Then

```text
t is lonely exactly when L(t) >= 1/n and R(t) >= 1/n.
```

The usual `min_i ||v_i t||` is just `min(L(t),R(t))`.  Keeping both sides gives
a labelled predecessor/successor pair in a circular order.  That is where
tournament structure can enter.

## What the computation shows

`lrc_distance_tournament_lens_s22.py` samples the active `n=14,15,16` rows:
initial boundary systems, quotient ladders, and the S21 deterministic
coarse-row completions.  It records:

- the two stationary bracket gaps and their owner labels;
- nearest and second-nearest shells around the stationary runner;
- the nearest-neighbor digraph on all runners;
- a lex-completed circular tournament from the runner order;
- score histograms and cyclic triples through the stationary vertex.

The cleanest new datum is that the dyadic `n=16` quotient ladders preserve a
local stationary clearance:

```text
n16 d=2   margin +1/176
n16 d=4   margin +1/176
n16 d=8   margin +1/176
n16 d=16  margin +1/176
```

This happens while the time-domain gap ratios shrink by powers of two:

```text
1/33, 1/66, 1/132, 1/264.
```

So dyadic export is not simply "the witness gets smaller."  The safe interval
in time shrinks, but the local circular-order clearance of the stationary
runner stays fixed at the sampled midpoint.  The old ArchGap/debt product and
the new bracket-margin invariant are measuring different projections.

There is a smaller echo at `n=15`: the `d=3` and `d=5` ladders both have
bracket margin `+7/360`, even though their endpoint debts and products differ.
That supports treating the two-neighbor bracket as a separate feature in the
row/column search state.

## Why nearest-neighbor alone is not enough

The global nearest-neighbor digraph can fail to see the stationary witness.  In
many positive-gap rows, no moving runner has `0` as its nearest neighbor.  The
stationary runner is lonely because its two incident circular gaps are wide
enough, not because it is globally isolated in the unmarked nearest-neighbor
graph.

This reframes the user's "two nearest neighbors" suggestion:

```text
bad feature:  unmarked nearest-neighbor graph
good feature: predecessor/successor of the marked stationary vertex
```

The second-nearest shell is still useful as redundancy data, but the key object
is the marked bracket.

## Tournament reframe

At each time, the runner positions define a circular order.  Orient each pair
by clockwise displacement less than `1/2`, breaking collisions or half-turn ties
lexicographically if needed.  This gives a circular tournament completion.

The LRC condition is then a marked good-cut condition:

```text
stationary vertex has both incident order gaps at least 1/n.
```

As time moves, the circular order changes by arithmetic collision events.  A
counterexample would need every order state in this flow to put some labelled
runner inside one of the two forbidden arcs adjacent to the stationary vertex.
That is much closer to tournament good-cut/SCC language than the scalar
`min_i ||v_i t||` formulation.

The caution is that ordinary tournament summaries are too coarse.  The
`n=16 d=8` ladder has the same quasi-regular score histogram as the initial
`n=16` boundary row, but it sits in the positive-gap export branch.  Labels and
metric gap sizes must remain attached.

## Next attempt route

The next computation should turn the two-neighbor bracket into a finite state
machine:

1. Track the labelled predecessor/successor pair of `0`.
2. Add events where a bracket owner crosses `1/n` or `1-1/n`.
3. Add collision events where circular-order adjacency changes.
4. Attach the S21 coarse divisor invoices and endpoint-protection invoices to
   those states.
5. Search for a dual certificate saying the flow cannot keep both bracket rows
   blocked while also paying all forced gate and divisor rows.

That would merge the tournament order-flow picture with the existing
Hall/IP/Bruhat-Tits machinery instead of treating it as another analogy.
