---
source: oracle-2026-06-01-S26; codex-2026-06-01-S26b
status: exploratory result (what H measures on the circular tournament)
tags:
  - tournament-analysis
  - tournament-clock
  - lonely-runner
  - hamiltonian-paths
  - three-cycles
  - loneliness-meter
---

# H as a Loneliness Meter: what the Hamiltonian-path count actually measures

S24 noticed that for the half-turn tournament `T(t)` of a runner system, the
directed Hamiltonian-path count `H(T(t))` falls when the runners bunch and rises
when they spread.  This session makes that statement precise and corrects the
tempting overstatement.

The correct thesis is:

> `H` is an exact detector of the half-circle bunched state (`H=1`), and above
> that it is a circular-tournament spread-class meter.  It is not a scalar
> max-gap meter and not a pointwise monotone loneliness score.

The lift is the half-turn circular tournament: `i -> j` iff `j` lies clockwise
from `i` by less than half a turn, with base-path tie-breaking.

## 1. The one sharp reading: H = 1 iff an empty semicircle exists

Across the extended samples `n=5..9`, there were zero mismatches for:

```text
H = 1   iff   max circular gap > 1/2
```

in the nondegenerate open-semicircle case.  The repo's finite clock uses
base-path tie completion, so the boundary `max_gap = 1/2` is counted on the same
transitive side; the computation therefore audits the condition as
`max_gap >= 1/2`.

This is the canonized THM-374 reading.  It is the cleanest sense in which `H`
detects loneliness: the half-turn tournament sees the `1/2` gap exactly.  In
runner language, this is the `n=2` Lonely Runner threshold, not the finer `1/n`
threshold.

## 2. Above H = 1, H is not a scalar gap meter

For `H > 1`, the max-gap ranges overlap heavily.  The same rounded max gap can
land in several distinct `H` classes:

```text
n=5: 26/74 rounded max-gap buckets map to multiple H values
n=6: 30/73
n=7: 33/70
n=8: 34/63
n=9: 36/57
```

There are direct pointwise monotonicity failures.  In one `n=7` sample,
`H=105` occurred with `max_gap=0.1812`, while a higher `H=151` occurred with the
larger `max_gap=0.4881`.  So the safe statement is not "larger H means smaller
max gap."  The safe statement is that `H` is a finite cell/class invariant of
the half-turn circular tournament, and these cells are correlated with spread.

The correlations remain strong but imperfect:

```text
n=5: corr(H,max_gap)=-0.833, corr(H,gap_entropy)=0.715
n=6: corr(H,max_gap)=-0.769, corr(H,gap_entropy)=0.672
n=7: corr(H,max_gap)=-0.700, corr(H,gap_entropy)=0.616
n=8: corr(H,max_gap)=-0.652, corr(H,gap_entropy)=0.579
n=9: corr(H,max_gap)=-0.624, corr(H,gap_entropy)=0.569
```

This decay is useful: as `n` grows, a single integer `H` is still a meaningful
class feature, but the local geometry increasingly needs companion data.

## 3. What H sees beyond scores: cyclic arrangement

The number of directed 3-cycles has an exact geometric reading:

```text
#3cycles(T) = #{triples not contained in any semicircle}
             = binom(n,3) - sum_v binom(outdeg(v),2).
```

This is a spread-triple count.  At `n=5`, it lines up perfectly with the four
circular classes:

```text
#3cycles = 0,3,4,5  <->  H = 1,9,11,15.
```

But `H` is finer than the score sequence or the 3-cycle count.  The extended
probe again finds the known collisions:

```text
n=6: scores=(2,2,2,3,3,3), #3cycles=8 -> H in {41,45}
n=7: scores=(2,2,3,3,3,4,4), #3cycles=12 -> H in {123,137}
```

So `H` is not just measuring degree spread.  It sees cyclic arrangement inside
the same score/3-cycle profile.

## 4. Two-neighbor data is genuinely extra

The script also records the static two-neighbor safety count at threshold `1/n`:
a point is counted when both adjacent circular gaps are at least `1/n`.  This is
the local datum that the half-turn tournament forgets.

Within the same `H` bucket, the safe count often ranges over several values.  For
example, in `n=7`, `H=105` has `safe@1/n` ranging from `0` to `5`.  Correlation
between `H` and this two-neighbor count is weak by `n=8,9` (`0.048`, `0.025`).

This is the operational lesson for LRC: `H` says whether the cloud is bunched or
spread at half-turn resolution; two-neighbor data says which runners have actual
local room at the `1/n` threshold.

## 5. n=14 LRC overlay: high H does not solve endpoint loneliness

The n=14 selected rows from the LRC tournament-clock overlay are a useful
stress test.  The clock vertices are the observer plus the 13 runners.

```text
row                 H          H/H0    max_gap    orig/th  safe@1/n score_width c3
n14 initial          24104937    1.000       1/14        1        14           1 112
n14 row-parent       22168229    0.920   481/3696    29/24         2           3 110
n14 gate             17826951    0.740   481/3696  4339/3696       1           3 106
n14 double-gate      17826951    0.740   481/3696  8035/7392       1           3 106
```

The hard rows are not bunched by the half-turn meter: `H` stays high, and
`max_gap` is well below `1/2`.  But the true LRC certificate is anchored at the
observer and at threshold `1/14`.  That is why HYP-1967's two-clock view is
necessary: half-turn spread/topology plus anchored endpoint protection.

## Bottom line

`H` is a good loneliness meter only after naming its resolution:

1. At half-turn resolution, `H=1` is exact: it detects open-semicircle bunching.
2. Above `H=1`, `H` is a circular tournament class meter correlated with spread.
3. It is not determined by max gap and not pointwise monotone in max gap.
4. It is finer than score and 3-cycle meters, because it sees cyclic arrangement.
5. LRC needs `H` paired with two-neighbor and endpoint-clock data.

## Artifacts

```text
04-computation/H_loneliness_meter_s26.py
05-knowledge/results/H_loneliness_meter_s26.out
05-knowledge/hypotheses/HYP-1970-h-half-turn-loneliness-meter.md
```
