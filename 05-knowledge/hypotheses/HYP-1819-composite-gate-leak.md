# HYP-1819: Composite gate leak

**Status:** EXPLORATORY proof-search hypothesis.

## Statement

Let `n=k+1` be composite.  A reduced Lonely Runner counterexample with `k`
speeds must contain at least one speed divisible by `n`, because otherwise
`t=1/n` is a lonely witness.

HYP-1819 says that this necessary `n`-gate is also destabilizing: once a
primitive speed set passes the gate, the endpoint-protection system either has
a positive gap or its protection core descends to a smaller speed set along
the divisible-by-`n` channel.

In particular, for the next composite targets:

```text
n=14  (fourteen runners, k=13)
n=15  (fifteen runners,  k=14)
```

the gate should create a leak rather than an all-protected core.

## Evidence

S364 cycles through:

```text
14-runner attack
15-runner idea
disproof-construction pressure
```

In two feedback laps:

- the initial segments `(1,...,13)` and `(1,...,14)` are boundary-only with
  unit skeletons;
- inserting the mandatory gate speed creates positive gaps;
- random high-coverage gated searches for `k=13` and `k=14` produce no open
  covers and no nonempty endpoint cores.

Representative data:

```text
n=14:
(1,2,3,4,5,7,8,9,10,11,12,13,14)
gap/thresh = 0.037879
peel_depth = 27
core_E = 0

n=15:
(1,2,3,4,5,6,7,8,9,10,11,12,13,15)
gap/thresh = 0.030303
peel_depth = 16
core_E = 0
```

The follow-up frontier script `lonely_runner_frontier_feedback_s364.py`
strengthened this in four ways:

- extended the exact `n=14` gated primitive box from `max_speed=16` to
  `max_speed=17`, finding `1820/1820` positive-gap cases and no full-measure
  candidates;
- ran the analogous exact `n=15` gated box through `max_speed=18`, finding
  `2380/2380` positive-gap cases and no full-measure candidates;
- printed the full 27-layer endpoint peel for the hardest `n=14` gate leak,
  which terminates at the empty core through subgroup moduli ending
  `... 2156, 728, 112, 140, 168, 196`;
- tested gate-overload and random gated disproof pressure; all survivors
  remained positive-gap with empty endpoint cores.

The same script also links this gate leak to HYP-1817: the exact `n=14`
coarse `r/14` obstruction is repaired on a prime grid, while the `n=15`
analogue already has a near-obstruction (`219/225` coarse pairs blocked)
resolved at denominator `17`.  This suggests the composite-gate descent and
micro-staircase routes are not separate tricks but two views of the same
finite boundary leak.

S373 pushed the disproof direction harder for `n=14` with speed-first
constructions.  It found no open-cover candidates, but did find a sharper
positive-gap near-miss:

```text
(1,7,14,21,28,35,49,56,63,70,77,84,91)
forbidden_length = 142/143
max_gap/thresh = 0.005411
unprotected endpoints = 84
```

This strengthens the "gate leak" picture: quotient ladders can almost cover
the circle, but the price is a large exposed endpoint layer rather than an
all-protected core.

## Why It Matters

The current published frontier proves the reduced conjecture through `k=12`.
The next cases `k=13` and `k=14` have composite denominators `14` and `15`.
Prime-polynomial ansatz methods do not directly explain this structure, but
CRT gates do.

If the composite gate always leaks or descends, then the next frontier may be
more approachable than a generic `k=13` finite check suggests.

## Disproof Search Consequence

A random dense speed set is probably a bad disproof model: it creates many
overlaps, but also many endpoint leaks.  A plausible disproof construction
would need a deliberately engineered endpoint-protection cycle.  Future
searches should generate protection cycles first and solve for speeds second.

## See

- `04-computation/lonely_runner_composite_gate_feedback_s364.py`
- `05-knowledge/results/lonely_runner_composite_gate_feedback_s364.out`
- `04-computation/lonely_runner_frontier_feedback_s364.py`
- `05-knowledge/results/lonely_runner_frontier_feedback_s364.out`
- `04-computation/lonely_runner_14_disproof_hunt_s373.py`
- `05-knowledge/results/lonely_runner_14_disproof_hunt_s373.out`
- `07-reflections/lonely-runner-composite-gate-feedback-s364.md`
- `07-reflections/lonely-runner-fourteen-disproof-hunt-s373.md`
- HYP-1816
- HYP-1828
