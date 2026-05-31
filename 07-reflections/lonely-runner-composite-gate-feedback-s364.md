---
source: codex-2026-05-31 S364
status: exploratory
tags:
  - lonely-runner
  - fourteen-runners
  - fifteen-runners
  - disproof-search
  - crt
  - endpoint-protection
---

# Lonely Runner Feedback Loop: 14, 15, and Disproof Pressure

The requested rule was:

```text
attack 14 runners -> dead end -> invent/work 15-runner route
15 dead end       -> invent/work disproof construction route
disproof dead end -> cycle back
```

I implemented that as:

```text
04-computation/lonely_runner_composite_gate_feedback_s364.py
05-knowledge/results/lonely_runner_composite_gate_feedback_s364.out
```

The script does exact forbidden-interval scoring first, then runs the heavier
S362 endpoint/core peeling only on dangerous survivors.

A follow-up pass sharpened the same loop with:

```text
04-computation/lonely_runner_frontier_feedback_s364.py
05-knowledge/results/lonely_runner_frontier_feedback_s364.out
```

This second script keeps the feedback loop explicit, but adds exact gated
boxes, the complete late-peel ledger for the hardest `n=14` leak, a
gate-overload disproof family, and a small `n=15` micro-staircase analogue.

## Lap Pattern

### 14-runner route

Reduced case:

```text
k=13, n=k+1=14.
```

The initial segment is tight:

```text
(1,2,...,13)
boundary_only
unit_skeleton=True
unprotected=6
core_E=0
```

But any counterexample must pass the `14`-gate, i.e. contain a speed divisible
by `14`.  The best local gated survivor was:

```text
(1,2,3,4,5,7,8,9,10,11,12,13,14)
gap/thresh = 0.037879
unprotected = 8
peel_depth = 27
core_E = 0
```

This is not a counterexample.  It is a long-leaking endpoint system.  The
creative pivot is to study the late peel layers, not just the gap.

The frontier pass extended the exact gated scan to `max_speed=17`:

```text
k=13, n=14, max_speed=17
pass_14_gate = 1820
positive_gap = 1820
full_measure = 0
open_cover = 0
```

The same hardest example has a 27-layer peel.  The late tail is especially
suggestive:

```text
remaining_E: 12 -> 8 -> 6 -> 4 -> 2 -> 0
subgroup:    2156, 728, 112, 140, 168, 196
```

So the dead end is not a vague failure: it is a finite draining process with a
small terminal denominator pattern.

### 15-runner route

Reduced case:

```text
k=14, n=k+1=15=3*5.
```

The initial segment is again tight:

```text
(1,2,...,14)
boundary_only
unit_skeleton=True
unprotected=8
core_E=0
```

The first gated replacement is:

```text
(1,2,...,13,15)
positive_gap
forbidden_length = 95717/96525
max_gap = 1/495
gap/thresh = 0.030303
unprotected = 12
peel_depth = 16
core_E = 0
```

This looks even cleaner than the `n=14` case.  The composite gate `15=3*5`
also leaks immediately.

The exact gated box through `max_speed=18` confirms the pattern:

```text
k=14, n=15, max_speed=18
pass_15_gate = 2380
positive_gap = 2380
full_measure = 0
open_cover = 0
```

This is now a real possible route for the 15-runner case: split the proof
between the `3`, `5`, and `15` quotient gates, then run the same
divisible-channel descent that the `14` case asks for.

### Disproof route

The disproof search sampled high-coverage gated sets for both `k=13` and
`k=14`, ranking by small positive gap.  The best random survivors had gap
ratios around `0.038-0.050`, many endpoints, many unprotected endpoints, and
empty cores.  No open cover appeared.

The important failure mode is not "almost all endpoints protected."  The
sampled sets have hundreds of unprotected endpoints.  They are globally dense
but locally leaky.

## New Common Pattern

The old unit-skeleton theorem says:

```text
initial segment (1,...,n-1) has witnesses at units a/n.
```

The new feedback-loop observation is:

```text
forcing the denominator speed n into the set breaks the unit skeleton and
creates a positive gap.
```

For the two next composite denominators:

```text
n=14: (1,2,3,4,5,7,8,9,10,11,12,13,14)
      gap/thresh = 0.037879

n=15: (1,2,3,4,5,6,7,8,9,10,11,12,13,15)
      gap/thresh = 0.030303
```

This suggests a lemma family:

```text
Composite gate leak.

Start with the tight initial segment.  Remove one speed and insert n.  If the
result passes the necessary n-gate, the endpoint system has an explicit
positive gap or an empty protection core.  More generally, any primitive set
that passes the n-gate projects to a smaller divisible-channel core; if that
projection is impossible, the system leaks.
```

## Possible 15-Runner Idea

For `n=15`, use the two prime channels separately:

```text
t=1/3  forces a speed divisible by 3,
t=1/5  forces a speed divisible by 5,
t=1/15 forces a speed divisible by 15.
```

A hypothetical core must protect the unit endpoints modulo `15`, but a speed
`15w` acts trivially on the first quotient layer and only reappears after
division by `15`.  This may make the projection cleaner than `n=14`: the
unit group has `phi(15)=8` endpoints and two independent prime gates.

## Possible Disproof Construction Idea

If a disproof exists, the feedback-loop data says it probably will not look
like a random dense gated set.  Random dense sets make many overlaps but also
many endpoint leaks.

The more plausible disproof architecture would be a deliberately engineered
protection cycle:

```text
endpoint e_i of speed v_i is protected by speed v_{i+1},
and the protecting interval endpoints are themselves protected cyclically.
```

So the next disproof search should generate endpoint-protection cycles
directly, then ask whether any speed set realizes them, rather than sample
speed sets first.

The follow-up gate-overload family makes this sharper.  Adding more multiples
of the gate does not produce a cover; it usually increases unprotected
endpoints:

```text
n=14 gate multiples 1..6: unprotected 24,28,40,52,96,148
n=15 gate multiples 1..6: unprotected 12,32,36,60,88,124
```

So "more quotient protection" is not monotone in the counterexample direction.
It protects the unit layer while creating many higher-denominator leaks.

## Micro-Staircase Add-On

The frontier pass also reconnects this feedback loop to the S363
micro-staircase route.

For `n=14`, the known coarse obstruction still blocks all `196` coarse
`r/14` candidates, but a denominator-17 prime-grid cell repairs it:

```text
bins=(0,1,2,3,4,4,5,6,7,8,9,9,10)
cell=[9/154,5/84), width=1/924
```

For `n=15`, a random/structured search found a near-obstruction blocking
`219/225` coarse `r/15` candidates and resolving at denominator `17`.  That is
not as strong as the exact `n=14` obstruction, but it is enough to suggest a
new 15-runner idea:

```text
use a 3x5 micro-staircase, not only the quotient gates.
```

The next useful script should classify cells, not speed sets: rank a candidate
by the pair

```text
(endpoint peel depth, minimum resolving cell width).
```

## Next Computation

1. Print all late peel layers for the `n=14` long-leak example with
   `peel_depth=27`.
2. Prove or refute the explicit gate-replacement formula for
   `(1,...,n-2,n)`.
3. Build a reverse endpoint-cycle generator: choose desired protection edges
   in the quotient and solve the corresponding divisibility inequalities for
   speeds.
4. For `n=15`, separate the `3`-gate and `5`-gate residues during peeling.
5. Classify `n=14` and `n=15` micro-staircase cells by resolving width.
