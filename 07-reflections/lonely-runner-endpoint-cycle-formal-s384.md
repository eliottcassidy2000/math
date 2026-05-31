# Lonely Runner Endpoint-Cycle Formal Session S384

**Session:** codex-2026-05-31-S384  
**Script:** `04-computation/lonely_runner_endpoint_cycle_formal_s384.py`  
**Stored output:** `05-knowledge/results/lonely_runner_endpoint_cycle_formal_s384.out`  
**Theorem:** THM-365

This session tried to formalize the next small piece of the endpoint-incidence
picture from S379.

The useful reduction is:

```text
counterexample
=> full all-protected endpoint system
=> nonempty protection core
=> directed endpoint-protection cycle
=> labelled arithmetic cycle satisfying THM-357 inequalities
```

The proof is elementary but clarifying.  In a nonempty core, every remaining
endpoint has a protector interval whose boundary endpoints are also in the
core.  So every endpoint has an outgoing arrow to another core endpoint.  A
finite directed graph with positive outdegree everywhere has a directed cycle.

## The False Friend

The negative side mattered more than I expected: abstract circular-arc topology
is not enough.

On `Z/3Z`, the three arcs

```text
(0 -> 2), (1 -> 0), (2 -> 1)
```

already cover all cells and protect every endpoint.  The S384 search finds
abstract all-protected covers for `q=3,...,9`, including restricted short-arc
versions.

So a pure endpoint-core theorem cannot prove LRC.  The core/cycle object needs
its arithmetic labels:

```text
owner speed u
protector speed p
center m
sign eps
strict inequality |p*(n*m+eps)-a*n*u| < u
```

That is the actual shape.  The endpoint graph is the skeleton; the integer
labels are the muscle.

## LRC Audit

The sampled LRC systems all peel to empty terminal core:

```text
initial n=14:       coreE=0, first removed=6
n14 seven-ladder:  coreE=0, first removed=84
n14 single-gate:   coreE=0, first removed=24
initial n=15:       coreE=0, first removed=8
n15 3x5 ladder:    coreE=0, first removed=150
```

This gives a clean interpretation of the seven-ladder.  It has a tiny visible
gap, but its endpoint-cycle attempt dies immediately with `84` unprotected
endpoints.  It is a near-miss in horizontal gap width, not a near-cycle in
the endpoint core.

## New Search Rule

For disproof construction, speed-first feels backwards.  A better order is:

1. Generate a directed endpoint-protection cycle shape.
2. Attach owner/protector labels and strict integer inequalities.
3. Solve realizability by integer speeds.
4. Only then audit full interval measure.

This turns HYP-1828 into a concrete solver architecture.  It also explains why
product-sum and torsion coordinates keep reappearing: they may be exactly the
short labelled-cycle attempts that almost close a quotient layer before
exposing a higher one.

## New Questions

1. What is the smallest labelled endpoint cycle that satisfies all local
   integer inequalities for `n=14`?
2. Can a cycle slack potential prove that every primitive labelled cycle leaks
   either a positive gap or an unprotected endpoint?
3. Does every cycle touching a unit endpoint force a speed divisible by `n`,
   then descend to a smaller quotient layer?
4. Can the S367 eight alpha stencils be recovered as minimal broken labelled
   cycles, rather than as micro-staircase cells?
5. In an LRC-TDA extractor, should we record the largest preterminal directed
   cycle just before peeling kills it?

The headline is:

```text
LRC is not "no protected endpoint cycles."
It is "no arithmetic labelled protected endpoint cycles that also balance measure."
```
