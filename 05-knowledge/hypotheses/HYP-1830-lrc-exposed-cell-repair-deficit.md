# HYP-1830: LRC exposed-cell repair deficit

**Status:** EXPLORATORY; supported by exact S375 cell ledgers for `n=14,15`.

## Claim

Near the scalar-ramp obstruction in the composite Lonely Runner micro-staircase
systems, exposed witness cells cannot be repaired locally without opening a
larger package of new exposed cells.

For `n=14`, the basic half-turn puncture

```text
(i=6, delta=7)
```

opens `56` cells.  The only one-step change that closes all `56` is reverting
to the scalar ramp and hence returns to zero misses.  The best non-reverting
repair move also covers all `56` old misses, but creates `308` total misses.

For `n=15`, the basic `5`-layer puncture

```text
(i=6, delta=5)
```

opens `120` cells.  The best non-reverting repair moves cover only `60` old
misses and create at least `340` total misses.

## Evidence

`04-computation/lonely_runner_feedback_loop_s375.py` cycles through the
user-requested three lanes:

1. `n=14` far-from-scalar vector search.
2. `n=15` forced analogue via the `5`-layer torsion shell.
3. disproof-construction pressure by local speed-set search.

The exact `n=14` half-turn support shell has best missed counts:

```text
support 1:  56
support 2: 112
support 3: 126
support 4: 182
support 5: 168
```

The exact `n=15` `{5,10}`-layer support shell has best missed counts:

```text
support 1: 120
support 2: 220
support 3: 280
support 4: 290
```

Constrained far-from-scalar local search also failed to find blockers.  For
`n=14`, scalar distance at least `3` or `5` still bottomed out at `112` missed
cells.  For `n=15`, scalar distance at least `3` bottomed out at `230`.

The disproof lane did not produce open covers.  The best local `n=14` gated
speed sets retained a positive gap ratio `0.037879` with empty endpoint core;
the best local `n=15` set retained positive gap ratio `0.030303`, also with
empty endpoint core.

## Interpretation

The scalar-puncture moat from HYP-1827 may be a visible face of a stronger
repair-deficit principle: exposed cells have a Hall-like imbalance.  A move can
cover some current misses only by withdrawing protection from a larger,
structured family of cells.

This reframes the next proof target.  Instead of searching for near-scalar
blockers, prove that every non-reverting local repair has positive exposure
defect.  Then combine that with scalar-gauge quotienting and endpoint descent
to force any full blocker either back onto the scalar ramp or into a far-field
region where endpoint cores peel away.

## Test Plan

1. Turn the one-step repair ledgers into exact symbolic inequalities for the
   eight `n=14` alpha stencils.
2. Extend from one-step repairs to two-step repair sets by matching old misses
   against newly exposed cells.
3. Prove the restricted torsion-shell lower bounds for `n=14` half-turn supports
   and `n=15` `{5,10}` supports.
4. Use the repair-deficit score as a branch-and-bound lower bound in the
   normalized `n=14` quotient search.
5. Couple repair deficit to endpoint-protection descent for speed-set
   counterexample pressure.

## Sources

- `04-computation/lonely_runner_feedback_loop_s375.py`.
- `05-knowledge/results/lonely_runner_feedback_loop_s375.out`.
- `07-reflections/lonely-runner-feedback-loop-s375.md`.
- HYP-1818, HYP-1823, and HYP-1827.
