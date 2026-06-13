---
id: HYP-1843
status: EXPLORATORY
source: codex-2026-05-31-S385
related:
  - HYP-1812
  - HYP-1828
  - HYP-1836
  - HYP-1842
---

# HYP-1843: LRC near-disproof danger is critical-radius surplus, not gap width

## Statement

Visible complement-gap width at the conjectural threshold is a poor measure of
how close a speed set is to a Lonely Runner counterexample.  The better danger
metric is the surplus of the max-min loneliness landscape:

```text
critical_radius / (1/n) - 1
```

together with endpoint-core size.  A tiny open gap can be a steep safe valley,
not a near-counterexample plateau.

## Evidence

S379 already showed that the `n=14` seven-ladder has a tiny visible gap but is
not close in height:

```text
gap/th = 0.005411
exact critical/th = 1.217391
unprotected endpoints = 84
```

S385 found the same behavior in the composite-denominator debt model:

```text
n14 14-ladder debt:
  gap/th = 0.002706
  gap-probe critical/th = 1.208333
  unprotected endpoints = 168
  coreE = 0
```

The `n14 single-gate` example has a larger visible gap but a smaller exact
surplus:

```text
gap/th = 0.045455
exact critical/th = 1.076923
unprotected endpoints = 24
```

So "smaller gap" and "more dangerous counterexample candidate" are different
orders.

## Predictions

1. Ranking disproof candidates by max-gap width alone will repeatedly select
   steep-valley false positives.
2. Serious candidates must drive `critical/th` toward `1` and create a
   nonzero all-protected endpoint core simultaneously.
3. Fejer/Riesz kernel pressure should correlate better with critical-radius
   surplus than with raw threshold-gap width.
4. A useful search objective is a two-coordinate Morse score:

```text
(critical_radius_surplus, endpoint_core_size)
```

with max-gap width only a secondary tie-breaker.

## Test Plan

1. Add exact or bounded critical-radius estimates to all S373-S385 near-
   disproof scans.
2. Compare rank ordering by `gap/th`, `critical/th`, endpoint debt, and peel
   depth.
3. Search for examples with small critical surplus but nontrivial endpoint
   protection; these are more meaningful than tiny-gap ladders.
4. Try to express S385's surplus as a kernel certificate in the style of
   HYP-1812.

## Sources

- `04-computation/lonely_runner_shape_questions_s379.py`
- `05-knowledge/results/lonely_runner_shape_questions_s379.out`
- `04-computation/lonely_runner_hypothesis_noise_s385.py`
- `05-knowledge/results/lonely_runner_hypothesis_noise_s385.out`
- `07-reflections/lonely-runner-hypothesis-noise-s385.md`
- HYP-1812
- HYP-1836
