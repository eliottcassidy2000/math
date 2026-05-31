---
id: HYP-1837
status: EXPLORATORY
source: codex-2026-05-31-S380
related:
  - HYP-1823
  - HYP-1828
  - HYP-1831
  - HYP-1832
  - HYP-1833
  - THM-360
---

# HYP-1837: Fourteen-runner composite-denominator disproof attempts pay endpoint-transfer debt

## Statement

At `n=14`, protecting the composite-denominator leak layers that arise from
the scalar-puncture moat can shrink open complement gaps, but it cannot by
itself produce an open-cover counterexample.  The cost appears as increased
endpoint-transfer debt: more unprotected boundary endpoints on descendant
denominators.

Equivalently, the denominator-14 anomaly is real, but its disproof pressure is
balanced by endpoint leakage.

## Evidence

S380 built a target layer from:

- unit endpoints `a/14`;
- known `98/182` leak points such as `9/98`, `29/182`, and `15/182`;
- largest exact gaps of the seven-ladder and single-14-gate seeds.

The resulting layer had `187` rational targets.  The top target-protecting
speeds were mostly multiples of `14`, confirming that THM-360's unit-layer
gate pressure persists into the composite leak layer.

The best new near-disproof was:

```text
speeds=(1,14,28,42,56,70,98,112,126,140,154,168,182)
```

This protects every designed target point and improves the seven-ladder's
gap ratio from `0.005411` to `0.002706`.  However, its exact audit is still
positive-gap:

```text
forbidden_length = 142/143
max_gap = 5/25872
unprotected endpoints = 168
first unprotected = 15/196
```

Across the `90` best exact-audited candidates:

```text
open_cover_candidates = 0
boundary_only_candidates = 1
positive_gap_candidates = 89
```

The only boundary-only object was the initial segment, and it still exposes
the unit boundary.  The composite constructions moved leakage to denominators
such as `2156,2352,2548` and `1078,1176,1274`.

## Why It Matters

HYP-1831 identified scalar moat curvature and quotient-ladder pressure as
recursive LRC metrics.  HYP-1837 says these metrics are coupled: improving one
channel in a counterexample search can worsen endpoint closure.

For the fourteen-runner frontier, this reframes the disproof route.  A speed
set with many `14`-gates can look extremely close to an open cover, but the
same gates create a larger endpoint-transfer obligation.

## Test Plan

1. Build the finite endpoint-protection graph on the unit plus `98/182` leak
   orbit.
2. Search protection graphs first, before assigning speeds.
3. Use the integer protection criterion from THM-360/S360 to determine whether
   any 13-speed set can realize a fully protected graph.
4. Prove a transfer-debt lower bound for the `14`-multiple ladder family
   `(1) union {14q : q != 6}` and its one/two-speed perturbations.

## Sources

- `04-computation/lonely_runner_14_composite_denominator_disproof_s380.py`
- `05-knowledge/results/lonely_runner_14_composite_denominator_disproof_s380.out`
- `07-reflections/lonely-runner-composite-denominator-disproof-s380.md`
- HYP-1828
- HYP-1831
- HYP-1832
- THM-360
