---
id: HYP-1844
status: EXPLORATORY
source: codex-2026-05-31-S385
related:
  - HYP-1828
  - HYP-1831
  - HYP-1832
  - HYP-1837
  - HYP-1842
---

# HYP-1844: LRC quotient protection obeys a debt-export law

## Statement

When a Lonely Runner construction uses divisor gates to protect a lower
quotient layer, it does not simply move monotonically toward a counterexample.
It exports debt to higher endpoint layers.

Informally:

```text
protecting a chosen quotient layer shrinks visible gaps,
but creates new unprotected endpoint obligations on descendant denominators.
```

This is the endpoint-transfer version of the composite-denominator anomaly at
`n=14`.

## Evidence

S380 found that the `14`-multiple ladder protects the designed unit/`98/182`
target layer and halves the seven-ladder gap ratio, but doubles endpoint debt:

```text
n14 seven-ladder:
  gap/th = 0.005411
  unprotected = 84
  low-layer debt = 36

n14 14-ladder debt:
  gap/th = 0.002706
  unprotected = 168
  low-layer debt = 72
```

S385 adds the incidence view:

```text
n14 seven-ladder:
  privateE = 710
  coreE = 0

n14 14-ladder debt:
  privateE = 1420
  coreE = 0
```

So the composite-gate construction improves the horizontal gap statistic while
duplicating exposed endpoint obligations.  The exact one-swap provocations
around both ladder families preserve the same positive-gap and unprotected-
endpoint pattern, suggesting a stable debt channel rather than an accidental
bad choice of one speed.

## Predictions

1. Any family that systematically protects `n`, `n*d`, or scalar-puncture
   leak layers will show an exported endpoint-debt term on descendant
   denominators.
2. For `n=14`, low-layer debt should satisfy a lower bound proportional to the
   number of protected `14/98/182` target orbits unless the construction
   collapses back to the initial segment.
3. Disproof searches should solve the endpoint-debt ledger before speed
   selection; otherwise they will keep finding tiny-gap positive covers.
4. The same export law should appear at `n=15` with order-3 CRT layers, but
   the debt denominator tree will differ from the `2*7` case.

## Test Plan

1. Track endpoint debt by denominator layer in all S373-S385 speed families.
2. Build a transfer matrix from protected low-layer targets to newly exposed
   descendant endpoints.
3. Prove a lower bound for the `14`-multiple ladder and its one/two-swap
   perturbations.
4. Repeat the debt ledger for `n=15` order-3 families to separate universal
   export from `14`-specific behavior.

## Sources

- `04-computation/lonely_runner_14_composite_denominator_disproof_s380.py`
- `05-knowledge/results/lonely_runner_14_composite_denominator_disproof_s380.out`
- `04-computation/lonely_runner_hypothesis_noise_s385.py`
- `05-knowledge/results/lonely_runner_hypothesis_noise_s385.out`
- `07-reflections/lonely-runner-hypothesis-noise-s385.md`
- HYP-1837
