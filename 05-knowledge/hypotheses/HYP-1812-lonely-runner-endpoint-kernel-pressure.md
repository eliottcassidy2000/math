---
id: HYP-1812
status: EXPLORATORY
source: codex-2026-05-30-S360
related:
  - THM-357
  - HYP-1802
  - HYP-1805
  - HYP-1808
  - HYP-1809
  - HYP-1810
  - HYP-1811
---

# HYP-1812: Lonely Runner Endpoint Kernel Pressure

## Statement

The all-protected endpoint certificate from THM-357 should be impossible
because it creates too much local endpoint pressure in a nonnegative
Fourier/Riesz test.

Let `V` be a primitive `k`-speed set, `n=k+1`, and

```text
Q = n * lcm(V).
```

Discretize the time circle at `Z/QZ`, and let

```text
S_v = { t : ||v t|| >= 1/n }
P_V(t) = product_v 1_{S_v}(t)
```

with closed safety at endpoints.  THM-357 says that a counterexample is
exactly a full-measure forbidden cover whose every endpoint is strictly
protected.  The hypothesis is that for every such all-protected endpoint graph
there exists a nonnegative trigonometric kernel `K` on `Z/QZ` such that

```text
sum_t K(t) P_V(t) > 0.
```

This would contradict `P_V=0` and produce a lonely witness.

Equivalently: full endpoint protection should force a spectral pressure
imbalance.  A Fejer- or Riesz-type kernel centered near a minimally protected
endpoint should amplify the closed safe product enough that the all-protected
open cover cannot remain invisible.

## Why This Is Plausible

THM-357 turns LRC into a finite incidence problem.  The recent finite-checking
papers turn LRC into finite ansatz/lifting/sieving problems.  Jensen's mixed
threshold paper supplies Fourier formulas for safe/unsafe indicators, while
Bedert's Riesz-product paper improves global loneliness bounds by building a
weighted product test.  These are all signs that the missing object may be a
nonnegative kernel certificate rather than a raw interval-cover argument.

The Fejer/root-sign work in HYP-1808/HYP-1809 gives the local model: a chamber
of signs can cast a sharply positive character shadow.  For LRC, the signs are
not tournament root signs but endpoint safe/unsafe signs on `Z/QZ`.

## Predictions

1. Boundary-only tight examples have an obvious endpoint-centered kernel:
   center at the first unprotected endpoint, usually `1/(k+1)`.
2. Near-tight positive-gap examples should require wider kernels whose mass
   sees a short interval, not just one endpoint.
3. Any hypothetical all-protected full-measure endpoint graph should exhibit
   abnormal low-frequency pressure in the safe-product Fourier coefficients.
4. Jensen's unequal-threshold formulas should identify the first obstruction
   through two-speed intersection terms before high-order products are needed.
5. The sieve/improper tuples from the 2025-2026 finite-checking papers should
   correspond to endpoint graphs where small kernels fail but larger quotient
   lifts restore pressure.

## Test Plan

1. Build the exact closed safe array `P_V` on `Z/QZ` for scanned tight and
   near-tight examples.
2. Search over Fejer kernels, shifted Fejer kernels, and short Riesz products
   for certificates with `sum K P_V > 0`.
3. Compare successful kernels against endpoint indegree histograms from
   `lonely_runner_endpoint_protection_s360.py`.
4. Force artificial all-protected interval-cover toys, if possible, and see
   what kernel features fail.
5. Translate the strongest certificate into the mixed-threshold Fourier
   language so it can be compared to Jensen's formulas.

## Sources

- THM-357.
- `04-computation/lonely_runner_endpoint_protection_s360.py`.
- `05-knowledge/results/lonely_runner_endpoint_protection_s360.out`.
- Alathea Jensen, `arXiv:2605.27941`.
- Benjamin Bedert, `arXiv:2511.16636`.
- Terence Tao, `arXiv:1701.02048`.
