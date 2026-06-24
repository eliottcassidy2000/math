---
id: HYP-2982
title: LRC14 analytic sieve packet weights and Goldbach smoothing guardrails
status: STUB / analytic-number-theory transfer lane reserved; computation pending
source: codex-2026-06-24-S163
related:
  - HYP-2981
  - HYP-2979
  - HYP-2978
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2901
  - HYP-2899
  - HYP-2679
  - HYP-1953
  - HYP-1963
  - THM-572
  - OPEN-Q-108
---

# HYP-2982: LRC14 Analytic Sieve Packet Weights

This stub reserves the prompt's analytic-number-theory lane: sums over primes,
`sum mu(n)`, `sum mu(n)/n`, `sum mu(n)^2/phi(n)`, large-sieve and circle-method
weights, upper-bound quadratic/Selberg sieve packets, exponential sums,
Goldbach smoothing choices, saddle-point / explicit-formula terms, and the
repo's Kaczynski/Kaczorowski boundary/exceptional-set threads.

Working thesis:

```text
The LRC14 proof should treat analytic-sieve weights as labelled packet
certificates, not scalar density estimates.  A smoothing or sieve quotient is
allowed only when it declares its kernel, its exceptional-set boundary, and the
packet labels it forgets.
```

Initial anchors:

- HYP-2899/HYP-2900 already keep denominator Mobius/totient data on the product
  ledger `Div(D) x B_r`.
- HYP-2901 shows no fixed finite denominator basis closes LRC14; prime powers
  and exact-period packets remain load-bearing.
- HYP-2978/HYP-2979 show divisor/Ramanujan packets are useful only when
  endpoint, Haar, route, or state-lift labels are retained.
- HYP-2981 supplies the Fejer interval certificate target, where smoothing
  choices become exact interval proof obligations.
- HYP-2679 and the Kaczynski/Bagemihl boundary-function reflections warn that
  approach regions and boundary values must be carried separately.

Pending computation: build a finite arithmetic atlas for the sums
`M(x)=sum_{n<=x}mu(n)`, `sum mu(n)/n`, `G(x)=sum mu(n)^2/phi(n)`, prime sums,
and a toy Selberg/large-sieve packet hierarchy, then compare the resulting
retention tournament to the current LRC14 packet stack.

