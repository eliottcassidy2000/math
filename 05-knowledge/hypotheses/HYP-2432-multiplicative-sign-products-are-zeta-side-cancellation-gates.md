# HYP-2432 - Multiplicative sign products are the zeta-side cancellation gates

**Status:** OPEN synthesis.
**Source:** codex-2026-06-11-P3.
**Companions:** HYP-2424, HYP-2428, HYP-2431, OPEN-Q-065, T782.

## Statement

The ordinary product atlas

```text
prod (1-q^m)^{b_m}
```

is the power-series shadow of Euler products for zeta and L-functions. Choosing
`b_m` from primes, Mobius, Liouville, or random multiplicative signs converts
multiplicative cancellation into a partition-function gate.

This suggests a zeta-side analogue of the random pentagonal problem:

- eta asks whether a sparse additive support/sign law keeps zeros out of the disk;
- zeta-side products ask which multiplicative sign carriers keep ordinary product
  coefficients from generic dense growth.

## Evidence and Warning

The pilot atlas shows:

- `prime_only` has sparse exponents but almost dense ghost support;
- `mobius` and `liouville` have dense ghost supports but small coefficient maxima
  in the window relative to random all-integer signs;
- scalar slope tournaments are again transitive, so slope alone is too coarse.

This is not a Riemann-hypothesis claim. It is a transfer of method: keep
exponent, ghost, and coefficient carriers separate before making any analytic
claim.

## Next Tests

1. Replace ordinary `q^m` products by Dirichlet truncations `p^{-s}` and compare
   zero/partial-sum statistics.
2. Add completely multiplicative random signs and true character schedules.
3. Build a two-observable Tournament Analysis:

```text
(coefficient slope, ghost irregularity)
```

or

```text
(Dirichlet zero pressure, ordinary coefficient leakage).
```

The expected useful cycles are where multiplicative cancellation and ordinary
partition cancellation disagree.
