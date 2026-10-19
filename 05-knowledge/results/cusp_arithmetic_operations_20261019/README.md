# Smooth arithmetic operations at the quadratic cusp

19 October 2026. An elementary transfer theorem built on the companion
next-order cusp phase and the prior every-target rough four-form sieve.

## Results

- A universal real-analytic parameter translation implements
  N(epsilon) -> N(epsilon)+q away from a shrinking integer-phase guard.
  Its formula is independent of the escape radius R and its endpoint
  constant C_R; no escape-threshold table is used.
- A real-analytic involution implements the exact complementary sum m
  on the guarded middle-third window. It depends on C_R and uses the
  essential offset m-1.
- The maps have exact partial-action and reflection relations throughout
  their common domains. Their explicit flow is a flow in parameter space,
  not an embedding of the original quadratic recurrence.
- The removed relative measure is O_R(log M/M^2), both for phase measure
  and ordinary parameter measure on comparable-scale windows.
- Bounded arithmetic statistics have an unguarded integral-to-sum error
  O_R(||F||_infinity log M/M). The guarded, calibrated integral divided by
  1-2eta is **exactly** the discrete sum. Positivity equivalence is claimed
  only for that guarded integral.
- For every sufficiently large target m, the prior rough sieve gives
  ordinary epsilon-measure at least c_R/(m^2(log m)^4) with exact
  complementary escape integers and all four neighboring forms 6a±1,
  6(m-a)±1 having at most 37 prime factors each, with multiplicity.

## Files and reproduction

- `cusp_arithmetic_operations_20261019.md`: complete transfer proof,
  exact domains, rounding conventions, arithmetic application, and scope.
- `verify_cusp_arithmetic_operations_20261019.py`: exact rational controls,
  numerical analytic identities, and directed-interval escape controls.
- `cusp_arithmetic_operations_20261019.out`: captured standard run.
- `cusp_arithmetic_operations_20261019.optimized.out`: captured optimized
  run, byte-identical to the standard output.
- `SOURCE-DEPENDENCIES.json`: content hashes of the inherited proof inputs.
- `review/`: independent proof audit, independent verifier, and replay receipts.
- `PUBLICATION_RECEIPT.json`: publication-copy transformations and replay checks.
- `MANIFEST.sha256`: integrity hashes of this packet, excluding itself.

Run:

    python verify_cusp_arithmetic_operations_20261019.py
    python -O verify_cusp_arithmetic_operations_20261019.py

Requires Python 3 and mpmath; tested with Python 3.12.14 and mpmath 1.3.0.
Each run takes about one second here. The script reads no input files,
uses no network, and writes no files. It uses explicit checks rather than
Python assertions, so optimized mode does not remove them.

The verifier passes 2,556 exact rational rounding/guard controls, 26 exact
piecewise integral controls, 84 numerical analytic-map/density controls,
and 120 actual translation plus 45 actual reflection pair controls using
directed interval arithmetic. The latter use the inherited certified
enclosure of C_2 and 150 distinct first-escape evaluations. An exact hostile
example demonstrates why positive unguarded integrals can be false
positives for the intended arithmetic predicate.

## Scope and review

The transfer proof is complete from its explicitly identified inputs.
The [independent audit](review/independent_cusp_operations_review.md)
returned **PASS**. No historical novelty or priority claim is made.
The phase theorem supplies the shrinking error;
the classical sieve supplies the almost-prime abundance. Neither is
replaced by a finite computation here.

K_R and the asymptotic onset are existential, not numerically calibrated.
Finite certified samples do not establish their values. Replacing the
exact phase constant by an approximation requires an additional guard for
its error. The ordinary epsilon measure is weighted and is not claimed to
be uniform arithmetic counting measure.

Twin-prime, balanced Goldbach, and four-prime twin-center integrals are
identified but not estimated from below. No prime conjecture, digit law,
Collatz claim, or global analytic embedding follows from this packet.
