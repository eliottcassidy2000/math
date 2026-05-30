# HYP-1811: LRC endpoint protection has no all-protected core

**Status:** EXPLORATORY proof-search hypothesis.

## Statement

For every primitive reduced Lonely Runner speed set, the endpoint-protection
incidence graph has an unprotected endpoint after finite peeling.  Equivalently,
there is no nonempty all-protected endpoint core capable of certifying a full
open forbidden cover.

## Evidence

S359 computes endpoint-protection graphs for known tight examples and near-tight
small positive-gap examples.  Known tight examples are not close to all-protected:
their unprotected endpoints form the unit-boundary skeleton.  Near-tight
positive-gap examples can have high protected ratios, but still retain explicit
unprotected endpoints.

## Why It Matters

A counterexample to the Lonely Runner Conjecture is exactly a full open cover
of the time circle, and any full open cover must protect every forbidden
endpoint.  Therefore a no-core/leaf-peeling theorem for endpoint protection
would prove the conjecture.

This mirrors endpoint-transfer collision hypergraphs in the tournament repo:
private endpoints are witnesses; all-protected collisions require an incidence
rank/peeling argument.

## Next Tests

- Implement explicit protection-core peeling.
- Search bounded boxes for nonempty all-protected subcores, not only full
  counterexamples.
- Track which speeds protect unit-boundary endpoints and whether protection
  forces new unprotected endpoints elsewhere.

## See

- `04-computation/lonely_runner_endpoint_protection_s359.py`
- `07-reflections/lonely-runner-endpoint-protection-s359.md`
