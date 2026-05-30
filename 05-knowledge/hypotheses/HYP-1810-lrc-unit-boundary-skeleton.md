# HYP-1810: LRC tight examples have a unit-boundary skeleton

**Status:** EXPLORATORY; verified on known/scanned tight examples in S359.

## Statement

For a reduced `k`-speed Lonely Runner tight example with threshold
`1/(k+1)`, the unprotected forbidden endpoints are

```text
{0} union {a/(k+1) : gcd(a,k+1)=1}.
```

The nonzero endpoints in this set are exactly the boundary lonely witnesses.

## Evidence

`lonely_runner_endpoint_protection_s359.py` supports this by exact endpoint
protection records for:

- initial segments with `k=3,4,5,7`;
- sporadic tight examples `{1,3,4,7}`, `{1,3,4,5,9}`,
  `{1,4,5,6,7,11,13}`, and `{1,2,3,4,5,7,12}`.

The stored output lists the unprotected samples; in each listed tight example
they match this unit-boundary skeleton.

## Why It Matters

This is sharper than saying `1/(k+1)` is a witness.  It says the entire tight
boundary collapses to the unit group of `Z/(k+1)Z`.  A counterexample must
therefore protect this unit-boundary skeleton, which requires at least one
speed divisible by `k+1` and may force a divisibility descent.

## See

- `04-computation/lonely_runner_endpoint_protection_s359.py`
- `05-knowledge/results/lonely_runner_endpoint_protection_s359.out`
- `07-reflections/lonely-runner-endpoint-protection-s359.md`
