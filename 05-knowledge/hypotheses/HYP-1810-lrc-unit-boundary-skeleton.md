# HYP-1810: LRC tight examples have a unit-boundary skeleton

**Status:** EXPLORATORY; verified on known/scanned tight examples in S359.

## Statement

For a reduced `k`-speed Lonely Runner tight example with threshold
`1/(k+1)`, the unprotected forbidden endpoints are the nonzero unit residues

```text
{a/(k+1) : 1 <= a <= k and gcd(a,k+1)=1}.
```

The point `0` belongs to the standard Dirichlet-pigeonhole orbit, but it is a
forbidden center rather than a forbidden endpoint.  The endpoint skeleton is
therefore the nonzero unit group of `Z/(k+1)Z`.

## Evidence

`lonely_runner_endpoint_protection_s359.py` supports this by exact endpoint
protection records for:

- initial segments with `k=3,4,5,7`;
- sporadic tight examples `{1,3,4,7}`, `{1,3,4,5,9}`,
  `{1,4,5,6,7,11,13}`, and `{1,2,3,4,5,7,12}`.

The stored output lists the unprotected samples; in each listed tight example
they match this unit-boundary skeleton.

THM-358 proves the statement for the initial-segment family
`{1,2,...,k}`.  It is exactly the equality case of Dirichlet's pigeonhole
approximation argument: the `k+1` points `0,t,2t,...,kt` must either have a
gap smaller than `1/(k+1)` or form the regular `(k+1)`-gon.

THM-360 proves the first divisibility filter: a unit endpoint `a/(k+1)` can
only be strictly protected by a speed divisible by `k+1`.  Thus any
full-open-cover counterexample must contain at least one such speed.

`lonely_runner_bohr_descent_s362.py` verifies this theorem computationally
through `n=36` and rechecks the known/scanned tight examples.  In the inherited
full-measure primitive boxes, every case has the unit skeleton.

## S602 Additive-Chain Update

`lrc_p0_collapse_additive_chains_s602.py` isolates the unit skeleton as the
operational `p0`-collapse predicate and confirms that the visible quotient is
larger than AP.  In targeted primitive boxes through the known `n=8` sporadics,
the collapsed rows are:

```text
AP rows for n=4,5,6,7,8;
(1,3,4,7);
(1,3,4,5,9);
(1,2,3,4,5,7,12);
(1,4,5,6,7,11,13).
```

All have boundary witnesses exactly `{a/n : gcd(a,n)=1}` and witness lcm `n`.
All are two-seed addition chains.  Thus the unit-boundary skeleton is not an
AP-only phenomenon; AP is the all-lower member of a larger additive-chain shell
family.  HYP-2153 records the new classification subproblem.

## Why It Matters

This is sharper than saying `1/(k+1)` is a witness.  It says the entire tight
boundary collapses to the unit group of `Z/(k+1)Z`.  A counterexample must
therefore protect this unit-boundary skeleton, which requires at least one
speed divisible by `k+1` and may force a divisibility descent.

## See

- `04-computation/lonely_runner_endpoint_protection_s359.py`
- `05-knowledge/results/lonely_runner_endpoint_protection_s359.out`
- `04-computation/lonely_runner_bohr_descent_s362.py`
- `05-knowledge/results/lonely_runner_bohr_descent_s362.out`
- `04-computation/lrc_p0_collapse_additive_chains_s602.py`
- `05-knowledge/results/lrc_p0_collapse_additive_chains_s602.out`
- `07-reflections/lonely-runner-endpoint-protection-s359.md`
- HYP-2153
- THM-358
- THM-360
