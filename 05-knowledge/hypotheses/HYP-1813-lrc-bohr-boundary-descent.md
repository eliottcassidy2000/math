---
id: HYP-1813
status: EXPLORATORY
source: codex-2026-05-30-S361
related:
  - THM-357
  - THM-358
  - THM-359
  - THM-360
  - HYP-1794
  - HYP-1802
  - HYP-1810
  - HYP-1811
  - HYP-1812
---

# HYP-1813: Lonely Runner Bohr-Boundary Descent

## Statement

Every hypothetical all-protected Lonely Runner endpoint core admits a strict
descent in the finite Bohr-boundary quotient.

More concretely, let `V` be a primitive `k`-speed set, put `n=k+1`, and

```text
Q = n * lcm(V).
```

If the forbidden arcs

```text
||v t|| < 1/n,  v in V
```

cover the circle with every endpoint strictly protected, then the protected
endpoint incidence graph contains one of the following:

1. a forced common divisibility class that can be divided out of the endpoint
   quotient; or
2. a structured progression/subquotient whose induced endpoint graph is
   smaller and still all-protected.

Iterating the descent must eventually expose an unprotected endpoint.  By
THM-357, that endpoint is a boundary lonely witness.

## Why This Is Plausible

The LRC safe time is an anti-approximation time:

```text
||v t|| >= 1/n for every v.
```

Thus a counterexample is not a failure of ordinary Dirichlet approximation.  It
is a finite inhomogeneous Bohr-cover whose boundary has no exposed point.
THM-357 makes the dangerous case exact: full forbidden measure plus all
endpoints protected.

Known tight examples do not resemble this dangerous case.  HYP-1810 says their
unprotected endpoints collapse to the unit-boundary skeleton

```text
{a/n : 1 <= a <= n-1 and gcd(a,n)=1}.
```

So a counterexample must first protect the unit residues.  That protection
requires arithmetic coincidences in `Q`; the descent guess is that those
coincidences cannot close indefinitely.  They either reveal a quotient factor
that should be divided away, or concentrate the endpoint set into a smaller
Bohr/GAP-like substructure where the same endpoint-protection problem repeats.

## Relation To Existing LRC Threads

HYP-1811 is the incidence version: no all-protected endpoint core should exist.

HYP-1812 is the analytic version: an all-protected endpoint core should create
detectable pressure for a nonnegative Fejer/Riesz safe-product kernel.

This hypothesis is the Diophantine approximation version: an all-protected core
should force a smaller finite Bohr-boundary quotient, giving a descent rather
than a one-shot contradiction.

These are three faces of the same proposed obstruction:

```text
peeling leaf       = exposed endpoint in the incidence graph
kernel pressure    = exposed safe mass in Fourier language
Bohr descent       = exposed quotient after denominator refinement
```

## Predictions

1. Endpoint-core peeling in bounded primitive boxes should always return the
   empty core.
2. Artificial nonempty cores, if constructed, should have visible quotient
   factor or generalized-arithmetic-progression structure.
3. Speeds divisible by `n=k+1` should be the first descent class because they
   are the natural protectors of the unit-boundary skeleton.
4. Near-tight positive-gap examples should fail at a definite quotient layer:
   either already at `Z/nZ`, or at a higher denominator introduced by a small
   set of protector speeds.
5. Fejer/Riesz safe-product certificates from HYP-1812 should localize around
   the same quotient layer selected by the descent.

## S362 Formalization

THM-358 proves the first descent model exactly.  For initial-segment speeds
`{1,...,n-1}`, the safe set is precisely

```text
{a/n : gcd(a,n)=1}.
```

This is the equality case of Dirichlet's pigeonhole theorem: unless the
points `0,t,2t,...,(n-1)t` form a regular `n`-gon, some difference gives
`||v t|| < 1/n`.  Thus the unit-boundary skeleton is not just an observed
artifact; it is the rigid equality case of the basic Diophantine lemma.

THM-359 formalizes endpoint/interval core peeling as a greatest-fixed-core
algorithm.  Starting with all forbidden intervals and endpoints, repeatedly
remove endpoints with no remaining protector and then remove intervals whose
boundary has been removed.  The terminal pair is the largest protection core
for that incidence system.

`lonely_runner_bohr_descent_s362.py` implements this exact finite system.  It
verifies THM-358 through `n=36`, checks the known tight examples, and finds
empty terminal cores in all inherited full-measure primitive-box cases.  The
first removed layer in tight examples is always the unit quotient layer
`unit_mod_n`.

THM-360 proves the first quotient filter inside that layer: a unit endpoint
`a/n` can only be protected by a speed divisible by `n`.  Hence every
full-open-cover counterexample must contain at least one speed divisible by
`k+1`.

## Test Plan

1. Extend endpoint-core peeling beyond the inherited S360 primitive boxes.
2. For every peeled layer, record the gcd of endpoint residues, the generated
   subgroup of `Z/QZ`, and the smallest quotient on which the layer remains
   distinguishable.
3. Compare tight examples, near-tight positive-gap examples, and random sets by
   the first quotient layer where unprotected endpoints appear.
4. Search protected endpoint sets for short generalized arithmetic progressions
   and for common-denominator collapses.
5. If a nonempty core appears in artificial covers, quotient by its detected
   structure and check whether protection survives.

## Sources

- THM-357.
- THM-358.
- THM-359.
- THM-360.
- HYP-1802, HYP-1810, HYP-1811, HYP-1812.
- `07-reflections/diophantine-approximation-lonely-runner-s361.md`.
- `04-computation/lonely_runner_bohr_descent_s362.py`.
- `05-knowledge/results/lonely_runner_bohr_descent_s362.out`.
- `04-computation/lonely_runner_endpoint_protection_s360.py`.
- `05-knowledge/results/lonely_runner_endpoint_protection_s360.out`.
