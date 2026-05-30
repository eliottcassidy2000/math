# Lonely Runner Endpoint Protection

**Session:** codex-2026-05-30-S359
**Status:** exploratory progress; exact rational probe added

## Starting Point

The previous Lonely Runner pass reframed the reduced conjecture as a finite
open-cover problem.  For `k` nonzero speeds `V` and threshold `1/(k+1)`, each
speed forbids

```text
F_v = {t in R/Z : ||v t|| < 1/(k+1)}.
```

The conjecture says the open union `F(V)=union_v F_v` is never all of the
circle.

There are three cases:

```text
positive gap       -> open interval of lonely witnesses
boundary-only      -> tight examples; witnesses are endpoints
full open cover    -> counterexample
```

The new object is the endpoint-protection graph.

## Endpoint Protection

Every forbidden interval endpoint has exact rational form

```text
t = (m +/- 1/(k+1))/v.
```

Call an endpoint protected if it lies strictly inside another forbidden
interval.  Then:

```text
full open cover  => every endpoint protected.
```

Conversely, if total forbidden length is `1` and some endpoint is unprotected,
then that endpoint is a boundary lonely witness.  This makes a counterexample
a finite incidence object:

```text
endpoints -> intervals that strictly protect them.
```

The script `04-computation/lonely_runner_endpoint_protection_s359.py` computes
this incidence graph exactly over `Fraction`.

Stored output:

```text
05-knowledge/results/lonely_runner_endpoint_protection_s359.out
```

## Main Observations

### 1. Known tight examples are not close to all-protected

For all standard tight examples checked, there are explicit unprotected
endpoints.  Examples:

```text
(1,2,3,4):          5 unprotected endpoints
(1,3,4,7):          5 unprotected endpoints
(1,3,4,5,9):        3 unprotected endpoints
(1,4,5,6,7,11,13):  5 unprotected endpoints
```

So the known tight stratum is not an almost-counterexample in the protection
graph.  It is a rigid boundary skeleton.

### 2. The boundary skeleton is the unit group modulo `k+1`

In every known tight example checked, the unprotected endpoints are exactly

```text
{0} union {a/(k+1) : gcd(a,k+1)=1}.
```

The lonely witnesses are the nonzero unit residues `a/(k+1)`.  The endpoint
`0` is unprotected but not lonely.

This holds for:

```text
initial k=3,4,5,7
sporadic {1,3,4,7}
sporadic {1,3,4,5,9}
sporadic {1,4,5,6,7,11,13}
sporadic {1,2,3,4,5,7,12}
```

This is much sharper than "the first witness is `1/(k+1)`."  Tight examples
appear to collapse to the quotient `Z/(k+1)Z`, and the surviving boundary
set is the unit group.

### 3. The first counterexample filter is quotient-divisibility

If no speed is divisible by `k+1`, then `t=1/(k+1)` is immediately lonely:

```text
v not == 0 mod (k+1)
=> ||v/(k+1)|| >= 1/(k+1).
```

So every counterexample must contain at least one speed divisible by `k+1`.
This is elementary, but it is now part of the protection-graph story: a
counterexample must first protect the unit-boundary skeleton.

### 4. Near-tight positive-gap examples still have many unprotected endpoints

The near-tight scans rank positive-gap sets by `max_gap/threshold`.
Even the tightest small positive-gap examples retain visible unprotected
endpoint sets.  For example:

```text
k=6, speeds=(1,4,5,6,7,11)
gap/threshold = 0.023810
protected_ratio = 0.923077
unprotected endpoints = 5
```

This is the closest thing in the scan to "almost all protected," and it still
has a small explicit survivor set.

## Emerging Proof Shape

The endpoint-protection formulation suggests a possible proof strategy:

```text
1. Any counterexample must contain a speed divisible by n=k+1.
2. Divisible speeds try to protect the unit-boundary skeleton a/n.
3. Protection of unit endpoints forces additional endpoint exposures elsewhere.
4. Iterate/peel until an unprotected endpoint remains.
```

This is analogous to endpoint-transfer private pivots in the tournament repo:

```text
unprotected endpoint = private witness
all-protected core   = collision block
proof                = leaf-peeling / no-core theorem
```

## New Conjectural Handles

### Unit-boundary skeleton

Tight instances may be exactly those where the protection graph collapses to
the unit-boundary skeleton.  This is not literally a classification yet, but
it is the strongest small-data pattern.

### Protection-peeling

A full open cover would require a nonempty all-protected endpoint core.  The
next finite theorem should try to prove no primitive speed set has such a core,
or at least that every all-protected core forces a smaller counterexample after
dividing a common quotient.

### Divisibility descent

Speeds divisible by `k+1` are necessary for a counterexample.  If many are
divisible, divide them by `k+1` and use lower-dimensional LRC on that subset;
if few are divisible, the non-divisible speeds may be unable to protect all
unit endpoints.  This is the same split as the earlier failed induction, but
now stated in endpoint-protection language.

## Next Experiments

- Optimize endpoint protection so larger mod-filter scans are feasible.
- For every endpoint, record the first protecting speed and build a directed
  protection hypergraph.
- Implement a leaf-peeling algorithm: repeatedly remove unprotected endpoints
  and intervals that only protect removed endpoints.
- Search directly for nonempty all-protected endpoint cores in bounded boxes.
- Test whether boundary-only examples beyond the known ones always have the
  unit-boundary skeleton.

