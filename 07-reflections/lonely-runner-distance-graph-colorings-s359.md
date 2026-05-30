---
source: codex-2026-05-30 S359
status: exploratory
tags:
  - lonely-runner
  - distance-graphs
  - circular-coloring
  - endpoint-protection
  - quotient-gaps
---

# Lonely Runner as Distance-Graph Coloring

Let `D={d_1,...,d_k}` be a primitive positive speed set.  Form the integer
distance graph

```text
G(Z,D): vertices are integers, and x~y iff |x-y| in D.
```

A time `t in R/Z` gives a regular circular coloring

```text
c_t(x) = tx mod 1.
```

For an edge of distance `d`, the circular separation is `||dt||`.  Therefore
`c_t` is a circular `(k+1)`-coloring exactly when

```text
||d_i t|| >= 1/(k+1)  for every d_i in D.
```

So the reduced Lonely Runner Conjecture is equivalent to:

```text
Every k-distance graph G(Z,D) admits a multiplier circular (k+1)-coloring.
```

This is slightly stronger in language than merely bounding the ordinary
chromatic number.  The coloring is not arbitrary; it is a group homomorphism
from `Z` to the circle.

## Quotient-Gap Translation

For fixed `D`, each distance forbids the open color arc

```text
||d_i t|| < 1/(k+1).
```

The set of bad multipliers is a finite union of rational open intervals in
the multiplier circle.  A good multiplier is a regular circular coloring.
The existing S356/S357 scans split the multiplier circle into:

```text
positive gap     = interval of regular colorings
boundary-only    = tight coloring, all witnesses at forbidden endpoints
open cover       = no regular circular (k+1)-coloring
```

The LRC is exactly the claim that the third case never occurs.

## Endpoint-Protection Coloring Graph

S359 turns the open-cover obstruction into an endpoint incidence graph.  A
forbidden endpoint is protected if another forbidden interval strictly covers
it.  In coloring language:

```text
endpoint       = a tight multiplier where some edge-distance is exactly
                 separated by 1/(k+1)
unprotected    = a valid boundary circular coloring
protected      = another edge-distance destroys that boundary coloring
open cover     = every endpoint is protected
```

Thus the next finite object is not the distance graph `G(Z,D)` itself, but the
boundary graph of multiplier colorings:

```text
vertices: rational endpoint multipliers
incidence: endpoint e is strictly forbidden by distance d_j
counterexample certificate: full measure plus positive protection at every
endpoint
```

The integer test is the one in HYP-1802.  For an endpoint

```text
e = ((k+1)m + eps) / ((k+1)d_i),  eps in {-1,+1},
```

distance `d_j` protects it exactly when some integer `a` satisfies

```text
| d_j * ((k+1)m + eps) - a * (k+1) * d_i | < d_i.
```

That makes the coloring obstruction a finite arithmetic hypergraph in the
quotient `(k+1) lcm(D)`.

## S359 Evidence

I added and ran:

```text
04-computation/lonely_runner_endpoint_protection_s359.py
05-knowledge/results/lonely_runner_endpoint_protection_s359.out
```

The run is exact over `Fraction` arithmetic.  Main observations:

- Known tight examples are not close to counterexamples in the endpoint-core
  sense: they all have unprotected endpoints, hence explicit boundary
  circular colorings.
- The unprotected samples are highly structured.  For the initial and known
  sporadic tight examples they include the expected quotient points such as
  `1/(k+1)`.
- Near-tight positive-gap examples through the default boxes
  `k=4,max=18`, `k=5,max=16`, `k=6,max=14` can have high protection ratios,
  but still leave several unprotected endpoints.
- The best scanned `k=6` near-tight positive example,
  `(1,4,5,6,7,11)`, has `gap/thresh = 0.023810` and
  `protected_ratio = 0.923077`, but still has `5` unprotected endpoints.
- The sporadic `n=8` tight examples have protected ratios about `0.92-0.94`,
  but still exactly `5` unprotected endpoints in this ledger.

## The Clean Proof Target

The distance-graph coloring version of HYP-1802 is:

```text
Endpoint-core obstruction.

For every primitive k-distance set D, if the bad multiplier arcs have full
measure, then at least one forbidden endpoint is unprotected.  Equivalently,
every full-measure regular-coloring obstruction has a boundary circular
(k+1)-coloring.
```

This suggests a peeling proof: repeatedly delete unprotected endpoints or
distances that cannot be part of an all-protected core.  A hypothetical LRC
counterexample would have to survive the peel as a nonempty all-protected
endpoint core.

## Immediate Filters

The trivial coloring multiplier is already useful:

```text
t = 1/(k+1).
```

If no distance in `D` is divisible by `k+1`, then this is a valid circular
`(k+1)`-coloring.  Hence every counterexample must contain a distance
divisible by `k+1`.

In distance-graph terms:

```text
Failure of the residue coloring x -> x/(k+1) mod 1 is mandatory.
```

The external finite-checking strategies that force divisibility by many
primes can be read as iterated failures of such quotient colorings.

## Next Session

1. Add a `--scan-mod-filter` CLI option to S359 for heavier runs outside the
   default stored result.
2. Compute endpoint-core peel depth: repeatedly remove endpoints with
   protection zero and distances whose intervals only protect already-removed
   endpoints.
3. Compare the surviving core, if any, against residue coloring failures
   modulo `k+1` and small prime quotients.
4. Record the distance-graph language in HYP-1802 so future LRC sessions use
   "regular circular coloring" as the default mental model.
