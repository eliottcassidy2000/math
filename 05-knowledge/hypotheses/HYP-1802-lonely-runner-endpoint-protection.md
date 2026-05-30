---
id: HYP-1802
status: EXPLORATORY
source: codex-2026-05-30 S357
related:
  - HYP-1794
  - THM-355
---

# HYP-1802: Lonely Runner Endpoint-Protection Obstruction

## Statement

Let `V={v_1,...,v_k}` be a primitive integer speed set and let
`delta=1/(k+1)`.  Pull back the open forbidden arcs

```text
||v_i t|| < delta
```

to the time circle.  A counterexample to the reduced Lonely Runner
Conjecture is equivalent to the following finite endpoint-protection
certificate:

1. the forbidden interval union has measure `1`; and
2. every forbidden endpoint is strictly covered by at least one other
   forbidden interval.

The hypothesis is that no primitive `V` admits such a certificate.

Equivalently, whenever the forbidden union has full measure, at least one
endpoint is unprotected, and that endpoint is a boundary lonely witness.

## Distance-Graph Coloring Form

Let `G(Z,V)` be the integer distance graph with vertices `Z` and edges `x~y`
whenever `|x-y| in V`.  A multiplier `t in R/Z` gives the regular circular
coloring

```text
c_t(x) = tx mod 1.
```

This is a circular `(k+1)`-coloring exactly when

```text
||v_i t|| >= 1/(k+1)  for all i.
```

Thus HYP-1802 says that a failure of regular circular `(k+1)`-colorability
would require a full-measure bad-multiplier cover whose every boundary
coloring is strictly killed by another distance.  The conjectural obstruction
is an impossible all-protected endpoint core.

## Finite Quotient Form

Every endpoint lies in

```text
Q(V) = (k+1) * lcm(V).
```

For an endpoint

```text
e = ((k+1)m + eps) / ((k+1)v_i),  eps in {-1,+1},
```

speed `v_j` protects `e` exactly when some integer `a` satisfies

```text
| v_j * ((k+1)m + eps) - a * (k+1) * v_i | < v_i.
```

Thus endpoint protection is a finite arithmetic hypergraph on the boundary
quotient.

## Evidence

`04-computation/lonely_runner_tight_scan_s357.py` exhaustively scanned small
primitive boxes:

```text
k=3, max_speed=24
k=4, max_speed=24
k=5, max_speed=20
k=6, max_speed=16
k=7, max_speed=14
```

No open-cover counterexample appeared.  The boundary-only instances were rare
and exactly recovered the standard small tight examples in the scanned boxes:

```text
(1,2,3)
(1,2,3,4), (1,3,4,7)
(1,2,3,4,5), (1,3,4,5,9)
(1,2,3,4,5,6)
(1,2,3,4,5,6,7),
(1,2,3,4,5,7,12),
(1,4,5,6,7,11,13)
```

In all these examples the first boundary witness is `1/(k+1)` and the
boundary quotient collapses to `k+1`.

`04-computation/lonely_runner_endpoint_protection_s359.py` adds the
endpoint-protection ledger.  In the stored default run, all known tight
examples have unprotected endpoints, and near-tight positive-gap examples in
the boxes `k=4,max=18`, `k=5,max=16`, and `k=6,max=14` still leave explicit
unprotected boundary colorings.  The tight sporadic `n=8` examples have high
protected ratios (`0.92-0.94`) but still have `5` unprotected endpoints.

## Why It Matters

HYP-1794's positive-gap/boundary split is a finite-open-cover trichotomy.  The
nontrivial problem is to rule out the third case: a full open cover of the time
circle.  Endpoint protection is the finite certificate that such a cover would
need.

This aligns with the current finite-checking frontier.  Rosenfeld and
Sungkawichai-Trakulthongchai use finite ansatz times, lifting, projection, and
  eventual properness to force divisibility by many primes.  HYP-1802 asks for a
local boundary obstruction that would make the improper residue tuples
impossible without exhaustive lifting.

## Test Plan

1. Extend the S357 scan to build the endpoint-protection graph.
2. Store unprotected endpoint residue classes for every tight example in small
   boxes.
3. Search for mandatory local motifs in hypothetical all-protected endpoint
   graphs.
4. Compare these motifs against improper residue tuples in the recent
   finite-checking papers.
5. Translate the endpoint-protection obstruction into safe-product Fourier /
   Walsh language.

## Sources

- `04-computation/lonely_runner_tight_scan_s357.py`
- `05-knowledge/results/lonely_runner_tight_scan_s357.out`
- `07-reflections/lonely-runner-tight-stratum-s357.md`
- `04-computation/lonely_runner_endpoint_protection_s359.py`
- `05-knowledge/results/lonely_runner_endpoint_protection_s359.out`
- `07-reflections/lonely-runner-distance-graph-colorings-s359.md`
- `arXiv:2604.23906`
- `arXiv:2605.27941`
