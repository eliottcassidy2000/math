---
source: codex-2026-05-30 S357
status: exploratory
tags:
  - lonely-runner
  - quotient-gaps
  - tight-instances
  - endpoint-protection
  - mixed-thresholds
---

# Lonely Runner: Tight Stratum and Endpoint Protection

## Current External Status

As of 2026-05-30, the Lonely Runner Conjecture is not solved in full.  The
current frontier moved quickly in 2025-2026:

- Rosenfeld proved the eight-runner case by finite checking
  (`arXiv:2509.14111`).
- Trakulthongchai proved the nine- and ten-runner cases.
- Sungkawichai and Trakulthongchai proved the eleven-, twelve-, and
  thirteen-runner cases (`arXiv:2604.23906`, submitted 2026-04-26).

In the reduced notation, the problem asks: for every primitive set

```text
V = {v_1, ..., v_k}
```

of distinct nonzero integer speeds, does there exist `t in R/Z` such that

```text
||v_i t|| >= 1/(k+1)  for every i?
```

The 2026 finite-checking papers are very close to this repository's quotient
language: they work with finite ansatz times, improper residue tuples, lifting,
projection, equivalence under coordinate permutation/sign/scaling, and
eventual properness.  Their Section 7 explicitly asks for a uniform ansatz
bound for every non-tight coprime speed tuple.

The newest adjacent paper I found this session, Jensen's `arXiv:2605.27941`
(submitted 2026-05-27), studies mixed thresholds.  It develops safe/unsafe
indicator functions, Fourier expansions, arithmetic-progression summation, and
exact two-speed intersection formulas.  This is the same analytic object as
the repo's gap-measure/safe-product view.

Sources:

- https://arxiv.org/abs/2604.23906
- https://arxiv.org/abs/2605.27941
- https://arxiv.org/abs/2509.14111

## Correction To HYP-1794

HYP-1794 said a lonely witness is either a positive quotient gap or a boundary
residue.  That is true, but it is also mostly topology.

For fixed `V`, define the open forbidden set

```text
F(V) = union_i {t in R/Z : ||v_i t|| < 1/(k+1)}.
```

Because `F(V)` is a finite union of open intervals with rational endpoints,
exactly one of the following happens:

1. `F(V)` has measure `< 1`.  Then there is a positive gap, hence an interval
   of lonely times.
2. `F(V)` has measure `1` but is not all of `R/Z`.  Then every witness is a
   boundary point.  This is the tight stratum.
3. `F(V) = R/Z`.  This is exactly a counterexample.

So the hard target is not the positive-gap/boundary split.  The hard target is
to rule out case 3, or to classify case 2 tightly enough that case 3 cannot
occur.

This reframes HYP-1794:

```text
LRC_k = no primitive k-speed set creates a full open cover of the time circle.
```

## Endpoint-Protection Reformulation

Every forbidden endpoint has the form

```text
t = (m +/- 1/(k+1)) / v_i
```

modulo 1.  Equivalently, all endpoints lie in

```text
Q(V) = (k+1) * lcm(V).
```

Call an endpoint protected if it lies strictly inside another forbidden
interval.  Then:

```text
V is a counterexample
iff
F(V) has measure 1 and every endpoint is protected.
```

If `F(V)` has measure 1 and some endpoint is unprotected, that endpoint is a
lonely boundary witness.  This gives a finite arithmetic obstruction graph:

```text
vertices: forbidden endpoints (i, m, sign)
edge e -> j: endpoint e lies strictly inside the forbidden interval for speed v_j
counterexample certificate: every endpoint has positive protection indegree
```

The membership test is a pure integer inequality.  For

```text
e = ((k+1)m + eps) / ((k+1)v_i),  eps in {-1, +1},
```

speed `v_j` protects `e` exactly when some integer `a` satisfies

```text
| v_j * ((k+1)m + eps) - a * (k+1) * v_i | < v_i.
```

This is the cleanest repo-native next object: endpoint protection is a finite
quotient hypergraph, and a proof of LRC would be a theorem that the
full-measure protection hypergraph cannot have all endpoints protected.

## S357 Exact Scan

I added:

```text
04-computation/lonely_runner_tight_scan_s357.py
05-knowledge/results/lonely_runner_tight_scan_s357.out
```

The scan enumerates primitive speed sets in small boxes and separates:

```text
positive_gap
boundary_only
open_cover
```

It also checks the known sporadic tight examples.

Main exact findings:

- No open-cover counterexample appeared in any scanned primitive box.
- Boundary-only instances are rare:
  - `k=3, max_speed=24`: 1 boundary-only set, namely `(1,2,3)`.
  - `k=4, max_speed=24`: 2 boundary-only sets, `(1,2,3,4)` and `(1,3,4,7)`.
  - `k=5, max_speed=20`: 2 boundary-only sets, `(1,2,3,4,5)` and `(1,3,4,5,9)`.
  - `k=6, max_speed=16`: 1 boundary-only set, `(1,2,3,4,5,6)`.
  - `k=7, max_speed=14`: 3 boundary-only sets, the initial segment and the two
    known `n=8` sporadic examples.
- For all known/scanned tight examples, the first boundary witness is
  `1/(k+1)`, even when the speeds are not an initial segment.  S360 later
  clarified that the full endpoint graph still lives at
  `Q(V)=(k+1)lcm(V)`.

This exactly matches the external tight examples:

```text
n=5: {0,1,3,4,7}
n=6: {0,1,3,4,5,9}
n=8: {0,1,4,5,6,7,11,13}, {0,1,2,3,4,5,7,12}
```

## Simple Divisibility Lemma

Let `n = k+1`.  If no speed in `V` is divisible by `n`, then `t = 1/n` is a
lonely witness:

```text
v_i mod n != 0
=> ||v_i/n|| >= 1/n.
```

Therefore every counterexample to `LRC_k` must contain at least one speed
divisible by `k+1`.

This is elementary, but it is the first quotient filter.  The recent finite
checking papers amplify the same idea: if every residue tuple is eventually
proper modulo a prime `p`, then any counterexample must have a speed divisible
by `p`; enough such primes contradict a finite upper bound on counterexample
size.

## Failed Induction Attempt That Still Points Somewhere

Suppose `n = k+1`, and let `D` be the speeds divisible by `n`.

If `|D|` is large, one can try:

1. divide `D` by `n`;
2. use the lower-dimensional LRC to choose `tau` making those divided speeds
   safe;
3. test the `n` times

```text
t_a = (a + tau) / n,  a = 0, ..., n-1.
```

The divisible speeds stay safe across all `a`.  Each non-divisible speed bans
only a small number of `a` when it is coprime to `n`, and a bounded number
depending on `gcd(v_i,n)` in general.

This gives a quick proof in high-divisibility cases, but not enough for the
general conjecture.  It is essentially the same shape as Jensen's mixed
threshold arithmetic-progression formula and the Barajas-Serra prime filtering
style reductions.

## Best Attack From Here

The positive-gap cases are probably not the problem.  They have an interval
of witnesses, and the recent papers suggest non-tight tuples may admit bounded
finite ansatz witnesses.

The target should be:

```text
Full-measure endpoint-protection theorem.

No primitive speed set V can make F(V) have measure 1 while protecting every
forbidden endpoint.
```

Concrete next steps:

1. Extend `lonely_runner_tight_scan_s357.py` with an endpoint-protection graph.
2. For all scanned tight examples, record the unprotected endpoints and their
   residue classes in `Q(V)`.
3. Search for local protection motifs that would be required in an open-cover
   counterexample.
4. Compare those motifs with Rosenfeld/Sungkawichai/Trakulthongchai improper
   residue tuples: "improper" should be exactly "no ansatz endpoint survives"
   at finite resolution.
5. Translate safe-product Fourier formulas into the repo's Walsh/inclusion
   algebra.  A counterexample is the extreme case where the safe product
   vanishes at every switch point and has zero integral.

I do not see a complete proof today.  The session did narrow the problem:
to solve LRC from this workspace, we should stop treating the continuous
circle as continuous and attack the finite endpoint-protection hypergraph.

## S360 Formalization Update

THM-357 now formalizes the finite-open-cover trichotomy used above.  It proves
that a counterexample is exactly a full-measure forbidden union in which every
endpoint is strictly protected.  The new script
`04-computation/lonely_runner_endpoint_protection_s360.py` builds this graph
exactly and verifies the integer protection criterion.

The main correction to carry forward is terminological: the first witness in
tight examples often appears at `1/(k+1)`, but the full endpoint graph lives at
`Q(V)=(k+1)lcm(V)`.
