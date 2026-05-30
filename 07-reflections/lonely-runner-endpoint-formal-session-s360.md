---
source: codex-2026-05-30-S360
status: formalized-plus-exploratory
tags:
  - lonely-runner
  - endpoint-protection
  - quotient-gaps
  - fejer-kernel
  - riesz-products
  - finite-checking
---

# Lonely Runner Endpoint Protection: Formal Session S360

## What Got Formalized

The most solid local result near the Lonely Runner Conjecture is now
THM-357.  It says that for a speed set `V`, the pulled-back forbidden set

```text
F(V) = union_v {t : ||v t|| < 1/(k+1)}
```

falls into an exact finite trichotomy:

1. `measure(F(V)) < 1`: there is a positive interval of lonely times.
2. `measure(F(V)) = 1` and some forbidden endpoint is unprotected: that
   endpoint is a boundary lonely witness.
3. `F(V)=R/Z`: this is exactly a full-measure all-protected endpoint graph.

The third case is a counterexample certificate.  So the local theorem does
not solve LRC, but it removes the continuous fog.  LRC becomes:

```text
no primitive V produces a full-measure endpoint graph in which every endpoint
is strictly protected.
```

That is the nearest honest formalization of the repo's quotient-gap idea.
THM-355 says empty finite fibers are zero rows and columns in transport.
THM-357 says LRC witnesses are finite quotient gaps too: either open cells
survive, or boundary endpoints survive.  A counterexample must erase both.

## The Endpoint Graph Computation

I added:

```text
04-computation/lonely_runner_endpoint_protection_s360.py
05-knowledge/results/lonely_runner_endpoint_protection_s360.out
```

The script builds the endpoint-protection graph exactly.  For an endpoint

```text
e = (n*m + eps)/(n*v_i),  eps in {-1,+1},
```

it checks protection by `v_j` using both the direct inequality

```text
||v_j e|| < 1/n
```

and the finite integer criterion

```text
|v_j*(n*m + eps) - a*n*v_i| < v_i.
```

The two tests agreed in every checked relation.  The bounded primitive scans
from S357 again found no open-cover candidate.  Every full-measure case was
boundary-only:

```text
k=3, max_speed=24: 1 full-measure case, 0 open covers
k=4, max_speed=24: 2 full-measure cases, 0 open covers
k=5, max_speed=20: 2 full-measure cases, 0 open covers
k=6, max_speed=16: 1 full-measure case, 0 open covers
k=7, max_speed=14: 3 full-measure cases, 0 open covers
```

A useful correction: S357's "quotient collapse to `k+1`" belongs to the
first visible boundary witness and its residue pattern.  The full endpoint
graph lives at

```text
Q(V) = (k+1) * lcm(V),
```

which can be large even for tight examples.

## External Frontier Checked

As of 2026-05-30, the current public frontier is finite and computational:

- Rosenfeld proved the eight-runner case in `arXiv:2509.14111`.
- Trakulthongchai proved nine and ten runners in `arXiv:2511.22427`.
- Sungkawichai and Trakulthongchai proved the eleven-, twelve-, and
  thirteen-runner cases, stated as `k in {10,11,12}` nonzero speeds, in
  `arXiv:2604.23906`.

The nearby analytic frontier is also relevant:

- Tao proved a finite-checking reduction and a small improvement over the
  trivial bound in `arXiv:1701.02048`.
- Bedert uses Riesz products to improve the global loneliness lower bound in
  `arXiv:2511.16636`.
- Jensen studies mixed thresholds and gives Fourier formulas for safe/unsafe
  indicators in `arXiv:2605.27941`.

The repo's endpoint theorem is not competing with these.  It supplies a small
exact object that can talk to both sides: finite endpoint incidence for the
computer-assisted line, and safe-product Fourier kernels for the analytic
line.

## Feedback Loop: Where The New Ideas Landed

First pass: the obvious theorem was the open-cover trichotomy.  That is now
THM-357.  It is elementary, but it is the right kind of elementary: the proof
target becomes finite and checkable.

Second pass: the endpoint graph has two layers that should not be conflated.
The first lonely witness often sits at `1/(k+1)`, but the protection graph's
native quotient is `Q=(k+1)lcm(V)`.  So "collapse" is a residue phenomenon,
not a statement that all endpoints live in the small quotient.

Third pass: the finite-checking papers' sieves and improper tuples may be
viewed as ways of showing that projected endpoint graphs cannot stay
all-protected after enough prime/residue lifts.  In this language, a sieve is
a forced unprotection mechanism.

Fourth pass: the Fejer/Riesz tangent suggests a proof technology.  Define the
closed safe product

```text
P_V(t) = product_v 1_{||v t|| >= 1/(k+1)}.
```

LRC asks whether `P_V` is positive somewhere.  THM-357 says that in the only
dangerous case, all open cells and all endpoints have been erased.  HYP-1812
asks whether a nonnegative Fejer/Riesz kernel on `Z/QZ` can always force

```text
sum_t K(t) P_V(t) > 0.
```

That would turn endpoint protection into a kernel-pressure contradiction.

Fifth pass: root-sign language gives a warning.  In circulant tournaments,
Fejer is the character shadow of one additive root-sign chamber.  For LRC,
endpoint protection is also a chamber system, but its signs live on rational
time residues and speed slabs, not on cyclic pairs `{d,-d}`.  The transfer is
not literal.  The right analogue is probably "which endpoint-protection
chambers have a positive low-frequency safe-product shadow?"

## Concrete Next Moves

1. Add a finite `Z/QZ` safe-product extractor for `P_V`.
2. Search for shifted Fejer kernels and short Riesz products that certify
   scanned tight and near-tight examples.
3. Compare successful kernels with endpoint indegree histograms.
4. Pull Jensen's two-speed intersection formula into exact rational tests for
   endpoint pairs.
5. Try to encode the 2025-2026 sieve language as "forced unprotected endpoint
   after a prime quotient lift."

The live conjectural sentence is small now:

```text
Every full-measure LRC endpoint-protection graph has an unprotected endpoint.
```

The formal theorem says this sentence is exactly where the remaining Lonely
Runner difficulty sits.
