---
id: THM-1100
title: THE EXTENDED RATIONAL-POINT SIEVE — exact residue criterion and useful finite atlas, but the proposed uniform primitive-family route is refuted by single-coordinate divisor loading
status: CORRECTED / UNIFORM ROUTE REFUTED by the pre-existing THM-566 and HYP-2876; the residue criterion, named-family certificates, and reported sample kill rates remain exact. THM-1098 strengthens the obstruction to actually-lonely primitive covering rows and proves the necessary logarithmic height cost
source: opus-2026-07-17-S373 (owner: work the surviving tools toward a new route)
depends_on: [THM-1035/1040 (the classical seven-moduli sieve, the q ≤ 14 case), THM-1050 (common-dilation invariance), THM-566 (primitive covering rows defeat every fixed denominator cap), HYP-2876 (arbitrary finite denominator bases are atlases, not closures), THM-1098 (explicit lonely obstruction and height cost), THM-1105 (first unblocked modulus / arithmetic-position experiment), THM-1110 (sharp forbidden window and gcd-stratified correction), THM-1095 (the retired ledger route), MISTAKE-154/156 (sampling guardrails)]
scripts: 04-computation/extended_sieve_opus_S373.py, sieve_adversarial_opus_S373.py, largest_gap_opus_S373.py -> 05-knowledge/results/
---

# THM-1100 — an exact finite atlas, not a uniform route

## Correction after the canon audit

The proposed final statement in the original version of this note was already
refuted by THM-566, proved on 2026-06-22.  For every `B`, that theorem constructs
the primitive covering row

```text
{1,2,...,11,13,84*lcm(1,...,B)}
```

whose last runner is at the observer at every rational point of reduced
denominator at most `B`.  HYP-2876 gives the same obstruction for an arbitrary
fixed finite list of denominators.  Thus common-gcd normalization does **not**
rescue a uniform cap: it removes common dilation, but it does not remove the
divisibility loaded into one coordinate.

THM-1098 sharpens the correction.  The obstructing rows can be chosen primitive,
covering, and provably lonely, with an explicit later rational witness.  At
height `H=lcm(1,...,B)` they force `q_min>B`, hence an unavoidable
`(1-o(1)) log H` lower cost.  It also proves rigorously that a fixed rational
address can lie in lonely components whose lengths tend to zero.

Everything below about the modular criterion and the displayed finite samples
is retained.  What is withdrawn is the inference that those samples support an
absolute `Q0` on primitive families.  Opus-S374 independently caught the same
regression and recorded in THM-1105 the useful sampled position law for `q0`,
the first modulus dividing no speed.  Its reported `96.8%` agreement between
`q0` and the least witness denominator is sample evidence, not a uniform law.

THM-1095 retired the uniform Bonferroni ledger. This file builds a route
from what survived, centred on the tool with the best kill rate.

## The observation the classical sieve under-uses

t = p/q is lonely iff, for every speed v,

> **min(vp mod q, q − (vp mod q)) · 14 ≥ q**

This depends **only on v mod q**. So "does a lonely rational of denominator
≤ Q exist" is a condition on the speeds mod lcm(1..Q) — **finite and
arithmetic**, not analytic. That is a different kind of problem from the
uniformity that blocked the ledger, which is the point of the pivot.

The classical seven-moduli sieve is exactly the case q ≤ 14, where the
band degenerates to "q divides no speed" (for q ≤ 14 every nonzero residue
passes). For q > 14 the band is a genuine but WIDE constraint — a random
speed fails with probability only ~1/7 — and there are φ(q) numerators to
try per q.

## The kill rate

| Q | random 13-families | the hard stratum (blocks all seven) |
|---|---|---|
| 14 | 86.5% | 0/54 |
| 20 | 99.5% | 52/54 |
| **30** | **100%** | **54/54** |

Every named adversarial family falls, and cheaply: the THM-1055 primitive
failure at 7/15, the THM-1060 L=31 family at 3/23, the AP d=8 at 1/2, the
sum-free {1,3,…,25} at 1/2, and {1,…,13} at 1/14 — the last with equality,
which is exactly right for the extremal family.

## Common dilation is only the first obstruction

`k·{1,…,13}` needs denominator `14, 28, 14, 14, 98` for
`k = 1,2,3,5,7`.  Thus common dilation inflates this example only along the
`2,7`-part of `k`, and THM-1050 correctly reduces common dilation to primitive
representatives.

The original inference that the route should therefore work on primitive
families was false.  THM-566's fixed core contains `1`, so it is primitive for
every `B`, while one far runner is divisible by every denominator through `B`.
Primitivity controls `gcd(S)`; bounded-denominator loneliness is controlled by
the much richer collection of divisor packets carried by the individual
coordinates.

## What the sampled search established

The minimal denominator over primitive families measured:

| method | max found |
|---|---|
| 600 random families | 25 |
| light adversarial hill-climb | 32 |
| harder hill-climb (45 restarts, 700 steps) | **39** |

Each increase in search effort raised the sampled maximum.  The exact data are
useful as a description of that bank, but the deterministic lcm family from
THM-566 was outside the sampling distribution and has arbitrarily large least
denominator.  In particular, “nothing sampled exceeds `39`” has no global
content.  The `q = 39` row

> V = {31, 32, 36, 74, 145, 210, 231, 260, 304, 459, 500, 552, 616}

was verified exhaustively: no lonely rational with `q < 39`, and `gcd(V)=1`.
That per-row certificate remains valid.

## The largest-gap reformulation, now rigorously separated

The natural control is “the largest uncovered gap is bounded below by
`L0>0`”, since an interval of controlled length meets a bounded rational grid.
The S373 measurements already warned against it: largest gaps were around
`0.001` while actual least denominators were around `39`.

THM-1098 makes the separation exact.  A runner of speed `M` cuts every lonely
component down to length at most `6/(7M)`.  Nevertheless the primitive covering
families

```text
{1,...,11,13,84(41k+1)}
```

all have the same lonely point `17/41`, while their largest component lengths
tend to zero.  **The gaps are situated, not large:** arithmetic address and
geometric thickness are independent coordinates.  The old sample is retained
as discovery evidence, but no longer bears the proof claim.

## Correct surviving route

The absolute-`Q0` conjecture is **refuted**, even for primitive covering rows
that are independently known to be lonely.  For a *fixed* `Q`, the residue
test modulo `lcm(1,...,Q)` remains exact and useful; it is an atlas, not a
closure.

A correct sieve theorem must be adaptive or height-dependent.  At minimum it
must account for the denominator-loading depth

```text
kappa(S) = max {B : every q <= B divides at least one speed of S},
```

Every rational witness has reduced denominator strictly larger than
`kappa(S)`.  THM-1098 shows `kappa` can be `(1-o(1)) log(max S)` on primitive
covering rows.  After deleting those divisibility-killed packets, the open
arithmetic task is to prove that some remaining numerator band is live.  For
`q>14`, “no speed is divisible by `q`” is necessary but no longer sufficient;
the thirteen danger bands can still cover all numerator residues.  This is the
precise residual that the finite residue computations measure.  THM-1110 now
gives the sharp one-modulus language: its `q=15` row refutes the extended
lemma, and the correct numerator union bound is stratified by `gcd(v,q)`.
At `q=90`, three gcd-five residues already cover all unit numerators, so even
the proposed eleven-speed unit count cannot be applied to arbitrary nonzero
residues.  The remaining content is necessarily cross-modulus compatibility,
not another single-`q` cardinality estimate.

Equivalently, `kappa(S)=q0(S)-1` for THM-1105's first modulus `q0` dividing
no speed.
