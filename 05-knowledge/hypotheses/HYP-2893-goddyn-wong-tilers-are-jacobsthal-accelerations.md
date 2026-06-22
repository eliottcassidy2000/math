---
id: HYP-2893
status: STRUCTURAL SIGNAL / family classification route; selected rows exact-verified
source: codex-2026-06-22-S106
tags: [lrc, lrc14, goddyn-wong, tight-tilers, sporadic, jacobsthal, acceleration, exact-tiling, tournament-analysis]
related:
  - THM-560
  - THM-523
  - HYP-2138
  - HYP-+2888
  - HYP-2889
  - HYP-2890
  - OPEN-Q-108
results:
  - 04-computation/lrc_goddyn_wong_sporadic_tilers_codex_s106.py
  - 05-knowledge/results/lrc_goddyn_wong_sporadic_tilers_codex_s106.out
external:
  - https://terrytao.wordpress.com/2017/01/10/some-remarks-on-the-lonely-runner-conjecture/
  - https://arxiv.org/pdf/2409.20160
---

# HYP-2893: Goddyn-Wong tilers are Jacobsthal accelerations of the AP tiler

The Goddyn-Wong LRC14 row

```text
{1,2,3,4,5,6,7,8,9,10,11,13,24}
```

is not an arbitrary sporadic.  It is the `n=13` instance of a tail
acceleration of the AP tiler `{1,...,n}`:

```text
replace v by 2v, where v=12.
```

The structural criterion, quoted in Tao's LRC remarks, is:

```text
{1,...,n} with v replaced by 2v remains tight
if gcd(v,j)>1 for every j in [n-v+1, 2n-2v+1].
```

So the live object is not "non-AP chaos."  It is a Jacobsthal-style interval of
consecutive nonunits modulo `v`.  For LRC14, `n=13`, `v=12`, and the interval
is `[2,3]`; both entries are nonunits modulo `12`.

## Exact audit

Script:

```text
04-computation/lrc_goddyn_wong_sporadic_tilers_codex_s106.py
```

Stored output:

```text
05-knowledge/results/lrc_goddyn_wong_sporadic_tilers_codex_s106.out
```

The script exact-checks selected rows by lower-envelope maximization of

```text
M(S)=max_t min_{s in S} ||s t||.
```

It also computes exact safe measure for moderate selected rows.  Verified:

```text
n=7,  q=8:  {1,2,3,4,5,7,12}              M=1/8,  safe=0
n=13, q=14: {1,2,3,4,5,6,7,8,9,10,11,13,24} M=1/14, safe=0
n=19, q=20: {1..17,19,36}                 M=1/20, safe=0
n=32, q=33: {1..29,31,32,60}              M=1/33, safe=0
n=73, q=74: {1..69,71,73,140,144}         M=1/74
```

The `n=73` row is the first audited multi-acceleration example:

```text
v=70: [n-v+1,2n-2v+1]=[4,7], gcds [2,5,2,7]
v=72: [2,3], gcds [2,3]
```

This matches the longer Goddyn-Wong example listed in Tao's comments:
`...,69,71,73,140,144`.

## Family law

The simplest infinite subfamily is

```text
n ≡ 1 mod 6,     S_n = {1,...,n-2,n,2n-2}.
```

Here the accelerated speed is `v=n-1`, and the nonunit interval is always
`[2,3]`.  The condition is exactly:

```text
2 | (n-1) and 3 | (n-1).
```

For LRC14, `n=13`, this gives `{1,...,11,13,24}`.

The broader family is:

```text
choose tail speeds v>n/2 such that [n-v+1,2n-2v+1] is a run of nonunits mod v;
replace each selected v by 2v.
```

Empirically, simultaneous accelerations can coexist at least in the tested
Goddyn-Wong row `n=73` with `v=70,72`.

## Relation to THM-560

THM-560 proves that the difference-closed exact tilers are exactly AP dilates
`d*{1,...,13}`.  HYP-2893 explains what replaces difference-closure in the
Goddyn-Wong branch:

```text
difference-closed AP tiler:
  all pairwise difference speeds are present, forcing equal spacing.

Goddyn-Wong tiler:
  a missing tail speed v is replaced by 2v because the local obstruction window
  [n-v+1,2n-2v+1] has no units modulo v.
```

So the residual classification is not "all non-difference-closed tilers."  It
is closer to a finite/recursive Jacobsthal-block atlas in the AP tail.  This is
directly useful for LRC14: the known sporadic `{1..11,13,24}` should be routed
as a single accelerated-tail atom, not as a generic wide or additive-energy
anomaly.

## Boundary profile for LRC14

At the q-grid points `a/14`, the immediate blockers are still AP-like:

```text
a= 1: left=[1],  right=[13]
a= 3: left=[5],  right=[9]
a= 5: left=[3],  right=[11]
a= 9: left=[11], right=[3]
a=11: left=[9],  right=[5]
a=13: left=[13], right=[1]
```

The accelerated runner `24` is not one of these first-order blockers.  Its job
is to cover the intervals that would have opened when speed `12` was removed.
That is the key conceptual correction: Goddyn-Wong is a second-order petal
repair, while AP/difference-closed tilers are first-order equal-spacing
systems.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
jacobsthal_window,
single_acceleration,
multi_acceleration,
difference_closed_AP,
arbitrary_sporadic_search.
```

Observable: retained structure for classifying exact tilers.  Gauge orients
toward the carrier with the most proof content and least brute-force search.
Hamiltonian path:

```text
jacobsthal_window
  > single_acceleration
  > multi_acceleration
  > difference_closed_AP
  > arbitrary_sporadic_search.
```

Assumption challenged: exact tilers are not exhausted by AP-like
difference-closure, but they also are not amorphous sporadics.  The useful
quotient preserves the AP tail, the accelerated speed, and the nonunit
interval; it destroys detailed speed magnitude away from the tail and therefore
must be coupled to THM-560 and the LRC14 residual/cap machinery.
