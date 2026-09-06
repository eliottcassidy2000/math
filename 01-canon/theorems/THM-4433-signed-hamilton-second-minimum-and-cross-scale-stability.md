---
id: THM-4433
title: "Signed Hamilton second minimum and cross-scale stability"
status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED
source: overnight-hexagon-sep05 second research wave
---

# THM-4433 -- Signed Hamilton second minimum and cross-scale stability

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**

For an edge signing of K_n let c_k be the number of negative unoriented
simple k-cycles, and e_(n,k)=(n-2)!/(n-k)! the single-edge value. Work
modulo vertex switching. For odd k the zero class is balanced; for even k
the balanced and antibalanced classes are zero.

## Cross-scale first minimum and uniform stability

For every n>=9 and every3<=k<=n, the first positive value of c_k is exactly
e_(n,k), attained only by a single negative edge, together with its global
negative for even k. The
[all-k proof and exact stratum counts](../../05-knowledge/results/overnight_hexagon_sep05_crossscale.md)
extend the transposition mechanism of THM-4427 without a K8 census premise.

Writing d2(n,k) for the next distinct positive value, for n>=10,

```text
lambda_n <= d2(n,k)/e_(n,k)
 <= min(2-2/(n-2), 2-4(k-3)/((n-2)(n-3))),
lambda_n=2(n^2-9n+22)/((n-2)(n-3)).
```

For n=9 the uniform lower constant is36/35. Thus the normalized second
value tends uniformly to two across all cycle lengths. For0<epsilon<1
and n>3+8/epsilon, a positive value at most(2-epsilon)e_(n,k) is forced
to be a first-minimum signing.

A vertex transposition produces a negative K_(2,r). Its normalized negative
k-cycle count is at least its normalized Hamilton count; the exact excess
retains q=r(n-2-r). This makes the Hamilton scale the worst obstruction.
It does not make the two candidate equality shapes scale-independent:
their difference changes sign at n=2k-3.

## Exact Hamilton second minimum

For every n>=16,

```text
d2(n,n)=2(n-4)(n-3)!.
```

Equality consists exactly of two disjoint negative edges, and their global
negative when n is even. There are3*binom(n,4) labelled switching classes
for odd n and6*binom(n,4) for even n, respectively one and two relabelling
orbits. The second signed-cycle operator gap is twice this value.

The [full classification, matching-integral positivity and equality proof](../../05-knowledge/results/overnight_hexagon_sep05_hamilton_second.md)
are part of this theorem. Transposition rigidity reduces signings to
matchings or one of P3,P4,C3,C4,C5. The exact signed matching sum is

```text
W_n(m)=1/2 sum_j binom(m,j)(-4)^j(n-j-1)!
      =1/2 integral_0^infinity exp(-t)t^(n-m-1)(t-4)^m dt.
```

Its positivity and contraction recurrence order the matching weights;
the five remaining motifs have explicit strict excess. Small hostile
values W_6(3)=-4 and W_5(1)=0 are retained, not extrapolated away.

The [earlier stability proof](../../05-knowledge/results/overnight_hexagon_sep05_hamilton_stability.md),
exact-second15,347-gate audit, and cross-scale33,510-gate audit are mutually
consistent and independently reviewed. Exact second values for every k
at lower orders or at general non-Hamilton scale are not asserted here.
