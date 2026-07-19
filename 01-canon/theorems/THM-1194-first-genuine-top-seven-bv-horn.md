---
id: THM-1194
title: The first genuine top-seven BV horn at v7=78
status: PROVED FINITE-EXACT.  No ordered thirteen-speed LRC(14) counterexample has v7=78 and v7>=13v6.  Distinctness fixes the lower core to {1,...,6}; an exact rational BV density on G_11(78) rules out coverage by all six faster combs.  The failed a=79 discovery attempt is telemetry only; every v7>=79 case remains open
source: codex-2026-07-18-S77 top-seven target adjustment
depends_on: [THM-1176, THM-1182]
related: [HYP-7715]
script: 04-computation/lrc14_top7_v78_bv_horn_thm1194_verify.py
certificate: 05-knowledge/results/lrc14_top7_v78_bv_horn_thm1194_certificate.json
output: 05-knowledge/results/lrc14_top7_v78_bv_horn_thm1194.out
---

# THM-1194 -- the first genuine top-seven BV horn

Let

```text
0<v_1<...<v_13                                           (1)
```

be distinct integers, and put

```text
M(V)=sup_t min_i ||v_i t||.                              (2)
```

> **Theorem.**  If
>
> ```text
> v_7=78,                    v_7>=13v_6,                 (3)
> ```
>
> then `M(V)>=1/14`.  In particular, this horn contains no LRC(14)
> counterexample.

The theorem is the first direct application of a phase-specific BV needle to
the actual top-seven deletion.  It is a single certified horn, not an interval
of carrier values.

## 1. Distinctness freezes the lower core

Equation (3) gives `v_6<=6`.  Six distinct positive integers lie below or at
`v_6`, so also `v_6>=6`.  Therefore

```text
v_6=6,              {v_1,...,v_6}={1,2,3,4,5,6}.       (4)
```

Every speed in this core is safe at radius `1/14` throughout

```text
I=[11/84,13/84].                                       (5)
```

The exact minimum distances on `I`, in speed order, are

```text
11/84, 11/42, 11/28, 8/21, 19/84, 1/14.               (6)
```

Thus each is at least `1/14`.  Notice that equality for speed six at the
right endpoint is sufficient for the non-strict lonely-runner conclusion.

## 2. One complete 78-slow gap lies inside the core interval

For `a=78` and `k=11`, the complete closed slow gap is

```text
G=G_11(78)=[155/1092,167/1092].                        (7)
```

Since

```text
11/84=143/1092 < 155/1092,
167/1092 < 169/1092=13/84,                              (8)
```

we have `G subset I`, with left and right margins `1/91` and `1/546`.
The carrier itself obeys

```text
||78t||>=1/14                    for every t in G.      (9)
```

## 3. The one-row exact BV certificate

Write

```text
Dbar_b={t: ||bt||<=1/14}.                              (10)
```

The certificate stores a rational probability step density `f` on `200`
equal bins of `G`.  It uses the low-frequency cutoff

```text
B=140*78=10920.                                        (11)
```

Exact rational integration proves

```text
integral_(Dbar_b) f < 1/6          for 79<=b<=10920.    (12)
```

The largest load in (12) occurs at `b=324` and is

```text
13498696201763/80999999992386
 =1/6-651898484/40499999996193.                        (13)
```

The density has total variation

```text
V=2121404250001400/499999999953.                       (14)
```

THM-1182 proves the BV tail inequality

```text
integral_(Dbar_b) f <=1/7+3V/(49b).                    (15)
```

At the first untested speed `b=10921`, the right side is

```text
2123224416495771/12741166665468997
 =1/6-1820166494371/76446999992813982.                 (16)
```

It decreases with `b`.  Hence

```text
integral_(Dbar_b) f<1/6               for every b>78.  (17)
```

The independent verifier recomputes (12) from

```text
F(y)=floor(y)/7+min(frac(y),1/7),                      (18)
|[u,v] intersect Dbar_b|
 =[F(bv+1/14)-F(bu+1/14)]/b,                          (19)
```

rather than replaying the discovery code's clipped teeth.

## 4. Contradiction to a top-seven cover

The six speeds `v_8,...,v_13` are all greater than `78`.  By (17) and the
union bound,

```text
integral_(union_(i=8)^13 Dbar_(v_i)) f < 6(1/6)=1.     (20)
```

Because `f` is a probability density supported on `G`, there is a
`t in G` outside all six closed faster danger combs.  At this time:

- the six lower speeds are safe by `G subset I` and (6);
- speed `v_7=78` is safe by (9);
- the six faster speeds satisfy `||v_i t||>1/14` by construction.

Therefore `min_i ||v_i t||>=1/14`, proving the theorem.  Using closed teeth in
(10) makes the certificate stronger than the strict-open convention needed
for LRC.

## 5. Carrier and tournament audit

As in THM-1182, the faithful carrier is the weighted incidence operator
between the `200` bins of `G` and the frequency obligations `b>78`, together
with the BV tail norm.  The runner-order tournament preserves only speed
order.  A load tournament, oriented by the sign of

```text
integral_(Dbar_b)f-integral_(Dbar_c)f                 (21)
```

and tied by speed label, is transitive: for any selected six speeds it has
score histogram `0,1,2,3,4,5`, no directed cycles, singleton SCCs, and one
tie Hamiltonian path.  It loses the exact load margins and `1/b` tail, so it
cannot replace the BV incidence certificate.

The challenged assumption is again that the carrier must be a runner or an
arc.  Here it is a dual probability density supported inside one slow gap.

## 6. Exact replay and the open next carrier

The independent repository verifier checks the certificate row, all `10842`
low speeds, the variation and tail, the inclusion `G subset I`, and all six
core minima.  Artifact hashes and measured runtime are recorded in the output
file.

The continuation search at

```text
a=79,                      k=11                         (22)
```

did **not** produce a passing BV row.  This is failure of the present
discovery search, not an infeasibility certificate and not a cover witness.
No assertion is made for `v_7=79`, for `80<=v_7<=90`, or for any larger
carrier.  Those horns remain open.
