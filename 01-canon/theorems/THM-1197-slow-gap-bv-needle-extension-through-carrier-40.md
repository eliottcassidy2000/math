---
id: THM-1197
title: The bounded-variation needle closure extends exactly from slow carrier 30 through carrier 40
status: PROVED FINITE-EXACT for every integer carrier 31<=a<=40, every phase, and all faster integer speeds.  Together with THM-1182 this closes every carrier 1<=a<=40.  The original a<=30 certificate is preserved byte-for-byte; this theorem has a separate 180-row certificate and two exact replays.  Carriers a>=41 and a uniform all-a density-selection theorem remain open
source: codex-2026-07-18-S76 BV extension and normalized-orbit audit
depends_on: [THM-1182, THM-1176]
related: [THM-1178, THM-1179, HYP-7715]
script: 04-computation/lrc14_slow_gap_bv_needle_extension_thm1197.py
independent_verifier: 04-computation/lrc14_slow_gap_bv_needle_extension_thm1197_verify.py
certificate: 05-knowledge/results/lrc14_slow_gap_bv_needle_extension_thm1197_certificate.json
output: 05-knowledge/results/lrc14_slow_gap_bv_needle_extension_thm1197.out
---

# THM-1197 -- BV needle closure through carrier 40

Keep the closed danger comb and complete closed slow gap

```text
Dbar_b={t in R: ||bt||<=1/14},
G_(a,k)=[(14k+1)/(14a),(14k+13)/(14a)].               (1)
```

> **Theorem.**  For every integer `31<=a<=40` and every `k in Z`, there is
> an absolutely continuous probability measure `mu_(a,k)`, supported on
> `G_(a,k)`, for which
>
> ```text
> mu_(a,k)(Dbar_b)<1/6                                (2)
> ```
>
> for every integer `b>a`.  Consequently no six faster danger combs cover
> the complete slow gap.  Repeated speeds are allowed.

The union-bound conclusion is immediate:

```text
mu(union_(i=1)^6 Dbar_(b_i))
 <=sum_(i=1)^6 mu(Dbar_(b_i))<1=mu(G_(a,k)).          (3)
```

With THM-1182, this proves the same statement uniformly for every
`1<=a<=40`.  Closing both gap and teeth only strengthens the result needed
for the strict-open LRC formulation.

## 1. Supplier and separate exact certificate

THM-1182 proves the sharp elementary BV estimate

```text
|integral f(t)1_(Dbar_b)(t)dt-1/7|
 <=3 TV_R(f)/(49b).                                   (4)
```

For completeness, the constant is the sup norm of the centered periodic
primitive of `1_(||x||<=1/14)-1/7`.  Its slopes are `6/7` and `-1/7`, its
peak-to-trough range is `6/49`, and hence its centered norm is `3/49`.
Stieltjes integration by parts after the substitution `x=bt` gives (4).
The global variation includes both jumps from the step density to zero at
the support endpoints.  Thus there is no missing boundary term.

The new certificate has `180` rows: one for every reflection representative

```text
31<=a<=40,             0<=k<=floor((a-1)/2).          (5)
```

Each row stores integer masses on `200` equal bins of `G_(a,k)`, chooses a
cutoff `B=Ma`, and checks exactly

```text
mu(Dbar_b)<1/6                    (a<b<=B),
1/7+3 TV_R(f)/(49(B+1))<1/6.                  (6)
```

Equation (4) supplies all `b>B`.  Translation by one and reflection
`k -> a-1-k (mod a)` supply every omitted phase exactly as in THM-1182.

The multiplier ledger is

```text
M       12   16   20   28   44   60   76
rows    77   25   25   34   17    1    1.             (7)
```

The two long rows are `(a,k,M)=(39,5,60)` and `(39,6,76)`.  This only records
which witnesses the stated discovery schedule found; it does **not** assert
that smaller multipliers are mathematically impossible.

An exact argmax replay clarifies the functional role of these long cutoffs.
Their worst low-frequency loads occur already at

```text
(a,k,M)=(39,5,60):  b=178=4*39+22,
(a,k,M)=(39,6,76):  b=273=7*39.                      (8)
```

These are far below the respective cutoffs `2340` and `2964`.  The large
`M` is therefore not evidence that a newly added high comb is dangerous; it
finances the larger BV norm demanded by the low-frequency separation.  In
the unit-gap coordinate, the exact variations are

```text
TV(g_(39,5))=9994269230600/499999999957  =19.988538...,
TV(g_(39,6))=5063307692240/199999999981  =25.316538.... (9)
```

This is the clean `H`-drift functional form: the low block purchases a
variation cost, while the centered primitive charges that cost as
`TV(g)/(14M)` in the asymptotic tail.  Increasing `M` trades a longer exact
ledger for a linearly larger admissible variation budget.

The worst exact low-frequency and tail loads are respectively

```text
18495875250722/110999999989455
 =1/6-8249495041/221999999978910,                    (10)

122774157111133/736749999920431
 =1/6-105057253633/4420499999522586.                 (11)
```

The generator replay uses the clipped-tooth formula from THM-1182.  The
independent pure-standard-library replay instead uses

```text
F(y)=floor(y)/7+min(frac(y),1/7),
|[u,v] intersect Dbar_b|
 =[F(bv+1/14)-F(bu+1/14)]/b.                        (12)
```

Thus the LP is discovery-only: both proof replays consume only stored
integers and rational arithmetic.  The original THM-1182 certificate is not
regenerated or modified.

The completed artifact hashes are

```text
generator_sha256  = c9b29d532c923468d22c974e39c476a1453ef3a501415e655a478994190f77ec
verifier_sha256   = aa4a0c04832fb5c902a77f2898cb75444d5709c412d71611a4c93044e9d7cc46
certificate_sha256 = 0de62d755647f85951c1c36bbb33558ea33b238c0c380d5a7135a4da2b207ba4
certificate_bytes = 437709
generator_replay  = 347.03 s
independent_replay = 502.29 s
```

## 2. The normalized orbit is the underlying object

The fixed-size description clarifies what an all-carrier theorem must see.
Write

```text
t=(k+u)/a,                 J=[1/14,13/14],
b=qa+s,                    0<=s<a,
alpha=s/a,                 phi_s={sk/a}.             (13)
```

Then, modulo an integer,

```text
bt=(q+alpha)u+phi_s,                                  (14)
```

so the exact low-frequency obligation is

```text
||(q+alpha)u+phi_s||<=1/14.                           (15)
```

For fixed `(a,k)`, all phases therefore lie on the cyclic torus orbit

```text
C_(a,k)={(s/a,{sk/a}):0<=s<a} subset (R/Z)^2,        (16)
```

and increasing `q` produces a vertical frequency ladder over the same orbit
point.  This is the precise toothpick self-similarity: the phase stalk is
fixed while each new rung adds one to the frequency.  It is more faithful
than treating the integers `b` as an unstructured list.

If `psi(u)du` is a probability measure on `J` and
`f(t)=a psi(at-k)`, then

```text
TV_R(f)=a TV_R(psi).                                  (17)
```

At cutoff `B=Ma`, the BV tail closes exactly when

```text
TV_R(psi)<7(Ma+1)/(18a)=7M/18+7/(18a).               (18)
```

Equivalently, in the unit-gap coordinate
`x=(7/6)(u-1/14)`, the probability density `g` must satisfy

```text
TV_R(g)<M/3+1/(3a).                                  (19)
```

Equations (15)--(19) are a compact parametric restatement of every row in
both BV certificates.  They expose the remaining all-`a` problem as a
**bounded-variation density-selection problem over cyclic torus orbits**:
select `psi_(a,k)` whose loads stay below `1/6` on the finite vertical ladder
and whose variation grows slowly enough for (18).  A density depending only
on the real ratio `b/a` discards `phi_s` and is not a faithful quotient.

## 3. Tournament and challenged-carrier audit

The proof object is the weighted bipartite incidence matrix between needle
bins and obligations `(q,s)` in (15), augmented by the BV norm.  The cyclic
orbit points, rather than runners or arcs, are the natural vertices on the
obligation side.  This quotient preserves the exact predicate (2), phase
arithmetic, the ladder, and the tail budget.  It discards endpoint-owner
chronology, which this mass proof does not use.

For a tournament diagnostic, orient two obligations toward the one with
larger certified load and break equality lexicographically by `(q,s)`.  This
gauge is a total order: for `m` selected obligations its score histogram is
`0,1,...,m-1`, it has no directed cycle, every SCC is a singleton, and it has
one Hamiltonian path.  The tournament preserves load order only.  It destroys
the rational margins, cyclic phase, and variation functional, so no
tournament fingerprint alone implies (3).

The challenged assumption is that a tournament vertex should be a runner,
gap, or tooth.  Here the useful carriers are (i) a cyclic phase-frequency
obligation and (ii) a smeared Kakeya needle bin.  Fixed circle sections,
section boundaries, wall events, residues, Fourier modes, and endpoint
owners were also considered; none retain both (15) and (18) as economically.

## 4. Honest frontier

This theorem is an exact ten-carrier extension, not an all-carrier proof.
Carriers `a>=41` remain open.  More importantly, the row-specific witnesses
do not yet give a uniform density-selection law on the cyclic orbit (16).
The analytic estimate (4) is still a paper lemma rather than a Lean theorem.
The next structural target is to control orbit arithmetic -- especially
large common divisors and short phase cycles -- while bounding the variation
needed to clear each vertical ladder.

The final repository artifact hashes are

```text
generator    c9b29d532c923468d22c974e39c476a1453ef3a501415e655a478994190f77ec
independent  aa4a0c04832fb5c902a77f2898cb75444d5709c412d71611a4c93044e9d7cc46
certificate  0de62d755647f85951c1c36bbb33558ea33b238c0c380d5a7135a4da2b207ba4
```
