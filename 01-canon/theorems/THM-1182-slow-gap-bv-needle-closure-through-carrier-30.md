---
id: THM-1182
title: Bounded-variation needle closure of every six-comb slow gap through carrier 30
status: PROVED FINITE-EXACT for every integer carrier 1<=a<=30, every phase, and all faster integer speeds.  A rational step probability density gives a strict 1/6 load cap for every faster comb, so even six repeated closed danger combs cannot cover a complete closed slow gap.  The BV tail lemma is proved on paper here but is not yet in Lean; carriers a>=31 remain open
source: codex-2026-07-18-S75 slow-gap exact dual/BV session
depends_on: [THM-1094, THM-1176]
related: [THM-1178, THM-1179, HYP-7715]
script: 04-computation/lrc14_slow_gap_bv_needle_thm1182.py
certificate: 05-knowledge/results/lrc14_slow_gap_bv_needle_thm1182_certificate.json
output: 05-knowledge/results/lrc14_slow_gap_bv_needle_thm1182.out
---

# THM-1182 -- bounded-variation needle closure through carrier 30

At radius `1/14`, use the closed danger comb

```text
Dbar_b={t in R: ||bt||<=1/14}                          (1)
```

and the complete closed slow gap

```text
G_(a,k)=[(14k+1)/(14a),(14k+13)/(14a)].               (2)
```

Closing the danger teeth only makes coverage easier.  The theorem below is
therefore stronger than the strict-open statement needed by LRC.

> **Theorem.**  Let `1<=a<=30` and `k in Z`.  There is an absolutely
> continuous probability measure `mu_(a,k)=f_(a,k) dt`, supported on
> `G_(a,k)`, such that
>
> ```text
> mu_(a,k)(Dbar_b)<1/6                                 (3)
> ```
>
> for every integer `b>a`.  Consequently no six faster integer danger combs
> cover `G_(a,k)`.  Distinctness is unnecessary: the conclusion still holds
> if speeds are repeated.

Indeed, for any six speeds `b_i>a`, the union bound and (3) give

```text
mu(union_i Dbar_(b_i))<6(1/6)=1=mu(G_(a,k)).          (4)
```

Thus the union misses a positive-`mu` subset of the gap.

## 1. The BV high-frequency supplier

Let

```text
chi(x)=1_(||x||<=1/14),       h=chi-1/7.              (5)
```

The mean-zero periodic function `h` has a continuous one-periodic primitive
`H` with

```text
||H||_infinity=3/49.                                  (6)
```

To see the constant, start a primitive at one danger endpoint.  On the
danger arc, of length `1/7`, its slope is `6/7`, so it rises by `6/49`.
On the complementary safe arc it falls by the same amount.  Centering its
range gives (6).

Let `f` be a compactly supported probability density of bounded variation on
the real line, where `TV_R(f)` includes the jumps from zero at the support
endpoints.  Stieltjes integration by parts gives

```text
integral f(t)h(bt)dt=-(1/b) integral H(bt) df(t).      (7)
```

Hence, for every positive integer `b`,

```text
| integral_(Dbar_b) f(t)dt-1/7 |
  <=3 TV_R(f)/(49b).                                  (8)
```

Danger endpoints have Lebesgue measure zero, so (8) is unchanged for open,
closed, or half-open teeth.

## 2. Exact low frequencies plus the analytic tail

For every reflected phase representative in Section 3, the certificate gives
a nonnegative rational step density on `N=200` equal bins of (2).  If the bin
length is

```text
ell=6/(7aN)                                            (9)
```

and the stored integer weights are `w_0,...,w_(N-1)`, then bin `j` has mass

```text
w_j/W,                 W=sum_j w_j,                   (10)
```

and density `w_j/(W ell)`.  Its global variation is exactly

```text
V=[w_0+sum_(j=1)^(N-1)|w_j-w_(j-1)|+w_(N-1)]/(W ell). (11)
```

Each row also chooses

```text
B=Ma,                M in {12,16,20,28}.              (12)
```

Direct rational interval integration verifies

```text
mu(Dbar_b)<1/6                for a<b<=B,              (13)
1/7+3V/[49(B+1)]<1/6.                                (14)
```

For `b>B`, (8) decreases with `b`, so (14) proves (3).  This is the
low-frequency/high-frequency split that turns the computational dual witness
into an all-speed theorem for each certified carrier.

For an independent exact integration formula, put

```text
F(y)=floor(y)/7+min(frac(y),1/7).                      (15)
```

For a bin `[u,v]`,

```text
|[u,v] intersect Dbar_b|
 =[F(bv+1/14)-F(bu+1/14)]/b.                          (16)
```

The independent verifier uses (16), whereas the generator uses explicit
clipped teeth.

## 3. Phase symmetry and completeness

Integer translation gives

```text
G_(a,k+a)=G_(a,k)+1,                                  (17)
```

and every integer comb is invariant under translation by one.  Reflection
`t -> -t` sends

```text
G_(a,k) -> G_(a,-k-1).                                (18)
```

Thus, modulo `a`, phase reflection is `k -> a-1-k`.  The certificate contains

```text
0<=k<=floor((a-1)/2).                                  (19)
```

for a total of `240` rows.  If `a` is odd,
`k=(a-1)/2` is the unique fixed point and is included.  If `a` is even there
is no fixed point.  Therefore (17)--(19) cover every integer phase.

## 4. Exact margins and replay

Of the `240` rows, `213` use `M=12`, `17` use `M=16`, `6` use `M=20`, and
`4` use `M=28`.  The largest certified low-frequency load is

```text
2166600522321/12999999998713
 =1/6-396864787/77999999992278.                        (20)
```

The largest certified tail bound is

```text
3791124999661/22749999996997
 =1/6-3249999031/136499999981982.                      (21)
```

The LP is used only to discover weights.  The stored integers and rational
interval identities are the proof certificate.  From the final repository
paths, a same-script exact replay passed under `python -O` in `176.97 s`; an
independent pure-standard-library replay using (16) passed in `191.26 s` under
concurrent machine load.

Artifact hashes after the final repository replay are recorded in the output
file and below:

```text
generator_sha256 = fc954b992f561783bd2601382afa35c220740bd4c402e559dbe3ab4c1730564b
verifier_sha256  = fac941a6bb011fbb6c65cfdbaeeee5a62df61f6aafc946530b83ba6bf7d85233
certificate_sha256 = 840bb70a2cf7611eec0670a76be7ae8a8b9c506b035ed4253d602df28d226457
certificate_bytes = 581289
```

## 5. Carrier and tournament audit

The faithful proof object is not a runner tournament.  It is the weighted
bipartite incidence operator

```text
A_(b,j)=|Dbar_b intersect I_j|/|I_j|                  (22)
```

between `200` needle-bin vertices and low-frequency proof obligations,
together with the BV tail norm.  Its dual bin weights preserve exactly the
mass obstruction (4); they intentionally discard endpoint-owner chronology.

If a tournament is imposed, use the pairwise observable

```text
O(b,c)=mu(Dbar_b)-mu(Dbar_c),                          (23)
```

orient toward the larger load, and break a tie by the speed label.  On any
selected six speeds this tournament is transitive: score histogram
`0,1,2,3,4,5`, no directed cycles, six singleton SCCs, and one tie Hamiltonian
path.  It preserves only load order.  It destroys the load magnitudes, the
variation budget, and the `1/b` tail in (8), so it cannot prove (4).

The challenged carrier assumption is therefore decisive: the useful carrier
is a dual probability density, or a field of smeared Kakeya needles, rather
than a runner, tooth, or arc.  Runners, slow gaps, wall events, residues,
Fourier modes, and proof obligations were all considered; the bin/obligation
incidence is the quotient that retains the needed predicate.

## 6. Remaining frontier

This closes the entire slow-carrier branch `a<=30`, uniformly over arbitrary
faster integer speeds.  It does **not** prove universal slow-gap noncoverage
for `a>=31`, and therefore does not by itself prove LRC(14).  The finite
certificate is exact.  The analytic BV integration-by-parts lemma above is a
paper proof and has not yet been formalized in Lean.
