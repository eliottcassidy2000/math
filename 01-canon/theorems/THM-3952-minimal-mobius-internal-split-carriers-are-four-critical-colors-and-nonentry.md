---
id: THM-3952
title: "Minimal Mobius internal-split carriers are four critical colors and do not enter A2"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  In the unit-debt degree-one-ratio specialization of THM-3950, the two
  assigned factors generate the normalization ring k[t] exactly when the
  source point at infinity is one of the four simple critical points of the
  universal cubic ratio map phi. The four resulting affine gauges have an
  explicit linear pencil member. At every other infinity color the best
  degree drop is (3,2), which violates the Abhyankar--Moh degree-divisibility
  theorem for a polynomial embedding of A1. THM-3951 excludes a
  same-function-field Keller atlas for all four critical gauges: at least two
  of the three clean graph/residual incidences remain affine. The C2^4
  Galois-closure ledger has genus one; replacing it by a polynomial cubic
  creates the C3 infinity packet and genus-zero ledger routed by THM-3941.
source: jc-zero-debt-lift / post-THM-3951 minimal-carrier classification, 2026-08-24
audit: >
  SELF-HOSTILE EXACT CANDIDATE. The companion verifies the critical
  polynomial and all four double-fibre factorizations, the four normalized
  representatives and their linear carrier coordinates, the noncritical
  (3,2) hostile cusp, the finite clean-incidence table, and all three
  Riemann--Hurwitz ledgers. The Abhyankar--Moh implication and the use of
  THM-3951 are proved in text and await independent hostile audit.
depends_on:
  - THM-3950-a1-internal-split-denominator-debt-and-equianharmonic-shadow
  - THM-3951-affine-plane-boundary-incidence-forest-and-equianharmonic-survivor-nonentry
related:
  - THM-3941-all-degree-centered-cubic-pole-carrier-routing
  - THM-3947-scalar-weighted-repeated-square-split-trichotomy
external:
  - "Abhyankar--Moh degree divisibility for polynomial embeddings of A1 in A2."
script: 04-computation/jc2_minimal_mobius_internal_split_carriers_thm3952.py
output: 05-knowledge/results/jc2_minimal_mobius_internal_split_carriers_thm3952.out
script_sha256: 4d61c821d69dfe6014c935d3a55750179299eb397806cc6176c406e7545b9be4
output_sha256: fd7580214f925cc64340eb6813d959e3f55e5a041352fa4728b2c8d622513ffa
semantic_sha256: f0cfdf815f49674ba40ed68740fe030abe5ae6fbf6782780cf23153a1b405e57
hash_basis: raw LF bytes
---

# THM-3952 -- the minimal Mobius carrier has exactly four affine colors

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**
Work over an algebraically closed field `k` of characteristic zero. Fix

```text
omega^2+omega+1=0,       delta=omega-omega^2,       delta^2=-3.       (1)
```

This theorem classifies the smallest polynomial carrier inside THM-3950's
nonconstant-ratio packet and then applies THM-3951. Here “minimal” has two
precise meanings:

```text
deg(R/S:P1 -> P1)=1,                    c in k*.                         (2)
```

Higher-degree ratios and nonunit extra debt `c` are not classified here.

## 1. The carrier and its polynomial-embedding gate

Let

```text
R=r1*t+r0,              S=s1*t+s0,
r1*s0-r0*s1 != 0.                                             (3)
```

After absorbing the nonzero constant `c` and the reciprocal scalar `lambda`
into target-coordinate scalings, THM-3950 gives

```text
a=(S-R)(S-omega R)^2,
b=(S+R)(S+omega^2 R)^2.                                  (4)
```

The ratio `x=R/S` is a Mobius coordinate on the projective normalization.
The affine carrier is a polynomially embedded line exactly when

```text
k[a(t),b(t)]=k[t].                                      (5)
```

Put

```text
phi(x)=((1-x)(1-omega*x)^2)/((1+x)(1+omega^2*x)^2).     (6)
```

Then `a/b=phi(R/S)`. The exact Wronskian is

```text
Nphi'(x)Dphi(x)-Nphi(x)Dphi'(x)
 =-2delta*x*(x+omega)*(x-omega^2).                      (7)
```

The fourth critical point is `x=infinity`; all four have ramification index
two. Their fibre factorizations are

```text
phi^(-1)(1):       a-b=R^2(R+delta S),
phi^(-1)(-omega):  a+omega b=(1-omega)S^2(3R+delta S)/3,
phi^(-1)(0):       a=(S-R)(S-omega R)^2,
phi^(-1)(infinity):b=(S+R)(S+omega^2 R)^2.              (8)
```

Thus the critical source colors are

```text
K={0,infinity,-omega,omega^2}.                           (9)
```

They are the four colors of one compactified cover. They remain distinct
affine boundary gauges because an affine change of `t` fixes its point at
infinity.

## 2. Critical infinity is necessary and sufficient

Homogenize `(3)` in `[T:Z]`. The two entries `(a,b)` are homogeneous cubics
in `(T,Z)`. Let

```text
x_infinity=[r1:s1].                                     (10)
```

The unique member of the target pencil that vanishes at `[T:Z]=[1:0]` is
the fibre of `phi` through `x_infinity`. If its multiplicity there is `e`,
then this pencil member contains exactly `Z^e` and has affine degree `3-e`.
There is no triple fibre by `(7)-(8)`.

If `x_infinity` is critical, `e=2`; after an invertible linear target change
the carrier degrees are `(3,1)`. The degree-one member belongs to `k[a,b]`,
so `(5)` follows.

If `x_infinity` is not critical, `e=1`; the normalized degree pair is

```text
(deg f,deg g)=(3,2).                                    (11)
```

Suppose `(5)` nevertheless held. The classical Abhyankar--Moh degree
theorem for an embedding `A1 -> A2` would force `deg f | deg g` or
`deg g | deg f`. Neither `3|2` nor `2|3`, a contradiction. Therefore

```text
k[a,b]=k[t]  iff  x_infinity in K.                      (12)
```

This uses Abhyankar--Moh only in the noncritical direction. The critical
direction is constructive and needs no embedding theorem.

Geometrically, the basepoint-free projective carrier has degree three. It
cannot map three-to-one to a line: that would put a cube in the pencil and
give a triple fibre, excluded by `(7)-(8)`. It is therefore birational onto
a rational cubic. At a critical infinity it has the local orders `(2,3)`,
so its unique genus-paying cusp is at infinity and the affine curve is
`A1`. At a noncritical infinity that point is smooth; the genus-paying
singularity remains affine. This is the projective mechanism behind
`(11)-(12)`.

## 3. The four normalized representatives

Up to an affine change of `t` and a common scaling of `(R,S)`, a Mobius map
with prescribed finite limit `xi` at infinity is `x=xi+1/t`; the infinite
limit is `x=t`. The four critical rows are therefore exactly

| `x_infinity` | `(R,S)` | degree-one member | finite clean colors |
|---|---|---|---|
| `0` | `(1,t)` | `a-b=delta*t+1` | `U=0,V=0` |
| `infinity` | `(t,1)` | `a+omega*b=((3-delta)t+1+delta)/2` | `R=0,U=0,V=0` |
| `-omega` | `(1-omega*t,t)` | `b=delta*t+omega` | `R=0,V=0` |
| `omega^2` | `(1+omega^2*t,t)` | `a=-delta*t-omega^2` | `R=0,U=0` |

Here, as in THM-3951,

```text
U=S+omega^2 R,                    V=S-omega R.            (13)
```

The table is an exhaustive list of affine boundary gauges, not four new
degree-three covers. As a hostile noncritical control, choose

```text
R=t+1,                  S=t,                  x_infinity=1.             (14)
```

Then `(deg a,deg b)=(2,3)`. Both derivatives vanish at

```text
t=(-3+delta)/6,                                           (15)
```

so the missing normalization point is visible as an affine cusp; this row
does not satisfy `(5)`.

## 4. All four minimal carriers fail the affine-plane atlas gate

For every row in the table, `c` is a unit and `R,S` are coprime. Hence

```text
gcd(c,RUV)=1.                                             (16)
```

THM-3951 applies to the associated finite normal cubic surface. Its graph
ramification prime and irreducible equianharmonic residual meet cleanly over
the three source colors

```text
x=0,                  x=-omega,                  x=omega^2,             (17)
```

corresponding respectively to `R=0,U=0,V=0`. Only one source color can be
the single point `t=infinity`. The last column of the table displays the
two or three clean incidences that remain affine in each gauge.

If the cubic function field admitted source coordinates making the target
map Keller, Zariski Main would realize that source as an `A2` open in the
finite normalization. Both ramification primes would be boundary curves,
but two distinct boundary curves of an affine-plane open cannot meet at two
distinct smooth points by THM-3951's incidence-forest lemma. Therefore none
of the four minimal Mobius carriers admits a same-function-field planar
Keller atlas.

## 5. Why the polynomial-carrier transition changes the genus debt

The cover `(6)` has four distinct simple branch values

```text
{0,1,-omega,infinity}.                                  (18)
```

Its Galois group is `S3`, its normal closure has degree six, and every branch
inertia is `C2`. Riemann--Hurwitz gives

```text
2g-2=-12+4*3=0,                    g=1.                 (19)
```

This is THM-3950's equianharmonic shadow. By contrast, a cubic polynomial
with two distinct finite critical points has one totally ramified infinity
and two distinct simple finite critical values. Its `S3` closure has

```text
inertia C3+C2+C2,
2g-2=-12+4+3+3=-2,                 g=0.                 (20)
```

The two finite critical values cannot collide in characteristic zero: after
centering the polynomial as `u^3-3rho^2*u+q`, their difference is
`-4rho^3`, nonzero for distinct critical points.

If the critical points coalesce, the polynomial is affine-equivalent to
`u^3`; because `mu_3 subset k`, its normal closure has degree three, not six,
and the correct cyclic ledger is

```text
inertia C3+C3,
2g-2=-6+2+2=-2,                    g=0.                 (21)
```

Thus, among connected tame degree-three covers, escaping the `C2^4`
elliptic ledger requires creating an index-three inertia point. Moving that
point and its value to infinity is precisely the rational-to-polynomial
carrier transition. THM-3941 then routes the resulting polynomial carrier.
This comparison is a design map, not a construction of a Keller atlas.

## 6. Scope and reproduction

The theorem closes exactly the unit-debt, degree-one-ratio, standard
one-factor internal assignment. It does not classify nonunit `c` sharing a
factor with `RUV`, higher-degree ratios, distributions across several cube
factors, non-`A1` primary branches, or nonmaximal cubic orders. In
particular it does not prove `JC(2)` and constructs no Keller map.

Run

```bash
python3 04-computation/jc2_minimal_mobius_internal_split_carriers_thm3952.py
python3 -O 04-computation/jc2_minimal_mobius_internal_split_carriers_thm3952.py
```

for the exact companion. It includes all four representatives, the
noncritical hostile cusp, the clean-incidence count, and the three
Riemann--Hurwitz controls.
