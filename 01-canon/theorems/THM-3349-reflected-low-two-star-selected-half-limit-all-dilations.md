---
id: THM-3349
title: "Reflected low-two-star selected half-limit at every dilation"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT.  On all 649 upper-median bodies
  remaining behind the robust-edge block, every one of the 11,856 ordered
  primitive rays on which all nine canonical two-star edges are low has a
  deterministic maximum-homogeneous-limit edge which closes at every common
  integer dilation.  An analytic midpoint/shear inequality proves the tail;
  only 42 of 140,711 selected channels require any finite head beyond scale
  two, and 66 exact checks close them.  This includes the current 561-body
  residual atlas, but it is not the whole reflected branch or LRC(14).
source: root/lrc-math-2026-08-12
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
related:
  - THM-3200-fixed-lrc-channel-cleared-overlap-quasipolynomial-and-mass-recurrence-boundary
  - THM-3211-uniform-lrc-channel-limit-bernoulli-cubic-and-sharp-floor
selector_script: 04-computation/lrc14_j7_reflected_low_two_star_limit_selector_scout_20260812.py
selector_output: 05-knowledge/results/lrc14_j7_reflected_low_two_star_limit_selector_scout_20260812.out
selector_script_sha256: cd8b08087f0f7e1e0c0c7d0be673629c7c2702c170c5c1e771e1d76df1d3cd1c
selector_output_sha256: cda9854e4222408a4a3f4e73218e9044126a992187db6047274614a6f01c98c6
script: 04-computation/lrc14_j7_reflected_low_two_star_half_limit_all_dilations_thm3349.py
output: 05-knowledge/results/lrc14_j7_reflected_low_two_star_half_limit_all_dilations_thm3349.out
script_sha256: 6271192c37cdfe6360da993a487005fa648ff1de910539f7e17819464e146b79
output_sha256: 14838fc79bc4b6062bafddd5d36300bbee3e8d11eac5e2ed867026edc03b9d71
hash_basis: LF-normalized bytes
---

# THM-3349 -- selected low-two-star half-limit at every dilation

**PROVED + FINITE-EXACT + VERIFIED-EXACT.**

## 1. Statement

Consider the `649` six-label bodies left after THM-2941's robust-edge block.
At the deterministic upper-median body-safe cell, let the boundary label and
the first other label be the two centres of the canonical nine-edge
two-star.  THM-2941's low-phase graph leaves exactly `11,856` ordered
primitive rays for which every one of those nine edges is low.

On each body/ray assignment, calculate the exact homogeneous common-dilation
limit on the nine edges.  Select the largest limit, breaking exact ties by the
lexicographically largest edge.  If the selected raw levels are

```text
(p,q)=d(P,Q),                       gcd(P,Q)=1,            (1)
```

then, at every common integer scale `s>=2`, its exact reflected overlap
satisfies

```text
I_s >= I_infinity/2.                                      (2)
```

The exact selector census also proves, on all

```text
649*11,856=7,694,544                                      (3)
```

assignments,

```text
I_infinity>debt(1),                  I_1>debt(1).         (4)
```

Every singleton debt term has the form

```text
e/[7(spL-e)].                                             (5)
```

For `s>=2`, its denominator is strictly larger than `s(pL-e)`, so

```text
debt(s)<debt(1)/s<=debt(1)/2<I_infinity/2<=I_s.           (6)
```

Equations `(4)` and `(6)` prove that every selected ray closes at every
integer common dilation `s>=1`.

This theorem includes the current `561` residual bodies.  It closes the
all-nine-low primitive ray bank on those bodies; it does not address an
assignment with a high edge, disconnected-low component scales, arbitrary
reflected levels, physical entry into this sufficient family, or LRC(14).

## 2. Exact indicator and slab midpoint

Fix one selected edge, body ruler `L`, body-safe cell `j`, and labels `e,f`.
Put

```text
R=ej mod L,              S=fj mod L,              h=sd,
chi(x)=1 if ||x||<=1/14 and 0 otherwise,
F(x,t)=chi(Px-(R+et)/L) chi(Qx-(S+ft)/L).                 (7)
```

The exact physical overlap is

```text
I_s=integral_0^1 F(ht,t)dt.                              (8)
```

Write `t=(k+u)/h`.  Unit periodicity in the first coordinate gives

```text
I_s=(1/h) sum_(k=0)^(h-1)
     integral_0^1 F(u,(k+u)/h)du.                        (9)
```

Define the static fibre and its composite midpoint sum by

```text
H(t)=integral_0^1 F(x,t)dx,
J_h=(1/h) sum_(k=0)^(h-1) H((k+1/2)/h).                  (10)
```

We compare `I_s` with `J_h`, then `J_h` with `I_infinity=integral H`.

## 3. Sharp shear bound inside each slab

For the first clause set `epsilon=e/(hL)` and

```text
lambda=1-epsilon/P.
```

Relative to a slab midpoint, the exact clause is `g(v(u))` while the frozen
clause is `g(u)`, where

```text
g(u)=chi(Pu-alpha),
v(u)=1/2+lambda(u-1/2).                                  (11)
```

The map `v` contracts toward `1/2`.  For a jump `c` of `g`, the `u`-set
whose segment from `u` to `v(u)` crosses `c` has length at most

```text
(epsilon/P)|c-1/2|/lambda.                               (12)
```

The `2P` jumps form two translated `P`-grids.  For one translated grid, the
convex, one-periodic distance-to-half sum is maximized at a grid breakpoint,
where it equals

```text
P/4                              if P is even,
(P^2+1)/(4P)                     if P is odd.             (13)
```

Put

```text
gamma_P=1/2                              (P even),
gamma_P=(P^2+1)/(2P^2)                   (P odd).          (14)
```

Summing `(12)` over the two jump grids bounds the first-clause disagreement
by

```text
gamma_P eP/(hLP-e).                                       (15)
```

The second clause similarly contributes `gamma_Q fQ/(hLQ-f)`.  The binary
product inequality `|ab-cd|<=|a-c|+|b-d|` therefore proves

```text
|I_s-J_h| <= gamma_P eP/(hLP-e)+gamma_Q fQ/(hLQ-f).       (16)
```

Endpoint conventions affect only null sets.

## 4. Determinant-controlled midpoint error

Let

```text
C=Qe-Pf,                         D=QR-PS.                 (17)
```

The primitive-fibre law is

```text
H(t)=Phi((-D-Ct)/L),
Phi(z)=[T_A(z)-T_B(z)]/(PQ),
A=(P+Q)/14,                     B=|P-Q|/14,              (18)
```

where `T_a(z)=sum_n(a-|z+n|)_+` is the periodic tent.  Every selected edge
is low, so `P+Q<=7` and `A<=1/2`.  Hence `Phi` is a periodic trapezoid whose
derivative is `0` or `+/-1/(PQ)` and whose derivative has total variation
`4/(PQ)` per circle.

In this atlas `|C|<=14(P+Q)<=98<L`, since all labels are at most `14` and
every ruler is at least `168`.  The path in `(18)` traverses less than one
circle.  Therefore

```text
TV(H')<=4|C|/(LPQ).                                       (19)
```

For a continuous piecewise-affine function, composite midpoint error on a
mesh of width `1/h` is at most `TV(H')/(8h^2)`.  To see this, decompose it
into an affine term plus slope-jump hinges `(t-a)_+`.  Midpoint quadrature is
exact on the affine term, and a hinge with slope jump `delta` contributes at
most `|delta|/(8h^2)`.  Thus

```text
|J_h-I_infinity|<=|C|/(2h^2LPQ).                          (20)
```

Combining `(16)` and `(20)`, with `h=sd`, gives the exact analytic bound

```text
|I_s-I_infinity| <= B_s,                                  (21)
B_s=gamma_P eP/(sdLP-e)
   +gamma_Q fQ/(sdLQ-f)
   +|C|/[2(sd)^2LPQ].                                     (22)
```

All three summands decrease in `s`.  Therefore the first integer `S>=2`
such that `B_S<=I_infinity/2` certifies `(2)` for the whole tail `s>=S`.

## 5. Exact selector census and finite head

The exact selector reconstructs all `(3)` assignments.  The exact limit is
computed after primitive reduction `(1)`; the degenerate case `C=0` is
handled by direct static interval overlap, since the Bernoulli quotient is
invalid there.  Every assignment satisfies `(4)`.  The unique weakest row is

```text
body=(2,4,6,9,11,12), L=5544, j=2970,
ray=(12,6,4,18,3,9), edge=(0,5),
I_infinity=1/3168,
debt(1)=20120188378425768804/99637927188239340409565,
I_infinity-debt(1)
 =35897170405386504838493/
  315652953332342230417501920>0.                          (23)
```

There are `140,711` distinct selected body/edge/raw-channel groups.  Their
tail-threshold histogram is

```text
S       2       3   4  5  6  9
count 140669    29   7  4  1  1.                          (24)
```

Only `42` groups and `66` physical scales lie before their analytic tails.
The exact floor-moment evaluator checks all of them with zero failures.  The
weakest is

```text
body=(2,4,6,9,11,12), edge=(0,5), raw pair=(12,9), s=2,
I_2=174636/553172005,
I_2-I_infinity/2=553321691/3504897823680>0.               (25)
```

The unique maximum threshold `S=9` occurs at

```text
body=(1,2,3,4,6,12), L=168, j=90,
edge=(1,5), raw pair=(6,1), I_infinity=29/1470.           (26)
```

The finite-head semantic digest is
`fd0b48f6975a38ca271106c75f97b9ea9a43d6348ede29289a2438b860d544b0`.

## 6. Hostile control and exact scope

The selector is load-bearing.  At the same canonical upper-median cell for
`H=(1,2,3,4,6,12)`, the low lane `(e,f;P,Q)=(1,3;5,1)` has

```text
I_2=0<I_infinity/2,                 I_infinity=1/5880.    (27)
```

Thus `(2)` is false for arbitrary low channels.  The theorem succeeds
because the exact maximum-limit selector chooses a compatible edge from the
nine-edge two-star on every atlas ray.
