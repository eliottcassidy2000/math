---
id: THM-2283
title: "Mixed rank-two safe-torus floor and scalar rank-three harvest"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every saturated
  rank-two relation torus on one guard and eight ordinary scalar
  coordinates has mixed safe-box Haar mass at least 72/16807. The proof
  classifies every support-at-most-three relation pattern; its only
  zero-union-bound case is repaired by an exact mean-one danger identity
  and one Haar pair. A squared-Fejer comparison, with an elementary
  Machin-series upper bound for pi, then proves that every live scalar
  five-unit/three-blocker cover has relation rank at least three by
  coefficient height 3540. The fixed original-row section has rank at
  least three by height 7080. No scalar profile is excluded, no owner or
  Fourier atom is selected, and LRC(14) remains open.
source: codex-2026-07-25-mixed-safe-torus-rank-three
depends_on:
  - THM-2185-rank-two-safe-cube-floor-and-height-500-rank-three-harvest
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2274-mixed-scalar-relative-rank-harvest-and-adaptive-pair-crossing
related:
  - THM-2270-simultaneous-balanced-cut-relation-and-six-uniform-orientation
  - THM-2278-two-shallow-proper-root-spectrum-and-gap-ancestry-activation
  - THM-2282-thirteen-adic-saturation-and-unit-anchored-minor
script: 04-computation/lrc14_mixed_rank_two_safe_torus_rank_three_thm2283.py
output: 05-knowledge/results/lrc14_mixed_rank_two_safe_torus_rank_three_thm2283.out
script_sha256: d1615934a277a493273beea05e361dbd44cd77371ba98d1b7ae484d9811ea7c9
output_sha256: 60a1d9f52cd7f77fe9f3a26f269c0bf12d43880a6a617de93be68920f8ce0db6
hash_basis: working-tree bytes (LF)
---

# THM-2283 -- mixed safe-torus floor and scalar rank-three harvest

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2185 proves that a saturated rank-two relation torus cannot avoid the
equal-length safe cube. The scalar quotient of the live LRC(14) problem is
not an equal-length cube: its guard-safe interval has mass `5/7`, while
each of the eight ordinary safe intervals has mass `6/7`.

This theorem proves the mixed analogue and applies it directly after
THM-2274:

```text
scalar null safe event
  + scalar relation rank at least two by height 2116
  + every saturated rank-two scalar torus has safe mass at least 72/16807
  + a degree-3540 positive Fourier comparison
  -> scalar relation rank at least three by height 3540.         (1)
```

The new difficulty is a six-coordinate sparse core containing the guard.
Its danger masses sum to one, so the union bound is exactly zero. The
repair is not another termwise Fourier estimate. The mean-one identity

```text
Pr(N=0)=E[(N-1)_+]                                  (2)
```

turns one unavoidable Haar pair into positive safe mass.

## 1. The mixed relation torus

Let

```text
w=(w_0,w_1,...,w_8) in Z^9,
w_i!=0                         for every i,           (3)

Lambda(w)={m in Z^9:m.w=0}.                          (4)
```

Let `L subset Lambda(w)` be saturated of rank two and put

```text
K_L={x in (R/Z)^9:l.x=0 mod 1 for every l in L}.     (5)
```

Saturation gives

```text
ann(K_L)=L.                                          (6)
```

For a coordinate set `T`, character duality therefore gives

```text
ann(pi_T(K_L))=L intersection Z^T.                   (7)
```

In particular, `pi_T(K_L)` is the full coordinate torus exactly when the
right side is zero. No member of `L` has support one, because

```text
c e_i in Lambda(w)  ->  c w_i=0,                    (8)
```

contrary to (3).

On `K_L`, define the guard danger and eight ordinary dangers by

```text
D_0={x:||x_0||<=1/7},             measure(D_0)=2/7,

D_i={x:||x_i||<1/14},             measure(D_i)=1/7,
                                      1<=i<=8.        (9)
```

Endpoint choices are immaterial. Write

```text
N(x)=#{i:x in D_i}.                                  (10)
```

The desired mixed safe box is the atom `{N=0}`.

## 2. The cubic danger minorant

Use the THM-2185 polynomial

```text
P(n)=-(n-1)(n-3)(n-4)/12

    =1-n+(5/6)binomial(n,2)-(1/2)binomial(n,3).      (11)
```

For every integer `0<=n<=9`,

```text
P(n)<=1_(n=0).                                       (12)
```

Indeed, `P(0)=1`, it vanishes at `1,3,4`, and it is negative at `2`
and at every `n>=5`. Hence

```text
measure_(K_L)(N=0)>=E_(K_L) P(N).                    (13)
```

If every projection onto at most three coordinates is Haar, the first
three factorial moments are independent. With

```text
p_0=2/7,             p=1/7,                          (14)
```

their elementary symmetric sums are

```text
sum_i p_i=10/7,

sum_(i<j)p_i p_j
 =8p_0p+binomial(8,2)p^2
 =44/49,

sum_(i<j<k)p_i p_j p_k
 =binomial(8,2)p_0p^2+binomial(8,3)p^3
 =16/49.                                             (15)
```

Substitution in (11) gives the mixed three-wise baseline

```text
B=1-10/7+(5/6)(44/49)-(1/2)(16/49)
 =23/147.                                            (16)
```

We now retain a positive floor when sparse relations spoil some of these
moments.

## 3. Complete sparse-relation classification

Call a nonzero member of `L` **sparse** if its support has size at most
three. Since `rank(L)=2` and support one is impossible, exactly one of the
following cases occurs.

### Case A: no sparse relation

Every projection onto at most three coordinates is Haar by (7), so

```text
measure_(K_L)(N=0)>=23/147.                          (17)
```

### Case B: one rational sparse line of support three

Let its primitive generator have support `A`, with `|A|=3`, and put

```text
q=measure_(K_L)(intersection_(i in A)D_i).           (18)
```

Every pair projection is Haar, and every triple except `A` is Haar.
The triple intersection in (18) is contained in each of its pair
intersections, so in both possible label types

```text
0<=q<=1/49.                                          (19)
```

Only the cubic factorial moment changes.

If `A` contains three ordinary coordinates, its independent triple mass is
`1/343`, and

```text
E P(N)
 =B-(1/2)(q-1/343)
 >=23/147-(1/2)(1/49-1/343)
 =152/1029.                                          (20)
```

If `A` contains the guard and two ordinary coordinates, its independent
triple mass is `2/343`, and

```text
E P(N)
 >=23/147-(1/2)(1/49-2/343)
 =307/2058.                                          (21)
```

### Case C: one rational sparse line of support two

Let `A={i,j}` be its support, let

```text
q=measure(D_i intersection D_j),
delta=q-p_i p_j.                                     (22)
```

For every `k` outside `A`, equation (7) identifies the projection onto
`A union {k}` with

```text
K_(L|A) times R/Z.                                   (23)
```

Thus the pair factorial moment changes by `delta`, and each triple
containing `A` changes by `p_k delta`.

If `A` is an ordinary pair, the sum of the seven outside danger masses is

```text
p_0+6p=8/7.
```

The coefficient of `delta` in (11) is therefore

```text
5/6-(1/2)(8/7)=11/42.                                (24)
```

Since `q>=0`, one has `delta>=-1/49`, and hence

```text
E P(N)>=23/147-(11/42)(1/49)
       =311/2058.                                    (25)
```

If `A` contains the guard, the seven outside coordinates are ordinary.
The correction coefficient and the lower bound for `delta` are

```text
5/6-(1/2)(7/7)=1/3,
delta>=-(2/7)(1/7)=-2/49.                            (26)
```

Therefore

```text
E P(N)>=23/147-(1/3)(2/49)
       =1/7.                                         (27)
```

### Case D: two independent sparse relations

Choose two independent sparse relations `a,b` and put

```text
U=supp(a) union supp(b),          m=|U|.             (28)
```

They span `L` over `Q`, so every relation in `L` is supported on `U`.
As in THM-2185,

```text
3<=m<=6,                                             (29)
```

and the torus factors as

```text
K_L=K_(L|U) times (R/Z)^(9-m).                       (30)
```

If the guard is outside `U`, the union bound on the core and exact
independence on the free coordinates give

```text
measure(N=0)
 >=(1-m/7)(5/7)(6/7)^(8-m).                         (31)
```

For `3<=m<=6`, the smallest displayed value is

```text
(1/7)(5/7)(6/7)^2=180/2401.                         (32)
```

If the guard lies in `U` and `3<=m<=5`, the corresponding floor is

```text
[1-2/7-(m-1)/7](6/7)^(9-m)
 =[(6-m)/7](6/7)^(9-m).                             (33)
```

Its smallest value is

```text
(1/7)(6/7)^4=1296/16807.                            (34)
```

Every bound (17), (20), (21), (25), (27), (32), and (34) is larger than
the final target below.

## 4. The mean-one six-coordinate core

The only unresolved case of Section 3 is

```text
m=6,                 guard in U.                    (35)
```

Let `N_U` be the danger count on the six-coordinate core. Its mean is
exactly

```text
E N_U=2/7+5(1/7)=1.                                  (36)
```

For any nonnegative integer random variable of mean one,

```text
Pr(N_U=0)=E[(N_U-1)_+].                              (37)
```

Here `0<=N_U<=6`, and pointwise

```text
(N_U-1)_+ >=(1/3)binomial(N_U,2).                    (38)
```

At least one pair projection of the core is Haar. Otherwise, for a fixed
core coordinate `i`, equation (7) would put a nonzero relation supported
on `{i,j}` into `L` for every other core coordinate `j`. These five
relations are linearly independent: the `j`th one has a unique nonzero
`j`-coordinate. This contradicts `rank(L)=2`.

The danger intersection of a Haar pair is either

```text
(1/7)^2=1/49
```

or

```text
(2/7)(1/7)=2/49.                                    (39)
```

Consequently

```text
E binomial(N_U,2)
 =sum_(i<j)measure(D_i intersection D_j)
 >=1/49.                                             (40)
```

Equations (37)--(40) give the core floor

```text
Pr(N_U=0)>=1/147.                                    (41)
```

The three free coordinates are ordinary, so (30) finally gives

```text
measure_(K_L)(N=0)
 >=(1/147)(6/7)^3
 =72/16807.                                          (42)
```

Together with Section 3 this proves:

> **Mixed rank-two safe-torus floor.** For every nine-coordinate integer
> row satisfying (3) and every saturated rank-two
> `L subset Lambda(w)`,
>
> ```text
> measure_(K_L)(
>   {||x_0||>1/7}
>   intersection
>   intersection_(i=1)^8{||x_i||>=1/14}
> )>=72/16807.                                       (43)
> ```

The constant is a universal lower bound, not a sharpness claim.

## 5. An exact Jackson certificate at degree 3540

Use THM-2185's normalized squared-Fejer kernel

```text
J_N=F_N^2/integral F_N^2
```

and, for any circle interval `I`, put

```text
q_(N,I)=J_N*1_I.                                     (44)
```

Then

```text
0<=q_(N,I)<=1,
degree(q_(N,I))<=2N-2,
||q_(N,I)-1_I||_1<=eta_N,                           (45)
```

where

```text
eta_N
 =1/2-[4/(pi^2 C_0)]
       sum_(1<=k<=2N-3, k odd) C_k/k^2,              (46)

C_0=N(2N^2+1)/3,

C_k=(4N^3-6Nk^2+2N+3k^3-3k)/6,       0<=k<=N,

C_k=((2N-k)^3-(2N-k))/6,              N<k<=2N-2.     (47)
```

For completeness, (46) is the exact distance-moment estimate. Translation
of a circle interval by `y` changes its indicator in `L^1` by at most
`2||y||`. The Fourier series of `2||y||` is

```text
1/2-(4/pi^2)sum_(k>=1, k odd)cos(2pi k y)/k^2.       (48)
```

Integrating (48) against `J_N`, whose normalized positive coefficients are
`C_k/C_0`, gives (46). Thus the same error bound applies to the guard and
ordinary intervals.

### 5.1 A rational upper bound for pi

The standard Machin identity can be certified without importing a decimal
approximation. Put

```text
a=arctan(1/5),             b=arctan(1/239).          (49)
```

The tangent addition formulas give

```text
tan(2a)=5/12,
tan(4a)=120/119,
tan(4a-b)=1.                                         (50)
```

Direct alternating-series bounds place `4a-b` strictly between `0` and
`1`. The elementary bound `pi>2` puts this interval inside the principal
tangent branch `(0,pi/2)`, so

```text
pi=16a-4b.                                           (51)
```

For `0<x<1`, the alternating arctangent series has an upper partial sum
when its last index is even and a lower partial sum when its last index is
odd. Therefore

```text
a<sum_(j=0)^6 (-1)^j/[5^(2j+1)(2j+1)],

b>sum_(j=0)^1 (-1)^j/[239^(2j+1)(2j+1)].            (52)
```

Equations (51)--(52) give

```text
pi
 <471661273023004128472
   /150134446131591796875

 <104348/33215.                                      (53)
```

The second strict margin is exactly

```text
464759401292/1565687795372314453125>0.               (54)
```

### 5.2 The accepted bandwidth

All `C_k` in (47) are positive. Substituting the strict upper bound (53)
for `pi` in (46) therefore gives a strict rational upper bound for
`eta_N`.

At

```text
N=1771,                  H=2N-2=3540,                (55)
```

the exact finite odd-mode sum gives

```text
eta_1771<23794/100000000.                            (56)
```

Since the two nine-factor product telescopes below cost `18 eta_1771`,

```text
72/16807-18(23794/100000000)
 =424089/420175000000
 >1/1000000
 >0.                                                  (57)
```

The exact scan shows that `N=1770` is not accepted by this rational-upper-
bound ledger. This is a boundary for the displayed certificate, not an
optimality theorem for Jackson smoothing or relation height.

As an independent coarser control, the classical bound

```text
pi<355/113
```

first clears the same ledger at `N=1772`, giving fallback heights `3542`
and `7084`. The companion reconstructs every load-bearing Jackson
coefficient by direct convolution at both bandwidths.

## 6. Scalar rank three

Use the live scalar row and cover from THM-2274:

```text
w_*=(H,q_1,...,q_5,c_1,c_2,c_3) in Z_(>0)^9,        (58)
```

with null mixed safe event

```text
measure(
  {||Ht||>1/7}
  intersection
  intersection_(i=1)^5{||q_i t||>=1/14}
  intersection
  intersection_(j=1)^3{||c_j t||>=1/14}
)=0.                                                  (59)
```

Write

```text
W_K^*
 =span_Q{m in Z^9:m.w_*=0, ||m||_infinity<=K}.       (60)
```

THM-2274 proves

```text
dim_Q W_2116^*>=2.                                   (61)
```

Choose two independent relations from (61), and let `L` be the saturation
of their integer span inside `Z^9`. Then

```text
L subset Lambda(w_*),             rank(L)=2.         (62)
```

Suppose for contradiction that

```text
dim_Q W_3540^*=2.                                    (63)
```

For the guard and eight ordinary intervals in (59), use the degree-`3540`
polynomials from (44). The one-parameter orbit integral of their product
retains exactly the Fourier indices

```text
m in Lambda(w_*),              ||m||_infinity<=3540.
```

By (61)--(63), every such index lies in the rational plane of `L`;
saturation makes it an element of `L`. Conversely every product frequency
in `L` survives Haar integration over `K_L`. Hence the two polynomial
integrals are exactly equal:

```text
I_orbit=I_(K_L).                                      (64)
```

Every coordinate character on both averaging spaces pushes Haar to Haar.
On the orbit this uses `w_i!=0`; on `K_L` it uses the absence of
support-one relations from (8). The nine-factor product telescope in (45)
therefore gives

```text
|I_orbit-0|<=9 eta_1771,                             (65)

|I_(K_L)-measure_(K_L)(mixed safe box)|
 <=9 eta_1771.                                       (66)
```

Equations (42), (64)--(66) would imply

```text
72/16807<=18 eta_1771,                               (67)
```

contradicting (56)--(57). Thus

```text
dim_Q W_3540^*>=3.                                   (68)
```

This conclusion holds for every one of the `165` live scalar profiles,
including the `45` boundary/repeated profiles not covered by the interior
THM-2266 atlas.

## 7. Fixed-section lift

THM-2203 lifts scalar relation coefficients to the fixed nine-coordinate
section of the original thirteen-speed row by

```text
(r_H,r_1,...,r_8)
 |->(2r_H,r_1,...,r_8),                              (69)
```

with zeros on the four forgotten coordinates. The map is integral and
injective, so it preserves three-dimensional relation rank. It doubles
height by at most two. Equation (68) therefore gives:

```text
fixed-section original relation rank>=3
by coefficient height 7080.                          (70)
```

This is quotient-faithful. THM-2185's ambient height-`500` rank-three
theorem alone cannot imply (68), because arbitrary ambient relations need
not vanish on the four coordinates forgotten by the scalar projection.

## 8. Scope and exact referee

The theorem raises the bounded relation rank of every live scalar cover:

```text
THM-2274: scalar rank>=2 by height 2116
THM-2283: scalar rank>=3 by height 3540.              (71)
```

It does not exclude a scalar profile. Rank three does not bound the speed
magnitudes, choose a blocker owner, preserve expiration ancestry, or
identify one exact integer Fourier atom. The torus floor (42) belongs to
the ambient rank-two relation torus; equality of the finite Fourier
integrals is the contradiction device, not an assertion that the original
one-parameter orbit is dense in that torus.

The connection and loss ledger is

```text
source:
  THM-2274's quotient-faithful scalar rank-two packet and the null
  mixed safe event;

target:
  a third bounded relation on the same scalar coordinate section;

map:
  prove a uniform mixed safe-box floor on every saturated rank-two
  torus, then compare its Jackson polynomial integral with the orbit;

preserved:
  all nine scalar labels, the guard/ordinary interval distinction,
  rational relation rank, and the fixed-section lift;

destroyed:
  the particular two starting relations, sparse-core labels, Fourier
  frequency, root sheet, owner, and post-expiration ancestry;

cheapest hostile probes:
  a support-three line containing the guard, a guard support-two line,
  the six-coordinate mean-one core, N=1770 in the sharp rational-pi
  ledger, and the four-coordinate scalar quotient loss;

needed sidecar:
  a residue-nondegenerate anchored rank-three minor and a theorem that
  lands its residue selector on one exact energetic integer frequency. (72)
```

The companion uses exact integer and `Fraction` arithmetic. It:

1. checks the cubic minorant for every count `0,...,9`;
2. verifies all Case A--C floors and the complete Case D union table;
3. checks (36)--(42), including every pointwise count inequality;
4. certifies the tangent identity, alternating-series directions, and
   exact fractions in (53)--(54);
5. reconstructs the squared-Fejer coefficients by triangular convolution;
6. scans the sharp and coarse rational-pi ledgers through their first
   accepted bandwidths; and
7. verifies the heights `3540` and `7080`.

Normal and optimized Python executions reproduce the stored transcript
byte for byte. LRC(14) remains open. QED.
