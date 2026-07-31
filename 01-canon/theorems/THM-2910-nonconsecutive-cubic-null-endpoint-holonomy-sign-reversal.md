---
id: THM-2910
title: "Nonconsecutive cubic-null endpoint-holonomy sign reversal"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The ordered
  four-slot block charts with supports
  (0,1,2,5) and (0,1,2,6) each contain exactly one positive plane on
  which the binary quadratic moment divides the binary cubic moment.
  Both planes retain atomic adjacent-tensor positivity and strictly
  positive local response-row TP3, but their quartic endpoint
  determinants have opposite signs.  Thus even simultaneous cubic
  nullity does not admit a support-uniform endpoint orientation.  An
  independently exact translated pair (2,3,4,8)/(3,4,5,9) shows that
  even the fixed gap pattern (1,1,4) does not supply the orientation.
source: root/nonconsecutive-cubic-null-sign-reversal-2026-07-29
audit: >
  An independent hostile audit reconstructed both factorial moment systems,
  cubic branch selectors, quotient-ring endpoint remainders and four TP3
  polynomials; verified the a0=0 coefficient-algebra scope; replayed normal,
  optimized and stored output byte-for-byte; and matched all LF hashes.
  A separate exact translated-support audit verified the same-gap
  (2,3,4,8)/(3,4,5,9) sign reversal by rational isolation and modular
  coprimality.
depends_on:
  - THM-2853-gamma-adjacent-tensor-cycle-weighted-positivity
  - THM-2872-four-slot-shared-multipole-quartic-norm-and-response-secant-reduction
related:
  - THM-2879-all-shift-cubic-null-endpoint-holonomy-exit
  - THM-2906-atomic-tp3-does-not-orient-mixed-endpoint-holonomy
  - THM-2909-nonconsecutive-response-row-tp3-arc
script: 04-computation/gmc_nonconsecutive_cubic_null_holonomy_reversal_thm2910.py
output: 05-knowledge/results/gmc_nonconsecutive_cubic_null_holonomy_reversal_thm2910.out
script_sha256: bdc738fb0b964d502596784ec3b2d747b18572bfbfb70186b52c254cd948e1bf
output_sha256: bc3c597fb0982ffd0dd8552babe64528b3da0cef9315cf2ba000a728443326c3
hash_basis: LF-normalized bytes
---

# THM-2910 -- nonconsecutive cubic-null holonomy reverses sign

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Put

```text
L(s^k)=k!,                    f_j=s^j/j!,
d_j=f_(j+1)-f_j.                                      (1)
```

For `r=5,6`, consider the ordered three-block chart

```text
B_1=f_1-f_0=d_0,
B_2=f_2-f_1=d_1,
B_3^(r)=f_r-f_2=sum_(j=2)^(r-1)d_j,

U_r=B_1+xB_3^(r),              V_r=B_2+yB_3^(r).      (2)
```

There is exactly one pair `(x_r,y_r)` in the positive quadrant for
which the binary quadratic moment of `U_r+zV_r` divides its binary
cubic moment.  These two positive cubic-null planes satisfy

```text
support (0,1,2,5):             J_5(x_5,y_5)<0,
support (0,1,2,6):             J_6(x_6,y_6)>0,        (3)
```

where `J_r` is the quartic endpoint coefficient determinant of
THM-2872.  Both local response-row `(1,2,3)` determinants are strictly
positive throughout their positive parameter rays.  Consequently

```text
atomic adjacent-tensor positivity
 + separate local response TP3
 + simultaneous cubic nullity

does not determine a support-uniform sign of J.        (4)
```

The repaired general target is sign-blind: exclude a common zero of
the cubic eliminant and the endpoint remainder, support by support, or
find an additional support-sensitive orientation coordinate.

## 1. Cubic-null equations

Write

```text
g_0=L(U^2),              g_1=L(UV),             g_2=L(V^2),
t_0=L(U^3),              t_1=L(U^2V),
t_2=L(UV^2),             t_3=L(V^3).                   (5)
```

The quadratic moment

```text
q(z)=g_0+2g_1z+g_2z^2                              (6)
```

divides the cubic moment exactly when

```text
I_1=3t_1g_0g_2-t_3g_0^2-2t_0g_1g_2=0,
I_2=3t_2g_0g_2-2t_3g_1g_0-t_0g_2^2=0.              (7)
```

For each support in `(2)`, the exact subresultant sequence in `x` has
degree profile

```text
(4,3,2,1,0).                                          (8)
```

Up to a nonzero rational scalar, its last member factors as

```text
Res_x(I_1,I_2)=G_r(y)^2 P_r(y),                       (9)
```

where

```text
G_5=108y^2+12y+1,              disc(G_5)=-288,
G_6=437y^2+18y+1,              disc(G_6)=-1424,       (10)
```

and `P_r` is primitive of degree fifteen.  Thus `G_r` has no real zero.
Each `P_r` has exactly one coefficient sign change.  Exact endpoint
evaluation gives

```text
P_5(1/84)<0<P_5(1/83),
P_6(1/161)<0<P_6(1/159).                              (11)
```

Descartes' rule and `(11)` therefore give one and only one positive
root `y_r`, with

```text
1/84 <y_5<1/83,
1/161<y_6<1/159.                                      (12)
```

The degree-one subresultant, after removing its harmless `G_r`
content, is

```text
A_r(y)x-N_r(y).                                       (13)
```

Both `A_r,N_r` have degree ten.  Bernstein certificates on the
intervals `(12)` prove `A_r,N_r>0` and

```text
1/37 <x_5=N_5(y_5)/A_5(y_5)<1/35,
1/105<x_6=N_6(y_6)/A_6(y_6)<1/104.                   (14)
```

Exact evaluation in `QQ[y]/(P_r)` verifies both equations `(7)`.
The two polynomials in `(2)` are linearly independent, so positivity
of the factorial inner product makes their quadratic Gram determinant
strictly positive.  Equations `(12)--(14)` hence define genuine
positive nondegenerate cubic-null planes.

## 2. Opposite quartic endpoint signs

Let

```text
A_i=L(U^(4-i)V^i),                         0<=i<=4.   (15)
```

The chart-free endpoint coefficient determinant from THM-2872 is

```text
J=(2A_1g_0-A_0g_1)g_2^2
 -(2A_3g_2-A_4g_1)g_0^2.                            (16)
```

There is an important typing boundary here.  Since `B_1=f_1-f_0`,
one has `U_r(0)=-1`; the derivative/response-secant interpretation in
THM-2872 Section 4 is not available.  Only the coefficient algebra of
THM-2872 Section 2 is used: if `q` divided the quartic moment, then its
two linear remainder coefficients would vanish, and in particular

```text
J=0.                                                   (17)
```

Because `A_r` is invertible modulo `P_r`, put

```text
xbar_r=N_r A_r^(-1) in QQ[y]/(P_r),
K_r=J_r(xbar_r,y) mod P_r.                            (18)
```

Horner evaluation in the quotient ring gives `deg K_r=14`.  Direct
Bernstein certificates on the same rational intervals prove

```text
K_5(y)<0             for 1/84<=y<=1/83,
K_6(y)>0             for 1/161<=y<=1/159.             (19)
```

For an independent coprimality check, after primitive integral
normalization,

```text
Res(P_5,K_5)=5 mod 11,
Res(P_6,K_6)=4 mod 11.                                (20)
```

Thus no branch collision is hidden in `(19)`, and evaluating at the
unique roots of `P_r` proves `(3)`.  In particular, neither positive
cubic-null plane can also be a quartic multipole line.

## 3. Local TP3 survives on both branches

For any positive adjacent cone `W`, define

```text
Delta_1(W)=det [
 L(d_1W)   L(d_1W^2)   L(d_1W^3)
 L(d_2W)   L(d_2W^2)   L(d_2W^3)
 L(d_3W)   L(d_3W^2)   L(d_3W^3)
].                                                    (21)
```

Literal factorial expansion gives

```text
Delta_1(U_5)=12 p_(5,U)(x),
Delta_1(V_5)=12 p_(5,V)(y),
Delta_1(U_6)=12 p_(6,U)(x),
Delta_1(V_6)=12 p_(6,V)(y),                           (22)
```

where the coefficient vectors, in increasing degree order, are

```text
p_(5,U):
 (1,855,558525,261898305,24354608535,
  1685007161865,6073568246025),

p_(5,V):
 (2551,670995,87912270,6362772150,191150808345,
  2766944215260,6073568246025),

p_(6,U):
 (1,2056,4918729,11672252448,2264524229763,
  535438981677081,2991445280149032),

p_(6,V):
 (2551,1844344,845012299,249433808247,19258405287303,
  873935140816599,2991445280149032).                  (23)
```

Every coefficient in `(23)` is positive.  Hence all four local minors
are positive for positive parameters, including at the cubic-null
points.  THM-2853 simultaneously supplies positivity of the underlying
mixed adjacent-difference tensors.  This places the sign reversal
strictly beyond the ambient hostile in THM-2906, whose two witness
cells did not satisfy `(7)`.

It also explains why THM-2909 is not enough to orient the endpoint:
even exact response-row curvature is an intrapolynomial order, whereas
`J` compares the two ends of a mixed binary quotient.

## 4. Translation does not repair the orientation

The failure is not only a comparison between two different gap triples.
Apply the same construction to the translated supports

```text
(2,3,4,8),                 (3,4,5,9).                 (24)
```

Both have ordered block-length pattern `(1,1,4)`, and the second is the
unit translate of the first.  A separate exact quotient-ring audit finds
one positive cubic-null branch in each chart, with

```text
support (2,3,4,8):                 J>0,
support (3,4,5,9):                 J<0.               (25)
```

The exact isolating boxes include

```text
13476/1442287 <y_(2348)<10591/1133516,
7992889016/10^12 <x_(2348)<7992889020/10^12,

17523/1757477 <y_(3459)<7961/798452,
8033892855/10^12 <x_(3459)<8033892860/10^12.          (26)
```

The reduced endpoint polynomials are coprime to their cubic eliminants
modulo `7` and `13`, respectively.  Thus neither the ordered block
architecture nor its gap pattern or translation class orients `J`.
The base depth is retained information.

This identifies the first failed implication:

```text
positive ordered blocks + cubic nullity + local curvature
 + fixed gap pattern
      -/-> one coherent scalar endpoint orientation.  (27)
```

The loss occurs when the two-component quartic remainder
`(R_0,R_1)` is compressed to the scalar endpoint determinant `J`.
A transportable general certificate must instead retain the full
remainder vector, or equivalently a maximal-minor/Pluecker atlas, and
test its simultaneous vanishing sign-blindly.

## 5. What this changes

THM-2879 proves a negative pseudo-remainder for its consecutive
all-shift family.  Equation `(3)` shows that its sign is not a generic
consequence of positivity plus cubic nullity.  A universal syzygy with
a fixed sign on all ordered four-slot supports is impossible.

The surviving program is:

```text
fixed support A
  -> cubic branch eliminant P_A
  -> endpoint remainder K_A
  -> prove gcd(P_A,K_A)=1, without prescribing sign(K_A). (28)
```

This is precisely the sign-blind resultant strategy used by the
consecutive projective atlas.  The theorem proves only the two high
shared charts `(2)`.  It does not classify other real branches, middle
or low pivot charts, boundaries, arbitrary four-slot supports, SFC(4),
or any new case of GMC.

## 6. Exact companion

The exact companion:

1. reconstructs all quadratic, cubic and quartic tensors twice, by
   indexed factorial expansion and literal polynomial evaluation;
2. verifies the subresultant profile and factorization `(9)--(10)`;
3. proves the unique positive roots and rational boxes `(11)--(14)`;
4. evaluates `I_1,I_2,J` in the two quotient rings;
5. proves the signed degree-fourteen Bernstein certificates `(19)`;
6. checks the two mod-11 resultants `(20)`; and
7. reconstructs and locks all four TP3 polynomials in `(23)`.

Every truth-bearing gate uses an explicit `require`, so optimized mode
performs the same checks.  Run

```text
python3 04-computation/gmc_nonconsecutive_cubic_null_holonomy_reversal_thm2910.py
python3 -O 04-computation/gmc_nonconsecutive_cubic_null_holonomy_reversal_thm2910.py
```

Both modes byte-match

```text
05-knowledge/results/gmc_nonconsecutive_cubic_null_holonomy_reversal_thm2910.out.
```

An independent hostile audit reproduced the two moment systems, branch
polynomials, positive selectors, opposite endpoint signs, local TP3
polynomials and mod-11 coprimality values.  It also checked the
chart-free `a_0=0` use of THM-2872, replayed ordinary and optimized modes
against the stored transcript, and matched the declared LF hashes.
It separately certified the translated same-gap pair `(24)--(26)`,
which is not needed for the two headline computations.

**QED.**
