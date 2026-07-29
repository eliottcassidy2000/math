---
id: THM-2879
title: "All-shift cubic-null endpoint-holonomy exit"
status: >
  PROVED + VERIFIED-EXACT.  On the unique positive cubic-null plane
  U_n=d_n+x_n d_(n+2), V_n=d_(n+1)+y_n d_(n+2) of THM-2855, the quartic
  endpoint-holonomy determinant is strictly negative for every integer
  n>=1.  Exact substitution through the positive linear subresultant and
  pseudo-division by the degree-fifteen branch eliminant leaves a
  degree-fourteen polynomial whose 14,745 coefficients are all negative.
  Hence beta*kappa_U<alpha*kappa_V, so no member of this infinite shifted
  family is a shared cubic--quartic multipole line.  A projective rechart
  also makes the shared-middle positive wedge cubic-empty; among strict
  cone-cutting consecutive wedges only the shared-low, equivalently
  high-chart X>0,Y<0, sector remains.  Cone boundaries and the general
  mixed four-slot branch remain open.
source: codex/mixed-secant-2026-07-29
depends_on:
  - THM-2855-shifted-positive-cone-transverse-asymptotic-family
  - THM-2872-four-slot-shared-multipole-quartic-norm-and-response-secant-reduction
related:
  - THM-2846-arbitrary-positive-cone-moment-three-transverse-boundary
  - THM-2873-two-ray-factorial-response-tp3-curvature
  - THM-2875-all-order-initial-pascal-ladder-and-support-square-range
script: 04-computation/gmc_all_shift_cubic_null_endpoint_holonomy_thm2879.py
output: 05-knowledge/results/gmc_all_shift_cubic_null_endpoint_holonomy_thm2879.out
script_sha256: 7aae88a219689a910229e69fd03cc7a28ae53ef671d8297a274f237ee8a7476e
output_sha256: b0cba47d600dcbb290fd960ca2967a9ca592f3e51703ad1df3aea75fce39619e
hash_basis: LF-normalized bytes
---

# THM-2879 -- all-shift cubic-null endpoint-holonomy exit

**PROVED + VERIFIED-EXACT.**

Put

```text
L(s^q)=q!,                  f_j=s^j/j!,
d_j=f_(j+1)-f_j.                                      (1)
```

THM-2855 proves that for every integer `n>=1` there is exactly one positive
pair

```text
x_n>0,                    y_n>0                       (2)
```

such that, for

```text
U_n=d_n+x_n d_(n+2),
V_n=d_(n+1)+y_n d_(n+2),                              (3)
```

the binary quadratic moment divides the binary cubic moment:

```text
L((U_n+zV_n)^2) divides L((U_n+zV_n)^3).              (4)
```

Write

```text
g_0=L(U_n^2),        g_1=L(U_nV_n),       g_2=L(V_n^2),
A_i=L(U_n^(4-i)V_n^i),                    0<=i<=4.    (5)
```

The endpoint-holonomy determinant of THM-2872 is

```text
J_n=
 (2A_1g_0-A_0g_1)g_2^2
 -(2A_3g_2-A_4g_1)g_0^2.                             (6)
```

Then

```text
J_n<0                                                   (7)
```

for every integer `n>=1`.  In the response coordinates of THM-2872,

```text
J_n=(2/3)g_0^2g_2^2
        (beta kappa_U-alpha kappa_V),                 (8)
```

so `(7)` is exactly

```text
beta kappa_U<alpha kappa_V.                           (9)
```

Thus the two endpoint determinations of a possible quartic quotient never
agree.  No plane in the entire shifted family `(3)` is a shared
cubic--quartic multipole line.

## 1. Positive common tensor normalization

For `k in {2,3,4}`, put

```text
M_k(n)=[k(n+2)]!/(n+2)!^k>0.                         (10)
```

For offsets `delta_i in {0,1,2}`, direct expansion of the adjacent
differences gives the rational function

```text
tau_(delta_1,...,delta_k)(n)
 =1/M_k(n) L(prod_i d_(n+delta_i))

 =sum_(epsilon_i in {0,1})
   (-1)^(k-sum epsilon_i)
   [kn+sum delta_i+sum epsilon_i]!
   (n+2)!^k
  /(
    [kn+2k]!
    prod_i[n+delta_i+epsilon_i]!
   ).                                                  (11)
```

Every denominator in `(11)` is positive for `n>=1`.  Replacing every
order-`k` tensor by `tau` divides all tensors of that order by the same
positive number.  Hence it preserves the cubic divisibility equations.
If `J_n^#` is the determinant `(6)` in the normalized tensors, then

```text
J_n=M_4(n)M_2(n)^3 J_n^#.                             (12)
```

Thus `J_n` and `J_n^#` have the same sign.  The exact companion checks
`351` literal `2^k` tensor expansions and three independent instances of
the scale law `(12)`.

## 2. The exact positive cubic branch

Let `P_1(n,x,y),P_2(n,x,y)` be the cleared cubic remainders.  Their
subresultant degrees in `x` are

```text
4,3,2,1,0.                                            (13)
```

The resultant factors as

```text
Res_x(P_1,P_2)
 =C(n) G_n(y)^2 P_n(y),                               (14)
```

up to a nonzero rational integer, where `C(n)` has no zero for `n>=1`,

```text
G_n(y)
 =(n+2)+2(2n+3)y+2(2n+3)y^2,                         (15)
```

and `P_n` has degree fifteen.  The discriminant of `(15)` is

```text
-4(2n+3)<0,                                           (16)
```

so `G_n` contributes no real branch.

Orient `P_n` so that its leading coefficient

```text
lambda(n)=LC_y(P_n)                                   (17)
```

is positive.  In fact every coefficient of `lambda` as a polynomial in
`n` is positive.  The coefficient-threshold staircase proved in THM-2855
and replayed by the present companion gives exactly one sign change in
the descending coefficient list of `P_n` for every integer `n>=1`.
Its leading coefficient is positive and its constant coefficient is
negative.  Hence `P_n` has exactly one positive root `y_n`.

After its positive common content is removed, the linear subresultant is

```text
mathcal A(n,y)x-mathcal N(n,y).                       (18)
```

The two coefficient dictionaries have respectively

```text
251 and 253 terms,                                    (19)
```

and every coefficient in both is positive.  Therefore at the positive
root,

```text
x_n=mathcal N(n,y_n)/mathcal A(n,y_n)>0.              (20)
```

This is the exact THM-2855 selector, now retained as the coordinate through
which the mixed endpoint statistic is transported.

## 3. The universal negative pseudo-remainder

In normalized tensors, write

```text
J_n^#=E(n,x,y)/D(n),                                  (21)
```

where

```text
D(n)=
 1024(n+3)(2n+1)^2(2n+3)^4
 (4n+1)(4n+3)(4n+5)(4n+7)>0,                         (22)
```

and `E` has degree five in `x`.  Clear the positive selector denominator:

```text
F(n,y)
 =mathcal A(n,y)^5
  E(n,mathcal N(n,y)/mathcal A(n,y),y).               (23)
```

This is an ordinary polynomial with

```text
deg_n F=121,        deg_y F=55,        terms(F)=6820. (24)
```

Pseudo-divide `F` by `P_n` in the variable `y`.  Since

```text
55-15+1=41,
```

the standard pseudo-division identity is

```text
lambda(n)^41 F(n,y)
 =Q(n,y)P_n(y)+K(n,y),             deg_y K<=14.       (25)
```

Exact expansion gives

```text
deg_n K=982,         deg_y K=14,
terms(K)=14745.                                         (26)
```

Every one of the `14,745` coefficients of `K` as a polynomial in `(n,y)`
is strictly negative.  The canonical coefficient-dictionary digest is

```text
c1ed1aa0ff682a7226f67c752aceb7bb
4924708a2973126fe903c62c686d67a2.                     (27)
```

Consequently,

```text
K(n,y)<0                         for n>=1, y>0.        (28)
```

At `(n,y_n)`, equation `(25)` and `P_n(y_n)=0` give

```text
lambda(n)^41F(n,y_n)=K(n,y_n)<0.                      (29)
```

All factors cleared in `(17)`, `(22)`, and `(23)` are positive.  Therefore

```text
J_n^#<0,
```

and `(12)` proves `(7)`.

The sign certificate is a mixed endpoint coordinate.  Separate TP3
curvature of the two response curves in THM-2873, and the initial Pascal
ladders of THM-2875, do not compare those two curves.  The selector plus
the reduced pseudo-remainder does.

## 4. Shared-line exclusion

For `n>=1`, both polynomials `(3)` vanish at zero.  Their positive span is
exactly the nonnegative wedge of their plane: the `d_n` coefficient is the
coefficient of `U_n`, and the `d_(n+1)` coefficient is the coefficient of
`V_n`.  Thus `U_n,V_n` are the positive extreme rays required in
THM-2872's response-secant chart.

Equation `(4)` supplies its cubic hypothesis.  THM-2872 says that quartic
divisibility first requires endpoint equality

```text
beta kappa_U=alpha kappa_V.                           (30)
```

Equation `(9)` excludes `(30)` strictly.  The midpoint defect is never
reached.  Therefore

```text
L((U_n+zV_n)^2) does not divide
L((U_n+zV_n)^4)                                      (31)
```

for every positive cubic-null pair `(x_n,y_n)`, and either complex root of
the quadratic has nonzero fourth factorial moment.

At `n=1`, this recovers the canonical THM-2846/2873 endpoint exit without
an interval enclosure.  More importantly, it fills THM-2855's finite-depth
gap between the depth-one control and its eventual asymptotic exit.

## 5. Three strict wedge charts

There is a useful projective consequence beyond the family `(3)`.  Put

```text
e_1=d_n,             e_2=d_(n+1),          e_3=d_(n+2).
```

A two-plane whose normal has strict mixed signs cuts the positive
`(e_1,e_2,e_3)` cone in a two-ray wedge.  After scaling the rays, exactly
one of the three coordinates is shared.  The three charts are

```text
high:    c_3=x c_1+y c_2,
middle:  c_2=x c_1+y c_3,
low:     c_1=x c_2+y c_3,                 x,y>0.      (32)
```

The high chart is exactly `(3)` and is excluded from the shared quartic
line by `(7)`.

Suppose a middle-chart plane were cubic-divisible.  Since `y>0`, its
equation can be solved as

```text
c_3=(-x/y)c_1+(1/y)c_2.                               (33)
```

This is the high chart with

```text
X=-x/y<0,                    Y=1/y>0.                 (34)
```

This is an invertible reparametrization of the same binary plane, so it
preserves divisibility of its quadratic moment into its cubic moment.
Section 2 classifies every cubic-divisible high-chart plane with `Y>0`:
the resultant has no positive-`Y` root outside `P_n`, and the linear
subresultant forces the exact selector

```text
X=mathcal N(n,Y)/mathcal A(n,Y)>0.                    (35)
```

Equations `(34)--(35)` contradict one another.  Thus the strict
shared-middle wedge is cubic-empty, before the quartic is tested.

Finally, the low chart becomes

```text
c_3=(1/y)c_1-(x/y)c_2,                                (36)
```

so in high coordinates it is exactly

```text
X=1/y>0,                     Y=-x/y<0.                (37)
```

Therefore the sole strict cone-cutting orientation not discharged by the
high exit or the middle rechart is the high-chart sector

```text
X>0,                         Y<0.                     (38)
```

This is a reduction, not an emptiness theorem for `(38)`.  The assumptions
`x,y>0` are load-bearing: zero coefficients are cone-boundary charts and
are not covered.  Cone-avoiding planes are also outside `(32)`.

## 6. Sharp boundary and hostile controls

The integer lower bound is exact for the branch statement.  At `n=0`, all
sixteen coefficients of `P_0(y)` are positive, so

```text
#{y>0:P_0(y)=0}=0.                                    (39)
```

There is no bottom positive cubic-null branch to which `(7)` could be
extended.

Independent fixed-depth calculations at

```text
n=1,2,8                                               (40)
```

start from the literal factorial tensors rather than `(11)`.  In all
three cases they recover the degree profile `(13)`, one positive
degree-fifteen branch, positive degree-ten selector polynomials, and a
fifteen-term quartic pseudo-remainder with every coefficient negative.
These controls detect normalization, selector-sign, and pseudo-division
orientation errors.

The theorem does **not** prove that every cone-cutting four-slot cubic line
has negative endpoint holonomy.  It closes the infinite shared-high-tooth
family `(3)`, eliminates the strict middle wedge by recharting, and leaves
the low sector `(38)`.  Cone-avoiding planes, zero-coordinate boundaries,
the `a_0=0` derivative boundary, the midpoint wall after a possible endpoint
equality, and the general mixed shared line remain open.  No claim about
SFC(4), the full Gaussian Moment Conjecture, or the planar Jacobian
conjecture follows.

## 7. Exact verification

Run

```text
python3 04-computation/gmc_all_shift_cubic_null_endpoint_holonomy_thm2879.py
python3 -O 04-computation/gmc_all_shift_cubic_null_endpoint_holonomy_thm2879.py
```

Both modes byte-match

```text
05-knowledge/results/gmc_all_shift_cubic_null_endpoint_holonomy_thm2879.out.
```

The companion verifies:

1. all normalized/direct tensors and the positive determinant scale;
2. the complete resultant profile, positive quadratic factor, and
   degree-fifteen sign staircase;
3. the sharp `n=0` no-branch boundary;
4. the positive linear selector and all clearing denominators;
5. the pseudo-division exponent and every coefficient of `(26)`; and
6. three independent literal-depth endpoint controls.

**QED.**
