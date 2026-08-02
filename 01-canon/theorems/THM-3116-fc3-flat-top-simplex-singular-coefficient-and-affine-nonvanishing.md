---
id: THM-3116
title: "FC(3) flat-top simplex coefficient: exact Gamma/Mittag--Leffler shifts, a radial limit, and affine arithmetic nonvanishing"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.  Let
  L_3(x^a y^b z^c)=a!b!c!, let f have exact degree D, and suppose its top
  form is a(x+y+z)^D with a nonzero.  If A_1 is the restriction of f_(D-1)
  to the coordinate two-simplex, then
  L_3(f^m)/(a^m Gamma(Dm+3)) tends uniformly through the radial-projective
  formula to int_Delta exp(A_1/(Da)) dA.  Thus a flat-top FC(3)
  counterexample forces this simplex exponential integral to vanish.  After
  algebraic specialization it cannot vanish when A_1/(Da) is affine:
  exact divided-difference formulas plus Lindemann--Weierstrass handle every
  confluent case. The same calculation proves that an algebraic affine Q has
  int_Delta exp(Q)=1/2 iff Q=0, sharpening THM-3039's forced-level bridge.
  Hence a flat-top counterexample must have projective
  subleading degree at least two; in particular none has exact degree at
  most two.  This does NOT exclude nonflat leading forms and does not prove
  FC(3).  Exact controls verify the r^2 Jacobian, 125 Dirichlet monomials,
  Gamma(Dm+3), E_(D,3)/E_(D,4), seven direct layer expansions, affine
  formulas, and two hostile cancellation examples, identically under normal
  and -O execution.
source: codex-2026-08-02-fc3-simplex
depends_on: []
related:
  - THM-1510
  - THM-2017
  - THM-3018
  - THM-3039
external:
  - "Riesz--Markov--Kakutani representation theorem."
  - "octonion/mathematics/fc, revised research draft (not peer reviewed; only the radial and flat-symbol calculations are independently recovered here)."
script: 04-computation/fc3_flat_top_simplex_asymptotic_thm3116.py
output: 05-knowledge/results/fc3_flat_top_simplex_asymptotic_thm3116.out
---

# THM-3116 — the FC(3) flat-top simplex coefficient

## 1. Statement and scope

Write

```text
L_3(x^alpha y^beta z^gamma)=alpha! beta! gamma!,
S=x+y+z,
Delta={(u,v): u>=0, v>=0, u+v<=1},
w=1-u-v.
```

Let `f` have exact total degree `D>=1`, with homogeneous decomposition
`f=f_D+f_(D-1)+...+f_0`, and put

```text
A_j(u,v)=f_(D-j)(u,v,w),       0<=j<=D.                         (1)
```

Assume the leading form is flat:

```text
f_D=a S^D,                    a in C^*.                         (2)
```

Then

```text
 L_3(f^m)/(a^m Gamma(Dm+3))
      --> c_0(f):=int_Delta exp(A_1(u,v)/(D a)) du dv.          (3)
```

The convergence before the final simplex integration is uniform in
`(u,v) in Delta`.  Consequently, if every positive factorial moment of `f`
vanishes, then

```text
int_Delta exp(A_1/(D a)) du dv=0.                               (4)
```

After algebraic specialization preserving (2), (4) is impossible whenever
`A_1/(Da)` is affine on `Delta`.  Therefore every flat-top FC(3)
counterexample must have

```text
deg_Delta A_1 >= 2.                                             (5)
```

In particular there is no flat-top counterexample of exact degree `D=1` or
`D=2`.  This is a closure of the flat branch only.  No assertion is made
about a nonconstant projective leading form, so this theorem is not FC(3).

The closest proved mechanism is the one-variable factorial saddle of
THM-1510 and its uniform mixed-factorial use in THM-2017.  The new point is
that the simplex variable is retained and the full lower-layer multinomial
tail is controlled uniformly, giving the exact two-dimensional coefficient
rather than a pointwise heuristic.

## 2. The exact radial--projective and Mittag--Leffler formulas

The factorial functional has the positive integral representation

```text
L_3(h)=int_[0,infinity)^3 h(x,y,z) exp(-S) dx dy dz.             (6)
```

Use

```text
x=ru,  y=rv,  z=r(1-u-v).
```

The Jacobian is exactly `r^2`.  Hence, with

```text
P(u,v,r)=f(ru,rv,rw)=sum_(j=0)^D A_j(u,v) r^(D-j),              (7)
I_m=L_3(f^m)/Gamma(Dm+3),
```

one has

```text
I_m = 1/Gamma(Dm+3)
      int_Delta int_0^infinity exp(-r) r^2 P(u,v,r)^m dr du dv. (8)
```

For `E_(D,beta)(z)=sum_(m>=0) z^m/Gamma(Dm+beta)`, termwise
integration for sufficiently small `|s|` gives the exact shifts

```text
G(s):=sum_(m>=0) I_m s^m
 =int_Delta int_0^infinity exp(-r) r^2 E_(D,3)(sP) dr du dv,     (9)

J(s):=sum_(m>=0) I_m s^m/(Dm+3)
 =int_Delta int_0^infinity exp(-r) r^2 E_(D,4)(sP) dr du dv.    (10)
```

Thus the three-variable shifts are

```text
Gamma(Dm+3),  E_(D,3),  E_(D,4),                               (11)
```

not the two-variable `Gamma(Dm+2), E_(D,2), E_(D,3)` shifts.  This is the
precise dimension-three extension of the radial bookkeeping in the supplied
external FC(2) draft; no other part of that draft is used here.

As a normalization check, `I_0=1/Gamma(3)=1/2`, the coordinate area of
`Delta`, rather than `1`.

The bookkeeping is dimension-free. For `N` variables, the same substitution
has Jacobian `r^(N-1)` and gives

```text
Gamma(Dm+N),  E_(D,N),  E_(D,N+1),
L_N(f^m)/(a^m Gamma(Dm+N))
   --> int_(Delta_(N-1)) exp(A_1/(Da)) dV.                       (11a)
```

The proof below uses only that `N` is fixed. The specialization `N=3` is
where the affine arithmetic of section 4 becomes a triangle statement.

## 3. Exact finite expansion and proof of the flat limit

Divide by `a` and replace `f` by `a^(-1)f`; it is enough to prove (3) when
`a=1`.  Then `A_0=1`.  Expanding (8), let

```text
k_0+...+k_D=m,       J=sum_(j=1)^D j k_j.
```

The radial integral is still nonsingular because `0<=J<=Dm`, and gives the
exact finite identity

```text
I_m = int_Delta sum_(k_0+...+k_D=m)
        [m!/(k_0!...k_D!)] prod_(j=1)^D A_j(u,v)^(k_j)
        Gamma(Dm+3-J)/Gamma(Dm+3) du dv.                        (12)
```

Let `M` bound every `|A_j|` on the compact simplex.  Split (12) at
`J=Dm/2`.

For `J>Dm/2`, the absolute multinomial mass is at most `(1+DM)^m`, while

```text
Gamma(Dm+3-J)/Gamma(Dm+3)
 <= Gamma(floor(Dm/2)+3)/Gamma(Dm+3).                            (13)
```

Stirling makes (13), even after the exponential multinomial factor,
superexponentially small.

For `J<=Dm/2`, every factor in the descending Gamma quotient is at least
`Dm/2`, so

```text
Gamma(Dm+3-J)/Gamma(Dm+3) <= (2/(Dm))^J.                        (14)
```

Also `m!/(m-K)!<=m^K`, where `K=sum_(j>=1)k_j`.  Therefore the absolute
sum in this range is dominated by the product exponential with parameters

```text
M (2/D)^j m^(1-j),       1<=j<=D.                               (15)
```

Every parameter with `j>=2` tends to zero; the total contribution of tuples
using any lower layer `A_j`, `j>=2`, is `O(1/m)`.  For the remaining tuples
only `k_1=k` is nonzero, and for each fixed `k`, uniformly on `Delta`,

```text
binom(m,k) Gamma(Dm+3-k)/Gamma(Dm+3) A_1^k
       --> (A_1/D)^k/k!.                                        (16)
```

The `j=1` part of (15) supplies an exponential-series majorant independent
of `m`; its tail is uniformly summable. Let `H_m(u,v)` denote the finite sum
inside the simplex integral in (12). Dominated convergence in `k`, then
uniform convergence on `Delta`, proves

```text
H_m(u,v) --> exp(A_1(u,v)/D) uniformly,
I_m --> int_Delta exp(A_1/D) du dv                              (17)
```

after integration.  Restoring `a` gives (3).

This also identifies the first generating singular coefficient without any
sectorial analytic-continuation assumption.  In the normalization `a=1`,
the coefficients of `G` tend to `c_0(f)`, so the Abel theorem gives

```text
lim_(s up to 1) (1-s)G(s)=c_0(f),
lim_(lambda up to 1) D(1-lambda)G(lambda^D)=c_0(f).              (18)
```

Thus, when `c_0(f)!=0`, the first real-radial singular term is
`c_0(f)/(D(1-lambda))`.  The stronger sectorial `O(log(1-lambda))`
expansion claimed in the external draft is not needed and is not imported.
If every positive factorial moment vanishes then `I_m=0` for `m>=1`, so
(17) directly yields (4).

## 4. Arithmetic nonvanishing for an affine subleading layer

Suppose, after algebraic specialization, that

```text
q(u,v):=A_1(u,v)/(Da)=c+alpha u+beta v,
alpha,beta,c in Qbar.                                           (19)
```

The factor `exp(c)` is nonzero, so consider `c=0`.  If
`alpha`, `beta`, and `alpha-beta` are all nonzero, direct integration gives
the second divided difference

```text
int_Delta exp(alpha u+beta v) du dv
 = 1/(alpha beta)
   + exp(alpha)/(alpha(alpha-beta))
   + exp(beta)/(beta(beta-alpha)).                              (20)
```

The three coefficients are nonzero algebraic numbers and the exponents
`0,alpha,beta` are distinct algebraic numbers.  Lindemann--Weierstrass
therefore makes (20) nonzero.

Every confluent case is equally explicit:

```text
beta=0, alpha!=0:
  int_Delta exp(alpha u) = [exp(alpha)-1-alpha]/alpha^2;         (21)

alpha=beta!=0:
  int_Delta exp(alpha(u+v))
      = [(alpha-1)exp(alpha)+1]/alpha^2.                         (22)
```

Equation (21) cannot vanish because it would make `exp(alpha)` algebraic.
Equation (22) has value `1` at `alpha=1`; away from `1`, vanishing would
again make `exp(alpha)` algebraic.  The case `alpha=0` is symmetric, and a
constant `q` gives `exp(c)/2`.  This proves affine nonvanishing.

There is also an exact affine conclusion for THM-3039's homogeneous
forced-level bridge. Write an affine algebraic phase by its barycentric
vertex values:

```text
Q=lambda_0 w+lambda_1 u+lambda_2 v,      lambda_i in Qbar.
```

Then

```text
int_Delta exp(Q) du dv = exp[lambda_0,lambda_1,lambda_2],         (23)
```

the second divided difference of the exponential. If the three nodes are
distinct, (23) is

```text
sum_(i=0)^2 exp(lambda_i)/prod_(j!=i)(lambda_i-lambda_j).         (24)
```

Every displayed coefficient is nonzero algebraic. Adding the proposed value
`-1/2` as the coefficient of `exp(0)`, Lindemann--Weierstrass excludes
`(24)=1/2`; a node equal to zero only combines the constant coefficient and
leaves the other nonzero exponential coefficients.

If exactly two nodes equal `mu` and the third is `mu+delta`, `delta!=0`,
simplex symmetry and confluence give

```text
exp(mu) [exp(delta)-1-delta]/delta^2.                             (25)
```

Equality with `1/2` is a linear algebraic relation among
`exp(mu+delta), exp(mu), exp(0)`. When these exponents are distinct,
Lindemann--Weierstrass excludes it. If `mu=0`, it would make
`exp(delta)=1+delta+delta^2/2` algebraic. If `mu+delta=0`, it would make
`exp(-delta)` algebraic unless `delta=-1`, where the two sides are `1` and
`1/2`. Thus the confluent nonconstant case is also impossible. Finally, if
all nodes equal `lambda`, the integral is `exp(lambda)/2`, which equals
`1/2` for algebraic `lambda` only when `lambda=0`. Therefore

```text
Q affine over Qbar and int_Delta exp(Q)=1/2   iff   Q=0.          (26)
```

This does not prove the nonlinear exponential-period condition in THM-3039,
but it completely removes its affine locus.

Algebraic specialization is legitimate and load-bearing.  If a complex
counterexample satisfies (2), scale it so the leading form is exactly
`S^D`.  Parameterize only the lower coefficients (and, when relevant, the
linear subspace on which `A_1` is affine).  All equations `L_3(f^m)=0` have
rational polynomial coefficients.  Their ideal is proper at the assumed
complex point and is finitely generated by Noetherianity, so the weak
Nullstellensatz supplies a point over `Qbar` in the same flat/affine
subspace.  Applying (19)--(22) contradicts (4).

For `D=2`, `A_1` is automatically the restriction of a linear form, and for
`D=1` it is constant.  This proves the final assertion of section 1.  More
generally, for `D>=2`, an affine `A_1` is equivalent to
`S^(D-2)` dividing `f_(D-1)`; a flat counterexample must avoid that divisor.

## 5. Riesz--Markov: the exact measure statement and the cancellation gap

The positive functional

```text
ell(h)=int_Delta h du dv,       h in C(Delta),                   (27)
```

has norm and total mass `1/2`.  Riesz--Markov--Kakutani represents it by the
unique regular positive Borel measure `mu_Delta`, here coordinate area.  For
the polynomial `q=A_1/(Da)`, let `nu=q_*mu_Delta`.  Then

```text
c_0(f)=int_C exp(z) dnu(z),       nu(C)=1/2.                     (28)
```

This is an exact pushforward statement, but positivity of `nu` does not make
the complex integrand positive.  Equation `c_0=0` implies the necessary
geometric condition

```text
0 in conv(exp(q(Delta))).                                        (29)
```

It is impossible if `q` is real-valued, or more generally if the oscillation
of `Im q` is strictly less than `pi`: after one rotation, `exp(q(Delta))`
lies in an open half-plane.  Condition (29) is not sufficient for the fixed
pushforward measure.

Likewise, vanishing of the holomorphic moments `int z^m dnu` does not let
Riesz--Markov identify `nu` with a point mass.  The test algebra contains
`z` but not `conjugate(z)` and is not self-adjoint; uniform circle and disc
measures are the standard hostile controls.  A conjugate-moment, support, or
arithmetic sidecar is still needed.

## 6. Two sharp hostiles and the external-draft audit

The one-dimensional algebraic exponential-integral theorem asserted in the
supplied FC(2) draft cannot be extended to simplex periods by replacing an
interval with a triangle.  Already the nonconstant algebraic affine phase

```text
q(u,v)=u+v
```

has

```text
int_Delta exp(q) du dv = int_0^1 s exp(s) ds = 1.                (30)
```

Thus blanket **transcendence** of nonconstant algebraic simplex exponential
periods is false.  Formula (30) does not threaten (4), because its value is
nonzero.

Over unrestricted complex coefficients even nonvanishing fails.  Choose a
nonprincipal Lambert branch and put

```text
alpha=1+W_k(-exp(-1)),       k notin {0,-1}.
```

Then `alpha!=0` and

```text
(alpha-1)exp(alpha)+1=0,
int_Delta exp(alpha(u+v)) du dv=0.                               (31)
```

Hermite--Lindemann shows that this nonzero `alpha` is transcendental.  This
is the exact failure boundary: complex measure positivity and the flat
coefficient alone permit cancellation; algebraic specialization removes
the affine cancellation locus.

The external draft is a self-described research draft, not a peer-reviewed
or formally verified dependency.  Independently recovered here are only:

1. its radial Gamma/Mittag--Leffler bookkeeping (now with the exact
   dimension-three shifts (11)); and
2. its flat first coefficient (now proved by the finite expansion (12),
   without its nonflat constellation or arithmetic `E`-function chains).

No claim from its nonflat-top exclusion or exponential-integral proof is
used.  Equations (30)--(31) precisely identify why its final one-variable
step does not automatically become an FC(3) step.

## 7. Verification

Run

```bash
python3 04-computation/fc3_flat_top_simplex_asymptotic_thm3116.py
python3 -O 04-computation/fc3_flat_top_simplex_asymptotic_thm3116.py
```

The frozen controls verify:

* the `r^2` Jacobian and `125` exact Dirichlet/Gamma monomials;
* `54` coefficient checks for `Gamma(Dm+3)`, `E_(D,3)`, and `E_(D,4)`;
* the direct factorial functional against (12) for `m=0,...,6` on a
  nonradial quadratic with every lower layer present;
* all generic and confluent affine formulas;
* the exact affine forced-level law (26), algebraic-period hostile (30),
  Lambert cancellation (31), and numerical convergence to (3).

The two execution modes are byte-identical.  QED in the stated flat scope.

```text
source sha256 = 91282b41aeaf3712dc89de47c22724ca17719d98d938cc48828e263a48c530d9
output sha256 = 68adfea560b4fd728524ac22127001899902a1899f2c5edc66195bc4471bea0d
```
