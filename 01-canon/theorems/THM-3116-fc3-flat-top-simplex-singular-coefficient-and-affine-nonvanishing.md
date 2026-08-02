---
id: THM-3116
title: "FC(3) flat-top simplex coefficient: exact radial limit and algebraic one-coordinate quadratic exclusion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Let
  L_3(x^a y^b z^c)=a!b!c!, let f have exact degree D, and suppose its top
  form is a(x+y+z)^D with a nonzero.  If A_1 is the restriction of f_(D-1)
  to the coordinate two-simplex, then
  L_3(f^m)/(a^m Gamma(Dm+3)) tends uniformly through the radial-projective
  formula to int_Delta exp(A_1/(Da)) dA.  Thus a flat-top FC(3)
  counterexample forces this simplex exponential integral to vanish.  After
  algebraic specialization it cannot vanish when A_1/(Da) is affine:
  exact divided-difference formulas plus Lindemann--Weierstrass handle every
  confluent case.  More generally it cannot vanish whenever
  A_1/(Da)=P(1-lambda_i) for a barycentric coordinate lambda_i and an
  algebraic polynomial P of degree at most two.  The generic quadratic case
  uses an exact first-order E-function system, a complete functional-linear-
  independence audit, and Beukers Corollary 1.4 at the ordinary point z=1.
  The same families have int_Delta exp(Q)=1/2 iff Q=0, sharpening
  THM-3039's forced-level bridge.
  Hence a flat-top counterexample must have projective
  subleading degree at least two; in particular none has exact degree at
  most two.  This does NOT exclude nonflat leading forms and does not prove
  FC(3).  Exact controls verify the r^2 Jacobian, 125 Dirichlet monomials,
  Gamma(Dm+3), E_(D,3)/E_(D,4), seven direct layer expansions, affine and
  one-coordinate quadratic formulas, the quadratic E-function recurrence,
  and two hostile cancellation examples, identically under normal and -O
  execution.
source: codex-2026-08-02-fc3-simplex
depends_on: []
related:
  - THM-1510
  - THM-2017
  - THM-3018
  - THM-3039
external:
  - "Riesz--Markov--Kakutani representation theorem."
  - "F. Beukers, A refined version of the Siegel--Shidlovskii theorem, Annals of Mathematics 163 (2006), 369--379, Corollary 1.4: https://annals.math.princeton.edu/wp-content/uploads/annals-v163-n1-p08.pdf"
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

There is a further exact exclusion inside the first remaining degree.  Put
`lambda_0=w, lambda_1=u, lambda_2=v`.  For no `i` can the algebraically
specialized first projective layer have

```text
A_1/(Da)=P(1-lambda_i),       P in Qbar[s], deg P<=2.           (5a)
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

### One-coordinate quadratics: the elementary boundary family

There is a small nonlinear family where the simplex density itself supplies
the antiderivative missing in the general quadratic problem.  Let
`lambda_i` be any barycentric coordinate and set `s=1-lambda_i`.  A slice of
the triangle at fixed `s` has coordinate length `s`, so for every continuous
`Phi`,

```text
int_Delta Phi(1-lambda_i) du dv = int_0^1 s Phi(s) ds.          (26a)
```

Consequently, for algebraic `b,c`,

```text
int_Delta exp(c+b(1-lambda_i)^2) du dv
 = exp(c)(exp(b)-1)/(2b),                         b!=0,
 = exp(c)/2,                                      b=0.          (26b)
```

The value is nonzero: for nonzero algebraic `b`, Hermite--Lindemann excludes
`exp(b)=1`.  It equals the forced level `1/2` only when `b=c=0`.  Indeed, for
`b!=0`, equality would give the algebraic linear relation

```text
exp(b+c)-exp(c)-b exp(0)=0.                                   (26c)
```

Lindemann--Weierstrass excludes it when the three exponents are distinct.
The only collisions are `c=0` and `b+c=0`; each would make respectively
`exp(b)` or `exp(-b)` algebraic.  The case `b=0` reduces to
`exp(c)=1`, hence `c=0`.

Algebraicity is sharp here.  Taking `b=-2 pi i` and `c=2 pi i` makes
`exp(b)=1` and the period in (26b) zero.  With `s=1-u`, this is the
polynomial phase `2 pi i(2u-u^2)` underlying THM-3039's canonical
nonalgebraic moment hostile.

The mechanism is dimension-free and explains why the exponent is a square
for FC(3).  On `Delta_(N-1)`, the pushforward of coordinate volume by
`s=1-lambda_i` has density `s^(N-2)/(N-2)!`.  Hence

```text
int_(Delta_(N-1)) exp(c+b(1-lambda_i)^(N-1)) dV
 = exp(c)(exp(b)-1)/(b (N-1)!),                 b!=0,           (26d)
```

with the continuous value `exp(c)/(N-1)!` at `b=0`.  The same argument gives
nonvanishing and equality with the simplex volume only at `b=c=0` over
`Qbar`.  For `N=3`, (26d) is exactly (26b).

This square family is not an arbitrary lucky example; it is the maximal
genuinely quadratic case detected by a polynomial antiderivative.  For a
one-coordinate phase `q=P(s)`, such a certificate has the form

```text
d/ds [R(s) exp(P(s))]=s exp(P(s)),
R'(s)+P'(s)R(s)=s.                                             (26e)
```

A degree comparison classifies polynomial solutions `R`.  Constant `P`
admits `R=s^2/2`.  If `P=Bs+C`, `B!=0`, then
`R=s/B-1/B^2`.  If `P=As^2+Bs+C`, `A!=0`, the highest-degree term forces
`R` to be constant, after which (26e) forces `B=0` and `R=1/(2A)`.
For `deg P>=3`, the highest term of `P'R` cannot be cancelled by `R'`, so no
nonzero polynomial solution exists.  Thus (26b) is exactly the last
one-coordinate polynomial phase reducible to endpoints by this elementary
mechanism.  The generic quadratic `As^2+Bs+C` with `AB!=0` is the first
non-elementary case; the following argument clears it by retaining its full
parameter `E`-function rather than asking for an elementary antiderivative.

### The generic quadratic: exact system and `E`-function arithmetic

Let `A,B,C in Qbar` with `AB!=0`, and define

```text
q(s)=A s^2+B s+C,
H(z)=int_0^1 exp(z q(s)) ds,
I(z)=int_0^1 s exp(z q(s)) ds,
lambda=A+B+C.                                                  (26f)
```

Integrating `d(exp(zq(s)))/ds` gives the exact endpoint identity

```text
2Az I(z)+Bz H(z)=exp(lambda z)-exp(Cz).                         (26g)
```

A second integration by parts, now using
`d(s exp(zq(s)))/ds`, eliminates the second `s`-moment and gives

```text
4Az H'(z)+[2A+(B^2-4AC)z]H(z)
       =(2A+B)exp(lambda z)-B exp(Cz).                          (26h)
```

The function `H` is an `E`-function, not merely holonomic.  Indeed,

```text
H(z)=sum_(m>=0) h_m z^m/m!,       h_m=int_0^1 q(s)^m ds.        (26i)
```

All `h_m` lie in one number field.  After clearing fixed denominators of
`A,B,C`, expansion of `q^m` shows that a common denominator for
`h_0,...,h_m` divides a fixed exponential factor times
`lcm(1,...,2m+1)`, whose logarithm is `O(m)`.  Every conjugate of `h_m` is
bounded exponentially in `m`.  Equation (26h) supplies the required linear
differential system; eliminating its two endpoint exponentials, which each
satisfy a first-order equation, gives a homogeneous scalar differential
equation over `Qbar[z]` for `H`.  These are exactly the arithmetic, height,
and differential conditions in the definition of an `E`-function.

It remains to audit functional linear independence before specializing.
Take the distinct functions among

```text
H(z), exp(lambda z), exp(Cz), 1.                               (26j)
```

Distinct algebraic exponentials are linearly independent over `Qbar(z)`.
Suppose a relation involving `H` existed and solve it for `H` as a rational
linear combination of the distinct exponentials.  Substitute it into
(26h), and inspect the coefficient of the exponential whose exponent is
`C`.  Its rational coefficient `R` would have to satisfy

```text
zR'(z)+(1/2+kappa z)R(z)=d,
kappa=B^2/(4A),
d=-B/(4A)                 if lambda!=C,
d=1/2                     if lambda=C.                         (26k)
```

Here `kappa!=0` and `d!=0`.  The second line includes the combined source
when `lambda=C`, equivalently `A+B=0`: `(2A+B)-B=2A` before division by
`4A`.  Equation (26k) has no rational solution.  A pole of order `n` at a
nonzero point produces an uncancellable pole of order `n+1` in `zR'`.  At
zero, its leading coefficient is `(-n+1/2)r_(-n)`, also nonzero.  Thus `R`
has no finite poles and is a polynomial.  But `kappa zR` raises the degree,
so no nonzero polynomial can equal the nonzero constant `d`; `R=0` also
fails.

This single coefficient audit covers every collision.  If `C=0`, then
`exp(Cz)=1` and its source is unchanged.  If `lambda=0`, only that endpoint
exponential merges with `1`.  If `lambda=C`, the source becomes `2A`; and
when `lambda=C=0`, both endpoint exponentials and `1` merge but that same
nonzero combined source remains.  Consequently the distinct functions in
(26j) are linearly independent over `Qbar(z)`.

They form a first-order rational system obtained from (26h) and the three
exponential equations; after duplicates are removed, its common denominator
is `T(z)=z`.  Beukers, *A refined version of the Siegel--Shidlovskii
theorem*, Corollary 1.4 (Annals 163 (2006), p. 371), therefore applies at
the algebraic point `xi=1`, because `xi T(xi)!=0`.  It makes the values of
the distinct functions in (26j) linearly independent over `Qbar`.

Now put `z=1` in (26g).  If `I(1)=0`, it becomes

```text
B H(1)-exp(lambda)+exp(C)=0;                                  (26l)
```

if `I(1)=1/2`, it becomes

```text
B H(1)-exp(lambda)+exp(C)+A=0.                                (26m)
```

Both are nontrivial even after every endpoint collision because the
coefficient `B` of `H(1)` is nonzero.  Corollary 1.4 excludes both.  The
axes were already settled: `A=0` is affine, while `B=0` is (26b).  Therefore,
for every barycentric coordinate `lambda_i` and every algebraic polynomial
`P` of degree at most two,

```text
int_Delta exp(P(1-lambda_i)) du dv !=0,
int_Delta exp(P(1-lambda_i)) du dv =1/2   iff   P=0.            (26n)
```

Algebraic specialization is legitimate and load-bearing.  If a complex
counterexample satisfies (2), scale it so the leading form is exactly
`S^D`.  Parameterize only the lower coefficients (and, when relevant, the
linear subspace on which `A_1` is affine or lies in one fixed span
`<1,1-lambda_i,(1-lambda_i)^2>`).  All equations `L_3(f^m)=0` have
rational polynomial coefficients.  Their ideal is proper at the assumed
complex point and is finitely generated by Noetherianity, so the weak
Nullstellensatz supplies a point over `Qbar` in the same constrained
subspace.  Applying (19)--(22) or (26n) contradicts (4).

For `D=2`, `A_1` is automatically the restriction of a linear form, and for
`D=1` it is constant.  This proves the final assertion of section 1.  More
generally, for `D>=2`, an affine `A_1` is equivalent to
`S^(D-2)` dividing `f_(D-1)`; a flat counterexample must avoid that divisor.
If its first projective layer is quadratic, it must also avoid each of the
three one-coordinate subspaces (5a).

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

This hostile also identifies the exact missing coordinate in a naive
`E`-function repair.  Introduce the parameter integral

```text
F(z)=int_Delta exp(z(u+v)) du dv
    =((z-1)exp(z)+1)/z^2
    =sum_(m>=0) z^m/(m!(m+2)).                                  (30a)
```

It is an `E`-function: its algebraic coefficients `1/(m+2)` have
exponentially bounded heights and common denominators, and it satisfies

```text
z F''(z)+(3-z)F'(z)-2F(z)=0.                                   (30b)
```

Nevertheless `F(1)=1`, because there is already the exact functional
relation

```text
z^2 F(z)-(z-1)exp(z)-1=0.                                      (30c)
```

The value relation `F(1)-1=0` is precisely the specialization of (30c): the
coefficient of `exp(z)` vanishes at `z=1`.  This is fully consistent with
Beukers' lifting theorem.  Therefore holonomicity or `E`-function status of
a general simplex parameter integral cannot by itself supply the needed
arithmetic obstruction.  One must control its complete functional-relation
module (or prove the relevant independence) and audit coefficient
degeneration at the evaluation point.  Affine phases are especially hostile:
their divided-difference formulas make such exponential-polynomial relations
systematic, not exceptional.

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
* all generic and confluent affine formulas and the derivative-aligned
  quadratic boundary formula (26b);
* eight exact coefficient recurrences for the quadratic endpoint identity
  (26g), nine for its `E`-function ODE (26h), and the distinct/colliding
  endpoint sources used in (26k);
* the exact affine forced-level law (26), algebraic-period hostile (30), its
  exact `E`-function relation and ODE (30a)--(30c), Lambert cancellation
  (31), and numerical convergence to (3).

The two execution modes are byte-identical.  QED in the stated flat scope.

An independent hostile audit on 2026-08-02 checked the flat-limit
multinomial/Gamma domination, every affine confluent
Lindemann--Weierstrass case, the derivative-aligned quadratic slice and
forced-level collision cases, the generic quadratic `E`-function coefficient
bounds and ODE, every endpoint-exponent collision pattern, the rational
pole/degree obstruction, the ordinary-point hypotheses of Beukers Corollary
1.4, the affine `E`-function ODE/relation, and the polynomial-antiderivative
classification.  It found no mathematical defect.
The audit rechecked the displayed exact identities and proof mechanisms but
did **not** supply a second independent implementation of the verification
script; `VERIFIED-EXACT` continues to refer to the frozen implementation
listed below.

```text
source sha256 = 02687e5c59e320d2a6be6d2f82610fa7130ed804e4350552c7658405b887bd8f
output sha256 = ce3ed343b6c6817c1d62f7a61f971bf8d792836706d60fbf4121335b320261fb
```
