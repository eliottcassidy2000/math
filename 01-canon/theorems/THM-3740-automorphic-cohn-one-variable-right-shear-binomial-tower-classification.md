---
id: THM-3740
title: "Automorphic Cohn one-variable right-shear all-exposure classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For an arbitrary
  diagonal constant right factor followed by one arbitrary polynomial upper
  right shear, every compatible lower-row closure is classified; dually,
  every compatible upper-row closure after one lower right shear is
  classified.  A lower closure forces v_X=0 and an upper closure forces
  u_T=0.  The remaining one-variable parameter integrates uniquely to a
  nonlinear source translation, and closure occurs exactly on the paired
  integer-depth binomial slopes of THM-3734.  Exceptional negative transport
  resonances are empty.  Every survivor is source-automorphic to a THM-3734
  potential and has no polynomial mate.  On the resulting one-variable right
  boundary, the crossed exposed rows are also completely classified: only
  the zero-parameter diagonal tower, two constant Broughton sheets, and one
  nonconstant depth-one sheet survive in each orientation.
source: root + jc_sparse_direct_search / 2026-08-22
audit: >
  PASS.  An independent audit rederived both compatible and both crossed raw
  PDEs, every mixed-variable top-degree case, the p=m=1 edges, both negative
  transport-resonance exclusions, both connection conjugacies and gradients,
  every zero/constant/depth-one boundary, and the arbitrary-mate transfers.
  Normal, optimized, and frozen output agree; script/output/semantic hashes
  and CHECKS=1171 match.
depends_on:
  - THM-3721-automorphic-cohn-one-right-shear-nonentry
  - THM-3734-automorphic-cohn-diagonal-binomial-divided-power-towers
  - THM-3736-automorphic-cohn-complete-constant-sl2-polynomial-exposure-classification
related: []
script: 04-computation/jc2_automorphic_cohn_one_variable_right_towers_thm3740.py
output: 05-knowledge/results/jc2_automorphic_cohn_one_variable_right_towers_thm3740.out
script_sha256: c800df5287182640b44dddebb07335408cd1166c09bc8e854e08c44f3e734df3
output_sha256: e16347cfd6deb288df1847ef3a42644afe334b8756106ee2e23c69cfe4dd6764
semantic_sha256: e3fc3fd4f99bda017d984a968d2d543d8292f779a0745a99681285cdb7b01cff
hash_basis: raw LF bytes
---

# THM-3740 -- one-variable right shears transport the whole binomial tower

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  THM-3721 found a
one-variable survivor behind the first right-shear Cohn cell, but only its
depth-one Broughton form was then visible.  THM-3734 reveals that it is the
first member of a complete divided-power inheritance tower.  The right shear
is not an extra deformation of the potential: it is an integrable connection
whose primitive translates one source coordinate.

Work over a characteristic-zero field containing the displayed square roots.
Put

```text
M_0=[4T^2,2XT-1;1+2XT,X^2],
D_a=diag(a,a^{-1}),                  A=a^2!=0,
E_+(v)=[1,v;0,1],                    E_-(u)=[1,0;u,1].  (1)
```

Use `curl(P,Q)=partial_T P-partial_X Q`.  In each matrix below write its rows
as `alpha,beta`.

## 1. Complete compatible-exposure classification

### Upper right shear, lower exposed row

Let

```text
N_+=M_0D_aE_+(v),                    v in k[X,T].       (2)
```

There is a polynomial `h` for which `beta+h alpha` is closed if and only if
there is an integer `r>=1` such that

```text
A=(r+1)/(2r),                        v=v(T).            (3)
```

Write `v(T)=sum_(j>=0) v_jT^j` and define the unique polynomial

```text
B(T)=sum_(j>=0) v_j T^(j+1)/[j+1+r/(r+1)],            (4)
Y=X+B(T),
h=B(T)/(2T)+Y^2 H_r^L(YT),                            (5)
H_r^L(z)=[(1+2z/r)^r-1-2z]/(4z^2).
```

All quotients in `(4)--(5)` are polynomial.  The closed row is the gradient
of

```text
Q=aY Phi_r^L(YT),
Phi_r^L(z)=r[(1+2z/r)^(r+1)-1]/[2(r+1)z].             (6)
```

### Lower right shear, upper exposed row

Let

```text
N_-=M_0D_aE_-(u),                    u in k[X,T].       (7)
```

There is a polynomial `h` for which `alpha+h beta` is closed if and only if
there is an integer `r>=1` such that

```text
A=r/[2(r+1)],                        u=u(X).            (8)
```

For `u(X)=sum_(j>=0)u_jX^j`, put

```text
C(X)=sum_(j>=0) u_j X^(j+1)/[j+1+r/(r+1)],            (9)
V=T+C(X),
h=2C(X)/X+V^2 H_r^U(XV),                              (10)
H_r^U(z)=[-(1-2z/r)^r+1-2z]/z^2.
```

Then the closed row is the gradient of

```text
Q=a^{-1}V Phi_r^U(XV),
Phi_r^U(z)=-r[1-(1-2z/r)^(r+1)]/[2(r+1)z].            (11)
```

Neither potential `(6)` nor `(11)` has a polynomial Jacobian mate.  This
classifies the compatible lower exposure after `E_+` and upper exposure after
`E_-` for arbitrary bivariate right parameters.  The next section also
classifies both crossed exposures once the right parameter is one-variable.
Crossed exposures with genuinely bivariate right data and words with two
nonconstant right factors remain outside the theorem.

## 2. Both crossed one-variable exposures

The compatible gate naturally leaves `v=v(T)` and `u=u(X)`.  On this whole
one-variable boundary, the crossed rows have the following complete list.

For `M_0D_aE_+(v(T))`, the upper row `alpha+h beta` is closed exactly as
follows:

```text
v=0:                 the full diagonal upper tower of THM-3734;
deg v>=1:            A=1/4, h=0;
v=c in k*:           (A,h)=(1/4,0) or (1,3/c).        (C1)
```

For `M_0D_aE_-(u(X))`, the lower row `beta+h alpha` is closed exactly as
follows:

```text
u=0:                 the full diagonal lower tower of THM-3734;
deg u>=1:            A=1, h=0;
u=c in k*:           (A,h)=(1,0) or (1/4,3/(4c)).     (C2)
```

To prove the nonconstant cases, their exact crossed equations are

```text
U_A(h)-Av[(1+2XT)h_X+2Th]=0,                          (C3)
L_A(h)+u[(2XT-1)h_T+2Xh]=0.                           (C4)
```

If `p=deg_Tv>=1` and nonzero `h` has highest `T` coefficient `c_N(X)`, the
unique degree-`p+N+1` coefficient of `(C3)` is

```text
-2A lead(v)[Xc_N'+c_N].                               (C5)
```

It would force `c_N=C/X`, impossible for a nonzero polynomial.  Hence `h=0`,
and `(C3)` reduces to `2(4A-1)T=0`, giving `A=1/4`.  The dual top coefficient
in `(C4)` is

```text
2 lead(u)[Tc_N'+c_N],                                 (C6)
```

so `h=0` and `2(A-1)X=0`, giving `A=1`.  When the right parameter is zero,
THM-3734 applies; when it is a nonzero constant, THM-3736 rules out
nonconstant `h` and the two raw THM-3726 equations give exactly the two
displayed sheets.

The nonconstant `h=0` potentials are still source images of depth-one
Broughton components.  At `A=1/4`,

```text
Q=a^{-1}(XT^2-T)+4a integral(T^2v(T) dT)
 =a^{-1}[(X+S(T))T^2-T],
S(T)=T^{-2} integral(T^2v(T) dT) in Tk[T].             (C7)
```

At `A=1`, the dual potential is

```text
Q=a(X+X^2T)+a^{-1} integral(X^2u(X) dX)
 =a[X+X^2(T+R(X))],
R(X)=X^{-2} integral(X^2u(X) dX) in Xk[X].             (C8)
```

The source changes in `(C7)--(C8)` have Jacobian one, so these potentials
have no mates by THM-3734.  The constant sheets have no mates by THM-3726.

## 3. The full two-variable compatible closure equations

Before any specialization of `u,v`, direct differentiation gives the lower
equation for `(2)`:

```text
4AT^2h_T-(2XT-1+4AT^2v)h_X
 +[2(4A-1)T-4AT^2v_X]h
 +2(A-1)X-A(1+2XT)v_X-2ATv=0.                         (12)
```

Equivalently, if `L_A(h)` denotes the diagonal lower equation `(6)` of
THM-3734, then

```text
L_A(h)-A partial_X[v(1+2XT+4T^2h)]=0.                 (13)
```

The upper equation for `(7)` is

```text
A(1+2XT)h_T-X^2h_X+2(A-1)Xh+2(4A-1)T
 +partial_T[u(2XT-1+X^2h)]=0.                         (14)
```

These identities retain arbitrary bivariate right parameters and every
zero/constant boundary.

## 4. Why lower closure forces `v_X=0`

Suppose `p=deg_X v>=1`; let `m=deg_Xh` when `h!=0`.  Coefficient comparison
in `(12)` exhausts the following cases.

- If `p,m>=2`, the `X^(p+m-1)` coefficient has the unique term

  ```text
  -4AT^2(p+m) lead_X(v)lead_X(h),                      (15)
  ```

  which is nonzero.

- If `p>=2,m=1`, the `X^p` coefficient forces
  `lead_X(h)=-1/(2T)`, not a polynomial.  If `m=0` (including `h=0`), the
  degree-`p` forcing term is uniquely nonzero.

- If `p=1,m>=2`, write `s(T)=lead_X(v)` and
  `c(T)=lead_X(h)`.  The `X^m` coefficient is

  ```text
  2AT c'+(4A-1-m)c-2AT(m+1)s c=0.                    (16)
  ```

  The last product has strictly larger `T`-degree than every other nonzero
  term, so `(16)` is impossible.

- If `p=m=1`, set `zeta(T)=1+2T lead_X(h)`.  The `X` coefficient becomes

  ```text
  T zeta'=[2T lead_X(v)-(A-1)/A]zeta.                 (17)
  ```

  Since `zeta(0)=1`, evaluation at zero forces `A=1`.  At `A=1`, the
  nonzero product `2 lead_X(v)zeta` has degree larger than `zeta'`, again
  impossible.  For `p=1,m=0`, the `X` coefficient says
  `2AT lead_X(v)=A-1`, which is equally impossible.

Therefore every lower survivor has `v_X=0`, with no degree bound.

## 5. The typed dual gate forces `u_T=0`

Now suppose `p=deg_Tu>=1` in `(14)` and put `m=deg_Th`.  The same conclusion
does not rely on an informal symmetry; the typed cases are:

- `p,m>=2`: the top term in `X^2 partial_T(uh)` is unique;
- `p>=2,m=1`: `lead_T(h)=-2/X`, impossible;
- `p>=2,m=0`: the forcing degree is unmatched;
- `p=1,m>=2`: for `s(X)=lead_T(u)` and `c(X)=lead_T(h)`, the top coefficient
  is

  ```text
  -X^2c'+2[A(m+1)-1]Xc+(m+1)X^2s c=0,                (18)
  ```

  whose product term has uniquely highest `X`-degree;
- `p=m=1`: the `T` coefficient is

  ```text
  -X^2c'+(4A-2)Xc+2sX^2c+2(4A-1)+4Xs=0.             (19)
  ```

  Evaluation at `X=0` forces `A=1/4`, after which the product term is still
  uniquely highest.  The `p=1,m=0` coefficient comparison also fails.

Hence an upper survivor necessarily has `u_T=0`.

## 6. The right shear is an integrable source connection

Take `v=v(T)`.  The triangular Euler map

```text
B |-> B'+B/(2AT)                                      (20)
```

sends `T^(j+1)` to `[j+1+1/(2A)]T^j`.  Away from the exceptional slopes

```text
A=-1/[2(j+1)],                                        (21)
```

it therefore has a unique inverse from `k[T]` to `Tk[T]`.  Let `B` be that
inverse and set

```text
Y=X+B,                     h=B/(2T)+g(Y,T).            (22)
```

Exact substitution turns `(12)` into the diagonal lower equation for `g`.
One can see the cancellation directly: the unsheared second component differs
from the `Y`-diagonal component by

```text
-B K(YT)/(2aT),                                       (23)
```

while the right shear contributes `avK(YT)` and the source derivative of
`aY Phi(YT)` contributes `aB'K(YT)`.  Their equality is exactly `(20)`.

No closure hides at an exceptional slope `(21)`.  If `m=deg_Xh>=2`, its
leading coefficient `c(T)` obeys

```text
2ATc'+(4A-1-m)c=0.                                    (24)
```

A polynomial monomial `T^j` in `(24)` forces

```text
A=(m+1)/[2(j+2)],                                     (25)
```

which cannot equal the negative value `(21)`.  The cases `m<=1` force
`A=1`, also nonexceptional.  Thus every actual closure admits `(22)`, and
THM-3734 yields precisely `(3)--(6)`.  At those slopes,
`1/(2A)=r/(r+1)`, giving the denominator in `(4)`.

Dually, for `u=u(X)` the connection equation is

```text
C'+2AC/X=u,                       C in Xk[X].           (26)
```

The exceptional eigenvalues are `A=-(j+1)/2`.  If `m=deg_Th>=2`, the top
coefficient equation

```text
-Xc'+2[A(m+1)-1]c=0                                (27)
```

would force `A=(q+2)/[2(m+1)]` for a monomial `X^q`, never exceptional;
`m<=1` forces `A=1/4`.  Therefore `(26)` is available for every survivor.
With

```text
V=T+C,                         h=2C/X+g(X,V),          (28)
```

equation `(14)` becomes the diagonal upper equation, proving `(8)--(11)`.

Both source translations have unit Jacobian:

```text
J(Y,T)=1,                       J(X,V)=1.               (29)
```

THM-3734 classifies the diagonal equations and excludes arbitrary polynomial
mates.  Source automorphisms preserve mate existence, so every survivor here
also has no mate.  **QED.**

## 7. The inheritance ladder and remaining route

The formulas glue all earlier boundaries exactly.

- `v=0` or `u=0` gives THM-3734.
- For constant `v_0`,
  `B=[2A/(2A+1)]v_0T`, which is the upper-triangular gauge of THM-3736 with
  matrix entry `b=av_0`.  Constant `u_0` gives
  `C=u_0X/(2A+1)`, the lower-triangular gauge with `c=a^{-1}u_0`.
- At `r=1`, the divided-power correction `H_1` vanishes, but `B/(2T)` or
  `2C/X` can remain nonconstant.  This is exactly the full one-variable
  Broughton survivor of THM-3721 and its dual.

Thus the first counterexample-directed escape beyond the constant orbit is
real but rigid: a one-variable right polynomial changes only the source
section of the same cyclotomic potential.  To escape the no-mate class, the
right word must use at least two genuinely interacting polynomial factors or
a different non-elementary core; mere bivariate dependence in this single
right factor is killed by `(15)--(19)`.

Reproduce with

```bash
python3 -B 04-computation/jc2_automorphic_cohn_one_variable_right_towers_thm3740.py
python3 -B -O 04-computation/jc2_automorphic_cohn_one_variable_right_towers_thm3740.py
```

The companion derives `(13)--(14)` symbolically; checks four right profiles,
both raw gradients, both connection defects, and both source translations
through depth eight; rejects 120 mixed-variable low-degree systems; and
checks 48 exceptional negative-resonance systems plus 384 exact leading-ODE
separations.  It also derives both crossed PDEs and verifies 22 positive
crossed boundary controls.  The displayed degree arguments carry the
unbounded proof.
