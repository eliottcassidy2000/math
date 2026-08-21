---
id: THM-3607
title: "Russell-cylinder rank-two linear-projection formal rigidity and degree-seven gate"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; PENDING INDEPENDENT AUDIT.
  For every polynomial cylinder graph w=h(x,q), no rank-two linear projection
  of the three polynomial coordinates (B,C,S) has a nonzero constant ordinary
  Jacobian.  This conclusion does not require the graph to retain the
  THM-3561 collision.  Linear algebra reduces the only hard branch to
  L=B+rho C, M=S+tau C.  Locally along the x=0 arm, a constant Jacobian gives
  a first-order equation whose unique formal solution is
  h=t-rho-ac; it has a genuine simple pole on D=0.  Independently, ordinary
  leading forms for the (B,S) branch force any hypothetical graph of degree
  at least three to have degree divisible by seven, while a complete
  full-collision-preserving degree-at-most-two calculation ends in -81/2.
source: root / cylinder_mixed_projection mixed-projection wildcard, 2026-08-21
audit: >
  PENDING INDEPENDENT AUDIT.  The exact companion verifies the Russell graph
  identities, rank-two row reduction, the fixed-C and (B,C) faces, the
  determinant-to-PDE sign, the explicit arm-completion inverse, formal
  recurrence, ordinary initial-form bracket and degree-seven arithmetic, the
  complete quadratic collision elimination, and the sharp punctured graph.
  Normal and optimized replays are byte-identical to the stored output.  The
  all-degree conclusion is proof-driven and is not promoted before a separate
  hostile proof audit.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
related:
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
script: 04-computation/jc2_russell_cylinder_mixed_projection_thm3607.py
output: 05-knowledge/results/jc2_russell_cylinder_mixed_projection_thm3607.out
script_sha256: 59214827b8361378f2b8952c2cb0850584bca3305e6a3d21d376e2cd21542ebb
output_sha256: 4ca192d1047b42ada9e672e8660e4635c32fb4a69d2c5cb69c26f3fc903b87b2
hash_basis: raw LF bytes
---

# THM-3607 -- Russell-cylinder rank-two linear-projection formal rigidity

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; PENDING INDEPENDENT
AUDIT.**  The theorem closes every rank-two linear projection using the
three coordinates `(B,C,S)` on a polynomial graph.  It does not rule out an
implicit non-graph `A2`, a nonlinear projection, or a projection involving
the remaining cylinder coordinate `Y`.

All rings and derivatives below are over `C`.

## 0. Statement and preservation/loss ledger

Retain the THM-3561 functions

```text
D=1+x^2q,
a=q/D^2,
c=xD(D+2),
b=(D-1)(D+2)^2=ac^2,
e=q(D+3)=a(b+4).                                      (1)
```

Thus `Jac_(x,q)(a,c)=-3` on `D!=0`.  For an arbitrary polynomial
`h in C[x,q]`, put

```text
B_h=b+ch,
S_h=((b+2)(e+3h^2)+ch(3e+h^2))/8.                    (2)
```

These are the `B` coordinate and the stable `S` coordinate of the Russell
cylinder formula after restricting to the graph `w=h`; also put `C=c`.
For every matrix `Lambda in Mat_(2x3)(C)` of row rank two, define

```text
(F_Lambda,G_Lambda)^T=Lambda(B_h,C,S_h)^T.
            Jac_(x,q)(F_Lambda,G_Lambda) notin C*.     (3)
```

In fact `(3)` holds without any collision hypothesis.  Therefore it holds
in particular when

```text
h(0,-3/4)=h(1,-3)=h(-1,-3),                           (4)
```

which is exactly the condition that the polynomial graph retain the **full
four-coordinate stable** triple collision.  For the projected `(B,S)` pair
alone, equality of the three squares `h(r_i)^2` is enough; Section 5 uses
the stronger full-stable condition `(4)` explicitly.

The information ledger is

| item | retained or lost |
|---|---|
| source object | the polynomial graph `w=h(x,q)` in the stabilized THM-3561 source |
| target map | an arbitrary rank-two linear projection of `(B_h,C,S_h)` |
| retained | polynomiality of both outputs; condition `(4)` retains the full stable collision |
| discarded | at least one of `B,C,S` and always `Y`; stable three-volume does not imply planar area |
| tested predicate | a nonzero constant ordinary two-dimensional Jacobian |
| decisive sidecar | the formal completion of the `x=0` arm in the etale `(a,c)` coordinates |
| sharp survivor | on the mixed branch, `h=t-rho-ac` gives Jacobian `3t` and retains the collision |

The last row is not a polynomial counterexample: it has the same `D=0` pole
as the original rational Keller pair.

## 1. Rank-two reduction and the three projection faces

Work first in the function field on `D!=0` and define

```text
A=ac+h.                                                (5)
```

Direct use of `b=ac^2` and `e=a(b+4)` in `(2)` gives

```text
B_h=cA,
S_h=a+3A^2/4+cA^3/8.                                  (6)
```

The point is that `A` need not be polynomial even though `(2)` is
polynomial.

Let `W` be the two-dimensional row space of a rank-two projection of
`(B,C,S)`.  Its intersection with `span(B,C)` is nonzero; choose a line in
that intersection and take a first output on it.

### 1.1 The fixed-`C` line

If the line is `C`, a nonzero constant Jacobian would give a polynomial
Jacobian mate for

```text
c=xp(u),       u=x^2q,       p(u)=(u+1)(u+3).          (7)
```

No such mate exists.  Give `x,q` weights `1,-2`.  The weight-zero part of
`Jac(F,c)` can only come from

```text
F_(-2)=x^(-2)R(u),                  R in u C[u].        (8)
```

For a homogeneous expression `x^mR(u)`, direct differentiation gives

```text
Jac(x^mR(u),c)=x^(m+2)(mRp'-R'p).                      (9)
```

Thus `Jac(F,c)=lambda in C*` would imply

```text
pR'+2p'R=-lambda,
(p^2R)'=-lambda p.                                    (10)
```

Writing `I=u^3/3+2u^2+3u`, integration gives
`p^2R=-lambda I+K_0`.  Evaluation at `u=-3` forces `K_0=0`, while at
`u=-1` the right side is `4lambda/3`, a contradiction.

### 1.2 The `(B,C)` plane

Otherwise normalize the first output to

```text
L=B_h+rho C,                         rho in C.          (11)
```

If the second output also has zero `S` coefficient, the pair is
`GL_2(C)`-equivalent to `(B_h,C)`.  But

```text
Jac(B_h,c)=Jac(b+ch,c)=c(-3c+Jac(h,c)),                (12)
```

which is divisible by the nonconstant polynomial `c` and cannot be a
nonzero constant.

### 1.3 The genuinely mixed line

It remains to consider a second output with nonzero `S` coefficient.  Row
reduction using `L` puts it in the unique normalized form

```text
M=S_h+tau C,                          tau in C.         (13)
```

Put

```text
Z=A+rho,                 A=Z-rho,
K_rho(Z,c)=Z(3A/2+3cA^2/8)-cA^3/8.                   (14)
```

Then `L=cZ`.  Regarding `Z` locally as a function of the etale coordinates
`(a,c)`, differentiation of `(6),(13)` gives

```text
det partial(L,M)/partial(a,c)
 =-(Z+c partial_c Z+(K_rho(Z,c)-tau c) partial_a Z).   (15)
```

Since `Jac_(x,q)(a,c)=-3`, the ordinary source determinant is

```text
Jac_(x,q)(L,M)
 =3(Z+c partial_c Z+(K_rho(Z,c)-tau c) partial_a Z).   (16)
```

This fixes the sign and reduces the only hard projection face to a scalar
transport equation.  Additive constants in the two outputs do not change
the argument.

## 2. The `x=0` arm forces constant boundary data

Suppose for contradiction that the left side of `(16)` is the nonzero
constant `3t`, where `t in C*`.  On the divisor `x=0`, formulas `(1)` give

```text
D=1,             a=q,             c=0,
Z(a,0)=f(a):=h(0,a)+rho.                              (17)
```

Restricting `(16)` to this divisor yields

```text
f+(3/2)(f-rho)f f'=t.                                 (18)
```

Here `f in C[a]`.  If `deg f=n>=1`, the derivative term in `(18)` has
degree `3n-1>n` and nonzero leading coefficient.  It cannot be cancelled by
`f` or by the constant right side.  Hence `f` is constant, and `(18)` then
forces `f=t`.  This uses no collision hypothesis.

## 3. Formal-arm uniqueness forces the punctured solution

The completion used next has an explicit inverse, not merely a pointwise
etale description.  In `C[a][[c]]`, there is a unique `R=1+O(c^2)` solving

```text
ac^2=(R-1)(R+2)^2,                                    (19)
```

because the derivative of the right side at `R=1` is `9`.  Then

```text
x=c/(R(R+2)),              q=aR^2                     (20)
```

satisfy `x^2q=R-1` and `D=R`.  Conversely `(a,c)` are given by `(1)`, so
`(19),(20)` prove the formal line-completion isomorphism

```text
C[a][[c]]  ~=  C[q][[x]].                             (21)
```

In particular `Z` has an expansion in `C[a][[c]]`.  Set `W=Z-t`.  The
boundary result says `c` divides `W`, and `(16)` becomes

```text
W+c partial_c W
 +(K_rho(t+W,c)-tau c) partial_a W=0.                 (22)
```

If `W` were nonzero, let

```text
W=z_N(a)c^N+O(c^(N+1)),       N>=1,
z_N in C[a] nonzero.                                   (23)
```

The coefficient of `c^N` in `(22)` is

```text
(N+1)z_N+(3/2)t(t-rho)z_N'=0.                         (24)
```

The `tau` term starts only in order `N+1`.  In `(24)`, the first term has
degree `deg z_N`, while the derivative term has strictly smaller degree (or
is zero).  Thus `(24)` has no nonzero polynomial solution.  This contradicts
`(23)`, so `Z=t` in the completion.  The localized source domain embeds in
its `x`-adic completion by Krull intersection, hence `Z=t` in the function
field.  Equations `(5),(14)` now force

```text
h=t-rho-ac=t-rho-xq(D+2)/D.                           (25)
```

But `D` does not divide `xq(D+2)`: modulo `D`, one has `x^2q=-1`, so `x`
is a unit and the numerator is `-2/x`, nonzero.  Thus `(25)` has a genuine
simple pole along `D=0` and cannot equal a polynomial `h`.  Together with
Sections 1.1--1.2 this proves `(3)` for every rank-two linear projection.

## 4. Independent ordinary initial-form gate

There is a useful global degree warning independent of the formal proof.
Use ordinary total degree and write

```text
c_7=x^5q^2,                 h_d=top homogeneous part of h.           (26)
```

If `d=deg h>=3`, the terms `rho c` and `tau c` are lower order, so the
unique top forms in the normalized mixed branch `(11),(13)` are

```text
L_(d+7)=c_7 h_d,
M_(3d+7)=c_7 h_d^3/8.                                  (27)
```

There is no tie at `d=3`: the competing degrees for `S_h` are `2d+9`,
`d+11`, and `13`, all strictly below `3d+7`.  The top Jacobian is

```text
-c_7 h_d^3 Jac(h_d,c_7)/4.                            (28)
```

Therefore a constant Jacobian would require `Jac(h_d,c_7)=0`.  If

```text
h_d=sum_(i=0)^d alpha_i x^i q^(d-i),                  (29)
```

then the coefficient attached to the `i`th monomial in that bracket is

```text
(7i-5d) alpha_i.                                      (30)
```

Distinct `i` give distinct monomials, so a nonzero `h_d` can survive only
when `i=5d/7`.  Consequently

```text
7 divides d,             h_d=kappa (x^5q^2)^(d/7).    (31)
```

Thus ordinary leading degree alone eliminates every `d>=3` not divisible
by seven.  Degree seven is the first apparent boundary; Sections 2--3 show
that its cancellation cannot extend to a polynomial solution.

## 5. Complete collision-preserving quadratic hostile

For an independent finite boundary calculation, specialize to the original
`rho=tau=0` pair `(B_h,S_h)` and write

```text
h=A_2 x^2+B_2 xq+C_2 q^2+D_2 x+E_2 q+F_2.            (32)
```

Condition `(4)` is equivalent to

```text
D_2=3B_2,                 16A_2+135C_2-36E_2=0.      (33)
```

After imposing `(33)`, three coefficients of
`Jac_(x,q)(B_h,S_h)` give a complete contradiction.  In order, they are

```text
[x^17 q^3] = -6561(15C_2-4E_2)^4/65536,
[x^13 q^7] = 3(B_2+1)^4/4                 after E_2=15C_2/4,
[x^13 q^3] = -81/2                        after B_2=-1.               (34)
```

The first line forces `E_2=15C_2/4`, which with `(33)` forces `A_2=0`.
The second forces `B_2=-1`, and the last is impossible.  Hence no
collision-preserving graph of total degree at most two can work.  Constants
and affine-linear graphs are included by setting the quadratic coefficients
to zero.

## 6. Sharp punctured endpoint

The formal solution rejected only by polynomiality is genuinely sharp on
every normalized mixed branch.  For `t in C*` and arbitrary `rho,tau in C`,
take on the punctured `(a,c)` plane

```text
h=t-rho-ac,             A=t-rho,             Z=t.     (35)
```

Then

```text
L=tc,
M=a+3(t-rho)^2/4+c((t-rho)^3/8+tau),
Jac_(x,q)(L,M)=3t.                                     (36)
```

At the three THM-3561 collision points, `ac=0`, so the stable graph value is
the same `t-rho` and `(36)` preserves the full collision.  At `rho=tau=0`
this is the `mu=0` endpoint outside the successful arm-pair family recorded
in THM-3605: the arm coordinate `A` becomes constant, but the mixed `S`
projection restores a Darboux pair.  Formula `(25)` shows exactly why it
does not extend across the source puncture.

## 7. Scope and cheapest controls

The theorem closes polynomial graphs over the fixed source `A2_(x,q)` for
every rank-two linear projection of `(B,C,S)`.  It leaves **OPEN**:

- implicitly embedded planes in the cylinder which are not graphs over
  `(x,q)`;
- nonlinear projections of `(B,C,S)` and projections involving `Y`;
- other cylinder isomorphisms or a change of the fixed source coordinates.

The cheapest positive mixed-branch control is `(35)`, which satisfies the
desired Jacobian and collision but exposes the pole.  On the easy faces,
the visible controls are the fixed-`c` weight ODE `(10)` and the factor `c`
in `(12)`.  The degree-seven initial form `h=x^5q^2` is the first control
that defeats only the coarse top-degree gate; formal-arm rigidity still
kills it.

## 8. Exact companion and reproduction

The deterministic companion verifies all displayed algebraic identities and
coefficient rows, finite initial-form controls, and the sharp punctured
endpoint.  The arbitrary-degree conclusion comes from the polynomial degree
argument in `(18)` and the first-nonzero-coefficient proof `(23)--(24)`.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_mixed_projection_thm3607.py
python3 -O 04-computation/jc2_russell_cylinder_mixed_projection_thm3607.py
```

Both streams must be byte-identical to
`05-knowledge/results/jc2_russell_cylinder_mixed_projection_thm3607.out`.
