---
id: THM-3611
title: "Russell-cylinder central-transverse nonlinear first-coordinate rigidity"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT
  CONTROLS; NOT YET PROVED CANON.  Let P(T,C) and K(U,C) be arbitrary
  polynomials.  If P_T(0,0) is nonzero, then no polynomial source graph followed by
  L=P(B,C), M=S+K(L,C) has nonzero constant ordinary Jacobian.  The central
  formal arm forces the localized quotient (P(B,C)-P(0,0))/C to be constant;
  after clearing C, B is algebraic over C(C).  The birational coordinate
  v=x^2q proves C(x,q)=C(C,v), hence B lies in C(C); polynomial intersection
  then puts B in C[C], contradicting its values 0 and -4 over two components
  of C=0.  The arm-collapsing coordinate C+B(B+4) has an independent exact
  square-identity proof.  Vanishing P_T(0,0), first outputs involving S or Y,
  and implicit non-graph planes remain open.
source: root/nonlinear_cylinder_shears maximal nonlinear-coordinate envelope, 2026-08-21
audit: PENDING -- provisional package requires independent hostile audit
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3608-russell-cylinder-nonlinear-target-shear-rigidity
related:
  - THM-3607-russell-cylinder-mixed-projection-degree-seven-gate
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
script: 04-computation/jc2_russell_cylinder_arm_separating_nonlinear_first_coordinate_thm3611.py
output: 05-knowledge/results/jc2_russell_cylinder_arm_separating_nonlinear_first_coordinate_thm3611.out
script_sha256: ac0a85269a0eb621907d769c95a3577ee6e23a16fbb31c300a997a2a9c88fd69
output_sha256: 56a302aaa9f67ebe345d76b0c24348b37212f097a4b72b063f90afcaee9794ce
hash_basis: raw LF bytes
---

# THM-3611 -- Russell-cylinder central-transverse nonlinear first-coordinate rigidity

**RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT.**  This is a proposed
strict extension of THM-3608 from first outputs affine in `B` over `C[C]` to
arbitrary polynomials in `(B,C)` transverse to `B` at the central arm.  It is not part of the proved dependency
graph until an independent audit promotes it.  All rings and derivatives are
over `C`.

## 0. Statement and quotient/descent ledger

Retain the THM-3561 functions

```text
D=1+x^2q,
a=q/D^2,
c=xD(D+2),
b=(D-1)(D+2)^2=ac^2,
e=q(D+3)=a(b+4),                    Jac_(x,q)(a,c)=-3.  (1)
```

For an arbitrary polynomial graph `h in C[x,q]`, set

```text
A=ac+h,
B_h=b+ch=cA,
S_h=a+3A^2/4+cA^3/8.                                  (2)
```

The last two functions are polynomial in `(x,q)`.  Let

```text
P(T,c) in C[T,c],       K(U,c) in C[U,c],
p_0=P(0,0),             alpha=P_T(0,0),
alpha!=0,                                                   (3)

L_h=P(B_h,c),           M_h=S_h+K(L_h,c).              (4)
```

Then the provisional theorem is

```text
Jac_(x,q)(L_h,M_h) notin C*.                           (5)
```

No collision hypothesis is used, so `(5)` remains true under the full-stable
collision condition on `h`.  The data ledger is

| item | retained or lost |
|---|---|
| source | a polynomial graph in the stabilized THM-3561 source |
| target pair | `(P(B,C), S+K(P(B,C),C))` |
| retained | polynomiality of both outputs; any full stable collision already retained by the graph |
| tested predicate | a nonzero constant ordinary planar Jacobian |
| local sidecar | `P_T(0,0)!=0`, which makes the central-arm first coordinate formally invertible |
| global sidecar | `C(x,q)=C(c,v)` for `v=x^2q`, relative algebraic closure, and `C[x,q] intersect C(c)=C[c]` |
| controlled loss | the normalized quotient `Z` below is generally localized, not a global polynomial |
| sharp local survivor | the formal algebraic branch `A=psi(t,c)`, giving determinant `3t` |

The quotient warning is load-bearing.  The proof never claims that
`(P(B_h,c)-p_0)/c` lies in `C[x,q]`.  It is represented by a polynomial in
the **formal arm variables** `(A,c)`, hence lies in the `D`-localized source
and its central completion.  Only after formal rigidity is it multiplied by
`c` to obtain a global polynomial identity; an algebraic-descent argument,
not division by `c`, then supplies the global contradiction.

## 1. The formal-arm polynomial and its inverse

Define

```text
Xi(A,c)=(P(cA,c)-p_0)/c.                               (6)
```

Every nonconstant monomial of `P(cA,c)-p_0` contains `c`, so

```text
Xi in C[A,c].                                          (7)
```

Put

```text
r=P_c(0,0),          Z=Xi(A,c),          L_h=p_0+cZ.   (8)
```

Taylor expansion at `c=0` gives

```text
Xi(A,0)=alpha A+r.                                     (9)
```

Since `alpha!=0`, formal inversion of `(9)` gives a unique

```text
A=psi(Z,c) in C[Z][[c]],
psi(Z,0)=(Z-r)/alpha.                                  (10)
```

This needs no analytic implicit-function theorem.  Comparing successive
`c` coefficients in `Xi(A,c)=Z`, the coefficient of each new unknown is the
nonzero constant `alpha`; induction puts every coefficient of `psi` in
`C[Z]`.

On the source, one may also write

```text
Z=(P(B_h,c)-p_0)/c,                                   (11)
```

but `(11)` is a quotient in the function field.  Its regular central-arm
meaning is `(6)`, with `A=ac+h` in `C[x,q,D^(-1)]`.  Condition
`alpha!=0` does not make this quotient a global polynomial.  Treating `Z`
as one would silently erase components of `c=0`; Section 4 instead clears
the quotient before making any global deduction.

## 2. Exact transport equation for arbitrary `P`

Put

```text
Q(A,c)=3A/2+3cA^2/8.                                  (12)
```

Regard `A=psi(Z,c)` in the central completion.  Then

```text
A_a=psi_Z Z_a,        A_c=psi_Z Z_c+psi_c.             (13)
```

In `K(U,c)`, let `K_c` mean the derivative in the second argument with `U`
fixed.  Differentiating `(2),(4),(8)`, all `K_U` terms cancel and gives

```text
det partial(L_h,M_h)/partial(a,c)
 =-(Z+cZ_c+V(Z,c)Z_a),                                (14)

V(Z,c)=Z Q(A,c) psi_Z
       -c(A^3/8+Q(A,c)psi_c+K_c(p_0+cZ,c)).            (15)
```

Since `Jac_(x,q)(a,c)=-3`, the source determinant is

```text
Jac_(x,q)(L_h,M_h)=3(Z+cZ_c+V(Z,c)Z_a).               (16)
```

For `P(T,c)=u(c)T+F(c)`, equations `(6),(10),(15)` specialize exactly to
the THM-3608 transport equation.  Thus no sign or normalization changes at
the affine boundary.

## 3. Boundary degree and formal uniqueness

Suppose the left side of `(16)` is the nonzero constant `3t`.  On the
central arm `x=0`, one has `a=q,c=0`, and

```text
f(a)=Z(a,0)=alpha h(0,a)+r in C[a].                   (17)
```

Equations `(10),(12),(15)` yield

```text
V(f,0)=3f(f-r)/(2alpha^2).                             (18)
```

Thus `(16)` restricts to

```text
f+3f(f-r)f'/(2alpha^2)=t.                             (19)
```

If `deg f=d>=1`, the derivative term in `(19)` has degree `3d-1>d` and
nonzero leading coefficient.  Hence `f` is constant and `(19)` forces

```text
f=t.                                                   (20)
```

Use the explicit formal-arm isomorphism from THM-3607/3608,

```text
C[a][[c]] ~= C[q][[x]],                               (21)
```

and write `W=Z-t`.  If

```text
W=z_N(a)c^N+O(c^(N+1)),       N>=1,       z_N!=0,     (22)
```

then the first nonzero coefficient in `(16)/3-t=0` is

```text
(N+1)z_N+3t(t-r)z_N'/(2alpha^2)=0.                   (23)
```

The derivative term has smaller `a`-degree than the first term, so `(23)`
has no nonzero polynomial solution.  Therefore

```text
Z=t                                                       (24)
```

in the central completion.  The ring `C[x,q,D^(-1)]` embeds in its `x`-adic
completion because `D` is an `x`-adic unit and the domain is separated.
Thus `(24)` holds in the localized source function field.

## 4. Clear the quotient, descend algebraically, then compare arms

Multiply `(11),(24)` by `c`.  This produces the global polynomial identity

```text
P(B_h,c)=p_0+tc.                                       (25)
```

Because `alpha=P_T(0,0)` is nonzero, the polynomial
`P(T,c)-p_0-tc` is nonconstant in `T`.  Equation `(25)` therefore makes
`B_h` algebraic over `C(c)`.

Put

```text
v=x^2q,              p(v)=(v+1)(v+3),
c=xp(v).                                                 (26)
```

There is a birational change of function-field coordinates

```text
C(x,q)=C(x,v)=C(c,v),
x=c/p(v),              q=v p(v)^2/c^2.                (27)
```

A ground field is relatively algebraically closed in a one-variable pure
transcendental extension: every nonconstant element of `C(c)(v)` is
transcendental over `C(c)`.  Indeed, if `R=f(v)/g(v)` is nonconstant, then
`v` is algebraic over `C(c)(R)` through `f(T)-R g(T)`; if `R` were algebraic
over `C(c)`, transitivity would make `v` algebraic, a contradiction.  Since
`B_h` is algebraic over `C(c)`, `(27)`
forces

```text
B_h in C(c).                                           (28)
```

The polynomial intersection is exact:

```text
C[x,q] intersect C(c)=C[c].                           (29)
```

Indeed, write an element of the left side as `f(c)/g(c)` with coprime
`f,g in C[c]`.  If `g` is nonconstant, choose a root `gamma` of `g`.
Then `f(gamma)!=0`, while the polynomial identity
`f(c)=g(c)H(x,q)` would make the nonconstant polynomial `c-gamma` divide
`f(c)=f(gamma)+(c-gamma)J(c)`, hence divide the nonzero constant
`f(gamma)`, impossible.  Thus `g` is constant and `(29)` follows.

Consequently `B_h=s(c)` for some `s in C[c]`.  But `B_h=b+ch` has two
incompatible values above `c=0`:

```text
(x,q)=(0,0):       D=1, c=0, b=0,       so s(0)=0,
(x,q)=(1,-1):      D=0, c=0, b=-4,      so s(0)=-4.   (30)
```

This contradiction proves `(5)`.  No comparison of `P(0,0)` with
`P(-4,0)` is needed.  The mechanism is instead: central transversality
forces `(25)`, relative algebraic closure descends `B_h` to `C(c)`, and the
two arms obstruct polynomial descent to `C[c]`.

## 5. Sharp local formal-algebraic endpoint

The formal conclusion is exact.  For any `t in C*`, let

```text
A_t(c)=psi(t,c) in C[[c]],          h_t=A_t(c)-ac.     (31)
```

Then in the central formal arm

```text
L_t=p_0+tc,
M_t=a+3A_t(c)^2/4+cA_t(c)^3/8+K(p_0+tc,c),
det partial(L_t,M_t)/partial(a,c)=-t,
Jac_(x,q)(L_t,M_t)=3t.                                (32)
```

This is a local formal-algebraic control, not a polynomial graph.  The
algebraic-descent and two-arm contradiction in Section 4 are exactly what
prevent it from globalizing.

For the simplest genuinely nonlinear example

```text
P(T,c)=T+T^2,                 p_0=r=0, alpha=1,        (33)
```

one has `Xi(A,c)=A+cA^2` and

```text
A_t(c)=(-1+sqrt(1+4ct))/(2c)
      =t-ct^2+2c^2t^3-5c^3t^4+14c^4t^5-... .         (34)
```

The control solves the central equation exactly, while Section 4 prevents
any polynomial graph from realizing it globally.

## 6. The arm-collapsing coordinate: independent square control

A sharp test of the Section 4 descent mechanism is the arm-identifying
polynomial coordinate

```text
P_*(T,c)=c+T(T+4).                                    (35)
```

It is a genuine polynomial coordinate: `(T,c) -> (T,P_*)` has inverse
`c=P_*-T(T+4)`.  Moreover

```text
P_{*,T}(0,0)=4,
P_*(-4,0)=P_*(0,0)=0.                                 (36)
```

Thus central formal rigidity applies although the first coordinate itself
identifies `B=0` and `B=-4`.  The general algebraic descent already closes
this case.  There is also an independent direct proof.  On a graph, set

```text
Y_h=ce+(2b+4)h+ch^2.                                  (37)
```

The Danielewski relation gives

```text
B_h(B_h+4)=cY_h,
L_*=c(1+Y_h).                                         (38)
```

If `(L_*,S_h+K(L_*,c))` had constant Jacobian, Sections 2--3 would force
`1+Y_h=t`, hence `Y_h=k` for a constant `k`.  But `(37)` and
`c^2e=b(b+4)` imply the exact square identity

```text
(ch+b+2)^2=4+cY_h=4+kc.                               (39)
```

If `k!=0`, the right side of `(39)` has ordinary total degree seven, because
the leading form of `c=xD(D+2)` is `x^5q^2`; a nonzero polynomial square has
even total degree.  Hence `k=0`.  Then `(39)` forces

```text
ch=-b       or       ch=-(b+4).                       (40)
```

The first identity is impossible modulo `D`, where `c=0,b=-4`; the second is
impossible modulo `x`, where `c=b=0` but `b+4=4`.  Thus `(35)` is independently
closed for every polynomial graph and every `K`.

The same proof closes the affine family

```text
P(T,c)=p_0+nu*c+mu*T(T+4),           mu!=0,           (41)
```

because its normalized first output is `Z=nu+mu Y_h`.

## 7. Sharp scope and open exits

The provisional result closes every nonlinear `P(B,C)` with
`P_T(0,0)!=0`; no arm-separation hypothesis remains.  It includes first
coordinates which are not polynomial coordinates and coordinates produced
by arbitrary target automorphisms whenever their selected `P(B,C)` is
central-transverse.

The hypotheses mark genuine method boundaries:

- if `P_T(0,0)=0`, the central first coordinate is not formally invertible;
  THM-3608 gives the exact nonconstant formal Keller hostile
  `P(T,c)=cT+c`;
- first outputs involving `S` or `Y`, a nonconstant coefficient of `S` in
  the second output, implicit non-graph planes, other cylinder
  isomorphisms, and changes of the fixed source coordinates remain **OPEN**.

## 8. Exact companion and reproduction

The deterministic companion verifies the Russell arm formulas, polynomial
`Xi` divisions and central jets for nonlinear `P`, the abstract `psi`
transport determinant and `K_U` cancellation, direct polynomial-source
rows, boundary and first-coefficient degree gates, the nonlinear formal
endpoint and its series, the controlled quotient identity, the two-arm
evaluations, the birational `(c,v)` inverse, relative-algebraic/intersection
controls, and the independent square no-go for `(35),(41)`.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_arm_separating_nonlinear_first_coordinate_thm3611.py
python3 -O 04-computation/jc2_russell_cylinder_arm_separating_nonlinear_first_coordinate_thm3611.py
```

Both streams must be byte-identical to
`05-knowledge/results/jc2_russell_cylinder_arm_separating_nonlinear_first_coordinate_thm3611.out`.
