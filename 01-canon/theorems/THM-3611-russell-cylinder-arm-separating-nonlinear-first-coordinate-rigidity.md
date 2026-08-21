---
id: THM-3611
title: "Russell-cylinder arm-separating nonlinear first-coordinate rigidity"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT
  CONTROLS; NOT YET PROVED CANON.  Let P(T,C) and K(U,C) be arbitrary
  polynomials.  If P_T(0,0) is nonzero and P(-4,0) differs from P(0,0),
  then no polynomial source graph followed by
  L=P(B,C), M=S+K(L,C) has nonzero constant ordinary Jacobian.  The central
  formal arm forces the localized quotient (P(B,C)-P(0,0))/C to be constant;
  after clearing C, the hostile D=0 arm contradicts the separation condition.
  The arm-collapsing polynomial coordinate C+B(B+4), which fails separation,
  is nevertheless closed separately by an exact square identity.  Vanishing
  P_T(0,0), general arm-identifying P, first outputs involving S or Y, and
  implicit non-graph planes remain open.
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
script_sha256: 7ffe6fd4088dbdf57f3857871d28ac0abb21eadedc979a5f385d706dad3ea47c
output_sha256: dfac53f79b0b4a266f1af7d4da2b1839995170da78e6230e06dce0515b922eb2
hash_basis: raw LF bytes
---

# THM-3611 -- Russell-cylinder arm-separating nonlinear first-coordinate rigidity

**RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT.**  This is a proposed
strict extension of THM-3608 from first outputs affine in `B` over `C[C]` to
arbitrary polynomials in `(B,C)`.  It is not part of the proved dependency
graph until an independent audit promotes it.  All rings and derivatives are
over `C`.

## 0. Statement and quotient-loss ledger

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
alpha!=0,               P(-4,0)!=p_0,                 (3)

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
| global sidecar | `P(-4,0)!=P(0,0)`, which keeps the hostile `B=-4` arm distinct from `B=0` |
| controlled loss | the normalized quotient `Z` below is generally localized, not a global polynomial |
| sharp local survivor | the formal algebraic branch `A=psi(t,c)`, giving determinant `3t` |

The quotient warning is load-bearing.  The proof never claims that
`(P(B_h,c)-p_0)/c` lies in `C[x,q]`.  It is represented by a polynomial in
the **formal arm variables** `(A,c)`, hence lies in the `D`-localized source
and its central completion.  Only after formal rigidity is it multiplied by
`c` to obtain a global polynomial identity.

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
`P(-4,0)!=p_0` in fact guarantees that the numerator of `(11)` is not
globally divisible by `c`; treating `Z` as a global polynomial would erase
the very hostile arm used in Section 4.

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

## 4. Clear the quotient, then separate the hostile arm

Multiply `(11),(24)` by `c`.  This produces the global polynomial identity

```text
P(B_h,c)=p_0+tc.                                       (25)
```

Now evaluate `(25)` on the prime divisor `D=0`.  There

```text
c=0,             B_h=b+ch=-4,                         (26)
```

because `b=(D-1)(D+2)^2`.  Equation `(25)` would force

```text
P(-4,0)=p_0=P(0,0),                                   (27)
```

contrary to `(3)`.  This proves `(5)`.  The mechanism has two independent
parts: `alpha!=0` supplies central formal rigidity, and `(27)` supplies the
global arm obstruction.

## 5. Sharp local formal-algebraic endpoint

The formal conclusion is exact.  For any `t in C*`, let

```text
A_t(c)=psi(t,c) in C[[c]],          h_t=A_t(c)-ac.     (28)
```

Then in the central formal arm

```text
L_t=p_0+tc,
M_t=a+3A_t(c)^2/4+cA_t(c)^3/8+K(p_0+tc,c),
det partial(L_t,M_t)/partial(a,c)=-t,
Jac_(x,q)(L_t,M_t)=3t.                                (29)
```

This is a local formal-algebraic control, not a polynomial graph.  The arm
separation in `(3)` is exactly what prevents it from globalizing.

For the simplest genuinely nonlinear example

```text
P(T,c)=T+T^2,                 p_0=r=0, alpha=1,        (30)
```

one has `Xi(A,c)=A+cA^2` and

```text
A_t(c)=(-1+sqrt(1+4ct))/(2c)
      =t-ct^2+2c^2t^3-5c^3t^4+14c^4t^5-... .         (31)
```

Here `P(-4,0)=12`, so the local solution cannot cross the hostile arm.

## 6. The arm-collapsing coordinate: method hostile, separately closed

The cheapest failure of the separation hypothesis is

```text
P_*(T,c)=c+T(T+4).                                    (32)
```

It is a genuine polynomial coordinate: `(T,c) -> (T,P_*)` has inverse
`c=P_*-T(T+4)`.  Moreover

```text
P_{*,T}(0,0)=4,
P_*(-4,0)=P_*(0,0)=0.                                 (33)
```

Thus central formal rigidity still applies, but the final arm-separation
step does not.  On a graph, set

```text
Y_h=ce+(2b+4)h+ch^2.                                  (34)
```

The Danielewski relation gives

```text
B_h(B_h+4)=cY_h,
L_*=c(1+Y_h).                                         (35)
```

If `(L_*,S_h+K(L_*,c))` had constant Jacobian, Sections 2--3 would force
`1+Y_h=t`, hence `Y_h=k` for a constant `k`.  But `(34)` and
`c^2e=b(b+4)` imply the exact square identity

```text
(ch+b+2)^2=4+cY_h=4+kc.                               (36)
```

If `k!=0`, the right side of `(36)` has ordinary total degree seven, because
the leading form of `c=xD(D+2)` is `x^5q^2`; a nonzero polynomial square has
even total degree.  Hence `k=0`.  Then `(36)` forces

```text
ch=-b       or       ch=-(b+4).                       (37)
```

The first identity is impossible modulo `D`, where `c=0,b=-4`; the second is
impossible modulo `x`, where `c=b=0` but `b+4=4`.  Thus `(32)` is separately
closed for every polynomial graph and every `K`.

The same proof closes the affine family

```text
P(T,c)=p_0+nu*c+mu*T(T+4),           mu!=0,           (38)
```

because its normalized first output is `Z=nu+mu Y_h`.

## 7. Sharp scope and open exits

The provisional result closes arbitrary nonlinear `P(B,C)` whenever the
central derivative and hostile-arm separation in `(3)` both hold, plus the
collapsed-arm family `(38)`.  It includes first coordinates which are not
polynomial coordinates and coordinates produced by arbitrary target
automorphisms whenever their selected `P(B,C)` meets these tests.

The hypotheses mark genuine method boundaries:

- if `P_T(0,0)=0`, the central first coordinate is not formally invertible;
  THM-3608 gives the exact nonconstant formal Keller hostile
  `P(T,c)=cT+c`;
- if `P(-4,0)=P(0,0)`, the quotient-clearing identity does not distinguish
  the two arms; family `(38)` is closed by its special `Y` sidecar, but a
  general arm-identifying `P` remains **OPEN** here;
- first outputs involving `S` or `Y`, a nonconstant coefficient of `S` in
  the second output, implicit non-graph planes, other cylinder
  isomorphisms, and changes of the fixed source coordinates remain **OPEN**.

## 8. Exact companion and reproduction

The deterministic companion verifies the Russell arm formulas, polynomial
`Xi` divisions and central jets for nonlinear `P`, the abstract `psi`
transport determinant and `K_U` cancellation, direct polynomial-source
rows, boundary and first-coefficient degree gates, the nonlinear formal
endpoint and its series, the controlled quotient identity, the hostile-arm
evaluation, and the complete square no-go for `(32),(38)`.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_arm_separating_nonlinear_first_coordinate_thm3611.py
python3 -O 04-computation/jc2_russell_cylinder_arm_separating_nonlinear_first_coordinate_thm3611.py
```

Both streams must be byte-identical to
`05-knowledge/results/jc2_russell_cylinder_arm_separating_nonlinear_first_coordinate_thm3611.out`.
