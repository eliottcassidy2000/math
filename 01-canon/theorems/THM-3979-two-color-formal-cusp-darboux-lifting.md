---
id: THM-3979
title: "Two-color formal cusp Darboux lifting"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT. On the
  minimal height-two completion of THM-3973, the x-adic completion splits
  into the disjoint boundary-cusp and interior-arm colors. The cusp color
  admits an exact formal Darboux pair with boundary value (y^2,y^3): after
  Y=y-x/2, one transverse quadrature and a formal square root put the source
  volume in the standard form X dX wedge dY. The arm color is ordinary
  Darboux. Hensel idempotents glue the two pairs, so the simultaneous
  boundary invoices lift to every finite order. The cusp jet operator has
  image (y) and one-dimensional cokernel, but its preceding kernel jet pays
  the next constant row by the nonzero scalar 3(r+1)L; hence there is no
  hidden formal compatibility obstruction. An explicit polynomial control
  reaches order five on both colors. This is only a formal/all-finite-order
  result: termination in B_2, global critical-point avoidance, finiteness,
  and a Keller pair remain OPEN.
source: jc-zero-debt-lift / post-THM-3973 simultaneous boundary-jet lane, 2026-08-24
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
related:
  - THM-3977-simultaneous-cusp-arm-family-critical-resultant
  - THM-3978-linear-seam-submersion-rational-mate-pole-obstruction
script: 04-computation/jc2_two_color_formal_cusp_lifting_thm3979.py
output: 05-knowledge/results/jc2_two_color_formal_cusp_lifting_thm3979.out
script_sha256: 9b6513f06a55c44fc1d5d677c4d48c3caa316ed7bcfd8e75a398ce47bab70a79
output_sha256: 09dbf58e85137a73aba69837c8fe1e55b54afbcf8a7a703be5b522a811000dd4
semantic_sha256: a325289714ede662300ee6fff5369641bff8c92dc9e8508775200089dce7b6a4
hash_basis: raw LF bytes
---

# THM-3979 -- the simultaneous cusp/arm invoices are formally unobstructed

**PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.** Work over
an algebraically closed field `k` of characteristic zero and choose

```text
6L^2=1.                                                        (1)
```

On `A2_(x,t)` put

```text
z=1+x^2t,             p=zt,             y=xzt^2,
B=B_2=k[x,z,p,y],      X=Spec(B).                              (2)
```

Thus `X` is the minimal exact-volume completion of THM-3973. Its principal
divisor `V(x)` has two disjoint colors

```text
D=V(x,z,p) isomorphic to A1_y,
L_1=V(x,z-1,y) isomorphic to A1_p,                              (3)
```

where `D=X minus A2_(x,t)` and `L_1` lies in the source plane. Let
`Bhat` denote the `x`-adic completion. There are elements

```text
Ahat, Chat in Bhat                                             (4)
```

such that

```text
dAhat wedge dChat = dx wedge dt,                               (5)
(Ahat,Chat)|D=(y^2,y^3),
(Ahat,Chat)|L_1=(p,0).                                         (6)
```

Consequently the simultaneous cuspidal-boundary and smooth-source-arm
requirements have no formal obstruction at any order. More precisely, for
every `N>=1` there are `A_N,C_N in B` whose completed Jacobian ratio satisfies

```text
(dA_N wedge dC_N)/(dx wedge dt) = 1 mod x^N Bhat.              (7)
```

The quotient in `(7)` is regular: on `D` both numerator and denominator
have the same simple zero, while on `L_1` the source volume is a generator.
The theorem makes **no** assertion that one compatible sequence terminates
in `B`, that the resulting target map is finite, or that a polynomial
Jacobian pair exists.

## 1. The completed divisor splits into two colors

The determinantal relations of THM-3973 specialize at height two to

```text
z(z-1)=x^2p,             p(z-1)=xy,             zy=xp^2.       (8)
```

Modulo `x`, they give

```text
B/(x) = k[y] times k[p],                                      (9)
```

with the first factor supported at `z=0` and the second at `z=1`. Since
idempotents lift uniquely in an adically complete ring,

```text
Bhat = Bhat_D times Bhat_L.                                   (10)
```

This splitting is also explicit. Put

```text
q=2z-1,                   W=(1+4x^2p)^(1/2),
e_L=(1+q/W)/2,            e_D=(1-q/W)/2.                      (11)
```

The square root is the unique one congruent to one modulo `x`. Equation
`q^2=1+4x^2p` gives

```text
e_D^2=e_D,       e_L^2=e_L,       e_De_L=0,       e_D+e_L=1.  (12)
```

Thus data on the two formal colors can be constructed independently and
then glued without cross terms; in particular `de_D=de_L=0`.

## 2. Exact cusp normal form on the boundary color

On `Bhat_D`, equations `(8)` imply

```text
z(z-1)^2=x^3y,                  p=xy/(z-1).                    (13)
```

The derivative with respect to `z` of the left side of `(13)` is one at
`z=0`. Formal implicit-function theory therefore gives a unique

```text
z=zeta(x,y) in x^3 k[y][[x]],
Bhat_D isomorphic to k[y][[x]].                               (14)
```

Since `t=(z-1)/x^2` in the fraction field, differentiation of `(13)` gives

```text
zeta_y=x^3/((zeta-1)(3zeta-1)),
dx wedge dt=x w(x,y) dx wedge dy,                             (15)
w=1/((zeta-1)(3zeta-1)) in k[y][[x]]^*,       w(0,y)=1.       (16)
```

Now center the cusp parameter by

```text
Y=y-x/2.                                                       (17)
```

At fixed `Y`, define `X=x+O(x^2)` by the single quadrature

```text
X^2=2 integral_0^x s w(s,Y+s/2) ds.                           (18)
```

The right side is `x^2` times a unit congruent to one, so the required
formal square root exists uniquely. Differentiating `(18)` at fixed `Y`
gives

```text
X X_x=x w(x,y),             X dX wedge dY=dx wedge dt.        (19)
```

The standard cusp pair

```text
A_D=Y^2+2LX,                 C_D=Y^3+3LXY                    (20)
```

then has

```text
det partial(A_D,C_D)/partial(X,Y)=6L^2X=X.                    (21)
```

Equations `(19)--(21)` prove `dA_D wedge dC_D=dx wedge dt`
exactly, not merely to bounded order. At `x=0`, `(20)` restricts to
`(y^2,y^3)`. Its first jets are

```text
A_D=y^2+x(2L-y)+x^2/4                         mod x^4,
C_D=y^3+x(3Ly-3y^2/2)+x^2(3y/4-3L/2)-x^3/8  mod x^4.         (22)
```

These are precisely the corrected cusp-side jets that generated the
finite polynomial experiment below.

## 3. Exact Darboux normal form on the source-arm color

On `Bhat_L`, write `v=z-1`. The equation

```text
v(1+v)=x^2p                                                     (23)
```

has a unique solution `v in x^2k[p][[x]]`; hence

```text
Bhat_L isomorphic to k[p][[x]],             t=p/z.             (24)
```

The pair

```text
A_L=t+2Lx,                         C_L=-x                      (25)
```

satisfies `dA_L wedge dC_L=dx wedge dt` and restricts to `(p,0)`
on `L_1`. Gluing `(20)` and `(25)` with `(11)` proves `(4)--(6)`.

## 4. Jet anatomy: one cokernel, paid one order later

The coordinate construction already proves all-order solvability. The
linearized jet calculation explains the apparent obstruction seen during
direct polynomial correction.

On the arm color, start from `(25)`. Adding `x^r f(t)` to `C_L` changes the
coefficient of `x^(r-1)` in the Jacobian by

```text
-r f(t).                                                       (26)
```

Thus the arm jet operator is surjective in every positive order.

On the cusp color, write

```text
A=y^2+x a_1(y)+...,          C=y^3+x c_1(y)+... .             (27)
```

The order-zero determinant equation forces the first relation below. The
constant part of the order-one equation is
`(3/2)a_1(0)^2=1`; choose the sign

```text
c_1=(3/2)y a_1,                  a_1(0)=2L!=0.                (28)
```

At order `r>=2`, new jets `x^r(a,c)` have leading response

```text
T_r(a,c)=r y(3ya-2c).                                          (29)
```

Hence

```text
image(T_r)=(y),       coker(T_r)=k,
ker(T_r)={(a,(3/2)ya):a in k[y]}.                              (30)
```

The missing constant is not an obstruction. If the kernel jet in `(30)`
is added at order `r`, its contribution to the constant part of the next
Jacobian row is

```text
(3/2)(r+1)a_1(0)a(0)=3(r+1)L a(0).                            (31)
```

This scalar is nonzero in characteristic zero. One chooses `a(0)` to pay
the next constant row, and the next jet kills the remaining multiple of
`y`. Induction proves the same all-order conclusion as `(18)--(21)` and
identifies the delayed payment mechanism exactly.

## 5. From the formal pair to every finite algebraic jet

For each `M`, completion gives

```text
B/(x^M) isomorphic to Bhat/(x^M).                              (32)
```

Reduce `(Ahat,Chat)` modulo `x^(N+2)` and choose arbitrary lifts in `B`.
Differentiation loses at most one `x`-adic order. On `L_1` the volume is a
unit; on `D` equation `(16)` says it is `x` times a unit. Therefore the two
extra function jets give the relative congruence `(7)`. This proves that
every finite system of the simultaneous boundary-Jacobian equations is
solvable by actual elements of `B`.

This inverse-limit statement is deliberately weaker than algebraization.
The degrees in `p` and `y` of direct corrections can grow, and no uniform
degree bound or finite-support termination is supplied by completion.

## 6. Exact finite polynomial control

The companion starts from the corrected global cusp/arm seed and applies
the color-localized corrections `z f(p)` and `(1-z)g(y)`. It verifies an
explicit pair `A^(4),C^(4) in B` for which the Jacobian error begins at
order five on both colors. The first surviving rows are

```text
on L_1:
 -(2L/15)(324p^5+675p^4+330p^3-135p^2+51p-4)x^5,             (33)

on D:
 (7/24)(144Ly^4-16Ly^2-27L-288y^5+224y^3)x^5.                (34)
```

At the preceding stage the cusp row had a nonzero constant. The exact
kernel payment was

```text
delta A=(5L/6)x^5(1-z),                                      (35)
```

with its coordinated `C` jet; it raised both colors by two orders. This is
a hostile control against mistaking the cokernel in `(30)` for a genuine
formal obstruction. The complete finite formulas are frozen in the script,
where symbolic substitution into both formal charts checks `(33)--(35)`.

## 7. Consequence for counterexample design

THM-3977 proves that several natural finite polynomial truncations of these
same invoices still have affine critical points. The present theorem shows
why adding more boundary jets alone cannot settle that global question:
the formal problem is exactly soluble. Any decisive obstruction must see
failure of algebraization/termination, a generic-fibre differential, an
interior critical locus, or finiteness of the target map. Conversely, a
positive construction must compress the infinite formal pair into finite
support without reintroducing an interior critical point.

## 8. Exact verification

Run

```bash
python3 04-computation/jc2_two_color_formal_cusp_lifting_thm3979.py
python3 -O 04-computation/jc2_two_color_formal_cusp_lifting_thm3979.py
```

The two outputs must byte-match. The companion checks the determinantal
relations, both completed charts, explicit Hensel idempotents, the formal
implicit series and volume unit, the quadrature/square-root identity, the
standard cusp determinant, both jet operators and delayed scalar, and the
explicit order-five polynomial control.
