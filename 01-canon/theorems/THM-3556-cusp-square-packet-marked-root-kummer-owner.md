---
id: THM-3556
title: "Cusp-square packet as a marked-root/Kummer inverse-cubic owner"
status: >
  PROVED + VERIFIED-EXACT + FINITE-EXACT POSITIVE PACKET.  For every
  polynomial U(v,y), the cusp-square packet (L,T,U,S) is the coefficient
  packet of an inverse cubic L X^3+T X+2U that factors into the marked root
  X=-1/v and a quadratic Kummer pair.  Its discriminant is -4 L S^2, and
  the marked root escapes projectively at v=0.  None of the six natural
  two-coordinate projections is Keller.  One explicit four-coordinate
  packet is nevertheless everywhere immersive, while no constant-linear
  projection has unit Jacobian.  Nonlinear decomposable projection remains
  open; no planar counterexample is claimed.
source: kps-s188
depends_on: []
related:
  - THM-3535-fixed-keller-full-wreath-and-all-level-linear-primitivity
  - THM-3554-punctured-kummer-collision-surface-normal-form
  - THM-3555-catalan-thickening-universal-cubic-root-cover
companion: 04-computation/cusp_square_marked_root_kummer_owner_kps_s188.py
output: 05-knowledge/results/cusp_square_marked_root_kummer_owner_kps_s188.out
script_sha256: e5bf45488fb2e8008bad4fe4693f04d8d7dc25bf5b265f9778e2385111dcc971
output_sha256: 670ac8f38a05eef09da1d836a704cfffdef5c10fa61f6334d7a08932213cc450
hash_basis: LF-normalized bytes
---

# THM-3556 -- cusp-square packet as a marked-root/Kummer inverse-cubic owner

**PROVED + VERIFIED-EXACT + FINITE-EXACT POSITIVE PACKET.**  The
four-coordinate cusp packet is not merely a discriminant identity.  It is an
explicit `1+2` inverse-cubic owner: one marked root and one quadratic Kummer
pair.  The marked root reaches infinity exactly where the leading cubic
coefficient vanishes.

The field has characteristic zero.

## 1. The packet and its inverse cubic

Let `U=U(v,y)` be any polynomial and define

```text
T=y^2-6vU,
S=y^3-9vUy,
L=v^2(8vU-y^2).                                       (1)
```

Then the polynomial identity

```text
S^2=T^3+27LU^2                                        (2)
```

holds.  Associate to `(1)` the cubic

```text
E(X)=LX^3+TX+2U.                                      (3)
```

It has the exact factorization

```text
E(X)=(vX+1)
     [v(8vU-y^2)X^2+(y^2-8vU)X+2U].                  (4)
```

For `v!=0`, the first factor gives the marked root

```text
X_0=-1/v.                                              (5)
```

The factorization is better read projectively.  The homogenized linear
factor is `vX+Z`, so its root is `[-1:v]`.  At `v=0` it becomes `[1:0]`, the
point at infinity.  Thus the marked root genuinely escapes when the cubic
leading coefficient `L`, which contains `v^2`, vanishes.

## 2. Kummer pair and discriminant-square law

The discriminant of the quadratic factor in `(4)` is

```text
Delta_2=y^2(y^2-8vU)=-y^2 L/v^2             in k[v,v^-1,y]. (6)
```

After removing the displayed square `(y/v)^2`, its quadratic square class is
`-L`.  The resultant of the linear and quadratic factors is

```text
Res=-2v(y^2-9vU).                                     (7)
```

Since `S=y(y^2-9vU)`, the product-discriminant formula gives

```text
disc_X(E)=Delta_2 Res^2=-4LS^2.                        (8)
```

The same result follows directly:

```text
disc(LX^3+TX+2U)
  =-4LT^3-108L^2U^2
  =-4L(T^3+27LU^2)
  =-4LS^2.                                             (9)
```

Thus `(2)` is exactly the discriminant-square owner equation.  The sole odd
discriminant factor is `L`; `S^2` records projection collisions between
finite roots.  This is the same typed split used by the fixed Keller inverse
cubic in THM-3535, but `(1)` gives a flexible two-source-parameter packet
rather than asserting that its coefficients already descend to a planar
Keller target.

## 3. Every natural two-coordinate projection is obstructed

Write `Jac` for `Jac_(v,y)`.  Exact differentiation gives

```text
Jac(T,S)=54vU(U+vU_v),                                 (10)
v | Jac(L,T), Jac(L,U), Jac(L,S),                      (11)
Jac(T,U)=-2(yU_v+3UU_y),                              (12)
Jac(U,S)|_(v=0)
 =3y[y U_v(0,y)+3U(0,y)(d/dy)U(0,y)].                 (13)
```

Equations `(10)`, `(11)`, and `(13)` visibly exclude a nonzero constant.
For `(12)`, suppose `U=sum_(i=0)^d a_i(y)v^i`.  If `d>=1`, let `j` be the
largest index with `a_j'!=0`.  The highest `v`-degree involving a derivative
in `UU_y` has coefficient `a_d a_j'`, with no competing term; hence it must
vanish.  Descending gives `a_i'=0` for every `i`.  Equation `(12)` then
reduces to `-2yU_v`, which cannot be a nonzero constant.  If `d=0`, a
nonconstant `U(y)` makes `UU_y` have degree `2 deg(U)-1>0`, while constant
`U` gives zero.  This handles the final pair.

Therefore none of

```text
(L,T), (L,U), (L,S), (T,U), (T,S), (U,S)               (14)
```

is a planar Keller map, for any polynomial `U`.

## 4. An everywhere-immersive positive packet

The obstruction in `(14)` is not a rank defect of the full packet.  Take

```text
U_*=1+y-y^2/2-(3/2)vy(y-3)                             (15)
```

and define `Z=(L,T,U_*,S):A^2->A^4` by `(1)`.  Let `M_ij` be its six
two-by-two Jacobian minors.  Exact Groebner reduction over `Q` gives

```text
(M_12,M_13,M_14,M_23,M_24,M_34)=(1) in Q[v,y].         (16)
```

Hence the differential of `Z` has rank two at every geometric affine point:
`Z` is everywhere immersive.  A second exact linear solve gives

```text
sum_(i<j) c_ij M_ij != 1                               (17)
```

for every six-tuple of constants `c_ij`.  Since the Jacobian of any constant-
linear projection `A^4->A^2` is such a constant combination, no
constant-linear projection of this packet is Keller.

The companion checks `(1)`--`(17)` with exact symbolic arithmetic and keeps
all truth gates active under optimized execution.  The Groebner and linear
systems are finite exact calculations, not numerical evidence.

## 5. The nonlinear projection equation

The unit ideal `(16)` says an arbitrary source-polynomial combination of the
six minors can equal one.  It does **not** yet produce a legal projection.
For polynomials `A,B in k[Z_1,Z_2,Z_3,Z_4]`, the chain rule is

```text
Jac(A(Z),B(Z))
 =sum_(i<j) [A_iB_j-A_jB_i](Z) M_ij.                  (18)
```

The six coefficient functions in `(18)` are highly constrained:

1. they must descend through the packet subalgebra `k[L,T,U_*,S]`;
2. they are Pluecker coordinates of the decomposable two-form `dA wedge dB`,
   so they obey

   ```text
   c_12 c_34-c_13 c_24+c_14 c_23=0;                  (19)
   ```

3. they satisfy the differential integrability conditions coming from an
   exact decomposable form; and
4. the resulting pair must retain a nontrivial fiber and the index/nonproper
   sidecars, not merely have formal rank two.

This converts the open nonlinear-projection search into a typed syzygy
problem: find a descending, integrable, decomposable Bezout certificate for
the minor ideal.  Constant coefficients fail by `(17)`; arbitrary source
coefficients exist by `(16)`; the gap between them is now exactly stated.

## 6. Counterexample architecture and scope

The packet combines the two boundary models from THM-3554 and THM-3555:

```text
marked factor vX+1       -> one root, escaping at v=0;
quadratic factor         -> Kummer pair with square class -L;
disc=-4LS^2              -> odd L is the prospective infinity owner;
four-coordinate packet  -> no differential rank defect.              (20)
```

A successful nonlinear projection in `(18)` would still require a full
audit; it could simply recover source coordinates and be an automorphism.
Conversely, failure in a bounded projection degree would not exclude other
packets or planar counterexamples.  What is proved is the owner
factorization, all natural projection no-gos, the explicit immersive packet,
and the exact nonlinear syzygy target.  No `JC(2)` conclusion is claimed.
