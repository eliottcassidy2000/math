---
id: THM-2985
title: "Multiparameter normal map and arc-factor separation"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED. For a square formal matrix of
  special-fibre corank r, the degree-r determinant initial form is, up to a
  unit, the determinant of the canonical kernel-to-cokernel normal map.
  Along a chosen formal arc, THM-2960's Smith bars are exactly the lifetimes
  of liftable special-fibre kernel classes. Arc contact, ambient tangent-cone
  multiplicity, and ambient divisibility are rigorously separated.
source: codex-multiparameter-normal-map-2026-07-30
depends_on:
  - THM-2960-local-smith-jet-fitting-barcode-and-negative-depth-chamber-atlas
related:
  - THM-2983-coloured-macaulay-wall-normal-cone-excess-contact-atlas
audit: >
  An independent hostile audit rederived the block-Schur initial form and
  the liftable-kernel/Smith identities, checked the tangent-arc factor
  counterexample, and required the no-common-kernel singular-pencil boundary
  now included below. No computation is used.
---

# THM-2985 -- multiparameter normal map and arc-factor separation

**PROVED + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2960 proves the one-parameter DVR Smith-barcode identity. This theorem
adds the missing multiparameter layer and separates three quantities that
coincide only in the transverse case: special-fibre corank, ambient
normal-cone degree, and contact order along a chosen arc.

## 1. The ambient first normal map

Let `k` be a field and

```text
R=k[[x_1,...,x_s]],       m=(x_1,...,x_s),
A in Mat_n(R),            A_0=A mod m,
corank(A_0)=r.                                             (1)
```

Put `K=ker(A_0)` and `C=coker(A_0)`. The linear part of `A` induces the
canonical map

```text
nu_1: K -> C tensor_k m/m^2.                              (2)
```

After choosing bases of the `r`-dimensional spaces `K,C`, regard `nu_1` as
an `r by r` matrix of linear forms. Then

```text
ord_m(det A) >= r,
in_m^r(det A)=u det(nu_1),             u in k^*.          (3)
```

The scalar `u` depends on complementary bases, but the projective tangent
cone does not. In particular,

```text
ord_m(det A)=r  iff  det(nu_1) is not the zero polynomial. (4)
```

A nonzero common right kernel of all coefficient matrices in `nu_1`, or a
nonzero common left kernel, is a sufficient certificate of excess contact.
It is not necessary. For example,

```text
nu_1 = [ 0   x   y]
       [-x   0   z]
       [-y  -z   0]                                      (5)
```

has identically zero determinant over every field, while the intersection
of the kernels of its three coefficient matrices is zero (and likewise on
the left). Thus the lawful invariant is the singular matrix space itself,
not only its common kernels.

## 2. Arcwise Smith lifetimes

Choose a local `k`-algebra map `gamma:R->k[[t]]` through the origin and
suppose `det A_gamma` is nonzero. Smith reduction over the DVR gives positive
integers `alpha_1,...,alpha_r` such that, up to invertible row and column
operations,

```text
A_gamma ~ diag(1,...,1,t^alpha_1,...,t^alpha_r),
ord_t(det A_gamma)=sum_i alpha_i.                         (6)
```

For `j>=1`, define the `k`-subspace

```text
K_j(gamma)={v_0 in K : some v(t)=v_0 mod t satisfies
                         A_gamma v(t)=0 mod t^j}.          (7)
```

Then

```text
dim_k K_j(gamma)=#{i:alpha_i>=j},
ord_t(det A_gamma)=sum_(j>=1) dim_k K_j(gamma),
ord_t(det A_gamma)-r=sum_(j>=2) dim_k K_j(gamma).         (8)
```

Thus excess arc contact is exactly the total lifetime beyond the first
layer of the liftable kernel classes. This is THM-2960's Smith barcode in
special-fibre-kernel form; `K_j/K_(j+1)` records the death grade.

## 3. Proof

Constant invertible row and column operations put `A_0` in block form

```text
A_0 = [ I_(n-r)  0 ] .                                   (9)
      [    0      0 ]
```

Write

```text
A = [ I+B  C ],             B,C,D,E in m.
    [  D    E ]                                          (10)
```

The block `I+B` is invertible. Block Gaussian elimination gives

```text
det A=det(I+B) det S,
S=E-D(I+B)^(-1)C.                                        (11)
```

The matrix `S` lies in `m`, and its linear part is exactly `(2)` because
the correction in `(11)` lies in `m^2`. Every term of `det S` has degree at
least `r`, and its degree-`r` part is the determinant of that linear part.
Since `det(I+B)` has nonzero constant term, `(3)` and `(4)` follow.

For `(8)`, reduce the Smith column matrix modulo `t`; its last `r` columns
identify the special-fibre kernel with the positive-valuation coordinate
lines. On a line with entry `t^alpha`, a nonzero residue vector lifts to a
solution modulo `t^j` exactly when `alpha>=j`. This proves the first identity.
The other two follow by summing

```text
alpha=sum_(j>=1) 1_(alpha>=j).                            (12)
```

Although Smith bases are not canonical, the subspaces `(7)` are, so their
dimensions and the death multiset are independent of the reduction.

## 4. Arc contact is not ambient divisibility

A factor or order seen after pullback to one arc cannot in general be
divided in `R`. Such division requires an ambient ideal or divisorial claim,
for example

```text
det A in (F^e) in R                                      (13)
```

for a specified prime divisor `F`. Take

```text
A(x,y)=diag(1,x+y^2),       F=y,
gamma(x)=0,                 gamma(y)=t.                  (14)
```

Then

```text
ord_t gamma(det A)=2=2 ord_t gamma(F),
det A=x+y^2 notin (F^2).                                 (15)
```

The extra order comes from tangency to the ambient tangent cone `x=0`, not
from an ambient square factor. Conversely, an arc contained in `F=0` pulls
`F^e` back to zero and has infinite rather than order-`e` contact. Hence one
arc order is neither necessary nor sufficient for `(13)`.

## 5. Scope and use

If `det(nu_1)=0`, the next lawful object is a Schur-complement/kernel flag or
an iterated normal cone after a direction is chosen. There need not be a
canonical higher ambient map before that choice. A Gram or Hessian zero mode
is weaker still: an isotropic pairing can be singular while the exact
presentation is injective, as in the repo's `S3`-equivariant SNC geometric
hostile.

For the width-six coloured Macaulay application, `(3)` explains the exact
half-wall and root-one tangent pencils, `(8)` is the correct language for
higher-contact integer walls, and `(14)` explains why a wall order on the
diagonal cannot simply be divided from the coloured determinant. Those
numerical applications belong to THM-2983 and are not proved here.

QED.
