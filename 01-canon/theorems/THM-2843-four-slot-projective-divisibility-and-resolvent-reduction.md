---
id: THM-2843
title: "Four-slot projective divisibility and resolvent reduction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every
  four-slot consecutive moment window either already has a
  positive-dimensional common-zero counterexample in its first three
  forms, or reduces to the determinant of multiplication by the fourth
  form on a finite complete intersection of length
  d(d+1)(d+2).  The determinant is a nonnegative real norm and is positive
  exactly when that cell is closed.  In the first window, common nullity is
  equivalent to a real two-plane on which the quadratic divides both the
  cubic and quartic.  In the reduced case, the six points form three
  conjugate pairs, but an exact hostile shows that this 2-by-3 count does
  not carry the standard holomorphic C2*C3 modular action.
source: root/four-slot-projective-resolvent-reduction-2026-07-28
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
related:
  - THM-2812-consecutive-three-slot-factorial-moment-six-detection
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2846-arbitrary-positive-cone-moment-three-transverse-boundary
script: 04-computation/gmc_four_slot_projective_resolvent_thm2843.py
output: 05-knowledge/results/gmc_four_slot_projective_resolvent_thm2843.out
script_sha256: 4832a9e4cda1608473e3a4bcfcb880a7fa3f1c6db47b7e939b4b5183cb9549aa
output_sha256: 919b9492b2036f5f558372c0240a212e054baf3ca2228bcfe53766151b9a002d
hash_basis: LF-normalized bytes
---

# THM-2843 -- four-slot projective divisibility and resolvent reduction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Fix

```text
0<=a_0<a_1<a_2<a_3,
e_j(s)=s^(a_j)/a_j!,
H_X(s)=sum_(j=0)^3 X_j e_j(s),                         (1)
```

and define the real homogeneous moment forms

```text
M_m(X)=L(H_X^m)
 =sum_(i_0+...+i_3=m)
   multinomial(m;i_0,...,i_3)
   (sum_j i_j a_j)!/product_j(a_j!)^(i_j)
   product_j X_j^(i_j).                                (2)
```

For a consecutive window beginning at `d=k+1>=1`, put

```text
G_0=M_d,       G_1=M_(d+1),
G_2=M_(d+2),   G_3=M_(d+3).                            (3)
```

This theorem gives an exact finite-algebra reduction of four-slot
nullity.  It does not prove that the final determinant is always nonzero.

## 1. Complete-intersection dichotomy

Let

```text
Z=V_(P^3_C)(G_0,G_1,G_2).                              (4)
```

Exactly one of the following holds.

1. `Z` has positive dimension.  Then

   ```text
   Z intersect V(G_3) != empty,                        (5)
   ```

   so the support/window already gives a common projective zero of all
   four moment forms.

2. `Z` is zero-dimensional.  Then `G_0,G_1,G_2` form a homogeneous
   regular sequence and

   ```text
   length(Z)=D=d(d+1)(d+2).                            (6)
   ```

For `(5)`, a positive-degree hypersurface meets every
positive-dimensional projective component unless it contains that
component, in which case the conclusion is immediate.  In the second
branch, the ideal has height three in the Cohen--Macaulay ring
`C[X_0,...,X_3]`; three generators of height three form a regular
sequence, and Bezout gives `(6)`.

Choose a real linear form `ell` that misses the finite support of `Z`,
and dehomogenize:

```text
A_ell=
 R[X_0,...,X_3]/(G_0,G_1,G_2,ell-1),       dim_R A_ell=D. (7)
```

Let `g` be the image of `G_3`.  Then

```text
N_ell=det_R(m_g:A_ell -> A_ell)                        (8)
```

satisfies

```text
N_ell!=0
 iff Z intersect V(G_3)=empty.                         (9)
```

Indeed, multiplication by `g` is invertible in the Artinian algebra
exactly when `g` survives in every residue field.

## 2. The determinant is a nonnegative conjugate norm

One of the three consecutive degrees `d,d+1,d+2` is even.  For every
nonzero real `X`,

```text
M_(2r)(X)
 =integral_0^infinity H_X(s)^(2r)e^(-s) ds
 >0.                                                   (10)
```

Consequently `Z(R)=empty`.  After complexification, the local factors of
`A_ell` occur in conjugate pairs.  Hence

```text
N_ell
 =product_({P,Pbar})
   |det_C(m_g on (A_ell tensor C)_P)|^2
 >=0.                                                  (11)
```

Combining `(9)--(11)`,

```text
N_ell>0
 iff the four-slot support/window has no common zero.  (12)
```

The first three window dimensions are

```text
k=0: D=6,             k=1: D=24,             k=2: D=60. (13)
```

Thus four-slot SFC reduces to strict positivity of explicit real norms.
The reduction does not prove those norms nonzero.

## 3. First-window moving-plane divisibility

Now take `k=0`.  Work in the three-dimensional real mean-zero space

```text
W_0={H:L(H)=0}.                                        (14)
```

With

```text
U_1=e_1-e_0,       U_2=e_2-e_1,       U_3=e_3-e_2,
H=xU_1+yU_2+zU_3,                                    (15)
```

put

```text
Q=L(H^2),             C=L(H^3),             F=L(H^4). (16)
```

The ternary quadratic `Q` is positive-definite over `R`.  It does not
divide `C`: otherwise restricting to `z=0` would make the corresponding
three-slot binary quadratic divide its cubic, contradicting THM-2824.
Thus `Q,C` form a complete intersection of length six in `P^2_C`.

There is a common projective zero of `Q,C,F` if and only if there is a
real two-plane `E subset W_0` such that

```text
Q|_E divides C|_E,             Q|_E divides F|_E.      (17)
```

For the forward direction, write a common point as `P=u+iv`.  Positive
definiteness makes `u,v` independent.  On
`E=span_R(u,v)`, the positive binary quadratic `Q|_E` has precisely the
two projective roots `[P],[Pbar]`; the real forms `C,F` vanish at both,
so the quadratic divides each restriction.  Conversely, either complex
root of a divisibility plane gives a common zero.

Equivalently, if `E=ker(lambda)`, there are real forms of the displayed
degrees satisfying

```text
C=Q B_1+lambda A_2,
F=Q B_2+lambda A_3.                                   (18)
```

When `Q intersect C` is reduced, its six points form three conjugate
pairs and hence three unordered real chord planes.  The norm `(11)` is
the product of the three squared fourth-moment obstructions.  There is no
canonical cyclic ordering of the three planes.

## 4. Binary sextic/octet resultant

A conic parametrization

```text
nu:P^1_C -> V(Q)
```

turns

```text
c_6=(C o nu)/2,              f_8=(F o nu)/3            (19)
```

into normalized binary forms of degrees six and eight.  The nonzero
scalar divisions do not change their zero divisors.  Then

```text
Q=C=F=0 has a projective solution
 iff Res(c_6,f_8)=0.                                   (20)
```

For the consecutive support `{0,1,2,3}`,

```text
Q=x^2+2xy+2xz+2y^2+6yz+6z^2
 =(x+y+z)^2+(y+2z)^2+z^2.                             (21)
```

Put `X=x+y+z`, `Y=y+2z`, `Z=z`, and parametrize

```text
X=u^2-v^2,            Y=i(u^2+v^2),            Z=2uv.
```

Exact substitution gives

```text
c_6=
 (-5-2i)u^6+(-36+18i)u^5v+(27+126i)u^4v^2
 +152u^3v^3+(-27+126i)u^2v^4
 +(-36-18i)uv^5+(5-2i)v^6,                            (22)

f_8=
 (44-96i)u^8+(-992-1264i)u^7v
 +(-10368+4160i)u^6v^2+(6624+39216i)u^5v^3
 +51048u^4v^4+(-6624+39216i)u^3v^5
 +(-10368-4160i)u^2v^6+(992-1264i)uv^7
 +(44+96i)v^8.                                        (23)
```

Their resultant is

```text
208741470184115575361509867388928
 =2^36 3^7 67 11702701 1771410437
 !=0.                                                  (24)
```

The leading coefficient of `c_6` is nonzero, so this affine resultant
does not lose the point `[1:0]`.  The sextic is squarefree, with

```text
Disc(c_6)
 =6291863054003994624
 =2^20 3^12 31^3 379.                                 (24a)
```

There is also an independent quotient certificate.  The face `z=0` has
resultant `116y^6`, so the complete intersection lies in the chart
`z=1`.  Direct multiplication by `F` in its six-dimensional real
quotient has determinant

```text
9070189700378194715889733632/707281>0.                 (24b)
```

Thus the first-window cell on `{0,1,2,3}` is closed exactly, by both the
binary resultant and the Artin-algebra norm.

## 5. The `C_2*C_3` modular action is not intrinsic

On the conic coordinate, take the standard projective generators

```text
sigma[u:v]=[-v:u],             tau[u:v]=[-v:u+v],
sigma^2=tau^3=1.                                      (25)
```

The exact sextic discriminant is

```text
Disc(c_6)=6291863054003994624!=0,
```

so its residual divisor is reduced.  Matching the `u^6` coefficient of a
possible scalar invariance forces

```text
lambda=(-21+20i)/29.
```

But exact coefficient extraction gives

```text
[u^5v](c_6 o sigma-lambda c_6)=(648+1620i)/29!=0,
[u^5v](c_6 o tau-lambda c_6)  =(1518+1272i)/29!=0.    (26)
```

Neither free factor preserves the conic--cubic residual divisor.
The genuine pairing involution is

```text
j[u:v]=[-conj(v):conj(u)],                             (27)
```

which is fixed-point-free and anti-holomorphic.  The honest `2 x 3`
structure is therefore conjugation within each of three unordered pairs.
It is not a canonical holomorphic `C_2*C_3` action.  Producing a cyclic
`C_3` would require an extra orientation or discriminant-square sidecar.

As a finite hostile control, reduction with `i=4` in `F_17` gives a
trivial stabilizer of `c_6` inside all `4896` elements of
`PGL_2(F_17)`.  This is not promoted to a characteristic-zero
automorphism theorem.  It only reinforces the typed conclusion: the
pairing involution is intrinsic, while an order-three action needs an
additional sidecar.

## 6. Facewise three-slot control is insufficient

The abstract real forms

```text
Q=x^2+y^2+z^2,
C=2x^3+2y^3+(x+y)z^2,
F=(x-y)x^3                                           (28)
```

have the common point

```text
[x:y:z]=[1:1:i sqrt(2)],                              (29)
```

while `Q,C` have no common point on any coordinate plane.  This is not a
factorial counterexample.  It is the sharp logical hostile to a proof
that checks only three-slot faces: the missing invariant is the moving
real divisibility plane `(17)`.

The theorem therefore supplies a finite algebra, a positive norm, and a
precise moving-plane target.  It does not close general four-slot SFC,
construct a modular action, or turn the `2 x 3` count into symmetry.

## 7. Exact companion

The exact companion checks:

1. the general four-slot moment formula on `24` deterministic cells and
   positive-definiteness of four first-window Gram forms;
2. the explicit quadratic, sextic, octic, resultant, factorization, and
   squarefree sextic on `{0,1,2,3}`;
3. the six-dimensional quotient multiplication norm by both reduction
   and an independent direct matrix determinant;
4. exact moving-plane ideal identities;
5. conjugate, nonreduced, and real-point controls for the norm sign;
6. the hidden-plane facewise hostile;
7. the exact `sigma` and `tau` coefficient defects; and
8. a finite `PGL_2(F_17)` control whose sextic stabilizer is trivial.

It explicitly does not claim the full `378`-cell atlas.  All truth-bearing
gates are explicit exceptions and all arithmetic is exact.  Reproduce
with

```text
python 04-computation/gmc_four_slot_projective_resolvent_thm2843.py
python -O 04-computation/gmc_four_slot_projective_resolvent_thm2843.py
```

Both modes byte-match the stored transcript.

## 8. Independent hostile audit

An independent audit rederived the complete-intersection dichotomy, the
real dehomogenization and nonreduced conjugate-local norm, the
moving-plane iff and global ideal lift, and the explicit modular
nonaction.  It caught and repaired the only exact normalization defect:
the displayed sextic and octic are `(C o nu)/2` and `(F o nu)/3`, not the
raw pullbacks.  It also enforced the reduced-case qualifier on the
three-pair language and replayed normal and optimized companions against
the stored transcript and declared LF hashes.

**QED.**
