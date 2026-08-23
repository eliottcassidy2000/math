---
id: THM-3837
title: "Quadratic corrections cannot realize the first comaximal line-hyperbola bichromatic fibre"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Fix a
  cubic spectral root alpha.  In the n=1 disconnected ansatz u=x,
  v=xy-1, f=uv, k=x-v+fq, h=alpha k+f, and
  d=-b(alpha)x^2/alpha^2+fD, no total-degree-at-most-two profiles q,D
  satisfy the THM-3827 selector equation.  First normal contact leaves four
  scalar parameters.  Second contact kills the two genuinely quadratic D
  directions, then the line and hyperbola force two incompatible values of
  d_x, with exact resultant -5817545.  No second row or Jacobian
  counterexample is constructed.
source: jc_sparse_direct_search / THM-3830 disconnected-fibre successor, 2026-08-23
audit: >
  INDEPENDENTLY HOSTILE-AUDITED by root on 2026-08-23 after the earlier
  affine packet was independently audited by jc_quartic_c3_construct.  The
  audit rederived the full Laurent first-contact comparison, verified that
  its six equations leave exactly d_x,d_xx,d_xy,q_xy, checked the universal
  second Taylor coefficient against both normal charts, and confirmed that
  the line Y and hyperbola X^8 coefficients kill d_xy and d_xx before the
  incompatible unit-multiplied L_u/L_v equations.  It also checked that the
  resultant argument works over every characteristic-zero field containing
  any root alpha and that the statement remains confined to the displayed
  selector cell.  The deterministic companion, corrected at the hyperbola
  derivative k_f=q-1/x, checks comaximality, an explicit determinant-one
  completion, the CRT sign boundary, the complete quadratic first-contact
  normal form, the exact universal second Taylor coefficient, the forced
  d_xy=d_xx=0 peel, the two terminal linear forms, and all unit resultants.
  Normal and optimized replay agree with the frozen transcript.
depends_on:
  - THM-3830-coordinate-cross-bichromatic-split-nonentry
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
related:
  - THM-3831-intrinsic-spectral-pencil-fibre-atlas-and-forced-cubic-two-arm-hit
  - THM-3832-nonlinear-cubic-root-ratio-triangular-birational-chart
  - THM-3836-cubic-factor-cofactor-darboux-packet
script: 04-computation/jc2_comaximal_line_hyperbola_affine_contact_thm3837.py
output: 05-knowledge/results/jc2_comaximal_line_hyperbola_affine_contact_thm3837.out
script_sha256: 74b8dd003383be9737bccf263ec76a276d8713bf07f76260cb5f4c336652c560
output_sha256: 4c446eb2f9c43bb787c0d520fe3a2bc82081f3088be86340fe44f5de6108861f
semantic_sha256: 6ece9fdb67e87d101d673126f289824b5131a5da993cdd3b4569a178c87bec77
hash_basis: raw LF bytes
---

# THM-3837 -- the first disconnected fibre fails through quadratic contact

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over a field `K` of characteristic zero containing a root `alpha` of

```text
r(z)=3z^3+7z^2+1.                                          (1)
```

Put

```text
Q(z)=7z^2+3,
a(z)=Q(z)r(z),
b(z)=(z+1)(2z+1)(3z-1),
B=b(alpha),                 rho=r'(alpha).                 (2)
```

For profiles of total degree at most two

```text
q=q_0+q_x x+q_y y+q_xx x^2+q_xy xy+q_yy y^2,
D=d_0+d_x x+d_y y+d_xx x^2+d_xy xy+d_yy y^2.             (3)
```

define

```text
u=x,             v=xy-1,             f=uv,
k=x-v+fq,        h=alpha k+f,
mu=-B/alpha^2,   d=mu x^2+fD.                              (4)
```

Use the THM-3827 homogeneous factors

```text
A_5=(7h^2+3k^2)(3h^3+7h^2k+k^3),
B_3=(h+k)(2h+k)(3h-k).                                    (5)
```

Then the selector equation

```text
A_5=d(kB_3+h^2d)                                          (6)
```

has no solution with profiles `(3)`.

This is an all-coefficient statement over `K`, but its support scope is
exactly `n=1` and total degree at most two.  It does not close cubic or
arbitrary profiles.

## 1. The first CRT-compatible disconnected fibre

The factors `u,v` are comaximal:

```text
(xy-1)-yx=-1.                                             (7)
```

Thus `f=0` is the disjoint union of the line `u=0` and the hyperbola
`v=0`.  The row in `(4)` has the prescribed boundary units

```text
k|_(u=0)=1,                   k|_(v=0)=x.                  (8)
```

It is genuinely unimodular.  Indeed, set

```text
E_0=1+xy(y-1),
m=E_0q-(y-1)^2,                C=E_0+alpha m.              (9)
```

The identity

```text
E_0(x-v)+f(y-1)^2=1                                           (10)
```

gives

```text
Ck-mh=E_0k-mf=1.                                           (11)
```

The selector also has exactly the intended CRT data.  On `u=0`, `d=0`.
On `v=0`, one has `k=x`, `h=alpha x`, and

```text
d=mu x^2=-B k^2/alpha^2,
kB_3+h^2d=x^4(B+alpha^2mu)=0.                              (12)
```

Hence `(6)` holds identically on both components before any interior
coefficient is solved.  This is the minimal escape from THM-3830's
intersecting coordinate cross.

The cubic identities needed below are

```text
B=-Q(alpha)=alpha rho,
mu=-9alpha-14,                                             (13)
Res(r,Q)=1615,                    disc(r)=-1615.            (14)
```

In particular, every displayed division is by a nonzero element.

## 2. First normal contact

Let

```text
E=A_5-d(kB_3+h^2d).                                       (15)
```

There are two independent formal charts along `f=0`.

On the selected line put

```text
x=epsilon,                 y=Y.                            (16)
```

The coefficient of `epsilon` in `E` is

```text
-Q(alpha)(d_yY+d_0+rho).                                  (17)
```

On the hyperbola put

```text
v=epsilon,             x=X,             y=(1+epsilon)/X.  (18)
```

Here `f=Xepsilon`, and the easy-to-miss derivative is

```text
partial k/partial f at epsilon=0 = q|_(v=0)-1/X,           (19)
```

not merely `q|_(v=0)`.  Direct expansion of the coefficient of `epsilon`
in `(15)` is equivalently the Laurent identity

```text
D|_(v=0)=rho+2B/alpha^2
 -(X/alpha^2)(2B q|_(v=0)+b'(alpha)-2B/alpha).             (19a)
```

Comparison in `K[X,X^-1]`, combined with `(17)`, forces exactly

```text
q_xx=q_yy=0,                   d_y=d_yy=0,
d_0=-rho,

d_xx=-(2B/alpha^2)q_x,
d_xy=2rho+(2B/alpha^2)(1-q_y),

q_0+q_xy=(-alpha^2 d_x-b'(alpha)+2B/alpha)/(2B).           (20)
```

Thus first contact is consistent and leaves the four parameters
`d_x,d_xx,d_xy,q_xy` free.  In particular, the obstruction is not the CRT
boundary itself and is not a spurious leading-degree failure.  Setting the
quadratic coefficients to zero recovers the earlier affine relations
`q_x=0,q_y=1+alpha,d_y=0,d_0=-rho`.

## 3. The incompatible second contacts

Substitute `(20)` before taking the next coefficients.  The coefficient of
`Y` in the `epsilon^2` term of the line chart `(16)` is

```text
-Q(alpha)d_xy,                                             (21)
```

so `d_xy=0`.  The coefficient of `X^8` in the `epsilon^2` term of the
hyperbola chart `(18)` is

```text
(alpha^2/4)d_xx^2,                                        (22)
```

so `d_xx=0`.  These conclusions use `Q(alpha)alpha!=0` and no choice of a
square root.  The remaining constant-in-`Y` coefficient in the line chart
is

```text
-Q(alpha)L_u,
L_u=d_x+18alpha+38+7alpha^2.                              (23)
```

In the hyperbola chart `(18)`, the coefficient of `X^5` in the
`epsilon^2` coefficient is

```text
-(1/3)U(alpha)L_v,
U(z)=49z^2-9z+7,
L_v=d_x+325alpha-55+147alpha^2.                           (24)
```

The second multiplier is also a unit:

```text
Res(r,U)=14535.                                           (25)
```

Consequently `(6)` would require `L_u=L_v=0`.  Their difference is

```text
L_v-L_u=140alpha^2+307alpha-93.                           (26)
```

But

```text
Res_z(r(z),140z^2+307z-93)=-5817545!=0.                  (27)
```

This contradiction proves the theorem.

## 4. Exact scope

The quadratic coefficient `q_xy` occurs in another hyperbola second-contact
coefficient, but it cannot change `(21)--(27)`.  Thus no saturation,
generic-nonzero assumption, or finite-field inference is hidden in the
closure.

The theorem leaves open cubic and higher profiles, other boundary units
`k|_(v=0)=x^n`, other comaximal component pairs, the remaining intrinsic `D`
and second-row laws, and the planar Keller equation.  No Jacobian
counterexample is claimed.  **QED.**
