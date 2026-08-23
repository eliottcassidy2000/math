---
id: THM-3837
title: "Affine corrections cannot realize the first comaximal line-hyperbola bichromatic fibre"
status: >
  PROOF CANDIDATE + VERIFIED-EXACT, PENDING INDEPENDENT HOSTILE AUDIT.  Fix a
  cubic spectral root alpha.  In the n=1 disconnected ansatz u=x,
  v=xy-1, f=uv, k=x-v+fq, h=alpha k+f, and
  d=-b(alpha)x^2/alpha^2+fD, no affine profiles q,D satisfy the THM-3827
  selector equation.  First normal contact leaves one scalar parameter;
  second contact on the line and hyperbola forces two incompatible values,
  with exact resultant -5817545.  Quadratic profiles remain OPEN over
  characteristic zero; a stated mod-11 two-contact scout is FINITE-EXACT
  only.  No second row or Jacobian counterexample is constructed.
source: jc_sparse_direct_search / THM-3830 disconnected-fibre successor, 2026-08-23
audit: >
  Self-audited after correcting the hyperbola derivative to
  k_f=q-1/x.  The deterministic companion checks comaximality, an explicit
  determinant-one completion, the CRT sign boundary, both first-contact
  expansions, the complete five-equation first-jet normal form, the two
  second-contact linear forms and all unit resultants.  It separately checks
  the quadratic mod-11 scout without promoting it to characteristic zero.
  Normal and optimized replay agree with the frozen transcript.  Independent
  hostile audit is pending.
depends_on:
  - THM-3830-coordinate-cross-bichromatic-split-nonentry
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
related:
  - THM-3831-intrinsic-spectral-pencil-fibre-atlas-and-forced-cubic-two-arm-hit
  - THM-3832-nonlinear-cubic-root-ratio-triangular-birational-chart
script: 04-computation/jc2_comaximal_line_hyperbola_affine_contact_thm3837.py
output: 05-knowledge/results/jc2_comaximal_line_hyperbola_affine_contact_thm3837.out
script_sha256: 002eeaecbed2133cf2b593dcff0173d95766e09e24e83f4dc50e739404b88ea2
output_sha256: 6084bdd2d0e74d5b62499c7964118ddfb22fb508fb09ccf5df0c6624dc75e42f
semantic_sha256: 3f345677c44bb03fa61e660f7e3db26af85cb5407fee65e1213c099bb29b268c
hash_basis: raw LF bytes
---

# THM-3837 -- the first disconnected fibre fails at second contact

**PROOF CANDIDATE + VERIFIED-EXACT, PENDING INDEPENDENT HOSTILE AUDIT.**
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

For affine profiles

```text
q=q_0+q_x x+q_y y,             D=d_0+d_x x+d_y y          (3)
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
exactly `n=1` and affine `q,D`.  It does not close quadratic or arbitrary
profiles.

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
q_x=0,                      q_y=1+alpha,
d_y=0,                      d_0=-rho,

q_0=(-alpha^2 d_x-b'(alpha)+2B/alpha)/(2B).                (20)
```

Thus first contact is consistent and leaves `d_x` free.  In particular,
the obstruction is not the CRT boundary itself and is not a spurious
leading-degree failure.

## 3. The incompatible second contacts

Substitute `(20)` before taking the next coefficients.  The constant-in-`Y`
part of the `epsilon^2` coefficient in the line chart `(16)` is

```text
-Q(alpha)L_u,
L_u=d_x+18alpha+38+7alpha^2.                              (21)
```

In the hyperbola chart `(18)`, the coefficient of `X^5` in the
`epsilon^2` coefficient is

```text
-(1/3)U(alpha)L_v,
U(z)=49z^2-9z+7,
L_v=d_x+325alpha-55+147alpha^2.                           (22)
```

The second multiplier is also a unit:

```text
Res(r,U)=14535.                                           (23)
```

Consequently `(6)` would require `L_u=L_v=0`.  Their difference is

```text
L_v-L_u=140alpha^2+307alpha-93.                           (24)
```

But

```text
Res_z(r(z),140z^2+307z-93)=-5817545!=0.                  (25)
```

This contradiction proves the theorem.

## 4. Quadratic boundary and exact scope

At the good specialization `K=F_11`, `alpha=1`, an exact coefficient-ideal
calculation for all total-degree-at-most-two `q,D` gives a unit ideal after
the first two contacts.  This is a **FINITE-EXACT scout only**: a unit ideal
after specialization does not by itself imply that the characteristic-zero
ideal is a unit.  Therefore quadratic `q,D` remain **OPEN** here.

The theorem also leaves open other boundary units `k|_(v=0)=x^n`, other
comaximal component pairs, the remaining intrinsic `D` and second-row laws,
and the planar Keller equation.  No Jacobian counterexample is claimed.
**QED, pending independent hostile audit.**
