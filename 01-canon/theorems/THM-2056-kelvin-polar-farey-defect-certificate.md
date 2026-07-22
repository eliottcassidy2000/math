---
id: THM-2056
title: "Kelvin-polar linearization and the Farey defect certificate for rank-two LRC gates"
status: >
  PROVED. Kelvin inversion sends the quadratic-versus-linear THM-2053 gate to
  membership in one fixed rational polar polygon. On every acute unimodular
  owner cone, one exact defect inequality certifies every interior lattice
  point, even when its boundary generators are themselves uncertified. The
  resulting finite Farey fan is a checkable determinant-gate certificate, not
  LRC(14) by itself.
source: codex-2026-07-21-LRC-kelvin-farey
script: 04-computation/lrc_kelvin_farey_scaled_core_codex_20260721.py
result: 05-knowledge/results/lrc_kelvin_farey_scaled_core_codex_20260721.out
script_sha256: 9f4d2780c6d9d24dfcdf8e99ae1a1fe327ad98e97e21044effd95769c0e4c309
result_sha256: 35c25efa343bf489f79c9593379378699483c254915f9d9c50a080d119c9b520
hash_basis: normalized repository blobs (LF)
depends_on:
  - THM-2053
  - THM-2055
related:
  - HYP-8871
  - MISTAKE-225
---

# THM-2056 -- Kelvin-polar and Farey-defect certificate

Keep the notation of THM-2055. Thus

```text
c_i=(u_i,z_i),                 d=(a,b),
K=conv{+c_i,-c_i},             R(a,b)=(-b,a),
D(d)=max_i |a z_i-b u_i|=h_K(Rd).                         (1)
```

Assume that the columns span `R^2`, and let

```text
B_D={x:D(x)<=1}=R^(-1)K^o.                                (2)
```

The THM-2053 sufficient gate is

```text
D(d)<=||d||^2/91.                                         (3)
```

## 1. Kelvin inversion makes the gate linear

For `d!=0`, define the Kelvin inversion

```text
I(d)=d/||d||^2.                                           (4)
```

**Theorem.** Gate (3) is equivalent to

```text
I(d) in P_K:=(1/91)B_D=(1/91)R^(-1)K^o.                  (5)
```

*Proof.* The support norm is positively homogeneous, so

```text
D(I(d))=D(d)/||d||^2.
```

Thus (3) is equivalent to `D(I(d))<=1/91`, which is exactly
(5). QED.

The tangent circles of THM-2055 have therefore not introduced a nonlinear
arithmetic object. They are the inverse images of the facets of one centrally
symmetric rational polygon. Explicitly, a signed hull vertex `p` gives the
facet

```text
p dot R x <= 1/91.                                        (6)
```

Its strict violation by `x=I(d)` is

```text
||d||^2 < 91 p dot R d,                                   (7)
```

which completes the square to THM-2055's tangent disk.

For primitive `d=(a,b)`, put `q=a^2+b^2`. Then

```text
I(d)=(a/q,b/q),
a/q+i b/q=1/(a-i b).                                      (8)
```

Both fractions have exact common denominator `q`, because
`gcd(a,b)=1`. Hence the residual is exactly the set of reciprocal visible
Gaussian lattice points lying outside `P_K`. This is a genuine arithmetic
address, but it is Gaussian reciprocity against a polygon, not a Heegner form
or a class group.

## 2. The Farey defect identity

Fix a signed hull owner `p=(r,s)` and put `w_p=(s,-r)`. On its owner cone
`C_p`, THM-2055 gives `D(d)=w_p dot d`. Define

```text
F_p(d)=||d||^2-91 w_p dot d,
A_p(d)=max(0,-F_p(d)).                                    (9)
```

Thus `F_p(d)>=0` is precisely the determinant certificate, and `A_p(d)` is
the nonnegative gate defect.

Let `u,v` be primitive generators of an acute unimodular cone contained in
`C_p`:

```text
|det(u,v)|=1,             u dot v>=0.                     (10)
```

Every lattice point in its relative interior is uniquely `m u+n v` with
positive integers `m,n`. Direct expansion gives the exact identity

```text
F_p(m u+n v)
 =m F_p(u)+n F_p(v)
  +(m^2-m)||u||^2+(n^2-n)||v||^2+2mn u dot v.             (11)
```

**Defect certificate.** If

```text
2 u dot v >= A_p(u)+A_p(v),                               (12)
```

then every lattice point in the relative interior of the cone is certified
by (3). Only the two primitive boundary generators can remain unresolved.

*Proof.* Write `A=A_p(u)` and `B=A_p(v)`. The two square terms in (11) are
nonnegative, while `F_p(u)>=-A` and `F_p(v)>=-B`. Therefore

```text
F_p(m u+n v)>=-mA-nB+2mn u dot v
             >=-mA-nB+mn(A+B)>=0,                        (13)
```

because `m,n>=1`. QED.

This contains the safe-endpoint rule as the case `A=B=0`. More usefully, it
can isolate a bad ray. For example, with `w_p=(1,0)`,

```text
u=(1,0):   F_p(u)=-90,
v=(91,1):  F_p(v)=1,
2 u dot v=182>=90.
```

Every interior lattice point of `cone(u,v)` is certified although `u` is not.

The acuteness hypothesis is load-bearing. With the same `w_p`, take
`u=(91,1)` and `v=(-90,-1)`. The pair is unimodular and both endpoints are
safe, but `u dot v<0` and `u+v=(1,0)` is not safe.

## 3. Finite regular-fan certificate

The uncertified primitive set in an owner cone is finite, because

```text
F_p(d)<0  implies  ||d||<91||w_p||.                       (14)
```

There is a finite acute unimodular fan in which every uncertified primitive
direction is a ray and every two-dimensional cone satisfies (12).

To see this, first insert the finitely many bad primitive rays and split every
angular sector to width at most `pi/2`. Around a bad primitive generator `u`,
choose a Farey neighbor `v_0` with `|det(u,v_0)|=1`. The neighbors
`v_k=v_0+k u` approach the ray of `u`, remain unimodular with it, and satisfy

```text
u dot v_k -> +infinity,             F_p(v_k)->+infinity.
```

Thus a sufficiently large `k` makes `v_k` safe and makes (12) hold. Do this
on both sides of every bad ray. The remaining rational gaps contain no bad
primitive directions; regularize them by the ordinary Euclidean/Farey
subdivision. Their ray generators are safe, and acuteness makes (12)
automatic. This proves existence of the finite certificate.

The certificate does not magically shrink the set of failures: a bad
primitive direction must still occur as a labelled ray. Its value is that the
entire infinite complement is discharged by finitely many determinant and
dot-product checks, with an exact slot on which to attach a pair-sum or Euler
certificate.

## 4. One-tail control

For `u=(1,...,13)` and `z=e_12`, THM-2055 gives

```text
D(a,b)=max(13|b|,|a-12b|).
```

The polar polygon in (5) is the parallelogram

```text
13|y|<=1/91,              |x-12y|<=1/91.                 (15)
```

For `d=(1,b)`, the reciprocal point is

```text
(x,y)=(1/(1+b^2),b/(1+b^2)).                             (16)
```

The first facet in (15) is exactly

```text
b^2-1183|b|+1>=0,
```

so every integral `|b|>=1183` is safe. The points `(1,0)` and `(1,2)` lie
outside (15), yet the corresponding LRC rows have exact maxima `1/14` and
`1/12`. This is the mandatory scope control: outside the polar polygon means
only that THM-2053 did not certify the row.

## 5. Predicate and loss ledger

Kelvin inversion and the defect fan preserve:

- the exact determinant-gate value and safe/fail predicate;
- the signed hull owner, primitive parameter, and rational slope;
- a finite proof that every nonlisted primitive direction is gate-safe.

They destroy the non-hull runner labels, phase height, pair-sum clock, and
boundary owner word. Thus the natural proof vertices are reciprocal lattice
points and Farey cones, not runners or quadratic-form classes. THM-2057 gives
the first successful sidecar: a modular clock/binding argument closes the
entire one-tail plane that (15) leaves with a large residual.
