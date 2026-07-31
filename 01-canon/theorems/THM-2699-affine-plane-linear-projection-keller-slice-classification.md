---
id: THM-2699
title: "Affine-plane and linear-projection Keller slice classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every affine
  source plane and affine-linear rank-two target
  projection of the THM-2696 S3<S4 reflection quotient whose induced planar
  Jacobian is a nonzero constant belongs to one of three explicit projective
  normal forms.  Each survivor is a triangular polynomial automorphism of
  A2, with a displayed polynomial inverse.  Thus this entire linear-slice
  family contains no planar Jacobian counterexample.  Nonlinear source
  surfaces, nonlinear target projections, arbitrary S4 resolvents, JC(2),
  and DC(2) remain open.
source: root/reflection-quotient-affine-slice-2026-07-28
depends_on:
  - THM-2696-reflection-completed-s4-relative-different-and-coordinate-invariant-jacobian-gate
related:
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2690-normal-crossing-cyclic-cubic-resolvent-exclusion-and-reflection-completion-boundary
script: 04-computation/jacobian_s4_affine_plane_linear_projection_thm2699.py
output: 05-knowledge/results/jacobian_s4_affine_plane_linear_projection_thm2699.out
script_sha256: d53a81c30700198ed103235dd700e94fc33979d8d87ef7401853b416312301ab
output_sha256: 709f54877515e08a9a31d97951ed7ce5381c94388b287588e0cb4b2f2bddcad2
hash_basis: LF-normalized bytes
---

# THM-2699 -- every affine-linear Keller slice is triangular

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2696 proves that the three-dimensional reflection quotient itself has
the intrinsic nonconstant different `s_1s_2-s_3`.  Its nonlinear level
`s_1s_2-s_3=c` removes that factor on the source but normalizes a singular
target surface.  The next simplest escape is now completely classified:

```text
affine source plane + affine-linear target projection + constant Jacobian.
```

There are constant-Jacobian slices, but none is exotic.  Every one is an
explicit triangular automorphism of the plane.

## 1. Coordinate-free planar Jacobian

Write the THM-2696 quotient as

```text
Phi(x,y,z)=(A,B,d)=(x^2-2y, y^2-2xz, z).                 (1)
```

Let

```text
Pi_(n,delta)={n_1 x+n_2 y+n_3 z=delta},                 (2)
```

where `n!=0`, and let an affine-linear target projection have independent
linear rows `r_1,r_2`.  Its translation is irrelevant to the Jacobian.  Put

```text
ell=r_1 cross r_2 !=0.                                  (3)
```

The pair `(n,delta)` is projective under simultaneous rescaling, while `ell`
is projective by itself.  Choose oriented affine coordinates
on `(2)` whose tangent vectors `p,q` satisfy `p cross q=n`; a different
choice merely multiplies the planar Jacobian by a nonzero constant.  The
cross-product/adjugate identity gives

```text
Jac(Pi -> A2)
 =ell dot (D Phi p cross D Phi q)
 =n^T adj(D Phi) ell.                                   (4)
```

Here

```text
adj(D Phi)=
 [2y  2   4x]
 [2z  2x  4x^2]
 [0   0   4xy-4z],                                      (5)
```

so `(4)` is the polynomial

```text
K=4 ell_3 n_2 x^2+4 ell_3 n_3 xy
 +(2 ell_2 n_2+4 ell_3 n_1)x+2 ell_1 n_1 y
 +(2 ell_1 n_2-4 ell_3 n_3)z+2 ell_2 n_1.              (6)
```

Thus the classification is a divisibility problem: `K-c` must lie in the
linear ideal `(n_1x+n_2y+n_3z-delta)` for some `c!=0`.

## 2. Exhaustive coefficient classification

First suppose `n_3!=0` and eliminate `z`.  The `xy` coefficient in `(6)` is
`4 ell_3 n_3`, hence

```text
ell_3=0.                                                 (7)
```

The remaining nonconstant coefficients are, up to the common nonzero
factor `2/n_3`,

```text
x: n_2(ell_2 n_3-ell_1 n_1),
y: ell_1(n_1 n_3-n_2^2).                               (8)
```

If `n_2=0`, a nonzero constant forces `n_1!=0` and
`ell` proportional to `(0,1,0)`.  This is Family I below.  If `n_2!=0`,
then `ell_1` cannot vanish, and `(8)` forces

```text
n_1 n_3=n_2^2,
ell_2 n_3=ell_1 n_1,
ell_3=0.                                                 (9)
```

This is Family III.

Now suppose `n_3=0`.  If `n_2!=0`, eliminate `y`.  Successively the `z`,
`x^2`, and `x` coefficients force

```text
ell_1=ell_3=ell_2=0,                                    (10)
```

contrary to `(3)`.  If `n_2=0`, then `n_1!=0`, the plane is `x=delta/n_1`,
and constancy forces only `ell_1=0`.  The constant is nonzero exactly under
the Family II condition below.  Equations `(7)`--`(10)` leave no fourth
case.

## 3. The three survivors and their polynomial inverses

Affine source-coordinate changes and invertible affine target-coordinate
changes preserve the property under study.  After such changes, the three
families are as follows.

### Family I: the `(A,d)` triangular slice

```text
n=(n_1,0,n_3),       n_1 n_3!=0,
ell proportional to (0,1,0).                            (11)
```

Use source coordinates `(d,y)` and target coordinates `(A,d)`.  Then

```text
x=(delta-n_3 d)/n_1,
y=(x^2-A)/2.                                            (12)
```

The planar Jacobian is `2`, and `(12)` is a polynomial inverse.

### Family II: a constant-`s_1` triangular slice

```text
n=(n_1,0,0),       x_0=delta/n_1,
ell=(0,ell_2,ell_3),
ell_2+2 ell_3 x_0!=0.                                   (13)
```

Use source coordinates `(y,z)` and target coordinates

```text
A=x_0^2-2y,
C=ell_3 B-ell_2 d
 =ell_3 y^2-(2 ell_3 x_0+ell_2)z.                       (14)
```

The planar Jacobian and inverse are

```text
Jac=2(ell_2+2 ell_3 x_0),
y=(x_0^2-A)/2,
z=(ell_3 y^2-C)/(2 ell_3 x_0+ell_2).                    (15)
```

### Family III: the rank-one-normal conic

Every nonzero solution of `(9)` can be represented as

```text
n=(a^2,ab,b^2),       ab!=0,
ell proportional to (b^2,a^2,0).                       (16)
```

Use source coordinates `(y,z)` and target coordinates

```text
d=z,
C=a^2 A-b^2 B.                                         (17)
```

The nonzero-constant condition is

```text
a^3+b delta!=0.                                         (18)
```

Indeed

```text
Jac=2(a^3+b delta)/a,                                   (19)

y=(delta^2-b^4 d^2-a^2 C)/(2a(a^3+b delta)),
x=(delta-ab y-b^2 d)/a^2.                               (20)
```

Equations `(20)` give the polynomial inverse.  The equalities in `(15)` and
`(18)` are sharp: when their displayed denominators vanish, the corresponding
planar Jacobian is zero rather than a hidden fourth family.

## 4. Consequence and exact boundary

Every affine-plane/affine-linear-projection Keller slice of `(1)` is one of
Families I--III and hence is a polynomial automorphism of `A2`.  Therefore

```text
the complete affine-linear slice family contains no JC(2) counterexample. (21)
```

This is stronger than checking coordinate projections one at a time: the
target kernel `ell` ranges over the entire projective dual plane, while the
source normal `n` ranges over every affine plane.

The affine hypothesis is load-bearing.  THM-2696's sharp level

```text
z=xy-c                                                       (22)
```

is nonlinear and therefore is not a missing case of this classification.
Likewise, a polynomial target projection can have a nonconstant differential
kernel and is not represented by one fixed `ell`.  The theorem does not
classify nonlinear source surfaces, nonlinear target maps, arbitrary `S4`
resolvents, `JC(2)`, or `DC(2)`.

## 5. Exact companion

Run

```text
python 04-computation/jacobian_s4_affine_plane_linear_projection_thm2699.py
python -O 04-computation/jacobian_s4_affine_plane_linear_projection_thm2699.py
```

Both modes must byte-match

```text
05-knowledge/results/jacobian_s4_affine_plane_linear_projection_thm2699.out.
```

The companion uses explicit `require` checks.  It derives `(5)`--`(10)`,
substitutes all three inverses, checks their Jacobians, and independently
enumerates every raw nonzero pair `(n,ell)` and every `delta` over `F_5` and
`F_7`.  The exact censuses are

```text
F_5: 976 =320+400+256,
F_7: 4572=1512+1764+1296,                               (23)
```

for Families I, II, and III respectively, with no mismatch between the
coefficient test and the three-family predicate.  The finite-field census is
an evidence cross-check; the characteristic-zero proof is the coefficient
elimination above.

An independent hostile audit rederived the projective typing and every branch
of the `n_3/n_2` split, checked the three target-row kernels and polynomial
inverses, and obtained the closed census formulas

```text
I=p(p-1)^3,       II=p^2(p-1)^2,       III=(p-1)^4.       (24)
```

It independently replayed both modes against the stored ten-line transcript
and matched both declared LF-normalized hashes.  The audit also confirmed
that `(n,delta)` must be rescaled jointly and that the nonlinear-source and
nonlinear-target exclusions are genuine scope boundaries.

QED.

