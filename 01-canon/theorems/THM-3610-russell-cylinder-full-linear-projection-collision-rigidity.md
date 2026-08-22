---
id: THM-3610
title: "Russell-cylinder full linear-projection collision rigidity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a
  polynomial cylinder graph retaining the full THM-3561 triple
  collision, no rank-two linear projection of all four Russell coordinates
  (B,C,Y,S) has nonzero constant ordinary Jacobian.  Projection planes that
  meet span(B,C) are excluded without the collision hypothesis by the
  fixed-C gate or a formal-arm transport equation.  A plane transverse to
  span(B,C) is injective on the source line x=0; Gwozdziewicz's one-line
  theorem would make a Keller projection an automorphism, contradicting the
  retained collision.  No nonlinear projection, implicit non-graph source
  plane, or JC(2) counterexample is claimed.
source: root / omitted-Y full-cylinder projection niche, 2026-08-21
audit: >
  PASS.  An independent hostile audit rederived the complete row-space
  trichotomy, the Y-only boundary ODE, the S+sigmaY determinant and formal
  recurrence, and the transverse-line reconstruction.  It checked
  Gwozdziewicz's one-line theorem against the primary source with every
  hypothesis typed, and verified the retained collision contradiction.
  Normal and optimized runs are byte-identical to the stored 15,028-gate
  transcript; the AST has no assertion gates, and documentation and diff
  checks pass.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3607-russell-cylinder-mixed-projection-degree-seven-gate
related:
  - THM-3589-danielewski-central-arm-every-line-and-kummer-trace-darboux-gates
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
  - THM-3608-russell-cylinder-nonlinear-target-shear-rigidity
external:
  - "Gwozdziewicz, Injectivity on one line, arXiv:alg-geom/9305008, Theorem 1.1."
script: 04-computation/jc2_russell_cylinder_full_linear_projection_thm3610.py
output: 05-knowledge/results/jc2_russell_cylinder_full_linear_projection_thm3610.out
script_sha256: d83e4452a94ff1a48f43e3feda7380e56b269385503bf9de8e6fca82bc3b0094
output_sha256: c96249605b3899195483da25e983c0776966c5a4693dc112dc373943bb7c8a65
hash_basis: raw LF bytes
---

# THM-3610 -- Russell-cylinder full linear-projection collision rigidity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The result
below closes every rank-two *linear* projection of the four
Russell-cylinder coordinates once the full stable collision is retained.
It does not close nonlinear target projections or implicit non-graph source
surfaces.

All rings and derivatives are over `C`.

## 0. Statement and ledger

Use the THM-3561 functions

```text
D=1+x^2q,
a=q/D^2,
c=xD(D+2),
b=ac^2,
e=a(b+4),                       Jac_(x,q)(a,c)=-3.       (1)
```

For an arbitrary polynomial graph `w=h(x,q)`, the Russell cylinder map of
THM-3605 becomes

```text
B=b+ch,
C=c,
Y=ce+(2b+4)h+ch^2,
S=((b+2)(e+3h^2)+ch(3e+h^2))/8.                        (2)
```

Let

```text
r_0=(0,-3/4),             r_+=(1,-3),
r_-=(-1,-3).                                                (3)
```

These are the three distinct source points in the THM-3561 collision.  If

```text
h(r_0)=h(r_+)=h(r_-),                                      (4)
```

then the full four-coordinate stable map `(B,C,Y,S)` has the same value at
all three points.  The claim is that for every rank-two matrix
`Lambda in Mat_(2x4)(C)`,

```text
(F_Lambda,G_Lambda)^T=Lambda(B,C,Y,S)^T,
Jac_(x,q)(F_Lambda,G_Lambda) notin C*.                    (5)
```

The source / target / loss ledger is:

| item | content |
|---|---|
| source | a polynomial graph in the stabilized source `A3_(x,q,w)` |
| target | an arbitrary rank-two linear projection of `(B,C,Y,S)` |
| retained | polynomial outputs and, by `(4)`, the full stable collision |
| discarded | two target directions and the stable three-form |
| tested predicate | nonzero constant planar Jacobian |
| decisive sidecars | formal completion of `x=0`; injectivity on that line |
| remaining exits | nonlinear projection, implicit non-graph surface, or a projection not retaining `(4)` |

## 1. Arm coordinates and the row-space trichotomy

On `D!=0`, put

```text
A=ac+h.                                                     (6)
```

Substitution in `(2)` gives the four remarkably small expressions

```text
B=cA,
C=c,
Y=4A+cA^2,
S=a+3A^2/4+cA^3/8.                                        (7)
```

Let

```text
V_0=span_C(B,C),              U=span_C(Y,S),
W=rowspan_C(Lambda).                                       (8)
```

Since `V_0 direct-sum U` is the four-dimensional coordinate space, the
integer

```text
r=dim_C(W intersect V_0)                                  (9)
```

is `0`, `1`, or `2`.

If `r=2`, then `W=V_0`; the pair is linearly equivalent to `(B,C)`, and
THM-3607 gives

```text
Jac(B,C)=c(-3c+Jac(h,c)),                                  (10)
```

which cannot lie in `C*`.

Suppose `r=1`.  A line in `V_0` is either `C` or, after scaling,
`B+rho C`.  If it is `C`, a constant determinant would give a polynomial
Jacobian mate for `c`; the fixed-`C` theorem in THM-3607 excludes this for
every polynomial second output.  The remaining line is treated in the next
two sections.

Finally, if `r=0`, projection `W -> U` is an isomorphism.  After a target
`GL_2(C)` change of basis, the outputs have the form

```text
F=Y+alpha B+beta C,
G=S+gamma B+delta C.                                      (11)
```

Section 4 kills precisely this transverse branch using `(4)`.

## 2. One `Y` direction and no `S` direction

Take the intersection generator and second row in the normalized form

```text
L=B+rho C=cZ,              Z=A+rho,
M=Y+tau C=g(Z,c),                                          (12)
```

where

```text
g(Z,c)=4(Z-rho)+c(Z-rho)^2+tau c.                         (13)
```

Regard `Z` as a function of the etale arm coordinates `(a,c)`.  A direct
differentiation gives

```text
det partial(L,M)/partial(a,c)
   =(c g_c-Z g_Z) Z_a,
Jac_(x,q)(L,M)=3(Z g_Z-c g_c)Z_a.                         (14)
```

On the line `x=0`, one has `a=q,c=0`; write `Z(a,0)=f(a)`.  If the second
quantity in `(14)` were the nonzero constant `3t`, then

```text
4 f f'=t.                                                 (15)
```

No polynomial `f` solves `(15)` for `t!=0`: a nonconstant `f` makes the
left side have degree `2 deg(f)-1`, while a constant `f` makes it zero.
Thus this branch is empty without using the collision.

## 3. Both `Y` and `S` directions

The last `r=1` branch can be normalized as

```text
L=B+rho C=cZ,
M=S+sigma Y+tau C=a+g(Z,c),                              (16)
```

with

```text
g(Z,c)=3(Z-rho)^2/4+c(Z-rho)^3/8
       +sigma(4(Z-rho)+c(Z-rho)^2)+tau c.                (17)
```

The same determinant cancellation used in THM-3607 gives

```text
Jac_(x,q)(L,M)
 =3[Z+cZ_c+(Z g_Z-c g_c)Z_a].                            (18)
```

Assume that `(18)` equals `3t`, `t!=0`.  On `x=0`, with
`f(a)=Z(a,0)`, equation `(18)` is

```text
f+f[3(f-rho)/2+4sigma]f'=t.                              (19)
```

If `f` is nonconstant, the term `(3/2)f^2f'` has degree
`3 deg(f)-1` and has no competing term of that degree.  Hence `f` is
constant, and `(19)` forces `f=t`.

The formal completion from THM-3607 is explicit.  In `C[a][[c]]`, the
unique series `R=1+O(c^2)` satisfying

```text
ac^2=(R-1)(R+2)^2                                        (20)
```

gives

```text
x=c/[R(R+2)],                    q=aR^2.                  (21)
```

Thus `Z` has a well-defined expansion in `C[a][[c]]`.  If

```text
Z-t=z_N(a)c^N+O(c^(N+1)),          N>=1,                 (22)
```

with `z_N!=0`, the first nonzero coefficient of `(18)` is

```text
(N+1)z_N+t[3(t-rho)/2+4sigma]z_N'=0.                    (23)
```

The first term has degree `deg z_N`; the derivative term has smaller
degree.  Hence `(23)` has no nonzero polynomial solution.  Therefore
`Z=t` in the completion and in the source function field.  Equations
`(6),(12)` force

```text
h=t-rho-ac=t-rho-xq(D+2)/D.                              (24)
```

Modulo `D`, the numerator `xq(D+2)` is nonzero, so `(24)` has a genuine
simple pole along `D=0`.  It cannot be a polynomial graph.  This empties
the final `r=1` branch, again without using `(4)`.

## 4. Transverse projections are injective on one line

It remains to treat `r=0`, so the outputs have form `(11)`.  On the affine
source line

```text
ell={x=0},             a=q, c=B=C=0,
A=h(0,q)=:phi(q).                                         (25)
```

Equations `(7),(11)` restrict to

```text
(F,G)|ell=(4phi(q), q+3phi(q)^2/4).                       (26)
```

This map `A1->A2` is injective: equality of the first coordinates gives
equality of `phi`, after which equality of the second gives equality of
`q`.

If `(F,G)` had nonzero constant Jacobian, Gwozdziewicz's cited one-line
theorem would make it an automorphism of `A2`, hence injective everywhere.
But at the three points `(3)`, one has

```text
a=-3/4,        b=c=B=C=0,
Y=4eta,        S=-3/4+3eta^2/4,
eta=h(r_0)=h(r_+)=h(r_-).                                (27)
```

Thus `(F,G)` has the same value at three distinct points, contradicting
injectivity.  This proves `(5)`.

The collision hypothesis is load-bearing only here.  Without `(4)`, the
`r=0` transverse branch remains open; Sections 2--3 still close every
projection plane meeting `span(B,C)`.

## 5. Boundaries and no-overclaim clause

1. The theorem is about **linear** target projections.  Polynomial target
   shears are the separate THM-3608 lane.
2. It is about polynomial source graphs `w=h(x,q)`.  An implicit embedded
   `A2` in the cylinder is not a graph and remains outside the statement.
3. The condition `(4)` retains the full four-coordinate collision.  A
   projected pair may collide under weaker conditions, but those are not
   used here.
4. Two out-of-scope boundary controls show that the cylinder itself contains
   Keller planes: on the non-graph plane `x=0`,
   `(Y,S)=(4w,q+3w^2/4)` is triangular in `(q,w)`; on the non-graph plane
   `q=0`, the nonlinear inverse cylinder coordinate paired with `C` is also
   a constant-Jacobian plane.  Neither plane retains the triple collision.
   These controls do not settle the `r=0` polynomial-graph branch without
   condition `(4)`.
5. No Darboux pair on the Danielewski target and no counterexample to
   `JC(2)` is constructed.

## 6. Exact companion contract

The companion must check, without truth-bearing `assert` statements:

- the four arm-coordinate identities `(7)` and `Jac(a,c)=-3`;
- the normalized determinant identities `(14)` and `(18)` for symbolic
  `rho,sigma,tau`;
- the boundary equations `(15),(19)` and the first formal coefficient
  `(23)` across an explicit degree/order grid;
- the exact three collision points, their common `(a,b,c,e)` data, and
  `(27)` under condition `(4)`;
- injectivity reconstruction from `(26)`;
- a finite exact Grassmannian row-space census confirming the `r=0,1,2`
  trichotomy and both normalized `r=1` branches;
- active hostile controls showing that `(Y,S)` on `x=0` and the inverse
  arm coordinate on `q=0` are Keller on collision-losing non-graph planes.

Passing that companion verifies the displayed algebra and the finite
row-space classification.  The all-degree conclusions use the polynomial
degree and formal-completion proofs above; no finite computation is treated
as a proof of `JC(2)`.
