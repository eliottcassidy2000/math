---
id: THM-3696
title: "Three-branch conductor and graded modules of the y=0 collision ring"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the Danielewski envelope
  D=C[b,c,e]/(c^2e-(1-b^2)), the collision ring is the graded subalgebra
  R=C[e,ce,bc].  Its weight-zero ring is the ordinary triple curve
  R_0=C[b(1-b^2),b^2(1-b^2)], whose normalization C[b] glues b=-1,0,1
  by equality of values and one exact derivative law.  Its conductor is
  b^2(1-b^2)^2C[b].  Every graded coefficient module of R is classified
  explicitly.  The global conductor R:D is zero, R is not Poisson-closed,
  and scalar brackets have a synchronized three-branch law: only weights
  (-2,1) survive at b=+-1, with equal values on both arms, while only
  (-2,1) and (-1,0) survive at b=0.  This constrains but does not close the
  remaining 2x5 and 3x4 Darboux cells or JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  PASS.  The independent audit checked the exact kernel and normalization,
  ordinary triple tangents, value-and-jet iff reconstruction, conductor
  converse, every semigroup boundary and Hermite module, the zero global
  conductor and birational nonfinite consequence, the non-Poisson witness,
  and all three-branch scalar evaluations.  It confirmed that the finite
  hostile window is only a control and that the written proofs are all-degree.
depends_on:
  - THM-3686-y0-collision-normalization-and-bracket-anatomy
related:
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
  - THM-3592-universal-exponent-two-three-by-three-weight-darboux-nonentry
  - THM-3695-y0-collision-ring-danielewski-embedding-and-seven-piece-floor
script: 04-computation/jacobian_y0_collision_three_branch_conductor_thm3696.py
output: 05-knowledge/results/jacobian_y0_collision_three_branch_conductor_thm3696.out
script_sha256: 2c5312c8d5696e8b2e648bbab9bc538b5aa10c4a46bce6b63f8c3cf09fc7ff84
output_sha256: 8e6bdac38643d420595bc55842dc63069867d04599cfb55bbb463a34af5d0cbb
hash_basis: LF-normalized bytes
---

# THM-3696 -- the collision ring is a three-branch gluing inside its envelope

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

All rings are over `C`.  Put

```text
h=1-b^2,
D=C[b,c,e]/(c^2e-h),              wt(b,c,e)=(0,1,-2),  (1)
R=C[e,ce,bc] subset D.                                      (2)
```

The harmless nonzero scalars in `A=3e`, `B=2ce`, `C=bc` have been removed.
THM-3686 identifies `(2)` with its `y=0` collision ring on the source chart.
For each integer `r`, define its coefficient module

```text
M_r={p(b) in C[b] : c^r p(b) in R}.                    (3)
```

For negative `r`, the expression in `(3)` is read in `D[c^-1]`; the stated
membership makes it regular in `D`.

The inclusion `(2)` is graded, and the Poisson bracket on `D` agrees with
the source bracket on pairs from `R`.  It is important that `R` is **not** a
Poisson subalgebra: Section 5 gives an explicit failed closure bracket.

## 1. The weight-zero curve and its normalization

Set

```text
X=bh=b-b^3,                         Y=b^2h=b^2-b^4.     (4)
```

A monomial `e^i(ce)^j(bc)^k` has weight and Laurent form

```text
r=k-j-2i,
e^i(ce)^j(bc)^k=c^r b^k h^(i+j).                       (5)
```

At weight zero, `k=j+2i`, so `(5)` is `X^jY^i`.  Hence

```text
R_0=M_0=C[X,Y] subset C[b].                             (6)
```

The parametrization satisfies

```text
X^4-X^2Y+Y^3=0.                                        (7)
```

This is the exact kernel relation.  Indeed, `b=Y/X`, and `b^3-b+X` is
irreducible over `C(X)` (a rational root would be a polynomial divisor of
`X`, and none works).  Thus `C(X,Y)=C(b)` has degree three over `C(X)`, while
`(7)` is cubic in `Y`.  Also `b` is integral over `R_0` by `b^3-b+X=0`.
Therefore `C[b]` is the normalization.  The only normalization points over
`(X,Y)=(0,0)` are

```text
b=-1,0,1.                                               (8)
```

The degree-three initial form of `(7)` is

```text
Y(Y-X)(Y+X),                                            (9)
```

so this is an ordinary triple point with three distinct tangent lines.

## 2. Exact value-and-jet gluing

For every `f in C[b]`,

```text
f in R_0
iff
f(-1)=f(0)=f(1),
f'(-1)+4f'(0)+f'(1)=0.                                (10)
```

Necessity follows by differentiating a polynomial in `(X,Y)` along the
three branches.  For sufficiency, let the common value be `alpha`, put

```text
beta=f'(0),              gamma=(f'(-1)+2beta)/2.       (11)
```

The polynomial

```text
f-(alpha+beta X+gamma Y)                               (12)
```

has a double zero at each point in `(8)`, by `(10)`.  It is therefore
divisible by

```text
X^2=b^2h^2.                                            (13)
```

Now `C[b]` is generated over `R_0` by `1,b,b^2`, and

```text
X^2 in R_0,              bX^2=XY in R_0,
b^2X^2=Y^2 in R_0.                                     (14)
```

Thus `X^2C[b] subset R_0`, proving sufficiency.

The same argument determines the conductor exactly:

```text
cond_(C[b]/R_0)=X^2C[b]=b^2(1-b^2)^2C[b].             (15)
```

Indeed `(14)` proves containment.  Conversely, if `fC[b] subset R_0`, apply
the value conditions in `(10)` to `fg` while prescribing the three values
of `g` independently.  This forces `f(-1)=f(0)=f(1)=0`.  Apply the derivative
condition to `fg`; independent values of `g` then force all three derivatives
of `f` to vanish.  Hence `X^2|f`.

## 3. Complete graded-module classification

The exact answer for `(3)` is

```text
r>=0:
  M_r=b^r R_0;

r=-1:
  M_-1=h(R_0+bR_0)
      =h{q in C[b]: q(-1)-2q(0)+q(1)=0};

r=-2m, m>=1:
  M_-2m=h^m(R_0+hR_0)
        =h^m{q in C[b]: q(-1)=q(1)};

r=-(2m+1), m>=1:
  M_-(2m+1)=h^(m+1)C[b].                              (16)
```

Here is the all-degree semigroup proof.  If `r>=0`, equation `(5)` gives
`k=r+j+2i`, and hence

```text
b^k h^(i+j)=b^r X^jY^i.                               (17)
```

For `r=-R<0`, the allowed lattice points satisfy `2i+j>=R`.  Reducing by
the two `R_0` moves `(i,j)->(i+1,j)` and `(i,j)->(i,j+1)` leaves precisely
the following minimal normalized generators:

```text
R=2m:    h^m times {1,h,...,h^m};
R=2m+1:  h^(m+1) times {1,h,...,h^m,b}.                (18)
```

For even `R`, the module `R_0+hR_0` is exactly the first hyperplane in

```text
R_0+hR_0={q:q(-1)=q(1)}.                               (19)
```

It already contains all higher powers of `h`.  For `R=1`, direct Hermite
reduction modulo the conductor gives

```text
R_0+bR_0={q:q(-1)-2q(0)+q(1)=0}.                      (20)
```

For odd `R>=3`, the normalized generators include `1,b,h`; because
`b^2=1-h` and `C[b]=R_0{1,b,b^2}`, they generate all of `C[b]`.  Equations
`(17)--(20)` prove every line of `(16)`, including the exceptional boundary
weight `-1`.

The two Hermite hyperplanes are exact, not just necessary.  Modulo `X^2`,
the five jets of

```text
1,X,Y,h,hX                                               (21)
```

span `q(-1)=q(1)`, while the five jets of

```text
1,X,Y,b,bY                                               (22)
```

span `q(-1)-2q(0)+q(1)=0`.  The companion performs both rank and annihilator
checks over `Q`.

## 4. The full extension has zero conductor

Although every coefficient module is controlled by the finite conductor
`(15)`, the full graded extension has

```text
(R:D)={a in D:aD subset R}=0.                           (23)
```

If `a!=0`, choose a nonzero homogeneous coefficient `c^r p(b)` of `a`.
For all sufficiently large `N`, membership of `c^Na` in `R` would force,
by the positive line of `(16)`,

```text
b^(r+N) divides p(b).                                  (24)
```

This is impossible for fixed nonzero `p` as `N` grows.  Hence no nonzero
element multiplies all of `D` back into `R`.  The extension is birational:
`c=(ce)/e` and `b=(bc)/c` in the common fraction field.  A finite birational
extension of these noetherian domains would have nonzero conductor after
clearing denominators of finitely many generators.  Thus `(23)` also proves
that `R subset D` is not finite.

## 5. Non-Poisson closure and the synchronized scalar law

In `D`, with the Danielewski bracket,

```text
{3e,2ce}=-12be.                                        (25)
```

The coefficient of `be=c^-2(bh)` at weight `-2` is `h q` with `q=b`.
It violates `q(-1)=q(1)` in `(16)`, so `be notin R`.  Thus `(25)` is an
explicit witness that `R` is graded but not Poisson-closed.  All Darboux
questions in `R` use the ambient bracket in `D`, equivalently the source
bracket of THM-3686.

For complementary homogeneous weights, after exchanging the two entries if
necessary, write

```text
c^-R p(b),                 c^(R-1)q(b),
J_R(p,q)=R p q'+(R-1)p'q.                              (26)
```

The module table `(16)` gives the exact three-branch evaluation law:

```text
R=1:
  p=hf, q in R_0;
  J_1(+-1)=0,                 J_1(0)=f(0)q'(0);

R=2:
  p=hf, f(-1)=f(1),          q=bg, g in R_0;
  J_2(-1)=J_2(1)=-2f(1)g(1),
  J_2(0)=2f(0)g(0);

R>=3:
  J_R(-1)=J_R(0)=J_R(1)=0.                              (27)
```

Reversing the bracket changes all signs but none of the vanishing or equality
statements.  Consequently:

1. at the two Danielewski arms `b=+-1`, a scalar row can be supplied only by
   `(-2,1)` addresses;
2. every individual `(-2,1)` address contributes the same value on both
   arms, so arm alternation is impossible in `R`; and
3. at the third glued point `b=0`, only `(-2,1)` and `(-1,0)` addresses can
   contribute.

This is an extra incidence sidecar beyond the universal Danielewski support
theorems.  It acts directly on the first inherited open cells `2x5` and
`3x4`, but no claim here closes either cell or proves `JC(2)`.

## 6. Exact reproduction

Run

```bash
python3 -B 04-computation/jacobian_y0_collision_three_branch_conductor_thm3696.py
python3 -B -O 04-computation/jacobian_y0_collision_three_branch_conductor_thm3696.py
```

Both modes byte-match the stored transcript.  The exact companion checks the
plane relation, normalization jets, conductor, all module boundary cases on
a hostile degree window, non-Poisson witness, and scalar evaluations.  The
semigroup and arbitrary-degree conductor arguments are the proofs above.

**QED.**
