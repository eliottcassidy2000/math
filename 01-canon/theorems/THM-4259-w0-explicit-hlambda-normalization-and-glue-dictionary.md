---
id: THM-4259
title: "W=0 explicit H-lambda normalization and glue dictionary"
status: >
  PROVED RELATIVE TO THM-4230/4241 + VERIFIED-EXACT. After compatible
  choices of the explicit hidden-map and H-lambda parameters, the visible
  and hidden projections are H_lambda+H_lambda o iota=-v and
  H_lambda-H_lambda o iota=omega^2 f+(omega^2-omega)g. Hence the concrete
  glue generator h=H_lambda+v-omega^2 g satisfies
  2h=v+omega^2 f+g and Th=omega^2 h-omega f, exactly matching the normalized
  full-Hom basis used by THM-4241/4249/4253/4258. This removes the basis
  normalization obstruction to the 1,512-incidence observer evaluation, but
  performs none of those evaluations. W=0, M=12, seam entry, JC(2), and
  DC(2) remain OPEN.
source: codex-padic-density-cartier-20260826
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4241-w0-hidden-cm-saturation-and-visible-hidden-index-four-gluing
related:
  - THM-4249-w0-cyclic-projector-missing-eigenline-attachment-squeeze
  - THM-4258-w0-three-sample-attachment-recurrence-and-two-torsion-sidecar
script: 04-computation/jc23_w0_explicit_hlambda_glue_dictionary_thm4259.py
output: 05-knowledge/results/jc23_w0_explicit_hlambda_glue_dictionary_thm4259.out
script_sha256: 0661009e9e17510ab632faa267d5e5b69e1e0e2b4b2718e675a52cf1f8a44be4
output_sha256: 012b14c8f144d01b540a4970a4e65a35375d39e71af334a923550aa1a92ee29f
hash_basis: raw LF bytes
audit: >
  PASS. Characteristic-zero Groebner reductions verify the elliptic-addition
  and pullback-differential identities. THM-4241 places the hidden projection
  in L with degree 12; the present coercive bound and companion enumerate the
  resulting complete 24-vector shell, prove its reductions modulo 397 are
  distinct, and solve a nonsingular differential system to identify one
  unique lift. This is finite-candidate identification, not a modular-zero
  lift. Normal, optimized, and fixed-hash-seed outputs byte-match.
---

# THM-4259 -- `W=0` explicit `H_lambda` normalization

**PROVED RELATIVE TO THM-4230/4241 + VERIFIED-EXACT. NO ATTACHMENT
INCIDENCE IS EXCLUDED HERE.**

## 1. Statement and compatible parameters

Retain

```text
C_0:x^6+y^4=1,                   O=Z[omega],
iota:(x,y)->(x,-y),              tau:(x,y)->(zeta_6 x,zeta_4 y),
omega=zeta_12^4,                 kappa=zeta_12^5.       (1)
```

THM-4230 gives the degree-six hidden map `f=f_a` and `g=Tf`, where

```text
a^4-2a^3-2a+1=0,
s^6(2a^3+3a^2-1)=4,
t=(1+y^2)/x^3,

X_f=(s^2/2)x(t^2-a^2)/t,
Y_f=(s^3/2)y(t^2+a^3)/t.                              (2)
```

THM-4241 gives a degree-four glue witness

```text
H_lambda=(X_H,Y_H),
X_H=alpha(y-lambda)/x^2,
Y_H=epsilon Q_lambda(y)/x^3,
Q_lambda(y)=y^2-(r lambda/2)y+3lambda^2,               (3)
```

with `r^2-12r+24=0`, `lambda^4(r-9)=1`,
`alpha^3=r lambda`, and `epsilon^2=-1`.

Choose the compatible branches

```text
r=4+2(a+a^-1),
alpha=-4/(lambda(8-r)),
epsilon=zeta_12^3.                                    (4)
```

Here “choose” has a precise arithmetic meaning.  Fix the coefficient-field
place over `397` at which

```text
(zeta_12,a,s,r,lambda,alpha,epsilon)
  =(157,161,27,363,30,92,334) mod 397.
```

The successive equations for `zeta_12,a,s,lambda,alpha` are simple at this
point; lift this triangular branch by Hensel's lemma, define `r,epsilon` by
`(4)`, and embed the resulting characteristic-zero number field in `C`.
Thus Section 4 reduces these chosen coefficients, not an unrelated
finite-field point.

Dividing the quartic in `(2)` by `a^2` proves
`(a+a^-1)^2-2(a+a^-1)-2=0`, so the first equation in `(4)` satisfies the
required quadratic for `r`. The other relations imply exactly

```text
alpha^3=r lambda,
alpha^2-r lambda^2+6lambda^2=0.                        (5)
```

> **Theorem.** In the full Hom lattice of THM-4241,
>
> ```text
> H_lambda+H_lambda composed_with iota=-v,              (6)
> H_lambda-H_lambda composed_with iota
>       =omega^2 f+(omega^2-omega)g.                    (7)
> ```
>
> Therefore
>
> ```text
> h=H_lambda+v-omega^2 g                                (8)
> ```
>
> is the concrete normalized glue generator:
>
> ```text
> 2h=v+omega^2 f+g,               Th=omega^2 h-omega f. (9)
> ```

We use the same letter for a curve morphism `C_0->E_0` and its induced class
in `Hom(J(C_0),E_0)`.  Equations `(6)--(9)` are Hom-class identities; `(6)`
is stronger and holds pointwise for the displayed representatives.  A
Hom-class identity need not fix the translation constants of arbitrary
curve-map representatives.

## 2. Visible projection by exact elliptic addition

Write `H=H_lambda` and `H_i=H_lambda composed_with iota`. On
`E_0:Y^2=X^3+1`, the addition slope is

```text
mu=(Y_(H_i)-Y_H)/(X_(H_i)-X_H)
  =-epsilon alpha^2/(2x).                              (10)
```

The resulting `X` coordinate is

```text
X_(H+H_i)=(2alpha lambda-alpha^4/4)x^-2=-x^-2,         (11)
```

where `(5)` and the relations in `(3)--(4)` give the final equality. The
remaining constant in the `Y` numerator is one half of
`alpha^2-alpha^3lambda+6lambda^2`, hence vanishes by `(5)`. Therefore

```text
Y_(H+H_i)=-epsilon y^2/x^3.                            (12)
```

The visible map of THM-4241 is

```text
v=(-x^-2,epsilon y^2/x^3).                             (13)
```

Negation on `E_0` fixes `X` and negates `Y`; equations `(11)--(13)` prove
`(6)` as an equality of Hom classes.

## 3. Hidden projection by differentials

Using `x^6=(1-y^2)(1+y^2)`, equation `(2)` becomes

```text
X_f=(s^2/2)x^-2((1-a^2)+(1+a^2)y^2),
Y_f=(s^3/2)y x^-3((1+a^3)+(1-a^3)y^2).                (14)
```

Put

```text
A_0=-a^3+2a^2+3,             B_0=a^3-2a^2-1.          (15)
```

Differentiating `(14)` with `6x^5dx+4y^3dy=0` and reducing by the quartic
in `(2)` gives

```text
eta_f=f^*(dX/Y)
     =(2/(3s))x^-5(A_0+B_0y^2)dy.                     (16)
```

Precomposition by `tau` multiplies `x^-5dy` by
`zeta_12^(-10+3)=zeta_12^5=kappa` and sends `y^2` to `-y^2`, hence

```text
eta_g=kappa(2/(3s))x^-5(A_0-B_0y^2)dy.                (17)
```

THM-4241 computes

```text
eta_H=(alpha/(3epsilon))x^-5
      (3+y^4-4lambda y^3)/Q_lambda(y) dy.              (18)
```

Since `iota` sends `y` and `dy` to their negatives, exact reduction by the
relations in `(3)` gives

```text
eta_(H-H_i)=(2alpha/(3epsilon))x^-5
             (y^2+lambda^2(r-9))dy.                   (19)
```

THM-4241 proves that `ell_H=H-H_i` lies in the saturated lattice
`L=O f direct-sum O g` and has `q(ell_H)=12`. Write

```text
ell_H=A f+B g,                  A,B in O.               (20)
```

Equating the constant and `y^2` coefficients in `(16)--(19)` gives the
two-by-two system

```text
A_0(A+kappa B)=s(alpha/epsilon)lambda^2(r-9),
B_0(A-kappa B)=s(alpha/epsilon).                        (21)
```

Equality of pullback differentials determines the Hom class: a homomorphism
of complex abelian varieties with zero differential is zero.

## 4. Finite-candidate identification at the good prime 397

The hidden Hermitian form is

```text
q(Af+Bg)=6N(A)+6N(B)+Tr(A conjugate(B)(-4-2omega)).     (22)
```

The off-diagonal entry in `(22)` has norm `12`.  Cauchy--Schwarz therefore
gives

```text
q(Af+Bg)>=(6-2sqrt(3))(N(A)+N(B)).                    (23)
```

If `q=12` and the integer `N(A)+N(B)` were at least five, the right side
would be at least `30-10sqrt(3)>12`.  Hence `N(A)+N(B)<=4`.  Since
`N(m+n omega)=(m-n/2)^2+3n^2/4` (and symmetrically with `m,n` exchanged),
this forces every Eisenstein coordinate into `[-2,2]`.  Exact enumeration of
that proved finite universe gives precisely `24` coefficient pairs with
`q=12`.  Thus `(20)` is known *before reduction* to be one member of a finite
explicit set.

Use the exact specialization

```text
p=397,
(zeta_12,a,s,r,lambda,alpha,epsilon)
 =(157,161,27,363,30,92,334),
omega=34,                         kappa=177.             (24)
```

The companion verifies that `397` is prime, every equation in `(1)--(5)`, all needed
denominators, exact order twelve for `zeta_12`, and the nonzero derivatives
of the triangular equations for `zeta_12,a,s,lambda,alpha`; `r` and
`epsilon` are then defined by `(4)`.  Hensel's lemma therefore supplies the
chosen characteristic-zero algebraic branch and its reduction map.  The
determinant of `(21)` is

```text
-2 kappa A_0 B_0=313 mod 397 !=0.                      (25)
```

Solving `(21)` gives

```text
(A,B)=(362,328) mod 397.                               (26)
```

All `24` degree-twelve coefficient pairs have distinct reductions in
`F_397^2`. Since

```text
omega^2=362,
omega^2-omega=328 mod 397,                             (27)
```

the unique lift of `(26)` is

```text
(A,B)=(omega^2,omega^2-omega).                         (28)
```

This proves `(7)`. The logical direction is crucial: this is not “an
identity vanishes modulo one prime, therefore it vanishes in characteristic
zero.” The characteristic-zero theorem first confines `(A,B)` to 24 exact
possibilities, and reduction is injective on that set. The good-prime system
then identifies which possibility occurs.  More precisely, the proof reduces
the characteristic-zero coefficient equations `(21)` at the fixed
coefficient-field place.  It does not infer equality of Hom classes from an
equality seen only on a special fibre.  The prime is larger than three and
all displayed parameter denominators are units, so the explicit formulas
also have good coefficient reduction, but no specialization-injectivity
claim is needed.

## 5. Concrete basis and cyclic action

Adding `(6)` and `(7)` gives

```text
2H_lambda=-v+omega^2 f+(omega^2-omega)g.               (29)
```

Using `omega^2+omega=-1`, definition `(8)` now gives the first identity in
`(9)`. It also differs from `H_lambda` by an element of
`D=V direct-sum L`, so it represents the same nontrivial glue class.

Apply `T` to `2h=v+omega^2f+g` and use

```text
Tv=omega^2v,                   Tf=g,       Tg=-omega f. (30)
```

For completeness, these characters are visible in the displayed formulas.
The square of `tau` fixes `t` and sends `(x,y)` to `(omega x,-y)`, so `(2)`
gives `T^2f=-omega f`; with `g=Tf`, this is `Tg=-omega f`.  Formula `(13)`
under `tau` has its `X` coordinate multiplied by `omega^2` and its `Y`
coordinate fixed, which is exactly `Tv=omega^2v` on `E_0`.

The result equals `2(omega^2h-omega f)`. The Hom lattice is torsion-free, so
multiplication by two is injective and the second identity in `(9)` follows.
Thus `(8)` is not merely some half-sum: it is the concrete basis vector with
the exact relation and `T` matrix used in the current residual censuses.

For every residual row `m=b f+c g+d h`, equation `(8)` gives the following
Hom-class formula:

```text
m=d H_lambda+d v+b f+(c-d omega^2)g.                   (31)
```

Every curve-map representative on the right has an explicit rational
formula.  On the degree-zero divisor `[tau Q-Q]`, translations cancel:
an induced Hom class represented by `M:C_0->E_0` evaluates as
`M(tau Q)-M(Q)`.  Hence THM-4258's translation-normalized observer

```text
F_m(Q)=m([tau Q-Q])                                    (32)
```

can now be evaluated without an unidentified two-torsion half or basis
change.

## 6. Scope and next computation

This theorem removes the normalization obstruction, but it does not perform
the `1,512` map-ratio incidence evaluations, prove either marked-ratio set
empty, close `W=0` or exact `M=12`, prove seam entry, `JC(2)`, or `DC(2)`.
The next exact computation is to combine `(24),(31),(32)` with each
map-specific torsion ratio and the canonical twelve-node orbit; THM-4258
reduces the uniform direct table to `4,536` group-value rows.

```bash
python3 -B 04-computation/jc23_w0_explicit_hlambda_glue_dictionary_thm4259.py
python3 -B -O 04-computation/jc23_w0_explicit_hlambda_glue_dictionary_thm4259.py
PYTHONHASHSEED=4259 python3 -B \
  04-computation/jc23_w0_explicit_hlambda_glue_dictionary_thm4259.py
```

All three streams byte-match the frozen output. **QED.**
