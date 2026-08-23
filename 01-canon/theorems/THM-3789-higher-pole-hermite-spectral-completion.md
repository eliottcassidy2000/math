---
id: THM-3789
title: "Higher-pole Hermite spectral completion"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.  For every
  squarefree polynomial q of degree d>=2 whose normalized primitive Q of q^2
  has pairwise distinct values at 0 and the roots of q, the rational Keller
  pair F=xq(t), P=Q(t)/F^3, t=x^2(1+xy), has an explicit smooth Hermite
  completion.  The resulting quasi-finite etale source map misses exactly d
  codimension-two points, and its coordinate ring is exactly
  k[x,y] intersect k(F,P).  The generic source degree is 2d+1.  Degree d=1
  is the sharp failure: a whole divisor is omitted and x^3 is an additional
  target-field polynomial, leading to THM-3785 instead.  The completed
  surface has units k* and Picard group Z^d.  More decisively, its moving
  axis arm gives a nonzero algebraic de Rham class, so it has no Darboux pair
  at all.  The entire distinct-critical higher-pole family is therefore
  closed as a planar Jacobian-counterexample lane.
source: jc_quartic_c3_construct / higher-pole Hermite-completion lane, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The proof checks every UFD valuation,
  boundary restriction, residual-sheet count, Poisson sign, etale minor,
  image stratum, height-one DVR descent, moving-arm Cech--de Rham residue,
  and the degree-one hostile boundary.
  The exact companion verifies the formal gradient packet, exact quadratic
  and cubic positive controls, all boundary constants, critical residuals,
  and the degree-one missing-divisor witness.  Normal and optimized runs
  byte-match the frozen transcript.  Independent hostile audit remains due.
depends_on: []
related:
  - THM-3782-simple-pole-spectral-danielewski-completion-and-target-field-gate
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
script: 04-computation/jc2_higher_pole_hermite_spectral_completion_thm3789.py
output: 05-knowledge/results/jc2_higher_pole_hermite_spectral_completion_thm3789.out
script_sha256: 131c23b59985fe86bf60356eef186be074a91250b3fe9bdc8fe978c8763c1bcc
output_sha256: a530de749cc4a29329d76a4d1481ecc75ae22670faf7314bb531288e34c0ded1
semantic_sha256: 58fde6b52fe4578249bc3413eed29dbf12a753077a39faeaacd0c8b7f1000e5a
hash_basis: raw LF bytes
---

# THM-3789 -- higher poles admit an all-degree Hermite completion

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.**  This
theorem identifies a new infinite family of smooth polynomial target-field
completions.  It does not give a polynomial Keller pair or a planar Jacobian
counterexample.  Instead, it computes the complete polynomial part of an
all-degree family of rational Keller target fields and isolates the exact
degree-one failure that produced the exceptional cubic pseudo-plane of
THM-3785.

Let `k` be an algebraically closed field of characteristic zero.  In
`k[x,y]`, put

```text
t=x^2(1+xy).                                                   (1)
```

Let `q in k[T]` be squarefree of degree `d>=2`, with `q(0)!=0`, and let

```text
Q(T)=integral from 0 to T of q(s)^2 ds.                         (2)
```

Write the distinct roots of `q` as `alpha_1,...,alpha_d`, and set

```text
lambda_i=Q(alpha_i).                                           (3)
```

The load-bearing spectral hypothesis is

```text
0,lambda_1,...,lambda_d are pairwise distinct.                 (4)
```

Define

```text
F=xq(t),                 W=Q(t),                 P=W/F^3,
Psi(S)=product_(i=1)^d (S-lambda_i).                            (5)
```

Then the two rational functions

```text
H=W Psi(W)/F^2,
Y=[H-Psi(W)]/F                                             (6)
```

are source polynomials.  They satisfy

```text
F^2H=W Psi(W),
H=Psi(W)+FY,
F^3Y=(W-F^2)Psi(W).                                           (7)
```

Consequently the homomorphism

```text
B=k[f,w,z]/(f^3z-(w-f^2)Psi(w))  -->  k[x,y],
(f,w,z) |-> (F,W,Y)                                             (8)
```

is injective.  If `S=Spec B`, the induced map

```text
phi:A2_(x,y) -> S                                               (9)
```

is quasi-finite and etale, with exact geometric image

```text
S \ {p_1,...,p_d},              p_i=(0,lambda_i,0).             (10)
```

The missing set has codimension two.  The complete polynomial intersection
and the generic field degree are

```text
k[x,y] intersect k(F,P)
 = B
 = k[F,W,Y],

[k(x,y):k(F,P)]=2d+1.                                         (11)
```

Thus `(8)` is not merely a convenient polynomial subring: it is the maximal
polynomial observable carried by the rational Keller target field.  Moreover,

```text
B*=k*,                 Pic(S)=Z^d.                              (11a)
```

The inverse symplectic form has nonzero algebraic de Rham class:

```text
[omega]!=0 in H^2_dR(S).                                       (11b)
```

Consequently there are no `A,C in B` with `{A,C}=1`.

## 1. The rational Keller seed

For the Jacobian convention `{A,C}=A_x C_y-A_y C_x`, direct differentiation
of `(1)` gives

```text
{F,t}=x^3q(t).                                                  (12)
```

Since `Q'=q^2`, this implies

```text
{F,W}=q(t)^2 {F,t}=x^3q(t)^3=F^3,
{F,P}={F,W/F^3}=1.                                             (13)
```

Thus `(F,P)` is a rational Keller pair with a third-order pole along
`F=0`.  Also

```text
k(x,y)=k(F,t),
x=F/q(t),                 y=[t-x^2]/x^3.                       (14)
```

In particular `F,t` are algebraically independent.

## 2. Why the two Hermite quotients are polynomial

Factor

```text
q(T)=c_q product_i (T-alpha_i).                                (15)
```

The height-one source components of `F=0` are the axis `A=(x)` and

```text
D_i=(t-alpha_i).                                                (16)
```

Because `alpha_i!=0`, the polynomial `t-alpha_i` is primitive and linear in
`y`, with coprime coefficients `x^3` and `x^2-alpha_i`; hence every `D_i` is
prime.  These components and `A` are pairwise coprime.

At the axis, `t=x^2(1+xy)` and `Q'(0)=q(0)^2!=0`, so

```text
ord_A(F)=1,              ord_A(W)=2,              ord_A(Psi(W))=0.
                                                                    (17)
```

At `D_i`, write `u=t-alpha_i`.  Squarefreeness gives

```text
q(t)=q'(alpha_i)u+O(u^2),
W-lambda_i=q'(alpha_i)^2 u^3/3+O(u^4).                         (18)
```

Hypothesis `(4)` makes `W` and all other spectral factors units there.
Therefore

```text
ord_(D_i)(F)=1,          ord_(D_i)(W Psi(W))=3.                 (19)
```

Equations `(17),(19)` show in the UFD `k[x,y]` that `F^2` divides
`W Psi(W)`.  Hence `H` is polynomial.  Its boundary values are

```text
H mod A=Psi(0),                 ord_(D_i)(H)=1.                 (20)
```

The first equality follows by dividing the leading terms
`W=q(0)^2x^2(1+xy)+O(x^4)` and `F^2=q(0)^2x^2+O(x^4)`.
Thus `H-Psi(W)` vanishes on `A`; on each `D_i`, both terms vanish.  The
pairwise-coprime factorization `(15),(16)` now shows that `F` divides this
difference.  This proves polynomiality of `Y` and all identities `(7)`.

The relation in `(8)` is irreducible: viewed as a primitive polynomial linear
in `z`, its two coefficients `f^3` and `(w-f^2)Psi(w)` are coprime.  Its image
has transcendence degree two by `(14)`, so `(7)` is the complete relation and
`(8)` is injective.

## 3. Smoothness and the Jacobian Poisson packet

Let

```text
D=f^3z-(w-f^2)Psi(w).                                         (21)
```

The gradient is

```text
D_z=f^3,
D_f=3f^2z+2fPsi(w),
D_w=-Psi(w)-(w-f^2)Psi'(w).                                  (22)
```

If `f!=0`, then `D_z!=0`.  If `f=0`, equation `D=0` forces
`w=0` or `w=lambda_i`.  At these arms,

```text
D_w(0,0,z)=-Psi(0)!=0,
D_w(0,lambda_i,z)=-lambda_i Psi'(lambda_i)!=0.                (23)
```

Thus `S` is smooth and therefore normal.  Differentiating `(7)`, or checking
directly in the source, gives the full Poisson packet

```text
{F,W}=F^3,
{F,Y}=Psi(W)+(W-F^2)Psi'(W),
{W,Y}=3F^2Y+2FPsi(W).                                        (24)
```

These three minors are respectively `D_z,-D_w,D_f`; hence they generate the
unit ideal on `S`.  The source map has full differential rank everywhere.

## 4. Exact boundary packets and exact image

The first nonzero axis jet and the critical-arm restrictions are

```text
Y(0,y)=Psi(0)y/q(0),                                          (25)

Y |_(D_i)
 = lambda_i Psi'(lambda_i)/(3q'(alpha_i)x^3),
y=(alpha_i-x^2)/x^3,                 x in G_m.                 (26)
```

Every coefficient in `(25),(26)` is nonzero by `(4)` and squarefreeness.
Thus the axis maps onto the whole arm `(f,w)=(0,0)`.  Each `D_i` maps onto
the punctured critical arm

```text
(f,w)=(0,lambda_i),                 z!=0,                      (27)
```

and misses exactly `p_i`.

It remains to check that no off-boundary point is lost.  Let `f!=0` and fix a
root `beta` of

```text
Q(T)-w=0.                                                       (28)
```

If `w` is not a critical value, every root has `q(beta)!=0`.  If
`w=lambda_i`, then `alpha_i` is an exact triple root by `(18)`.  The polynomial
in `(28)` has degree `2d+1`, so its residual factor has degree

```text
2d-2>=2.                                                        (29)
```

By `(4)`, no residual root is another critical point.  Hence in all cases a
root with `q(beta)!=0` exists.  Put

```text
x=f/q(beta),                 y=[beta-x^2]/x^3.                 (30)
```

This reconstructs the target point; `z` is then forced by `(21)` because
`f!=0`.  Equations `(25)--(30)` prove the exact image `(10)`.  All fibres are
finite.  Full differential rank makes `phi` unramified; quasi-finiteness and
miracle flatness between the smooth equidimensional surfaces make it flat.
Therefore `phi` is etale.

## 5. Maximal intersection and degree

From `(14)` and the standard degree of a polynomial function-field map,

```text
[k(x,y):k(F,W)]
 =[k(t):k(Q(t))]
 =deg Q=2d+1.                                                  (31)
```

Since `W=F^3P`, the fields `k(F,W)` and `k(F,P)` agree.  Also
`Frac(B)=k(F,W)` by `(21)`.

Now take

```text
G in k[x,y] intersect k(F,W).                                  (32)
```

The omitted locus `(10)` has codimension two, so the generic point of every
height-one prime `q` of the normal surface `S` lies in the etale image.  Choose
a height-one source prime `p` above it.  Etaleness gives ramification index
one, and therefore

```text
ord_p(G)=ord_q(G).                                              (33)
```

The left side is nonnegative because `G` is a source polynomial.  Thus `G`
has nonnegative order at every height-one prime of `B`.  The Krull
intersection of the height-one DVRs of the normal domain `B` gives `G in B`.
The reverse inclusion was established in Sections 1--2, proving `(11)`.

### 5.1 Units, Picard group, and the birational preliminary

The reduced divisor `f=0` is the disjoint-at-the-generic-point union of the
`d+1` arm primes

```text
L_0=(f,w),                 L_i=(f,w-lambda_i),
div(f)=L_0+L_1+...+L_d.                                      (33a)
```

After inverting `f`, equation `(21)` eliminates `z` and gives

```text
B_f=k[f,f^(-1),w],                                           (33b)
```

a UFD with unit group `k* f^Z`.  Nagata localization therefore says that
`Cl(B)` is generated by the arm classes.  If a principal divisor is supported
only on the arms, its defining rational function is a unit in `(33b)`, hence
is `c f^m`; its arm-order vector is `m(1,...,1)`.  Thus `(33a)` is the only
relation:

```text
B*=k*,
Cl(B)=Z^(d+1)/Z(1,...,1)=Z^d.                                (33c)
```

Smoothness identifies `Pic(S)=Cl(B)`.

Suppose now that `{A,C}=1` in `B` and `k(A,C)=Frac(B)`.  The symplectic
Poisson packet `(24)` makes `S -> A2_(A,C)` etale and birational, hence an
open immersion.  A divisorial complement would give a nonconstant unit on
`S`, contradicting `(33c)`.  A codimension-two complement has global ring
`k[A,C]` by normality; since `S` is affine, the open immersion must then be
all of `A2`, contradicting `Pic(S)=Z^d`.  Therefore a Darboux pair, if one
existed, would have target-field degree at least two.  The next argument
removes even that possibility.

### 5.2 The moving-arm Cech--de Rham obstruction

Factor the right side of `(21)` using the `d+1` polynomial root sections

```text
beta_0(f)=f^2,                 beta_i(f)=lambda_i (1<=i<=d),
(w-f^2)Psi(w)=product_(i=0)^d (w-beta_i(f)).                    (34)
```

For each arm, retain it and delete the others:

```text
U_i=S \ union_(j!=i) L_j,
a_i=[w-beta_i(f)]/f^3.                                        (35)
```

Put `Psi_i(f,w)=product_(j!=i)(w-beta_j(f))`.  On `U_i`, the two
presentations

```text
a_i=[w-beta_i(f)]/f^3=z/Psi_i(f,w)                             (36)
```

show that `a_i` is regular.  They give mutually inverse plane coordinates

```text
U_i ~= A2_(f,a_i),
w=beta_i(f)+f^3a_i,
z=a_i Psi_i(f,beta_i(f)+f^3a_i).                               (37)
```

The `U_i` cover `S`, and every intersection of at least two distinct charts
is

```text
D(f) ~= G_m,f x A1_w.                                         (38)
```

On `D(f)`, the inverse form of `(24)` is

```text
omega=df wedge dw/f^3.                                        (39)
```

Equation `(37)` extends it regularly over every chart and gives local
primitives

```text
omega=df wedge da_i=d eta_i,                 eta_i=-a_i df.    (40)
```

On an overlap,

```text
eta_j-eta_i=[beta_j(f)-beta_i(f)]f^(-3)df.                     (41)
```

The coefficient of `df/f` in `(41)` is

```text
rho_j-rho_i,                 rho_i=[f^2]beta_i(f),
(rho_0,rho_1,...,rho_d)=(1,0,...,0).                           (42)
```

This residue is the whole obstruction, not a heuristic local pole.  The
cover `(37),(38)` has the same Cech--de Rham rows as the arm-plane cover in
THM-3600:

```text
H^1_dR(U_i)=0,
H^1_dR(U_i intersection U_j)=k[df/f].                          (43)
```

The edge residues `(rho_j-rho_i)` satisfy the triple-overlap cocycle law.
The resulting degree-two class is the vertex vector `(rho_i)` modulo the
constant vector, so

```text
H^2_dR(S) ~= k^(d+1)/k(1,...,1),
[omega]=[(1,0,...,0)]!=0.                                    (44)
```

This explicitly distinguishes the global algebraic de Rham class from the
fact that `(40)` supplies a primitive on every individual plane chart.

Finally, if `{A,C}=1` in `B`, then the symplectic identity gives

```text
dA wedge dC=omega=d(A dC),                                    (45)
```

contradicting `(44)`.  Hence the completed surface has no Darboux pair in
any degree or support.

## 6. The degree-one hostile is the sharp boundary

Let instead `d=1`, so after scaling

```text
q(T)=a(T-alpha),                 alpha!=0.                      (46)
```

At its sole critical value `lambda=Q(alpha)`, one has exactly

```text
Q(T)-lambda=a^2(T-alpha)^3/3.                                  (47)
```

There is no residual factor in `(29)`.  Consequently the curve

```text
w=lambda,                 z=0,                 f!=0             (48)
```

on the putative completion has no source preimage: its only possible `t` is
`alpha`, where `q(t)=0`, incompatible with `f!=0`.  The image now misses a
divisor, so the height-one descent in Section 5 is invalid.  The missing
polynomial is visible explicitly:

```text
x^3=F^3/[3a(W-lambda)] in k(F,W).                              (49)
```

For `q=T+c`, retaining this Kummer coordinate and its first jet gives exactly
the cubic pseudo-plane maximal intersection of THM-3785, whose symplectic
form is exact.  Thus the transition `d=1` to `d>=2` is structural: in degree
one the putative surface still has the moving-arm class `(44)`, but it is not
the maximal intersection because it misses a divisor.  The required Kummer
filling changes the surface and kills that class.  In degree at least two,
the residual sheets in `(29)` repair codimension-one coverage while the
moving-arm residue survives and becomes a decisive global obstruction.

## 7. Exact quadratic control

For the first genuinely completed case

```text
q(T)=T^2+c,                 c!=0,                               (50)
Q(T)=T^5/5+2cT^3/3+c^2T,
Psi(W)=W^2+64c^5/225.                                         (51)
```

The two critical values are opposite and nonzero.  The surface and Poisson
packet specialize to

```text
F^3Y=(W-F^2)(W^2+64c^5/225),                                  (52)

{F,W}=F^3,
{F,Y}=3W^2-2F^2W+64c^5/225,
{W,Y}=3F^2Y+2F(W^2+64c^5/225).                                (53)
```

At a root `alpha^2=-c`, the boundary value is

```text
Y=64alpha^9/(675x^3).                                         (54)
```

The exact companion verifies polynomiality of `(6)`, `(52)--(54)`, the
abstract gradient signs, the quadratic and a cubic positive control, every
critical residual degree, the exact image packets, the moving-arm Cech
residue vector, and the degree-one hostile witness.  Normal and optimized
executions byte-match the frozen transcript.

## 8. Scope and failure boundaries

The theorem needs all of the following:

- squarefreeness of `q`, which gives the exact cubic contacts `(18)`;
- `q(0)!=0`, which separates the axis from the critical components;
- distinctness of `0,lambda_1,...,lambda_d`, which keeps all arms smooth and
  prevents two critical components from sharing a target packet; and
- `d>=2`, which supplies the residual good sheet `(29)`.

Repeated roots or colliding critical values require a confluent, higher-jet
Hermite completion and are **OPEN**.  For the distinct-critical family proved
here, however, `(44),(45)` close every Darboux support and target-field degree.
No planar Jacobian counterexample can arise from this completed family.
**QED, conditional only on independent hostile audit.**
