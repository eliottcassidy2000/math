---
id: THM-3797
title: "Confluent quadratic Hermite jet-completion no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The unique
  quadratic collision Q(0)=Q(1) in the third-order Hermite carrier has an
  explicit two-step polynomial jet completion.  Its naive hypersurface is
  singular and its normalization is still nonfinal; adjoining M and then N
  produces a smooth three-arm affine surface.  The source A2 maps to it
  quasi-finitely and etale, missing exactly two points, so this surface is the
  complete polynomial intersection with the rational Keller target field.
  Nevertheless its moving-root residue vector is (1,0,0) modulo constants,
  so its symplectic form is not exact and it has no polynomial Darboux pair.
source: jc_sparse_direct_search / confluent higher-pole Hermite lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The collision locus,
  residual noncritical roots, UFD boundary orders, quotient signs, exact
  boundary constants, normalization and unique A1 point, six-relation
  saturation, reduced three-arm fibre, chart inverses, Poisson signs,
  two-point image complement, etale valuation descent, degree-five field
  extension, and Cech--de Rham residue were independently rederived.  A
  separate saturation/Jacobian calculation confirms that B0[u] is normal
  with only the stated A1 singularity.  Normal and optimized executions
  byte-match the frozen transcript, all hashes match, and documentation
  checks pass.
related:
  - THM-3789-higher-pole-hermite-spectral-completion
  - THM-3791-moving-root-danielewski-resonant-jet-de-rham-law
script: 04-computation/jc2_confluent_quadratic_hermite_scout.py
output: 05-knowledge/results/jc2_confluent_quadratic_hermite_scout.out
script_sha256: 36abf389045ac12baa684f1837b35638d90ca54ff17c616a9fe228790641cc83
output_sha256: 91e91ae33d009dd4ab74844cff7721e88f6eb99fbea8b52e2085b90635676f4f
semantic_sha256: f2c7f95b532e3be227a1bef021a14376d8c47de447de9090658d63c2c893d6b8
hash_basis: raw LF bytes
---

# THM-3797 -- the confluent quadratic carrier needs two jets and still has no Darboux pair

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
closes the first ramified spectral boundary left outside THM-3789 and
THM-3791.  The collision does create new polynomial observables, and stopping
at the normalization loses the decisive geometry.  After the complete
two-step affine modification, however, the same moving-axis residue survives.

Let `k` be an algebraically closed field of characteristic zero.  Choose

```text
b=(5+sqrt(-15))/20  or  (5-sqrt(-15))/20,
10b^2-5b+1=0,                                                (1)
```

and put

```text
t=x^2(1+xy),
q(T)=(T-1)(T-b),
Q(T)=integral_0^T q(s)^2 ds
    =T(T-1)^3(2T+1-5b)/10.                                  (2)
```

Then

```text
Q(0)=Q(1)=0,
ell=Q(b)=b(b-1)^3(1-3b)/10 !=0.                             (3)
```

Define the rational Keller seed

```text
F=xq(t),                 W=Q(t),                 P=W/F^3.    (4)
```

The following two quotients and the singular-completion coordinate are source
polynomials:

```text
M=W(W-ell)/F^2,
N=M(M+ell)/F,
Z=(W-F^2)M/F.                                                (5)
```

Let

```text
R=k[F,W,M,N,Z] subset k[x,y],                 S=Spec R.       (6)
```

Then `S` is smooth, and the source morphism

```text
phi:A2_(x,y) -> S                                               (7)
```

is quasi-finite and etale.  Its image is the complement of exactly two closed
points.  Moreover

```text
k[x,y] intersection k(F,P)=R,
[k(x,y):k(F,P)]=5.                                           (8)
```

The inverse symplectic form on `S` has class

```text
[omega]=[(1,0,0)] !=0 in H^2_dR(S)=k^3/k(1,1,1).             (9)
```

Consequently no `A,C in R` satisfy `{A,C}=1`.  In particular no polynomial
Keller pair obtained inside the rational target field `k(F,P)` can turn this
confluent carrier into a planar Jacobian counterexample.

## 1. The collision locus and the rational Keller seed

The condition `Q(1)=0` before imposing `(1)` is

```text
Q(1)=(10b^2-5b+1)/30.                                       (10)
```

Thus `(1)` is the complete normalized quadratic collision locus.  Its two
points are conjugate, `b` is different from `0,1`, and `ell` is nonzero.  The
other critical fibre factors as

```text
Q(T)-ell=(T-b)^3 R_b(T),                                     (11)
R_b(T)=T^2/5+(b-5)T/10-3(5b-11)/100.
```

Modulo `(1)`,

```text
Res_T(R_b,q)=-(35b-13)/2000 !=0.                             (12)
```

Thus the residual quadratic over `W=ell` contains no critical root.  Over
`W=0`, besides the triple root `T=1`, the simple roots `T=0` and
`T=(5b-1)/2` are also noncritical.

Exactly as in the higher-pole seed calculation,

```text
{F,t}=x^3q(t),                 {F,W}=F^3,                     (13)
```

and therefore

```text
{F,P}=1.                                                        (14)
```

Also

```text
k(x,y)=k(F,t),
x=F/q(t),                 y=(t-x^2)/x^3.                     (15)
```

## 2. The complete quotient ladder

The reduced components of `F=0` are

```text
A=(x),                 D_1=(t-1),                 D_b=(t-b). (16)
```

They are pairwise coprime primes.  The relevant orders are

```text
                 A       D_1       D_b
ord(F)           1        1         1
ord(W)           2        3         0
ord(W-ell)       0        0         3.                       (17)
```

Consequently `M` in `(5)` is polynomial.  More explicitly, `(11)` gives

```text
M=(1+xy)q(t)(2t+1-5b)R_b(t)/10.                              (18)
```

Its exact boundary values are

```text
M|A=-ell,                 M|D_1=0,                 M|D_b=0.  (19)
```

Thus `M+ell` vanishes on `A`, while `M` vanishes on `D_1,D_b`.
The UFD factorization `(16)` proves that `F` divides `M(M+ell)`, so `N` is
polynomial.  Similarly `W-F^2` has sufficient axis order and `M` vanishes on
both critical components, proving that `Z` is polynomial.

There is a useful completeness statement already at this stage.  The three
boundary points in the `(w,m)`-plane are

```text
(0,-ell),                 (0,0),                 (ell,0).    (20)
```

Their reduced vanishing ideal is

```text
(w(w-ell), wm, m(m+ell)).                                   (21)
```

After division by `f`, the first two generators give `fM` and

```text
wM/f=Z+fM,                                                    (22)
```

while the third gives `N`.  Hence `N` is the only genuinely new first jet
after adjoining `M`; it is not an arbitrary guessed quotient.

## 3. Why normalization alone is not the answer

If one retains only the old Hermite coordinate `Z`, the naive completion is

```text
B_0=k[f,w,z]/(f^3z-(w-f^2)w(w-ell)).                          (23)
```

Its singular locus is the whole line

```text
Sing(Spec B_0)=V(f,w).                                        (24)
```

Put

```text
u=w(w-ell)/f.                                                  (25)
```

This element is integral over `B_0`: it satisfies

```text
fu=w(w-ell),
f^2(z+u)=wu,
u^2=f(w-ell)(z+u).                                            (26)
```

The last equation is monic in `u`.  The finite birational ring `B_0[u]` is
regular away from

```text
p=(f,w,u,z)=(0,0,0,0).                                       (27)
```

At `p`, `w-ell` is a unit and the first equation eliminates `w`; with
`v=(w-ell)(z+u)`, the completed local equation is

```text
u^2=fv.                                                        (28)
```

This is the normal `A_1` surface singularity.  Thus `B_0[u]` is exactly the
normalization of `B_0`, but it is not the final smooth completion.

Since `u=fM`, adjoining `M=u/f` is the next affine modification.  The ring
`B_0[M]` is smooth.  Indeed its two special-fibre lines are the collision line
`(f,w,z)=(0,0,0)` and the critical line `(f,w,m)=(0,ell,0)`; the saturated
relations have respectively the constant Jacobian minors `-ell^2` and
`ell^2` there, while `f!=0` is a plane.  But the source map is still not
quasi-finite: the entire axis maps to `(f,w,m,z)=(0,0,-ell,0)`, and the entire
component `D_1` maps to `(0,0,0,0)` on its exceptional line.  The second jet
`N` in `(5)` is exactly what restores the missing tangent coordinate on both
source components.

## 4. Exact presentation and the three smooth charts

Use lower-case letters for the five abstract generators in `(6)`.  The exact
kernel is

```text
I=(
 f^2m-w(w-ell),
 fn-m(m+ell),
 fz-(w-f^2)m,
 (m+ell)(fm+z)-nw,
 (w-ell)(fm+z)-fm^2,
 nw(w-ell)-fm^2(m+ell)
).                                                            (29)
```

Equivalently,

```text
I=(f^2m-w(w-ell), fn-m(m+ell), fz-(w-f^2)m):f^infinity.       (30)
```

The exact companion computes both containments in `(29),(30)`, not just a
one-way list of consequences.  In particular

```text
R=k[f,w,m,n,z]/I,                 R_f=k[f,f^(-1),w].          (31)
```

The special fibre is exactly the disjoint union at its generic points of
three affine lines

```text
L_A: (f,w,m,z)=(0,0,-ell,0),                 coordinate n,
L_1: (f,w,m,z)=(0,0,0,0),                    coordinate n,
L_b: (f,w,m,n-z)=(0,ell,0,0),                coordinate z.   (32)
```

There are already three independent relation gradients on every line in
`(32)`: using the variable triples `(w,m,z)`, `(w,m,z)`, and `(w,m,n)`, one
gets triangular minors with diagonal entries respectively

```text
(ell,ell,-ell),                 (ell,-ell,ell),
(-ell,-ell,-ell).                                               (32a)
```

Off `f=0`, equation `(31)` is a plane.  Hence the surface defined by `(29)`
is smooth and normal before any chartwise divisibility is invoked.

The following formulas give a complete plane atlas.  In each display `a` is
free and all omitted variables are the displayed polynomials in `(f,a)`.

For the moving-axis arm, put `s=1+fa`:

```text
w=f^2s,
m=s(f^2s-ell),
n=m(fs^2-ell*a),
z=f^2 a m.                                                     (33)
```

For the root `1` arm, put `r=f^3a-ell`:

```text
w=f^3a,
m=f a r,
n=a r(ell+f a r),
z=f^2 a r(fa-1).                                              (34)
```

For the root `b` arm, put `r=ell+f^3a`:

```text
w=r,
m=f a r,
n=a r(ell+f a r),
z=a r(ell+f^3a-f^2).                                         (35)
```

Substitution into all six relations `(29)` gives zero.  Conversely, when
`f!=0` the inverse coordinates are

```text
a_A=(w-f^2)/f^3,                 a_1=w/f^3,
a_b=(w-ell)/f^3.                                                (36)
```

At `f=0`, formulas `(33)--(35)` restrict respectively to

```text
n=ell^2 a,                 n=-ell^2 a,
n=z=ell^2 a.                                                  (37)
```

On the open retaining `L_A`, the relations show successively that `w` is
divisible by `f^2`, `m+ell` by `f`, and `w-f^2` by `f^3`; on the open
retaining `L_1`, they show `m` is divisible by `f` and `w` by `f^3`; on the
open retaining `L_b`, they show `m` is divisible by `f` and `w-ell` by
`f^3`.  Normality therefore makes the three rational functions `(36)`
regular on their stated opens.  Formulas `(33)--(37)` are their polynomial
inverses.  Thus each chart is `A2_(f,a)`, retains exactly its named arm, and
the three charts cover `S`.  This also rules out hidden components in the
saturation computation.  Every pairwise or triple overlap is `D(f)`.

## 5. Poisson packet, exact image, and etaleness

On `f!=0`, transport the bracket from `(13)` by `{f,w}=f^3`.  Differentiating
the rational definitions of `m,n,z` gives

```text
{f,m}=f(2w-ell),
{f,n}=(2w-ell)(2m+ell),
{f,z}=ell*f^2-2ell*w-2f^2w+3w^2.                             (38)
```

These formulas are polynomial and therefore hold on all of `S`.  On the
three arms in `(32)`, the appropriate local minors are

```text
{f,n}|L_A=ell^2,
{f,n}|L_1=-ell^2,
{f,z}|L_b=ell^2.                                              (39)
```

Hence the bracket is everywhere symplectic and its inverse is

```text
omega=df wedge dw/f^3                                         (40)
```

on the dense chart `f!=0`.

The exact source-boundary packets are

```text
x=0:
 (F,W,M,Z,N)=(0,0,-ell,0,(ell^2/b)y),

t=1:
 (F,W,M,Z,N)=(0,0,0,0,kappa/x^3),

t=b:
 (F,W,M,Z,N)=(0,ell,0,kappa/x^3,kappa/x^3),                  (41)
```

where

```text
kappa=ell^2/[3(b-1)] !=0.                                    (42)
```

Thus the axis fills `L_A`, while `D_1` and `D_b` fill the punctured lines
`L_1` and `L_b`.  The only missing points are

```text
p_1=(0,0,0,0,0),                 p_b=(0,ell,0,0,0),           (43)
```

in the coordinate order `(f,w,m,n,z)`.

No point with `f!=0` is lost.  If `w` is not a critical value, every root of
`Q(T)-w` is noncritical.  For `w=0`, the root `T=0` has `q(0)=b!=0`.  For
`w=ell`, `(11),(12)` supply a residual root with `q(T)!=0`.  Given such a
root `tau`, reconstruct

```text
x=f/q(tau),                 y=(tau-x^2)/x^3.                  (44)
```

Equations `(41)--(44)` prove

```text
image(phi)=S\{p_1,p_b}.                                       (45)
```

All fibres are finite.  Off `F=0`, the minor `{F,W}=F^3` is nonzero; on the
three boundary components the pullbacks of the units `(39)` are nonzero.
Thus `phi` has full differential rank everywhere.  Quasi-finiteness and
miracle flatness between the smooth equidimensional surfaces make `phi`
etale.

## 6. The maximal polynomial intersection

The field equalities `(15)` give

```text
Frac(R)=k(F,W)=k(F,P),
[k(x,y):k(F,W)]=deg Q=5.                                    (46)
```

Take `G in k[x,y] intersection k(F,W)`.  The complement `(43)` has
codimension two, so the generic point of every height-one prime of the normal
surface `S` lies in the etale image.  Choose a source height-one prime above
it.  Etaleness gives ramification index one, hence the target valuation of
`G` equals a nonnegative source valuation.  The intersection of the
height-one DVRs of `R` now gives `G in R`.  The reverse inclusion follows
from Section 2, proving `(8)`.

## 7. The residue that survives confluence

The three plane coordinates in `(36)` have root jets

```text
beta_A(f)=f^2,                 beta_1(f)=0,
beta_b(f)=ell.                                                   (47)
```

On each chart,

```text
omega=df wedge da_i=d(-a_i df).                               (48)
```

On an overlap,

```text
-a_j df+a_i df=[beta_j(f)-beta_i(f)]f^(-3)df.                 (49)
```

The only logarithmic coefficient is the coefficient of `f^2` in `(47)`:

```text
rho=([f^2]beta_A,[f^2]beta_1,[f^2]beta_b)=(1,0,0).            (50)
```

The charts are affine planes and every multiple intersection is
`G_m,f x A1_w`.  The same explicit truncated-simplex calculation as in
THM-3791 therefore gives

```text
H^1_dR(S)=0,
H^2_dR(S)=k^3/k(1,1,1),
[omega]=[(1,0,0)]!=0.                                        (51)
```

This is a direct calculation on the final smooth model; the squarefree
special-fibre hypothesis of THM-3791 is not being silently applied to the
singular hypersurface `(23)`.

Finally, if `{A,C}=1` in `R`, then

```text
dA wedge dC=omega=d(A dC),                                   (52)
```

contradicting `(51)`.  This proves the all-support no-go.

## 8. Scope and the structural lesson

The two conjugate values in `(1)` are both covered.  This theorem does not
claim that every collision pattern in every degree has the same affine
modification tree.  It proves the complete quadratic `Q(0)=Q(1)` boundary.

Three stopping points must not be conflated:

1. the naive one-equation completion is nonnormal along a line;
2. its normalization is normal but retains an `A_1` point;
3. adjoining `M` is smooth but still collapses source divisors.

Only after the second jet `N` does one recover the smooth quasi-finite etale
atlas needed for the exact intersection.  Confluence therefore creates real
new polynomial structure, but it does not erase the moving-axis coefficient:
the obstruction moves from a squarefree hypersurface to an iterated
jet-completion tree.  **QED.**
