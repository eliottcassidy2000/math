---
id: THM-3973
title: "Exact-volume simple-cubic determinantal affine-plane completions"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT. For every
  n>=2, the determinantal ring generated inside k[x,t] by
  z=1+x^n t, p=zt, and y=x^(n-1)zt^2 is a smooth normal affine surface X_n
  containing A2_(x,t) with one A1 boundary D, scalar units, Cl(X_n)=Z[D],
  and div(dx wedge dt)=D. Unlike THM-3971, the source volume is globally
  exact: beta=-dz/((n-1)x^(n-1)) is regular and d beta=dx wedge dt. Thus
  X_n passes the boundary, class, simple-cubic canonical, and exact-volume
  completion invoices simultaneously. Every higher top-generator variant
  collapses to the same B_n. The constant has bracket length at most two,
  and an exact rational one-bracket compression has only the divisors x and
  2z-1 in its denominators. Homogeneous pairs, each canonical generator,
  every generator-linear pair, and 7033 sparse low-filtration rows have no
  polynomial mate. A polynomial Darboux pair and finite cubic target map
  remain OPEN; no Jacobian counterexample is claimed.
source: jc-zero-debt-lift + jc-degree6-one-place + root / post-THM-3971 residue-cancellation design, 2026-08-24
depends_on:
  - THM-3922-affine-plane-open-boundary-basis-class-group-obstruction
related:
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
  - THM-3968-canonical-vector-different-affine-plane-boundary-obstruction
  - THM-3971-canonical-debt-determinantal-affine-plane-completion
  - THM-3972-simple-collision-affine-p-graph-blowup-normalization
script: 04-computation/jc2_exact_volume_simple_cubic_completion_thm3973.py
output: 05-knowledge/results/jc2_exact_volume_simple_cubic_completion_thm3973.out
script_sha256: c0b9f024a33e629a49d7f816df4dd26f03cc93ef36a8fab67c8bca4febd9269b
output_sha256: 8c0d2e6ab2d3a2695c0a6ecb184c31ed9aab22c830d89974ff2c00282c7670d5
semantic_sha256: 6674f41354cb93ced1324905338268d55dba40242bc789fb889aeeca099bda9b
hash_basis: raw LF bytes
---

# THM-3973 -- a simple-cubic completion passes the exact-volume gate

**PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.** Work over
an algebraically closed field `k` of characteristic zero. Fix `n>=2` and,
inside `k[x,t]`, put

```text
z=1+x^n t,             p=zt,             y=x^(n-1)zt^2,
B_n=k[x,z,p,y],         X_n=Spec(B_n).                     (1)
```

Then `X_n` is a smooth normal affine surface with

```text
U=D(x) union D(z) isomorphic to A2_(x,t),
D=X_n minus U=V(x,z) isomorphic to A1_y,                  (2)
B_n^*=k^*,                Cl(X_n)=Z[D].                    (3)
```

The source volume has the simple-cubic canonical divisor but zero de Rham
class:

```text
eta=dx wedge dt,          div_Xn(eta)=D,                  (4)
beta=-dz/((n-1)x^(n-1)) in Gamma(X_n,Omega^1),
d beta=eta.                                                   (5)
```

Thus its two completion coordinates are

```text
kappa=ord_D(eta)=1,              rho=[eta]=0.             (6)
```

This is the first family in the present search that simultaneously passes
the affine-plane boundary, primitive class-group basis, simple `(2,1)`
cubic canonical-different, and exact-volume gates. It is a positive
completion **passport**, not yet a counterexample: no polynomial pair
`A,C in B_n` with `J_(x,t)(A,C)=1`, and no finite target map, is constructed.

## 1. Exact determinantal presentation

Let

```text
M_n=[[z,p,x],[x^(n-1)p,y,z-1]].                           (7)
```

Then

```text
B_n isomorphic to k[x,z,p,y]/I_2(M_n),                    (8)
```

where the three minors are

```text
zy-x^(n-1)p^2,
z(z-1)-x^n p,
p(z-1)-xy.                                                (9)
```

The substitutions `(1)` annihilate `(9)`, so `(8)` has a natural
surjection onto `B_n`. To prove that no extra component or saturation is
hidden, localize first at `z`. Setting `t=p/z`, equations `(9)` become

```text
z=1+x^n t,             y=x^(n-1)zt^2,                    (10)
```

and the localized quotient is

```text
k[x,t,(1+x^n t)^-1].                                     (11)
```

On `D(z-1)`, eliminate `p=xy/(z-1)`. The remaining equation is the domain

```text
k[x,z,y,(z-1)^-1]/(z(z-1)^2-x^(n+1)y).                  (12)
```

The first minor follows from the hypersurface equation in `(12)`. Since
`D(z)` and `D(z-1)` cover the whole quotient, `(11)--(12)` prove the global
isomorphism `(8)`, including primality of the determinantal ideal.

The fraction field is `k(x,t)`, since `t=(z-1)/x^n`, so `X_n` is an
integral finite-type surface.

## 2. The affine-plane open and smooth boundary

The two localizations of `B_n` are

```text
B_n[x^-1]=k[x,x^-1,t],
B_n[z^-1]=k[x,t,z^-1].                                   (13)
```

They glue to `A2_(x,t)`: the two polynomials `x` and `z=1+x^nt` have no
common source zero. Their union is therefore exactly the open `U` in `(2)`.

Modulo `(x,z)`, the last relation in `(9)` gives `p=0`, while `y` remains
free. Hence

```text
B_n/(x,z)=k[y],                                           (14)
```

so `D` is one reduced affine line and there is no second divisorial
boundary component.

The open `U` is smooth. Near `D`, use the hypersurface chart `(12)` and put

```text
F=z(z-1)^2-x^(n+1)y.
```

Then

```text
F_z=(z-1)(3z-1),              F_z|D=1.                   (15)
```

Thus `X_n` is smooth along `D`, hence smooth and normal everywhere.
Applying THM-3922 to `(2)` proves `(3)`: the sole boundary class is a
primitive generator, not merely a nonzero class.

## 3. Simple canonical debt and exact source volume

Differentiate `y=x^(n-1)zt^2` at fixed `x`. Using `x^nt=z-1` gives

```text
partial y/partial t
 =x^(n-1)t(3z-1)
 =x^-1(z-1)(3z-1).                                       (16)
```

Therefore near `D`

```text
eta=dx wedge dt
 =x/((z-1)(3z-1)) dx wedge dy.                           (17)
```

The denominator is a unit on `D`, while `eta` is nowhere zero on `U`.
Consequently `(17)` proves the exact divisor row `(4)`. In the basis `(3)`,
`[K_Xn]=[D]`, precisely the canonical-different coordinate required by
THM-3968 for one simple cubic ramification prime.

Unlike the height-one-pole family in THM-3971, this volume has no logarithmic
residue. On `U`, define

```text
beta=-t dx-(1/(n-1))d(xt)
    =-dz/((n-1)x^(n-1)).                                  (18)
```

It is visibly regular in the source coordinates on `U`. Differentiating the
boundary equation gives

```text
(z-1)(3z-1)dz=(n+1)x^n y dx+x^(n+1)dy.                   (19)
```

After division by `x^(n-1)`, the right side of `(19)` still has factors
`x` and `x^2`. Thus the last expression in `(18)` is regular near `D` as
well, and the two descriptions glue to a global regular one-form on `X_n`.
Finally,

```text
d beta=x^-n dx wedge dz=dx wedge dt=eta.                 (20)
```

This proves `(5)--(6)`.

The Gysin row still gives `H^2_dR(X_n)=k`, because `U=A2` and `D=A1`.
Equation `(20)` says specifically that the Keller source volume is the zero
class, not that the whole cohomology group vanished. This is the exact
route around THM-3971's residue-one obstruction.

## 4. All higher top-generator variants are the same ring

The construction initially appears to have a second height parameter. For
`m>=2`, put

```text
y_m=x^(n(m-1)-1)zt^m,
B_(n,m)=k[x,z,p,y_m].                                     (21)
```

Let `u=z-1=x^nt`. In the fraction field,

```text
p=x^-n u(u+1),
y_m=x^(-n-1)u^m(u+1),                                    (22)
y_m=u^(m-2)y_2,
x^(n-1)p^2=(u+1)y_2.                                     (23)
```

Since `gcd(u+1,u^(m-2))=1`, choose `A_m,B_m in k[u]` with

```text
A_m(u)(u+1)+B_m(u)u^(m-2)=1.                             (24)
```

Multiplication by `y_2` and `(23)` give

```text
y_2=A_m(u)x^(n-1)p^2+B_m(u)y_m in B_(n,m).               (25)
```

The reverse inclusion follows from `y_m=u^(m-2)y_2`. Hence

```text
B_(n,m)=B_(n,2)=B_n                 for every m>=2.       (26)
```

The apparent height direction is pure generator redundancy. The true
positive family has the single pole-order parameter `n>=2`.

## 5. Exact graded pieces

Give the source ring the grading

```text
wt(x)=1,                wt(t)=-n,              wt(u)=0.   (27)
```

The generators in `(1)` are homogeneous: `p` has weight `-n`, `y` has
weight `-(n+1)`, and `z=1+u` has weight zero. For `r>=0`, one has

```text
(B_n)_r=x^r k[u].                                           (28)
```

For `k>=1`, the negative component is exactly

```text
(B_n)_(-k)
 =x^-k u^ceil(k/n)(u+1)^ceil(k/(n+1)) k[u].               (29)
```

To prove `(28)--(29)`, write a generator monomial as

```text
x^c p^a y^b f(u).
```

If its weight is `-k`, then

```text
c=na+(n+1)b-k>=0                                          (30)
```

and its coefficient after extracting `x^-k` is

```text
u^(a+2b)(u+1)^(a+b)f(u).                                  (31)
```

The minimum possible `u` exponent in `(31)` is `ceil(k/n)`, attained with
`b=0`; the minimum possible `u+1` exponent is `ceil(k/(n+1))`, attained
with `a=0`. Since `u` and `u+1` are coprime, the ideal generated by all
rows `(31)` is their product with exactly those two exponents. This proves
`(29)`. The nonnegative argument is the same extraction without a
denominator and proves `(28)`.

This module formula is the main finite support sidecar for future searches.
It records both mandatory colors `u=0` and `u=-1` in every negative-weight
coefficient.

## 6. Homogeneous and canonical-generator nonentry

The bracket has weight shift `n-1`. More precisely, for homogeneous rows

```text
A=x^r f(u),                 C=x^s g(u),                   (32)
```

direct differentiation gives

```text
J(A,C)=x^(r+s+n-1)(r f g'-s f'g).                         (33)
```

Suppose `(33)` were one. Then `r+s=1-n`. If both weights were negative,
`(29)` would make both `f` and `g` divisible by `u(u+1)`, and the last
factor in `(33)` would still be divisible by `u(u+1)`. Both weights cannot
be nonnegative. Thus, after swapping the rows, take `r>=0` and
`s=-(r+n-1)<0`.

If `r=0`, the expression is a multiple of `g`, hence not one. If `r>0`,
evaluate at both `u=0` and `u=-1`. For a nonzero constant, `g` must have a
simple zero at each address. Formula `(29)` forces

```text
ceil((r+n-1)/n)=ceil((r+n-1)/(n+1))=1,                   (34)
```

so `r=1` and `s=-n`. Write `g=u(u+1)G`. The leading term of

```text
f g'+n f'g                                               (35)
```

has coefficient `deg(g)+n deg(f)`, which is nonzero in characteristic zero,
and degree at least one. It cannot be the constant one. Therefore:

```text
no two homogeneous elements of B_n form a Darboux pair.  (36)
```

The four canonical generators also have no mate even against an arbitrary
inhomogeneous polynomial. Weight separation reduces the only potentially
constant component to the following exact rows.

1. A mate of `x` would satisfy `C_t=1`, hence
   `C=t+f(x)`. But `ord_D(t)=-n` while `f(x)` is regular, so this cannot lie
   in `B_n`.
2. For every `A in B_n`,

   ```text
   J(A,z)=x^(n-1)(xA_x-ntA_t),                            (37)
   ```

   which cannot equal one.
3. Only a weight-one component `A=xF(u)` can mate `p`; it would have to solve

   ```text
   n u(1+u)F'+(1+2u)F=1.                                 (38)
   ```

   For nonconstant `F` of degree `d`, the leading coefficient is multiplied
   by `nd+2`; a constant `F` also fails by its `u` coefficient.
4. Only `A=x^2F(u)` can mate `y`; its bracket is

   ```text
   u[(n+1)u(1+u)F'+2(2+3u)F],                            (39)
   ```

   which is divisible by `u`.

The cheapest boundary-address correction does not help. For
`lambda in k`,

```text
y+lambda p=p(x^(n-1)t+lambda).                            (39a)
```

When `lambda=0`, row 4 applies. When `lambda!=0`, the two factors in
`(39a)` vanish simultaneously at

```text
x=lambda^-1,             t=-lambda^n,             z=0.   (39b)
```

Hence both partial derivatives of `y+lambda p` vanish there. It cannot be
one member of a Darboux pair at any height `n`.

These are all-degree controls on the distinguished generators, not a
classification of arbitrary multigraded pairs.

## 7. Exact bracket length two and rational compression

Put

```text
g=u(u+1),                 w=1+2u=2z-1.                   (40)
```

The elementary Bezout identity

```text
w^2-4g=1                                                     (41)
```

turns the regular primitive into a two-bracket certificate. Since
`p=x^-n g`, differentiation gives

```text
x dp=-n p dx+x^(1-n)w dz.                                (42)
```

Using `(41)--(42)` in `(18)`, then subtracting the exact differential
`d(-wx p/(n-1))`, yields

```text
beta congruent -wp dx+(6/(n-1))xp dz       modulo exact forms. (43)
```

Taking differentials gives the exact polynomial identity

```text
1=J(-wp,x)+J((6/(n-1))xp,z).                             (44)
```

In fact the two summands are `1+6g` and `-6g`. Thus the constant has
Poisson-bracket length at most two on every `B_n`; the positive problem is
the compression of these two brackets to one.

That compression already exists rationally. Define

```text
P=g/((n-1)w^2),                 Q=w^3/x^(n-1).            (45)
```

Then

```text
dQ=(n-1)w^2x^-n[-w dx+(6/(n-1))x dz],
P dQ=-wp dx+(6/(n-1))xp dz,                               (46)
J(P,Q)=1.                                                  (47)
```

The exact remaining debts are visible: `P` has a pole on `w=0`, where
`g=-1/4`, and `Q` has a pole of order `n-1` along `x=0`. Neither belongs to
`B_n`. Any successful polynomial compression must absorb these two
divisors without restoring a boundary unit or a critical curve.

## 8. Minimal bounded search and finiteness protocol

There is an exact all-parameter submersion obstruction before the bounded
mate search. On `B_2`, every generator-linear function whose restriction to
`D=A1_y` is nonconstant can be normalized, up to a scalar and additive
constant, to

```text
C=y+A x+B z+C_0 p.                                       (47a)
```

Direct elimination of `t` from `C_x=C_t=0` gives

```text
Res_t(C_x,C_t)=3x^3 H(x),                                (47b)
```

where, from degree seven down to degree zero, the coefficients of `H` are

```text
3B^3,
9A^2-6ABC_0,
-B^2C_0,
2AC_0^2+3B^2,
4AC_0-3BC_0^2,
-4A+2BC_0,
C_0^3,
-C_0^2.                                                   (47c)
```

If `C_0!=0`, then `H` is nonconstant and `H(0)=-C_0^2!=0`, so it has a
nonzero root. If `C_0=0` and `A!=0`, divide the remaining `x^2` and use the
nonzero constant `-4A`; if `A=0,B!=0`, the remaining factor contains
`Bx^3+1`. The endpoint `A=B=C_0=0` is `C=y`, already killed in Section 6.
At every nonzero root, the leading `t` coefficients of `C_x,C_t` are
`3x^2,3x^3`, so the resultant zero is a genuine affine common root rather
than a root at infinity. Thus every row `(47a)` has a critical point and
cannot occur in a Darboux pair.

For the minimal member `n=2`, give `x,z,p,y` filtration one and let `F_N`
be the span of their monomials of total generator degree at most `N`. Exact
integer row reduction gives

```text
dim(F_0),...,dim(F_4)=1,5,14,28,47.                       (48)
```

First take arbitrary generator-linear rows

```text
A=a_1x+a_2z+a_3p+a_4y,
C=c_1x+c_2z+c_3p+c_4y.                                   (49)
```

Equating every coefficient of `J(A,C)-1` gives a reduced Groebner basis
`[1]`. Hence no pair in `(49)` exists over any characteristic-zero field.

For the larger linear-mate census, use the deterministic nonconstant basis

```text
y,p,z,x,y^2,py,p^2,zy,zp,z^2,xp,xz,x^2                  (50)
```

of `F_2`. For each basis row and every signed sum with support two, three,
or four, normalize the first sign to `+1` and solve the full exact linear
system

```text
A in F_2,                  C in F_4,        J(A,C)=1.     (51)
```

There are

```text
13+2*C(13,2)+4*C(13,3)+8*C(13,4)=7033                    (52)
```

tested `A` rows and no survivor. The calculation uses exact integer ranks;
there is no finite-field inference. It is only a hostile cutoff and is not
used to infer the open all-degree claim.

If a future pair survives, it becomes a planar counterexample only after
the usual two checks. First, it is already nonautomorphic if its entries lie
in `B_n`: since `ord_D(t)=-n`, one has `t notin B_n`, whereas an automorphic
pair would force `k[A,C]=k[x,t] subset B_n`. Second, to realize `X_n` as the
actual finite normal completion one must additionally prove

```text
B_n is finite over k[A,C].                                (53)
```

For a concrete survivor, `(53)` is an exact elimination/integrality test on
the four generators. The present census has no survivor, so no vacuous
finiteness claim is made.

Reproduce with

```bash
python3 04-computation/jc2_exact_volume_simple_cubic_completion_thm3973.py
python3 -O 04-computation/jc2_exact_volume_simple_cubic_completion_thm3973.py
```

Both runs must print `CHECKS=7415` and byte-match the frozen output.

## 9. Exact scope and next design gate

THM-3973 does **not** construct a Keller pair. It proves that the coarse
completion invoices which killed the preceding candidates are mutually
compatible on an explicit non-complete-intersection determinantal surface:

```text
one rational unibranch A1 boundary,
primitive free boundary class,
simple-cubic canonical vector,
exact source volume,
constant bracket length at most two.                       (54)
```

The sibling graph blowup in THM-3972 has a related abstract completion
passport but retains a ramification obstruction for its named cubic graph
map. No abstract isomorphism between that surface and `X_n` is asserted;
their marked pole/key data differ and require a separate comparison.

The next positive gate is now sharply localized: either compress `(44)`
polynomially, or deform the rational pair `(45)` so that the `x` and `w`
denominators are absorbed while `(2)--(6)` survive. Arbitrary multigraded
pairs and the finiteness row `(53)` remain **OPEN**, as does `JC(2)`.
**QED.**
