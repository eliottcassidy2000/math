---
id: THM-2127
title: "Full resonant trains, exact coprime two-face closure, and affine-root residues"
status: >
  PROVED. THM-2102's proper-power first-defect equation has a full resonant
  solution: each resonance seeds a new fractional-binomial train, but every
  later train has strictly higher h-adic order than the original one. Hence an
  exact two-weighted-face component f=h^m+A with h power-free and gcd(h,A)=1
  has a Jacobian mate only in the terminal case deg(h)+deg(A)=w_1+w_2; there
  h,A are polynomial coordinates and the pair is an automorphism. Separately,
  in polynomial coordinates (u,v), every f=u+P(A(u)v+B(u)) with deg P>=2 has
  a rational, hence polynomial, mate exactly when A is constant. Its negative
  direction is an all-order moving-pole residue obstruction. These close two
  proper-power subclasses but do not prove planar JC.
source: codex-2026-07-22-JC2-proper-power-all-order
depends_on:
  - THM-2102
related:
  - THM-2045
  - THM-2071
  - HYP-8950
---

# THM-2127 -- proper-power all-order obstructions

Write

```text
{F,G}=F_x G_y-F_y G_x.                               (1)
```

The theorem has two complementary parts. The first solves the complete
resonant recursion for a two-face component. The second changes to an exact
approximate-root coordinate and detects a nonzero residue after **all** face
orders have been assembled.

## 1. The full resonant fractional-binomial ladder

Let `w=(w_1,w_2)` be a positive integral weight, put

```text
W=w_1+w_2,                                            (2)
```

and let `h in C[x,y]` be a nonconstant, power-free, `w`-homogeneous
polynomial of weighted degree `d`. Power-free means that the gcd of the
multiplicities in the irreducible factorization of `h` is one.

Suppose `f,g in C[x,y]` satisfy

```text
{f,g}=1.                                              (3)
```

Assume their leading faces have a common primitive root,

```text
in_w(f)=h^m,                 in_w(g)=c h^n,            (4)
```

where `m>=2`, `n>=1`, and `c!=0`. Let `rho>0`, let `K>=1`, and assume that
through defect `K rho` the only lower face of `f` is

```text
A!=0,              deg_w(A)=a=m d-rho>0.             (5)
```

Precisely, every homogeneous component of

```text
f-h^m-A                                                (6)
```

has defect strictly greater than `K rho` from the top degree `m d`.

Put

```text
D=(m+n)d-W.                                           (7)
```

Assume only

```text
K rho<D.                                                (8)
```

For `j=0,...,K`, let `G_j` be the homogeneous component of `g` of defect
`j rho` from `n d`, with `G_j=0` if that degree is absent. Thus

```text
G_0=c h^n.                                            (9)
```

Put

```text
s=gcd(d,rho),              L=d/s,              e=rho/s. (10)
```

Then there are scalars `lambda_b`, with `lambda_0=c`, such that every face in
this defect chain has the exact rational-function expansion

```text
G_j=sum_(0<=b<=floor(j/L))
 lambda_b binom((n-b e)/m,j-b L)
 h^[n-b e-m(j-b L)] A^[j-b L],             0<=j<=K.  (11)
```

The scalar `lambda_b` first appears at the resonant defect `j=bL`. Equation
(11) lives in `C(x,y)`, but its left side is polynomial; polynomiality of the
whole resonant sum is therefore a necessary condition for a mate.

Conceptually, the `b`th train is the defect expansion of

```text
lambda_b f^((n-b e)/m).                               (11a)
```

The high-weight part of the mate is a truncated Puiseux polynomial in `f`.
Each resonance lowers the permitted seed exponent by `e`; coprimality will
make the original `b=0` train the unique deepest `h`-adic pole.

### The graded rational centralizer

We first record the rigidity that makes the ladder unique.

> **Centralizer lemma.** If `R in C(x,y)^*` is `w`-homogeneous of degree `e`
> and `{h,R}=0`, then `d|e` and
> ```text
> R=lambda h^(e/d)                                    (12)
> ```
> for some `lambda in C^*`.

To prove it, the bracket equation makes the gradients of `R` and `h`
parallel in the rational function field. Weighted Euler identities give

```text
w_1 x h_x+w_2 y h_y=d h,
w_1 x R_x+w_2 y R_y=e R.                              (13)
```

Eliminating the proportionality factor between the gradients yields, with
`dR,dh` denoting differentials,

```text
d*(dR/R)=e*(dh/h),
R^d=C h^e                                             (14)
```

for a nonzero constant `C`. If

```text
h=product_l p_l^(a_l),            gcd_l a_l=1,        (15)
```

the valuations in (14) say that `d` divides every `e a_l`, hence `d|e`.
Taking valuations once more gives (12), after absorbing a constant `d`th
root. This also covers negative `e` and negative powers of `h`.

### Face recursion

The bracket of homogeneous faces of defects `a,b` has weighted degree

```text
D-a-b.                                                (16)
```

For `1<=j<=K`, this degree is positive by (8), whereas the right side of (3)
has degree zero. The defect-`j rho` bracket component must therefore vanish.
The sparseness assumption (6) leaves exactly

```text
m h^(m-1){h,G_j}+{A,G_(j-1)}=0.                       (17)
```

For any rational exponent `alpha` and starting index `j_0`, the train

```text
binom(alpha,j-j_0)h^[m alpha-m(j-j_0)]A^[j-j_0]       (18)
```

solves (17) after its seed at `j=j_0`. This follows directly from

```text
m r binom(alpha,r)
 =(m alpha-m(r-1))binom(alpha,r-1).                   (19)
```

We now prove (11) by induction. Continue every previously introduced train
one step using (17)--(19), and subtract their sum from `G_j`. The difference
centralizes `h`. Its weighted degree is

```text
n d-j rho.                                            (20)
```

If `L` does not divide `j`, then `d` does not divide `j rho`, so the degree in
(20) is not divisible by `d`. The centralizer lemma makes the difference
zero. If `j=bL`, then

```text
n d-j rho=d(n-b e),                                   (21)
```

and the only possible difference is `lambda_b h^(n-b e)`. It is exactly the
new `b`th seed in (11). This proves the full resonant formula. No genericity
of the resonance scalars is assumed. QED.

### The original pole survives every resonance

Write

```text
n=q m+r,                 q>=0, 0<r<m,                 (22)
```

take `K=q+1`, retain `K rho<D`, and assume

```text
gcd(h,A)=1.                                           (23)
```

At `j=K`, the exponent of `h` in the `b`th summand of (11) is

```text
r-m+b(mL-e)=r-m+b a/s.                               (24)
```

The `b=0` coefficient `c binom(n/m,q+1)` is nonzero. Its `h`-adic valuation is
strictly smaller than every `b>0` valuation because `a/s>0`; coprimality (23)
means the powers of `A` add no `h`-adic order. Thus no resonant train can
cancel the pole

```text
c binom(n/m,q+1) A^(q+1)/h^(m-r).                    (25)
```

It follows that a polynomial mate forces

```text
(q+1)rho>=D.                                          (26)
```

This is the key all-order point: resonance changes the formula, but not its
lowest `h`-adic tooth.

### Exact two-face coprime closure

Assume now that the component has exactly two nonconstant weighted faces,

```text
f=h^m+A,                 gcd(h,A)=1,                  (27)
```

with the hypotheses above but no mate top face fixed in advance. If `f` has
a polynomial Jacobian mate, then `f` is a polynomial coordinate and every
Keller pair containing it is an automorphism.

To prove this, repeatedly subtract a scalar power of `f` from a mate whenever
its leading face is `c h^n` with `m|n`. This preserves the bracket and strictly
lowers weighted degree. The descent cannot end with

```text
m d+deg_w(g)<=W.                                      (28)
```

Strict inequality leaves no bracket component of weight zero. At equality,
the top equation would be `{h^m,in_w(g)}=1`, impossible because its left side
is divisible by the nonconstant polynomial `h^(m-1)`.

The reduced mate therefore has leading face `c h^n` with

```text
n=q m+r,                 0<r<m,                       (29)
```

and `D>0`. The pole result forces `(q+1)rho>=D`. Since `rho=md-a`,

```text
D-(q+1)rho=(q+1)a+r d-W<=0.                           (30)
```

In particular `a+d<=W`. If `a+d<W`, then the homogeneous bracket `{h,A}` has
negative weighted degree and is zero. The centralizer lemma would make the
nonconstant `A` a positive power of `h`, contradicting (27). Hence

```text
a+d=W,                    {h,A}=kappa in C^*.          (31)
```

Substituting `W=a+d` back into (30) gives

```text
q a+(r-1)d<=0.
```

Positivity forces `q=0` and `r=1`: every reduced mate has top exponent
`n=1`. Thus equality has already collapsed the face descent to its terminal
tame control.

THM-2113 makes `(h,A/kappa)` a polynomial coordinate pair. Moreover

```text
q_0=-h/kappa,                   {f,q_0}=1,             (32)
```

and `(f,q_0)` is explicitly invertible: recover `h=-kappa q_0`, then
`A=f-h^m`, and finally invert the coordinates `(h,A/kappa)`. Any other
polynomial mate is `q_0+Q(f)`. This proves the exact two-face closure. QED.

## 2. Complete affine-root residue classification

Let `(u,v)` be polynomial coordinates on the affine plane, normalized by

```text
{u,v}=1.                                              (33)
```

Choose

```text
A,B in C[u],                    A!=0,
P in C[z],                      deg P>=2,              (34)
```

and put

```text
h=A(u)v+B(u),
f=u+P(h).                                              (35)
```

Then the following are equivalent:

1. there is a rational function `g in C(u,v)` with `{f,g}=1`;
2. there is a polynomial `g in C[u,v]` with `{f,g}=1`;
3. `A` is a nonzero constant.

When `A=a in C^*`, an explicit mate is

```text
g_0=h/a,                                              (36)
```

and `(f,g_0)` is a polynomial automorphism. Every polynomial mate is

```text
g=g_0+Q(f),                         Q in C[f],         (37)
```

and every rational mate has the same form with `Q in C(f)`.

### The constant case

Since

```text
{f,h}={u,h}=A(u),                                     (38)
```

equation (36) is a mate when `A=a`. The inverse to `(u,v)->(f,h/a)` is
polynomial:

```text
h=a g_0,
u=f-P(a g_0),
v=[a g_0-B(u)]/a.                                     (39)
```

In these coordinates, `{f,g}=1` is simply partial differentiation of `g`
with respect to `g_0`. This proves (37) and its rational analogue.

### A nonconstant coefficient creates a moving nonzero residue

Assume now that `A` is nonconstant and suppose a rational mate exists. Make
the birational change

```text
t=f,                         z=h.                     (40)
```

It has rational inverse

```text
u=t-P(z),
v=[z-B(t-P(z))]/A(t-P(z)),                            (41)
```

so `C(u,v)=C(t,z)`. Moreover

```text
{t,z}=A(t-P(z)).                                      (42)
```

Writing the putative mate as `G(t,z)`, its Jacobian equation becomes

```text
partial G/partial z=1/A(t-P(z)).                      (43)
```

We show that the rational differential on the right has a nonzero residue,
which an exact rational derivative cannot have.

Choose a root `alpha` of `A`, of multiplicity `k`. Around `alpha`, write

```text
1/A(alpha+s)=sum_(j=1)^k c_j s^(-j)+holomorphic,
c_k!=0.                                               (44)
```

Put `w_0=t-alpha`. Over an algebraic closure of `C(t)`, choose a root `beta`
of

```text
P(beta)=w_0.                                          (45)
```

It is a simple root: a simultaneous equation `P'(beta)=0` would make `beta`
and `P(beta)` constant, contrary to the transcendence of `t`. Let

```text
z=phi(w),                    H(w)=phi'(w)              (46)
```

be a local inverse branch of `P` through this root. Since
`s=t-alpha-P(z)=w_0-w`, the residue of the differential in (43) at `z=beta`
is

```text
R(w_0)=sum_(j=1)^k
 (-1)^j c_j H^(j-1)(w_0)/(j-1)!.                     (47)
```

This algebraic function is not identically zero. Indeed, if `M=deg P>=2`, an
inverse branch at infinity has a Puiseux expansion

```text
phi(w)=lambda w^(1/M)(1+O(w^(-1/M))),
H^(j-1)(w)=C_j w^(1/M-j)(1+O(w^(-1/M))),              (48)
```

where every `C_j` is nonzero. Let `j_0` be the least index with
`c_(j_0)!=0`. Its term in (47) has the unique largest exponent
`1/M-j_0`, so it cannot cancel with the terms of larger index. Thus
`R(w_0)!=0` in the algebraic function field.

The derivative of a rational function has zero residue at every place: in a
local Laurent series, the coefficient of the inverse first power would have
to come from differentiating the constant term. Equation (47) therefore
contradicts (43). No rational mate exists when `A` is nonconstant, completing
the equivalence and the theorem. QED.

## 3. Scope and sharpness

The affine-root classification allows arbitrary `B`, arbitrary lower
coefficients of `P`, and arbitrary degree at least two. It contains the tame
control

```text
f=x+(x+y)^2
```

after a coordinate choice, and it excludes families such as

```text
f=u+(u v)^m,                    m>=2.                 (49)
```

The condition `deg P>=2` is essential for the **rational** statement. For
example, when `P(z)=z`, `A(u)=u^k`, `B=0`, and `k>=2`, equation (43) has the
rational primitive

```text
[t-z]^(1-k)/(k-1).                                    (50)
```

The theorem does not settle a general proper-power top face. The full ladder
requires `f` to have no additional face before the tested defect; its exact
closure assumes only the two coprime faces `h^m,A`. The residue theorem instead
requires the primitive root to be affine in one polynomial coordinate. What
the two arguments add is an exact division of labor: resonant trains cannot
cancel the original `h`-adic pole in the two-face branch, while moving residues
kill every rational continuation in the affine-root family at once.
