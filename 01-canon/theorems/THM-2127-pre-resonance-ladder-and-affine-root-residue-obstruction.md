---
id: THM-2127
title: "Full resonant trains, radical-first-face closure, and affine-root residues"
status: >
  PROVED. THM-2102's proper-power first-defect equation has a full resonant
  solution: each resonance seeds a new fractional-binomial train, but every
  nonspecial train has higher factor valuation than the original one. With
  arbitrary later weighted faces, if the first lower face A is not divisible
  by rad(h), a mate forces h,A to be coordinates, every later face to lie in
  C[h], and the pair to be an automorphism. Thus every hard branch satisfies
  rad(h)|A. Separately, in polynomial coordinates (u,v), every
  f=u+P(A(u)v+B(u)) with deg P>=2 has a rational, hence polynomial, mate
  exactly when A is constant. These close broad proper-power subclasses but
  do not prove planar JC.
source: codex-2026-07-22-JC2-proper-power-all-order
depends_on:
  - THM-2102
  - THM-2113
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

The theorem has three complementary parts. The first solves the complete
resonant recursion for a sparse initial defect chain. The second assembles
**all** later weighted faces and shows that one irreducible factor of the top
root still sees a unique deepest pole unless `rad(h)|A`. The third changes to
an exact approximate-root coordinate and detects a nonzero moving residue.

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

Notice that `a=md-rho=s(mL-e)`, so `a/s` is a positive integer.

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

> **Centralizer lemma.** If `R in C(x,y)^*` is `w`-homogeneous of degree `E`
> and `{h,R}=0`, then `d|E` and
> ```text
> R=lambda h^(E/d)                                    (12)
> ```
> for some `lambda in C^*`.

To prove it, the bracket equation makes the gradients of `R` and `h`
parallel in the rational function field. Weighted Euler identities give

```text
w_1 x h_x+w_2 y h_y=d h,
w_1 x R_x+w_2 y R_y=E R.                              (13)
```

Eliminating the proportionality factor between the gradients yields, with
`dR,dh` denoting differentials,

```text
d*(dR/R)=E*(dh/h),
R^d=C h^E                                             (14)
```

for a nonzero constant `C`. If

```text
h=product_l p_l^(a_l),            gcd_l a_l=1,        (15)
```

the valuations in (14) say that `d` divides every `E a_l`, hence `d|E`.
Taking valuations once more gives (12), after absorbing a constant `d`th
root. This also covers negative `E` and negative powers of `h`.

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

The `b=0` coefficient `c binom(n/m,q+1)` is nonzero. Fix any irreducible
factor `pi` of `h`. Coprimality (23) gives `ord_pi(A)=0`, while (24) shows
that the `pi`-adic valuations of the summands strictly increase with `b`
because `a/s>0`. Thus no resonant train can cancel the pole

```text
c binom(n/m,q+1) A^(q+1)/h^(m-r).                    (25)
```

It follows that a polynomial mate forces

```text
(q+1)rho>=D.                                          (26)
```

This is the key all-order point: resonance changes the formula, but not its
lowest `h`-adic tooth.

## 2. Arbitrary later faces and the radical-first-face closure

Now let

```text
f=h^m+A+R,                                             (F1)
```

where `h` and `m` are as above, `A` is the first lower homogeneous face,

```text
deg_w(A)=a=md-rho>0,                                  (F2)
```

and every face of `R` has degree strictly below `a`. Assume the sharp factor
condition

```text
rad(h) does not divide A.                             (F3)
```

If `f` has a polynomial Jacobian mate, then

```text
{h,A}=kappa in C^*,
f=h^m+A+Q(h)                         for Q in C[T],    (F4)
```

and `(f,-h/kappa)` is a polynomial automorphism pair. Thus every hard
proper-power branch must satisfy `rad(h)|A`.

### Mate reduction

Repeatedly subtract `c f^(n/m)` from a mate whenever its leading face is
`c h^n` with `m|n`. This preserves the bracket and strictly lowers integral
weighted degree. Whenever `md+deg_w(g)>W`, the positive top bracket with
`h^m` vanishes, so the polynomial centralizer lemma makes the new leading
face another `c h^n`. The finite descent cannot finish below weight sum,
where no bracket term has weight zero, or at equality, where it would demand

```text
{h^m,in_w(g)}=1,
```

impossible because the left side is divisible by `h^(m-1)`. We therefore
reach a reduced mate with

```text
in_w(g)=c h^n,          n=q m+r,       0<r<m,
D=(m+n)d-W>0,           K=q+1.                         (F5)
```

### The full-face formal train

Index every weighted face by defect:

```text
f=sum_(delta>=0) f_delta,       deg_w(f_delta)=md-delta,
g=sum_(t>=0) g_t,               deg_w(g_t)=nd-t,       (F6)
```

so `f_0=h^m`, `f_rho=A`, no `f_delta` occurs for `0<delta<rho`, and
`g_0=c h^n`. Introduce a bookkeeping variable `z` and put

```text
F(z)=sum_delta f_delta z^delta,
G(z)=sum_t g_t z^t.                                  (F7)
```

The weighted decomposition of `{f,g}=1` is exactly

```text
{F(z),G(z)}=z^D.                                      (F8)
```

For every integer `ell>=0`, define in `C(x,y)[[z]]`

```text
S_ell(z)=h^(n-ell)[F(z)/h^m]^((n-ell)/m).             (F9)
```

The binomial series is well-defined because the bracketed series has constant
term one, and `{F,S_ell}=0`. For every `T<D`, there are scalars `lambda_ell`,
with `lambda_0=c`, such that

```text
G(z)=sum_(ell d<=T) lambda_ell z^(ell d) S_ell(z)
                         modulo z^(T+1).              (F10)
```

This is the full resonant train with **no** sparseness assumption on the later
faces. To prove (F10), subtract all trains seeded below a defect `t` and let
`H_t` be the first residual coefficient. Since `t<D`, equations (F8)--(F9)
give `{h^m,H_t}=0`. The residual has degree `nd-t`; the graded centralizer
lemma makes it zero unless `t=ell d`, when it is
`lambda_ell h^(n-ell)`. Subtracting the corresponding train kills it and
continues the finite induction.

### The decisive tooth survives every later face

Suppose `K rho<D`. In the coefficient of `z^(K rho)`, the `ell=0` train has
the nonzero term

```text
c binom(n/m,K) A^K/h^(m-r).                           (F11)
```

Every term in every train is indexed by a seed `ell`, a number `N` of selected
lower faces, and defects `delta_1,...,delta_N`, with

```text
ell d+sum_i delta_i=K rho,             delta_i>=rho.  (F12)
```

The special term (F11) is exactly `ell=0`, `N=K`, and every
`delta_i=rho`. Every other term has `N<K`, and

```text
ell d<=(K-N)rho<(K-N)md,
ell<m(K-N).                                           (F13)
```

Its bare `h` exponent `n-ell-mN` is therefore strictly greater than
`n-mK=r-m`.

By (F3), choose an irreducible `pi|h` with `pi` not dividing `A`. The special
term has negative valuation `(r-m)ord_pi(h)`. Every other term has strictly
larger `pi`-valuation because all later faces are polynomials. Thus (F11) is
the unique deepest pole and cannot cancel, contradicting polynomiality of
`g_(K rho)`. Therefore

```text
K rho>=D.                                             (F14)
```

### Terminal closure

Substituting `K=q+1`, `rho=md-a`, and `n=qm+r` into (F14) gives

```text
(q+1)a+r d<=W,
a+d<=W.                                               (F15)
```

If `a+d<W`, then `{h,A}=0`; if equality held but the bracket still vanished,
the centralizer lemma would make `A` a positive rational power of `h`. Either
conclusion makes every irreducible factor of `h` divide `A`, contrary to
(F3). Hence `a+d=W` and `{h,A}=kappa!=0`.

THM-2113 makes `h` a coordinate. If `v_0` is a coordinate mate, then
`{h,A/kappa-v_0}=0`, so in coordinates `(h,v_0)` one has

```text
A/kappa=v_0+P(h).                                     (F16)
```

Thus `(h,A/kappa)` is a polynomial coordinate pair. Express every later
homogeneous face in those graded coordinates of degrees `d,a`. Its degree is
strictly below `a`, so it contains no positive power of `A/kappa` and belongs
to `C[h]`. Summing the faces gives (F4), and

```text
q_0=-h/kappa,                  {f,q_0}=1              (F17)
```

is explicitly an automorphism mate. QED.

### Sharp radical boundary

The hypothesis (F3) is sharp for this mechanism. With ordinary weights,

```text
f=y^4+y^2+x,                  g=y,                    (F18)
```

one has `{f,g}=1`, `h=y`, `m=4`, `A=y^2`, `rho=2`, and the later face `x`
has defect three. Here `rad(h)|A`. The reduced mate has `n=1`, `K=1`, and
`D=3`, so `K rho<D`; nevertheless

```text
F(z)^(1/4)=y+(1/4)y^(-1)z^2+O(z^3),                  (F19)
```

and the resonant `ell=2` seed contributes `-(1/4)y^(-1)z^2`, cancelling the
original pole exactly. The remaining hard object is therefore not “resonance”
alone but cancellation of factor-initial forms on the divisor `rad(h)|A`.

## 3. Complete affine-root residue classification

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
g=g_0+Q(f),                         Q in C[T],         (37)
```

and every rational mate has the same form with `Q in C(T)`.

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
and `P(beta)` constant, contrary to the transcendence of `t`. The polynomial
`P(z)-w` is irreducible over `C(w)` because the rational map `P^1_z->P^1_w`
has function-field degree `deg P`; hence the finite local inverse germ below
and an inverse Puiseux place at infinity lie on the same algebraic function
field. Let

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

Here `H^(j-1)` denotes the `(j-1)`st derivative with respect to `w`, not a
power.

This algebraic function is not identically zero. Indeed, choose a place above
`w=infinity` on the same algebraic inverse branch. If `M=deg P>=2`, it has a
Puiseux expansion

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

## 4. Scope and sharpness

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

The theorem does not settle a general proper-power top face. The arbitrary-
later-face closure leaves precisely the factor-initial locus `rad(h)|A`, where
the control (F18)--(F19) shows that a later resonant seed really can cancel the
first pole. The residue theorem instead requires the primitive root to be
affine in one polynomial coordinate. What the arguments add is an exact
division of labor: away from the factor-initial locus, no collection of later
faces or resonant seeds can cancel the unique deepest factor valuation; on the
affine-root locus, moving residues kill every rational continuation at once.
