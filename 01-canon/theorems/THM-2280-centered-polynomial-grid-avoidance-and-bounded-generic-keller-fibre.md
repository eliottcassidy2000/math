---
id: THM-2280
title: "Centered polynomial-grid avoidance and bounded generic Keller fibres"
status: >
  PROVED. Over a characteristic-zero field, finitely many nonzero
  multivariate polynomials can be avoided simultaneously on a centered
  integer grid whose coordinate radius is half the sum of their separate
  coordinate degrees, rounded up. The grid cardinality and centered radius
  are sharp for arbitrary polynomials with only those degree bounds. For a
  monic planar Keller map with y-degrees d=deg_y(P) and e=deg_y(Q), this
  produces an integer target (u,v), with |u|<=ceil(e/2) and
  |v|<=ceil(d/2), whose fibre has the exact generic cardinality N.
  Finite banks of nonzero target polynomials can be avoided at the summed
  bidegree bound. A scoped THM-2240 corollary likewise avoids a finite
  polynomial bank along a finite-dimensional grade-six gauge slice only
  when every pulled-back predicate is nonzero. This is bounded algebraic
  genericity, not JC(2), DC(2), or a continuation theorem.
source: codex-2026-07-25-centered-polynomial-avoidance
depends_on:
  - THM-2240-dc2-grade-response-gauge-is-not-a-continuation-state
  - THM-2241-monic-transverse-response-depth-and-resultant-nonproper-quotient
related:
  - THM-2260-positive-finite-fibre-capacity-and-thin-predicate-boundaries
  - THM-2270-simultaneous-balanced-cut-relation-and-six-uniform-orientation
---

# THM-2280 -- centered polynomial avoidance and a bounded Keller fibre

THM-2270 multiplies finitely many nonzero linear restrictions and evaluates
the product on a centered coefficient grid. The underlying interpolation
argument does not require linear factors. This theorem records the
separate-degree version and applies it where the target predicate is
actually polynomial.

The two applications have different logical force:

```text
monic Keller resultant:
  an unconditional bounded integer target with the generic fibre count;

DC grade-six gauge:
  a representative-level avoidance statement, conditional on every
  predicate remaining nonzero after pullback to the chosen gauge slice.
                                                                  (1)
```

The second clause is deliberately weaker. THM-2240 proves that the current
DC response is not a continuation state, so algebraic genericity of one
representative cannot be promoted to a quotient or termination statement.

## 1. Centered separate-degree avoidance

For an integer `n>=0`, define the centered consecutive set

```text
S(n)={-floor(n/2),-floor(n/2)+1,...,ceil(n/2)}.       (2)
```

It has exactly `n+1` elements and

```text
max{|z|:z in S(n)}=ceil(n/2).                         (3)
```

Let `K` be a field of characteristic zero. Take integers `r,s>=1` and
nonzero polynomials

```text
f_1,...,f_s in K[T_1,...,T_r].                       (4)
```

For every coordinate put

```text
a_(j,i)=deg_(T_i) f_j,
A_i=sum_(j=1)^s a_(j,i).                             (5)
```

Then there is an integer point

```text
t=(t_1,...,t_r) in Z^r                               (6)
```

such that

```text
f_j(t)!=0                         for every j,
|t_i|<=ceil(A_i/2)                for every i.       (7)
```

Here the integers are mapped into `K` through its prime field.

### Proof

Form

```text
F(T)=product_(j=1)^s f_j(T).                          (8)
```

The polynomial ring over a field is an integral domain, so `F` is nonzero.
Viewing the factors as one-variable polynomials in `T_i` over the integral
domain

```text
K[T_1,...,T_(i-1),T_(i+1),...,T_r]
```

shows that leading coefficients cannot cancel under multiplication. Hence

```text
deg_(T_i) F=A_i                                      (9)
```

for every `i`.

We use the elementary grid-vanishing lemma:

> If `G in K[T_1,...,T_r]` has
>
> ```text
> deg_(T_i)G<=n_i
> ```
>
> for every `i`, and `E_i subset K` has `n_i+1` distinct elements, then
> `G` cannot vanish on `E_1 x ... x E_r` unless `G=0`.

For completeness, prove it by induction on `r`. The one-variable statement
is the root bound. For the induction step write

```text
G=sum_(k=0)^(n_r) G_k(T_1,...,T_(r-1)) T_r^k.       (10)
```

Fix a point of `E_1 x ... x E_(r-1)`. If `G` vanishes on the full grid,
the resulting polynomial in `T_r` has `n_r+1` roots and degree at most
`n_r`, so it is zero. Thus every coefficient `G_k` vanishes on the
smaller grid. The induction hypothesis makes every `G_k` zero, proving
the lemma.

Apply the lemma to `F` and

```text
E_i=S(A_i).                                          (11)
```

Characteristic zero makes their integer images distinct. Since `F` is
nonzero, some point `t` of the product grid has `F(t)!=0`. Every factor in
(8) is then nonzero at `t`, and (3) gives (7). QED.

## 2. Exact grid boundary

Both numerical features of Section 1 are sharp for a theorem based only on
the separate coordinate degrees.

First, a degree-`n` polynomial can vanish at any prescribed `n` scalar
points:

```text
G(T_i)=product_(z in E)(T_i-z),          |E|=n.      (12)
```

Thus `n+1` points in coordinate `i` cannot uniformly be replaced by `n`.

Second, if `n+1` distinct integers all have absolute value at most `R`,
then

```text
n+1<=2R+1,
R>=ceil(n/2).                                        (13)
```

The set `S(n)` attains this radius. Therefore the bound in (3) is the
smallest possible centered integer radius for `n+1` sample points.

These are universal grid-certificate boundaries. They do not claim that a
particular polynomial lacks a smaller nonzero integer value.

## 3. Simultaneous polynomial banks

For later use, specialize Section 1 to two variables. Let

```text
g_1,...,g_s in K[U,V]\{0},

A=sum_j deg_U g_j,
B=sum_j deg_V g_j.                                  (14)
```

There is a point `(u,v) in Z^2` such that

```text
|u|<=ceil(A/2),
|v|<=ceil(B/2),
g_j(u,v)!=0                     for every j.         (15)
```

No irreducibility, coprimality, transverse-intersection, or real-positivity
hypothesis is hidden here. The load-bearing facts are exactly:

```text
each g_j is nonzero in the ambient polynomial ring;
the coefficient field has characteristic zero;
the sampling carrier contains the full centered product grid.      (16)
```

If a purported predicate is the zero polynomial, or becomes zero after
restriction to a smaller parameter space, the product argument supplies
nothing.

## 4. A bounded integer target for a monic Keller map

Let

```text
P,Q in C[x,y],
Jac(P,Q)=kappa in C*,

d=deg_y P>=1,
e=deg_y Q>=0,                                        (17)
```

and assume that `P` is monic in `y`. Define the THM-2241 resultant

```text
R_F(X;U,V)
 =Res_Y(P(X,Y)-U,Q(X,Y)-V).                          (18)
```

Let `N` be the generic degree of `F=(P,Q)` and put

```text
c(U,V)=[X^N]R_F(X;U,V).                              (19)
```

THM-2241 proves that `c` is nonzero and that, at every target point,

```text
c(u,v)!=0
 iff deg_X R_F(X;u,v)=N
 iff #F^(-1)(u,v)=N.                                (20)
```

The last cardinality is geometric cardinality. The Keller condition makes
every local intersection multiplicity one.

### 4.1 The resultant bidegree

The Sylvester determinant is homogeneous of degree `e` in the coefficient
block of

```text
P(X,Y)-U
```

and homogeneous of degree `d` in the coefficient block of

```text
Q(X,Y)-V.                                            (21)
```

This includes the boundary `e=0`, where no coefficient from the first block
is selected and the resultant is the `d`th power of the degree-zero second
polynomial.

The variable `U` occurs linearly in one coefficient of the first block,
and `V` occurs linearly in one coefficient of the second. Therefore

```text
deg_U R_F<=e,
deg_V R_F<=d.                                       (22)
```

Taking one `X`-coefficient cannot increase either target degree, so

```text
deg_U c<=e,
deg_V c<=d.                                         (23)
```

These bounds may be strict through cancellation; only the upper bounds are
used.

### 4.2 The bounded generic fibre

Apply Section 1 directly to the nonzero polynomial `c`, using the grids

```text
S(e) x S(d).                                         (24)
```

Equations (22)--(23) and the grid-vanishing lemma give an integer target

```text
(u,v) in Z^2,

|u|<=ceil(e/2),
|v|<=ceil(d/2),
c(u,v)!=0.                                           (25)
```

By (20),

```text
#F^(-1)(u,v)=N.                                      (26)
```

This is unconditional under (17). If `c` is constant, its actual
bidegree is `(0,0)` and `(u,v)=(0,0)` already works. THM-2241 then says the
map is proper and hence an automorphism. In the nonconstant case, (25)
merely finds one bounded lattice point off the nonproper curve `V(c)`.

### 4.3 Additional target predicates

Let

```text
h_1,...,h_m in C[U,V]\{0}.                           (27)
```

Apply Section 3 to the bank

```text
c,h_1,...,h_m.
```

There is an integer target with

```text
|u|<=ceil((e+sum_l deg_U h_l)/2),
|v|<=ceil((d+sum_l deg_V h_l)/2),                    (28)

c(u,v)h_1(u,v)...h_m(u,v)!=0.                       (29)
```

Thus the target has the exact generic fibre count `N` and avoids every
listed algebraic hypersurface simultaneously.

The hypotheses `h_l!=0` are ambient statements. If one is obtained by
restricting a polynomial from a larger parameter space, its restriction
must first be checked to remain nonzero.

## 5. Scoped DC representative-fibre corollary

We now apply the same lemma to a finite-dimensional slice of THM-2240's
grade-six correction fibre. This is not a Keller resultant statement and
does not descend to the grade-response quotient.

THM-2240's standard correction has

```text
Atilde_0=-(3/4)q(u-1),
Btilde_0=0,                                           (30)
```

and normalized next residual

```text
H_7^0=-(3/2)q(10u^2-36u+29).                        (31)
```

Fix `n>=0` and parameters `t=(t_0,...,t_n)`. On the integration-constant
axis take

```text
a_t(u)=sum_(i=0)^n t_i u^i,

Atilde_t=Atilde_0+a_t,
Btilde_t=0.                                          (32)
```

Because `a_t` is independent of `q`, equation (10) of THM-2240 gives

```text
d_6(Atilde_t,Btilde_t)
 =d_6(Atilde_0,Btilde_0).                            (33)
```

Every point of this parameter space therefore represents the same
grade-six cancellation. Its exact grade-seven residual is

```text
H_7(t)=H_7^0+L_a(a_t),                               (34)

L_a(a)=20(u-2)a+2(2u^2-6u+1)a'.                     (35)
```

In particular, every coefficient of `H_7(t)` is an affine-linear function
of `t`.

Let

```text
Phi_1,...,Phi_s                                       (36)
```

be polynomial predicates in the finitely many coefficients of the
polynomial `H_7(t)`. Define their pullbacks

```text
p_j(t)=Phi_j(H_7(t)) in Q[t_0,...,t_n].              (37)
```

Assume explicitly that

```text
p_j is not the zero polynomial             for every j.   (38)
```

Put

```text
A_i=sum_(j=1)^s deg_(t_i)p_j.                        (39)
```

Section 1 supplies an integer gauge representative satisfying

```text
|t_i|<=ceil(A_i/2),
Phi_j(H_7(t))!=0                    for every j.      (40)
```

If `Phi_j` has total degree at most `delta_j`, then affine-linearity in
(34) gives

```text
deg_(t_i)p_j<=delta_j.
```

Thus, with

```text
Delta=sum_j delta_j,
```

the coarser uniform bound

```text
|t_i|<=ceil(Delta/2)                                 (41)
```

is always valid.

The restriction (38) is load-bearing. THM-2240 proves that `L_a` is
injective, and already

```text
L_a(1)=20(u-2)!=0,                                   (42)
```

so this slice genuinely changes the next residual while preserving the
current response. But injectivity of `L_a` does not imply that an arbitrary
polynomial predicate is nonzero on its affine image. A predicate may vanish
on the entire `a`-axis. One must then enlarge the gauge slice, for example
to a finite-dimensional part of the `J` or `b` axes, and recheck the
pullback. THM-2240's explicit `C=1`, `J=q` splitter,

```text
Delta H_7=8q(4u^2-13u+13),                           (43)
```

is the cheapest control showing why different axes can retain different
information.

Consequently (40) says only:

```text
a finite bank of proper algebraic conditions cannot cover this
chosen representative slice when every pullback is nonzero.       (44)
```

It does not say that the predicates are continuation-invariant, that a
later residual vanishes, or that any correction ladder terminates.

## 6. Connection and loss ledger

### Monic Keller application

```text
source:
  the nonzero top-X resultant coefficient c(U,V) of THM-2241;

target:
  one bounded integer target carrying all N generic sheets;

map:
  use Sylvester multihomogeneity to bound the U,V degrees, then
  evaluate c and any auxiliary target polynomials on a centered grid;

preserved:
  simultaneous algebraic nonvanishing and the exact fibre count N;

destroyed:
  the embedded geometry of V(c), target-shear covariance as a
  coordinate-free bound, and every route from N to degree one;

needed sidecar:
  the resultant coefficient c and an independent argument forcing
  c to be constant or N=1;

cheapest hostile test:
  a nonconstant c may miss a small integer point while its zero curve
  still contains every nonproper target.                              (45)
```

### DC representative application

```text
source:
  a finite-dimensional affine slice of THM-2240's grade-six
  cancellation fibre;

target:
  one representative avoiding a finite next-residual predicate bank;

map:
  pull every predicate back to the gauge parameters and apply the
  centered grid only after verifying that no pullback is zero;

preserved:
  the current grade-six response and the listed representative-level
  nonvanishing predicates;

destroyed:
  quotient invariance, later residuals, state-dependent continuations,
  and termination data;

needed sidecar:
  a gauge-invariant continuation profile or explicit future residual
  data, plus a nonzero-pullback proof for every predicate;

cheapest hostile tests:
  L_a(1)=20(u-2) on the a-axis and the C=1 splitter (43).              (46)
```

## 7. Scope

The theorem is an exact finite algebraic selection result. It proves:

1. a sharp centered product-grid lemma for arbitrary finite polynomial
   banks over characteristic zero;
2. an unconditional bounded integer target with the generic fibre count
   for every monic planar Keller map; and
3. a conditional, representative-sensitive DC gauge-bank corollary.

It does not prove that the generic Keller degree is one, make the lattice
bound invariant under arbitrary complex coordinate changes, detect the
whole nonproper curve from one fibre, turn the DC response quotient into a
continuation state, or prove `JC(2)` or `DC(2)`. No companion computation is
needed: every bound follows from exact degree counting and the displayed
interpolation proof. QED.
