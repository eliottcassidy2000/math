---
id: THM-3783
title: "Quadratic tower etale-surface maximal polynomial observable"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the m=1
  boundary of the THM-3774 rational Keller tower, the complete polynomial
  intersection with the rational target field is the smooth surface
  r^3g=z^2-r/4.  The source affine plane maps onto it etale with generic
  degree two and one sheet escaping along its unique multiple arm.  Its
  symplectic form is exact, its units are constant, and Pic is Z/2.  Thus
  the exponent-one obstruction of THM-3779 genuinely disappears, but every
  Darboux pair has surface degree at least three by THM-3794 and would give
  a planar Keller map of even degree at least six.  Neither homogeneous
  nor two-by-two Euler-weight supports can do this; one named aligned
  two-by-three orientation also has a universal transport-factor no-go.
  Other two-by-three cells and arbitrary Darboux pairs remain open.
source: jc_sparse_direct_search / THM-3779 quadratic-boundary continuation, 2026-08-23
audit: >
  TWO INDEPENDENT HOSTILE AUDITS PASSED (root and jc_zero_debt_lift,
  2026-08-23).  The audits rederived
  the generator and field identities, smoothness, surjective nonfinite
  etaleness, height-one DVR descent for the full intersection, the unique
  doubled arm and Picard class, units, LND kernel, global Euler primitive,
  degree floor, homogeneous anatomy, complete two-by-two sign/factor cases,
  and the aligned two-by-three transport factor.  The Broughton dependencies
  and degree-two Galois gate were checked against their proved scopes.
  Normal and optimized runs byte-match the hardened frozen output; all
  three raw hashes and 94 active gates agree.
depends_on:
  - THM-3774-three-component-rational-keller-cover-tower
  - THM-3716-monomial-broughton-hamiltonian-obstruction-family
  - THM-1330-keller-monoid-exact-picture-inverse-jelonek-cusp-rule
  - THM-3794-constant-unit-surfaces-have-no-quadratic-etale-plane-map
related:
  - THM-3779-three-component-tower-maximal-danielewski-polynomial-observable
  - THM-3782-simple-pole-spectral-danielewski-completion-and-target-field-gate
script: 04-computation/jc2_quadratic_tower_etale_surface_maximal_observable_thm3783.py
output: 05-knowledge/results/jc2_quadratic_tower_etale_surface_maximal_observable_thm3783.out
script_sha256: 772e78bcaba0bf459d9842716fc86ffa4a27991941da40b7fc4a3b96dd761068
output_sha256: 515fa9c4e55dc117865410da4e92f74b0d50b4bba259a93f4c8622d9cadb7265
semantic_sha256: 98d73f9751663c7d1ed2bcb73335ee71d956be9173d59d8864225c2e8defd170
hash_basis: raw LF bytes
---

# THM-3783 -- the quadratic boundary is a smooth etale surface

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This
theorem computes the sharp `m=1` exception left open by THM-3779.  It does
not construct a planar Jacobian counterexample.  It does isolate a new exact
surface on which such a counterexample could live.

Work over an algebraically closed field `k` of characteristic zero.  In the
`m=1` member of THM-3774 put

```text
A=1+xy,
B=1+x^3A,
U=xAB,
P=(2B-1)/x.                                           (1)
```

Then

```text
J(U,P)=1,
t^2-Pt+U=0,                  [k(x,y):k(U,P)]=2,       (2)
```

where `t=x^2A`.  The rational pair is not polynomial because of the pole of
`P` on `x=0`.

## 1. Three forced polynomial observables

The discriminant identity in the exceptional quadratic row is

```text
P^2-4U=1/x^2.                                         (3)
```

It forces the first extra polynomial target-field observable

```text
r=x^2=1/(P^2-4U).                                     (4)
```

Two more convenient observables are

```text
z=x/2+x^5y=rP/2-r^2,
g=y+x^4y^2.                                           (5)
```

They satisfy

```text
r^3g=z^2-r/4.                                         (6)
```

Conversely,

```text
U=rg+2z+r^2,
P=2(z+r^2)/r.                                         (7)
```

Thus

```text
k(U,P)=k(r,z)=Frac k[r,z,g]/(r^3g-z^2+r/4).           (8)
```

The algebraic independence of `r,z` follows already from

```text
J(r,z)=2r^3!=0.                                       (9)
```

Hence `(6)` is the complete relation among the three observables.

For comparison, THM-3779's exponent-one observables remain polynomial but
are strict sub-observables.  If `V=UP` and `E=P(V-1)`, then in the new ring

```text
V=1+2gz+6r^2g+6rz+2r^3,                              (10)

E=g+3r+4r^2g^2+16rgz+24r^3g+16r^2z+4r^4,            (11)

UE=V(V-1).                                            (12)
```

The missing divisor in the old completion is exactly what makes `r,z,g`
visible.

## 2. The smooth surface and its exact source atlas

Put

```text
S=k[r,z,g]/(r^3g-z^2+r/4),              Y=Spec S.     (13)
```

The gradient of the defining equation is

```text
(3r^2g+1/4, -2z, r^3),                                 (14)
```

which has no common zero.  Thus `Y` is smooth and normal.  Equations
`(4),(5)` define a morphism

```text
phi:A2_(x,y) -> Y.                                    (15)
```

Its three Jacobian minors are

```text
J(r,z)=2r^3,
J(r,g)=4z,
J(z,g)=1/2+6r^2g.                                    (16)
```

They are twice the three signed partial derivatives in `(14)`, so they
generate the unit ideal.  Therefore `phi` is etale.

It is also surjective.  Given a geometric point `(rho,zeta,eta)` with
`rho!=0`, choose `xi` with `xi^2=rho` and set

```text
x=xi,
y=(zeta-xi/2)/xi^5.                                  (17)
```

Then

```text
y+x^4y^2=(zeta^2-rho/4)/rho^3=eta                    (18)
```

by `(6)`.  On `rho=0`, equation `(6)` forces `zeta=0`, and

```text
x=0,                         y=eta                    (19)
```

gives the required point.  Hence

```text
phi(A2)=Y.                                             (20)
```

The generic fibre has two points because `x^2=r` and `r` is not a square
in `k(r,z)`.  But the arm

```text
L=V(r,z) ~= A1_g                                      (21)
```

has the single source lift `x=0,y=g`.  Thus `phi` is a surjective etale
map of generic degree two which is not finite: its second sheet escapes at
infinity exactly along `L`.  This is the geometric information lost by
looking only at the quadratic field degree.

## 3. The exact maximal intersection

Inside `k(x,y)`, let

```text
R_1=k[x,y] intersection k(U,P).                       (22)
```

Equations `(4)--(6)` give `S subset R_1`.  Conversely, let `F in R_1` and
let `q` be any height-one prime of `S`.  Surjectivity and etaleness of
`phi` supply a height-one prime `p` of `k[x,y]` above `q`.  The induced DVR
extension is unramified, so

```text
ord_p(F)=ord_q(F).                                    (23)
```

The left side is nonnegative because `F` is a source polynomial.  Hence
`ord_q(F)>=0` for every height-one prime of the normal ring `S`.  The Krull
intersection of its height-one DVRs proves

```text
k[x,y] intersection k(U,P)
   =k[r,z,g]/(r^3g-z^2+r/4).                          (24)
```

This is the complete repair of THM-3779 at `m=1`, not a bounded generator
guess.

## 4. Units, Picard torsion, and one additive action

Localizing `(13)` at `r` eliminates `g`:

```text
S_r=k[r,r^(-1),z].                                    (25)
```

The only height-one prime above `r=0` is `L=(r,z)`.  At its generic point,
`z` is a uniformizer and `(6)` gives

```text
ord_L(r)=2,                       div(r)=2L.           (26)
```

Nagata's class-group sequence applied to `(25)` says that `Cl(S)` is
generated by `[L]`, and `(26)` gives `2[L]=0`.  The class is not zero: if
`L=div(f)`, then `f` becomes a unit in `(25)`, hence `f=c r^n`; its divisor
would be `2nL`, never `L`.  Smoothness identifies class group and Picard
group.  The same unit calculation gives

```text
S^*=k^*,                         Pic(Y)=Z/2.           (27)
```

There is also an explicit locally nilpotent derivation

```text
delta(r)=0,                  delta(z)=r^3,
delta(g)=2z.                                             (28)
```

It preserves `(6)` and is visibly locally nilpotent.  On `(25)` it is
`r^3 partial_z`, so

```text
ker_S(delta)=k[r].                                    (29)
```

Consequently the Makar--Limanov invariant satisfies

```text
ML(S) subset k[r].                                    (30)
```

Equality in `(30)` is not claimed here.

## 5. Poisson structure and the vanished de Rham obstruction

The source bracket induces on `S` exactly the packet `(16)`.  Its inverse
symplectic form is characterized on `r!=0` by

```text
omega=dr wedge dz/(2r^3),                 phi^*omega=dx wedge dy.       (31)
```

Unlike the exponent-one surface of THM-3779, this form is exact.  Give
`(r,z,g)` Euler weights

```text
wt(r)=2,                    wt(z)=1,       wt(g)=-4.   (32)
```

The relation has weight two and the Poisson bracket has degree `+3`.  The
Euler vector field

```text
Eul=2r partial_r+z partial_z-4g partial_g             (33)
```

satisfies `Lie_Eul(omega)=-3omega`.  Cartan's formula therefore gives the
global regular primitive

```text
alpha=-(1/3) i_Eul omega,                  d alpha=omega.               (34)
```

On the source,

```text
phi^*alpha=-(4y dx+x dy)/3,               d(phi^*alpha)=dx wedge dy.   (35)
```

Thus the logarithmic cohomology gate genuinely disappears at the quadratic
boundary.  Exactness is a construction permission, not a construction.

## 6. What any Darboux pair would have to accomplish

Suppose now that `k=C` and `F,G in S` satisfy `{F,G}=1`.  Their pullbacks
would be a polynomial Keller pair on the source.  Put

```text
d=[k(U,P):k(F,G)].                                    (36)
```

The bracket condition makes this a finite positive integer, and `(2),(8)`
give

```text
[k(x,y):k(F,G)]=2d.                                   (37)
```

If `d=1`, the source extension would be quadratic and hence Galois.
THM-1330's Galois Keller theorem would force invertibility and degree one,
a contradiction.  THM-3794 excludes `d=2` directly on this constant-unit
affine surface.  Therefore every hypothetical Darboux pair in `S` obeys

```text
d>=3,                 [k(x,y):k(F,G)]>=6 and even.    (38)
```

It must have surface degree at least three; both the birational and quadratic
target-change routes are completely closed.

The three most visible component candidates already fail separately:

1. `U` has the three-component residue obstruction of THM-3774;
2. `2z=x+2x^5y` becomes `X+X^5Y` after rescaling `Y=2y`;
3. `g=y+y^2x^4` is the transposed `m=2,n=4` monomial Broughton profile.

THM-3716 therefore excludes polynomial mates for both `z` and `g`.  A live
pair must mix the surface generators rather than select one of them.

## 7. Exact Euler-weight anatomy

Put

```text
s=1+2x^4y.                                            (39)
```

In the quadratic extension, the deck involution is

```text
(x,s) |-> (-x,-s),                                    (40)
```

and

```text
r=x^2,             z=xs/2,             g=(s^2-1)/(4x^4).               (41)
```

Every nonzero homogeneous element of `S` of weight `a in Z` has a unique
form

```text
x^a f(s),                         f in k[s],           (42)
```

where

```text
f(-s)=(-1)^a f(s),
(s^2-1)^ceil(-a/4) divides f(s)       when a<0.        (43)
```

Conversely `(43)` makes `(42)` a polynomial in `r,z,g`.  Necessity comes
from regularity at `x=0`; parity transports the required zero at `s=1` to
`s=-1`.  Sufficiency follows from `(41)` and the remaining exponent
`a+4ceil(-a/4) in {0,1,2,3}`.

For two homogeneous elements, direct differentiation gives

```text
J(x^a f(s),x^b h(s))
 =2x^(a+b+3)[a f h'-b f'h].                           (44)
```

### 7.1 No homogeneous Darboux pair

A scalar bracket requires `a+b=-3`.  If both weights are negative, `(43)`
makes both profiles vanish at `s=1`, so the bracket in `(44)` vanishes
there.  If `a>=0>b`, a nonzero value at `s=1` forces the negative profile
to have a simple zero and `a!=0`.  Since `b=-3-a`, this leaves only

```text
(a,b)=(1,-4)                       or (-4,1).           (45)
```

In the first case `f` is odd of degree at least one and `h` is divisible by
`s^2-1`.  The scalar equation is

```text
2(fh'+4f'h)=1.                                        (46)
```

Its leading coefficient is proportional to
`deg(h)+4deg(f)>0`, at degree `deg(f)+deg(h)-1>=2`.
This is impossible.  The other orientation is identical with the roles
reversed.  Hence there is no homogeneous Darboux pair.

### 7.2 Neither output may be homogeneous

If one output has one weight and the other has several, the different
weight components in `(44)` land in distinct output weights.  The unique
bucket of weight zero would itself be a homogeneous scalar bracket,
already excluded.  Thus every live output has at least two Euler weights.

### 7.3 The complete two-by-two no-go

Let the exact active weight supports be

```text
supp(F)={a,a+d},              supp(G)={b,b+e},         d,e>0.           (47)
```

If `d!=e`, all four pair sums are distinct, so a scalar again requires a
forbidden homogeneous bracket.  Hence `d=e`.  The scalar cannot occupy an
endpoint bucket for the same reason; it must occupy the double middle
bucket.  Therefore

```text
a+b+d=-3,                                             (48)

{F_a,G_b}=0,                 {F_(a+d),G_(b+d)}=0.     (49)
```

Commuting nonzero homogeneous profiles of opposite nonzero weights do not
exist: `(44)` would force positive divisor multiplicities to have opposite
sign.  A weight-zero member that commutes with a nonzero weight is constant.
These observations leave the following exhaustive regimes.

* For `d<=3`, all nonconstant profiles in the middle bracket vanish at
  `s=1`; any weight-zero coefficient in `(44)` is itself zero.  The scalar
  bucket vanishes at `s=1`.
* For `d>=4`, the zero/positive endpoint boundary reduces the middle scalar
  equation to `A fH'+B f'H=constant`, whose leading degree is positive
  because the negative profile contains `s^2-1`.
* In the remaining regime both low weights are negative and both top
  weights positive.  Write `A=-a,B=-b`, so

  ```text
  A+B=d+3,          C=a+d=B-3,          D=b+d=A-3.    (50)
  ```

  Endpoint commutation and unique factorization give

  ```text
  f=p^(A/delta),    h=lambda p^(B/delta),   delta=gcd(A,B),
  F=q^(C/epsilon),  H=mu q^(D/epsilon),     epsilon=gcd(C,D).           (51)
  ```

  The profile `p` is nonconstant by `(43)`.  If `A!=B`, the two middle
  summands have different degrees, with difference

  ```text
  (B-A)[deg(p)/delta+deg(q)/epsilon],                  (52)
  ```

  and the leading coefficient of the higher one is nonzero.  If `A=B`,
  the middle bracket is a nonzero scalar multiple of
  `(A-3)p'q+A p q'`, or zero; its degree is positive because
  `deg(p)>=2`.  Neither case is a nonzero constant.

Thus no Darboux pair has at most two Euler weights in each output.  The
first support shapes not closed by this theorem are `2 by 3`, `3 by 2`,
and `3 by 3`, rather than another sparse `2 by 2` decoration.

## 8. One aligned two-by-three orientation also peels

The first recurrent aligned placement has, for `d>=4`,

```text
supp(F)={1-d,1},
supp(G)={-4,d-4,2d-4},                                (53)
```

with the scalar in the bucket

```text
(1-d)+(d-4)=1+(-4)=-3.                                (54)
```

Write

```text
A=d-1,       C=A-3,       E=2A-2,       N=2A-3,
delta=gcd(A,4),       u=A/delta,         v=4/delta.   (55)
```

The bottom and top singleton buckets must commute.  After absorbing
nonzero scalars, their four profiles are

```text
F_(1-d)=p^u,              G_(-4)=lambda p^v,
F_1=q,                    G_(2d-4)=mu q^E.            (56)
```

Cancellation in the other double bucket is a first-order exact derivative
and gives the complete middle profile

```text
G_(d-4)=mu E p^u q^N+nu q^C.                          (57)
```

Substitution into the scalar bucket factors it as

```text
(q p'+delta p q') times
[
 -mu E N u p^(2u-1)q^(2A-4)
 -nu C u p^(u-1)q^(A-4)
 +lambda v p^(v-1)
].                                                     (58)
```

For `A=3`, the term carrying `C=0` is simply absent.  The first factor in
`(58)` is nonconstant: both `p` and the odd profile `q` are nonconstant, and
its leading coefficient is

```text
deg(p)+delta deg(q)>0.                                (59)
```

A product of polynomials cannot be a nonzero scalar when one factor is
nonconstant.  Thus `(53)` is empty.  The small seams `d=1,2,3` die directly
from `(43),(44)`; the companion separately checks the first integrated seam
`d=3`.

Literal output-swap gives the dual `3 by 2` orientation.  No claim is made
for every aligned placement or for irregular `2 by 3` cells; these are the
next exact coefficient systems to attack.

## 9. Controls and scope

The companion named in the metadata checks all identities through `(41)`,
the smooth/etale Groebner unit ideals, the rational deck involution, the
Euler primitive and locally nilpotent derivation, 32 reduced generator
monomials, and `(57),(58)` for `A=3,...,12`.  Normal and optimized executions
must byte-match the frozen transcript before promotion.

The theorem proves:

1. the complete maximal polynomial target-field observable for `m=1`;
2. its surjective nonfinite etale degree-two geometry;
3. its units, Picard torsion, Poisson packet, and exact symplectic primitive;
4. the even-degree-at-least-four floor for any Darboux pair;
5. the homogeneous and complete `2 by 2` weight-support no-go; and
6. the named aligned `2 by 3` transport-factor no-go.

It proves neither `ML(S)=k[r]`, nor a complete `2 by 3` classification, nor
the existence or nonexistence of an arbitrary Darboux pair on `Y`, nor a
planar Jacobian counterexample.  **QED.**
