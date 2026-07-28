---
id: THM-2621
title: "Planar degree-four inverse-spectral Keller congruence and sheet-defect pole ledger"
status: >
  PROVED + VERIFIED-SYMBOLIC + INDEPENDENTLY HOSTILE-AUDITED AFTER
  MISTAKE-301 SCOPE REPAIR.  A hypothetical
  planar Keller map of generic field degree four admits generic determinant-one
  source coordinates in which P is monic in y and x is primitive over
  K=C(P,Q).  If f is the monic inverse minimal quartic of x and y=b(x), then
  the Keller condition is exactly the coefficientwise congruence
  f_v b_u-f_u b_v = kappa^(-1) f_T modulo f.  The four reduced coefficient
  equations are explicit.  The resultant is R=c f, and on every prime
  Jelonek component D|c the coefficient indexed by the generic surviving
  fibre cardinality has full pole order nu_D(c).  Formal quartic depression
  is not symplectic without a dh wedge dy correction.  For every polynomial
  Keller realization the trace-Liouville form Tr(x dy)-4 kappa^(-1)u dv is
  exact: it is d Tr(H), where x dy-kappa^(-1)P dQ=dH on the affine source.
  Its affine poles are second-kind poles supported on the Jelonek boundary
  and every divisor residue is zero.  Branchwise Laurent cancellation, a
  localized power-sum potential, and in the D4 lane a relative opposite-pair
  potential survive before the full trace; these are explicit realization
  gates absent from the PDE alone.  At every D4 deck-pole prime, some affine
  determinant-one source gauge makes the anti-invariant potential polar
  there, although its differential remains exact with zero residue.  A
  rational D4 sheet-defect hostile permits
  arbitrary residue once polynomial realization is removed.  An exact
  punctured C4 hostile also satisfies the entire
  PDE with c=1 but has primitive-coordinate discriminant and b-denominator
  divisors, proving that the PDE alone does not encode affine polynomial
  realization.  This reduces but does not solve the D4/A4/S4 planar degree-four
  frontier, JC(2), or DC(2).
source: codex-2026-07-28-planar-degree-four-inverse-spectral
depends_on:
  - THM-2241-monic-transverse-response-depth-and-resultant-nonproper-quotient
related:
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
  - THM-2546-integral-coordinate-dichotomy-and-parity-lens-scope
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2612-d4-deck-pole-tax-and-depressed-resolvent-gcd-gate
  - THM-2627-d4-jelonek-quadratic-character-rank-and-component-gate
  - THM-2628-d4-opposite-pair-escape-and-deck-pole-census
  - MISTAKE-301
script: 04-computation/jacobian_planar_degree4_inverse_spectral_thm2621.py
output: 05-knowledge/results/jacobian_planar_degree4_inverse_spectral_thm2621.out
script_sha256: baf7034775833a26f814cbf82e71df679bdec8919e3834d0a38bbbd08d4234b2
output_sha256: d3e31d2ec9ca251d1e7bd4fd07e147f6335d388b4ce45fe7b72e042af73ea23c
hash_basis: working-tree bytes (LF)
---

# THM-2621 -- a degree-four Keller map is a marked inverse quartic with boundary poles

**PROVED + VERIFIED-SYMBOLIC + INDEPENDENTLY HOSTILE-AUDITED AFTER
MISTAKE-301 SCOPE REPAIR.**

The planar field-degree-four branch is open.  THM-2465 reduces every affine
point-cap instance to it, while THM-2598 leaves geometric monodromy
`D_4,A_4,S_4`.  Abstract quartic-resolvent identities do not remember the
second inverse coordinate or the symplectic connection.

This theorem supplies exactly that missing marked sidecar.  It replaces an
unrestricted polynomial coefficient hunt by a rational pair `(f,b)`: the
minimal polynomial of one primitive inverse coordinate and the degree-three
interpolant recovering the other.  The Keller equation becomes four
first-order coefficient equations, and every missing-sheet divisor forces an
exact pole in a specified coefficient of `f`.

## 1. Generic monic primitive coordinates

Let

```text
F=(P,Q): A^2_C -> A^2_C,

u=P,                 v=Q,

Jac_(x,y)(u,v)=kappa in C*,

[C(x,y):C(u,v)]=4.                                      (1)
```

There is a generic affine `SL_2(C)` source coordinate system `(x,y)` in
which

```text
P is monic in y,                  L=C(x,y)=K(x),

K=C(u,v).                                                (2)
```

Indeed, a linear form fails to be primitive only when two of the four generic
conjugates agree on it, a finite union of proper directions.  Monicity fails
only when the chosen `y` direction lies in the finite zero set of the top
homogeneous part of `P`.  These are nonempty Zariski-open conditions on the
source frame.  Scaling the two frame directions inversely preserves
determinant one and normalizes the leading coefficient to one, so `kappa` is
unchanged.

Let

```text
R(T;u,v)=Res_Y(P(T,Y)-u,Q(T,Y)-v).                       (3)
```

THM-2241 says that its generic degree is the generic fibre cardinality and
that its leading coefficient cuts the nonproper-value set.  Since `x` is
primitive and the generic degree is four,

```text
R(T;u,v)=c(u,v) f(T;u,v),

c=[T^4]R in C[u,v] minus {0},

f=T^4+a_3T^3+a_2T^2+a_1T+a_0 in K[T]                    (4)
```

with `f` the irreducible monic minimal polynomial of `x`.  The four generic
roots are distinct because the extension is separable.  Since `L=K(x)`, there
is a unique companion

```text
y=b(x),

b(T)=b_0+b_1T+b_2T^2+b_3T^3 in K[T].                    (5)
```

The pair `(f,b)` is the marked inverse spectral object.  The cubic resolvent
of `f` forgets `b`; that forgotten coordinate is precisely where the Keller
connection lives.

## 2. The Keller equation is one quartic congruence

All partial derivatives of `f` and `b` below are **coefficientwise at fixed
`T`**.  Differentiating `f(x;u,v)=0` gives

```text
x_u=-f_u/f_T,                     x_v=-f_v/f_T.          (6)
```

Meanwhile

```text
y_u=b_u+b_T x_u,                  y_v=b_v+b_T x_v.       (7)
```

The `b_T` terms cancel in the inverse Jacobian:

```text
det partial(x,y)/partial(u,v)
 =x_u y_v-x_v y_u
 =(f_v b_u-f_u b_v)/f_T.                                 (8)
```

Because `f_T(x)!=0` generically, (8) proves the equivalence

```text
Jac_(x,y)(u,v)=kappa

iff

f_v b_u-f_u b_v == kappa^(-1) f_T        modulo f.       (9)
```

This equivalence holds inside the existing function-field presentation.  An
abstract rational pair `(f,b)` satisfying (9) defines a generically
symplectic cover; it need not algebraize to polynomial coordinates on
`A^2`.

### 2.1 The four coefficient equations

Put

```text
D_k=sum_(i+j=k) ((a_i)_v (b_j)_u-(a_i)_u (b_j)_v),

0<=i,j<=3.                                               (10)
```

Reduction of `sum_(k=0)^6 D_k T^k` modulo `f` gives

```text
rho_0=D_0-a_0D_4+a_3a_0D_5+a_0(a_2-a_3^2)D_6,

rho_1=D_1-a_1D_4+(a_3a_1-a_0)D_5
      +(a_0a_3+a_1a_2-a_1a_3^2)D_6,

rho_2=D_2-a_2D_4+(a_3a_2-a_1)D_5
      +(-a_0+a_1a_3+a_2^2-a_2a_3^2)D_6,

rho_3=D_3-a_3D_4+(a_3^2-a_2)D_5
      +(-a_1+2a_2a_3-a_3^3)D_6.                         (11)
```

Thus (9) is exactly the four first-order rational PDEs

```text
(rho_0,rho_1,rho_2,rho_3)
 =kappa^(-1)(a_1,2a_2,3a_3,4).                          (12)
```

## 3. Formal depression loses the marked symplectic connection

The algebraic substitution

```text
t=x+h(u,v)                                               (13)
```

changes the inverse symplectic form by

```text
dt wedge dy
 =dx wedge dy+dh wedge dy.                               (14)
```

For the quartic (4), the depressed choice is `h=a_3/4`, because the cubic
coefficient of `f(T-h)` is `a_3-4h`.  Unless `dh wedge dy=0`, this is not a
Keller-preserving source-coordinate change.  It is target-dependent and
generally rational.

Thus depression is lawful for abstract Galois and resolvent calculations but
cannot be inserted into (9) while retaining the original right side.  The
missing term (14) explains why THM-2598's universal depressed-resolvent
identities alone do not constrain a polynomial Keller realization.

## 4. Every Jelonek component forces a full inverse-coefficient pole

Write the polynomial resultant as

```text
R=cT^4+r_3T^3+r_2T^2+r_1T+r_0,             r_j=c a_j.   (15)
```

Let `D` be an irreducible prime divisor of `c`, put

```text
e=nu_D(c)>0,                                               (16)
```

and reduce `R` at the generic point of `D`.  The reduction is nonzero:
otherwise the two specialized equations would share a curve, contradicting
the quasi-finiteness of an etale Keller map.  Define

```text
k_D=deg_T(R mod D)<4.                                     (17)
```

By THM-2241 this is the generic geometric fibre cardinality along `D`.
Coincident `x` values are counted with the sum of their simple point
multiplicities, so (17) remains valid when the primitive projection
specializes noninjectively.

By definition of the specialized degree,

```text
nu_D(r_(k_D))=0,

nu_D(r_j)>0                         for j>k_D.            (18)
```

Dividing by `c` gives the exact pole law

```text
nu_D(a_(k_D))=-e.                                        (19)
```

Equivalently, if `delta_D=4-k_D` sheets escape generically along `D`, then
the first surviving coefficient `a_(4-delta_D)` has the full pole order of
the resultant lead.  In particular a one-sheet defect forces

```text
nu_D(a_3)=-e,                    a_3=-Tr_(L/K)(x).        (20)
```

At least one normalized inverse coefficient therefore has a pole on every
Jelonek component.  If all `a_j` were polynomial, every coefficient of
`R=c f` would vanish modulo `D`, contradicting the nonzero specialization.

The same mechanism recovers THM-2546's sporadic cubic core

```text
R=L T^3+pT-2c.
```

At generic `L=0,p!=0`, one has `k_D=1`; the two-sheet defect appears as the
full pole `a_1=p/L`.  Its absent quadratic term comes from an already regular
depressed escape coordinate, not a target-dependent formal shift.

### 4.1 Three boundary ledgers must remain separate

The theorem distinguishes:

```text
c=0:              loss of finite sheets / Jelonek nonproperness;

disc_T(f)=0:      collision of primitive x-coordinates;

denominators(b):  failure of the chosen power-basis reconstruction.       (21)
```

These divisors may meet but are not interchangeable.  Equation (19) concerns
only primes of `c`; a pole of `b` need not signal a missing physical sheet.

## 5. A punctured `C_4` hostile satisfies the whole PDE

On `G_m x A^1` with coordinates `(Y,X)`, define

```text
u=Y^4,                     v=(X-Y)/(4Y^3).              (22)
```

Then

```text
Jac_(X,Y)(u,v)=-1.                                       (23)
```

This is a finite etale cyclic degree-four cover of `G_m x A^1`.  The
coordinate `X` is primitive and has minimal polynomial

```text
f(T)=T^4-16uvT^2-256u^3v^4+32u^2v^2-u,                 (24)
```

while

```text
Y=b(X),

b(T)=((1+48uv^2)T-4vT^3)/(1-256u^2v^4).                (25)
```

Exact elimination and reduction give

```text
Res_Y(Y^4-u,T-Y-4vY^3)=f(T),

f_v b_u-f_u b_v == -f_T                     modulo f,   (26)
```

as required by `kappa=-1`.

Here the resultant leading coefficient is

```text
c=1,                                                       (27)
```

but

```text
den(b)=1-256u^2v^4,

disc_T(f)=-256u^3(16uv^2-1)^2(16uv^2+1)^4.              (28)
```

The poles in (25) occur when the chosen primitive `X` projection ceases to
separate sheets; they are not Jelonek poles.  This exact hostile proves that
the PDE (9), even with polynomial coefficients of `f`, does not encode affine
polynomial realization.  In (22), `v` is Laurent rather than polynomial, the
source omits `Y=0`, and the base is not simply connected.  It is not a Keller
counterexample on `A^2` and does not contradict Campbell's exclusion of
polynomial `C_4` monodromy.

## 6. The trace--Liouville class vanishes for polynomial realizations

Let the trace run over the four inverse sheets and define the rational
one-form

```text
alpha=Tr_(L/K)(x dy).                                    (29)
```

On every branch, `dx wedge dy=kappa^(-1)du wedge dv`.
Exterior differentiation and trace therefore give

```text
d alpha=4 kappa^(-1) du wedge dv.                        (30)
```

Consequently

```text
omega_F=alpha-4 kappa^(-1)u dv                           (31)
```

is closed.  For an actual **polynomial** Keller map there is a stronger fact
that the first version of this theorem missed.  On the affine source put

```text
theta_F=x dy-kappa^(-1)P dQ.                            (31a)
```

The Keller equation gives

```text
d theta_F=dx wedge dy-kappa^(-1)dP wedge dQ=0.          (31b)
```

This is a closed polynomial one-form on `A^2`, hence it has a polynomial
primitive `H_F`, unique up to a constant.  Explicitly, if
`theta_F=A dx+B dy`, the radial homotopy formula

```text
H_F(x,y)=integral_0^1 [x A(tx,ty)+y B(tx,ty)] dt
```

is polynomial and has derivative `theta_F`.  Trace commutes with exterior
differentiation in the finite separable extension `L/K`, while
`Tr(P dQ)=4u dv`.  Therefore

```text
boxed: omega_F=d Phi_F,          Phi_F=Tr_(L/K)(H_F).   (31c)
```

The trace--Liouville class is identically zero and every divisor residue of
`omega_F` vanishes.  In particular, no hypothetical polynomial `D_4,A_4`, or
`S_4` Keller map can have a nonzero trace--Liouville residue.  Conversely, a
nonzero residue computed from an abstract rational pair `(f,b)` obstructs its
extension to a polynomial Keller map.

There is still an exact support and pole-order ledger.  Off the Jelonek set
the Keller map is finite etale.  Since `H_F` is regular upstairs, the
finite-algebra trace makes `Phi_F` regular there.  Thus every affine pole of
`Phi_F` or `omega_F` is supported on the Jelonek set.  If at its generic
point

```text
Phi_F=t^(-m)(a_0+a_1 t+...),        m>0,       a_0!=0,
```

then

```text
omega_F=d Phi_F=-m a_0 t^(-m-1)dt+t^(-m)d_D a_0+... .  (31d)
```

The normal pole order is exactly `m+1>=2` and the logarithmic residue is
zero: these are second-kind poles.

Before trace there is a sharper branchwise condition.  At every divisorial
valuation over a finite target divisor, `u,v` are regular and
`x dy=dH_F+kappa^(-1)u dv`; hence

```text
Res(x dy)=0.                                             (31e)
```

For a uniformizer `t` and Laurent expansions

```text
x=sum_i x_i t^i,                   y=sum_j y_j t^j,
```

this is the explicit constraint

```text
boxed: sum_j j x_(-j)y_j=0.                             (31f)
```

It is imposed separately on every normalization branch and retains
information that corestriction can cancel.

There is also a coefficient-only version that can be imposed directly on an
abstract pair `(f,b)`.  Put `p_m=Tr_(L/K)(x^m)`.  Newton's identities give

```text
p_1=-a_3,
p_2=a_3^2-2a_2,
p_3=-a_3^3+3a_3a_2-3a_1,
p_4=a_3^4-4a_3^2a_2+2a_2^2+4a_3a_1-4a_0.               (31g)
```

Since `y=sum_(j=0)^3 b_j x^j`, put

```text
E=sum_(j=0)^3 b_j p_(j+1)=Tr_(L/K)(xy),

C=sum_(j=0)^3 b_j/(j+1) dp_(j+1).                       (31h)
```

Direct differentiation on the four sheets gives

```text
alpha=sum_(j=0)^3 [p_(j+1) db_j+j/(j+1)b_j dp_(j+1)],

alpha=dE-C,

d alpha=sum_(j=0)^3 1/(j+1) dp_(j+1) wedge db_j
       =4 kappa^(-1)du wedge dv.                         (31i)
```

Combining this with (31c) gives

```text
Gamma_F:=C+4 kappa^(-1)u dv
        =d(E-Phi_F).                                     (31j)
```

This exactness has a support refinement.  Put

```text
V=A^2 minus A_F=D(c),                 U=F^(-1)(V).
```

Both `H_F` and `xy` are regular on `U`, and finite-algebra trace preserves
regular functions.  Therefore

```text
Phi_F,E,E-Phi_F in O(V)=C[u,v,1/c].                     (31k)
```

Thus the **entire affine polar principal part** of `alpha`, `omega_F`, and
`Gamma_F` cancels away from `c=0`.  Pure primitive-discriminant and pure
power-basis reconstruction divisors cannot support the aggregate trace form,
even when individual `a_i,b_j` have poles there.  At a Jelonek prime the
remaining principal part is exact and has zero logarithmic residue.

Because `u dv` is regular at every affine base prime, polynomial realization
also forces the explicit residue shadow

```text
boxed: Res_D sum_(j=0)^3 b_j/(j+1) dp_(j+1)=0
                           for every affine base prime D. (31l)
```

If `b_j=sum_n b_(j,n)pi^n` and
`p_(j+1)=sum_n p_(j+1,n)pi^n` in a uniformizer `pi`, this is

```text
sum_(j=0)^3 1/(j+1) sum_(r in Z)
       r b_(j,-r)p_(j+1,r)=0.                            (31m)
```

At a Jelonek component, (31l)--(31m) must coexist with the full inverse-coefficient
pole in (19).  It is a marked-companion constraint, not a consequence of the
pole order or cubic resolvent alone.

Exactness of `omega_F`, and its primitive modulo constants and exact affine
quadratic corrections, is invariant under determinant-one affine source and
target coordinate changes.  For the linear source part

```text
x'=ax+by,                  y'=cx+dy,          ad-bc=1,
```

one has the exact identity

```text
x' dy'-x dy
 =d(ac x^2/2+bc xy+bd y^2/2).                           (32)
```

The target Liouville term obeys the same identity, and translations add only
exact linear differentials.  Thus (31) is the exact trace--Liouville
differential of the multivalued inverse, independent of the chosen affine
symplectic frame.  The useful sidecar is its primitive modulo constants and
its branchwise splitting before trace, not a potentially nonzero de Rham
class.

For the hostile (22), direct trace gives

```text
alpha=4v du,

d alpha=-4du wedge dv=4 kappa^(-1)du wedge dv,

omega_F=4d(uv).                                          (33)
```

Here the source form has the explicit primitive

```text
H_F=uv+Y^2/2,
```

and `Phi_F=Tr(H_F)=4uv`, exactly as (31c) predicts.

In the `D_4` lane, let `tau` be THM-2612's nontrivial root-field involution.
Trace kills its entire anti-invariant channel, so retain instead

```text
G_tau=H_F-tau(H_F).                                     (33a)
```

Then, without any further hypothesis,

```text
dG_tau=x dy-tau(x)d tau(y),
tau(G_tau)=-G_tau,
Tr_(L/K)(G_tau)=0,
G_tau^2 in E_2=L^tau,       Norm_(L/E_2)(G_tau)=-G_tau^2. (33b)
```

For deck transformations `g`, the family `G_g=H_F-g(H_F)` obeys the action
cocycle identity `G_(gh)=G_g+g(G_h)`.  A shear `y -> y+s x` changes `G_tau`
by `(s/2)(x^2-tau(x)^2)` (and the transverse shear does the analogous job for
`y`).  This closes the gauge-existence step.  Let `mathfrak e` be a pole prime
of `tau`.  The affine functions `x,y` are regular at `mathfrak e`, while at least one of
`tau(x),tau(y)` is polar by THM-2612.  If `tau(x)` is polar, then
`x^2-tau(x)^2` is polar, and in the quotient of the local field by its
valuation ring at `mathfrak e` at most one scalar `s` can cancel the principal
part of `G_tau`.  Thus all but at most one shear parameter make the new
`G_tau` polar at `mathfrak e`.  If only `tau(y)` is polar, use
`x -> x+s y`.  If one also wants to retain the primitive/monic chart of
Section 1, exclude its finite set of bad shear parameters as well; the ground
field is infinite.  Consequently,

```text
for every deck-pole prime mathfrak e, some affine SL_2 source gauge
makes G_tau polar at mathfrak e.                         (33b')
```

This is still an exact-potential pole, not a residue obstruction:
`dG_tau` has zero residue at every divisor.

In the notation of THM-2627, `L=E=Omega^H`, and the intermediate called
`E_2=L^tau` here is exactly `M_deck=Omega^J=E^tau`, not the distinct
matching/discriminant field `M_Delta=Omega^V`.  Moreover
when `G_tau!=0`, its squareclass in `M_deck` describes the upper quadratic
extension `E/M_deck`; the pole-forcing gauge above guarantees this nonzero
condition.  It does not say that this squareclass
generates `M_deck/K` or has the deck-character parity vector.  The square may
even lie in `K`.  This is a typed pole probe, not a monodromy exclusion.

There is a complementary invariant-channel gate before the quadratic
`M/K` trace.  Put

```text
M=L^tau,                       [L:M]=[M:K]=2,

p_m^tau=x^m+tau(x)^m,

C_tau=sum_(j=0)^3 b_j/(j+1) d p_(j+1)^tau.              (33c)
```

Taking the relative trace of `theta_F=dH_F` gives

```text
Gamma_tau:=C_tau+2 kappa^(-1)u dv
 =d Tr_(L/M)(xy-H_F)                         in Omega^1_M. (33d)
```

On the finite-etale quotient `W=U/<tau>`, the pair potential on the right is
regular.  Thus every divisor of `M` has zero residue **before** its two
`M/K` conjugates can cancel, and its affine polar support is confined to the
boundary omitted from `W`.  This is strictly finer than the base gate
(31j)--(31m).  Together, (33b) and (33d) retain both the anti-invariant action
potential and the invariant opposite-pair trace.

THM-2612's deck-coordinate pole is compatible with these exactness laws.  A
coordinate pullback can have a pole while the associated trace--Liouville
differential has only exact higher-order principal parts and zero logarithmic
residue.  The pole invoice and the residue invoice must not be identified.

The polynomial hypothesis in (31a) is load-bearing.  On the punctured rational
family

```text
u=Y^4,                    v=(X-h(Y))/(4Y^3),
```

one still has `kappa=-1`, but now

```text
theta_F=X dY+u dv=d(uv)+h(Y)dY.                         (33e)
```

Taking `h(Y)=Y-4lambda/Y` gives

```text
omega_F=4d(uv)-4lambda du/u,
```

so arbitrary traced residue appears as soon as polynomial extension is
removed.  Taking

```text
h(Y)=1/(Y-1)-1/(Y-i)
```

instead gives opposite branch residues over `u=1`, while
`Norm((Y-1)/(Y-i))=(1-u)/(i^4-u)=1`, so their trace vanishes.  Thus zero
traced residue does not recover branchwise polynomial exactness.  These are
punctured rational controls, not polynomial maps of `A^2`.

### 6.1 A `D_4` sheet-defect hostile with arbitrary residue

The polynomial hypothesis remains load-bearing even after fixing `D_4`
monodromy and THM-2621's exact sheet-defect pole.  On a rational symplectic
chart put

```text
u=s^4+s^2,                         v=t/(4s^3+2s),

X=1/(s-1),
Y=-t(s-1)^2+lambda(s-1).                              (34a)
```

Then

```text
dX wedge dY=ds wedge dt=du wedge dv.                   (34b)
```

The primitive coordinate `X` has eliminant

```text
R=(2-u)X^4+6X^3+7X^2+4X+1=(2-u)f.                     (34c)
```

At `u=2`, `R` specializes to a cubic, so one sheet escapes and
`a_3=6/(2-u)` has the full pole prescribed by (19).  The extension
`C(s)/C(s^4+s^2)` has root involution `s -> -s`; its normal closure is
obtained by adjoining `sqrt(-1-s^2)` and has degree eight.  Hence its
geometric monodromy is `D_4`.  The squared-pair resolvent is

```text
U(U^2+2U+1+4u).                                         (34d)
```

The companion `Y=b(X)` has polynomial base coefficients

```text
b_0=2(-2lambda+10uv-3v),
b_1=-7lambda+26uv-26v,
b_2=2(-3lambda+11uv-16v),
b_3=-(u-2)(-lambda+4uv-6v).                             (34e)
```

Exact reduction gives

```text
f_v b_u-f_u b_v == f_X,                                 (34f)

disc_X(R)=-16u(4u+1)^2,
disc_X(f)=-16u(4u+1)^2/(u-2)^6.                         (34g)
```

Thus `u=2` is neither a raw discriminant collision nor a companion
denominator.  Nevertheless the coefficient trace is

```text
omega_F=d(-(20u+4)v)+lambda dlog(2-u),

Res_(2-u)(omega_F)=lambda.                              (34h)
```

The deck involution itself is

```text
tau(X)=-X/(2X+1),

tau(Y)=(2X+1)[2lambda-(2X+1)Y].                         (34i)
```

It fixes `u,v`, has `tau^2=1` and Jacobian one, and its pole lies over
`u=2`.  The parameter `lambda` is arbitrary.  Thus the quartic PDE,
one-sheet coefficient pole, polynomial companion coefficients, and `D_4`
monodromy do not imply even the base residue gate.  At `lambda=0`, the deck
pole remains while every logarithmic residue vanishes, proving that
second-kind data can carry the pole invoice and that zero residue is not a
`D_4` exclusion.

### 6.2 Base exactness does not imply opposite-pair exactness

A second `D_4` hostile proves that (33d) is strictly finer than the full base
trace.  Put

```text
u=s^4+s^2,                    v=t/(4s^3+2s),

x=s,                          y=t-6s/(u-2).              (35a)
```

Again `dx wedge dy=du wedge dv`.  Now

```text
f=T^4+T^2-u,

b(T)=4vT^3+[2v-6/(u-2)]T,                               (35b)

f_v b_u-f_u b_v == f_T modulo f.                        (35c)
```

The quartic is monic, so `c=1`, and

```text
disc_T(f)=-16u(4u+1)^2.                                 (35d)
```

Thus `u=2` is a pure companion-denominator divisor, not a Jelonek or
primitive-discriminant divisor.  The full base trace is nevertheless exact:

```text
omega_F
 =[12v-12/(u-2)^2]du+(12u+4)dv
 =d([12u+4]v+12/(u-2)).                                 (35e)
```

So its de Rham class and residues vanish by cancellation.  Its pole at
`u=2` violates the stronger support clause (31k), since `c=1` would require
the physical trace potential to be polynomial.

The intermediate quadratic field is `M=C(r,v)` with

```text
r=s^2,                         u=r^2+r.                  (35f)
```

After removing the exact term `d(2uv)`, the pair form (33d) is

```text
Gamma_tau-d(2uv)
 =-2 dlog g(r),                g(r)=(r-1)/(r+2).         (35g)
```

It has residues `-2,+2` at `r=1,-2`.  The nontrivial `M/K` conjugation is
`r -> -1-r` and sends `g` to `g^(-1)`, so these two classes cancel in the
base trace:

```text
g(r)g(-1-r)=1.                                             (35h)
```

Hence base exactness cannot replace opposite-pair exactness.  Polynomial
realization requires two refinements of the base class which must remain
separate:

```text
                         base de Rham/residue cancellation
                              /                    \
 affine support in C[u,v,1/c]     D4 pair exactness before M/K trace.
                                                               (35i)
```

This hostile fails both upper refinements, so it proves that each contains
information absent from the base trace; it does not order the two refinements
against each other.  The hostile is rational, not a polynomial Keller map.

## 7. Resolvent interface and exact residual

THM-2598's matching cubic remembers the three pairings of the four roots and
quotients by `V_4`.  It forgets

```text
the attached values y_i=b(x_i),
the coefficientwise symplectic congruence (9),
the trace--Liouville primitive, branch splitting, and localized coefficient gate (31c)--(31m),
the anti-invariant action potential and invariant opposite-pair gate (33a)--(33d),
and which normalization branches remain in the affine source.             (34)
```

Thus `(b,(9),(31c)--(31m),(33a)--(33d))` is the missing marked-origin sidecar for transferring any
grade-three resolvent anatomy back to a quartic Keller map.  In the `D_4`
lane, it must coexist with THM-2612's birational deck-involution pole and the
present/omitted boundary ownership of THM-2598.  These are compatible
invoices, not yet a contradiction.

The degree-four planar frontier is now the following typed classification
problem: classify rational pairs `(f,b)` satisfying (9), the pole law (19),
the polynomial exactness, support, branch, and pair laws (31c)--(33d), monodromy
`D_4,A_4`, or `S_4`, and extension of their total space to `A^2` with
`u,v in C[x,y]`.  This is strictly smaller than an unrestricted planar
coefficient search, but no lane is excluded here.

## 8. Exact companion and scope

Run

```text
python3 04-computation/jacobian_planar_degree4_inverse_spectral_thm2621.py
python3 -O 04-computation/jacobian_planar_degree4_inverse_spectral_thm2621.py
```

Both modes byte-match the stored transcript.  The exact symbolic companion
checks:

1. all four universal remainder formulas (11);
2. the Jacobian, elimination, companion recovery, and PDE congruence for the
   punctured `C_4` hostile;
3. the exact discriminant and reconstruction-denominator factorization;
4. twelve sheet-defect valuation controls covering `k=0,1,2,3` and pole
   orders `1,2,3`;
5. the hostile trace-curvature sign and explicit source primitive;
6. the universal source-defect curvature cancellation under the Keller
   equation;
7. the arbitrary traced-residue and norm-one branch-cancellation punctured
   hostiles;
8. the `D_4` sheet-defect hostile's eliminant, companion, PDE, discriminants,
   deck involution, and arbitrary residue;
9. the strict opposite-pair hostile's exact base trace, pure companion
   divisor, conjugate pair residues, and norm-one cancellation; and
10. the determinant-one Liouville exact-differential identity.

The generic-coordinate intersection, resultant factorization, valuation law,
and trace argument are proofs above, not extrapolations from the fixtures.
An independent hostile audit rederived the PDE sign and coefficient formulas,
checked the generic-coordinate and quasi-finite specialization arguments,
verified the `C_4` model and trace signs independently, and enforced the three
separate boundary ledgers in (21).  A subsequent hostile audit caught the
closed-versus-exact gap and supplied (31a)--(31f); MISTAKE-301 records the
failed implication and repair.

No polynomial degree-four Keller map, monodromy lane exclusion, JC(2), DC(2),
or GMC-to-JC interface follows.  In particular, a nonzero divisor residue of
`omega_F` is impossible rather than merely unproved; useful residue data must
be retained branchwise or pairwise before trace as in (31e)--(31f) and
(33d).  No local mechanism forcing the support and exactness gates is proved
here.

QED.
