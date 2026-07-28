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
  is not symplectic without a dh wedge dy correction.  After MISTAKE-301's
  scope repair, the trace-Liouville form
  Tr(x dy)-4 kappa^(-1)u dv is not merely closed: every genuine polynomial
  Keller realization makes it rationally exact, so all of its divisor
  residues vanish.  Its companion-coordinate formula gives an explicit
  residue-cancellation gate absent from the quartic PDE alone.  A rational
  `D4` hostile satisfies the PDE and one-sheet pole law with arbitrary
  residue `lambda`, proving the gate independent even after monodromy is
  fixed.  An exact punctured C4 hostile also satisfies the entire
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
  - MISTAKE-301
script: 04-computation/jacobian_planar_degree4_inverse_spectral_thm2621.py
output: 05-knowledge/results/jacobian_planar_degree4_inverse_spectral_thm2621.out
script_sha256: ff340be47192826e18140b9a2130ddd84bcbebb7c29b4f857696fd05727af74a
output_sha256: c596499c2c859ed1a8f47b55f104f4265ae6f9952db01802ae11119cb325c8af
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

is closed.  For an abstract rational pair `(f,b)` this is all that follows
from the coefficientwise PDE.  A genuine polynomial Keller realization has
strictly more structure.

On the polynomial source put

```text
theta=x dy-kappa^(-1)u dv.                              (31a)
```

It is a polynomial one-form and

```text
d theta=dx wedge dy-kappa^(-1)du wedge dv=0.             (31b)
```

The algebraic Poincare lemma on `A^2_C` gives a polynomial `H in C[x,y]`
with `theta=dH`.  Taking the degree-four field trace commutes with exterior
differentiation and gives the exact identity

```text
omega_F=Tr_(L/K)(theta)=d Tr_(L/K)(H).                    (31c)
```

Thus every genuine polynomial realization satisfies

```text
[omega_F]=0 in Omega^1_K/dK,

Res_D(omega_F)=0                         for every prime D. (31d)
```

The first promoted version allowed this class to be nonzero for a polynomial
`D_4,A_4`, or `S_4` map.  MISTAKE-301 records the repaired implication: a
nonzero class belongs only to an abstract inverse-spectral PDE solution and
is an obstruction to, rather than a possible invariant of, polynomial
realization.

### 6.1 Explicit coefficient residue gate

Let

```text
p_m=Tr_(L/K)(x^m).
```

Newton's identities give

```text
p_1=-a_3,
p_2=a_3^2-2a_2,
p_3=-a_3^3+3a_3a_2-3a_1,
p_4=a_3^4-4a_3^2a_2+2a_2^2+4a_3a_1-4a_0.                (31e)
```

Since `y=sum_(j=0)^3 b_j x^j`, put

```text
E=sum_(j=0)^3 b_j p_(j+1)=Tr_(L/K)(xy),

C=sum_(j=0)^3 b_j/(j+1) dp_(j+1).                        (31f)
```

Direct differentiation on the four sheets gives

```text
alpha
 =sum_(j=0)^3 [p_(j+1) db_j
               +j/(j+1)b_j dp_(j+1)],

alpha=dE-C,

d alpha=sum_(j=0)^3 1/(j+1) dp_(j+1) wedge db_j
       =4 kappa^(-1)du wedge dv.                          (31g)
```

Combining this with (31c) gives the stronger exact coefficient identity

```text
Gamma_F:=C+4 kappa^(-1)u dv
        =d(E-Tr_(L/K)(H)).                                (31h)
```

Thus `Gamma_F`, not merely its exterior derivative, must be exact.  There is
also a sharp support constraint.  Put

```text
V=A^2 minus A_F=D(c),                 U=F^(-1)(V).
```

The map `U -> V` is finite etale.  Since `H` and `xy` are regular on `U`,
finite-algebra trace gives

```text
G=Tr_(L/K)(H),       E=Tr_(L/K)(xy) in O(V)=C[u,v,1/c],

omega_F=dG,          Gamma_F=d(E-G).                       (31h')
```

Consequently the entire polar principal parts of `alpha`, `omega_F`, and
`Gamma_F` cancel at every affine prime not dividing `c`.  In particular a
pure primitive-discriminant divisor or a pure reconstruction-denominator
divisor cannot support the aggregate trace form, even though individual
`a_i,b_j` terms may have poles there.  Only Jelonek primes can support its
remaining principal parts, and those parts are exact with zero logarithmic
residue.

Because `u dv` is regular at every affine base prime, (31h) imposes the
explicit cancellation law

```text
Res_D sum_(j=0)^3 b_j/(j+1) dp_(j+1)=0
                               for every affine prime D. (31i)
```

More concretely, choose a uniformizer `pi` at `D` and write

```text
b_j=sum_n b_(j,n)pi^n,             p_(j+1)=sum_n p_(j+1,n)pi^n.
```

Then (31i) is the matched-Laurent equation

```text
sum_(j=0)^3 1/(j+1) sum_(r in Z)
       r b_(j,-r) p_(j+1,r)=0.                            (31j)
```

Only finitely many pairs contribute to the residue.  Individual summands in
(31j) depend on the primitive coordinate and companion presentation; their
sum is the lawful exact-class constraint.

At a Jelonek component, (31i)--(31j) must coexist with the full coefficient pole in
(19).  It is a new marked-companion constraint, not a consequence of the
cubic resolvent or of the pole order alone.

### 6.2 The `D_4` gate exists before the quadratic base trace

In the `D_4` lane let `tau` be THM-2612's unique source-deck involution and
put

```text
M=L^tau,                       [L:M]=[M:K]=2,

p_m^tau=Tr_(L/M)(x^m)=x^m+tau(x)^m,

C_tau=sum_(j=0)^3 b_j/(j+1) d p_(j+1)^tau.               (32a)
```

Taking the relative trace of `theta=dH` and repeating the calculation above
gives

```text
Gamma_tau:=C_tau+2 kappa^(-1)u dv
 =d(Tr_(L/M)(xy-H))                         in Omega^1_M. (32b)
```

Hence every divisor of the intermediate quadratic field has zero residue
**before** its two `M/K` conjugates can cancel in the base trace.  This is
strictly finer than (31i).  More precisely, on the finite-etale quotient

```text
W=U/<tau> -> V,
```

the pair potential `Tr_(L/M)(xy-H)` is regular.  Thus its polar support is
confined to the divisors omitted from `W`, not to arbitrary resolvent or
power-basis denominators.  The anti-invariant part is exact as well:

```text
x dy-tau(x)d(tau(y))=d(H-tau(H)).                         (32c)
```

Equations (32b)--(32c) are directly testable against THM-2612's functional
formula for `tau(x)` and `tau(y)=b(tau(x))`.  They also separate two boundary
notions: THM-2612 forces a pole in a deck-coordinate pullback, while (32b)
forces the associated relative trace--Liouville form to have zero
logarithmic residue.  An exact rational differential can have higher-order
principal parts, so these facts are compatible; at the trace--Liouville
level the unavoidable deck pole can appear only through second-kind data.

The exact zero class and the gates (31h), (32b) are invariant under affine determinant-one source and
target coordinate changes.  For the linear source
part

```text
x'=ax+by,                  y'=cx+dy,          ad-bc=1,
```

one has the exact identity

```text
x' dy'-x dy
 =d(ac x^2/2+bc xy+bd y^2/2).                           (32d)
```

The target Liouville term obeys the same identity, and translations add only
exact linear differentials.  Thus the realization gate (31c)--(32c) is
independent of the chosen affine symplectic frame.

For the hostile (22), direct trace gives

```text
alpha=4v du,

d alpha=-4du wedge dv=4 kappa^(-1)du wedge dv,

omega_F=4d(uv).                                          (33)
```

Its trace--Liouville class is zero, as every genuine polynomial realization
must be.  This hostile still separates discriminant and reconstruction poles;
it does not test the new residue obstruction because its class already
vanishes.

### 6.3 A `D_4` PDE hostile with arbitrary sheet-defect residue

The zero-residue condition is not hidden inside (9), the pole law (19), or
even abstract `D_4` monodromy.  On a rational symplectic chart with source
parameters `(s,t)`, put

```text
u=s^4+s^2,                         v=t/(4s^3+2s),

X=1/(s-1),
Y=-t(s-1)^2+lambda(s-1).                              (33a)
```

Direct differentiation gives

```text
dX wedge dY=ds wedge dt=du wedge dv.                   (33b)
```

The primitive coordinate `X` has polynomial eliminant

```text
R=(2-u)X^4+6X^3+7X^2+4X+1=(2-u)f.                     (33c)
```

At `u=2`, the specialized degree is three, so exactly one sheet escapes and
`a_3=6/(2-u)` has the full pole required by (19).  The extension

```text
C(u) subset C(s),                 u=s^4+s^2,             (33d)
```

has `D_4` geometric monodromy: the root-field involution is `s -> -s`, and
adjoining a square root of `-1-s^2` makes the normal closure degree eight.
Equivalently, the squared-pair resolvent is

```text
U(U^2+2U+1+4u),                                          (33e)
```

with one rational root and a generically irreducible quadratic factor.

The companion `Y=b(X)` has **polynomial** base coefficients

```text
b_0=2(-2lambda+10uv-3v),
b_1=-7lambda+26uv-26v,
b_2=2(-3lambda+11uv-16v),
b_3=-(u-2)(-lambda+4uv-6v).                              (33f)
```

Exact reduction modulo (33c) gives

```text
f_v b_u-f_u b_v == f_X,                                  (33g)
```

while

```text
disc_X(R)=-16u(4u+1)^2,
disc_X(f)=-16u(4u+1)^2/(u-2)^6.                          (33h)
```

Thus the sheet-defect divisor is neither a raw discriminant collision nor a
`b`-denominator divisor.  Formula (31f) now gives

```text
alpha=[-20v-lambda/(2-u)]du-4(4u+1)dv,

omega_F=d(-(20u+4)v)+lambda dlog(2-u),

Res_(2-u)(omega_F)=lambda.                               (33i)
```

On this hostile the deck involution has the explicit affine-coordinate form

```text
tau(X)=-X/(2X+1),

tau(Y)=(2X+1)[2lambda-(2X+1)Y].                          (33j)
```

It fixes `u,v`, has `tau^2=1` and Jacobian one, and its pole `2X+1=0` lies
over `u=2`.  The relative pair trace already carries residue `lambda`; the
base trace is the `lambda dlog(2-u)` term in (33i).  At `lambda=0` the deck
pole remains while every logarithmic residue vanishes, giving the sharp
second-kind control predicted by (32b).

The parameter `lambda` is arbitrary.  Therefore the PDE, the exact
one-sheet pole law, polynomial companion coefficients, and `D_4` monodromy
still permit any trace--Liouville residue.  For `lambda!=0`, (31d) rules out
a polynomial `A^2` realization immediately.  For `lambda=0`, the residue
obstruction vanishes but the rational chart (33a) remains nonpolynomial, so
the gate is necessary and not sufficient.

## 7. Resolvent interface and exact residual

THM-2598's matching cubic remembers the three pairings of the four roots and
quotients by `V_4`.  It forgets

```text
the attached values y_i=b(x_i),
the coefficientwise symplectic congruence (9),
the trace--Liouville exactness and coefficient residue gates (31c)--(32c),
and which normalization branches remain in the affine source.             (34)
```

Thus `(b,(9),(31c)--(32c))` is the missing marked-origin sidecar for transferring any
grade-three resolvent anatomy back to a quartic Keller map.  In the `D_4`
lane, it must coexist with THM-2612's birational deck-involution pole and the
present/omitted boundary ownership of THM-2598.  These are compatible
invoices, not yet a contradiction.

The degree-four planar frontier is now the following typed classification
problem: classify rational pairs `(f,b)` satisfying (9), the pole law (19),
the exactness/residue gates (31c)--(32c),
monodromy `D_4,A_4`, or `S_4`, and extension of their total space to `A^2`
with `u,v in C[x,y]`.  This is strictly smaller than an unrestricted planar
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
5. the hostile trace-curvature sign; and
6. the determinant-one Liouville exact-differential identity; and
7. the `D_4` hostile's symplectic Jacobians, eliminant, companion recovery,
   PDE congruence, one-sheet coefficient pole, discriminants, polynomial
   companion coefficients, power-sum trace form, deck involution, and
   arbitrary residue.

The generic-coordinate intersection, resultant factorization, valuation law,
and trace argument are proofs above, not extrapolations from the fixtures.
An independent hostile audit rederived the PDE sign and coefficient formulas,
checked the generic-coordinate and quasi-finite specialization arguments,
verified the `C_4` model and trace signs independently, and enforced the three
separate boundary ledgers in (21).

No polynomial degree-four Keller map, monodromy lane exclusion, local
mechanism forcing the cancellations in (31i)--(32b), JC(2), DC(2), or GMC-to-JC
interface follows.

QED.
