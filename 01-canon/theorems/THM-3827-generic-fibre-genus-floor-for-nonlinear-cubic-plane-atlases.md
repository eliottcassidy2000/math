---
id: THM-3827
title: "Dual generic-fibre genus floors for nonlinear cubic plane atlases"
status: >
  DUAL GENUS FLOORS AND THE SIX-MEMBER REDUCIBILITY GATE PROVED + CITED +
  VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  THE SECOND G_M ARM AND
  BICHROMATIC PASSPORT ARE PROVISIONAL PROOF CANDIDATES + VERIFIED-EXACT,
  PENDING INDEPENDENT HOSTILE AUDIT.  For the generative closed-polynomial
  factors h=p(g) and k=q(ell) of both pulled-back THM-3811 row functions, the
  smooth projective geometric generic fibres of g and ell each have genus at
  least three.  At equality they are respectively isomorphic after base
  change to explicit degree-eight and degree-seven THM-3822 hyperelliptic
  sidecars, with two versus one points at infinity.  The new sharpenings give
  V_U(k)=G_m as well as V_U(h)=G_m, and say that either h splits or one fixed
  spectral fibre splits into components carrying both square-root signs.
  Existence and uniqueness up to affine change of each generative factor,
  and the equivalence with relative algebraic closedness, are cited from
  Arzhantsev--Petravchuk rather than reproved.  No Jacobian counterexample is
  constructed.
source: jc_quartic_c3_construct / generic-fibre genus reframe, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS for Sections 1--6 (root /
  jc-cohn-boundary, 2026-08-23).  The audit checked the cited source directly:
  Arzhantsev--Petravchuk Proposition 1 and Corollary 1 give existence and
  affine uniqueness of a generative closed polynomial, while Lemma 3 is
  exactly the required relative-algebraic-closure equivalence.  Independently,
  it reconstructed the regular generic curve, the function-field injection
  from a nonconstant sidecar coordinate, proper extension and
  Riemann--Hurwitz, the denominator-valuation proof of R intersect K(g)=K[g],
  and the reduced special-fibre contradiction component by component.  For
  the five-slope strengthening it separately checked algebraic independence
  of h,k from the intrinsic D-quadratic, the second-order h-adic lift,
  pairwise coprimality of the five pencil members, and both surviving subset
  degrees.  No quantifier or scope repair was needed.  INDEPENDENT HOSTILE
  AUDIT PASS for the dual genus floor (jc_sparse_direct_search, 2026-08-23):
  it checked survival of the degree-seven discriminant after q(ell), the
  squarefree genus-three/one-infinity sidecar, both algebraic and
  transcendental branches, Riemann--Hurwitz, equality, normal and optimized
  replay, and the frozen hashes at commit 85eb0017c6.  The second intrinsic
  G_m arm and the strengthened component-sign passport are self-audited and
  await an independent pass.  The exact companion verifies the
  monic degree-eight sidecar, a squarefree hostile fibre, the full generic
  discriminant, stability under a nonconstant Stein composition, the
  genus-three count, every excluded Riemann--Hurwitz genus, and reduced-arm
  controls, the opposite degree-seven sidecar and its full discriminant, the
  complete k-zero Laurent parametrization and original cubic laws, the second
  completed square, the component-sign equation, the five distinct pencil
  slopes, and both terminal subset obstructions.  Normal and optimized replay
  agree with the frozen output.
depends_on:
  - THM-3822-nonlinear-cubic-plane-atlas-sl2-and-punctured-arm-gate
related:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
citation:
  - "Arzhantsev--Petravchuk, Closed and Irreducible Polynomials in Several Variables, arXiv:math/0608157v2, Proposition 1, Lemma 3, and Corollary 1."
script: 04-computation/jc2_nonlinear_cubic_atlas_generic_fibre_genus_floor_thm3827.py
output: 05-knowledge/results/jc2_nonlinear_cubic_atlas_generic_fibre_genus_floor_thm3827.out
script_sha256: dd6025305fd479c6118d57409e390ffa9751e15da533f465b50a998c59d45722
output_sha256: a59b2c098e2cd115adb6da15dd711a4e6e2752c1d15b8ed5378bcc589012e86d
semantic_sha256: d46b81bd79d27d7b4396f8f934dd773e694c0eb1a6384ffe03381e6934db7a29
hash_basis: raw LF bytes
---

# THM-3827 -- both row fibrations of a plane atlas need genus at least three

**DUAL GENUS FLOORS AND SIX-MEMBER REDUCIBILITY PROVED + CITED +
VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED; SECOND G_M ARM AND
BICHROMATIC PASSPORT PROVISIONAL PROOF CANDIDATES + VERIFIED-EXACT, PENDING
INDEPENDENT HOSTILE AUDIT.**  Work over an algebraically closed field `K` of
characteristic zero.  Put

```text
R=K[x,y],                 L=K(x,y).                              (1)
```

Recall the THM-3822 polynomial

```text
H(T,Z)=
 Z^8-14T^2Z^6+12T^2Z^5-12T^3Z^5+49T^4Z^4
 +112T^4Z^3+84T^5Z^3+36T^5Z^2+36T^6Z^2
 +196T^6Z+84T^7.                                               (2)
```

## 1. The function-field genus gate

Let `g in R` be a nonconstant **closed polynomial**, by which this theorem
means exactly

```text
K(g) is relatively algebraically closed in L.                   (3)
```

Let `Gamma_g` be the smooth projective curve over `K(g)` with function
field `L`, and let `gamma` be its geometric genus.  Condition `(3)` and
characteristic zero make `L/K(g)` regular, so `Gamma_g` is geometrically
integral.

Let `p in K[S]` be nonconstant and set

```text
h=p(g).                                                         (4)
```

Suppose there are `k,w in R` with

```text
H(h,k)=w^2.                                                     (5)
```

Then exactly one of the following holds:

```text
(a) k in K[g]; or
(b) gamma>=3, and over an algebraic closure of K(g) there is a
    nonconstant map Gamma_g -> C_h of degree d satisfying

        2gamma-2 = 4d + ramification_degree,                    (6)

    where C_h is the genus-three curve W^2=H(p(g),Z).           (7)
```

In particular, if `gamma<=2`, then `k in K[g]`.  If `gamma=3` and
`k notin K[g]`, `(6)` forces `d=1` and zero ramification, so the two smooth
projective geometric curves are isomorphic.

## 2. The sidecar really has genus three after every Stein composition

As a polynomial in `Z`, `H(T,Z)` is monic of degree eight.  Its discriminant
is

```text
Disc_Z H(T,Z)=19144454963200 T^44 G(T),                          (8)

G(T)=
 466560000T^8-96264205056T^7-180219652272T^6
 -35458022604096T^5-85186647741177T^4-616061164068T^3
 -497726825264T^2+2712579840T-678144960.                        (9)
```

This is nonzero.  A smaller hostile certificate is the fibre `T=-1`:

```text
H(-1,Z)=
 Z^8-14Z^6+24Z^5+49Z^4+28Z^3+196Z-84,                          (10)
```

whose gcd with its derivative is one.  Since `g` is transcendental over
`K` and `p` is nonconstant, `p(g)` is transcendental.  Therefore `(8)`
remains nonzero after `T=p(g)`.  The right side of `(7)` is a squarefree
polynomial of degree eight, so `C_h` is geometrically integral and has
genus `(8-2)/2=3`.

## 3. Proof of the dichotomy

Put `F=K(g)`.  Equation `(5)` gives an `L`-point

```text
(Z,W)=(k,w)                                                       (11)
```

on the affine part of `C_h`.  If `k` is algebraic over `F`, relative
algebraic closedness `(3)` gives `k in F`.

Otherwise `k` is transcendental over `F`.  Then `(11)` induces an injection
of the function field `F(C_h)` into `L`: a nonzero kernel in the affine
curve ring would have finite residue field over `F` and would make `k`
algebraic.  Thus `(11)` is a nonconstant rational map

```text
Gamma_g --> C_h.                                                (12)
```

After extending constants to an algebraic closure, properness extends
`(12)` across every missing point.  In characteristic zero it is separable,
and Riemann--Hurwitz gives `(6)`.  Because `C_h` has genus three, its main
term is `d(2*3-2)=4d`.  This proves the genus floor and the genus-three
equality statement.

It remains to turn `k in F` into a polynomial statement.  The elementary
intersection needed here is

```text
R intersect K(g)=K[g].                                         (13)
```

Indeed, write an element of the left side as `a(g)/b(g)` with coprime
`a,b in K[S]`.  If `b` is nonconstant, choose a root `alpha` of `b`.
Then `a(alpha)!=0`.  For any irreducible factor `q` of `g-alpha`, the
height-one valuation `v_q` is positive on `b(g)` and zero on `a(g)`, a
contradiction to membership in `R`.  Hence `b` is constant, proving `(13)`
and part `(a)`.

## 4. Consequence for a plane atlas of the nonlinear cubic surface

Let

```text
psi:A2_(x,y) -> U                                               (14)
```

be a dominant etale morphism to the THM-3811 affine cubic surface, and
write `h,k` for the pullbacks of its intrinsic functions `A/D,omega/D`.
THM-3822 gives a polynomial `w` satisfying `(5)`.  It also proves

```text
V_U(h)=G_m,                                                     (15)
```

and says that every irreducible component of the etale base change
`V_A2(h)` carries `k` as a nonconstant unit.

Arzhantsev--Petravchuk, Proposition 1 and Corollary 1, prove that the actual
nonconstant pulled-back polynomial has a generative factorization

```text
h=p(g)                                                          (16)
```

with `g` closed; this factor is unique up to affine change.  Their Lemma 3
identifies closedness with the precise relative-field condition `(3)`.
Therefore the genus `gamma` of the primitive generic fibre of `g` satisfies

```text
gamma>=3.                                                       (17)
```

For suppose `gamma<=2`.  Section 1 and `(13)` give `k=q(g)` for some
`q in K[S]`.  Meanwhile `(15)` is a reduced divisor and etale base change
preserves reducedness, so `p(g)` is squarefree in `R`.  In particular `p`
has no repeated root; otherwise `(g-alpha)^2` would divide `p(g)`.

Choose any root `alpha` of `p`, and any irreducible component of
`g-alpha`.  It is an irreducible component of `V(h)`, but on it

```text
k=q(alpha),                                                     (18)
```

a scalar.  This contradicts the nonconstant-unit conclusion following
`(15)`.  Hence `(17)` holds.

## 5. Genus equality boundary

The theorem replaces the vague phrase "multi-ended arm" by a sharper
birational condition:

```text
Every closed-polynomial Stein coordinate underlying a plane atlas has
geometric generic-fibre genus at least three.                   (19)
```

Thus coordinate, rational, elliptic, and genus-two fibrations are all
excluded in one argument.  At the first surviving value `gamma=3`, the
generic fibre must be geometrically isomorphic to the explicit
hyperelliptic curve `(7)`.  If the map has degree `d>=2`, `(6)` already
forces `gamma>=2d+1>=5`.

## 6. A six-member reducible-and-bichromatic packet

The same square sidecar has a second exact completion.  Put

```text
A_5=(7h^2+3k^2)(3h^3+7h^2k+k^3),
B_3=(h+k)(2h+k)(3h-k).                                         (20)
```

Direct expansion gives

```text
H(h,k)=(kB_3)^2+4h^2A_5.                                      (21)
```

First note that the pulled-back `h,k` are algebraically independent.  In
the intrinsic ring, `h` is not a unit because `V_U(h)=G_m`.  If
`7h^2+3k^2=0`, algebraic closedness and the determinant-one law would give
`k=lambda h` and then `1=h(lambda C-m)`, a contradiction.  Thus the
THM-3822 quadratic in `D` makes `D` algebraic over `K(h,k)`.  Its linear
reconstruction formulas then put `C,m,A,omega,theta` in an algebraic
extension of `K(h,k)`.  Since `K(U)` has transcendence degree two, so does
`K(h,k)`.  Dominance of the plane atlas preserves this independence.

Consider the six pencil members

```text
h=0,                  h-alpha_i k=0  (1<=i<=5),                 (22)
```

where the `alpha_i` are the roots of

```text
a(z)=(7z^2+3)(3z^3+7z^2+1).                                   (23)
```

They are five distinct finite slopes: `disc(a)=353831803500!=0`.  The exact
surviving passport is stronger than mere reducibility:

```text
either h is reducible, or some h-alpha_i k is reducible and its prime
components carry both signs w=+kB_3 and w=-kB_3 generically.             (24)
```

Call the latter property **bichromatic**.  To prove `(24)`, suppose `h` is
irreducible.  Modulo `(h)`, the arm law makes `k` a unit and `B_3=-k^3`, so
`kB_3=-k^4` is a unit.  From `w^2=H(h,k)` and `(21)`, one of
`w-kB_3,w+kB_3` vanishes modulo `h`.  Changing the sign of `w` if necessary
and taking one more `h`-adic residue gives

```text
w-kB_3=2h^2d,                A_5=d(kB_3+h^2d)                  (25)
```

for some `d in K[x,y]`.  Indeed, after writing `w-kB_3=hc`, cancellation
in `(w-kB_3)(w+kB_3)=4h^2A_5` and reduction modulo `h` force `h|c`;
division by the unit `2` gives `(25)`.

The determinant-one identity makes `gcd(h,k)=1`.  Hence any two distinct
pencil members `h-alpha k`, including `h` and the three members occurring in
`kB_3`, are coprime after pullback.  The exact resultant in `(28)` below and
`A_5(1,0)=21` therefore show

```text
gcd(A_5,hkB_3)=1.                                               (26)
```

Let `r` now be a prime factor of one fixed member `h-alpha_i k`.  At its
generic point, `h` and `kB_3` are units, while `(21)` gives exactly one of
the signs `w=+kB_3` or `w=-kB_3`.  Taking the `r`-valuation in `(25)` shows
that `d` contains the full multiplicity of `r` in `A_5` in the plus case,
and contains no `r` in the minus case.

If all five fixed-slope fibres were monochromatic, every prime component of
each whole member would make the same choice.  Thus UFD factorization and
`(25)` would make `d`, up to a scalar, the product of a subset of the five
**whole members**, even when those members are reducible or nonreduced.  Let
the size of this subset be `r_0`.

The formal grading in the algebraically independent pair `(h,k)` turns
`A_5=d kB_3+h^2d^2` into degrees

```text
5 on the left,              r_0+4 and 2r_0+2 on the right.      (27)
```

Every `r_0` except `1,2` has an unmatched extreme degree.  For `r_0=1`, the
degree-five part would force `A_5=d kB_3`, but

```text
Res_z(A_5(z,1),(kB_3)(z,1))=-31298700!=0,                      (28)
```

and `A_5(1,0)=21`, so `kB_3` shares no projective linear factor with `A_5`.
For `r_0=2`, the degree-six part would force `h^2d=-kB_3`, impossible because
`gcd(h^2,kB_3)=1`.  Thus not all five fibres are monochromatic.  One contains
prime components of both signs, so it is reducible and bichromatic.  This
proves `(24)` and, in particular, the earlier reducible-fibre claim.

Thus the surviving construction lane is narrower than a generic high-genus
pencil: unless `h` itself splits, it must allocate components of both square-
root signs inside one fixed spectral member before solving the second-row and
Keller equations.

## 7. The opposite row has the same floor and opposite infinity parity

There is first an intrinsic geometric reason that the opposite row is not a
formal symmetry.  Set `k=0` in the actual THM-3811 ring `B`.  The three
reconstruction identities `(3),(5),(6)` of THM-3822 force

```text
h=3/(7C^2),       m=-7C^2/3,       D=7C^4/9,                    (29)
A=C^2/3,          omega=0,          theta=-7C^2/3.
```

Conversely these assignments satisfy all three original cubic multiplication
laws, reconstruct the displayed nonzero `D`, and hence extend from `S` to
`B=S[A/D,omega/D]`.  Moreover `hC^2=3/7` makes `C` a unit in `B/(k)`, and
the formulas `(29)` express every generator of that quotient in
`K[C,C^-1]`.  The two maps are inverse, so

```text
B/(k)=K[C,C^-1],                 V_U(k) isomorphic to G_m.       (30)
```

Thus an etale plane atlas pulls both intrinsic arms back etale.  On every
component of `V(h)`, the function `k` is a nonconstant unit by THM-3822; on
every component of `V(k)`, the function `C`, and hence
`h=3/(7C^2)`, is a nonconstant unit.  This is the geometric source of the
even/odd asymmetry below.  Nonconstancy here is automatic: an etale morphism
from a curve component to `G_m` cannot have zero-dimensional image.

Apply the same cited generative-factor theorem to the other pulled-back row
entry:

```text
k=q(ell),                                                       (31)
```

where `ell in R` is closed and `q in K[S]` is nonconstant.  As a polynomial
in its first argument, the sidecar has degree seven and constant leading
coefficient `84`.  Its exact discriminant is

```text
Disc_T H(T,Z)=-683730534400 Z^43 J(Z),                          (32)

J(Z)=6480000Z^8+1952190576Z^7-14515170152Z^6
     +76957426508Z^5+405669771962Z^4-55140029819Z^3
     -17308754768Z^2+1276289280Z+234420480.                    (33)
```

This is nonzero.  Since `q(ell)` is transcendental, the polynomial
`H(T,q(ell))` is squarefree of degree seven in `T`; its smooth projective
hyperelliptic curve has genus three and one geometric point at infinity.

Let `Gamma_ell` be the smooth projective generic fibre of `ell`.  The point
`(T,W)=(h,w)` on this odd sidecar gives exactly the dichotomy of Section 1:
either `h` is transcendental over `K(ell)` and hence induces a nonconstant
map from `Gamma_ell` to a genus-three curve, or `h` is algebraic over
`K(ell)`.  In the second case, closedness of `ell` gives `h in K(ell)`, and
the intersection argument `(13)` gives `h in K[ell]`.  But `(31)` then makes
`h,k` algebraically dependent, contradicting the independence proved in
Section 6.  Therefore

```text
genus(Gamma_ell)>=3.                                           (34)
```

There is also a componentwise proof of the algebraic branch that exactly
mirrors Section 4.  Reducedness of the intrinsic divisor `(30)` and etale
base change make `q(ell)` squarefree.  Choose a root `beta` of `q` and a
component of `ell-beta`.  If `h in K[ell]`, then `h` is scalar on that
component, contradicting the nonconstant-unit law following `(30)`.

At equality, Riemann--Hurwitz again forces the map to have degree one, so
the generic fibre is geometrically isomorphic to the odd sidecar.  The
degree-eight sidecar from Sections 1--5 has two geometric points at infinity,
whereas this degree-seven sidecar has one.  Thus a surviving atlas must carry
two primitive high-genus fibrations with opposite infinity parity, not merely
one complicated arm.

## 8. Exact scope

The generative-polynomial existence and uniqueness used here are **CITED**,
not reproved by the exact companion.  They make the independently audited
dual claims `(19)` and `(34)` unconditional for every dominant etale plane
atlas over the stated field.  The second intrinsic arm `(30)` and the
bichromatic strengthening `(24)` are new provisional refinements pending an
independent hostile pass.  No planar Jacobian counterexample is claimed.
**QED for the dual genus floors and reducibility gate; the two sharpenings
remain subject to independent hostile audit.**
