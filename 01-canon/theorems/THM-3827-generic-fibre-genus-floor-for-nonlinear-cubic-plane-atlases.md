---
id: THM-3827
title: "Generic-fibre genus floor for nonlinear cubic plane atlases"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
  AUDIT.  For the generative closed-polynomial factor h=p(g) of the
  pulled-back THM-3811 arm function, the smooth projective geometric generic
  fibre of g has genus at least three.  In genus three it must be isomorphic
  to the explicit THM-3822 hyperelliptic sidecar after base change.  Existence
  and uniqueness up to affine change of g, and the equivalence with relative
  algebraic closedness of K(g), are cited from Arzhantsev--Petravchuk rather
  than reproved.  Independently, among h=0 and five explicit fixed-slope
  members h-alpha_i*k=0 of the intrinsic pencil, at least one is reducible.
  No Jacobian counterexample is constructed.
source: jc_quartic_c3_construct / generic-fibre genus reframe, 2026-08-23
audit: >
  SELF-AUDITED proof candidate.  The 22-gate exact companion verifies the
  monic degree-eight sidecar, a squarefree hostile fibre, the full generic
  discriminant, stability under a nonconstant Stein composition, the
  genus-three count, every excluded Riemann--Hurwitz genus, and reduced-arm
  controls, the second completed square, the five distinct pencil slopes,
  and both terminal subset obstructions.  Normal and optimized replay agree
  with the frozen output.
  Independent audit of the cited generative-factor import, the
  relative-constant-field step, and the special-fibre quantifiers is still
  required before promotion to PROVED.
depends_on:
  - THM-3822-nonlinear-cubic-plane-atlas-sl2-and-punctured-arm-gate
related:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
citation:
  - "Arzhantsev--Petravchuk, Closed and Irreducible Polynomials in Several Variables, arXiv:math/0608157v2, Proposition 1, Lemma 3, and Corollary 1."
script: 04-computation/jc2_nonlinear_cubic_atlas_generic_fibre_genus_floor_thm3827.py
output: 05-knowledge/results/jc2_nonlinear_cubic_atlas_generic_fibre_genus_floor_thm3827.out
script_sha256: 82290ea8abca925913393c950e28dc920c20517ad8156a6aca5b74342efb3516
output_sha256: f25dcb0dc688e67e8ce5b8e3e5eba71da7df5bc7290ffe2210fd07e4f5435450
semantic_sha256: 42bfbc4e18afc0ff86b28847525f1c9f5504331abdf1061f323c79b0d31043bd
hash_basis: raw LF bytes
---

# THM-3827 -- a plane atlas needs generic-fibre genus at least three

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; PENDING INDEPENDENT
HOSTILE AUDIT.**  Work over an algebraically closed field `K` of
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

## 6. A six-member reducible-fibre packet

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

They are five distinct finite slopes: `disc(a)=353831803500!=0`.
At least one member in `(22)` is reducible in `K[x,y]`.

Suppose otherwise.  In particular `h` is irreducible.  Modulo `(h)`, the
arm law makes `k` a unit and `B_3=-k^3`, so `kB_3=-k^4` is a unit.  From
`w^2=H(h,k)` and `(21)`, one of `w-kB_3,w+kB_3` vanishes modulo `h`.
Changing the sign of `w` if necessary and taking one more `h`-adic residue
gives

```text
w-kB_3=2h^2d,                A_5=d(kB_3+h^2d)                  (24)
```

for some `d in K[x,y]`.  Indeed, after writing `w-kB_3=hc`, cancellation
in `(w-kB_3)(w+kB_3)=4h^2A_5` and reduction modulo `h` force `h|c`;
division by the unit `2` gives `(24)`.

Over `K`, equation `(23)` factors `A_5` into the five pairwise distinct
linear binary forms `h-alpha_i k`.  By the assumed irreducibility of their
pullbacks, UFD factorization in `K[x,y]` and `(24)` make `d`, up to a scalar,
the product of a subset of these five forms.  Let its size be `r`.

The formal grading in the algebraically independent pair `(h,k)` turns
`A_5=d kB_3+h^2d^2` into degrees

```text
5 on the left,                 r+4 and 2r+2 on the right.       (25)
```

Every `r` except `1,2` has an unmatched extreme degree.  For `r=1`, the
degree-five part would force `A_5=d kB_3`, but

```text
Res_z(A_5(z,1),(kB_3)(z,1))=-31298700!=0,                      (26)
```

and `A_5(1,0)=21`, so `kB_3` shares no projective linear factor with `A_5`.
For `r=2`, the degree-six part would force `h^2d=-kB_3`, impossible because
`gcd(h^2,kB_3)=1`.  This contradiction proves the reducible-fibre claim.

Thus the surviving construction lane is narrower than a generic high-genus
pencil: it must also allocate components at one of six fixed spectral
members before solving the second-row and Keller equations.

## 7. Exact scope

The generative-polynomial existence and uniqueness used here are **CITED**,
not reproved by the exact companion.  They make `(19)` unconditional for
every dominant etale plane atlas over the stated field, while the function-
field argument of Sections 1--4 is self-contained once `h=p(g)` and `(3)`
are supplied.  No planar Jacobian counterexample is claimed.  **QED, subject
to independent hostile audit.**
