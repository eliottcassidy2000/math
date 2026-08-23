---
id: THM-3827
title: "Generic-fibre genus floor for nonlinear cubic plane atlases"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
  AUDIT.  Let h=p(g) be the pulled-back THM-3811 arm function, where g is
  a closed polynomial: K(g) is relatively algebraically closed in K(x,y).
  If a dominant etale plane atlas exists, the smooth projective geometric
  generic fibre of g has genus at least three.  In genus three it must be
  isomorphic to the explicit THM-3822 hyperelliptic sidecar after base
  change.  The theorem assumes the displayed closed-polynomial factor; it
  does not assert a universal polynomial Stein-factor theorem and does not
  construct a Jacobian counterexample.
source: jc_quartic_c3_construct / generic-fibre genus reframe, 2026-08-23
audit: >
  SELF-AUDITED proof candidate.  The 15-gate exact companion verifies the
  monic degree-eight sidecar, a squarefree hostile fibre, the full generic
  discriminant, stability under a nonconstant Stein composition, the
  genus-three count, every excluded Riemann--Hurwitz genus, and reduced-arm
  controls.  Normal and optimized replay agree with the frozen output.
  Independent audit of the relative-constant-field and special-fibre
  quantifiers is still required before promotion to PROVED.
depends_on:
  - THM-3822-nonlinear-cubic-plane-atlas-sl2-and-punctured-arm-gate
related:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
script: 04-computation/jc2_nonlinear_cubic_atlas_generic_fibre_genus_floor_thm3827.py
output: 05-knowledge/results/jc2_nonlinear_cubic_atlas_generic_fibre_genus_floor_thm3827.out
script_sha256: 5e91767d38f4e7ff0f4f31b00293e61b1e9b5daee2af65482c25dd46f2c7cd31
output_sha256: 662991a4bcdd2e59488253127ccb262f286b293be61b8896be7fb1902398759f
semantic_sha256: fe01574799af6d0f72d21666698aaca59722261dbcbf8491ed758a6c0e196c7f
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

Assume the actual pulled-back polynomial admits a displayed factorization

```text
h=p(g)                                                          (16)
```

with `g` closed in the precise sense `(3)`.  Then the genus `gamma` of the
primitive generic fibre of `g` satisfies

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

## 5. Equality boundary and exact scope

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

The qualification in `(19)` is load-bearing.  This theorem does **not**
prove that every polynomial `h` has a polynomial decomposition `h=p(g)`
whose `g` satisfies `(3)`.  It applies whenever such a closed polynomial
factor is supplied, including the important case that `h` itself is
closed (`p(S)=S`).  No planar Jacobian counterexample is claimed.  **QED,
subject to independent hostile audit.**
