---
id: THM-3841
title: "A deleted ramification divisor with a three-puncture branch forbids plane domination"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  A reusable
  valuation lemma says that if a normal finite completion
  deletes a prime divisor E, then every dominant generically finite plane
  morphism to the complement forces the image of E to be a component of the
  composite's Jelonek set.  For the THM-3811 nonlinear cubic surface this
  image is the irreducible discriminant Delta, whose affine normalization is
  P1 minus three points.  Jelonek--Lason polynomial uniruledness forbids such
  a planar nonproperness component.  Hence the surface admits no dominant
  morphism from A2 at all, and in particular no Keller atlas.  The argument
  is all-degree and does not prove JC(2) outside this completion.
source: root + jc_quartic_c3_construct / deleted-divisor valuation and Jelonek lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / jc-cohn3709, 2026-08-23).  The
  audit rederived extension of the
  divisorial valuation through the finite source-field extension, uniqueness
  of centers on the separated completion, the valuative nonproperness
  implication, component typing, normalization lifting, and the
  polynomial-uniruledness contradiction.  It also repaired two presentation
  gaps: absolute irreducibility now follows geometrically from the
  basepoint-free degree-six birational parametrization, and arbitrary
  characteristic-zero ground fields are handled by explicit finite-coefficient
  descent to C rather than an unexpanded Lefschetz slogan.  The deterministic
  companion verifies rational irreducibility and squarefreeness of Delta, the discriminant and
  different norm, the exact rational parametrization and inverse, the
  basepoint-free projective normalization map, its two projective infinity
  points and three normalization punctures, and the negative valuation
  certificate.  Normal and optimized replay agree with the frozen transcript.
depends_on:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3836-cubic-factor-cofactor-darboux-packet
  - THM-3840-forced-cubic-two-arm-jelonek-passport
citation:
  - "Z. Jelonek and M. Lason, Quantitative properties of the non-properness set of a polynomial map, arXiv:1411.5011v2, Theorems 1.2 and 3.2."
  - "Nguyen Van Chau, Note on the Jacobian condition and the non-proper value set, Ann. Polon. Math. 84 (2004), 203--210, DOI 10.4064/ap84-3-2; arXiv:math/0305088."
script: 04-computation/jc2_deleted_ramification_three_puncture_jelonek_nonentry_thm3841.py
output: 05-knowledge/results/jc2_deleted_ramification_three_puncture_jelonek_nonentry_thm3841.out
script_sha256: f1e7e2a51633119b81efc1d13c875aa07a92222471c81e58e04d2523747a0c87
output_sha256: ea3dc03bd168768236e8cd036d49463709cc0e0ac0dc5b3ae6a649e8cb36f043
semantic_sha256: 9964f46981085706bb9dd2508c6a6c473ee09d2b9410eaf0c50a4c76e6440b29
hash_basis: raw LF bytes
---

# THM-3841 -- a deleted three-puncture branch cannot come from the plane

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `K` of characteristic zero.
Let `Xbar=Spec S`, `E`, and

```text
U=Xbar minus E
```

be the normal cubic completion, its irreducible ramification divisor, and the
affine etale surface of THM-3811.  There is no dominant morphism

```text
psi:A2_(x,y) -> U.                                               (1)
```

In particular, there is no dominant etale plane atlas and no planar Keller
map obtained by composing `(1)` with

```text
pi:U -> A2_(A,C).                                                (2)
```

This is stronger than a support or degree obstruction: it closes every
polynomial plane parametrization of this specific surface.  It does not
exclude another cubic completion and does not prove `JC(2)`.

## 1. Deleted-divisor valuations force Jelonek components

We first isolate the reusable mechanism.  Let

```text
pi_bar:Xbar -> A2
```

be a finite dominant morphism from a normal integral surface, let `E` be a
prime divisor whose image `Gamma=pi_bar(E)` is a divisor, and put
`U=Xbar minus E`.  Suppose that a dominant morphism

```text
psi:A2 -> U                                                     (3)
```

exists.  Then the composite `F=pi_bar psi` is a generically finite polynomial
map and

```text
Gamma is contained in S_F,                                     (4)
```

where `S_F` is its nonproperness set.

Indeed, put `L=K(Xbar)` and `M=K(A2_source)`.  Dominance in equal dimension
makes `M/L` finite.  Let `v_E` be the normalized divisorial valuation of
`L`.  A localization of the integral closure of its DVR in `M` gives an
extension `w` with

```text
w|_L=e v_E,                         e>=1.                        (5)
```

No such `w` has a center on the affine source.  If it had center `p`, then
the morphism `(3)` would make `v_E` centered at `psi(p)` in `U`.  But `v_E`
already has center the generic point of `E` on the separated surface `Xbar`,
and a valuation has at most one center on a separated scheme.  Since
`E` is disjoint from `U`, this is impossible.

The restriction of `w` to `K(A,C)` is centered at the generic point of
`Gamma`, because the finite morphism sends `E` dominantly to `Gamma`.  If
`F` were proper over a neighborhood of that generic point, the valuative
criterion would give `w` a source center, contradicting the preceding
paragraph.  This proves `(4)`.  Since the nonproperness set of a generically
finite polynomial plane map is closed and pure of dimension one when
nonempty, every irreducible `Gamma` in `(4)` is an irreducible component of
`S_F`.

The center argument is intrinsic; a regular function with a pole on `E` is
a useful certificate.  In the THM-3811 completion one has

```text
div_Xbar(A)=P_0+P_1,                div_Xbar(D)=E+P_0,
h=A/D in Gamma(U,O_U).                                          (6)
```

Therefore

```text
div_Xbar(h)=P_1-E,                  v_E(h)=-1.                   (7)
```

Every extension in `(5)` has `w(h)=-e<0`; it cannot be centered on `A2`
because `psi^*h` is a source polynomial.  This independently certifies the
only delicate center step in the present application.

## 2. The whole nonlinear discriminant is forced into `S_F`

For THM-3811, `pi_bar(E)` is the irreducible discriminant curve

```text
Delta=
 A(C+5A)(4C+19A)(3C-17A)
 +C^2(162A^3+126A^2C-4C^3)-27A^2C^4=0.                         (8)
```

The generic point of `E` does not lie over `A=0`: equivalently,
`div_Xbar(A)=P_0+P_1` has no `E` term, and directly

```text
Delta(0,C)=-4C^5.                                               (9)
```

The different `D` has order one on `E`, and its exact norm is

```text
Norm(D)=-A^2 Delta.                                             (10)
```

Thus `(6)--(7)` apply exactly at the generic ramification valuation.  The
deleted-divisor lemma gives

```text
V(Delta) is an irreducible component of S_F.                    (11)
```

This is the global propagation missing from THM-3840: its two explicitly
forced spectral endpoints are points of the component in `(11)`, but one
generic deleted-divisor valuation already forces the entire component.

## 3. Exact normalization and its three punctures

Put

```text
G(q)=3q^2+7,                    R(q)=q^3-7q-6.                  (12)
```

THM-3811 gives the birational parametrization

```text
A(q)=-2q^2R(q)/G(q)^2,          C(q)=-qR(q)/G(q).               (13)
```

For completeness, a rational inverse at the generic point is

```text
q=A(27A-9C^2+7C)/(2(C^2-21A^2)).                               (14)
```

Homogenize with `[u:v] in P1`,

```text
G_h=3u^2+7v^2,                  R_h=u^3-7uv^2-6v^3.
```

The map to the projective closure of `(8)` is

```text
[u:v] |->
[-2u^2R_hv : -uR_hG_h : G_h^2v^2].                            (15)
```

It has no base point.  At `v=0` its middle coordinate is nonzero, while at
`G_h=0` its first coordinate is nonzero because

```text
Res_q(G,R)=6460!=0.                                             (16)
```

The three coordinates in `(15)` are basepoint-free forms of degree six.
Consequently the pullback of a general line has degree six.  Equation `(14)`
makes `(15)` generically one-to-one onto its image, so that image is an
irreducible projective curve of degree six.  The degree-six homogenization of
`Delta` vanishes on the image; it must therefore be the image equation up to
a scalar.  This proves absolute irreducibility as well as showing that `(15)`
is the normalization map.  The affine target chart is the last coordinate
nonzero.
Its inverse image is therefore

```text
v G_h !=0.
```

Consequently the affine normalization of `V(Delta)` is exactly

```text
N=P1 minus ({v=0} union {G_h=0})
 =Spec K[q,G(q)^-1].                                           (17)
```

The quadratic `G` has discriminant `-84`, so `(17)` removes three distinct
normalization points: parameter infinity and the two roots of `G`.  The
projective curve itself has two distinct points at infinity:

```text
v=0       |-> [0:1:0],
G_h=0     |-> [1:0:0].                                         (18)
```

The two roots of `G_h` are separate normalization branches over the second
point.

## 4. Polynomial uniruledness gives the contradiction

Jelonek--Lason prove over `C` that the nonproperness set of a generically
finite polynomial map from affine space is covered by polynomial curves;
their Theorem 3.2 gives a degree bound.  For general algebraically closed
`K` of characteristic zero, all coefficients of a hypothetical `psi`, the
surface equations, and a nonzero Jacobian minor certifying dominance lie in
a finitely generated extension `K_0/Q`.  Embed `K_0` into `C` and base-change.
The same polynomial identities and nonzero minor give a dominant complex
morphism to the same deleted completion, so the complex contradiction below
already excludes the original `psi`.  We may therefore work over `C`.

Choose a general point of the component `(11)`, away from all other
components of `S_F`.  A nonconstant polynomial curve through it must
lie in `V(Delta)` and, since both are curves, is dominant.

The normal source `A1` lifts this map through the normalization `(17)`.  On
coordinate rings it would give

```text
K[q,G(q)^-1] -> K[t],                q |-> p(t).                 (19)
```

Dominance makes `p` nonconstant.  But `(19)` requires `G(p(t))` to be a unit
of `K[t]`, whereas

```text
deg G(p)=2 deg p>0.                                             (20)
```

This contradiction proves that `(3)`, hence `(1)`, cannot exist.

Over `C` there is an independent Keller-specific endpoint.  Nguyen Van Chau
proves that the nonempty nonproper-value set of a nonsingular polynomial map
`C^2 -> C^2` has one point at infinity.  Equation `(18)` gives two already on
the forced component, so a Keller composite is impossible before counting
the third normalization branch.

The proof closes the nonlinear THM-3811 cubic completion in all degree.  It
does not use THM-3838's degree-five floor, does not classify other cubic
completions, and does not prove the planar Jacobian conjecture.  The reusable
design rule is sharper: a deleted completion divisor can support a plane
atlas only if its image is a polynomially uniruled affine curve.  **QED.**
