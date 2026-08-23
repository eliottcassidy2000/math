---
id: THM-3843
title: "Russell-arm Darboux restrictions are immersed normalizations and forced Jelonek components"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every Darboux
  pair on the THM-3785 Russell pseudo-plane restricts
  to a finite birational normalization with nowhere-vanishing differential
  on its unique multiple arm.  The explicit
  cubic plane atlas identifies that arm with a source line.  If the arm map
  were injective, Gwozdziewicz's injective-line theorem would make the
  composite planar Keller map an automorphism, contradicting its field degree
  at least nine.  Hence two distinct smooth arm points have the same image,
  which is a singular point of a forced Jelonek component.  Existence of an
  arbitrary Darboux pair remains open.
source: root + jc_quartic_c3_construct / Russell pseudo-plane arm lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / jc-cohn3709, 2026-08-23).  The audit
  rederived the intrinsic arm jet,
  polynomial Luroth step, finiteness and normalization, source-line transfer,
  field-degree tower, Gwozdziewicz hypothesis, singularity conclusion, and
  Jelonek-component transfer.  This strengthens THM-3790 from differential
  nonvanishing and noninjectivity to an exact normalization statement.
  The companion verifies the Russell relation and Poisson packet, exact line
  to arm isomorphism, universal derivative Bezout identity, common-parameter
  derivative hostile, and a nodal birational-immersion positive control.  It
  also checks the unique pole place and its unique extension in the
  polynomial-Luroth step, the finite normalization inference, the exact
  field-degree tower, and an explicit Laurent escaping arc over every arm
  point for the Jelonek transfer.  Normal and optimized replay agree with the
  frozen transcript.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3790-cubic-pseudoplane-arm-nodal-immersion-gate
related:
  - THM-3787-cubic-pseudoplane-complete-low-support-darboux-nonentry
  - THM-3794-degree-two-etale-map-constant-unit-obstruction
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
citation:
  - "J. Gwozdziewicz, Injectivity on one line, arXiv:alg-geom/9305008, Theorem 1.1."
  - "Z. Jelonek and M. Lason, Quantitative properties of the non-properness set of a polynomial map, arXiv:1411.5011."
script: 04-computation/jc2_russell_arm_birational_immersion_self_identification_thm3843.py
output: 05-knowledge/results/jc2_russell_arm_birational_immersion_self_identification_thm3843.out
script_sha256: 0b32d9962c59a122e826afce3725984dd6c9a1d13db19dd46995e7abcb842b22
output_sha256: eea6f9a8ca95aeb4ee44bce32500722e242bade04a1aac2f71e2570b049d8d0f
semantic_sha256: a24dc6e0b848766ce63799ad0c879a463f1e93e8c2162d28694f6dbebad95c7a
hash_basis: raw LF bytes
---

# THM-3843 -- the Russell arm must self-identify

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `K` of characteristic zero,
fix `c in K*`, and let

```text
Y=Spec B,
B=K[r,z,e]/(r^2e-z^3+c^3r)                                    (1)
```

be the smooth Russell pseudo-plane of THM-3785.  Its Poisson packet is

```text
{r,z}=3r^2,             {r,e}=9z^2,
{z,e}=3c^3+6re.                                               (2)
```

Let `P,Q in B` be any Darboux pair,

```text
{P,Q}=lambda in K*.                                           (3)
```

On the unique multiple arm

```text
L=V(r,z)=A1_e                                                 (4)
```

write

```text
p(e)=P|_L,                    q(e)=Q|_L.                       (5)
```

Then the restriction

```text
gamma=(p,q):L -> A2_(P,Q)                                     (6)
```

is the finite normalization of its image, its differential is nowhere zero,
but it is not injective.  Hence there are
distinct `e_1,e_2 in K` with

```text
(p(e_1),q(e_1))=(p(e_2),q(e_2)).                              (7)
```

Their common image is a singular point with at least two smooth normalization
branches on the polynomial curve `Gamma=gamma(L)`.  Moreover `Gamma` is an
irreducible component of the nonproperness set of the planar Keller map
obtained from the THM-3785 cubic source atlas.

The theorem is all-degree and applies to every hypothetical Darboux pair on
this fixed surface.  It does not construct such a pair and does not close
the remaining noncentral multigraded support.

This is a strict strengthening of THM-3790: that theorem already proves the
nowhere-zero differential and noninjectivity and develops the first nodal
support gate.  The new step here is the polynomial-Luroth argument proving
that the arm itself is the normalization, followed by the forced
Jelonek-component transfer.

## 1. The Darboux equation is a derivative Bezout identity on the arm

Choose polynomial representatives of `P,Q` in `K[r,z,e]`.  Their `z`-jets
on `L` are intrinsic: changing a representative by a multiple of the
relation in `(1)` changes its `z`-derivative by a multiple of `-3z^2`, which
vanishes on `L`.  Put

```text
alpha(e)=P_z|_L,                    beta(e)=Q_z|_L.             (8)
```

The chain-rule expansion of `(2)` is

```text
{P,Q}=
 3r^2(P_rQ_z-P_zQ_r)
 +9z^2(P_rQ_e-P_eQ_r)
 +(3c^3+6re)(P_zQ_e-P_eQ_z).                                  (9)
```

Restricting `(9)` to `(4)` gives the exact one-variable identity

```text
3c^3[alpha q'-p' beta]=lambda.                                (10)
```

Thus `(p',q')=K[e]`.  In particular the tangent vector `(p',q')` never
vanishes, so `(6)` has everywhere-injective differential (it is a
*differential immersion*) and at least one of `p,q` is nonconstant.  We do
not use “immersion” here in the scheme-theoretic locally-closed-embedding
sense, which would be incompatible with the noninjectivity proved below.

## 2. Polynomial Luroth makes the arm map birational

Suppose, toward a contradiction, that

```text
K(p,q) is a proper subfield of K(e).                            (11)
```

Luroth's theorem gives `K(p,q)=K(s)` for a rational function `s(e)`.  The
polynomial form needed here follows from places.  Let `E=K(p,q)`.  Any pole
place of either `p` or `q` in `E` has a pole above it in `K(e)`, while both
polynomials have only the source place at infinity as a pole.  Since at
least one is nonconstant, `E` therefore has exactly one pole place
`infinity_0` for the pair, and source infinity is its unique extension to
`K(e)`.  A degree-one place of a rational function field can be chosen as
the pole of a Luroth coordinate.  Choose `s` with its only pole at
`infinity_0`.  Its pullback has no finite pole, so

```text
s in K[e].                                                      (12)
```

The same pole argument inside `K(s)` gives

```text
p=P_0(s),                    q=Q_0(s)                            (13)
```

for one-variable polynomials `P_0,Q_0`.  The degree of `s` equals
`[K(e):K(p,q)]`; under `(11)` it is at least two.  Characteristic zero then
makes `s'` a nonconstant common factor of

```text
p'=P_0'(s)s',                    q'=Q_0'(s)s',                  (14)
```

contradicting `(10)`.  Therefore

```text
K(p,q)=K(e).                                                     (15)
```

The map is also finite.  If, say, `p` is nonconstant, then `e` is integral
over `K[p]`: divide the relation `p(T)-p(e)=0` by the nonzero leading
coefficient of `p`.  Hence `K[e]` is finite over `K[p,q]`.  Equations
`(15)` and normality of `K[e]` show that `(6)` is exactly the finite
normalization map of its image curve `Gamma`.  Together with `(10)`, it is an
immersed normalization: finite and birational, with nowhere-zero
differential.

## 3. The source line forbids injectivity

THM-3785 supplies the surjective etale cubic atlas

```text
phi:A2_(x,y) -> Y,

r=x^3,
z=x(c+x^3y),
e=3c^2y+3cx^3y^2+x^6y^3.                                     (16)
```

On the source line `ell=V(x)`, formula `(16)` becomes

```text
ell -> L,                    y |-> e=3c^2y,                     (17)
```

an isomorphism.  The composite

```text
F=(P,Q) phi:A2 -> A2                                           (18)
```

has constant Jacobian `lambda`.  Put

```text
d=[Frac B:K(P,Q)].                                              (19)
```

THM-3785 proves `d>=3`, while the atlas field degree is three, so

```text
[K(x,y):K(P,Q)]=3d>=9.                                        (20)
```

If `(6)` were injective, then `(17)` would make the planar Keller map `(18)`
injective on the affine line `ell`.  Gwozdziewicz's Theorem 1.1 applies over
every algebraically closed characteristic-zero field and says that such a
Keller map is a polynomial automorphism.  This contradicts `(20)`.  Hence
`(7)` is forced.

Because `(6)` is the normalization map, a smooth point of `Gamma` has only
one normalization point above it.  The two distinct points in `(7)` therefore
have a singular common image.  Identity `(10)` says the normalization has
nonzero differential at both points, so this is a genuine multi-branch
self-identification, not a cuspidal failure of the arm differential.

## 4. The self-identified curve is a forced Jelonek component

The cubic atlas `(16)` is nonproper over every point of `L`: one Kummer lift
stays finite and the other two escape.  This can be seen without a properness
slogan.  Fix `(0,0,e_0) in L`, choose a nontrivial cube root of unity `xi`,
and put

```text
w(t)=xi c+[e_0/(3xi^2c^2)]t^3,
x=t,                    y=(w(t)-c)/t^3.                         (20a)
```

As `t->0`, the source Laurent arc escapes because `y` has a pole, while
`r=t^3 ->0`, `z=tw(t)->0`, and

```text
e=(w(t)^3-c^3)/t^3 -> e_0.                                    (20b)
```

Composing this arc with the regular morphism `(P,Q)` shows

```text
Gamma is contained in S_F.                                     (21)
```

The map `(6)` is finite, so `Gamma` is a closed irreducible curve.  The
nonproperness set of the generically finite polynomial map `(18)` is a curve
(or is empty); this is the standard Jelonek nonproperness theorem, and it is
nonempty by `(21)`.  Consequently the closed irreducible curve `Gamma` is an
irreducible Jelonek component.
Its polynomial normalization gives the required one place at infinity, but
Sections 1--3 show that the finite part must contain a singular
self-identification.

The field scope is literal.  Gwozdziewicz's Theorem 1.1 is stated over an
arbitrary algebraically closed field of characteristic zero.  For any
complex-geometric formulation of the nonproperness theorem, all coefficients
and identities above descend to a finitely generated characteristic-zero
subfield, which embeds in `C`; nonproperness and the dimensions of its closed
components commute with the resulting faithfully flat base change.  Thus the
Jelonek-component conclusion descends back to `K`; no uncountable-field
embedding is being assumed.

Thus the first surviving Russell-pseudo-plane counterexample cell has an
exact geometric obligation beyond degree and Euler support: its unique arm
must be drawn as the immersed normalization of a noninjective polynomial
curve inside the target Jelonek set.  **QED.**
