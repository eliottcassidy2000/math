---
id: THM-2696
title: "Reflection-completed S4 relative different and coordinate-invariant Jacobian gate"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For any proper inclusion H<G of finite complex reflection groups,
  the quotient between their polynomial invariant spaces has relative
  Jacobian J_G/J_H, the product of the new reflection-hyperplane factors;
  hence it is never Keller, in any polynomial basic invariants.  For the
  even-sign-change semidirect product V4 semidirect S3=S4, the quotient map
  from the S3 source quotient to the S4 target quotient has intrinsic Kahler different
  (s_1 s_2-s_3), equal to the missing plus-reflection divisor.  Its graph
  quartic and squared-pair resolvent have identical discriminants, and the
  pullback splits as the source S3 discriminant times this relative Jacobian
  square.  Every polynomial generator change, and even every dominant
  polynomial sandwich, retains the nonconstant factor.  The sharp slice
  s_1 s_2-s_3=c!=0 has source A2 and unit relative Jacobian, but its
  coefficient image is an explicit singular nonnormal surface whose
  normalization is that A2.  This excludes only the fixed
  reflection-completed quotient family and its polynomial coordinate
  variants; no general S4, JC(2), or DC(2) conclusion follows.
source: a4-resolvent-next-gate-scout-2026-07-28
depends_on:
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
related:
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-2690-normal-crossing-cyclic-cubic-resolvent-exclusion-and-reflection-completion-boundary
  - THM-2695-secondary-kummer-bockstein-picard-divisibility-spectrum-and-prime-alignment-boundary
script: 04-computation/jacobian_s4_reflection_relative_different_thm2696.py
output: 05-knowledge/results/jacobian_s4_reflection_relative_different_thm2696.out
script_sha256: 7158a8cc523c5f75b8bfa6a2bbf502dc31dd762df15aa519d0a9a488f84dbfaa
output_sha256: 2e55ad6ea545847f891fbc9cf04c99a90ffca5ea9f4818c9194b6627d2d2393f
hash_basis: LF-normalized bytes
---

# THM-2696 -- the reflection-completed quotient has an intrinsic different

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2655 exhibits the sharp abstract `S4` carrier

```text
V4 semidirect S3 = S4
```

on the toric quotient `d^2=abc`, but its induced quartic quotient map is not
Keller in the displayed elementary coordinates.  The first result here is
that this failure is not a bad choice of quotient coordinates.  The factor
is the intrinsic relative Kahler different: changing polynomial generators
can move its equation, but cannot remove its divisor.

The second result isolates the tempting planar escape.  On the level
`s_1s_2-s_3=c!=0`, the source really is a polynomial `A2` and the relative
Jacobian really is a unit.  Nevertheless its coefficient image is a
singular nonnormal surface, not a polynomial `A2` target.  This closes the
naive constant-different slice without pretending to close arbitrary planar
coefficient surfaces.

## 1. The two quotient coordinate systems

Let

```text
Z=A3_(x,y,z).
```

Let `V=V4` be the four even sign changes, let `H=S3` permute `x,y,z`, and put
`G=V semidirect H`, so `G` is the standard reflection representation of
`S4`.  The two polynomial invariant rings are

```text
S=C[x,y,z]^H=C[s_1,s_2,s_3],
T=C[x,y,z]^G=C[A,B,d],                                  (1)

s_1=x+y+z,
s_2=xy+xz+yz,
s_3=xyz,

A=x^2+y^2+z^2=s_1^2-2s_2,
B=x^2y^2+x^2z^2+y^2z^2=s_2^2-2s_1s_3,
d=xyz=s_3.                                              (2)
```

Thus the degree-four quotient map `F:Spec S -> Spec T` is

```text
F(s_1,s_2,s_3)
 =(s_1^2-2s_2, s_2^2-2s_1s_3, s_3).                    (3)
```

This notation is load-bearing.  The `s_i` are source quotient coordinates;
`A,B,d` are target invariants.  They are not two names for one coordinate
triple.

On the three nonzero characters of `V4`, the order-three subgroup cyclically
permutes the characters and the order-two reflection inverts that cycle.
In matrices over `F2` one may take

```text
C3: [[0,1],[1,1]],        C2: [[0,1],[1,0]].             (4)
```

They satisfy `C3^3=C2^2=1`, `C2 C3 C2=C3^-1`, and generate
`GL_2(F2)=S3`.  This is the rigorous `C3` rotation / `C2` reflection
completion behind `S4/V4=S3`; no free-product identification is being made.

## 2. The missing reflection divisor is the relative different

The special calculation belongs to a general reflection-quotient lemma.
Let `V` be a complex vector space and let

```text
H subset G subset GL(V)                                  (4a)
```

be finite complex reflection groups.  Choose polynomial basic invariants
`h_1,...,h_n` for `H` and `g_1,...,g_n` for `G`.  The inclusion
`C[V]^G subset C[V]^H` gives a finite polynomial map `Phi` satisfying

```text
(g_1,...,g_n)=Phi(h_1,...,h_n).                          (4b)
```

For a reflecting hyperplane `L`, let `e_W(L)` be the order of the pointwise
stabilizer of `L` in `W`, taking `e_W(L)=1` when `L` is not a reflecting
hyperplane of `W`.  The basic-invariant Jacobian formula is

```text
J_W=lambda_W product_L alpha_L^(e_W(L)-1),              (4c)
```

where `alpha_L=0` cuts out `L` and `lambda_W in C^*`.  One quick proof is
local: at the generic point of `L`, the inertia is cyclic of order `e_W(L)`
and the quotient is `t |-> t^e`, so the Jacobian vanishes to order `e-1`.
The total degree `sum_i(deg(g_i)-1)` equals the sum of these local orders,
leaving no additional factor.

Since `H_L subset G_L` are cyclic, `e_H(L)` divides `e_G(L)`.  Applying the
chain rule to `(4b)` gives

```text
Jac(Phi)(h(x))
 =lambda product_L alpha_L^(e_G(L)-e_H(L)).             (4d)
```

The right side is therefore a polynomial (and automatically `H`-invariant).
It is the pullback of the relative Kahler different.  It is constant exactly
when every pointwise stabilizer agrees.  In that case every reflection of
`G` belongs to `H`; since `G` is generated by reflections, `G=H`.
Consequently:

```text
a quotient induced by a proper inclusion of finite complex reflection
groups is never Keller.                                 (4e)
```

Replacing either set of basic invariants by any other polynomial generator
system multiplies `(4d)` only by nonzero constants.  The lemma concerns this
special finite reflection-quotient family, not arbitrary finite polynomial
maps.

### 2.1 The `S3 subset S4` ratio

Direct differentiation of `(3)` gives

```text
Jac(F)=4(s_1s_2-s_3).                                   (5)
```

Set

```text
J_0=s_1s_2-s_3=(x+y)(x+z)(y+z).                        (6)
```

Formula `(5)` is the `S3 subset S4` instance of `(4d)`.  The reflection arrangement of `H`
has anti-invariant

```text
D_H=(x-y)(x-z)(y-z),                                    (7)
```

whereas the full `G` arrangement has anti-invariant

```text
D_G=product_(i<j)(x_i^2-x_j^2)=D_H J_0.                (8)
```

The three minus hyperplanes belong to the source `S3` quotient.  The three
plus hyperplanes are precisely the new relative ramification introduced by
reflection completion.

The degree ledger is the same `2/3/4` mechanism in compressed form.  The
source basic invariants have degrees `(1,2,3)`, while the target basics
`(A,d,B)` have degrees `(2,3,4)`.  Hence the relative Jacobian degree is

```text
(2+3+4)-(1+2+3)=3,                                     (8a)
```

exactly the three added plus hyperplanes and the cubic degree of `J_0`.
Thus the quartic action, its cubic resolvent, and the binary `V4` kernel meet
in the reflection degree ledger rather than through a numerical analogy.

Because `S` is finite over `T` and both are smooth polynomial rings of the
same dimension, the zeroth Fitting ideal of relative differentials is

```text
Fitt_0 Omega_(S/T)=(Jac(F))=(J_0).                      (9)
```

This is the relative Kahler different and its codimension-one ramification
divisor.  Formula `(9)`, unlike a displayed determinant, is independent of
polynomial coordinates.

THM-2695 proves that the standard `Cl[2]` plane on the intermediate toric
quotient is `mu_4`-nonliftable.  That secondary coefficient Bockstein and the
relative different `(9)` are separate obstructions: this theorem neither
constructs nor computes a `mu_4` lift.

## 3. Quartic and cubic-resolvent discriminants agree exactly

The four coset values of the source root are

```text
x+y+z,   x-y-z,   -x+y-z,   -x-y+z.                    (10)
```

Their graph quartic is

```text
P(X)=X^4-2AX^2-8dX+(A^2-4B).                           (11)
```

The three squared pair sums are `4x^2,4y^2,4z^2`, so the corresponding
scaled cubic resolvent is

```text
C(U)=U^3-4AU^2+16BU-64d^2.                             (12)
```

Taking the products of pairwise root differences gives the exact identity

```text
Disc_X(P)=Disc_U(C)
 =4096 product_(i<j)(x_i^2-x_j^2)^2.                   (13)
```

For the normalized cubic

```text
Q(W)=W^3-AW^2+BW-d^2                                  (14)
```

and the source cubic discriminant

```text
Delta_H=Disc_T(T^3-s_1T^2+s_2T-s_3),                   (15)
```

equations `(6)`--`(8)` become

```text
Disc(Q)(F)=Delta_H J_0^2
          =Delta_H (Jac(F)/4)^2,                       (16)

Disc(P)(F)=256 Delta_H Jac(F)^2.                       (17)
```

Thus the classical equality `disc(quartic)=disc(resolvent cubic)` does not
erase the missing divisor.  It displays it twice: target discriminant
pullback equals the source reflection discriminant times the square of the
relative different.

## 4. Polynomial quotient coordinates cannot remove the factor

Let `u_1,u_2,u_3` be any polynomial generator system of `S`, and let
`v_1,v_2,v_3` be any polynomial generator system of `T`.  Each generator
system defines a polynomial automorphism of `A3`.  The Jacobian determinant
of a polynomial automorphism is a nonzero constant because the chain rule
with its polynomial inverse makes it a unit in a polynomial ring.

Consequently, in the new source and target coordinates,

```text
det(partial(v_1,v_2,v_3)/partial(u_1,u_2,u_3))
 =lambda J_0(s(u)),             lambda in C^*.          (18)
```

In particular, every polynomial quotient-coordinate realization of the
fixed inclusion `T subset S` has a nonconstant Jacobian.  The theorem is not
limited to elementary symmetric generators.

There is a slightly larger infinite no-go family.  If `g,h:A3->A3` are
dominant polynomial maps, then characteristic zero makes their Jacobians
nonzero, and

```text
Jac(h o F o g)
 =Jac(h)(F(g)) * 4J_0(g) * Jac(g).                      (19)
```

Dominance makes pullback injective, so `J_0(g)` remains a nonconstant
nonunit.  All factors in `(19)` are nonzero polynomials; their product cannot
be a nonzero constant.  Hence no dominant polynomial sandwich of this fixed
quotient map is Keller.  Such a sandwich need not preserve the original
`S4` cover, so `(19)` is stated as a non-Keller composition family, not as a
new monodromy classification.

The exact companion checks `(18)` on the nontrivial triangular changes

```text
u_1=s_1+s_2^2,       u_2=s_2+s_3,       u_3=s_3,
v_1=A+B^2,           v_2=B+d,           v_3=d,           (20)
```

and checks `(19)` after a dominant non-automorphic quadratic source map.

## 5. The sharp constant-different planar boundary

The finite inclusion `T subset S` is flat by miracle flatness: `S` is
Cohen--Macaulay and `T` is regular of the same dimension.  Localization at
`J_0` makes `(9)` a unit, so `F` is etale on its source etale locus
`D(J_0)`.  This is already a sharp categorical boundary:
`S[J_0^-1]` has the nonconstant unit `J_0`, hence cannot be a polynomial ring
over `C`, whose only units are constants.

A more dangerous escape is an individual nonzero level.  Fix `c in C^*` and
set

```text
J_0=c,
s_1=u,        s_2=v,        s_3=uv-c.                   (21)
```

The source level is exactly `A2_(u,v)`, and `(5)` restricts to the nonzero
constant `4c`.  Its coefficient map is

```text
A=u^2-2v,
B=v^2-2u^2v+2cu,
d=uv-c.                                                 (22)
```

Nevertheless the image is not an affine-plane target.  Define

```text
R_c(A,B,d)
 = 4A^3d^2-A^2B^2-18ABd^2+2ABc^2+4B^3
   +27d^4-18d^2c^2+8dc^3-c^4.                          (23)
```

Then `(22)` lands on `R_c=0`.  The map

```text
A2_(u,v) -> V(R_c) subset A3_(A,B,d)                   (24)
```

is finite and birational.  Finiteness is inherited from the finite invariant
ring inclusion `T subset S`.  For birationality, eliminate `v` from `(22)`:

```text
u^3-Au-2(d+c)=0,
Bu^2+2Adu+(d+c)(3d-c)=0.                                (25)
```

Their resultant is

```text
(d+c)^2 R_c,                                            (26)
```

and the preceding subresultant is linear in `u`, with coefficient

```text
L=4A^2d^2-AB^2+Bc^2-2Bcd-3Bd^2.                        (27)
```

Thus `u`, and then `v=(u^2-A)/2`, is rationally recovered on the dense open
`L!=0`.  Exact Euclidean reduction also gives

```text
gcd(R_c,(d+c)L)=1 in C(c)[A,B,d],                       (27a)
```

so no component is hidden wholly in the exceptional resultant locus.  The
same elimination shows that `(R_c)` is the contraction kernel.
Since `A2` is normal, `(24)` is the normalization of `V(R_c)`.

The target surface is singular for every `c!=0`.  Choose `t` with

```text
t^3=-8c^2.                                               (28)
```

At

```text
(A,B,d)=(-2t,t^2,-c),                                   (29)
```

the three partial derivatives of `R_c` are

```text
4t^2(8c^2+t^3),
4t(8c^2+t^3),
-8c(8c^2+t^3),                                          (30)
```

so they all vanish.  The nonnormality is visible without a criterion: the
two distinct points

```text
(u,v)=(0,t),                  (t^2/(2c),0)               (31)
```

have the same image `(29)`.  Therefore `V(R_c)` is not isomorphic to `A2`.

This is the sharp survivor and failure boundary:

```text
constant relative different on a polynomial source slice
does not imply a polynomial coefficient target.                         (32)
```

The next planar gate is to classify other polynomial surfaces transverse to
the relative different whose coefficient image is smooth, or to prove that
no compatible pair of polynomial source/target parameterizations can make
the quotient data Keller.  The present theorem closes only the obvious
constant-`J_0` slice.

## 6. Exact scope

The theorem proves an intrinsic non-Keller result for the fixed
reflection-completed quotient `C[A,B,d] subset C[s_1,s_2,s_3]`, every choice
of polynomial generators of those two rings, and the dominant polynomial
composition family `(19)`.  It also proves that the most direct planar
constant-different slice normalizes a singular target surface.

It does **not** classify arbitrary `S4` resolvent normalizations, show that
every equality of quartic and cubic discriminants has factorization `(16)`,
compute THM-2695's secondary `mu_4`-lifting class, or rule out a different
smooth planar coefficient surface.  It proves no general degree-four
point-cap theorem, `JC(2)`, or `DC(2)`.

## 7. Exact companion

Run

```text
python 04-computation/jacobian_s4_reflection_relative_different_thm2696.py
python -O 04-computation/jacobian_s4_reflection_relative_different_thm2696.py
```

Both commands reproduce, byte for byte after LF normalization,

```text
05-knowledge/results/jacobian_s4_reflection_relative_different_thm2696.out.
```

The script uses explicit `require` checks.  It verifies the `C2/C3`
standard-plane action, all root and coefficient formulas, both discriminant
normalizations, the relative Jacobian, one nontrivial generator change, one
dominant sandwich, the planar resultant and linear subresultant, the singular
gradient, and the two normalization preimages.
