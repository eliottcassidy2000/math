---
id: THM-3419
title: "Generic Kummer response and regular sector rank"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Let K have characteristic zero, d>=2,
  P=ax+b+g(x)z^d with a!=0, and C_P=K[x,z]/D_P(K[x,z]).  For nonzero g,
  put N=deg(rad(g)).  Every one of the d fiber-exponent sectors of C_P has
  K[P]-rank N, so the total generic response has dimension dN.  After a
  field extension containing the dth roots of unity, the generic de Rham
  response is N copies of the regular mu_d representation; the Hamiltonian
  response is its one-character twist and hence is regular as well.  Root
  multiplicities change compactification genus and puncture counts but cancel
  from every sector rank.  Constant nonzero g gives N=0 and C_P=0 integrally;
  g=0 has the same zero response but is a separate non-Kummer boundary.  This
  computes the response sectors left open in THM-3418 and adds no new Keller
  case or JC(2) conclusion.
source: root-2608-jc-generic-kummer-response-2026-08-15
depends_on:
  - THM-3418-one-monomial-nonlinear-fiber-keller-classification
related:
  - THM-3348-linear-z-generic-puncture-response-and-one-root-valuation
  - THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms
script: 04-computation/jc_generic_kummer_response_thm3419.py
output: 05-knowledge/results/jc_generic_kummer_response_thm3419.out
script_sha256: e60107df6d609d210c61ded6437f9a7e5323f029e271c094695d587989aafcd0
output_sha256: 1b194ef30b97810b1890dbd6d5f3bde8f7ba9d5228bc3f44d03b97cfd5049310
hash_basis: LF-normalized bytes
---

# THM-3419 -- generic Kummer response and regular sector rank

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement and inheritance

Let `K` be a field of characteristic zero, let `d>=2`, and take

```text
R=K[x,z],                   P=ax+b+g(x)z^d,              a in K*,
D=D_P=Jac(P,-),             C_P=R/D(R).                  (1)
```

The exponent of `z` modulo `d` gives the exact decomposition from
[THM-3418](THM-3418-one-monomial-nonlinear-fiber-keller-classification.md):

```text
C_P=direct_sum_(j in Z/dZ) C_j,                          (2)
```

where `C_j` is the target-weight-`j` quotient by the image of input weight
`j+1`.  Multiplication by `P` preserves every `C_j`.

Assume first that `g!=0`, and put

```text
N=deg(rad(g)),                F=K(P).                     (3)
```

Then

```text
rank_(K[P])(C_j)=dim_F(C_j tensor_(K[P]) F)=N             (4)
```

for every `j in Z/dZ`.  Consequently

```text
dim_F(C_P tensor_(K[P]) F)=dN.                            (5)
```

After extending `F` to contain a primitive `d`th root of unity, the natural
`mu_d` action `z |-> zeta z` on the generic de Rham cohomology is

```text
H^1_dR(X/F)=N copies of the regular mu_d representation.  (6)
```

The Hamiltonian cokernel is the one-character twist of `(6)`, hence is also
`N` copies of the regular representation.  Equation `(4)` is defined over
`K`; adjoining roots of unity is needed only to name the individual
characters.

This is the nonlinear-`z` analogue of THM-3348's generic puncture rank, but
the cover now has `d` cyclic sectors.  Its closest hostile is the integral
response in THM-3418: root multiplicity drives an infinite coefficient tower,
whereas the generic rank proved here remembers only the number of distinct
roots.

## 2. The generic fiber and the Kummer cover

Write `t` for the abstract generator of `K[P]~=K[t]`.  Flat localization gives

```text
A=F tensor_(K[P]) R
 ~=F[x,z]/(ax+b+g(x)z^d-t).                               (7)
```

The polynomials `g(x)` and `t-ax-b` are coprime in `F[x]`.  Choose `u,v` with

```text
u g+v(t-ax-b)=1.
```

In `A`, the second factor is `gz^d`, so `g(u+vz^d)=1`.  Thus `g` is a unit and

```text
A~=F[x,g^(-1),z]/(z^d-h(x)),
h(x)=(t-ax-b)/g(x).                                       (8)
```

Put

```text
U=Spec F[x,g^(-1)],                  beta=(t-b)/a.         (9)
```

Then `x:X=Spec A -> U` is a finite degree-`d` Kummer cover.  It is etale over
`U\{beta}` and has one totally ramified point

```text
p=(beta,0)                                                     (10)
```

over `beta`.  The point lies in `U` because `g(beta)!=0`.  The curve is
geometrically connected: over an algebraic closure, `h` has valuation one at
`x=beta`, so `z^d-h` is Eisenstein in the corresponding discrete valuation
ring.  It is smooth as well.  Explicitly,

```text
(1/a)P_x-(g'z/(adg))P_z=1 in A.                           (11)
```

This identity is also the load-bearing type check for the de Rham comparison.

## 3. Hamiltonian cokernel equals de Rham cohomology

Because `D` is `K[P]`-linear, flatness and right exactness give

```text
C_P tensor_(K[P]) F ~= A/D(A).                            (12)
```

Define

```text
Phi: Omega^1_(A/F) -> A,
Phi(dx)=-P_z,                   Phi(dz)=P_x.               (13)
```

The relation `d(P-t)=P_x dx+P_z dz` maps to zero, so `Phi` is well-defined.
The unimodular identity `(11)` makes it an isomorphism.  More explicitly,

```text
eta=(1/a)(dz+(g'z/(dg))dx),              Phi(eta)=1.       (14)
```

For every `q in A`,

```text
Phi(dq)=P_x q_z-P_z q_x=D(q).                            (15)
```

Therefore

```text
A/D(A) ~= Omega^1_(A/F)/dA=H^1_dR(X/F).                  (16)
```

There is one important grading shift.  The form `eta` has `mu_d` weight one,
whereas `Phi(eta)=1` has weight zero.  Equivalently, `Phi` sends de Rham weight
`r` to Hamiltonian target weight `r-1`.  Thus

```text
C_j tensor F ~= H^1_dR(X/F)_(j+1).                        (17)
```

This is why the representation calculation below must be typed before it is
applied to the sector colimits of THM-3418.

## 4. Euler characteristic and total rank

All coefficients in `(7)` descend to a finitely generated subfield `K_0` of
`K`.  Cokernels, the decomposition `(2)`, and algebraic differentials commute
with field extension.  Embed `K_0(t)` in `C`, sending `t` to a number
transcendental over the coefficient field.  Faithful base change therefore
reduces the dimension and character calculation to the corresponding smooth
complex affine curve.

The complex line `U` is the affine line with the `N` distinct roots of `g`
removed.  Removing `beta` as well gives

```text
chi_c(U\{beta})=1-(N+1)=-N.                               (18)
```

Over this set, `(8)` is an unramified cover of degree `d`.  Adding back the
single point `p` yields

```text
chi_c(X)=d(-N)+1=1-dN.                                   (19)
```

The connected smooth affine curve `X` has `H_c^0(X)=0` and
`H_c^2(X)~=C`.  Hence `(19)` gives

```text
dim H_c^1(X)=dN.                                          (20)
```

Poincare duality gives the same dimension for ordinary `H^1`, and algebraic
de Rham comparison transports it to `(16)`.  This proves `(5)`.

## 5. Fixed-point trace and the regular packet

Let `sigma_k` be a nonidentity element of `mu_d`.  A fixed point satisfies

```text
(zeta^k-1)z=0,
```

so `z=0`, and the equation of `X` then forces `x=beta`.  Thus `p` is the only
fixed point.  Since `P_x(p)=a`, `z` is a local parameter at `p`; in that
parameter `sigma_k` is multiplication by `zeta^k`.  Its real fixed-point
index is

```text
sign det_R(1-zeta^k)=+1.                                  (21)
```

The compact-support Lefschetz formula and orientation preservation give

```text
1=L_c(sigma_k)
 =-Tr(sigma_k|H_c^1(X))+Tr(sigma_k|H_c^2(X))
 =-Tr(sigma_k|H_c^1(X))+1.                               (22)
```

Therefore every nonidentity trace on `H_c^1` is zero, while the identity trace
is its dimension `dN`.  This is exactly the character of `N` copies of the
regular representation.  Poincare duality replaces it by the dual
representation on ordinary `H^1`; the regular representation is self-dual.
Equations `(17)` and `(20)--(22)` now prove `(4)` and `(6)`.

## 6. Independent compactification audit

Factor over an algebraic closure as

```text
g=c product_(i=1)^N (x-alpha_i)^(e_i),       r=sum_i e_i. (23)
```

On the base projective line, the valuations of `h` in `(8)` are

```text
1 at beta,             -e_i at alpha_i,       r-1 at infinity. (24)
```

They sum to zero, including the constant case `r=N=0`.  Put

```text
c_i=gcd(d,e_i),             c_infty=gcd(d,|r-1|),
gcd(d,0)=d.                                                (25)
```

Riemann--Hurwitz for the compactified Kummer cover gives

```text
2 genus(Xbar)-2
 =-2d+(d-1)+sum_i(d-c_i)+(d-c_infty).                     (26)
```

The affine curve deletes `sum_i c_i+c_infty` points: all points above the
roots of `g` and infinity.  Consequently

```text
dim H^1(X)
 =2 genus(Xbar)+(sum_i c_i+c_infty)-1
 =dN.                                                      (27)
```

Thus every multiplicity- and inertia-dependent term cancels.  Formula `(27)`
is an independent check of `(19)--(20)` and explains why the generic response
cannot see the repeated-root thickness that governs THM-3418's integral
nontermination.

## 7. Boundaries, hostiles, and JC scope

- If `g=c!=0` is constant, then `N=0` and the theorem gives zero generic
  rank.  More strongly, `K[x,z]=K[P,z]` and `D=a partial_z`, so `C_P=0`
  integrally.
- If `g=0`, the Kummer division in `(8)` is invalid and `rad(g)` is undefined.
  Nevertheless `P=ax+b`, again `D=a partial_z`, and `C_P=0` integrally.
- The hypothesis `a!=0` is sharp for the rank formula.  For
  `P=xz^2` one has `N=1`, but the generic fiber is `G_m`, whose `H^1` has
  dimension one rather than `dN=2`.
- At `d=1` there is one sector and the same cover argument gives rank `N`,
  exactly THM-3348.  The present theorem is stated for `d>=2` because its role
  is to close the nonlinear-fiber sector question from THM-3418.
- Simple, repeated, split, and nonsplit roots all satisfy the response theorem.
  Gradient unimodularity is not assumed.  For a Keller mate, THM-3418 first
  forces repeated roots and then excludes every nonconstant `g` by its
  residue-one recurrence.
- Nonzero `H^1` is not by itself a Keller obstruction: the mate equation asks
  whether the particular class of `1` vanishes in `C_P`.  This theorem computes
  the ambient generic response; THM-3418 already supplies the classification.

Accordingly, `(4)--(6)` are a Hamiltonian-response theorem, not a new case of
the Jacobian Conjecture.

## 8. Exact computational referee

The standard-library companion independently evaluates the valuation/Riemann--
Hurwitz route `(23)--(27)` and the affine-cover Euler route `(18)--(20)` on all
`3,751` packets with `2<=d<=12`, at most four distinct roots, and root
multiplicities from one through four, including every nonzero-constant
boundary.  It checks `26,257` regular-character traces, `26,257` sector ranks,
`108,779` Hamiltonian grading shifts, and `341` `d=1` regressions.  Named
controls include a simple root, a repeated root, two repeated roots, mixed
inertia gcds, and a multiplicity packet compatible with a nonsplit polynomial.
Normal and optimized runs are byte-identical.

Reproduce with

```text
python3 04-computation/jc_generic_kummer_response_thm3419.py
python3 -O 04-computation/jc_generic_kummer_response_thm3419.py
```

Artifact hashes are pinned in the frontmatter.
