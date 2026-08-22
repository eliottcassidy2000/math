---
id: THM-3319
title: "Exceptional-quadratic two-clutch algebraic etale persistence"
status: >
  PROVED by algebraic geometry + VERIFIED-EXACT tangent audit in the two
  declared accessory fields; INDEPENDENT IMMUTABLE AUDIT PENDING.  Release
  both affine clutch slopes in the THM-3306/3309 family by C=c+d*x and E'=k.
  The two linear-subresultant coefficients cut out an algebraic surface whose
  projection to the (d,k)-plane is etale at each degree-36 exceptional point.
  Its completion is the unique two-parameter formal continuation through
  (d,k)=(1,1).  The exact d- and k-tangents move both x and c, and their
  2-by-2 determinant is nonzero.  On a neighborhood in this etale germ the
  surviving quadratic subresultant remains a connected finite-etale C2 deck,
  the two cubics retain gcd exactly two, and the true gradient still vanishes.
  Thus the exceptional critical deck is not an isolated fixed-slice accident.
  This is a local algebraic etale germ only: no global component, rational
  section, mate, inverse, JC(2), or DC(2) consequence follows.
source: root/creative-synthesis-next/2026-08-03
depends_on:
  - THM-3306-affine-c-critical-section-square-discriminant-and-transverse-base-locus
  - THM-3309-exceptional-quadratic-deck-passport-and-gradient-unimodularity-obstruction
related:
  - THM-3289-affine-transverse-c0-e0-coupled-clutch-critical-no-go
  - THM-3312-exceptional-quadratic-trace-norm-and-cofactor-antidescent
  - THM-3318-hamiltonian-divergence-torsion-ladder-for-x-plus-xr-z
script: 04-computation/jc_exceptional_quadratic_two_clutch_formal_persistence_scout_20260803.py
output: 05-knowledge/results/jc_exceptional_quadratic_two_clutch_formal_persistence_scout_20260803.out
script_sha256: 4a578ecf29351027762d3fd1e48dadc0d5ba162478e97c91eafbeed55c774f38
output_sha256: 166b12fa4e81b8b0c4ba46fc532c220e4f3088ad89e2b2b67163034f61c7da28
hash_basis: LF-normalized bytes
---

# THM-3319 -- exceptional-quadratic two-clutch algebraic etale persistence

**PROVED by algebraic geometry + VERIFIED-EXACT tangent audit in both
accessory fields; INDEPENDENT IMMUTABLE AUDIT PENDING.**

## 1. Family and fixed residue deck

Retain either THM-3212 accessory response field `K_i`, its response pair

```text
V=4SDT^2/Gamma^2,       A=2SET/Gamma,
g=ST=rad(V),            2VA'-AV'=2V,                    (1)
```

and the simultaneous affine clutch family

```text
C=c+d x,                E=e_0+kx,
P=(Vz^2+z+C)^2+Az+E.                                    (2)
```

The constant `e_0` is gradient-inert.  On `V!=0`, put `y=Vz` and

```text
L=y^2+y+CV,
R_1=2L(2y+1)+VA,
R_2=V^3k+V^2y+L(-V'y+2V^2d).                            (3)
```

Let the degree-one member of the generic subresultant sequence be

```text
S_1=a(x,c,d,k)y+b(x,c,d,k),                              (4)
```

using the inherited boundary localization.  At `(d,k)=(1,1)`, THM-3306
gives

```text
a=a_0(x),                 b=b_0(x)+c b_1(x),
A_i=K_i[x]/(a_0),         c_*= -b_0/b_1 in A_i,          (5)
```

where `A_i/K_i` has degree `36` and

```text
a_0'(x)b_1(x) in A_i^*.                                  (6)
```

The preceding quadratic subresultant can be written on a neighbourhood of
`xi_i` as

```text
S_2=P_2y^2+Q_2y+R_2^(0).
```

Here `R_2^(0)` denotes the constant coefficient of the quadratic row; `R_2`
remains reserved for the localized gradient cubic in `(3)`.  At
`(d,k)=(1,1)`, `S_2` specializes to the separable nonsplit polynomial of
THM-3309.  Using the same symbols for its residue coefficients and putting
`t=-y`, its residue algebra is

```text
B_i=A_i[t]/(P_2t^2-Q_2t+R_2^(0)).                        (7)
```

It is a field of relative degree two, with discriminant a nonsquare unit.

## 2. Algebraic etale germ and its completion

Let `U_i` be the inherited boundary localization of
`Spec K_i[x,c,d,k]` on which the normalized coefficients `a,b` and the
subresultant rows are regular, and define

```text
Z_i=V(a,b) subset U_i,             pi_i:Z_i -> A^2_(d,k). (8)
```

Let `xi_i` be the closed point over `(d,k)=(1,1)` defined by `(5)`; its
residue field is `A_i`.  At `xi_i`,

```text
det partial_(x,c)(a,b)=a_0'(x)b_1(x) in A_i^*.           (9)
```

The relative Jacobian criterion therefore makes `pi_i` etale at `xi_i`.
After shrinking around that point, `Z_i` is an algebraic smooth surface and
its projection to the physical clutch plane is etale.  This is an algebraic
germ with a degree-36 closed point, not a rational section over `K_i` and not
a classification of a global irreducible component.

Set

```text
h_d=d-1,                    h_k=k-1,
R_i=A_i[[h_d,h_k]].                                        (10)
```

Completing the etale local germ at `xi_i` gives `R_i`.  Equivalently, there
are unique formal series

```text
x_i(h_d,h_k), c_i(h_d,h_k) in R_i                         (11)
```

with constant terms `x mod a_0,c_*` such that both coefficients in `(4)`
vanish identically:

```text
a(x_i,c_i,1+h_d,1+h_k)=0,
b(x_i,c_i,1+h_d,1+h_k)=0.                                (12)
```

The raw subresultant row may carry the inherited boundary factor.  That
factor is a unit at `(5)` because `g` is a unit in `A_i`, so it does not
change the local zero scheme or the Jacobian rank.  Equation `(6)` is both
the algebraic etale gate for `(8)` and the formal implicit-function gate for
`(12)`.

The first derivatives are determined in `A_i` as follows.  For
`p in {d,k}`, let `a_p,b_p` be the physical parameter derivatives, including
the chain term `x partial_C` for `p=d`.  Then

```text
dot x_p = -a_p/a_x,
dot c_p = -(b_x dot x_p+b_p)/b_c.                        (13)
```

In both accessory fields all four entries in `(13)` are nonzero, and

```text
dot x_d dot c_k-dot x_k dot c_d !=0 in A_i.              (14)
```

Thus the physical clutch plane moves the exceptional base point in two
independent directions; the continued deck is not merely the fixed object
under a redundant parameter change.

## 3. Persistence of the connected critical deck

On `Z_i`, the degree-one subresultant `S_1=ay+b` vanishes.  After shrinking
around `xi_i`, `P_2` and the leading coefficients of both original cubics are
units, while `S_2` remains nonzero of degree two.  The inherited subresultant
recurrence, with all row normalizations absorbed into a regular unit, gives

```text
Res_y(S_2,S_1)
 =P_2b^2-Q_2ab+R_2^(0)a^2
 =u Res_y(R_1,R_2),        u in O_(U_i,xi_i)^*.
```

Because `a=b=0` on `Z_i`, the cubic resultant therefore vanishes identically
there, not merely to first order.  Thus `S_2` is the last nonzero PRS row over
the function field of the germ, and

```text
deg_y gcd(R_1,R_2)=2                                      (15)
```

Equation `(12)` describes the completed base germ only.  Over that base,
adjoining `t=-y` subject to `S_2(-t)=0` gives the rank-two formal common-root
algebra.  The exact tangent calculation in Section 4 is a redundant
first-order check of the resultant identity, not the source of its all-order
vanishing.

The quadratic discriminant is a unit at `xi_i`, with the THM-3309 nonsquare
as residue.  Thus the quadratic algebra is finite etale of rank two over the
local ring of `Z_i` at `xi_i`.  Any idempotent reduces to an idempotent of its
special fibre `(7)`; because that fibre is a field, Nakayama's lemma rules out
a nontrivial lift.  Consequently the local algebraic cover is connected and
its completion is the connected formal `C_2` deck over `R_i`.  No branch is
selected over either germ.

The exact triangular gradient identities remain

```text
R_1=V P_z,
R_2=V^3P_x-(V'y/2)R_1.                                  (16)
```

Because `V` is a unit on the germ, both roots of `(15)` satisfy

```text
P_x=P_z=0.                                               (17)
```

Hence every `Jac(P,Q)` vanishes after pullback to the algebraic deck.
Gradient unimodularity still fails before a mate-integrability class can be
defined.  The finite-clutch factor `k-A'C` and the owner divisor `g` are units
at the residue point, so after shrinking the etale germ neither old wall is
acquired locally.

## 4. Exact tangent audit

The companion reconstructs the generic cubic subresultant sequence before
specialization.  For each accessory field it:

- rebuilds the degree-36 field and the nonsplit discriminant packet;
- differentiates `(4)` by the symbolic chain rule;
- independently recovers the same derivatives by exact central differences
  at parameter values `0` and `2`;
- solves `(13)` and verifies both implicit equations coordinatewise;
- proves `(14)` by a nonzero degree-35 residue element;
- checks all three quadratic-row coefficients remain nonzero at the residue;
- checks the resultant tangent in both directions; and
- audits zero assertions and zero floating literals.

The two determinant digests are

```text
4111  3f0a0d1d6385099c25f6fa42aaeb11be990d59dba2260c56f9c0c1395f3b63ec
3211  1aba7e3b235a75045ab90bbeaeb5e58c00ef8f0af2cc23854ead0ca4fb1b65d1
```

Normal and optimized modes reproduce the frozen output byte for byte.

## 5. Boundary and obstruction taxonomy

This result is a local algebraic etale germ near `(d,k)=(1,1)` and only in
the two named accessory fields.  It does not identify or normalize a global
irreducible component, give a `K_i`-rational section of the clutch plane,
classify the degree-119 residual, or control distant owner-wall intersections.

THM-3318 supplies the complementary mate obstruction: there the gradient is
unimodular, the divergence class is defined, and a special-fibre torsion
ladder makes it nonzero.  Here the moving critical deck keeps the gradient
nonunimodular, so the pipeline stops one gate earlier.  The shared word
“quadratic” gives no carrier map between those families.

No polynomial mate, inverse, `JC(2)`, or `DC(2)` consequence follows.

QED.
