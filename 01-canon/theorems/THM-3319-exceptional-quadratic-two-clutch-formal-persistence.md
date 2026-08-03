---
id: THM-3319
title: "Exceptional-quadratic two-clutch formal persistence"
status: >
  PROVED by formal algebra + VERIFIED-EXACT tangent audit in the two declared
  accessory fields; INDEPENDENT IMMUTABLE AUDIT PENDING.  Release both affine
  clutch slopes in the THM-3306/3309 family by C=c+d*x and E'=k.  Over each
  degree-36 base field A_i, the transverse linear-subresultant base ideal has
  a unique two-parameter formal continuation through (d,k)=(1,1).  The exact
  d- and k-tangents move both x and c, and their 2-by-2 determinant is nonzero.
  Along this formal surface the surviving quadratic subresultant remains a
  connected finite-etale C2 deck, the two cubics retain gcd exactly two, and
  the true gradient still vanishes on both branches.  Thus the exceptional
  critical deck is not an isolated fixed-slice accident.  This is local
  formal persistence only: no global component, mate, inverse, JC(2), or
  DC(2) consequence follows.
source: root/creative-synthesis-next/2026-08-03
depends_on:
  - THM-3306-affine-c-critical-section-square-discriminant-and-transverse-base-locus
  - THM-3309-exceptional-quadratic-deck-passport-and-gradient-unimodularity-obstruction
related:
  - THM-3289-affine-transverse-c0-e0-coupled-clutch-critical-no-go
  - THM-3318-hamiltonian-divergence-torsion-ladder-for-x-plus-xr-z
script: 04-computation/jc_exceptional_quadratic_two_clutch_formal_persistence_scout_20260803.py
output: 05-knowledge/results/jc_exceptional_quadratic_two_clutch_formal_persistence_scout_20260803.out
script_sha256: 26a6864a31c8834e1b842490127bc2c5467ae7dd98ba9cb775db2b0a1c62bcfd
output_sha256: e9c01067f13dafe679c73c9f150a8540c60589a7611533486771d3c96c3b6d1c
hash_basis: LF-normalized bytes
---

# THM-3319 -- exceptional-quadratic two-clutch formal persistence

**PROVED by formal algebra + VERIFIED-EXACT tangent audit in both accessory
fields; INDEPENDENT IMMUTABLE AUDIT PENDING.**

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

The preceding quadratic subresultant specializes to the separable nonsplit
polynomial of THM-3309.  Its residue algebra

```text
B_i=A_i[t]/(P_2t^2-Q_2t+R_2)                             (7)
```

is a field of relative degree two, with discriminant a nonsquare unit.

## 2. Two-parameter formal continuation

Set

```text
h_d=d-1,                    h_k=k-1,
R_i=A_i[[h_d,h_k]].                                          (8)
```

There are unique formal series

```text
x_i(h_d,h_k), c_i(h_d,h_k) in R_i                          (9)
```

with constant terms `x mod a_0,c_*` such that both coefficients in `(4)`
vanish identically:

```text
a(x_i,c_i,1+h_d,1+h_k)=0,
b(x_i,c_i,1+h_d,1+h_k)=0.                                (10)
```

The raw subresultant row may carry the inherited boundary factor.  That
factor is a unit at `(5)` because `g` is a unit in `A_i`, so it does not
change the zero scheme or the Jacobian rank in the completed local ring.
Equation `(6)` is exactly the formal implicit-function gate for `(10)`.

The first derivatives are determined in `A_i` as follows.  For
`p in {d,k}`, let `a_p,b_p` be the physical parameter derivatives, including
the chain term `x partial_C` for `p=d`.  Then

```text
dot x_p = -a_p/a_x,
dot c_p = -(b_x dot x_p+b_p)/b_c.                         (11)
```

In both accessory fields all four entries in `(11)` are nonzero, and

```text
dot x_d dot c_k-dot x_k dot c_d !=0 in A_i.              (12)
```

Thus the physical clutch plane moves the exceptional base point in two
independent directions; the continued deck is not merely the fixed object
under a redundant parameter change.

## 3. Persistence of the connected critical deck

Along `(10)`, the degree-one subresultant vanishes.  The degree-two row has a
unit leading coefficient at the closed point, hence remains degree two over
`R_i`; the two original cubics retain degree three.  The subresultant PRS over
the fraction field therefore makes that quadratic the last nonzero row, so

```text
deg_y gcd(R_1,R_2)=2                                      (13)
```

throughout the formal section.  The computation independently checks that
the resultant has zero tangent in both clutch directions.

The quadratic discriminant has the THM-3309 nonsquare unit as constant term,
so the lifted quadratic algebra is finite etale over `R_i`.  A finite-etale
idempotent over the complete local ring would reduce to an idempotent of its
special fibre.  Since `(7)` is a field, the lift has no nontrivial idempotent.
Consequently it is a connected rank-two cover whose special fibre is `B_i`.
This is the formal `C_2` deck; no branch is selected over `R_i`.

The exact triangular gradient identities remain

```text
R_1=V P_z,
R_2=V^3P_x-(V'y/2)R_1.                                  (14)
```

Because `V` is a unit on the section, both roots of `(13)` satisfy

```text
P_x=P_z=0.                                               (15)
```

Hence every `Jac(P,Q)` vanishes after pullback to the formal deck.  Gradient
unimodularity still fails before a mate-integrability class can be defined.
The finite-clutch factor `k-A'C` and the owner divisor `g` have unit constant
terms at the residue point, so the formal section does not acquire either
old wall locally.

## 4. Exact tangent audit

The companion reconstructs the generic cubic subresultant sequence before
specialization.  For each accessory field it:

- rebuilds the degree-36 field and the nonsplit discriminant packet;
- differentiates `(4)` by the symbolic chain rule;
- independently recovers the same derivatives by exact central differences
  at parameter values `0` and `2`;
- solves `(11)` and verifies both implicit equations coordinatewise;
- proves `(12)` by a nonzero degree-35 residue element;
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

This result is formal near `(d,k)=(1,1)` and only in the two named accessory
fields.  It does not algebraize the formal section, identify a global
irreducible component, classify the degree-119 residual, or control distant
owner-wall intersections.

THM-3318 supplies the complementary mate obstruction: there the gradient is
unimodular, the divergence class is defined, and a special-fibre torsion
ladder makes it nonzero.  Here the moving critical deck keeps the gradient
nonunimodular, so the pipeline stops one gate earlier.  The shared word
“quadratic” gives no carrier map between those families.

No polynomial mate, inverse, `JC(2)`, or `DC(2)` consequence follows.

QED.
