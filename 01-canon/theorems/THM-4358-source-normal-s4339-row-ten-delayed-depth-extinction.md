---
id: THM-4358
title: "Source-normal S4339 row-ten delayed-depth extinction"
status: >
  PROVED FINITE-ROW RELATIVE TO THM-4308/4315/4357 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. On THM-4357's two-dimensional source plane, the
  row-nine bracket and projected-depth gate is a squarefree quadratic curve:
  two conjugate affine alpha_11-lines over the algebraic closure, with an
  affine seven-dimensional terminal fibre. The row-ten scalar bracket leaves
  two apparent source points, but the full G_10 equation selects a tangent
  that violates an already-required row-nine P_2 functional by the fixed
  nonzero value 9854451712/1430375. Hence no point of this finite row-nine
  gate reaches row ten. No all-row lift, seam entry, Keller pair, JC(2), or
  DC(2) conclusion is asserted.
source: root + S4339 row-nine scout + clean-room referee / next-sharp session, 2026-09-02
depends_on:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
  - THM-4357-source-normal-row-eight-endpoint-pullback-stratification
related:
  - THM-4316-source-normal-row-ten-cubic-corner-extinction
  - THM-2693-odometer-skew-product-three-event-escape-and-uniform-delayed-depth-four-nilpotence
mistake_firewall:
  - MISTAKE-354
  - MISTAKE-522
  - MISTAKE-539
  - MISTAKE-540
primary_script: 04-computation/jc2_source_normal_s4339_row10_delayed_depth_extinction_thm4358.py
primary_output: 05-knowledge/results/jc2_source_normal_s4339_row10_delayed_depth_extinction_thm4358.out
primary_script_sha256: 258edfaeb0bff29a14f91ca34922604965da7820146fc030db3990e7decb1d5c
primary_output_sha256: 8b59909205798b9c0188ce5d9018e274c9fef218944ac4040568544163dba6fe
referee_script: 04-computation/jc2_source_normal_s4339_row10_delayed_depth_extinction_independent_referee_thm4358.py
referee_output: 05-knowledge/results/jc2_source_normal_s4339_row10_delayed_depth_extinction_independent_referee_thm4358.out
referee_script_sha256: d0aa3ed110ef668d993c8a0e81464b9c87ddf927df96fa2b2882b9dafb01ae94
referee_output_sha256: 6d31642ae6ed8753211b583e042775aef7ec25f1c11276a07f758393fa24485c
hash_basis: raw LF bytes
audit: >
  PASS. The 148-check primary reconstructs the canonical row-eight response,
  the general row-nine bracket obstruction, both row-nine depth modules, and
  the row-ten selection. An import-free clean-room implementation independently
  rebuilds the source rows, Student functionals, sparse projection matrices,
  quadratic quotient reduction, primitive P_2 annihilator, and incompatibility
  constant. Normal and optimized runs byte-match both frozen outputs.
---

# THM-4358 -- Source-normal S4339 row-ten delayed-depth extinction

**PROVED FINITE-ROW RELATIVE TO THM-4308/4315/4357 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THIS CLOSES ONE FINITE SOURCE-NORMAL SLICE THROUGH A
NECESSARY ROW-TEN TEST. IT IS NOT AN ALL-ROW LIFT, SEAM-ENTRY, KELLER-PAIR,
`JC(2)`, OR `DC(2)` THEOREM.**

## 1. Statement and inheritance

Work over an algebraically closed characteristic-zero field in THM-4308's
source-normal, residual-weight-at-most-twelve finite universe. THM-4357's
smallest named surviving row-eight source pullback is

```text
S_4339:
Phi=0,     xi_10=1563264759296/115860375,     beta_11=0,
eta,alpha_11 free,

Delta=896/15,       Theta=512/75,       K=-32/5,
upsilon_5=-731648/2025,                 zeta_3=0,
U=-5200877686784/1042743375,
W= 2539122786304/115860375,             Z=0.             (1)
```

The fixed `U,K,W,U+W` are nonzero and the cubic endpoint is squarefree by
THM-4357. The theorem here is:

> The exact row-nine bracket plus projected-`P_2/P_3` locus over `(1)` is
> the squarefree curve `F(eta)=0` in `(4)`. Fibrewise it has one selected
> row-eight terminal point and a new affine seven-dimensional row-nine
> terminal fibre. No point of that curve satisfies the necessary row-ten
> bracket equation while retaining row-nine projected depth. Therefore the
> entire finite slice `S_4339` dies before row ten.

The inheritance pass is:

- closest proved mechanism: THM-4315's general row-nine Student cokernel and
  rank-seven selection of the old terminal fibre;
- canonical hostile: the row-ten scalar condition below is nonempty, so a
  bracket-only extinction claim is false;
- corrected near miss: THM-4315 audited row-nine depth only on its cubic
  corner, so depth membership had to be rebuilt on `S_4339`;
- least-used sidecar: the old row-nine `P_2` left-null functional after the
  next source row has selected a unique tangent.

The live board was

```text
source curve | old-terminal selection | new tangent fibre | next-row cokernel
prior depth functional | specialization quotient | all-row firewall.       (2)
```

## 2. Row nine: a curve, not extinction

The literal source row is

```text
G_9=(20U+10W+4Z)x^6 +(10alpha_11+6beta_11)x^7
    +(5upsilon_5+4xi_10)x^8 +(eta+zeta_3)x^9.            (3)
```

THM-4315's named Hasse-gated polynomial `E9` restricts to

```text
E9|S_4339=-11F(eta)/102987,

F(eta)=2393096045625 eta^2-415832184456871936.           (4)
```

The literal primitive Student evaluation uses a different nonzero rational
normalization; only its zero locus is identified with `(4)`. Since

```text
eta^2=415832184456871936/2393096045625,                  (5)
```

and the reduced numerator and denominator are nonzero nonsquares, `F` is
squarefree and irreducible over `Q`. Thus `(4)` is two disjoint affine
`alpha_11`-lines over the algebraic closure, exchanged by the quadratic
Galois involution, and has no `Q`-rational point. A numeric label `1,2` for
the two lines would require a choice of square root; the invariant discrete
object is the unordered quadratic pair.

The seven-dimensional THM-4308 terminal tangent maps with rank seven under
`G_9`. Hence every source point of `(4)` selects a unique point in the old
row-eight affine fibre.

## 3. Fresh row-nine projected-depth audit

Rebuild the declared finite modules from monomials `x^a u^b p^c y^e` with
`a+b<=d` and source row at most nine. Their exact matrices are

```text
pi_9(P_2): 75 coordinates, 160 columns, rank 59, left nullity 16;
pi_9(P_3): 85 coordinates, 251 columns, rank 73, left nullity 12.           (6)
```

Write the ten new determinant-tangent coordinates as

```text
theta_9=theta9_0+theta9_1 x+...+theta9_9 x^9.            (7)
```

The combined 28 left-null residuals have tangent rank three, with pivots
`theta9_7,theta9_8,theta9_9`. After solving those pivots, every residual is
zero coefficientwise in

```text
Q[alpha_11,eta]/(F).                                    (8)
```

There is therefore no additional source equation at row-nine projected
depth. Fibrewise over `(4)`,

```text
J_8 ~= A^7,       image(J_9 -> J_8) is one point,
J_9 ~= A^7.                                                (9)
```

These are finite projection statements. They do not assert membership in
the infinite modules `P_2,P_3`.

## 4. Row ten: the next row exposes old missing depth

On `(1)`, the literal next source row is

```text
G_10=(2006771806208/13903245)x^8
     +5alpha_11 x^9
     +(1521403518976/115860375)x^10.                    (10)
```

The primitive row-ten Student functional is

```text
(46189,0,14586,0,15444,0,30888,0,99792,0,489888).       (11)
```

It annihilates every `theta_9` contribution. Modulo `F`, its scalar
compatibility condition is

```text
143 E10*/3128230125=0,

E10*=1010418330375 alpha_11 eta+619241095293435904.     (12)
```

This condition is genuinely nonempty. Since `eta!=0`, it fixes one
`alpha_11` on each of the two geometric components:

```text
alpha_11 eta=-619241095293435904/1010418330375.         (13)
```

The full `G_10` coefficient map has rank ten on the ten coordinates `(7)`.
Whenever `(12)` holds, it therefore selects a unique row-nine tangent. Put

```text
a_(n,r)=[x^r]A_n.
```

The primitive functional

```text
H_A=35a_(5,0)-20a_(6,2)+10a_(7,4)-4a_(8,6)+a_(9,8)    (14)
```

annihilates all 160 columns of `pi_9(P_2)`. On the `G_10`-selected tangent,
reduction modulo `F` gives

```text
H_A=N_H/14189651847000,

N_H=5052091651875 alpha_11 eta+3193963923683016704.    (15)
```

But the two necessary numerators obey the exact unit witness

```text
N_H-5E10*=97758447215837184 !=0.                        (16)
```

Thus `(12)` and row-nine `P_2` membership cannot hold together. Equivalently,
after imposing `(12)`,

```text
H_A=9854451712/1430375 !=0.                             (17)
```

This proves the theorem. No row-ten depth module is needed: the selected
tangent already fails the inherited row-nine projection.

## 5. Mechanism and hostile boundary

The mechanism differs from THM-4316's scalar-gcd extinction. A next-row
linear equation may have one scalar cokernel condition and, after that
condition vanishes, select a unique point of the preceding tangent space.
Compatibility with the next row does not imply that the selected point lies
in the prior depth subspace. Here `(14)` is the separating functional and
`(16)` is the exact certificate.

The hostile boundary is `(12)--(13)`: both apparent scalar survivors are
real geometric points over the algebraic closure. Dropping the prior-depth
sidecar would retain them. Conversely, computing `pi_10(P_d)` would add cost
after the contradiction has already occurred. The right operation order is

```text
next scalar cokernel -> selected old tangent -> prior-depth consumer.       (18)
```

This is analogous in shape, not implication, to THM-2693's delayed LRC
carrier: several coarse states can survive before a later consumer reads an
uncompressed sidecar and forces exact extinction.

## 6. Audit and scope

The primary performs 148 exact checks and reconstructs all inherited data
before specializing. The clean-room referee imports no primary code and
independently rebuilds the source and bracket rows, primitive Student
functionals, sparse `pi_8/pi_9` matrices, quotient arithmetic modulo `F`,
rank ledgers, the 160-column annihilation by `(14)`, and `(12)--(17)`. It also
checks the normalization distinction following `(4)`. Normal and optimized
runs byte-match both frozen outputs named in the frontmatter.

The result closes only `S_4339` in the fixed finite source-normal universe
through a necessary row-ten test. It does not prove source-normal or exact-
seam entry for an arbitrary hypothetical counterexample, an all-row `B_2`
lift, polynomial termination, a Keller pair, `JC(2)`, or `DC(2)`. All of
those global questions remain open.
