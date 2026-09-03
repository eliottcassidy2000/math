---
id: THM-4360
title: "Source-normal zeta-zero row-ten delayed-depth extinction"
status: >
  PROVED FINITE-ROW RELATIVE TO THM-4308/4315/4357 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. On the full three-dimensional row-eight source
  plane Phi=0 and xi_10=xi_Z, including the beta_11!=0 stratum and its
  beta_11=0 boundary, the row-nine bracket/depth gate has two reduced
  affine-plane components. The row-ten scalar bracket leaves two affine
  beta_11-lines, but its selected tangent violates an already-required
  row-nine P_2 functional by the fixed nonzero value
  9854451712/1430375. Hence the whole finite plane dies at row ten. The
  beta_11 direction moves the selected jet but is invisible to all final
  obstruction coordinates. No all-row lift, seam entry, Keller pair,
  JC(2), or DC(2) conclusion is asserted.
source: root + zeta-zero row-nine scout + clean-room referee / next-sharp session, 2026-09-02
depends_on:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
  - THM-4357-source-normal-row-eight-endpoint-pullback-stratification
related:
  - THM-4337-zeta-zero-exact-weight-twelve-endpoint-wall-extinction
  - THM-4358-source-normal-s4339-row-ten-delayed-depth-extinction
  - THM-4359-source-normal-row-eight-constructible-response-affine-modification
mistake_firewall:
  - MISTAKE-354
  - MISTAKE-522
  - MISTAKE-539
  - MISTAKE-540
primary_script: 04-computation/jc2_source_normal_zeta_zero_row10_delayed_depth_extinction_thm4360.py
primary_output: 05-knowledge/results/jc2_source_normal_zeta_zero_row10_delayed_depth_extinction_thm4360.out
primary_script_sha256: 9d2553e3bc56b63bd70d53f84a2d1eed867da3757fdbf3023434137bf012ee40
primary_output_sha256: c1727c0a7b77dc4df4720d38cfd106fe96df7e358ae8518c39b6e5631bb89dc5
independent_referee_script: 04-computation/jc2_source_normal_zeta_zero_row10_delayed_depth_extinction_independent_referee_thm4360.py
independent_referee_output: 05-knowledge/results/jc2_source_normal_zeta_zero_row10_delayed_depth_extinction_independent_referee_thm4360.out
independent_referee_script_sha256: d2eb8ba0a5859d351b57f422de33ffb8cef89aec3d0c51db5e9305a219f4a6d6
independent_referee_output_sha256: 6323c3e816ef1bf50bc9087489c077f5e3585b07131f58ede9e7320d525dc4fb
hash_basis: raw LF bytes
audit: >
  PASS. The primary keeps beta_11 symbolic throughout. A separate import-free
  Fraction/sparse-polynomial implementation independently reconstructs the
  source rows, bracket recursion, Student cokernels, projected-depth matrices,
  selected tangents, affine beta shear, and constructible geometry. Both
  202-check implementations pass normally and under -O with byte-identical
  LF output.
---

# THM-4360 -- Source-normal zeta-zero row-ten delayed-depth extinction

**PROVED FINITE-ROW RELATIVE TO THM-4308/4315/4357 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THIS CLOSES ONE THREE-DIMENSIONAL SOURCE-NORMAL PLANE
THROUGH A NECESSARY ROW-TEN TEST. IT IS NOT AN ALL-ROW LIFT, SEAM-ENTRY,
KELLER-PAIR, `JC(2)`, OR `DC(2)` THEOREM.**

## 1. Statement and inheritance

Work over an algebraically closed characteristic-zero field in THM-4308's
source-normal, residual-weight-at-most-twelve finite universe. Put

```text
xi_Z=1563264759296/115860375.
```

The row-eight pullback plane isolated by THM-4357 is

```text
Sbar=V(Phi,xi_10-xi_Z)=A^3_(eta,alpha_11,beta_11),

Delta=896/15,        Theta=512/75,       K=-32/5,
upsilon_5=-731648/2025,                  zeta_3=0,
U=-5200877686784/1042743375,
W= 2539122786304/115860375,              Z=0.             (1)
```

Its constructible decomposition relevant to the endpoint atlas is

```text
Sbar=D(beta_11) disjoint_union V(beta_11),
D(beta_11)=S_4337^o ~= A^2_(eta,alpha_11) x G_m,
V(beta_11)=S_4339 ~= A^2_(eta,alpha_11).                 (2)
```

THM-4358 closed the boundary `S_4339`. The present proof does not infer the
open stratum from that specialization: it retains `beta_11` symbolically and
proves the stronger statement on all of `Sbar`.

The inheritance pass was:

- closest proved mechanism: THM-4358's next-row scalar survivor followed by
  an inherited-depth consumer;
- canonical hostile: the scalar row-ten gate is nonempty even when
  `beta_11` is retained;
- corrected near miss: a parameter absent from the final obstruction need
  not be absent from the selected jet;
- least-used sidecar: the derivative of the selected tangent along the
  apparently silent parameter.

The live concept board was

```text
constructible stratum | closure boundary | moving affine address
Student cokernel | selected old tangent | prior-depth consumer
finite-versus-infinite firewall.                              (3)
```

## 2. Row nine: two planes, and beta really moves the jet

The literal source row is

```text
G_9=(20U+10W+4Z)x^6 +(10alpha_11+6beta_11)x^7
    +(5upsilon_5+4xi_10)x^8 +(eta+zeta_3)x^9.            (4)
```

THM-4315's named normalized gate polynomial restricts to

```text
E9_named=-11F(eta)/102987,
F(eta)=2393096045625 eta^2-415832184456871936.           (5)
```

The literal primitive Student evaluation is instead

```text
143F(eta)/56308142250.                                   (6)
```

Equations `(5)` and `(6)` are nonzero rational multiples and have the same
zero locus; they must not be identified as the same normalization. The
polynomial `F` is irreducible over `Q`, squarefree, and has nonzero roots.
Thus over the algebraic closure the row-nine source gate in `Sbar` has two
reduced `A^2_(alpha_11,beta_11)` components. On the actual open stratum each
becomes `A^1 x G_m`; their `beta_11=0` boundaries are exactly THM-4358's two
`alpha_11`-lines.

The seven-dimensional THM-4308 terminal tangent maps with rank seven under
`G_9`, selecting one old row-eight point over every source point of `F=0`.

Rebuilding the exact projected-depth universes gives

```text
pi_9(P_2): 75 coordinates, 160 columns, rank 59, nullity 16;
pi_9(P_3): 85 coordinates, 251 columns, rank 73, nullity 12.              (7)
```

On the ten new tangent coordinates `q_0,...,q_9`, the combined 28 left-null
residuals have rank three, with pivots `q_7,q_8,q_9`. Setting the other seven
coordinates to zero gives

```text
q_7=-317075357581312/347581125,
q_8=-(125alpha_11+100beta_11+198eta)/15,
q_9=-553237643264/23172075.                              (8)
```

Every residual then vanishes coefficientwise in
`Q[alpha_11,beta_11,eta]/(F)`. Projected depth adds no new source equation,
and the fresh row-nine terminal fibre is `A^7`. The `100beta_11` term in
`(8)` is the first exact warning that beta is a retained address rather than
an absent variable.

## 3. Row ten: scalar survivors and delayed extinction

The literal next source row is

```text
G_10=(15U+10W+6Z)x^8 +(5alpha_11+4beta_11)x^9
     +(upsilon_5+xi_10)x^10.                            (9)
```

The primitive row-ten Student row is

```text
(46189,0,14586,0,15444,0,30888,0,99792,0,489888).       (10)
```

Applied before quotient reduction, it gives

```text
143/84462213375 *
(27281294920125 alpha_11 eta
 +241702700608125 eta^2-25279541057221296128).          (11)
```

Only after `(11)` is formed, reduction modulo `F` yields

```text
143N_E/3128230125,
N_E=1010418330375 alpha_11 eta+619241095293435904.       (12)
```

The direct `4beta_11 x^9` contribution cancels exactly against the
beta-dependent predecessor selected by `G_9`; the derivative of `(11)` with
respect to `beta_11` is zero. Since `eta` is nonzero on `F=0`, the scalar
row-ten gate consists of two affine beta-lines in the closure:

```text
F(eta)=0,
alpha_11 eta=-619241095293435904/1010418330375.          (13)
```

Their intersections with `D(beta_11)` are two copies of `G_m`. Thus a claim
of scalar-cokernel extinction would be false.

The eleven coefficients of `G_10` have tangent rank ten and generic
augmented rank eleven. On `(12)` they select a unique row-nine tangent. When
this selected point is tested against all 28 already-required row-nine
depth functionals, precisely three residual expressions can remain before
imposing `(12)`:

```text
index 12: N_H/14189651847000,
index 14:  13N_E/153248239947600,
index 26: -13N_E/204330986596800,

N_H=5052091651875 alpha_11 eta+3193963923683016704.     (14)
```

Index 12 is the primitive `P_2` annihilator

```text
H_A=35a_(5,0)-20a_(6,2)+10a_(7,4)-4a_(8,6)+a_(9,8),   (15)
```

where `a_(n,r)=[x^r]A_n`. It annihilates all 160 columns of
`pi_9(P_2)`. The exact unit witness

```text
N_H-5N_E=97758447215837184 !=0                          (16)
```

shows that `(12)` and `(15)=0` cannot hold together. Indeed `(12)` forces

```text
H_A=9854451712/1430375 !=0.                             (17)
```

Therefore no scalar survivor retains the already-required row-nine
projected depth. This proves extinction on the whole closure `Sbar`, hence
also on both constructible pieces in `(2)`.

## 4. The invisible moving direction

Although the final obstruction is independent of `beta_11`, the selected
jet is not. Its exact derivative along the beta direction is

```text
d(q_0,...,q_9)/d beta_11 =
(2384/945,0,-504928/8505,0,9484/105,0,-568/7,0,-20/3,0). (18)
```

All 28 selected prior-depth residuals have zero derivative along `(18)`,
even before reduction modulo `F`. Thus an affine shear makes the finite gate
a product with the beta-line, while the actual jet still moves. The quotient
by beta is safe for these named obstruction coordinates only because `(18)`
lies in every consumer kernel; it is not a declaration that beta carries no
information.

This is the same controlled-forgetting pattern as THM-4359's constructible
specialization audit: identify the fibre direction, compute its image under
the next consumer, and quotient only after proving that image is zero.

## 5. Hostiles, geometry, and scope

The exact hostile bank is:

1. Dropping the prior-depth consumer retains the two genuine scalar-survivor
   affine lines, or two `G_m` components on the open stratum.
2. Increasing `beta_11` by one changes `q_8` in `(8)` by `-20/3` and changes
   five coordinates in `(18)`, despite fixing `E9`, `E10`, and `H_A`.
3. Replacing `D(beta_11)` by its closure adds precisely the two beta-zero
   boundary components; extinction is valid on both only because beta was
   retained in the calculation.
4. Squarefreeness of `F` and its nonzero constant term are needed for the
   two-component geometry and for division by `eta` in `(13)`.

The proof order is

```text
row-nine bracket -> selected row-eight point -> row-nine depth fibre
-> row-ten scalar cokernel -> selected row-nine tangent
-> consume already-required row-nine depth.                         (19)
```

Both exact implementations follow this order independently. This theorem
concerns finite projections only. It does not prove infinite `P_d`
membership, extension beyond row ten, polynomial termination, entry of an
arbitrary hypothetical counterexample into this source-normal slice, a
Keller pair, `JC(2)`, or `DC(2)`.

Reproduce from the repository root:

```text
python3 -B 04-computation/jc2_source_normal_zeta_zero_row10_delayed_depth_extinction_thm4360.py
python3 -B -O 04-computation/jc2_source_normal_zeta_zero_row10_delayed_depth_extinction_thm4360.py
python3 -B 04-computation/jc2_source_normal_zeta_zero_row10_delayed_depth_extinction_independent_referee_thm4360.py
python3 -B -O 04-computation/jc2_source_normal_zeta_zero_row10_delayed_depth_extinction_independent_referee_thm4360.py
```
