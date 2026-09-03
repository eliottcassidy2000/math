---
id: THM-4361
title: "Source-normal beta-zero row-ten joint-depth extinction"
status: >
  PROVED FINITE-ROW RELATIVE TO THM-4308/4315/4357 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. On the beta_11=Z=0, Phi*U*W*zeta_3!=0
  row-eight pullback, row nine leaves a nonempty two-parameter source locus
  with affine-seven terminal fibre. The row-ten bracket plus inherited depth
  cuts this to exactly six geometric source points, naturally three Phi^2
  roots times a sign sheet. Each separate row-ten P_2 or P_3 projection is
  compatible, but their combination is not: after P_2 selects the three
  visible terminal coordinates, a primitive tetrahedral P_3 functional is
  nonzero at all six points by an exact coprimality certificate. No all-row
  lift, seam entry, Keller pair, JC(2), or DC(2) conclusion is asserted.
source: root + beta-zero row-nine scout + clean-room referee / next-sharp session, 2026-09-02
depends_on:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
  - THM-4357-source-normal-row-eight-endpoint-pullback-stratification
related:
  - THM-4334-beta-zero-exact-weight-twelve-endpoint-wall-extinction
  - THM-4358-source-normal-s4339-row-ten-delayed-depth-extinction
  - THM-4359-source-normal-row-eight-constructible-response-affine-modification
  - THM-4360-source-normal-zeta-zero-row-ten-delayed-depth-extinction
  - THM-4362-tetrahedral-diagonal-annihilator-and-sharp-depth-threshold
mistake_firewall:
  - MISTAKE-354
  - MISTAKE-522
  - MISTAKE-539
  - MISTAKE-540
primary_script: 04-computation/jc2_source_normal_beta_zero_row10_joint_depth_extinction_thm4361.py
primary_output: 05-knowledge/results/jc2_source_normal_beta_zero_row10_joint_depth_extinction_thm4361.out
primary_script_sha256: e39e01dc1be019b95886b0965738d696c329042d7de3cbcedc4872955930c602
primary_output_sha256: 5c5f9474b545f6c9ff401912ef46d83ad4644c0319b0acf16b0f2f9c1e556925
independent_referee_script: 04-computation/jc2_source_normal_beta_zero_row10_joint_depth_extinction_independent_referee_thm4361.py
independent_referee_output: 05-knowledge/results/jc2_source_normal_beta_zero_row10_joint_depth_extinction_independent_referee_thm4361.out
independent_referee_script_sha256: 71d6c94bd419456a79420fa12ae38334231dc2a2b9403ddbd12584920fd4f569
independent_referee_output_sha256: 34ebf09468b0a833c762d90ea043807a9ea0d2f78f99d27cab069ecf31588f40
hash_basis: raw LF bytes
audit: >
  PASS WITH WORDING REPAIR. The primary and an import-free Fraction/sparse-
  polynomial reconstruction independently build both rows and all four depth
  matrices. The referee verifies that P_2 and P_3 are separately compatible
  and only jointly inconsistent; the P_3 functional becomes obstructive
  after P_2 terminal selection. Normal/-O/frozen LF streams match.
---

# THM-4361 -- Source-normal beta-zero row-ten joint-depth extinction

**PROVED FINITE-ROW RELATIVE TO THM-4308/4315/4357 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THIS CLOSES ONE FINITE BETA-ZERO SOURCE PULLBACK BY A
JOINT ROW-TEN DEPTH TEST. `P_3` ALONE IS COMPATIBLE. THIS IS NOT AN ALL-ROW
LIFT, SEAM-ENTRY, KELLER-PAIR, `JC(2)`, OR `DC(2)` THEOREM.**

## 1. Source graph and inheritance

Work over an algebraically closed characteristic-zero field in THM-4308's
source-normal, residual-weight-at-most-twelve finite universe. Write

```text
P=Phi,                 X=xi_10,                 Y=P^2.
```

THM-4357's pullback of the THM-4334 endpoint gate is

```text
Z=beta_11=0,                    P*U*W*zeta_3!=0,
zeta_3=-3P/2,

eta=(12506118074368-173745000P^2-926883000X)/(195463125P),
U=(475515904-109350X)/200475,
W=-(4343625P^2-17172000X+143826305024)/4009500.          (1)
```

In particular `P!=0`. The sign involution

```text
(P,eta,alpha_11,zeta_3) -> -(P,eta,alpha_11,zeta_3),
(X,U,W,Z)                ->  (X,U,W,Z)                  (2)
```

will organize the finite survivor set.

The inheritance pass was:

- closest proved mechanism: THM-4315's row-nine Student cokernel followed by
  exact projected-depth selection;
- canonical hostile: THM-4360's smaller zeta-zero plane dies through a prior-
  depth consumer, whereas the present open graph avoids the exceptional
  fibre of THM-4359's affine modification;
- corrected near miss: a `P_3` left-null functional does not make `P_3`
  inconsistent until the competing `P_2` equations select its coordinates;
- least-used sidecar: the two-sheet sign address over the invariant `Y=P^2`.

The live board was

```text
open source graph | row-nine affine fibre | row-ten bracket residual
P_2 terminal selection | P_3 diagonal functional | sign sheet
finite-versus-infinite firewall.                                     (3)
```

## 2. Row nine survives with an affine-seven fibre

The row-nine Student equation is linear in `alpha=alpha_11`. In the invariant
package `A=P^3 alpha`, it is exactly

```text
649499164991015625 A =
 4085005423617421875 Y^2
+24824465812575000000 YX
-278518552828793671680000 Y
-7302452813356500000 X^2
+197059040065115394048000 X
-1329425408965288765218095104.                         (4)
```

Thus row nine fixes `alpha` and imposes no further source equation. The
literal `G_9` map has rank seven on THM-4308's old affine-seven terminal
space and selects one old terminal point.

The exact projected-depth matrices are

```text
pi_9(P_2): 75 rows, 160 columns, rank 59, left nullity 16;
pi_9(P_3): 85 rows, 251 columns, rank 73, left nullity 12.              (5)
```

Their 28 null residuals have rank three on the ten new row-nine tangent
coordinates, with pivots `q_7,q_8,q_9`. They are compatible after `(1),(4)`;
the other seven coordinates remain free. Hence the row-nine gate has a
two-dimensional source locus `(P,X)` and an affine-seven terminal fibre.

For example,

```text
P=1, X=0,
eta=12505944329368/195463125,
alpha=-1329703923433112135272353229/649499164991015625   (6)
```

obeys the strict row-eight graph and the row-nine bracket/depth gate, but
fails the next `X` equation. It is a rational hostile to premature row-nine
extinction.

## 3. The row-ten bracket leaves six source points

On `(1)` the literal source row is

```text
G_10=(15U+10W)x^8+5alpha x^9+(upsilon_5+X)x^10.          (7)
```

Its map on the seven free row-nine tangent directions has rank seven. After
solving those directions, the degree-eight residual is

```text
-4(10125X+1928704)/61875,                               (8)
```

and, after setting `X=-1928704/10125`, the remaining residual is

```text
2q(P^2)/(94585080322265625P^4),                         (9)
```

where

```text
q(Y)=2779225183740234375Y^3
     -194721282033880320000000Y^2
     -1868800030080493839974400000Y
     -9659395340042262184105231777792.                  (10)
```

All other residuals vanish. Thus the exact row-ten bracket locus over the
row-nine depth fibre is

```text
X=-1928704/10125,                    q(P^2)=0.            (11)
```

The cubic is squarefree and has nonzero constant term. At `(11)`,

```text
U=225611776/91125,
W=-13(3375Y+114294784)/40500.                            (12)
```

The following identities modulo `31` prove respectively that `q` is
squarefree and that `q` has no common root with the `W` factor:

```text
(13Y+3)q+(6Y^2-11Y)q'=1,
10q+(12Y^2-12Y-15)(-4Y-15)=1.                           (13)
```

The reduced polynomials keep their displayed degrees, so `(13)` lifts the
corresponding coprimality conclusions to characteristic zero. Consequently
all strict conditions in `(1)` hold at every root. The three distinct,
nonzero `Y` roots each have two square roots `P`, giving exactly six
geometric source points.

Their remaining coordinates are

```text
eta=-8(16875P^2-1231806464)/(151875P),

alpha=(2466241171875P^4-171004979896320000P^2
       -825436857625187713024)/(392122265625P^3).        (14)
```

Equation `(2)` exchanges the two points over each root of `q`. The invariant
address is therefore a `Y` root plus a sign sheet; numbering the six points
`1,...,6` would add a noncanonical ordering and discard this structure.

## 4. Joint projected depth, not standalone P3

The row-ten depth matrices are

```text
pi_10(P_2): 88 rows, 193 columns, rank 68, left nullity 20;
pi_10(P_3): 99 rows, 304 columns, rank 83, left nullity 16.             (15)
```

On the eleven fresh row-ten terminal tangent coordinates, the `P_2`
residuals have rank three and uniquely fix the visible coordinates
`r_8,r_9,r_10`; all twenty `P_2` null equations then vanish. The primitive
functional

```text
H_C=56c_(5,0)-35c_(6,2)+20c_(7,4)
    -10c_(8,6)+4c_(9,8)-c_(10,10)                      (16)
```

lies in the left kernel of the full `pi_10(P_3)` matrix. Here
`c_(n,r)=[x^r]C_n`. After the required `P_2` terminal selection, exact
reconstruction gives the identity

```text
H_C=-P d(P^2)/132327846238037905244160
    +q(P^2)/(248114711696321072332800000P),              (17)

d(Y)=1482253431328125Y^2
     -103851350418069504000Y
     -1468662667625265243357184.                         (18)
```

In particular, the simplified value

```text
H_C=-P d(P^2)/132327846238037905244160                   (19)
```

is asserted only on the bracket locus `q(P^2)=0`. A second identity modulo
`31` is

```text
(8-12Y)q+(-6Y^2+4Y+15)d=1.                              (20)
```

Thus `gcd(q,d)=1`. Since `P!=0`, `(19)` is nonzero at every one of the six
points. No point satisfies both row-ten projected-depth conditions.

This order is essential. Before `P_2` selects `r_8,r_9,r_10`, the `P_3`
terminal matrix has rank equal to augmented rank, namely three, and is
compatible. The same is true for `P_2` alone. The combined terminal matrix
has rank three and augmented rank four. The precise mechanism is therefore

```text
P_2 selects shared coordinates -> P_3 annihilator reads them -> conflict, (21)
```

not standalone `P_3` extinction.

## 5. Why the coefficients are tetrahedral

The absolute coefficients in `(16)` are

```text
(56,35,20,10,4,1)
=(C(8,3),C(7,3),C(6,3),C(5,3),C(4,3),C(3,3)).           (22)
```

They are consecutive tetrahedral numbers with alternating signs along the
diagonal

```text
(n,r)=(5,0),(6,2),(7,4),(8,6),(9,8),(10,10).            (23)
```

THM-4362 proves that this is a universal finite-difference annihilator of
`pi_10(P_3)`, not post-hoc numerology. The role of the present calculation is
different: it shows that `P_2` selection places the candidate jet outside
that universal `P_3` hyperplane.

This also gives the exact cross-problem information map suggested by the
repo's tournament and LRC work:

```text
source: six labelled source points with terminal coordinates
quotient: three Y-roots
sidecar: sign sheet and P_2-selected terminal packet
consumer: H_C
destroyed without sidecar: the sign and the coordinate read by H_C.        (24)
```

No tournament, LRC, or figurate theorem is imported into the proof.

## 6. Hostiles, audit, and scope

The exact hostile bank includes:

1. the rational point `(6)`, proving row nine is genuinely nonempty;
2. either row-ten depth module by itself, proving the joint qualifier is
   necessary;
3. the off-locus term in `(17)`, preventing substitution of `(19)` before
   imposing `q=0`;
4. the characteristic-29 control

   ```text
   (P,X,eta,alpha,beta)=(1,7,8,8,0),
   (zeta_3,U,W)=(13,10,24),       q(1)=0,       d(1)=23,
   ```

   which satisfies the finite source/bracket equations and strict gates but
   retains the depth obstruction. It is an arithmetic control, not a
   characteristic-zero point of the theorem.

The primary and import-free referee reconstruct the source rows, Student
cokernels, all four depth matrices, ranks and augmented ranks, cubic and
quadratic coprimality certificates, sign involution, and finite-field hostile.
The theorem is confined to exact finite projections in the declared source-
normal universe. It proves no infinite `P_d` membership, extension past row
ten, polynomial termination, general seam entry, Keller-pair statement,
`JC(2)`, or `DC(2)`.

Reproduce from the repository root:

```text
python3 -B 04-computation/jc2_source_normal_beta_zero_row10_joint_depth_extinction_thm4361.py
python3 -B -O 04-computation/jc2_source_normal_beta_zero_row10_joint_depth_extinction_thm4361.py
python3 -B 04-computation/jc2_source_normal_beta_zero_row10_joint_depth_extinction_independent_referee_thm4361.py
python3 -B -O 04-computation/jc2_source_normal_beta_zero_row10_joint_depth_extinction_independent_referee_thm4361.py
```
