---
id: THM-4395
title: "Source-normal weight-fourteen row-ten global affine absorption"
status: >
  PROVED FINITE-ROW RELATIVE TO THM-4390 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. In the complete fixed THM-4308 source-normal residual-weight-at-
  most-fourteen family, the row-ten bracket and projected P_2/P_3 depth
  equations form a global affine graph: xi_10, beta_11, c42, and the inherited
  row-nine c14 are solved with constant denominators. The source projection is
  A^6 and every terminal fibre is A^8, including at Phi=0. Row eleven, higher
  weights, chart or seam entry, all-row lifting, Keller pairs, JC(2), and
  DC(2) remain OPEN.
source: root + independent referee / JC2 continuation session, 2026-09-03
depends_on:
  - THM-4390-source-normal-weight-fourteen-row-nine-face-absorption
related:
  - THM-4385-source-normal-row-ten-elliptic-sign-quotient
  - THM-4389-source-normal-weight-thirteen-row-ten-nonisotrivial-elliptic-pencil
  - THM-4390-source-normal-weight-fourteen-row-nine-face-absorption
mistake_firewall:
  - MISTAKE-287
  - MISTAKE-541
primary_script: 04-computation/jc2_source_normal_weight14_row10_global_affine_absorption_thm4395.py
primary_output: 05-knowledge/results/jc2_source_normal_weight14_row10_global_affine_absorption_thm4395.out
primary_script_sha256: 010aae7201abd404dcff4e7c81e4cc1947fbf9e975d6d00d04eb72c2a1e746a0
primary_output_sha256: 9d5686a958ad31aa07b4a2657102aa26ec70a06152a595d7ccea9b72d10e0ecb
independent_audit_script: 04-computation/jc2_source_normal_weight14_row10_global_affine_absorption_thm4395_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_source_normal_weight14_row10_global_affine_absorption_thm4395_independent_audit.out
independent_audit_script_sha256: 44bd447feac75c4d3a5132d57319f865264a6d295b469c3b09089e11ed0e1c36
independent_audit_output_sha256: d27c0492f7211be45d9b36a94be4ea0d213379cfb8bf0f6f312601d599f42da1
hash_basis: raw LF bytes
audit: >
  PASS. The primary imports THM-4390's proved helper implementation and adds
  29 theorem-specific checks. The independent audit imports no THM-4390 or
  THM-4395 code; it reconstructs the literal complete source from the audited
  THM-4308/4315 row operators and executes 80 checks. Both verify all bracket
  coefficients, all projected-depth equations, both boundary and dense
  controls, constant denominators, and the exact fibre dimensions. Normal,
  optimized, and hash-seeded streams byte-match the frozen LF outputs. No
  finite-field inference is used.
---

# THM-4395 -- source-normal weight-fourteen row-ten global affine absorption

**PROVED FINITE-ROW RELATIVE TO THM-4390 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. THIS DETERMINES THE COMPLETE FIXED SOURCE-NORMAL RESIDUAL-WEIGHT-AT-
MOST-FOURTEEN FAMILY ONLY THROUGH ROW TEN. ROW ELEVEN, HIGHER WEIGHTS, CHART
OR SEAM ENTRY, AN ALL-ROW LIFT, POLYNOMIAL TERMINATION, A KELLER PAIR,
`JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. The global row-ten graph

Work over a characteristic-zero field in the complete fixed source universe
of THM-4390:

```text
H_14=H_12+c51*p^5*y+c23*p^2*y^3
           +c70*p^7+c42*p^4*y^2+c14*p*y^4.              (1)
```

THM-4390 proves that rows through eight retain the old Hasse gate and row
nine globally solves `c14`. Continue the exact bracket recursion and both
projected depth modules through row ten. Then the entire source projection is
a global affine graph over the six coordinates

```text
(Phi,eta,alpha_11,c51,c23,c70).                          (2)
```

Equivalently, it is `A^6`. The four solved coordinates are

```text
(beta_11,c42,xi_10,c14),                                 (3)
```

with the triangular formulas in Sections 2--3. Every denominator is a
nonzero rational integer; no `Phi!=0`, `eta!=0`, or other localization is
used. In particular the description includes the complete `Phi=0` fibre.

Every point of this source graph has an affine-eight terminal row-ten fibre.

## 2. The two bracket equations

After the row-nine graph and depth selection, the row-ten bracket selector
has shape `11 x 7`, rank seven, and pivots `(0,1,2,3,4,5,6)`. Its first
independent condition is

```text
13365000*Phi^2+15035625*Phi*eta
+6014250*c42+50787000*c70+57672000*xi_10
-964604821504=0.                                         (4)
```

Thus

```text
xi_10=-(13365000*Phi^2+15035625*Phi*eta
        +6014250*c42+50787000*c70-964604821504)
       /57672000.                                        (5)
```

Substituting `(5)`, the second bracket condition is

```text
B_10=
-104916222000*Phi^2+122625090000*Phi*alpha_11
+19707603750*Phi*beta_11+20802470625*Phi*c51
-246422138625*Phi*eta+20802470625*alpha_11*eta
-89131914000*c42-194981256000*c70
+61312545000*eta^2+2707389207937024=0.                   (6)
```

Both certificates compare all coefficients in degrees `0,...,10`; equations
`(4),(6)` are the complete capped bracket residual, not chosen statistics.

## 3. Projected depth supplies beta, then bracket supplies c42

Before solving `(6)`, the row-ten projected-depth universes are

```text
pi_10(P_2): 88 x 193, rank 68;
pi_10(P_3): 99 x 304, rank 83.                           (7)
```

The joint terminal system has shape `36 x 11`, rank three, and pivots
`(8,9,10)`. It leaves exactly one condition:

```text
-91*Phi+15*beta_11+18*eta=0.                             (8)
```

Consequently

```text
beta_11=91*Phi/15-6*eta/5.                               (9)
```

After `(9)`, equation `(6)` becomes

```text
14643240750*Phi^2+122625090000*Phi*alpha_11
+20802470625*Phi*c51-270071263125*Phi*eta
+20802470625*alpha_11*eta-89131914000*c42
-194981256000*c70+61312545000*eta^2
+2707389207937024=0.                                    (10)
```

Its `c42` coefficient is the nonzero constant `-89131914000`; hence

```text
c42=(14643240750*Phi^2+122625090000*Phi*alpha_11
     +20802470625*Phi*c51-270071263125*Phi*eta
     +20802470625*alpha_11*eta-194981256000*c70
     +61312545000*eta^2+2707389207937024)
    /89131914000.                                       (11)
```

Equations `(11)`, `(5)`, and THM-4390's row-nine formula then determine
`c42`, `xi_10`, and `c14` in that order; `(9)` determines `beta_11`
independently. The certificates substitute the resulting four graphs back
into every literal bracket and depth equation. All vanish. Since the terminal
matrix in `(7)--(8)` has eleven variables and rank three, its remaining fibre
is `A^8`.

The two free weight-thirteen coordinates `c51,c23` play different roles:
`c51` occurs in `(11)`, while `c23` remains absent through row ten. Neither
is constrained.

## 4. Boundary and hostile controls

No coefficient used as a pivot in `(4),(8),(10)` contains a source variable.
At `Phi=0`, the same formulas specialize without degree drop. The independent
audit also chooses a separate rational point on that boundary and evaluates
every original bracket and projected-depth equation, not only the reduced
three-equation presentation.

A dense rational control is obtained from the free coordinates

```text
(Phi,eta,alpha_11,c51,c23,c70)=(1,2,3,5,7,11).          (12)
```

The graph gives

```text
beta_11=11/3,
c42=2705560867462399/89131914000,
xi_10=4294106984995001/316913472000,
c14=218733048242820581/2707381887750.                    (13)
```

All five weight-thirteen/fourteen coefficients are nonzero. At this point,
THM-4385's old elliptic polynomial and THM-4389's weight-thirteen carrier
polynomial have respective values

```text
-814948800613227916,
-814952337033234166.                                    (14)
```

Both are nonzero, while all row-ten equations in the enlarged universe
vanish. Thus the result is not a re-presentation of either capped carrier.

## 5. Mechanism and stopping boundary

Weight fourteen does not create a more complicated replacement for the old
row-ten elliptic geometry. Its even face is visible to THM-4328's stochastic
observer; successive visible coordinates are absorbed by response variables:
`c14` pays row nine, depth fixes `beta_11`, and `c42` pays the final row-ten
bracket condition. The total row-ten carrier is therefore rational and, in
fact, affine.

This mechanism cannot be extrapolated to row eleven. The still-free `c23`
channel first has later nonlinear interactions, and new depth conditions may
cut the affine base `(2)`. Nor is extinction monotone in the weight cap:
weight at least fifteen can enter rows already compiled. Chart/seam entry and
all global Jacobian-conjecture consequences remain open.

## 6. Reproduction

Artifacts:

```text
04-computation/jc2_source_normal_weight14_row10_global_affine_absorption_thm4395.py
05-knowledge/results/jc2_source_normal_weight14_row10_global_affine_absorption_thm4395.out
04-computation/jc2_source_normal_weight14_row10_global_affine_absorption_thm4395_independent_audit.py
05-knowledge/results/jc2_source_normal_weight14_row10_global_affine_absorption_thm4395_independent_audit.out
```

Replay from the repository root:

```text
python3 -B 04-computation/jc2_source_normal_weight14_row10_global_affine_absorption_thm4395.py
python3 -B -O 04-computation/jc2_source_normal_weight14_row10_global_affine_absorption_thm4395.py
python3 -B 04-computation/jc2_source_normal_weight14_row10_global_affine_absorption_thm4395_independent_audit.py
python3 -B -O 04-computation/jc2_source_normal_weight14_row10_global_affine_absorption_thm4395_independent_audit.py
```

The independent audit imports no THM-4390 or THM-4395 code. It rebuilds the
literal source from the audited THM-4308/4315 operators, uses a different
elimination organization, and checks both a `Phi=0` control and the dense
hostile `(12)--(14)`. No finite-field argument is used.
