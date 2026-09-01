---
id: THM-4316
title: "Source-normal row-ten cubic-corner extinction"
status: >
  PROVED RELATIVE TO THM-4312 AND THM-4315 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED.  Every source-normal residual-weight-at-most-twelve
  cubic-corner k=1 point surviving the exact row-nine bracket/depth gate
  fails the necessary row-ten bracket compatibility.  The row-nine locus is
  the quintic Q(Phi^2)=0; the row-ten Student--Stein obstruction is a cubic
  P10(Phi^2), and an explicit modular Bezout identity proves gcd(Q,P10)=1.
  No row-ten depth theorem is needed.  This closes the cubic-corner k=1
  lane in the fixed normalized finite-row universe, not the U=0 or Z=0
  walls, seam entry, JC(2), or DC(2).
source: root / planar-Jacobian stochastic-process bridge session, 2026-09-01
audit: >
  PASS.  The primary SymPy path expands literal G10, reconstructs the unique
  G9 transport, verifies the full row-nine tangent variation, derives the
  row-ten Student cokernel and obstruction, rebuilds Q, and checks exact and
  mod-7 coprimality.  A standard-library Fraction/sparse-Laurent audit
  reconstructs source rows G4--G10 and the bracket transport without
  importing either primary certificate or a CAS; it independently verifies
  the Laurent obstruction identity and exact/mod-23 Bezout certificates.
  Normal and optimized runs byte-match both frozen transcripts.
depends_on:
  - THM-4312-source-normal-cubic-corner-repeated-face-collapse
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
related:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
  - THM-3163-universal-finite-prefix-markov-realization-and-physical-sidecar-boundary
  - THM-4317-source-normal-k1-completed-local-resolution-green-kernel
primary_script: 04-computation/jc2_source_normal_row10_cubic_corner_extinction_thm4316.py
primary_output: 05-knowledge/results/jc2_source_normal_row10_cubic_corner_extinction_thm4316.out
primary_script_sha256: 64e9c368f756e22c79b893cd2c43cd52aaee1ab6bfd9bf11b5ae79913ca05200
primary_output_sha256: b59ea69a67285067d65dd524075c37dd824709a48e900de074b5a51e9eabf4fc
independent_audit_script: 04-computation/jc2_source_normal_student_stein_row10_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_source_normal_student_stein_row10_independent_audit.out
independent_audit_script_sha256: b50eaf6013526473d1e5e8e1b74cb93175bf5850d013f24970ea552594daae62
independent_audit_output_sha256: 801e80ca42cb86fd4543b97dbccd1fd66bb3af745800c1e17ac858d778552b05
hash_basis: raw LF bytes
---

# THM-4316 -- Source-normal row-ten cubic-corner extinction

**PROVED RELATIVE TO THM-4312 AND THM-4315 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.  NO CUBIC-CORNER `k=1` POINT IN THE EXACT ROW-NINE
GATE HAS A DEGREE-CAPPED BRACKET LIFT THROUGH ROW TEN.  THIS IS NOT A
SEAM-ENTRY, `JC(2)`, OR `DC(2)` THEOREM.**

## 1. Statement and inheritance

Retain the fixed characteristic-zero source-normal gauge, residual weight at
most twelve, THM-4308's row-eight bracket/depth gate, THM-4312's exact
cubic-corner `k=1` equations, and THM-4315's row-nine compatibility.  Put
`X=Phi^2`.

THM-4315 proves that the allowed row-nine source points are exactly the ten
points above the five distinct roots of its quintic `Q(X)`.  The theorem
here proves

```text
Q(X)=0             implies no bracket lift through G_10. (1)
```

The inheritance board is deliberately short:

- closest mechanism: THM-4315's Student--Stein row cokernel;
- hostile: its `F_19` row-nine survivor, which shows that small proposal
  probability is not extinction;
- sidecar: the prescribed literal source row `G_10`, without which the
  universal tangent operator selects no continuation;
- decisive test: exact coprimality of the row-nine and row-ten source
  equations.

## 2. Literal row ten and its tangent quotient

Direct source expansion gives

```text
G_10=(15U+10W+6Z)x^8
     +(5alpha_11+4beta_11)x^9
     +(upsilon_5+xi_10)x^10.                            (2)
```

On the cubic corner, `W=-2U`, `Z=U`, and `beta_11=-alpha_11`, so `(2)`
becomes

```text
G_10=Ux^8+alpha_11x^9+c_2x^10.                         (3)
```

THM-4315's row-nine equation selects the unique row-eight point which can
reach `G_9`.  Solve the row-nine determinant equation with a particular row
`v_9^0`; every other solution is

```text
v_9=v_9^0+theta_9 v_0',               deg(theta_9)<=9. (4)
```

Its entire effect on the next bracket prediction is

```text
D_10 theta_9
 =[(x^2+6)theta_9'-20x theta_9]/20.                    (5)
```

This map has rank ten and primitive cokernel

```text
(46189,0,14586,0,15444,0,30888,0,99792,0,489888).      (6)
```

Thus the obstruction below is independent of every row-nine tangent choice,
including the seven-dimensional projected-depth fibre.

## 3. The row-ten obstruction

THM-4315's `E9=0` response fixes `alpha_11` as a rational function of
`Phi`.  Substituting it and the corner response into `(2)`, then applying
`(6)` to `G_10` minus the particular bracket prediction, gives exactly

```text
E_10=
143 P_10(Phi^2)
-----------------------------------------,               (7)
2518243204518514160156250 Phi^4
```

where

```text
P_10(X)=
 68842335386673891964107421875 X^3
-114199708075156991490870528000000 X^2
-49516799750570385919467992383488000 X
-2955996966894649715961849999793382752256.              (8)
```

The cubic is primitive and squarefree.  The denominator in `(7)` is legal:
THM-4315 proves that `Q` is coprime to `X`, so every row-nine point has
`Phi!=0`.

## 4. Exact extinction by coprimality

Let `Q` be the row-nine quintic in THM-4315 `(14)`.  Exact Euclidean
division gives

```text
gcd(Q,P_10)=1.                                           (9)
```

A compact certificate is visible modulo seven:

```text
Q      =3X^5+2X^4+X^3-2X^2-3,
P_10   =3X^3-X^2+1,

(X^2+3X+3)Q
+(-X^4+3X^3-2X^2+2X+3)P_10=1          in F_7[X].        (10)
```

Every row-nine point satisfies `Q(X)=0`, whereas every row-ten bracket lift
would have to satisfy `P_10(X)=0` by `(7)`.  Equations `(9)--(10)` make
those requirements disjoint.  This proves `(1)`.

No row-ten depth projection is required: bracket compatibility is already a
necessary condition for any depth-filtered lift.  This is deterministic
zero offspring at every row-nine source point, not merely small or zero
measure under a chosen proposal.

For the stochastic hostile from THM-4315, `X=7` survives row nine over
`F_19`, but

```text
P_10(7)=15 !=0                         in F_19.           (11)
```

Thus the same control cleanly separates proposal rarity (`19^(-38)`) from
the actual next-row obstruction.

## 5. Consequence and scope

Combining THM-4312, THM-4315, and this theorem gives:

> Inside the fixed normalized, residual-weight-at-most-twelve cubic corner,
> the only row-eight repeated-face lane is `k=1`; its row-nine gate consists
> of ten source points; none reaches row ten.

This closes that cubic-corner lane in the declared finite-row universe.  It
does not address the exact-`M=12` endpoint walls `U=0` or `Z=0`, prove entry
into the source-normal seam from an arbitrary hypothetical counterexample,
establish polynomial termination, or imply `JC(2)` or `DC(2)`.

Reproduce with

```text
python3 -B 04-computation/jc2_source_normal_row10_cubic_corner_extinction_thm4316.py
python3 -B -O 04-computation/jc2_source_normal_row10_cubic_corner_extinction_thm4316.py
python3 -B 04-computation/jc2_source_normal_student_stein_row10_independent_audit.py
python3 -B -O 04-computation/jc2_source_normal_student_stein_row10_independent_audit.py
```

Both normal and optimized streams must byte-match their frozen outputs.

**QED.**
