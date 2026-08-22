---
id: THM-3680
title: "Russell-cylinder Q_dagger coupled stable lift"
status: >
  PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.  The degree-eight
  Q_dagger fold has an actual rational target-ring lift with J_0=1 and
  J_1=J_2=0.  A mod-1000003 rank calculation at restriction cutoff 369 selects
  a 947-by-947 pivot square; that square is rebuilt and solved over Q, and the
  resulting four actual target vectors make the full unselected rational J_1
  and J_2 polynomial residuals vanish.  The nearby cutoff-366 failure is only
  a modular control, not a characteristic-zero minimality theorem.  No J_3,
  J_4, Keller-pair, or JC(2) conclusion is claimed.
source: kps-s194 / Q_dagger coupled actual-lift continuation, 2026-08-21
depends_on:
  - THM-3677-russell-cylinder-degree-eight-fourth-debt-parabola
  - THM-3678-russell-cylinder-qdagger-actual-j0-lift
related:
  - THM-3642-russell-cylinder-zero-debt-actual-lift-and-fourth-stable-closure
script: 04-computation/jc2_russell_cylinder_qdagger_coupled_stable_lift_thm3680.py
output: 05-knowledge/results/jc2_russell_cylinder_qdagger_coupled_stable_lift_thm3680.out
script_sha256: 5a4d32258bdf42aae6c6dc2069eca909fa16a4f199daac6069e3d5fafb07bd34
output_sha256: fdc710ef599f3517971bc57f288858e35bec6fcdb95dcc1cc844592632f9ec2c
hash_basis: raw LF bytes
---

# THM-3680 -- Q_dagger coupled stable lift

**PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.**  The sharp
zero-second-and-fourth-debt candidate does not stop at formal retained data: it
has an actual target-ring lift through the first two stable equations.

All rings are over `C`; the certificate itself is rational.

## 1. Frozen predecessor and equations

Use the Russell cylinder, compiler, quadratic fold, and restriction maps from
THM-3678:

```text
R=C[b,c,e]/(c^2e-b(b+4)),

D=1+x^2q,
b=(D-1)(D+2)^2,       c=xD(D+2),       e=q(D+3),       (1)

(x,t) |-> (x,Q_dagger(x)+t^2,w=t).                    (2)
```

Write `gamma` for restriction to `q=Q_dagger(x)` and `delta` for one normal
`q` derivative followed by restriction.  THM-3678 freezes actual elements
`F_1,G_1 in R` such that, for `U=c,V=e`,

```text
J_0=c'G_1-F_1e'=1.                                    (3)
```

For readability, every target element in the formulas below denotes its
`gamma` restriction, and primes denote `x` derivatives.  Introduce actual
target elements `F_2,G_2,F_3,G_3`.  Direct expansion of the pulled-back source
Jacobian gives

```text
J_1=K_1-2F_2e'+2c'G_2,                                (4)

J_2=K_2+F_2'G_1-2F_2G_1'+2F_1'G_2-F_1G_2'
        -3F_3e'+3c'G_3,                               (5)
```

where the frozen known terms are

```text
K_1=2c'delta(e)+F_1'G_1-F_1G_1'-2delta(c)e',          (6)

K_2=3c'delta(G_1)+2F_1'delta(e)+delta(c)'G_1
       -F_1delta(e)'-2delta(c)G_1'-3delta(F_1)e'.      (7)
```

The companion independently expands and gates these full polynomial formulas;
the calculation is not performed only at the three retained points.

## 2. Modular pivot selection

The restriction-degree weights are `(30,21,18)`.  Let

```text
M_N={b^i c^j e^k:30i+21j+18k<=N}                      (8)
```

in nested `(i,j,k)` order, and use one copy of `M_N` for each of
`F_2,G_2,F_3,G_3`.  Reducing the coupled coefficient matrix `(4)`--`(5)`
modulo `1000003` gives

```text
 N    |M_N|  rows J_1+J_2  columns  operator rank  augmented rank
366     953      387+561      3812        941             942
369     973      390+564      3892        947             947. (9)
```

The first row is a hostile finite-field control.  It does **not** prove that a
rational cutoff-366 solution is impossible: reduction modulo one prime need
not preserve rational solvability when solution denominators contain that
prime.  Accordingly, no characteristic-zero minimality claim is attached to
`369`.

At cutoff `369`, modular RREF selects `947` operator columns.  Transposed RREF
then selects `947` independent rows.  The deterministic pivot hashes are

```text
columns 0aae9293...635c37,
rows    6f65529d...49315f.                             (10)
```

The companion gates the full hashes and every pivot.

## 3. Characteristic-zero lift

The selected `947 x 947` matrix and right-hand side are rebuilt from the
original rational polynomial coefficients, not lifted from their modular
values.  Exact FLINT elimination solves this square over `Q`.  Substitution of
the selected solution into **all** `390+564` coefficient rows then gives

```text
J_1=0,                         J_2=0 in Q[x].           (11)
```

Thus `(3)` and `(11)` are full rational polynomial identities.  The actual
target coefficient vectors and their restrictions are frozen by

```text
       target support  restriction degree  target hash       restriction hash
F_2         281               369           ab540e28...0d144e  f0cddc4c...8689ff
G_2         278               366           340a5c25...d7d7c66 5c46e799...df1536
F_3         281               369           e51511c6...2eb750  b4eb6e64...29b731
G_3         104               369           b42158e3...37479ae bf1ebe51...7164f3. (12)
```

The full 64-hex-digit hashes are pinned in the companion and stored transcript.
Because each vector is indexed by the raw target monomials `(8)`, it defines an
actual element of `R`; the proof is not merely a source-polynomial Bezout
calculation.

## 4. Why this is the live frontier

THM-3642 exhibited an actual `J_0=1,J_1=J_2=0` lift for `Q_6`, but its universal
fourth-stable identity immediately forced a nonzero fourth debt.  For
`Q_dagger`, THM-3677 proves instead that the retained second and fourth debts
both vanish.  Combining that fact with this theorem gives the first audited-
chain candidate in this family known simultaneously to

```text
have an actual target-ring lift through J_2,
and pass the universal retained J_4 debt test.          (13)
```

This is a survival result, not a counterexample.  The next equations require
actual `F_4,G_4,F_5,G_5` and the higher normal jets of the frozen vectors.
Neither `J_3=0` nor `J_4=0` is proved here, and no infinite algebraization or
Keller-pair conclusion follows.

## 5. Reproduction

From the repository root:

```bash
python -B 04-computation/jc2_russell_cylinder_qdagger_coupled_stable_lift_thm3680.py
python -O -B 04-computation/jc2_russell_cylinder_qdagger_coupled_stable_lift_thm3680.py
```

The companion pins and replays the exact THM-3678 predecessor before rebuilding
the modular skeleton, rational square, and full residuals.  Both modes must
return the stored transcript.
