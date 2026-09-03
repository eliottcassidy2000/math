---
id: THM-4390
title: "Source-normal weight-fourteen row-nine face absorption"
status: >
  PROVED FINITE-ROW RELATIVE TO THM-4308/4315 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. In the complete fixed THM-4308 source-normal
  residual-weight-at-most-fourteen family, rows through eight retain the old
  Hasse gate. The unique row-nine bracket equation is the old E9 equation
  minus one explicit nonzero linear coordinate on the complete weight-
  fourteen face; it globally solves c14, and projected P_2/P_3 depth is
  automatic. Row ten, higher weights, chart or seam entry, all-row lifting,
  Keller pairs, JC(2), and DC(2) remain OPEN.
source: root + independent referee / JC2 continuation session, 2026-09-03
depends_on:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
related:
  - THM-4328-seam-covariant-student-stein-face-visibility
  - THM-4388-source-normal-weight-thirteen-row-twelve-depth-extinction
mistake_firewall:
  - MISTAKE-541
primary_script: 04-computation/jc2_source_normal_weight14_row9_face_absorption_thm4390.py
primary_output: 05-knowledge/results/jc2_source_normal_weight14_row9_face_absorption_thm4390.out
primary_script_sha256: b87113fc8c4c83c19d1076377f3aa14c03b8b64f8d31ae9551281bc224c26526
primary_output_sha256: e820b44572f8966c278b0bd8b4c1939eb517c416e0061aa24cee1300da278bde
independent_audit_script: 04-computation/jc2_source_normal_weight14_row9_face_absorption_thm4390_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_source_normal_weight14_row9_face_absorption_thm4390_independent_audit.out
independent_audit_script_sha256: 80315b1fe71929d4590bb3fe2ca5b395ff2fd8ba02db44d9a4a89879eadd457a
independent_audit_output_sha256: f2a0434157f7ca3fd412f7d62d5320913bac233eb3395a654cdf267b0fd61ec4
hash_basis: raw LF bytes
audit: >
  PASS. The primary performs 66 checks. The independent certificate imports
  no THM-4390 code, reconstructs the two complete new faces, uses a separately
  organized bracket and depth elimination, and performs 51 checks. Both
  import the audited THM-4308/4315 row operators. Normal, optimized, and
  hash-seeded executions byte-match the frozen LF outputs. No finite-field
  inference is used.
---

# THM-4390 -- source-normal weight-fourteen row-nine face absorption

**PROVED FINITE-ROW RELATIVE TO THM-4308/4315 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THIS DETERMINES THE COMPLETE FIXED SOURCE-NORMAL
RESIDUAL-WEIGHT-AT-MOST-FOURTEEN FAMILY ONLY THROUGH ROW NINE. ROW TEN,
HIGHER WEIGHTS, CHART OR SEAM ENTRY, AN ALL-ROW LIFT, POLYNOMIAL TERMINATION,
A KELLER PAIR, `JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. Universe and conclusion

Work over a characteristic-zero field in THM-4308/4315's fixed normalized
source chart. Starting with its complete residual-weight-at-most-twelve
polynomial, adjoin both complete next faces:

```text
H_14=H_12+c51*p^5*y+c23*p^2*y^3
           +c70*p^7+c42*p^4*y^2+c14*p*y^4.              (1)
```

Indeed, the nonnegative solutions of `2a+3b=13` and `2a+3b=14` are

```text
13: (a,b)=(5,1),(2,3),
14: (a,b)=(7,0),(4,2),(1,4).                             (2)
```

After the inherited row-eight Hasse gate, the entire row-nine bracket and
projected-depth condition on this family is

```text
E_9-(3339765000*c70+898128000*c42+216513000*c14)=0,     (3)
```

where

```text
E_9=613527750*Phi^2-511211250*Phi*alpha_11
    -3154140000*Phi*eta-255605625*eta^2
    +6736896000*xi_10-46483785515008.                    (4)
```

The coefficient of `c14` in `(3)` is the nonzero constant `-216513000`.
Therefore `(3)` is a global affine graph for `c14`, without localizing at
`Phi`, `eta`, or any other source coordinate:

```text
c14=(E_9-3339765000*c70-898128000*c42)/216513000.       (5)
```

Thus THM-4315's old row-nine hypersurface is not a necessary equation after
the complete weight-fourteen face is restored. It is the zero-face slice of
the graph `(3)`.

## 2. Exact face flags and stochastic visibility

Under `p=t(1+x^2*t)` and `y=x*t*p`, the weight-fourteen row increments are

```text
Delta_14 G_n=
 [binom(7,n-7)c70+binom(6,n-8)c42+binom(5,n-9)c14]
 *x^(2n-14),                                      7<=n<=14. (6)
```

The first three labelled channels are

```text
(h_7,h_8,h_9)=
(c70, 7c70+c42, 21c70+6c42+c14).                        (7)
```

This is an integral unimodular coordinate system, with inverse

```text
(c70,c42,c14)=(h_7,h_8-7h_7,h_9-6h_8+21h_7).           (8)
```

In the fixed normalized chart, THM-4328's row-specific Student--Stein
observer on these three even powers is diagonal with entries

```text
(1,2/5,36/85).                                           (9)
```

It is lossless. By contrast, the two weight-thirteen channels contribute
odd powers, and all their direct Student coordinates vanish. The exact
reconstruction confirms that `c51,c23` do not occur in `(3)`.

This contrast explains the change but is not by itself the bracket proof.
At rows seven and eight the visible weight-fourteen data are first absorbed
by already selected response coordinates:

```text
partial W/partial c70=-13/6,
partial Z/partial c70=-494/81,
partial Z/partial c42=-13/18.                            (10)
```

The new `c14` channel first enters at row nine, where it absorbs the old
scalar equation as in `(3)`.

## 3. Rows through eight remain on the old Hasse gate

The primary and independent certificates rebuild the literal polynomial
`(1)` rather than importing a weight-thirteen or weight-fourteen scout. The
rows through eight give the same projected-depth source ideal as THM-4308:

```text
Delta=896/15,       Theta=512/75,       zeta_3=-3Phi/2. (11)
```

The exact row-eight universes and terminal system remain

```text
pi_8(P_2): 63 x 131, rank 51;
pi_8(P_3): 72 x 204, rank 63;
terminal: 21 x 9, rank 2, pivots (7,8).                  (12)
```

Equation `(11)` is stable because the earlier visible face directions change
the selected `W,Z` responses in `(10)`, not because the new face is absent.

## 4. Row-nine bracket and depth

After `(11)`, the row-nine bracket selector has shape `10 x 7`, rank seven,
and pivots `(0,1,2,3,4,5,6)`. Its unique consistency equation is exactly
`(3)`. Both certificates check coefficientwise that the bracket difference
has no term outside degrees `0,...,9`; `(3)` is the complete capped row
condition, not a scalar projection of an unchecked tail.

After imposing the graph `(5)`, the row-nine projected-depth universes are

```text
pi_9(P_2): 75 x 160, rank 59;
pi_9(P_3): 85 x 251, rank 73.                            (13)
```

The joint terminal system has shape `28 x 10`, rank three, and pivots
`(7,8,9)`. Every projected-depth residual vanishes after its selected
solution. Thus row-nine depth is automatic on the bracket graph.

For a rational positive/hostile control, take

```text
(Phi,eta,alpha_11,xi_10,c70,c42)
  =(1,2,3,5,7,11),
c14=-1056627491057/4920750.                              (14)
```

This point satisfies `(3)`, while its old `E_9` value is nonzero. It proves
directly that the inherited hypersurface `E_9=0` has been absorbed, not
retained under another presentation.

## 5. Stopping boundary and reproduction

THM-4388 proves extinction by row twelve only before the weight-fourteen
face is present. The new channels already change the selected row-seven and
row-eight responses and erase the old row-nine source equation. Nothing here
transports THM-4388's later coefficient graphs, 26-point bracket scheme, or
depth obstruction to `(1)`. The next sharp calculation is the complete
row-ten bracket and projected-depth fibre over the graph `(5)`.

Artifacts:

```text
04-computation/jc2_source_normal_weight14_row9_face_absorption_thm4390.py
05-knowledge/results/jc2_source_normal_weight14_row9_face_absorption_thm4390.out
04-computation/jc2_source_normal_weight14_row9_face_absorption_thm4390_independent_audit.py
05-knowledge/results/jc2_source_normal_weight14_row9_face_absorption_thm4390_independent_audit.out
```

Replay from the repository root:

```text
python3 -B 04-computation/jc2_source_normal_weight14_row9_face_absorption_thm4390.py
python3 -B -O 04-computation/jc2_source_normal_weight14_row9_face_absorption_thm4390.py
python3 -B 04-computation/jc2_source_normal_weight14_row9_face_absorption_thm4390_independent_audit.py
python3 -B -O 04-computation/jc2_source_normal_weight14_row9_face_absorption_thm4390_independent_audit.py
```

Both artifacts import the audited THM-4308/4315 row operators. The independent
audit imports no THM-4390 code and separately organizes the bracket and depth
eliminations. No finite-field argument is used.
