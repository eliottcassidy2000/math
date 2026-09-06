---
id: THM-4410
title: "Source-normal least-weight-twenty row-thirteen affine continuation"
status: >
  PROVED FINITE-ROW RELATIVE TO THM-4403 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. In the restricted source obtained by adjoining kappa20*p*y^6 to
  THM-4403's two-channel tail, row thirteen has one active source condition
  and kappa20 pays it through a nonzero constant pivot. The source projection
  remains the global A^4 over (Phi,eta,alpha11,c51), and every terminal fibre
  is A^9. The complete exact-t-valuation-thirteen monomial response has rank
  one, all odd channels vanish, and p*y^6 is its least-weight visible channel.
  This is not the complete weight-at-most-twenty family and proves no chart or
  seam entry, all-row lift, termination, Keller pair, JC(2), or DC(2) claim.
source: root + row13_late_scout + row13_independent_audit / JC2 and arXiv continuation session, 2026-09-03
audit: >
  PASS. The primary certificate extends the frozen THM-4403 state and checks
  the row-thirteen bracket, complete exact-valuation response diagonal,
  projected depths, affine graph, dense and boundary controls, and an
  off-graph hostile. A clean-room certificate imports only the audited
  THM-4308/4315 row operators, reconstructs the complete weight-at-most-
  fourteen prefix and rows eleven/twelve, and then independently derives row
  thirteen. Its 121 exact checks reproduce every matrix, pivot, residual,
  response ratio, graph, depth rank, and control. Normal, optimized, and
  fixed-hash-seed clean-room streams byte-match the frozen LF transcript. No
  finite-field inference is used.
depends_on:
  - THM-4403-source-normal-two-channel-weight-eighteen-row-twelve-affine-continuation
  - THM-4399-source-normal-row-eleven-late-weight-eighteen-response-absorption
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
related:
  - THM-4395-source-normal-weight-fourteen-row-ten-global-affine-absorption
  - THM-4328-seam-covariant-student-stein-face-visibility
  - THM-2044-explicit-rank-two-poisson-counterexample-by-symplectic-suspension
primary_script: 04-computation/jc2_source_normal_weight18_row13_continuation_scout_s616.py
primary_output: 05-knowledge/results/jc2_source_normal_weight18_row13_continuation_scout_s616.out
primary_script_sha256: d194dcba2203b2c90fdb452c82cfc192d55beaf5d4bac3e412c5fca7ac921eb3
primary_output_sha256: c447a33e54715acadd54c0ecc4c3aa2a7d462ddea334bfb408f6435c9e64c6cf
independent_audit_script: 04-computation/jc2_source_normal_weight18_row13_continuation_independent_audit_s616.py
independent_audit_output: 05-knowledge/results/jc2_source_normal_weight18_row13_continuation_independent_audit_s616.out
independent_audit_script_sha256: 9ae70b915847fa7842df8c0d9fa2260550459b7677ade40622522d0f02fa9c1d
independent_audit_output_sha256: cd36174f459caa0da5e660a68f4eba1ff0715f995178f9092b36410ffef38351
hash_basis: raw LF bytes
---

# THM-4410 -- a least-weight-twenty channel carries the restricted tail through row thirteen

**PROVED FINITE-ROW RELATIVE TO THM-4403 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. THIS THEOREM TREATS ONE SELECTED WEIGHT-TWENTY CHANNEL ON TOP OF THE
TWO SELECTED WEIGHT-EIGHTEEN CHANNELS OF THM-4403. IT DOES NOT TREAT THE
COMPLETE WEIGHT-FIFTEEN-THROUGH-TWENTY SOURCE, CHART OR SEAM ENTRY, AN
ALL-ROW LIFT, POLYNOMIAL TERMINATION, A KELLER PAIR, `JC(2)`, OR `DC(2)`.**

## 1. Restricted universe and first visibility

Work over a characteristic-zero field in the fixed source-normal chart of
THM-4395. Retain THM-4403's source

```text
H_tail=lambda18*p^3*y^4+nu18*y^6
```

over the complete residual-weight-at-most-fourteen Hamiltonian `H_14`, and
adjoin exactly

```text
kappa20*p*y^6.                                          (1)
```

Here `p=t(1+x^2*t)` and `y=x*t*p`. Hence

```text
[t^r](p*y^6)=0 for r<13,          [t^13](p*y^6)=x^6,   (2)
```

and its residual weight is `20`. It changes none of the equations through
row twelve. By THM-4403, before row thirteen the source is already the global
affine graph

```text
A^4_(Phi,eta,alpha11,c51),                              (3)
```

with constant-denominator solved coordinates

```text
(c23,c70,lambda18,nu18)                                (4)
```

and terminal fibre `A^9`.

## 2. The unique active row-thirteen equation

After substituting all graphs in `(4)`, the degree-capped row-thirteen
bracket matrix has shape `14 x 9`, rank nine, and pivots `(0,...,8)`. Its raw
cokernel has dimension five, but on this inherited restricted source only one
primitive source residual survives. With a fixed primitive normalization the
condition after adding `(1)` is

```text
R_13+C_13*kappa20=0,                                   (5)
C_13=146361421124462229507072000000000,                (6)
```

where the exact 29-term quartic is

```text
R_13=
-503004647918674504953858834703125*Phi^4
+281399664232361822810302158750000*Phi^3*alpha11
+25632935297351247505641297187500*Phi^3*c51
+168646406739873822576865026562500*Phi^3*eta
+997371179190601416579226275000000*Phi^2*alpha11^2
+240934798305225203288209537500000*Phi^2*alpha11*c51
+507884048762601333843180109687500*Phi^2*alpha11*eta
+14443391196586112356470459375000*Phi^2*c51^2
+60840212468567704280567156250000*Phi^2*c51*eta
+183371558796978972752388813750000*Phi^2*eta^2
+31103763798486192216057642285268992000*Phi^2
+240934798305225203288209537500000*Phi*alpha11^2*eta
+28886782393172224712940918750000*Phi*alpha11*c51*eta
+1058211391659169120859793431250000*Phi*alpha11*eta^2
+68958748086694688156504823062691840000*Phi*alpha11
+120467399152612601644104768750000*Phi*c51*eta^2
+8907407707584660309945846625812480000*Phi*c51
+241125556732625043168769406250000*Phi*eta^3
+9821967221983641126420067415654400000*Phi*eta
+14443391196586112356470459375000*alpha11^2*eta^2
+3921152972505592649000878080000000*alpha11^2
+5294669582455545986439782400000000*alpha11*c51
+120467399152612601644104768750000*alpha11*eta^3
+8905711617228736291829492932116480000*alpha11*eta
+272733666678685404521280000000000*c51^2
+8156495129025030884010270720000000*c51*eta
+249342794797650354144806568750000*eta^4
+34494752372782415324210996992081920000*eta^2
+534024156288269955770790987227002574995456.          (7)
```

The pivot `(6)` is a nonzero rational integer, so `(5)` gives the global
constant-denominator graph

```text
kappa20=-R_13/C_13.                                    (8)
```

There is no localization in `Phi`, `eta`, or another source parameter.
Substitution of `(8)` into every literal row-thirteen bracket coefficient
gives zero. Replacing the right side of `(8)` by `-R_13/C_13+1` leaves the
exact nonzero residual `C_13`, providing an off-graph hostile control.

## 3. Complete exact-valuation-thirteen response line

Every monomial in `p,y` of exact `t`-valuation thirteen is

```text
p^13, p^11*y, p^9*y^2, p^7*y^3,
p^5*y^4, p^3*y^5, p*y^6,                              (9)
```

with residual weights `(26,25,24,23,22,21,20)`. Their leading row is
`(1,x,...,x^6)`. In the active primitive quotient their exact responses are

```text
(545467333357370809042560000000000,
 0,
 130912160005768994170214400000000,
 0,
 102452994787123560654950400000000,
 0,
 146361421124462229507072000000000).                  (10)
```

Thus the response matrix has rank one. Relative to `p*y^6`, the ratios are

```text
(805/216,0,161/180,0,7/10,0,1).                       (11)
```

These are the normalized even moments of THM-4315's row-thirteen
Student--Stein law; parity kills every odd channel. Consequently `p*y^6` is
the least-weight visible monomial on the complete exact-valuation-thirteen
diagonal. This is a statement about that diagonal, not about the full
unrestricted source at weight twenty.

## 4. Projected depth and affine geometry

After `(8)`, the exact projected-depth universes are

```text
pi_13(P_2): 133 x 308, rank 97,
pi_13(P_3): 147 x 491, rank 117.                       (12)
```

Their joint terminal system has shape `66 x 14`, rank five, and pivots
`(9,10,11,12,13)`. It imposes no additional source condition and leaves nine
free tangent coordinates. Combining `(8)` with the inherited graphs proves

```text
source projection through row thirteen =
  A^4_(Phi,eta,alpha11,c51),
every terminal fibre = A^9.                            (13)
```

The certificates evaluate all original bracket and depth coefficients at a
dense point and at a `Phi=eta=0` boundary point. Both controls vanish after
the graph substitutions; no denominator changes or boundary rank assumption
is used.

## 5. Mechanism and boundary

Rows eleven through thirteen now exhibit the same useful late-channel rule:
a Hamiltonian monomial can be chosen by **first visibility**, so it preserves
the solved prefix and pays the next Student-visible compatibility through a
constant pivot. The rank-one response and parity zeros show that this is not
generic coefficient abundance. At row twelve an inherited lower-row
coordinate was still needed for a Student-invisible side condition; at row
thirteen no analogous side condition survives after projected depth.

The exact universe is only

```text
H_14+lambda18*p^3*y^4+nu18*y^6+kappa20*p*y^6.          (14)
```

It excludes all other residual-weight-fifteen-through-twenty terms. It gives
no converse from projected depth to full `B_2` membership, no compatible
infinite lift, no termination, and no route into a planar Keller chart. The
planar Jacobian and Dixmier conjectures remain open.

## 6. Reproduction

Replay from the repository root:

```text
python3 -B 04-computation/jc2_source_normal_weight18_row13_continuation_scout_s616.py
python3 -B -O 04-computation/jc2_source_normal_weight18_row13_continuation_scout_s616.py
PYTHONHASHSEED=9173 python3 -B 04-computation/jc2_source_normal_weight18_row13_continuation_scout_s616.py
python3 -B 04-computation/jc2_source_normal_weight18_row13_continuation_independent_audit_s616.py
python3 -B -O 04-computation/jc2_source_normal_weight18_row13_continuation_independent_audit_s616.py
PYTHONHASHSEED=314159 python3 -B 04-computation/jc2_source_normal_weight18_row13_continuation_independent_audit_s616.py
```

Each certificate is deterministic across its three modes. The primary
certificate performs `44` exact checks; the clean-room audit performs `121`.
Their raw LF hashes are recorded in the frontmatter.
