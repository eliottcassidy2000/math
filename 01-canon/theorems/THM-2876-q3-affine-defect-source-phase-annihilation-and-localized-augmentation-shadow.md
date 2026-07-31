---
id: THM-2876
title: "q3 affine-defect source-phase annihilation and localized augmentation shadow"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The nearest affine q3 endpoint-rectangle defect transports the entire
  26-sample Prony atlas by K=10(omega^6-1), preserving the normalized
  THM-2861 Hermitian edge.  Restoring the actual source coefficient
  contributes the inverse phase omega^7 and annihilates the physical
  two-ended boundary termwise on all ten address copies.  K meets the
  Bockstein generator only on the 13-localized cotangent augmentation
  line; it is a trivial-character vertex scalar and does not repair the
  flat q3/q11/q7 seam.  No row exclusion or LRC(14) proof follows.
source: root/lrc-semilinear-rectangle-bridge-2026-07-28
depends_on:
  - THM-2859-horn-collar-q0-hinge-minimal-v4-globalization-and-witt-endpoint-obstruction
  - THM-2868-two-origin-endpoint-projector-and-projective-kummer-carry-descent
  - THM-2874-endpoint-kummer-galois-clutch-and-bockstein-seam-transgression
related:
  - THM-2847-q3-q11-transverse-endpoint-edge-and-e3-realization-horn
  - THM-2861-endpoint-hermitian-edge-holonomy-and-semilinear-clutch-test
script: 04-computation/lrc14_q3_affine_defect_source_phase_annihilation_thm2876.py
output: 05-knowledge/results/lrc14_q3_affine_defect_source_phase_annihilation_thm2876.out
script_sha256: 5124b82985b44dcfe52b6d484661f83494a60cbd6c3b3f53120c9a8ef12d0071
output_sha256: 50bc1589bb1b1d8f452dc39a6d8af79f048df10f26113b9496041b9fe292fdfc
hash_basis: LF-normalized bytes
---

# THM-2876 -- q3 affine-defect source-phase annihilation and localized augmentation shadow

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This is not a physical current construction or an LRC(14) row exclusion.
The companion pins THM-2859, THM-2868 and both THM-2874 evidence layers,
reconstructs the q3 endpoint rectangles from the canonical evaluator, and
works in both certified endpoint fields.

## Verdict

The rank-one q3 address defect carries the complete THM-2868/2861
Hermitian edge **at endpoint-numerator level**, but it vanishes after the
actual source coefficient is restored.

More exactly, let `C_m^-` and `C_m^+` be the endpoint numerators at q3
before and after the physical `66h` displacement.  Then

```text
C_m^+ = omega^6 C_m^-                                      (1)
```

for every one of the 26 lawful raw samples.  The displacement changes
neither Prony node.  Pairing the positive and negative affine-defect fibres
therefore gives

```text
D_m^endpt
 =10(C_m^+-C_m^-)
 =K C_m^-,

K=10(omega^6-1)
  =10(omega^3-1)(omega^3+1) !=0.                         (2)
```

Every local Prony split and offset transport commutes with this scalar.
Thus the left branch remains `chi_3`, the right branch remains trivial, the
projective ratio is unchanged, and the normalized THM-2861 Hermitian edge
still has Fourier support `{0,3,10}`, phase

```text
E_r=omega^3 conjugate(E_r),
```

and thirteen distinct nonzero values.  The unnormalized edge is multiplied
by

```text
K conjugate(K)=100(2-omega^6-omega^7).                   (3)
```

But the actual source coefficient has the inverse physical phase:

```text
P^+=omega^7 P^-,             omega^7 omega^6=1.          (4)
```

Consequently the correctly paired full-current boundary is

```text
10(P^+C_m^+-P^-C_m^-)=0                                  (5)
```

on all 26 samples.  The surviving Hermitian edge is therefore an
endpoint-only coefficient shadow.  It is not a physical full-current edge.

After localization at 13, the scalar `K` lies on the same augmentation
line as THM-2874's `omega^3-1`; equivalently, they have the same cotangent
direction.  Integrally their principal ideals are not equal, because the
factor `10` is not a cyclotomic unit.  This localized contact is not a
Bockstein coupling.  `K` is a constant `F=Q(zeta_91)`-rational,
relative-character-zero vertex scalar.  It leaves the collapsed frequency
seam holonomy equal to one, while the true ancestry/macro-E3 seam has
holonomy `omega^3`.

## 1. The exact affine defect

At the two q3 states the full endpoint rectangles are

```text
E^- = A3 x {0,3,4,5,8,9,10,11,12},
E^+ = A3 x {0,5,6,7,8,9,10,11,12},

A3={0,1,...,9}.                                          (6)
```

Under

```text
g(a,b)=(a,2b),
```

the exact signed boundary is

```text
E^+-g(E^-)
 =1_A3 tensor (delta_12-delta_3).                        (7)
```

The negative target address `b=3` has source preimage `b=8`.
Consequently the two physical ends of `(7)` are

```text
positive:  (a,12) at step 68,     a in A3;
negative:  (a, 8) at step 2,      a in A3.              (8)
```

The literal endpoint builder verifies:

- at step 2, `b=12`, `b=3`, and `b=8` are all present;
- at step 68, `b=12` is present and `b=3` is absent;
- all ten positive addresses and all ten negative-preimage addresses
  restrict to the same respective single physical interval;
- every address in the global Boolean support complement restricts to the
  empty carrier.

This separates three superficially similar signed quantities:

```text
same start carrier:     10C_m^--10C_m^- =0;
same end carrier:       10C_m^+-0       =10C_m^+;
two-ended affine edge:  10C_m^+-10C_m^- =K C_m^-.       (9)
```

Only the last line is the two-ended nearest-affine boundary.

## 2. Why all 26 local Prony charts survive endpoint-only

Write the THM-2868 numerator as

```text
C_m=alpha_L lambda_L^m-alpha_R lambda_R^m.
```

For the displacement `Delta=66h`, exact exponent arithmetic gives

```text
26 R_dil Delta =0                 mod N,
12 R_dil Delta =1092(N/2366)      mod N,
xi^1092=omega^6.                                        (10)
```

Thus

```text
lambda_L^+=lambda_L^-,
lambda_R^+=lambda_R^-,
alpha_L^+=omega^6 alpha_L^-,
alpha_R^+=omega^6 alpha_R^-.                            (11)
```

This is much stronger than merely observing the same recurrence
numerically: translation changes both endpoint weights by one common scalar
and leaves both nodes literally fixed.  Therefore the original local
offsets

```text
(0,0,0,0,1,0,0,0,2,0,0,0,0)
```

remain lawful, and for every formal section `r`,

```text
U_r^D=K U_r,              V_r^D=K V_r.                  (12)
```

It follows immediately that

```text
U_(r+1)^D=omega^3 U_r^D,
V_(r+1)^D=V_r^D,
U_r^D/V_r^D=U_r/V_r.                                   (13)
```

The exact finite-field replay independently performs all 26 raw endpoint
evaluations, the thirteen two-sample splits, and both repaired offset
transports.  It then recomputes the Hermitian DFT and obtains

```text
support={0,3,10},       13 distinct nonzero edge values. (14)
```

The scalar in `(2)` belongs to `F` and is fixed by the relative carry
Galois group.  Hence it carries frequency character zero; the nontrivial
character in `(13)--(14)` still comes entirely from `U_r`.

## 3. Why the full current dies

THM-2859's exact physical phase audit uses the same displacement and the
same source coefficient `P` imported by THM-2868.  It gives

```text
source x-sweep phase:       omega^7,
target endpoint phase:      omega^6.                    (15)
```

Both values are checked directly in the two endpoint fields, and the source
before displacement equals THM-2868's pinned `COMMON_SOURCE`.  The target
address operation does not alter that source coordinate: the probe gates
the same `P^-` and `P^+` independently on each of the ten positive/negative
address pairs and evaluates every restriction rather than merely inserting
a multiplicity ten.  Combining `(1)` and `(15)` gives

```text
P^+C_m^+
 =omega^7 P^- omega^6 C_m^-
 =P^-C_m^-.
```

Thus every term in the positive ten-address fibre equals the corresponding
term in the transported negative fibre.  This is checked on all `10 x 26`
address/sample pairs in both exact endpoint fields.  Equation `(5)` is not
a Fourier or averaging cancellation; it is termwise physical phase
cancellation.

This is the decisive typing boundary.  Freezing `P=P^-` while transporting
only the endpoint numerator produces the nonzero scalar `K`.  Transporting
the actual source and target together produces zero.  A future physical
repair must therefore break or polarize this inverse-phase pairing
lawfully; another endpoint-address permutation cannot do it.

## 4. Exact relation—and nonrelation—to THM-2874

There is a genuine 13-localized algebraic contact:

```text
K/(omega^3-1)=10(omega^3+1).                            (16)
```

In the cotangent quotient

```text
Z[omega]/(13,(omega-1)^2),
```

the ratio in `(16)` is `20=7 mod13`.  Hence

```text
K = 7(omega^3-1)        up to first cotangent order.     (17)
```

So the endpoint affine defect and the true Cech defect occupy the same
one-dimensional augmentation tangent line.  Since
`10(omega^3+1)` is not an integral cyclotomic unit (the obstruction is the
factor `10`), this does **not** assert equality of their integral principal
ideals.  It is a useful 13-local normalization, not a realization.

The map types remain different:

```text
K:
  one constant endpoint-translation coefficient over q3;

THM-2874 Bockstein:
  failure of a global edge gauge around
  q3 -> q11 -> q7 versus q3 -> q7.
```

Multiplying every frequency vertex by `K` cancels from all edge ratios.  The
collapsed section `r=q-3` therefore retains

```text
frequency holonomy =1,
ancestry holonomy  =omega^3.                            (18)
```

The adjacent-frequency Hermitian identity
`E=omega^3 conjugate(E)` does not change `(18)`: its orientation is on the
`r -> r+1` frequency edge, not on the physical q3/q11/q7 macro-E3 triangle.
Moreover the full-current carrier to which a physical seam could attach is
zero by `(5)`.

## 5. The global endpoint complement is the hostile, not the repair

Each q3 rectangle has `90` present and `79` absent endpoint addresses.  On
its own physical interval every address in that 79-point complement has
empty restriction, so its raw coefficient sum is zero.  This Boolean
complement records where the endpoint carrier is absent; it does not create
a complementary physical carrier.

The macro truth is independently:

```text
q3 step 2:   E3,
q3 step 68:  E3,
q7:          not-E3.                                    (19)
```

Thus `(7)` never crosses the macro-E3 boundary.  The 79-address
endpoint-support complement is neither the primary THM-2874 deep-bit
all-nine-safe hostile nor the addendum's full-pattern complement at the
zero chart.  The latter independently supplies a q7 coefficient vertex,
but no q7 coefficient, uncancelled q11 origin, or common-support
Bockstein transporter follows from the affine defect `(7)`.

## 6. Cheapest remaining decisive experiment

The new obstruction removes the pure two-ended translation route.
THM-2874 now executes the coefficient-support half of the second option
below and proves that its formal seam is flat.  The smallest remaining
construction must therefore retain one of the following extra coordinates:

1. a source-fixed polarization that compares the two endpoint numerators
   without simultaneously applying the inverse source phase; or
2. a typed q7 macro-complement transition at the same source, target, and
   endpoint chart, together with one uncancelled q11 origin and retained
   `QA/QAB` ancestry.

For either candidate, the decisive test is no longer whether the
`{0,3,10}` edge exists—it already does endpoint-only.  The test is whether
the resulting **nonzero full-current** seam has holonomy `omega^3`.
Holonomy one is the exact flat hostile.

## Reproduction

```bash
python3 04-computation/lrc14_q3_affine_defect_source_phase_annihilation_thm2876.py
python3 -O 04-computation/lrc14_q3_affine_defect_source_phase_annihilation_thm2876.py
```

The companion contains no executable Python `assert` statements.  Normal
and optimized outputs are byte-identical.  SHA-256:

```text
script  5124b82985b44dcfe52b6d484661f83494a60cbd6c3b3f53120c9a8ef12d0071
output  50bc1589bb1b1d8f452dc39a6d8af79f048df10f26113b9496041b9fe292fdfc
```

The independent audit rederived the affine sign/preimage convention,
all 26 endpoint phases and Prony splits, the uniform common source
coefficient on every `10 x 26` address/sample pair, and the termwise
annihilation.  It also checked the `13`-localized—not integral-unit—scope
of the augmentation comparison.  Normal and optimized replays byte-match
the stored transcript, and the companion contains no executable Python
`assert`.
