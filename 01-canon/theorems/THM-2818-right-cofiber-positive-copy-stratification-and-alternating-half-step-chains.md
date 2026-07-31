---
id: THM-2818
title: "Right-cofiber positive-copy stratification and alternating half-step chains"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The complete
  THM-2771 right cofiber has exactly twenty-eight nonzero physical
  (clock,target) cells.  After the inherited delayed coefficient selector,
  every cell is a positive sum of equal-length, equal-weight,
  ancestry-identical-within-cell interval copies.  Twenty-five cells contain
  two live copies.  The three target-twelve cells contain respectively
  121,265,254 live copies in alternating half-step chains, explaining every
  formerly opaque multiplier without cancellation.  The coefficient
  Bockstein is additive copy by copy.  There are two global ancestry
  prototypes, and native-factor, carrier, endpoint, target-convolution, and
  root-deck failures prevent a typed common-atom cospan.  No row exclusion or
  LRC(14) conclusion follows.
source: root/right-cofiber-copy-stratification-2026-07-28
depends_on:
  - THM-2771-joint-c7-c13-right-wing-mixed-spectrum-and-commuting-square-no-go
  - THM-2791-full-arm-orbit-transfer-and-lower-central-chord
related:
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
  - THM-2751-root-zero-clutch-mayer-vietoris-wing-shear
  - THM-2754-diagonal-clock-81-label-root-zero-clutch-addendum
  - THM-2782-semantic-arm-right-wing-local-unit-and-endpoint-deck-boundary
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2809-sharp-labelwise-eleven-label-maximum-on-fourteen-marked-rails
  - THM-2814-projective-allocation-square-holonomy-and-idempotent-provenance-no-go
  - THM-2819-nonwrapping-endpoint-conjugacy-odometer-wrap-death-and-sharp-target-eleven-face
script: 04-computation/lrc14_right_cofiber_positive_copy_stratification_thm2818.py
output: 05-knowledge/results/lrc14_right_cofiber_positive_copy_stratification_thm2818.out
script_sha256: 85edac9bb03f1fef198343268f4fb1226bec122d45ded79a049f8fa9a73882a8
output_sha256: 225bd77c27d5972e7dad663e46be3c4c20e2b9449615018773e6822360356a33
secondary_script: 04-computation/lrc14_right_cofiber_positive_copy_stratification_independent_audit_thm2818.py
secondary_output: 05-knowledge/results/lrc14_right_cofiber_positive_copy_stratification_independent_audit_thm2818.out
secondary_script_sha256: 684b8f71b74c83396bb4ecb78e11084cb88b25d343be1f526e21641093e8ea42
secondary_output_sha256: 63e53308d8561749f21a2bc5a1d013b6446d69998ad1488c4bed4ed6e4d95c19
hash_basis: LF-normalized bytes
---

# THM-2818 -- the right-cofiber multipliers count positive physical copies

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2771 records a sparse 7 by 13 right-cofiber coefficient table whose
nonzero entries have primitive multipliers 2, 121, 265, and 254.  A
coefficient table alone does not say whether those numbers arise from
positive multiplicity, cancellation, unequal atoms, or a hidden quotient.
This theorem resolves that ambiguity at the physical interval level.

Every nonzero entry is a literal positive sum of equal copies.  The ordinary
multiplier 2 is a pair of copies.  The three larger multipliers count the live
positions of finite alternating chains.  This gives a copywise explanation
of the coefficient and its mod-13 Bockstein, while retaining the exact
physical data that still prevent a common-atom construction.

## 1. Physical universe and live-copy definition

Keep THM-2771's rail-eight right cofiber, the full semantic section, source
carry 12, translated target carry 6, and clocks e in F_7 and targets t in
F_13.  For each cell, form the source one-sided carrier A, the pulled target
carrier B, their common part M=A intersect B, and

~~~text
R_(e,t)=B minus M.                                         (1)
~~~

Decompose R_(e,t) into maximal weighted open intervals.  Apply the inherited
delayed terminal functional separately on the source interval and its target
translate.  On every interval in all 91 cells the resulting pair is either

~~~text
(C,C),       C=103478815440,                               (2a)

or

(0,0).                                                     (2b)
~~~

There is no unequal pair, second nonzero content, or negative coefficient.
Call an interval live in case (2a) and dead in case (2b).  This distinction
is essential: raw support and coefficient-live support agree in the
ordinary cells and differ at target twelve.

## 2. Complete cell census

Exactly twenty-eight cells are nonzero.  Their target supports are

~~~text
S_1={3,4,5,6,7,8,9,10,11,12},
S_2={2,3,4,5,6,7,8,9,12},
S_3={2,3,4,5,6,7,10,11,12}.                              (3)
~~~

The remaining clocks are empty.  Twenty-five of these cells have

~~~text
2 raw pieces = 2 live + 0 dead.                            (4)
~~~

The three exceptional rows are

| cell | raw | live | dead |
|---|---:|---:|---:|
| (1,12) | 241 | 121 | 120 |
| (2,12) | 528 | 265 | 263 |
| (3,12) | 506 | 254 | 252 |

Every piece has length

~~~text
L=26444880.                                                (5)
~~~

At clock one every piece has weight

~~~text
w_1=27581135604,                                           (6a)
~~~

and at clocks two and three every piece has weight

~~~text
w_2=w_3=27580222516.                                       (6b)
~~~

Thus equal coefficients in a cell come from equal physical interval pieces,
not merely from a post-integration coincidence.

## 3. Why the sporadic multipliers are 121, 265, and 254

Order each exceptional row by left endpoint and put

~~~text
h=T/(2*13^5)=401080680.                                    (7)
~~~

Start a new block exactly when consecutive left endpoints are not h apart.
Every within-block delayed word is

~~~text
1,0,1,0,...                                                (8)
~~~

from a live head, and every interblock left-endpoint jump is 50h.  The exact
block data (length,live,dead) are

~~~text
clock 1: (145,73,72), ( 96,48,48);

clock 2: (143,72,71), (289,145,144), (96,48,48);

clock 3: (143,72,71), (289,145,144), (74,37,37).           (9)
~~~

Consequently

~~~text
121=ceil(145/2)+ceil(96/2),

265=ceil(143/2)+ceil(289/2)+ceil(96/2),

254=ceil(143/2)+ceil(289/2)+ceil(74/2).                   (10)
~~~

Equations (8)--(10) are the mechanism behind the sporadic branch.  The
multiplier is the number of positive live positions in finite parity chains;
it is not a signed or Fourier cancellation.  The first cell at which raw
support is not already uniformly live is exactly (1,12).

## 4. Local ancestry identity and its sharp global boundary

For every cell, construct the literal THM-2791 U- and V-contributor walls.
The hull from the first to last right-cofiber piece lies strictly inside one
U chamber and one V chamber.  Therefore every raw piece in that cell has the
same ordered ancestry-label sets, including the dead pieces in the
exceptional chains.

There are three wall-free chamber pairs but only two label prototypes:

~~~text
clock 1:
  |U|=966606, |V|=28534, |U||V|=27581135604;

clocks 2 and 3:
  |U|=966574, |V|=28534, |U||V|=27580222516.              (11)
~~~

The first lexicographic change is the cell (2,2).  Relative to clock one, U
loses exactly 32 labels and gains none; V is unchanged.  The displayed
THM-2791 path

~~~text
(a,b,e')=(59162,26,56658)                                 (12)
~~~

survives in both prototypes.  Hence “ancestry-identical copies” is a
within-cell theorem, not a claim of one global prototype and not an endpoint
origin.

## 5. Coefficients and the copywise Bockstein

One live copy contributes the seven-clock vector

~~~text
(0,x_e,x_e,x_e,x_e,x_e,x_e),                             (13)
~~~

where

~~~text
x_1=w_1 C=2854063240791928925760,
x_2=x_3=w_2 C=2853968755527296447040.                    (14)
~~~

A dead copy contributes zero.  Therefore the full cell vector is exactly
the live count times (13), proving the four THM-2771 multipliers by positive
addition.

The common factorization is

~~~text
g_0=5905329039529920,

x_1=g_0*483303,        x_2=x_3=g_0*483287,

gcd(483303,483287)=1.                                    (15)
~~~

Both x_e and g_0 have 13-adic valuation one.  Define the copy Bockstein by
beta_e=(x_e/13) mod 13.  Then

~~~text
beta_1=9,             beta_2=beta_3=2,                   (16)

beta_(e,t)=m_(e,t) beta_e mod 13.                         (17)
~~~

Thus the ordinary totals are 5 at clock one and 4 at clocks two and three;
the three target-twelve totals are respectively 10, 10, and 1.  This proves
that the Bockstein is additive copy by copy in all twenty-eight cells.

## 6. The eleven-label support spine

The union of (3) is

~~~text
S_1 union S_2 union S_3={2,3,...,12}.                    (18)
~~~

Its clock-fibre multiplicity is

~~~text
m(t)=2 for t in {2,8,9,10,11},
m(t)=3 for t in {3,4,5,6,7,12},                          (19)

5*2+6*3=28.
~~~

This explains the cell count as a two-or-three-sheeted cover of one
eleven-label spine.  Set (18) coincides with the positive residual label set
proved by THM-2809.  That coincidence is exact support combinatorics only:
THM-2809 varies label configurations
against a full marked source, whereas (1) is a fixed physical right cofiber
after semantic and coefficient filters.  No common mechanism or map follows.

Target twelve is the only label supported at all three clocks whose raw/live
split is nontrivial.  This makes it the precise boundary at which a proposed
label-conjugacy or common-atom bridge must retain the finite chain ends rather
than only the multiplier.

## 7. Why equal copies still do not form a common physical cospan

The selected ordinary cell (e,t)=(1,4) gives a minimal hostile.  Let I be
the common interval and J_-,J_+ the two live right-cofiber copies:

~~~text
I   =(142004992589460,142005019034340),
J_- =(142004190428100,142004216872980),
J_+ =(142082000080020,142082026524900).                   (20)
~~~

In factor order

~~~text
(E3,clock,q1,q2,c2,c3),                                   (21)
~~~

the four masks (source native, source pulled, target native, target
adjacent-pulled) are

~~~text
I:   111111, 111111, 111111, 111111;
J_-: 011111, 111111, 111111, 011111;
J_+: 011101, 111111, 111111, 011101.                     (22)
~~~

Both cofiber copies fail native E3; J_+ also fails native c2.  Across the
thirteen carrier twists, I has delta-zero support on both sides, whereas
each J has empty source support and delta-zero target support.

The source/target endpoint-mask zero/one counts are

~~~text
I:   (88,81)/(88,81),
J_-: (169,0)/(88,81),
J_+: (79,90)/(70,99).                                    (23)
~~~

No nontrivial endpoint translation identifies I with either J; the identity
matches only the target mask of J_-.  These failures occur before any
Fourier transform.  THM-2806's independently audited target-twelve addendum
gives the complementary global warning: its selected target-twelve simplex
has empty common carrier before address restriction.

Therefore equal length, weight, content, local ancestry, and additive
Bockstein do not supply a common atom or a typed boundary map.  A successful
bridge must add, and prove preservation of, native factors, source carrier,
endpoint origin, target convolution, and root-deck transport.

## 8. Exact proof and independent audit

Run

~~~text
python3 04-computation/lrc14_right_cofiber_positive_copy_stratification_thm2818.py
python3 -O 04-computation/lrc14_right_cofiber_positive_copy_stratification_thm2818.py

python3 04-computation/lrc14_right_cofiber_positive_copy_stratification_independent_audit_thm2818.py
python3 -O 04-computation/lrc14_right_cofiber_positive_copy_stratification_independent_audit_thm2818.py
~~~

Each normal/optimized pair byte-matches its stored output; both scripts have
zero Python assertion nodes.  The primary reconstructs all twenty-eight
cells plus the selected factor/carrier/endpoint sidecars.  The independent
audit imports only four hash-pinned canonical constructors, scans all 91
cells rather than receiving the support table, discovers the exceptional
blocks and ancestry prototypes, and independently recomputes the coefficients,
primitive factors, valuations, and Bocksteins.

Exactly proved are the positive-copy decomposition, the alternating-chain
formulas, local ancestry identity with its two-prototype boundary, copywise
coefficient and Bockstein laws, the eleven-label support spine, and the
selected typed hostile.  Not proved are a common allocated atom, a chart-to-
cofiber map, a target convolution action, root/Cech transport, row exclusion,
or LRC(14).

**QED.**
