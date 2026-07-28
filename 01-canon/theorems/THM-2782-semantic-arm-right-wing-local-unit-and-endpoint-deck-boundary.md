---
id: THM-2782
title: "Semantic-arm right-wing local unit and endpoint-deck boundary"
status: >
  RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
  INDEPENDENT AUDIT.  Three explicit rail-eight right-wing packets retain
  the physical E3 -> Q_(3,{1,2}) target word, both relative-present cuts,
  root one, carry six, and exact positive weighted mass.  At clock one the
  three-arm
  contrast is an augmentation-quotient unit, and the intrinsic THM-2771
  Bockstein decoder has thirteen nonzero target coordinates.  The arms are
  a literal three-point +13 central-digit segment modulo 169, but the
  coefficient fails to descend around that quotient and the adjacent
  low-digit cospan alternates zero/nonzero endpoint coefficients.  No full
  physical H13 action, allocated THM-2625 endpoint, chart/root transition,
  row exclusion, or LRC(14) conclusion is proved.  Until independent
  promotion, no proved result may depend on this candidate.
source: root/semantic-arm-right-wing-attachment-2026-07-28
depends_on:
  - THM-2712-semantic-following-congruence-lock-and-address-coboundary-descent
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
  - THM-2751-root-zero-clutch-mayer-vietoris-wing-shear
  - THM-2771-joint-c7-c13-right-wing-mixed-spectrum-and-commuting-square-no-go
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
related:
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2697-filtered-affine-handoff-germ-category-and-base-signature-holotopy-boundary
script: 04-computation/lrc14_semantic_arm_right_wing_central_digit_thm2782.py
output: 05-knowledge/results/lrc14_semantic_arm_right_wing_central_digit_thm2782.out
script_sha256: 7fbc6bb1ec303ded98eaad6e5d8205eb3d247258ada32b6f9904fc439ebb11fb
output_sha256: 13f570d63f212808171cecdb4d8f9aa41884fbdc7ed571dbfe27122b412fadc4
hash_basis: LF-normalized bytes
---

# THM-2782 -- the right wing has a physical semantic central segment, but it does not descend

**RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
INDEPENDENT AUDIT.**

THM-2771 proves a complete coefficient Bockstein for the physical
`C_7 x C_13` right wing, but leaves its target convolution unattached to a
semantic endpoint.  THM-2779 independently identifies the minimal
`169`-point Heisenberg carrier, but has no physical semantic current on that
carrier.

This theorem constructs the missing **local semantic attachment**.  Three
explicit physical packets retain the actual terminal word, carry, root half,
and relative-present cuts.  Their three-arm target contrast is exactly the
augmentation unit anticipated by THM-2751, and THM-2771's intrinsic decoder
fills all thirteen target coordinates.  The arms also lie on one central
`+13` address segment modulo `169`.

The construction does not close the global carrier problem.  The coefficient
depends on a higher address digit and therefore does not descend to the
`169`-point quotient.  The adjacent low-digit cylinders glue as sets while
their endpoint coefficients alternate between zero and nonzero.  These two
exact hostiles isolate the remaining holotopy square.

## 1. Exact physical packet

Put

```text
p=13,                    R=p^6=4,826,809,
T=297,836,897,838,480,   tau=7/R,
z_*=46,873,542,509,301 / 100,360,982,066,072,
epsilon=1 / 100,360,982,066,072.                         (1)
```

For an integer address `n`, define

```text
x_n={z_*+7n/R},          y_n=x_n+7/R,
C_n=(y_n-epsilon,y_n+epsilon).                           (2)
```

For physical clock `e` and lawful target labels `(sigma,t)`, let

```text
A_(e,sigma,t)
 =E3 intersect clock_e
   intersect q1-safe(-sigma/13) intersect c2-safe(+sigma/13)
   intersect q2-safe(-t/13)     intersect c3-safe(+t/13),

B=T_tau^(-1)A,           M=A intersection B,
L=A\B,                   Rwing=B\A.                      (3)
```

The coefficient used below is obtained in the following order:

1. retain both relative-present complements at terminal clock `1`;
2. push `Rwing` to the target coordinate;
3. retain the `c3` root-one half;
4. restrict to the physical cylinder `C_n`;
5. take the actual delayed carry-six/root-one coefficient of
   `D^6 Q_(3,{1,2})`.

Every operation precedes integration.  Thus the result is not a
postprocessed support table.

Use the exact base addresses

```text
n_1=3,454,614,          n_2=1,393,828,          n_3=708,364.   (4)
```

For clock `e`, the three arms have addresses

```text
n_e, n_e+13, n_e+26.                                    (5)
```

All nine target cylinders in `(5)` retain:

```text
rail eight;
the physical present clock e;
E3 -> Q_(3,{1,2});
both relative-present complements;
lawful target labels;
target root one;
actual predecessor carry six.                            (6)
```

The first source cylinder has semantic record `(0,None)` and the other two
already have the `E3` record.  The common endpoint asserted here is the
**target** semantic endpoint.

## 2. The exact three-arm tables

After dividing each clock by its positive common content, the arm-by-target
tables are

```text
clock 1:
  0001111111111
  0000000000001
  0000000000001

clock 2:
  0011111111001
  0000000000001
  0000000000001

clock 3:
  0011111100111
  0000000000001
  0000000000001.                                        (7)
```

Columns are `t=0,...,12`.  Every displayed `1` is a literal positive
physical cell.  The raw contents are

```text
c_1=790,161,473,087,466,480,
c_2=c_3=790,135,314,376,327,920,                         (8)
```

with

```text
v_13(c_e)=1
```

for all three clocks.  Every positive cell at clock `e` has exact weighted
mass

```text
c_e/R.                                                    (9)
```

In reduced form these weighted masses are

```text
c_1/R=60,781,651,775,958,960/371,293,
c_2/R=c_3/R=60,779,639,567,409,840/371,293.              (10)
```

## 3. The clock-one contrast is an augmentation unit

At clock one the exceptional column `12` is common to all three arms.  Apply
the integral three-arm augmentation projection

```text
Q3(v)_i=3v_i-(v_0+v_1+v_2).                              (11)
```

It kills that common column.  The arm-zero/arm-one difference is

```text
W=z^3+z^4+...+z^11.                                      (12)
```

In `Z[C_13]`,

```text
W(z)(z^2+z^6+z^10)=delta_0+2N_13,
N_13=1+z+...+z^12.                                       (13)
```

Hence `W` is an integral unit after quotienting the norm mode.  This is the
positive inverse anticipated abstractly in THM-2751, now attached to actual
target-semantic right-wing cells.

The intrinsic coefficient Bockstein has scalar

```text
(c_1/13) mod13=2.                                        (14)
```

Convolution with THM-2771's decoder

```text
K_beta=(3,5,5,5,7,12,2,8,2,9,8,11,0)                    (15)
```

gives the arm-zero `Q3` vector

```text
(11,9,2,12,10,9,10,4,3,5,6,1,12) in F_13^13.           (16)
```

Every entry is nonzero.  The other two arm vectors are equal, and the sum
of all three projected arm vectors is zero pointwise.  Thus `(16)` is a
physical semantic-arm attachment of the THM-2771 **coefficient**
Bockstein.  The signed projection and decoder are not themselves a positive
Boolean factor insertion.

## 4. The arms form a central-digit segment

Write

```text
Omega=Z/169Z,              n=v+13w,       v,w in F_13.    (17)
```

After choosing the standard digits, the THM-2779 endpoint Heisenberg model
is

```text
X:n ->14n,                 (v,w)->(v,w+v),
Y:(v,w)->(v+1,w),          carry-suppressed low digit,
Z:n ->n+13,                (v,w)->(v,w+1),                (18)
```

and exact exhaustion gives

```text
XY=ZYX,                    [X,Y]=Z.                       (19)
```

All three base addresses satisfy

```text
n_e mod169=85=7+13*6.                                    (20)
```

Consequently the three arms `(5)` are

```text
(v,w)=(7,6),(7,7),(7,8),                                 (21)
```

the first three points of one central `Z`-direction segment.  All target
endpoints on this segment retain the same semantic word and carry.

This is an exact physical lift of a **three-point central-direction
segment** in the address quotient.  It is not yet an identification with
one allocated THM-2625 endpoint-origin orbit: that requires the missing
same-ancestry current intertwiner.

## 5. All central colours survive, but quotient descent fails

Along the thirteen consecutive physical lifts

```text
n_1+13j,                    j=0,...,12,                   (22)
```

all target cylinders retain the semantic word, carry six, and low digit
`v=7`.  In this lifted central order the primitive right-wing table is

```text
j=0:       0001111111111,
j=1,...,6: 0000000000001,
j=7,...,12:0000000000000.                                (23)
```

Thus target columns `3,...,11` have central profile `delta_0`, while column
`12` has profile `1_{0,...,6}`.  Every active raw target column has all
twelve primitive central characters.  For the displayed integral lift of
`K_beta`, every one of the thirteen decoded target columns is nonconstant,
so cyclotomic irreducibility gives all twelve primitive central characters
in every decoded column.

The coefficient nevertheless does not descend to `Omega`.  The next lift
`j=13` has the same residue `(v,w)=(7,6)` as `j=0`, but

```text
j=0:   (target-3,target-12)=(c_1,c_1),
j=13:  (target-3,target-12)=(0,0).                       (24)
```

The high address digit is therefore load-bearing.  Equations `(21)--(24)`
give a physical charged central segment, not a cyclically `Z`-equivariant
current on the `169` endpoint-origin model.

## 6. The adjacent-base cospan glues support, not current

The cylinders obey the literal identity

```text
target cylinder of base n = source cylinder of base n+1. (25)
```

The bases

```text
6714 -> 6715 -> 6716                                      (26)
```

have digits

```text
(6,9)->(7,9)->(8,9).                                     (27)
```

They form a low-digit `Y` segment, not a central `Z` segment.  On all three
arms the endpoint coefficients are

```text
edge 6714->6715, target label 3:
  target Rwing coefficient =0,
  source M coefficient     =c_1;

edge 6715->6716, target label 12:
  target Rwing coefficient =c_1,
  source L coefficient     =0.                           (28)
```

The positive side in each row has weighted mass `c_1/R`.  The open cylinders
glue
exactly, but the carry/root endpoint conventions alternate.  Hence `(25)`
is a sharp support-level Mayer--Vietoris cospan with no same-coefficient
transition.

The full cyclic successor is also not the carry-suppressed `Y` at the
low-digit wrap.  At `v=12` it pays one central `Z` correction, exactly as in
THM-2657 and THM-2779.  Full factor covariance under `X,Y,Z`, including the
closing high-digit edge, remains open.

## 7. Three label hostiles

The following coordinates are distinct:

| object | coordinate |
|---|---|
| physical semantic address | `n in Z/13^6Z` |
| lawful target packet | `(sigma,t) in F_13^2` |
| endpoint-origin model | `(v,w) in F_13^2` |

At clock one,

```text
(sigma,t,arm)=(0,3,0): coefficient c_1,
(1,3,1):               coefficient 0,
(0,4,1):               coefficient 0.                   (29)
```

Thus central arm translation is neither `sigma->sigma+1` nor `t->t+1`.

There is also a tempting false numerical identification:

```text
Omega low digit v=7,       predecessor carry=6=-7 mod13. (30)
```

THM-2542's marker `a` is a root-deck transition; predecessor carry six is
THM-2640's carry coordinate; `v` is the mod-`169` address digit; and the
selected packet has root one.  Congruence `(30)` supplies no map among these
typed labels.

## 8. Exact source, target, and remaining square

The result has the following invoice:

```text
source:
  the physical rail-eight right-wing interval packet on the cylinders (2);

target:
  one target-semantic E3 -> Q_(3,{1,2}) arm contrast with all target
  Bockstein residues and all lifted central characters;

map:
  physical three-arm restriction -> Q3 -> divide by 13 -> K_beta;

preserved:
  rail, physical clock, target word, both relative-present cuts, lawful
  target labels, root one, carry six, exact positive weighted mass, arm
  address,
  and the lifted central digit;

destroyed or absent:
  positivity after Q3/K_beta, high-digit descent, full H13 covariance,
  allocated THM-2625 endpoint origin, THM-2542 marker/root transition,
  and common source-owner endpoint current.                              (31)
```

The cheapest decisive next test is no longer another local support scan.
It is one of:

1. a high-digit filtration on which the charged profile `(23)` becomes a
   lawful graded `Z` current;
2. a common-current map that sends the physical address segment `(21)` to
   the endpoint-origin central cycle without identifying `(29)`;
3. a corrected `Y` cospan whose coefficient, not only its cylinder, glues
   across `(25)`.

## 9. Exact companion and scope

Run

```bash
python 04-computation/lrc14_semantic_arm_right_wing_central_digit_thm2782.py
python -O 04-computation/lrc14_semantic_arm_right_wing_central_digit_thm2782.py
```

Both modes byte-match

```text
05-knowledge/results/lrc14_semantic_arm_right_wing_central_digit_thm2782.out.
```

The companion uses exact integers and `Fraction`s.  It reconstructs the
physical interval carriers and delayed coefficients, verifies every label
in `(6)`, the tables `(7)`, contents and weighted masses `(8)--(10)`, the
group-ring
identity `(13)`, the decoded vector `(16)`, all `169` Heisenberg digit
identities, the full thirteen-lift table `(23)`, the descent hostile `(24)`,
both cospan mismatches `(28)`, and the target-label hostiles `(29)`.

No full physical Heisenberg action, quotient-descended current,
carrier-allocation square, THM-2542 root transition, common source-owner
endpoint, scalar-row exclusion, or LRC(14) conclusion is proved.

**Awaiting independent audit; not QED.**
