---
id: THM-4437
title: "LRC14 all-parity network reduction to three low circuits"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4414/4422 + FINITE-EXACT +
  INDEPENDENTLY AUDITED. Every primitive sorted distinct positive
  ternary-unit triple having no signed relation with coefficient magnitudes
  (1,1,1), (1,1,2), or (1,2,2) has every complete network projection at most
  6/77 and its minimum strictly below 6/77. The only individual equalities
  are named below. THM-4441/4444/4445 subsequently
  classify all low circuits. Arbitrary entry, synchronization, and LRC(14)
  remain open.
source: root mixed-parity continuation + independent referee, 2026-09-06
depends_on:
  - THM-4414-lrc14-six-separated-contact-capacity-collapse
  - THM-4422-lrc14-projection-deficit-and-beatty-row-reduction
related:
  - THM-4434-lrc14-universal-scale-three-network-projection-bound
  - THM-4441-lrc14-signed-122-sharp-ray-closure
  - THM-4444-lrc14-signed-112-sharp-one-ray-classification
  - THM-4445-lrc14-signed-111-sharp-obstruction-classification
coefficient_script: 04-computation/lrc14_all_parity_coefficient_box_thm4437.py
coefficient_output: 05-knowledge/results/lrc14_all_parity_coefficient_box_thm4437.out
coefficient_script_sha256: 36fde3f9e00b094d43927c4a0fd23a723370804649b8748d1171d00297580153
coefficient_output_sha256: 09927dd2733d2f5bfd5000bc6691ae37f00f4854a6b4795d85cc2f2f15c5f0cd
native_script: 04-computation/lrc14_all_parity_low_circuit_reduction_thm4437.cpp
native_output: 05-knowledge/results/lrc14_all_parity_low_circuit_reduction_thm4437.out
native_script_sha256: 0c6c1121caf75502680d03b307d97032669d0f8cdb8c59d7afd140dbc2cf3a03
native_output_sha256: 6707e02b8435ad497a0a16e35cfb4d0341ccc604541fd9336d4c09d9cda2d270
independent_analytic_script: 04-computation/lrc14_all_parity_analytic_chain_thm4437_independent.py
independent_analytic_script_sha256: 0caebf01c71e03e8f862e0e466b3a7c517359986e4a75c920664625afcb18782
independent_coefficient_script: 04-computation/lrc14_all_parity_coefficient_box_thm4437_independent.py
independent_coefficient_script_sha256: 837dd8d7bdebdc9cc7535ea7845d21064498b925134950269575b6067f98548a
independent_raw_script: 04-computation/lrc14_all_parity_raw_head_thm4437_independent.cpp
independent_raw_script_sha256: a931ecec254e5a2ae689d395a3c6b400608656aeebcebf0120440a291ba40c1b
independent_audit: 05-knowledge/results/lrc14_all_parity_low_circuit_reduction_thm4437.md
hash_basis: raw LF repository bytes
audit: >
  PASS. Root and a clean-room referee independently rebuilt the parity-free
  750-pattern coefficient box and the 9,720,930-row height-611 head. The
  referee used raw carrier enumeration, not the native sheet engine or the
  short-relation classifier. Ordinary and optimized bounded replays agree.
  The audit caught and repaired the proposed strict all-coordinate wording:
  three nonlow individual coordinates attain 6/77, while every nonlow
  minimum is strict.
---

# THM-4437 -- LRC14 all-parity network reduction to three low circuits

**PROVED ELEMENTARY RELATIVE TO THM-4414/4422 + FINITE-EXACT +
INDEPENDENTLY AUDITED.** This removes parity from the generic local theorem.
This theorem retains three low-circuit families and does not prove entry into
the scale-three chart. The concurrent [parity-free closure and exact threshold
classification](../../05-knowledge/results/overnight4_20260906_lrc_parityfree_native.md#5-combined-nonadditive-ceiling-and-exact-old-threshold-classification)
now combines this theorem with the proved additive, norm-four, and THM-4441
norm-five results. All local ceilings are closed; `LRC(14)` remains **OPEN**.

## 1. Statement

Let

```text
w=(a,b,c),       1<=a<b<c,       gcd(a,b,c)=1,       3 does not divide abc.
```

No parity condition is imposed. Retain THM-4414's complete carrier set

```text
Lambda(w)={C in Z^3:C dot w=0,
  every C_i nonzero mod 3,
  14|C_i|<3(w_j+w_k) for every {i,j,k}={1,2,3}},
```

and its three exact projections

```text
E_i(w)=sum_(C in Lambda(w)) min(3/(7c),
  [3(w_j+w_k)-14|C_i|]/(14w_jw_k)).                  (1)
```

Suppose `w` has no signed integer relation whose sorted absolute primitive
coefficient triple is one of

```text
(1,1,1),       (1,1,2),       (1,2,2).               (2)
```

Then

```text
E_i(w)<=6/77 for i=a,b,c,       and       min_i E_i(w)<6/77.  (3)
```

The complete individual equality locus in (3) is

```text
E_b(7,16,22)=6/77,
E_c(14,17,22)=6/77,
E_c(4,19,22)=6/77.                                   (4)
```

All other projections in the stated domain are strict. Thus the strict
minimum in (3), rather than strictness of every coordinate, is the sharp
selector statement.

For sorted positive distinct speeds, the first circuit in (2) is exactly the
additive family `a+b=c`. The second is the three norm-four identities from
THM-4434, and the third is the next signed circuit shell. The theorem is a
reduction *to* these families, not a claim about them.

## 2. Inheritance and the repaired boundary

The closest proved mechanism is THM-4434's relation-independent zonotope
slice, discrepancy estimate, carrier intercept, and planar short-relation
argument. Oddness entered that proof in the finite coefficient universe, the
even-`S` cutoff split, and the native speed head. All three restrictions are
removed here.

The canonical hostile is

```text
w=(2,5,7),       E=(22/245,6/49,1/10),
```

whose minimum `22/245` is greater than `6/77`; it lies on the `(1,1,1)`
circuit. Excluding the additive family alone is insufficient:
`(2,11,20)` lies on `(1,1,2)` and has minimum `11/140>6/77`.

The corrected near miss is the stronger proposed wording “every nonlow
projection is strict.” It fails first at `(7,16,22)`, where `E_b=6/77`.
The exact audit finds the other two equalities in (4), but no nonlow
projection above the target and no nonlow minimum equality. The missing
coordinate was which selector projection attains the ceiling; scalarizing
immediately to the minimum hid that harmless boundary.

## 3. The parity-free coefficient theorem

Let `v dot w=0` be a primitive nonzero relation, put

```text
S=||v||_1,       M=max_i |v_i|.
```

THM-4434 Sections 2--4 construct an even decreasing section-width function
for the image of the error cube under

```text
e -> (v dot e,(w cross e)_i/v_i).
```

Its integral is `9(a+b+c)/49`; exact owner residues and rectangle
discrepancy give the slope load

```text
F_v(w)/c <= 6/49+4/(7M).                              (5)
```

Nothing in this construction uses the parity of `w` or `S`. Therefore
`M>=19` gives

```text
F_v(w)/c <=142/931<15/98.                             (6)
```

For `M<=18`, the complete parity-free magnitude universe is

```text
0<=p_1<=p_2<=p_3<=18,
at least two p_i nonzero, gcd(p_1,p_2,p_3)=1,
at most one p_i zero mod 3,
p!=(0,1,1),
p not in {(1,1,1),(1,1,2),(1,2,2)}.                  (7)
```

There is deliberately no even-norm filter. The support-two exclusion would
force equal speeds. Two coefficient residues divisible by three would force
the third one to be divisible by three because every speed is a ternary unit,
contradicting primitivity.

The full box before the three residual exclusions has `750` patterns:
`701` of full support and `49` of support two. Exact rational polygon
clipping and a separate cube-edge construction agree on every full-support
sector; the support-two rectangle formula supplies the boundary check. Every
nonlow pattern satisfies

```text
F_v(w)<=15c/98,                                       (8)
```

with unique sharp coefficient pattern `(1,7,8)`. The three excluded maxima
are respectively `2/7`, `2/7`, and `3/14`, so the exclusions are structural,
not unused bookkeeping.

## 4. Carrier count and the height-611 reduction

The parity-free defect count is the same as in THM-4434. If the allowed
defect list is empty the carrier conclusion is immediate; otherwise its
strict open-interval count, combined with (8), gives

```text
N:=|Lambda(w)| < 15c/98+2S/7+4/3.                    (9)
```

Hence THM-4422's automatic gate `N<=2c/11`, which pays all three projections
strictly, holds when

```text
c >= (308/31)S+4312/93.                              (10)
```

The projected relation lattice and explicit hexagonal `l1` ball from
THM-4434 do not use parity. Planar Minkowski supplies a primitive relation

```text
S<4 sqrt(c/3).                                        (11)
```

If this relation has one of the patterns (2), the speed triple is excluded
by hypothesis. Otherwise (9)--(10) apply. For `S<=57`, the right side of
(10) is at most

```text
56980/93 < 613.
```

For `S>=58`, the difference

```text
3S^2/16-(308/31)S-4312/93
```

is positive at `58` and increasing, while (11) gives `c>3S^2/16`. Thus all
eligible rows with `c>=613` are strict. The only intervening integer `c=612`
is divisible by three and is outside the universe. The proof therefore
reduces exactly to `c<=611`.

## 5. Complete native head and independent raw audit

The head universe is

```text
1<=a<b<c<=611,       gcd(a,b,c)=1,       3 does not divide abc.  (12)
```

The maintained producer reuses THM-4434's audited six-sheet literal interval
engine but removes every parity filter. It computes all three projections
before classifying a row by brute signed permutations of (2). The exact
totals are

```text
rows                                                   9,720,930
low-circuit rows                                          67,534
rows outside the three low circuits                    9,653,396
minimum failures outside the low circuits                      0
outside-low individual projections greater than 6/77          0
outside-low individual projection equalities                   3. (13)
```

The three equalities are exactly (4). Every outside-low minimum is strict;
its largest head value is `3/70` at `(2,11,40)`.

The low-circuit cross-tab is diagnostic, not an extrapolation:

```text
pattern      rows       min>6/77       min=6/77
(1,1,1)     14,220        14,219              0
(1,1,2)     28,438             1              1
(1,2,2)     24,876             0              0.       (14)
```

The lone `(1,1,2)` failure is `(2,11,20)` and its equality is `(1,5,11)`.
The zero failures for `(1,2,2)` are finite evidence only; that all-height
family is not silently included in the theorem.

A clean-room referee independently enumerates every raw carrier from its
Diophantine equation and strict roofs. It imports neither the six-sheet
engine nor the short-relation classifier, and it reproduces (12)--(14), the
three equalities, every parity count, and the two hostile controls. The
coefficient referee independently rebuilds cube-edge sections and normalized
speed-polygon vertices. This supplies separate implementations of both
finite obligations; the full [audit report](../../05-knowledge/results/lrc14_all_parity_low_circuit_reduction_thm4437.md)
records the proof-line checks and frozen hashes.

## 6. Physical consumer and exact scope

THM-4414 gives the physical failure mass

```text
mu(F_w)<=min_i E_i(w)<6/77                             (15)
```

in the theorem's domain. Consequently, for every finite positive-speed body
`C` with

```text
G_C={y in R/Z: ||cy||>=1/14 for all c in C},
```

the same scale-three argument as THM-4434 gives

```text
mu(G_C)>=6/77  implies  G_(3C union {a,b,c}) nonempty. (16)
```

Common ternary-unit dilation of the tail is removed by the corresponding
circle covering map, so the physical conclusion also applies after primitive
reduction. It remains a sufficient local consumer, not an entry theorem.

This theorem isolates the generic mixed-parity obstruction at the three
lowest signed circuits. [THM-4441](THM-4441-lrc14-signed-122-sharp-ray-closure.md)
subsequently proves the sharper `min E<=46/665<6/77` on the whole `(1,2,2)`
family. Only `(1,1,1)` and `(1,1,2)` remain possible local hostiles. The
additive sharp `6/55` theorem does not restore the old target. THM-4444 then
reduces `(1,1,2)` to hostile `(2,11,20)` and boundary `(1,5,11)`. Arbitrary
entry, synchronization, and `LRC(14)` remain **OPEN**.

## 7. Reproduction

```powershell
python -B 04-computation/lrc14_all_parity_coefficient_box_thm4437.py
python -B -O 04-computation/lrc14_all_parity_coefficient_box_thm4437.py
g++ -std=c++20 -O3 -DNDEBUG 04-computation/lrc14_all_parity_low_circuit_reduction_thm4437.cpp -o thm4437
./thm4437 611
```

The coefficient output is byte-identical under normal and optimized Python.
Root also compiled ordinary and optimized native binaries and obtained
byte-identical height-79 transcripts before the full optimized head. All
mathematical comparisons are exact integer or rational comparisons.
