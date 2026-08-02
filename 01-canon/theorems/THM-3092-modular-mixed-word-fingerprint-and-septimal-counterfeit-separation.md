---
id: THM-3092
title: "Modular mixed-word fingerprint and septimal counterfeit separation"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For exact-order marked generators s^2=c^3=1, every epimorphism to
  the natural four-point S4 has ord(sc)=4, every epimorphism to AGL_2(F_3)
  has ord(sc)=8, and every epimorphism to PGL_2(F_8) has ord(sc)=7 or 9.
  The corresponding inner marking counts are 1,2,6.  Only the S4 word has a
  nonzero half-power in the affine translation kernel V4.  At degree nine,
  order 7 is the split C7 torus while order 9 is nonsplit; neither degree nor
  the C7 scalar fibre chooses a modular marking.  A sharp AGL_2(F_3) hostile
  shows that the mixed word still cannot detect the affine origin.  No
  quartic realization, Keller, tree, or LRC consequence is asserted.
source: root-modular-mixed-word-2026-08-01
depends_on:
  - THM-3090-affine-projective-prime-power-handshake-and-septimal-counterfeit
related:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks
  - THM-2775-modular-s4-to-weyl-d3-generator-frame-and-affine-parity-blindness
  - THM-3088-punctured-projective-direction-algebra-and-exceptional-parity-saturation
script: 04-computation/modular_mixed_word_fingerprint_thm3092.py
output: 05-knowledge/results/modular_mixed_word_fingerprint_thm3092.out
script_sha256: 29aa3f7c87c4c1a188517d0612766da11c9ad7a6a9c5e959320dfd329a9247e9
output_sha256: 68053b52a425db17227807d1c2c8ec867d8783a0d03752fee434a2933b4b93a3
hash_basis: LF-normalized bytes
---

# THM-3092 -- modular mixed-word fingerprint and septimal counterfeit separation

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. The marked question

Put

```text
Gamma=<s,c | s^2=c^3=1>=C2*C3=PSL_2(Z).                    (1)
```

THM-2768 proves that the four-point binary/ternary clutch is the spherical
marked quotient `Delta(2,3,4)=S4`.  THM-3090 proves that equality of the
natural affine/projective degrees has one further prime-power solution:

```text
AGL_2(F_3) on F_3^2,       PGL_2(F_8) on P1(F_8),          (2)
```

both of degree nine.  Degree forgets the first word in which the two free
factors interact.  For a marked pair `(s,c)`, define

```text
w=sc.                                                       (3)
```

The order of `cs` is the same, since `cs=s(sc)s`.  Replacing `c` by its
inverse is included in the complete marked atlas below.  Thus `(3)` is an
orientation-safe mixed-word test.

## 2. Complete exact pair atlas

In each natural faithful permutation group, range over every exact-order
pair

```text
ord(s)=2,                 ord(c)=3,                         (4)
```

and record `(ord(sc), |<s,c>|)`.  The exact tables are

```text
S4 (either four-point chart):
  #ord2=9, #ord3=8;
  (2,6):24, (3,12):24, (4,24):24.                          (5)

AGL_2(F_3):
  #ord2=45, #ord3=80;
  (2,6):576,
  (6,6):144, (6,18):720, (6,54):864,
  (8,48):432, (8,432):864.                                 (6)

PGL_2(F_8):
  #ord2=63, #ord3=56;
  (2,6):504, (7,504):1512, (9,504):1512.                  (7)
```

Consequently the epimorphism rows are exactly

```text
Gamma -> S4:          ord(w)=4,       24 marked epimorphisms;
Gamma -> AGL_2(F_3):  ord(w)=8,      864 marked epimorphisms;
Gamma -> PGL_2(F_8):  ord(w)=7 or 9, 3024 marked epimorphisms. (8)
```

The sets of possible generating-word orders in `(8)` are pairwise disjoint.
Thus the first mixed word separates the true four-point clutch from both
degree-nine counterfeits and separates the two counterfeit groups from each
other.  In triangle language, the three targets occur respectively through
marked quotients of

```text
Delta(2,3,4),       Delta(2,3,8),       Delta(2,3,7) or Delta(2,3,9). (9)
```

Only the first equality is the spherical presentation itself; the other
triangle groups are infinite and `(9)` asserts finite quotients, not
isomorphisms.

## 3. Marking multiplicity

All three target groups have trivial centre.  Simultaneous conjugation on a
generating pair therefore has trivial stabilizer: an element fixing both
generators centralizes the generated group.  Dividing the epimorphism counts
in `(8)` by the target orders gives

```text
S4:           24/24 = 1 inner marking class;
AGL_2(F_3):  864/432 = 2 inner marking classes;
PGL_2(F_8): 3024/504 = 6 inner marking classes,             (10)
```

with three projective classes of word order seven and three of word order
nine.  This makes the exceptional clutch sharper: after an identification
of the four natural points, its exact-order modular marking is unique up to
inner conjugacy.  Neither degree-nine action has such uniqueness.

Equation `(10)` does not choose a marking from geometric data.  It only
classifies markings once a quotient map from `(1)` has been supplied.

## 4. The half-face detects the quartic translation kernel

The mixed-word cycle types on the natural sets are

```text
S4 generating row:          (4);
AGL_2(F_3) generating row:  (1,8);
PGL_2(F_8) generating rows: (1,1,7) or (9).                (11)
```

In the affine four-point chart, every generating `w` has

```text
w^2 in V4\{1},                                                (12)
```

where `V4=F_2^2` is the translation subgroup.  In the projective four-point
chart the same square is a double transposition, hence lies in the canonical
normal matching `V4`.  This is exactly the half-face appearing in THM-2775:
the `(2,3,4)` face square is the nonzero quartic matching translation.

For every generating pair in `AGL_2(F_3)`, `w^4` has cycle type

```text
(1,2,2,2,2).                                                (13)
```

It is not a translation: every nonzero element of the ternary translation
subgroup has cycle type `(3,3,3)`.  The projective counterfeit has odd word
order, so it has no nontrivial half-face at all.  Hence among the entire
prime-power degree handshake, the condition

```text
some proper power of the first mixed word is a nonzero affine translation (14)
```

selects the `S4=V4 semidirect S3` clutch.

This remains finite-group anatomy.  A graph quartic or Keller cover must
still realize the marked quotient and place `(12)` as geometric monodromy or
inertia.  THM-2775's affine-parity boundary shows why the abstract half-face
does not prove that placement.

## 5. Two sharp counterfeits inside the counterfeit

The order-eight test does not recover the ternary affine origin.  Among the
`1296` pairs in `AGL_2(F_3)` with `ord(w)=8`, exactly

```text
432 generate a point stabilizer GL_2(F_3) of order 48,
864 generate the full affine group of order 432.              (15)
```

Both rows give the same word cycle type `(1,8)`.  The point-stabilizer row
has natural orbits `1+8`; the full row is transitive.  Thus even the complete
permutation of the mixed word cannot distinguish a split linear complement
from the full affine lift.  One must retain the joint pair action, an affine
origin class, or transitivity.  This is the degree-nine analogue of
THM-2595's warning that local free-factor gauges do not determine the global
affine cocycle.

There is a different split on the projective side.  An order-seven element
of `PGL_2(F_8)` lies in the split torus `F_8^*`: its natural cycle type is
`1+1+7`.  An order-nine element lies in the nonsplit norm-one torus and is a
`9`-cycle.  The two rows have equal size `1512`.  Therefore THM-3088's
residual scalar fibre

```text
F_8^*=C7                                                     (16)
```

is visible in one half of the modular markings but is not selected by the
degree-nine action.  Residual `C7` structure alone does not repair the
counterfeit or choose a Farey/Hurwitz marking.

## 6. Exact proof and evidence

The companion pins the LF hash of THM-3090's finite-field action constructor.
That constructor enumerates every matrix in `GL_2(F_2)`, `GL_2(F_3)`, and
`GL_2(F_8)`, every affine shift, and every scalar-normalized projective
matrix, producing faithful groups of orders `24,432,24,504`.

For every pair satisfying `(4)`, the companion:

1. computes `ord(sc)` and the exact natural cycle type;
2. closes `<s,c>` by breadth-first multiplication and records its order and
   point-orbit partition;
3. checks the full tables `(5)--(7)`;
4. enumerates each centre and proves the marking counts `(10)`;
5. constructs the binary and ternary translation subgroups and checks every
   half-power assertion `(12)--(13)`; and
6. retains the split-complement hostile `(15)`.

Run

```text
python 04-computation/modular_mixed_word_fingerprint_thm3092.py
python -O 04-computation/modular_mixed_word_fingerprint_thm3092.py
```

Both modes byte-match the stored transcript after LF normalization.  The
program uses explicit `require` gates and no truth-bearing Python assertions.

```text
PROVED IN THE CANDIDATE:
  the complete exact-order marked pair tables (5)--(7);
  the disjoint generating mixed-word spectra (8);
  inner marking counts 1,2,6;
  the unique S4 translation half-face;
  the AGL split-complement and PGL split/nonsplit hostiles.

NOT PROVED:
  independent hostile audit or promotion;
  a canonical modular marking on a quartic or physical carrier;
  a graph-quartic, resolvent, Jelonek, or Keller realization;
  a literal binary/ternary tree identification;
  a tournament, LRC(14), JC(2), or DC(2) consequence.          (17)
```

QED (candidate).
