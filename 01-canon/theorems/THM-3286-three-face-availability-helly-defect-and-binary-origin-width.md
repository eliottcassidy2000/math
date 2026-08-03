---
id: THM-3286
title: "Three-face availability Helly defect and binary origin width"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On the three named bank-I2 support faces (1,2), (1,3), and (2,3), the
  twenty-two THM-3249 response rows cover every nonreset physical state.
  The natural third face has a distinct pole profile.  Its pairwise
  availability intersections with either promoted face are never empty,
  but singleton (5) is the unique proper three-way Helly defect.  Together
  with the two inherited THM-3275 conflicts, this gives three and only three
  triple obstructions.  Their local branch width and the minimum fixed
  face-origin alphabet both equal two: one bit is necessary and sufficient,
  and the third face may share either binary class.  This is a finite
  response-bank theorem, not FC(3), SFC(3), GMC, or positivity.
source: root/creative-synthesis-recover/2026-08-03
audit: >
  The primary assertion-free exact reconstruction derives all three pole
  banks and resets from the product-Gamma source, constructs the physical
  universes and reset-directed neighbours by dual methods, rebuilds all
  8,397 response vectors, and exhausts every pair and triple availability
  intersection.  The independent hostile audit bypasses the primary
  implementation and THM-3249 availability tables: it starts from THM-3238's
  product-Gamma coefficients, independently reconstructs partitions, upsets,
  physical submultisets, neighbours, and all twenty-two responses, and
  exactly matches the three face digests, pair/triple loci, binary
  colourings, and dependency locus.  It additionally proves by exhaustion
  that THM-3278's canonical pair (16,17) has no global three-face extension.
  Both implementations in normal and optimized mode byte-match their stored
  transcripts.
depends_on:
  - THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge
  - THM-3275-unrestricted-twenty-two-row-face-blind-selector-obstruction
related:
  - THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary
  - THM-2292-common-catalytic-section-and-helly-calibration-nerve
script: 04-computation/fc3_three_support_face_availability_hypergraph_scout_20260803.py
output: 05-knowledge/results/fc3_three_support_face_availability_hypergraph_scout_20260803.out
script_sha256: c31eff8a10c6a6e1ab7e3ea6759388b02cbd9591695182f2d3756e08da38c8c3
output_sha256: a3e44f6d2eb0e26386e399e7d582e6eb67512125233eb1b49d2b62c4d0869e05
independent_audit_script: 04-computation/fc3_three_support_face_availability_hypergraph_independent_audit_20260803.py
independent_audit_output: 05-knowledge/results/fc3_three_support_face_availability_hypergraph_independent_audit_20260803.out
independent_audit_script_sha256: 371d6f12fccaec8a99eea853e8a119e44ce127404a40fa0d6f0bdc5ea1723547
independent_audit_output_sha256: 9ed7094c717855c0c0a0d61d0d067748e05aecb2051ec3d7d3a3c8ded53471e5
hash_basis: LF-normalized bytes
---

# THM-3286 -- three-face availability Helly defect and binary origin width

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This finite theorem is established by an exact primary reconstruction and a
separately implemented hostile audit from the product-Gamma coefficient
source.  Both normal and optimized replays byte-match their stored outputs.

## 1. Three typed faces and their availability sets

Retain the twenty-two lawful response rows of
[THM-3249](THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge.md).
For `ab in {12,13,23}`, evaluate them on the support-`(a,b)`, bank-`I2`
product-Gamma face.  The physical pole multiset `P_ab`, reset `Q_ab`, and
complete nonempty state-bank size are

```text
face  pole multiplicities                              reset                 states
F12   ((1,4),(2,3),(3,2),(4,1),(5,1))                 (1,2,2,3,4,5)           239
F13   ((1,4),(2,3),(3,2),(4,2),(5,2),(6,1),(7,1),     (1,3,3,4,5,6,7,8)      4319
       (8,1))
F23   ((1,2),(2,4),(3,1),(4,3),(5,3),(6,1),(7,1),     (2,3,4,4,5,6,7,8)      3839
       (8,1)).                                                               (1)
```

In particular, `F23` is the natural third edge of the support triangle, but
its pole multiplicity profile and its state-bank size differ from both
promoted faces.  It is not a literal relabelled copy of either bank.

Let `D_ab=S(P_ab)\{Q_ab}` be the nonreset decision universe.  For
`n in D_ab`, define

```text
A_ab(n)={i in {1,...,22}:
          some one-pole move n -> n' decreases multiset l1 distance to Q_ab
          and r_i^ab(n')>r_i^ab(n)}.                    (2)
```

Thus availability means a strict local response ascent along a
reset-directed physical move.  It does not mean positivity of a Gaussian
moment or of a common linear functional.

The audited exact reconstruction gives

```text
A_ab(n) != empty for every n in D_ab and every ab in {12,13,23}.            (3)
```

## 2. Pair universes and obstruction loci

For faces `f,g`, put

```text
U(f,g)=D_f intersect D_g,
H(f,g)={n in U(f,g):A_f(n) intersect A_g(n)=empty}.     (4)
```

The complete pair census is

```text
|U(12,13)|=238,   H(12,13)={(3,4,5),(1,3,4,5)},
|U(12,23)|= 94,   H(12,23)=empty,
|U(13,23)|=1726,  H(13,23)=empty.                       (5)
```

The first row exactly reproduces the all-twenty-two-row obstruction of
[THM-3275](THM-3275-unrestricted-twenty-two-row-face-blind-selector-obstruction.md).
The other two rows are new facewise tests: `F23` has a common available row
with each promoted face at every joint decision state.

## 3. The unique proper three-way Helly defect

The triple decision universe

```text
U(12,13,23)=D_12 intersect D_13 intersect D_23          (6)
```

has exactly `94` states.  Its total availability intersection is empty at
exactly

```text
(5),       (3,4,5),       (1,3,4,5).                   (7)
```

At singleton `(5)` the three availability sets are

```text
A_12={2,8,9,11,12,14,16,18,22},
A_13={3,4,5,6,7,10,12,13,17,19,20,21},
A_23={3,4,5,6,7,8,10,13,16,17,18,19,20,21,22}.         (8)
```

Their three pairwise intersections are

```text
A_12 intersect A_13={12},
A_12 intersect A_23={8,16,18,22},
A_13 intersect A_23={3,4,5,6,7,10,13,17,19,20,21},     (9)
```

while

```text
A_12 intersect A_13 intersect A_23=empty.              (10)
```

Therefore `(5)` is a proper three-way Helly defect.  It is unique: at the
other two states in `(7)`, the first pair intersection is already empty by
THM-3275.  Their full row sets are

```text
n=(3,4,5):
 A_12={2,5,8,9,11,12,14,16,18,22},
 A_13={3,4,6,7,10,13,17,19,20,21},
 A_23={3,4,5,6,7,10,12,13,17,19,21};

n=(1,3,4,5):
 A_12={2,5,9,11,12,14,16,18,22},
 A_13={3,4,6,7,10,13,17,19,20,21},
 A_23={1,2,...,22}.                                   (11)
```

This also changes the role of a known control without contradicting it:
`(5)` remains compatible for the pair `(F12,F13)` through row `12`, exactly
as THM-3275 states, but `F23` does not make row `12` available there.

## 4. Local branch width and the global origin bit

For any raw state `n`, let `I(n)` be the set of named faces on which `n` is a
decision state.  When `|I(n)|>=2`, define `lambda(n)` as the least number of
labels needed to partition `I(n)` into blocks `B` satisfying

```text
intersection_(f in B) A_f(n) != empty.                (12)
```

The exhaustive histogram is

```text
(|I(n)|,lambda(n)) : count
(2,1):1776,       (3,1):91,       (3,2):3.             (13)
```

Thus the maximum local branch number is two; no raw state requires three
origin values.

There is also a fixed global version.  A face labelling `c` is lawful when,
for every raw state and every label value, all active faces of that label
have a common available row.  One label is impossible by any state in `(7)`.
Two labels are sufficient, and in face order `(F12,F13,F23)` the complete
list of lawful binary assignments is

```text
(0,1,0), (0,1,1), (1,0,0), (1,0,1).                   (14)
```

Indeed, `F12` and `F13` must have opposite labels because of either inherited
conflict in `(5)`.  Once they are opposite, `F23` may share either class
because both new pair obstruction loci are empty.  Choosing a least common
row within each labelled block gives a literal lawful policy on all

```text
238+4318+3838=8394                                      (15)
```

face/state decisions.  Hence the minimum fixed origin alphabet has size two,
or one bit.

The minimum raw-vector dependency locus of that bit is exactly the three
states in `(7)`.  Dependence is necessary there because the total
intersection is empty.  Away from `(7)`, choose the least row in the total
intersection and use it for every active label, so dependence can be removed.

## 5. What THM-3278 does and does not supply

[THM-3278](THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary.md)
restricts the two THM-3275 conflict fibres to a twelve-row core and identifies
their small/full split with that core's unique bipartition.  It then relates
the resulting two-face cut coordinate to a canonical paired selector,
critical-group hostiles, a finite-index cyclic lattice, and an abstract
norm-fibre transfer.

None of those statements contains `F23`, the distinct pole bank in `(1)`, the
triple universe in `(6)`, or the proper missing simplex `(8)--(10)`.  In
particular, `(5)` was not an all-row two-face conflict at all.  Conversely,
the present theorem does not extend THM-3278's core bipartition, phase
orientation, critical-group, lattice, or transfer conclusions to `F23`.
Equation `(14)` is only a static availability colouring.  It is not a new
phase character or a physical chronological walk.

There is also a sharp controller boundary.  The hostile audit exhausts all
four optimal colourings in `(14)` under the restriction that the controller
may use only THM-3278's canonical row pair `(16,17)`: none is globally lawful
on the three faces.  Thus the canonical two-row controller does not extend
from the THM-3275 pair to this three-face bank.  The one-bit sufficiency in
Section 4 instead permits an unrestricted, state-dependent choice among all
twenty-two rows within each active colour class.

The new distinction is therefore precise:

```text
availability nerve:  has one genuine missing two-simplex;
witness-colouring width: remains binary.               (16)
```

## 6. Exact verification and independent hostile audit

The exact companion performs the following without cached cover tables:

1. it pins the promoted THM-3249 and THM-3275 scripts and outputs;
2. it derives all three pole banks and resets from the pinned product-Gamma
   coefficient source;
3. it constructs each nonempty physical submultiset bank both by Cartesian
   multiplicities and by filtered combinations with replacement;
4. it constructs reset-directed neighbours independently through multisets
   and padded count vectors;
5. it rebuilds the twenty-two response coordinates on all `8,397` states;
6. it checks facewise coverage, reproduces THM-3275's full atlas digest, and
   exhausts all pair, triple, local-labelling, and fixed-labelling cases; and
7. it contains no `assert` node, floating literal, randomness, or fitted
   recurrence.  Normal and optimized output byte-match the stored transcript.

Run

```text
python3 04-computation/fc3_three_support_face_availability_hypergraph_scout_20260803.py
python3 -O 04-computation/fc3_three_support_face_availability_hypergraph_scout_20260803.py
```

and compare LF-normalized bytes with the declared output and hashes.

The independent audit does not import the primary executable or any
THM-3249 availability table.  It reconstructs the three face banks from the
THM-3238 product-Gamma coefficient source, uses independent partition/upset,
physical-submultiset, and neighbour implementations, rebuilds all twenty-two
response coordinates, and then recomputes every claimed invariant.  It
matches the three face-availability digests, every pair and triple locus, the
local-width histogram, all four binary colourings, and the exact dependency
locus.  As a further hostile test, it finds no global three-face extension
using only the row pair `(16,17)`.

Run

```text
python3 04-computation/fc3_three_support_face_availability_hypergraph_independent_audit_20260803.py
python3 -O 04-computation/fc3_three_support_face_availability_hypergraph_independent_audit_20260803.py
```

and compare LF-normalized bytes with the independent output and hashes in the
frontmatter.  The audit executable pins the primary artifacts and source
dependencies, but deliberately does not pin this theorem file, so theorem
promotion creates no verification cycle.

## 7. Scope and non-consequences

The theorem concerns exactly the three named bank-`I2` physical faces, the
twenty-two inherited lawful rows, and strict local reset-directed
availability.  The availability quotient forgets response magnitude, which
neighbour realizes an ascent, chronological compatibility, history, and
moment positivity.  It proves no `FC(3)`, `SFC(3)`, Gaussian Moment
Conjecture, positive-functional, arbitrary-bank, or arbitrary-support
conclusion.

QED.
