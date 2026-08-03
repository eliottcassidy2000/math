---
id: THM-3275
title: "Unrestricted twenty-two-row face-blind selector obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Across the two
  promoted bank-I2 FC3 response faces, exactly two of the 238 shared
  joint-nonreset padded count vectors have disjoint availability sets even
  after all 22 lawful response rows are allowed. Hence no deterministic
  memoryless selector of the untagged count vector can serve both faces. One
  face-origin bit is necessary and sufficient, and its minimum dependency
  locus is exactly those two adjacent vectors. This is a local response-bank
  theorem, not FC(3), SFC(3), or a positive-functional theorem.
source: root/2026-08-03
audit: >
  The assertion-free exact companion pins THM-3244, THM-3249 and THM-3266,
  reconstructs both complete response banks, enumerates all 238 shared
  decision fibres and all 22 rows, freezes the full availability atlas, and
  verifies the localized tagged policy on 4,556 decisions. An independent
  hostile audit bypassed the stored cover tables and rebuilt all 4,558
  response vectors directly from the product-Gamma coefficient formula on
  four workers. It reproduced the universes, histogram, two conflicts,
  directed edge, minimalities, policy census and dependency locus. Normal,
  optimized and stored outputs agree byte-for-byte.
depends_on:
  - THM-3244-unique-reset-exposure-deletion-graph-nonmorse-boundary
  - THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge
  - THM-3266-all-common-two-face-row-pairs-require-one-origin-bit
related:
  - THM-3258-depth-two-affine-farkas-clutch-and-complete-reset-distance-gauge-no-go
script: 04-computation/fc3_unrestricted_twenty_two_row_face_blind_selector_no_go_20260803.py
output: 05-knowledge/results/fc3_unrestricted_twenty_two_row_face_blind_selector_no_go_20260803.out
script_sha256: 2ad9eacf8e893f881b5616672d5fde872a80612b780372dd2318d57f75b6ea30
output_sha256: 664b7677e89873dda26cd004878c134e621086095d4155d307b52213471baf63
hash_basis: LF-normalized bytes
---

# THM-3275 -- unrestricted twenty-two-row face-blind selector obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. The exact selector question

Retain the complete bank-I2 response faces

```text
F_12 = support-(1,2), with 239 states,
F_13 = support-(1,3), with 4,319 states,
```

from
[THM-3249](THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge.md)
and
[THM-3244](THM-3244-unique-reset-exposure-deletion-graph-nonmorse-boundary.md).
All 239 padded count vectors of `F_12` occur in `F_13`. The small reset has
no small-face decision obligation, leaving 238 shared vectors where both
faces must choose a strict one-pole ascent toward their respective reset.

For such a vector `n`, let

```text
A_12(n) = lawful rows among all 22 on F_12,
A_13(n) = lawful rows among all 22 on F_13.              (1)
```

A deterministic memoryless selector depending only on the untagged padded
count vector serves both faces if and only if

```text
A_12(n) intersect A_13(n) != empty                       (2)
```

for every one of the 238 shared decision vectors. This formulation allows an
arbitrary lookup function into all 22 rows; it imposes no fixed pair,
threshold tree, linear blend or reset-distance gauge.

## 2. Complete obstruction locus

Both availability sets in `(1)` are nonempty at every shared decision
vector. Their intersection is empty at exactly two:

```text
n=(3,4,5):
  A_12={2,5,8,9,11,12,14,16,18,22},
  A_13={3,4,6,7,10,13,17,19,20,21};

n=(1,3,4,5):
  A_12={2,5,9,11,12,14,16,18,22},
  A_13={3,4,6,7,10,13,17,19,20,21}.                    (3)
```

Consequently `(2)` fails and no untagged memoryless selector into all 22
rows exists.

The exact common-availability histogram is

```text
(intersection size, number of vectors)=
(0,2),(1,4),(2,2),(7,1),(9,8),(10,8),(11,11),(12,11),
(13,15),(14,10),(15,7),(16,11),(17,10),(18,16),
(19,13),(20,29),(21,52),(22,28).                        (4)
```

Thus 236 of 238 fibres are compatible and the obstruction is localized.

## 3. Boundary geometry and minimality

The two conflicts form exactly one directed edge in the small-face reset
order:

```text
(3,4,5) --insert 1 toward Q_12--> (1,3,4,5).            (5)
```

Among the complete conflict locus:

- `(3,4,5)` is the unique minimum-cardinality and unique
  multiset-inclusion-minimal conflict;
- `(1,3,4,5)` is uniquely closest to the small reset, at padded-count
  `l1` distance two; the first conflict has distance three;
- their distances from the full reset are respectively five and four.

This also repairs the interpretation of the six fixed-pair hostiles in
[THM-3266](THM-3266-all-common-two-face-row-pairs-require-one-origin-bit.md).
Their all-row common availability is

```text
(5)             -> {12},
(4,5)           -> {5},
(1,3,5)         -> {5},
(1,4,5)         -> {5},
(3,4,5)         -> empty,
(1,3,4,5)       -> empty.                              (6)
```

Rows 5 and 12 repair four fixed-pair conflicts, but no row repairs the final
two-state seam. In particular, singleton `(5)` is a positive control for the
unrestricted question, not its hostile witness.

## 4. Exact sidecar cost

Zero origin bits are impossible by either vector in `(3)`. One persistent
bit identifying which of the two faces supplied the state is sufficient:

1. on each of the 236 compatible shared vectors, choose the least common
   available row;
2. on either conflict, choose the least row available on the tagged face;
3. on a full-only vector, choose the least full-face available row.

This policy is lawful on exactly

```text
238 small-face decisions + 4,318 full-face decisions = 4,556.            (7)
```

On the 236 compatible shared vectors, the least-common-row census is

```text
row 1:184, row 2:25, row 3:21, row 5:5, row 12:1.        (8)
```

Every lawful tagged policy must depend on the origin at both conflicts,
because the two row sets there are disjoint. The minimum static origin width
is therefore one bit, and the minimum raw-vector dependency locus of that bit
is exactly the two vertices in `(5)`.

## 5. Mechanism and next object

The relevant invariant is not merely the global width of a face label. For
two response faces define the **origin-dependence locus**

```text
H(f,g)={n : A_f(n) intersect A_g(n) is empty}.           (9)
```

Here `H(12,13)` is a two-vertex, one-edge seam. With three or more faces, the
hypergraph of empty simultaneous intersections is the natural next object:
it can distinguish a sparse local binary sidecar from a genuinely
higher-valued transition defect. This structure is not a tournament; the
availability relation has ties, multiple lawful rows and empty intersections.

## 6. Scope and reproduction

The theorem concerns exactly the two named promoted bank-I2 faces and the 22
lawful response rows. Availability means a local strict response ascent, not
a positive Gaussian-moment functional. The theorem excludes arbitrary
deterministic memoryless functions of the padded counts, but not controllers
whose history already retains face origin. It proves no other-face closure,
no `FC(3)` or `SFC(3)`, and no Factorial Conjecture conclusion.

Run

```text
python3 04-computation/fc3_unrestricted_twenty_two_row_face_blind_selector_no_go_20260803.py
python3 -O 04-computation/fc3_unrestricted_twenty_two_row_face_blind_selector_no_go_20260803.py
```

and compare LF-normalized bytes with the declared output. The companion has
no assertion node or floating literal, pins six direct artifacts, freezes the
complete availability-atlas digest, and checks positive and hostile controls.

QED.
