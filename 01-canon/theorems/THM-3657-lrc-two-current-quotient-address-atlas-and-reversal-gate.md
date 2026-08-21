---
id: THM-3657
title: "LRC two-current quotient address atlas and reversal gate"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.
  In THM-3593's rank-two state/constant correction quotient, the 169
  two-current addresses split into 37 kernel rows, one 124-address reversal-
  fixed projective line, and eight singleton lines in four reversal pairs.
  The 17 source-current character rows occupy distinct lines disjoint from
  every two-current line.  Hence any linear cancellation of a source defect
  by two-current rows must use one of the eight exceptional addresses.  This
  is a static finite-field gate, not chronology, current, physical entry,
  characteristic-zero transport, row exclusion, or LRC(14).
source: kps-s191 / THM-3593 address-atlas continuation, 2026-08-21
depends_on:
  - THM-3593-lrc-common-a4-anova-graph-flag
related:
  - THM-3534-r5-middle-response-relative-cospan-and-twisted-h1-collapse
  - THM-3647-lrc-single-reversal-paired-branch-spectral-projector
  - THM-3654-lrc-fixed-branch-rigidity-eigencriterion
script: 04-computation/lrc_two_current_quotient_address_reversal_gate_thm3657.py
output: 05-knowledge/results/lrc_two_current_quotient_address_reversal_gate_thm3657.out
script_sha256: f0323550a039bd3c59bc3367a9a48503ff10db2b642e99e21d9499e984492ccd
output_sha256: fdebbbba161149edec8c86b9f382ede76153e85aa01da635f0900075fe751628
semantic_sha256: 792059590ad928cdb096763c6fd429b8f02b991593efef68964e26c57f7b7059
hash_basis: raw LF bytes
---

# THM-3657 -- the LRC correction quotient has eight exceptional addresses

**PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
AUDIT.**  The theorem is an address-level classification inside the static
rank-two quotient of THM-3593.  Its cancellation corollary permits arbitrary
field coefficients, but it does not assert that such coefficients or address
combinations are physically realizable.

## 1. The correction coordinate

Work over

```text
F=F_p,                 p=755373809845391722745761.     (1)
```

THM-3593 writes the common raw four-plane as a relation-led graph

```text
B=graph(U:R4 -> E2),              ker U=K2,            (2)
```

where `E2=(P_0+P_s)A` is the state/constant correction plane.  For any raw
row `a` in the common plane, put

```text
e(a)=(P_0+P_s)a in E2.                                (3)
```

Thus `e(a)=U(P_r a)`, and in particular

```text
e(a)=0  iff  P_r a lies in K2.                         (4)
```

The companion reconstructs both presentations of the common plane from the
pinned parents: the 17 source-current rows `a_q`, indexed by the inherited
character bank, and the 169 two-current rows `a_(r0,r1)`, indexed by
`F_13^2`.  It applies `(3)` row by row.  Both families span the same
two-plane `E2`; its canonical RREF digest is

```text
2d82d5adbb3c5a71035ab622e36e51e2a6a701a941a9f6f31ca308b8de13d17c. (5)
```

All coordinates below use that RREF basis.  Numerical slopes depend on this
choice, while zero rows, projective multiplicities, disjointness, reversal
conjugacy, and the cancellation gate do not.

## 2. The complete two-current atlas

The zero set is exactly

```text
Z = F13 x {6}
    union (F13\{12}) x {0}
    union (F13\{0})  x {12},             |Z|=37.       (6)
```

This is an equality, not a one-way vanishing test.  By `(4)`, these and only
these address rows have relation component in `K2`.

Put

```text
X = {(12,0),(0,12),(12,1),(0,11),
     (6,5),(6,7),(3,9),(9,3)}.                         (7)
```

Among the remaining 132 nonzero rows, all 124 labels outside `X` occupy one
projective line, while the eight labels in `X` occupy eight pairwise distinct
singleton lines.  Hence the complete nonzero class-size multiset is

```text
(124,1,1,1,1,1,1,1,1).                               (8)
```

The six boundary/middle-neighbor labels in `(7)` recover THM-3534's endpoint
residual, while `(3,9),(9,3)` are the crossed representatives of its marked
two-chamber location `Loc_(3,9)`.  This provenance is explanatory only: the
present theorem reconstructs and checks all 169 rows directly, and it does
not identify the two typed response spaces.

## 3. Reversal acts by an exact quotient involution

Let simultaneous point reversal act on address pairs by

```text
j(r0,r1)=(12-r0,12-r1).                               (9)
```

In the RREF coordinates of `(5)`, define

```text
b     =371578917865089240854253,
b^(-1)=555904782330327598552358,
J_E   =((0,b),(b^(-1),0)).                             (10)
```

Exact row-by-row reconstruction proves

```text
e(a_(j(r0,r1)))=e(a_(r0,r1)) J_E     for all 169 labels,
J_E^2=I.                                               (11)
```

The two eigendirections are

```text
L_+=F(1,b),
L_-=F(1,-b).                                           (12)
```

The 124-address class in `(8)` is exactly `L_+\{0}` at the projective level.
No two-current row realizes `L_-`.  The exceptional set is the union of the
four transposed orbits

```text
(12,0)<->(0,12),   (12,1)<->(0,11),
(6,5) <->(6,7),    (3,9) <->(9,3).                    (13)
```

Thus the apparent `37+124+8` count is governed by the same branch reversal
that fixes digit `6`; it is not merely a histogram coincidence.  THM-3647
finds the reversal's spectral collapse in the six-point address algebra,
whereas `(11)` is its newly reconstructed action on the independent
two-dimensional correction quotient.  No identification of those two typed
spaces is being made.

## 4. The multi-stranger exceptional-address gate

Here *multi-stranger* means only an arbitrary number of two-current address
rows in one linear relation.  It is not THM-518's unrelated large-speed
"stranger" terminology or its Weyl-decoupling statement.

Every one of the 17 source-current corrections is nonzero, and their 17
projective lines are pairwise distinct.  More strongly,

```text
{F e(a_q):q in the source bank}
  intersect
{F e(a_r):r in F13^2\Z} = empty.                       (14)
```

Let `q` be any source label and suppose that a linear cancellation in `E2`
has the form

```text
e(a_q)+sum_i lambda_i e(a_(r_i))=0.                   (15)
```

If every active label `r_i` lay outside `X`, then `(6)--(8)` would put every
summand on `L_+` or at zero.  Their sum would lie on `L_+`, forcing
`e(a_q)` onto `L_+`, contrary to `(14)`.  Therefore

```text
some i has lambda_i!=0 and r_i in X.                  (16)
```

This is the promised multi-stranger gate: regardless of how many generic
two-current rows are allowed, at least one of eight exceptional addresses is
necessary to absorb one source defect.  It converts a 169-address search into
an eight-address hitting obligation.  It does **not** prove that an
exceptional address is sufficient, nor that chronology exposes any member of
`X`, nor that arbitrary coefficients in `(15)` are legal.

## 5. Verification and strict boundary

Reproduce with

```bash
python3 -B 04-computation/lrc_two_current_quotient_address_reversal_gate_thm3657.py
python3 -B -O 04-computation/lrc_two_current_quotient_address_reversal_gate_thm3657.py
```

The assertion-free companion source-pins THM-3593, reconstructs both parent
tensors, verifies the common correction basis, classifies every label,
checks the exact covariance `(11)`, both eigendirections and all four
exceptional orbits, and proves `(14)` by two-dimensional rank tests.  The
semantic digest is

```text
792059590ad928cdb096763c6fd429b8f02b991593efef68964e26c57f7b7059. (17)
```

Everything here is over the pinned finite field and inside one static common
row space.  There is no chronology, current, physical-entry theorem,
characteristic-zero lift, admissible coefficient theorem, row exclusion, or
LRC(14) conclusion.  **QED.**
