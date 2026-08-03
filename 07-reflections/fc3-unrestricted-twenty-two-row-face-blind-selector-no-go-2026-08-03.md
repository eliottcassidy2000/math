# Unrestricted twenty-two-row face-blind selector no-go on two FC3 response faces

**Status: PROMOTED AS
[THM-3275](../01-canon/theorems/THM-3275-unrestricted-twenty-two-row-face-blind-selector-obstruction.md),
VERIFIED-EXACT, AND INDEPENDENTLY HOSTILE-AUDITED.** The assertion-free
companion reconstructs the two promoted bank-I2 response banks from pinned
THM-3249 and THM-3244 sources. An independent direct product-Gamma engine
recomputed all 4,558 response vectors without using the stored cover tables.
Normal and optimized replays byte-match the frozen output. This note does not
claim FC(3), SFC(3), another support/bank face, or a Gaussian-moment
functional.

## 1. Question inherited from THM-3266

THM-3266 proves that none of the 24 fixed row pairs covering both the
support-(1,2), bank-I2 and support-(1,3), bank-I2 faces admits a face-blind
count-vector selector. Its universal singleton witness `(5)` does not settle
the larger question because row 12 is available on both faces there.

The unrestricted question is exact. For a shared padded count vector `n`, let

```text
A_12(n) = rows among all 22 with a strict Q_12-directed one-pole ascent,
A_13(n) = rows among all 22 with a strict Q_13-directed one-pole ascent.
```

A deterministic selector depending only on `n` can serve both faces exactly
when

```text
A_12(n) intersect A_13(n) is nonempty                     (1)
```

at every vector where both faces have a nonreset decision. This reduces the
problem to an exact fibre-intersection census, without imposing a threshold
tree, fixed pair, linear blend, or lookup representation.

## 2. Exact obstruction locus

All 239 padded vectors of the small face occur in the 4,319-vector full face.
The small reset has no small-face decision obligation, so there are 238 shared
joint-nonreset inputs. Both facewise availability sets are nonempty at every
one. Their intersection is empty at exactly two states:

```text
n=(3,4,5):
  A_12={2,5,8,9,11,12,14,16,18,22},
  A_13={3,4,6,7,10,13,17,19,20,21};

n=(1,3,4,5):
  A_12={2,5,9,11,12,14,16,18,22},
  A_13={3,4,6,7,10,13,17,19,20,21}.                  (2)
```

Thus no function from the untagged padded count vector to any of the 22 rows
is lawful on both faces. This is the unrestricted statement deliberately left
open by THM-3266.

The complete intersection-size histogram on the 238 inputs is

```text
(size,number of vectors)=
(0,2),(1,4),(2,2),(7,1),(9,8),(10,8),(11,11),(12,11),
(13,15),(14,10),(15,7),(16,11),(17,10),(18,16),
(19,13),(20,29),(21,52),(22,28).                      (3)
```

So the obstruction is sharply localized: 236 of 238 shared decision vectors
admit a face-blind row, while precisely two do not.

## 3. All natural minimalities

The two hostile states form one short boundary chain:

```text
(3,4,5) --insert 1 toward Q_12--> (1,3,4,5).           (4)
```

They exhaust the conflict locus, and therefore give the following exact
classifications.

- `(3,4,5)` is the unique minimum-pole-cardinality witness and the unique
  multiset-inclusion-minimal witness.
- `(1,3,4,5)` is the unique conflict closest to the small reset, at count-vector
  `l1` distance 2; `(3,4,5)` has distance 3.
- Their full-reset distances are respectively 5 and 4.

The no-go is therefore not a diffuse incompatibility across the two boxes. It
is one Q-directed edge on which the entire 22-row atlas changes sides.

This also resolves the six-state intersection in THM-3266. The common-row
sets on those six states are

```text
(5)             -> {12},
(4,5)           -> {5},
(1,3,5)         -> {5},
(1,4,5)         -> {5},
(3,4,5)         -> empty,
(1,3,4,5)       -> empty.                             (5)
```

Rows 5 and 12 repair four states after leaving any fixed common pair. The
remaining adjacent two-state seam cannot be repaired even after allowing all
22 rows. The singleton `(5)` is a positive control, not the unrestricted
hostile witness.

## 4. Exact sidecar cost

Zero origin bits are impossible by either state in (2). One persistent bit
encoding which of the two faces supplied the state is sufficient: choose any
available row on the indicated face.

More sharply, the bit need influence the selector on exactly two raw vectors.
On each of the other 236 shared nonreset vectors, choose the least member of
the common availability set. On the two conflicts, choose the least row on the
tagged face. On full-only inputs, choose the least full-face row. This policy
is exact on

```text
238 small-face decisions + 4,318 full-face decisions = 4,556 cases. (6)
```

Its common face-blind row census on the 236 compatible vectors is

```text
row 1:184, row 2:25, row 3:21, row 5:5, row 12:1.      (7)
```

Every lawful selector must distinguish the face at both states in (2), since
the two availability sets are disjoint there. Hence the minimum static
face-origin width is one bit, and the minimum raw-vector dependency locus for
that bit is exactly the two-state seam.

## 5. Mechanism and reframe

The missing coordinate is not globally expensive. Forgetting the face leaves
almost the entire row atlas mergeable, then creates a rank-zero compatibility
hole on one directed edge. This suggests treating face origin as a sparse
transition defect rather than attaching an undifferentiated tag to the whole
state box.

For further face extensions, the useful invariant is therefore not only the
number of bits in a global face label. It is the **origin-dependence locus**

```text
H(f,g)={n : A_f(n) intersect A_g(n) is empty}.         (8)
```

Here `H(12,13)` has two vertices and one Q-directed edge. With three or more
faces, the hypergraph of empty simultaneous intersections may reveal whether
one categorical face label is wasteful, whether local binary sidecars suffice,
or whether a genuine higher transition cocycle appears. This is a concrete
next computation and is better typed than calling the present two-face split
a tournament: availability has ties, multiple lawful rows, and empty
intersections, so no intrinsic total orientation is present.

## 6. Scope and reproduction

The universe is exactly the two promoted complete bank-I2 faces and the 22
lawful THM-3238 response rows. “Face-blind” means an arbitrary deterministic
memoryless function of the padded count vector. The result is stronger than a
fixed-pair or threshold-tree no-go, but it does not exclude a controller whose
history itself retains origin, and it says nothing about other maintained
faces.

Availability is a local strict response ascent, not a global positive moment
functional. No FC(3), SFC(3), Gaussian Moment Conjecture, or Factorial
Conjecture conclusion follows.

Raw artifacts:

- `04-computation/fc3_unrestricted_twenty_two_row_face_blind_selector_no_go_20260803.py`
- `05-knowledge/results/fc3_unrestricted_twenty_two_row_face_blind_selector_no_go_20260803.out`

Replay both

```text
python3 04-computation/fc3_unrestricted_twenty_two_row_face_blind_selector_no_go_20260803.py
python3 -O 04-computation/fc3_unrestricted_twenty_two_row_face_blind_selector_no_go_20260803.py
```

and compare LF-normalized bytes with the frozen transcript. The companion has
no `assert` node or floating literal, pins six direct artifacts, reconstructs
both response banks, freezes the full availability-atlas digest, checks
positive and hostile controls, and verifies the localized one-bit selector on
all 4,556 tagged decision cases.
