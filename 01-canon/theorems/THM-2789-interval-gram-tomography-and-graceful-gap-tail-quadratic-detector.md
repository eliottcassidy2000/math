---
id: THM-2789
title: "Interval Gram tomography and graceful gap-tail quadratic detector"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For intervals on a finite line, every nonempty higher
  intersection equals the minimum pairwise Gram entry.  Boolean Mobius
  inversion therefore reconstructs the complete incidence-column
  multiset from the pairwise Gram matrix and ambient length.  Applied to a
  graceful gap matrix, this detects and counts every suffix-gap candidate
  from quadratic-size data.  It does not recover gap order or the marked
  attachment rank, and the reconstruction fails for arbitrary set systems
  and circular arcs without a retained cut.
source: root/interval-gram-gap-tail-2026-07-28
depends_on:
  - THM-2787-signed-path-sum-weyl-orbit-and-gap-tail-leaf-insertion
related:
  - THM-2383-polarized-complete-subcube-gram-tomography
  - THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge
  - THM-2786-binary-golomb-universal-edge-difference-separation-and-graceful-compression-boundary
script: 04-computation/interval_gram_gap_tail_tomography_thm2789.py
output: 05-knowledge/results/interval_gram_gap_tail_tomography_thm2789.out
script_sha256: 86694905982ceadcbfc58c2d1c7ceb6f3f56dad133901ec9d1722a7e2da11a58
output_sha256: dd9a8522ca3fba1fb38f94ce03851d1818351ab980981adf722e3fad98ac59c4
hash_basis: LF-normalized bytes
---

# THM-2789 -- interval Gram tomography and graceful gap-tail quadratic detector

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2787 turns one-leaf graceful extension into a suffix-column problem in
a consecutive-ones matrix.  A priori that still appears to require the
whole ordered matrix.  On a line, however, pairwise overlaps retain every
higher intersection.  Consequently a quadratic-size Gram matrix recovers
the complete multiset of columns exactly.

This is a strong compression and a sharp nonclosure.  It detects every
suffix column, but it forgets where that column occurs and which side of the
gap contains a marked parent.

## 1. Quantitative Helly identity on a finite line

Let `X={1,...,N}` with its usual order, and let

```text
I_1,...,I_m subset X
```

be intervals, possibly empty.  Write `C` for their incidence matrix and

```text
G=C C^T,                 G_ij=|I_i intersect I_j|.          (1)
```

For every nonempty `S subset {1,...,m}`, put

```text
h(S)=|intersection_(i in S) I_i|.                         (2)
```

Then

```text
h(S)=min_(i,j in S) G_ij.                                 (3)
```

To prove `(3)`, write a nonempty interval as `[l_i,r_i] intersect X`,
using `l_i>r_i` for the empty interval.  Choose `a in S` with maximal left
endpoint and `b in S` with minimal right endpoint.  Then

```text
intersection_(i in S) I_i = I_a intersect I_b.            (4)
```

If the total intersection is empty, this chosen pair is empty as well.
Every pair contains the total intersection, so `(4)` is the smallest
pairwise overlap.  This proves `(3)`, including singleton `S` through the
diagonal entry `G_ii`.

Equation `(3)` is stronger than ordinary one-dimensional Helly: it recovers
the exact size, not only nonemptiness, of every higher intersection.

## 2. Gram reconstruction of the column multiset

For each `x in X`, define its incidence column

```text
A_x={i:x in I_i} subset [m],
```

and let

```text
n(A)=#{x in X:A_x=A}.                                    (5)
```

Set `h(empty)=N`.  Counting columns containing a fixed set gives

```text
h(S)=sum_(A superset S) n(A).                             (6)
```

Boolean Möbius inversion and `(3)` now give the explicit reconstruction

```text
n(A)=sum_(B superset A) (-1)^(|B|-|A|) h(B),              (7)

h(empty)=N,
h(B)=min_(i,j in B) G_ij                  (B nonempty).   (8)
```

Therefore `(G,N)` determines the complete multiset

```text
{{A_x:x in X}}                                           (9)
```

with multiplicity.  In particular, two interval families on equally long
finite lines have the same Gram matrix if and only if they have the same
incidence-column multiset.

The phrase *quadratic detector* below refers to the size of the sufficient
statistic `G`, which has `m^2` entries.  Formula `(7)` is a Boolean Möbius
sum; no claim of an `O(m^2)` reconstruction algorithm is made.

## 3. Graceful gap-tail application

Fix a graceful labeling of an `m`-edge tree as in THM-2787.  Restrict its
gap matrix to the `m` internal cuts

```text
X={1,...,m}.
```

For the unique edge of difference `d`, put

```text
I_d={t in X:C_f(d,t)=1}.                                 (10)
```

THM-2787 gives that `I_d` is an interval of cardinality `d`.  Form the
gap Gram matrix

```text
G_de=sum_(t=1)^m C_f(d,t) C_f(e,t)
    =|I_d intersect I_e|.                                (11)
```

For `1<=k<=m`, let

```text
T_k={k,k+1,...,m} subset [m].                            (12)
```

The number of internal gaps whose column is the suffix `T_k` is exactly

```text
n(T_k)
 =sum_(B superset T_k)(-1)^(|B|-|T_k|)
    min_(i,j in B)G_ij.                                  (13)
```

Thus `(11)` detects and counts every internal Ferrers/suffix column required
by THM-2787.  The two exterior gaps `t=0,m+1` are the two all-zero columns.
There are no other all-zero columns: for `1<=t<=m`, the threshold set is a
nonempty proper vertex set, and connectedness forces its tree cut to be
nonempty.

Consequently, the existence and multiplicity of every suffix-cut candidate
is already present in quadratic-size overlap data.  This is exactly the
information forgotten by treating the edge intervals only through their
individual lengths `1,...,m`.

It is not the whole marked extension condition.  THM-2787 additionally
requires, at the actual gap position `t`,

```text
|t-iota_t(f(v))|=k.                                      (14)
```

The Gram matrix recovers the multiset of columns but not their linear order,
the position `t`, the marked vertex label `f(v)`, or the side of the
threshold cut containing `v`.  Those are the exact sidecars still needed
for a Graceful Tree induction.

## 4. Three sharp failure boundaries

### 4.1 Arbitrary set systems

On four points, take column multisets

```text
{12},{13},{23},empty
```

and

```text
{123},{1},{2},{3}.                                      (15)
```

Both have Gram matrix

```text
G=((2,1,1),(1,2,1),(1,1,2)),                            (16)
```

but their triple intersections have sizes zero and one.  Thus `(3)` and
the Gram reconstruction fail without interval convexity.

### 4.2 Ordered gaps

On the three-point line, take one interval of each length one, two, and
three.  The start vectors

```text
(0,0,0)       and       (1,0,0)                         (17)
```

have the same Gram matrix and the same column multiset

```text
{123},{23},{3},                                         (18)
```

but the first two columns occur in opposite order.  These two ordered
families are not related by reflecting the line.  Hence even inside the
interval category, Gram data cannot identify the gap position in `(14)`.

### 4.3 Circular arcs

On the six-cycle, let

```text
I_1={0,1,2},   I_2={2,3,4},   I_3={4,5,0}.              (19)
```

Every pair intersects in one point while the triple intersection is empty.
The arcs cover the entire circle, so there is no common omitted point at
which to cut them into line intervals.  A circular-arc or LRC transplant
therefore requires a retained anchor/cut; cyclic order alone is
insufficient.

## 5. Exact finite controls

The companion enumerates every family consisting of one interval of each
length `1,...,m` on an `m`-point line through `m=8`.  There are `m!` such
families.  In every case it verifies `(3)`, `(7)`, nonnegative reconstructed
multiplicities, and constancy of the reconstructed multiset on each Gram
fibre.

The numbers having at least one nonempty suffix column are

```text
m:             1   2   3   4    5    6     7      8
all families:  1   2   6  24  120  720  5040  40320
suffix exists: 1   2   6  24  120  718  4990  39744.     (20)
```

These are interval-system controls, not counts of realizable graceful tree
gap matrices.  They show both why the suffix condition is common in the
ambient consecutive-ones space and why it is not automatic from row lengths
alone.

Run

```bash
python 04-computation/interval_gram_gap_tail_tomography_thm2789.py
python -O 04-computation/interval_gram_gap_tail_tomography_thm2789.py
```

The script uses explicit exceptions and no truth-bearing assertions.

```text
PROVED HERE (candidate):
  quantitative interval Helly identity;
  complete incidence-column reconstruction from (G,N);
  exact graceful suffix-column count from the gap Gram matrix;
  the two exterior-zero-column boundary;
  arbitrary-set, ordered-column, and circular-arc hostiles;
  exact interval-family controls through m=8.

NOT PROVED:
  an O(m^2)-time reconstruction algorithm;
  recovery of the ordered gap matrix from G;
  the marked rank condition;
  the marked gap-tail extension conjecture;
  the Graceful Tree Conjecture;
  a circular/LRC application without a retained cut;
  a Keller, JC(2), DC(2), knot, or tournament consequence.             (21)
```

QED, conditional only on candidate status promotion after exact replay and
independent hostile audit.
