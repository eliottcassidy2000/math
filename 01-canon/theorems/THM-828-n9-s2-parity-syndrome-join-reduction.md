---
id: THM-828
title: The n=9 raw lower codec has exactly 58 reflection doubletons supported on a punctured rank-four mirror-defect cube
status: PROVED REDUCTION + FINITE-EXACT COMPLETE CANONICAL JOIN
source: codex-2026-07-15-S13/referee
depends_on: [THM-809, THM-811, THM-814, THM-818, THM-825]
related: [THM-801, THM-813, THM-839, HYP-6880]
verification:
  - 04-computation/mobius_cech_n8_node_atlas_export_codex_S13.py
  - 04-computation/mobius_cech_n9_exact_join_codex_S13.cpp
  - 05-knowledge/results/mobius_cech_n9_exact_join_codex_S13.out
  - 05-knowledge/results/mobius_cech_n9_exact_join_codex_S13.json
  - 05-knowledge/results/mobius_cech_n9_exact_join_witnesses_codex_S13.tsv
---

# THM-828 — the exact `n=9` join and its four-dimensional defect

Let `X_9^0` be the `2^27` size-nine staircase tilings whose apex bit is zero.
For `x in X_9^0`, let `A(x),B(x),C(x)` be its three size-eight B3 faces and
let `H_8` be THM-818's oriented lower observation

```text
H_8(y)=(nu_8(y),nu_8(kappa y),chi_8(y)).
```

Write `S2(x)` for the raw four-state reflection-pair histograms on layers
`tau=3,...,8`, followed by the one-count on the fixed `tau=9` layer.  The raw
lower address considered here is

```text
Lambda_9(x)=(H_8(Ax),H_8(Bx),H_8(Cx),chi_9(x),S2(x)).    (1)
```

Its complete partition on `X_9^0` is

```text
lines                         134,217,728
address cells                 134,217,670
singleton cells               134,217,612
doubleton cells                           58
collision excess                          58
maximum multiplicity                         2.          (2)
```

Every doubleton is exactly `{x,sigma x}`.  Thus equality of raw addresses
determines the grid-reflection orbit, although `Lambda_9` does not descend to
that quotient globally: only these 58 nontrivial reflection orbits have a
reflection-fixed address.  They are literal **false palindromes**—their
address is reflection-fixed while their tiling is not.

## 1. The raw-S2 difference code has dimension fourteen

For an equal-address pair `(u,v)`, put `D=u xor v`.  On a nonfixed layer,
each reflection pair is a letter of `G=F_2^2`, recording its low and high
bits.  Equality of the two four-state histograms implies that the XOR-sum of
the letter differences is zero.  Therefore `D` has even parity separately
on the low and high halves of each of the six layers `tau=3,...,8`.
The fixed layer has even `D` parity, and both apex bits are zero.  These are
fourteen independent equations in the 28 upper bits:

```text
12 half-layer parity equations + fixed parity + apex = 0. (3)
```

Consequently the admissible difference code has dimension `28-14=14` and
contains exactly `16,384` words.  In the repository's size-nine tile order,
one explicit basis is

```text
001010 4400000 000808 840000 020004 020400 1002000
1080000 010002 010200 100080 104000 200100 0208000.      (4)
```

The condition is existence-exact at this size.  If a layer has `s<=3`
reflection pairs and differences `delta_1,...,delta_s in G` with sum zero,
then, up to reordering, they have one of the forms

```text
s=1: (0),        s=2: (a,a),
s=3: (0,a,a) or (a,b,a+b).
```

These are respectively a fixed point, a transposition, or a three-cycle of
base letters in `G`, so one can choose the base word to make the two
histograms equal.  The fixed layer is the analogous binary statement.

This also gives a closed checksum for all ordered apex-zero raw-S2-equal
pairs before any face condition.  For layer sizes `1,1,2,2,3,3`, the sums of
squared multinomial coefficients are `4,4,28,28,256,256`; the three free
fixed positions contribute `sum_k binom(3,k)^2=20`.  Hence their number is

```text
4^2 * 28^2 * 256^2 * 20 = 16,441,671,680.               (5)
```

Equation (3) says only that some base realizes equality.  A particular base
glued from three metagraph faces must still pass literal `S2`; confusing
these two quantifiers would create many false survivors.

## 2. The three-face join factors exactly through two faces

For a size-eight face difference `d`, define the exact collision list

```text
L[d]={y : H_8(y)=H_8(y xor d)}.                          (6)
```

Every nonzero syndrome `D` in (3) has nonzero `A` and `C` restrictions.
Indeed, if `d_A=0`, the possible support is confined to the exposed top-row
cell in each constrained half-layer; the parity equations then kill those
bits one at a time, and fixed parity plus apex zero kill the remainder.
Thus `D=0`, a contradiction.  The proof for `d_C=0` is symmetric.

The `A` and `C` faces cover all 28 upper cells except the already fixed apex
and overlap in 15 cells.  Therefore an overlap-compatible pair

```text
x_A in L[d_A],       x_C in L[d_C]                       (7)
```

glues to a unique apex-zero upper base `u`; then `v=u xor D` and the `B`
face is forced.  It suffices to bucket `L[d_C]` by its 15-bit overlap, scan
one orientation `x_A < x_A xor d_A`, derive `B(u)`, and test (1) directly.
This is an exact two-face factorization of THM-818's ordered Cech join, not a
heuristic pruning rule.

The exact size-eight node atlas has `2^21` unsigned 16-bit entries and SHA256

```text
30debad3387a4ea0ef51108ea132115efda2ac2fcdfcc2c5c1d4d23155095835. (8)
```

It is deterministically exported from the certified 6,880-class classifier;
the atlas is scratch data, while its exporter and checksum are stored.
Multiple executions—including two C++ builds and an independently written
Python join—gave the same census.

## 3. Exact collision genealogy

The complete nonzero-syndrome genealogy is

| gate | supported `D` sectors | canonical candidate pairs |
|---|---:|---:|
| both `L[d_A]` and `L[d_C]` nonempty | 319 | — |
| literal `A/C` overlap match | 45 | 9,540 |
| forced `B` face satisfies `H_8` | 16 | 636 |
| upper colour agrees | 16 | 636 |
| raw `S2` through `tau=3,4` | 16 | 636 |
| through `tau=5` | — | 88 |
| through `tau=6` | — | 82 |
| through `tau=7` | — | 68 |
| through `tau=8` | 11 | 58 |
| fixed `tau=9` | 11 | 58 |

The underlying `R_8` audit is also reproduced: `5,997,416` relation rows,
`3,900,264` off-diagonal rows, `249,149` nonempty difference lists, and
maximum `|L[d]|=5,360`.  The largest `A/C` lists encountered have size
`4,104`; the largest overlap bucket has size 16.  The sorted pair list has
standard FNV-1a-64 digest `0x53b4b074be8ae851`.

The 116 endpoints in the 58 final pairs are distinct.  Completeness of the
join and endpoint disjointness prove (2): a fibre of size at least three
would contribute intersecting equal-key pairs.

## 4. The residual object is a punctured four-cube

Every final pair satisfies `v=sigma u`, so over `F_2`

```text
D=(1+sigma)u,       sigma D=D,                           (9)
```

and `D` vanishes on reflection-fixed cells.  Intersecting (9) with the
fourteen-dimensional syndrome code leaves a six-dimensional mirror-syndrome
space: the size-one layers `tau=3,4` vanish; `tau=5,6` contribute one even
direction each; and `tau=7,8` contribute two each.

The eleven actually occupied difference sectors span a four-dimensional
subspace.  One basis is

```text
192486, 8c2c0c, 11b4600, 4483414.                       (10)
```

All three face projections retain rank four:

| `D` | `d_A` | `d_B` | `d_C` |
|---|---|---|---|
| `192486` | `1a103` | `34a46` | `03249` |
| `8c2c0c` | `48306` | `98e0c` | `11858` |
| `11b4600` | `9e980` | `35300` | `2368c` |
| `4483414` | `10850a` | `50a14` | `89068` |

Thus no single face projection explains the dimension drop; it is a
compatibility phenomenon in the `H_8` kernel and literal layer data.

The 11 occupied sectors have fibre multiplicities

```text
192486:6  8c2c0c:2  95088a:2  11b4600:2  18e4e8a:4
1976a0c:4  4483414:4  4511092:2  4c41818:26
5c67a9e:2  5df5e18:4.                                  (11)
```

Their Hamming weights are `8,12,16` with pair counts `44,8,6`.  Four nonzero
points of the span are missing:

```text
1026286, 4dd3c9e, 54a5692, 5537214.                     (12)
```

In basis coordinates (10), these are `0101,1011,1101,1100`.  The first has
empty `A` and `C` collision lists.  The remaining three have exact gate
counts

```text
D         |L_A| |L_C|   AC     B   final S2
4dd3c9e       8     8     4     4       0
54a5692    1028  1028  7944   504       0
5537214      12    12     4     4       0.              (13)
```

So one puncture is imposed by metagraph-face structure and three by the
failure of the particular compatible bases to realize raw-S2 equality.  The
middle puncture alone accounts for 504 of all 636 B-compatible candidates;
it dominates the pre-S2 mass and its subsequent collapse.

Join occupied basis-neighbours in the four-cube.  The resulting syndrome
graph is connected, has 14 edges, diameter five, degree histogram
`{1:1,2:4,3:5,4:1}`, and coordinate-edge counts `(4,3,3,4)`.  This is an
objective recursive adjacency and ordering scaffold for the defect sectors;
it is not a claim that these coordinate edges are blue or black metagraph
edges.

## 5. One antisymmetric bit is necessary and sufficient

Number the positions in the low half of layer `tau` from left to right and
define the skew moment

```text
J_tau(x)=sum_p p (x_low(p)-x_high(p)),
J(x)=(J_3(x),...,J_8(x)).                                (14)
```

Reflection sends `J` to `-J`.  Every one of the 58 false-palindrome pairs has
nonzero `J`, and the first nonzero component is at `tau=5,6,7` for respectively
`38,12,8` pairs.  Therefore the single bit

```text
eta(x)=1 if the first nonzero component of J(x) is positive, else 0          (15)
```

separates every doubleton.  The default in the zero case is immaterial
because no collision endpoint has `J=0`.  Thus `(Lambda_9,eta)` is injective
on all `2^27` oriented lines.  At least one bit is information-theoretically
necessary because (1) is not injective, so this repair is minimal.

THM-825's full first positional moments also separate every pair.  Equation
(14) identifies the smaller datum actually missing: antisymmetric position.
A total skew sum still vanishes on 42 pairs, and a clock-weighted skew sum on
eight.  Thus neither of these two natural contractions suffices; the
lexicographic sign is the certified uniform choice used here.

## 6. Black-curvature disintegration boundary

All 58 pairs have `UABC=KKKK`.  Each of their `A,B,C` endpoint-node pairs is
a cross rather than a loop, and all 58 ordered six-node face signatures are
distinct.  Hence the failure is not concentrated at one repeated metagraph
node tuple.

More sharply, the full THM-811 tuple `(q_0,q_1,S,epsilon)` is identical at
the two endpoints of every collision and every pair lies on the Smith balance
wall `epsilon=0`.  Their `(q_0,q_1)` distribution is

```text
(0,0):13 (0,1):5 (0,2):8 (1,1):5 (1,2):2 (1,3):2
(2,3):6  (3,3):3 (3,4):2 (3,5):10 (4,5):1 (5,7):1.     (16)
```

Curvature successfully disintegrates black flow in THM-811/827, but it
cannot orient the two sheets exactly on its balance wall.  The missing datum
is chirality, not additional scalar curvature.

## Tournament Analysis and preservation boundary

Use the eleven occupied sectors as vertices.  The pairwise observable is the
tuple `(first nonzero skew layer, Hamming weight, fibre multiplicity)`.  A
structural switch prioritizes earlier skew detection and smaller support; an
empirical switch prioritizes larger fibre mass; hexadecimal `D` is the tie
Hamiltonian path.  Both gauges are transitive tournaments with score
histogram `0,...,10`, no directed triangle, eleven singleton SCCs, and one
Hamiltonian path.  Twenty-two edges flip between the gauges.  This orders the
four-cube sectors without pretending that their auxiliary tournament is the
merged metagraph itself.

The challenged assumption is that a node tuple, scalar curvature profile, or
raw layer count should be the recursive vertex.  The exact carrier here is a
reflection-defect sector together with its fibre of literal gluing witnesses.
It preserves the raw address kernel, face projections, upper colour, layer
histograms, and reflection action.  It destroys metric LRC gaps, owner labels,
wall-crossing chronology, and centered-continued-fraction phase.  Accordingly
this theorem is a complete `n=9` codec result, not a proof of LRC(14) and not
an all-size finite-state theorem.

The next structural test is to act on the rank-four space (10) by the
centered-CF coordinate-copy maps: decide whether it is invariant, splits into
THM-813 `Q` orbits, or maps into THM-814's positional-curvature kernel.  That
test now has only fifteen nonzero sectors rather than an ambient line scan. ∎
