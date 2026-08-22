---
id: THM-3462
title: "Five-point line-metric interval-sum preorder atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The ten labelled
  distances of five ordered line points realize exactly 477 positional total
  preorders.  Their tie census is
  114,162,124,36,5,14,16,3,2,1 in the stated block types; the projective
  arrangement f-vector is (20,125,218,114).  Reversal gives 241 classes with
  five fixed signatures; S5 gives no additional positional collapse and
  28,620 arbitrary labelled signatures.  Primitive gap representatives have
  sharp universal total 23 and coordinate height 10.  No turnpike, tournament,
  LRC, Rule 30, Keller, or Jacobian conclusion follows.
source: codex-2026-08-15-five-point-interval-sum-atlas
audit: >
  two algorithmically independent exact routes in one stdlib companion:
  global 4-by-4 Cramer vertex enumeration over 49 affine equations and
  recursive exact Fourier--Motzkin wall insertion; full normal/optimized/stored
  replay, AST/security, dependency/hash, symmetry, hostile, documentation, and
  scope gates; independent compressed-Cramer, Fourier--Motzkin-prefix,
  S5-orbit, representative, hostile-sampling, replay, hash, security,
  documentation, and routing audit PASS.
depends_on:
  - THM-3457-four-point-line-metric-turnpike-preorder-atlas
related:
  - THM-3454-fibonacci-selected-u-spine-farey-lorentz-isometry-and-one-tie-edge-order
script: 04-computation/five_point_line_metric_interval_sum_preorder_atlas_thm3462.py
output: 05-knowledge/results/five_point_line_metric_interval_sum_preorder_atlas_thm3462.out
script_sha256: ac057c7acce3794ba7772bf548485436eee6f228633581e37d1c762d89892a95
output_sha256: bf28b02ce29de0bd9c83fb47207b8d1b76999939fa0e49bfe9cda1acf950ed16
covector_sha256: 38c1e9f857ebfe1ab5473bcab15695500fa5de7293cba436924199b4f0024eff
signature_sha256: ee30ec8827977bc514ab4b7430fd6f71c6eb9464694ebae865cd6d922f4e0ca6
semantic_sha256: 90323a5738e47db6643fea438d3a6b7e3d744b6119568f1e35b2bd03b3297d88
hash_basis: LF-normalized bytes
---

# THM-3462 -- five-point line-metric interval-sum preorder atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This is a complete elementary classification of the labelled distance
preorders of five ordered points on a line.  It is not a literature-priority
claim and not a reconstruction theorem for an unlabelled distance multiset.

## 1. Typed setup and inheritance

Let

```text
x_0<x_1<x_2<x_3<x_4,
a=x_1-x_0,  b=x_2-x_1,  c=x_3-x_2,  d=x_4-x_3,          (1)
```

with `a,b,c,d>0`.  The ten labelled distances are the nonempty contiguous
interval sums:

| edge | indicator in `(a,b,c,d)` | distance |
|---|---:|---|
| `01` | `1000` | `a` |
| `12` | `0100` | `b` |
| `23` | `0010` | `c` |
| `34` | `0001` | `d` |
| `02` | `1100` | `a+b` |
| `13` | `0110` | `b+c` |
| `24` | `0011` | `c+d` |
| `03` | `1110` | `a+b+c` |
| `14` | `0111` | `b+c+d` |
| `04` | `1111` | `a+b+c+d` |

Increasing distance defines a total preorder on these ten edge labels.  A
strict chamber gives a transitive tournament on ten comparison vertices;
on a wall, equal distances remain tied rather than receiving an artificial
orientation.

The validity card is therefore:

| field | declaration |
|---|---|
| geometric vertices | the ordered line points `0<1<2<3<4` |
| comparison vertices | the ten labelled edges of `K5` |
| pairwise observable | comparison of their positive interval-sum lengths |
| orientation | smaller distance points to larger distance |
| ties | missing strict arcs, or bidirected arcs in the weak relation |
| preserved target | the complete labelled distance preorder |
| lost | magnitudes, common scale, translation, reflection gauge, and recurrence origin |

The closest proved mechanism is
[THM-3457](THM-3457-four-point-line-metric-turnpike-preorder-atlas.md),
which gives the corresponding `25`-preorder atlas for four points.  The two
canonical hostile chambers here are `(6,5,4,8)` and `(2,3,4,10)`: they force
the sharp representative bounds.  The corrected near miss is the general
arrangement warning in [MISTAKE-223](../MISTAKES.md): an Euler count or common
wall language does not transport another problem's predicate.  The retained
sidecar is the full fifteen-coordinate wall covector, not merely its number of
ties or its strict tournament shadow.

## 2. Exactly fifteen equality walls

There are `binom(10,2)=45` distance comparisons.  Two interval indicators are
coordinatewise comparable exactly when one interval contains the other, and
then positivity forces strict inequality.  An interval of length `L` contains
`L(L+1)/2-1` proper nonempty subintervals.  Hence the number of forced pairs is

```text
sum_(L=1)^4 (5-L)(L(L+1)/2-1) = 0+6+10+9 = 25.          (2)
```

The other `20` comparisons are mixed.  Cancelling the overlap of two
contiguous intervals turns equality into equality of two disjoint nonempty
contiguous gap blocks.  Conversely each such block equality occurs.  After
orienting the leftmost nonzero coefficient positively, the `20` mixed pairs
deduplicate to exactly these `15` walls:

```text
a=b,       a=b+c,       a=b+c+d,
a=c,       a=c+d,       a=d,
a+b=c,     a+b=c+d,     a+b=d,
b=c,       b=c+d,       b=d,
a+b+c=d,   b+c=d,       c=d.                            (3)
```

Some walls compare more than one pair of distances; this explains `20>15`.
There are no other equality walls, because every remaining comparison is one
of the strict containments counted in `(2)`.

Let `h_1,...,h_15` be the primitive integer forms in `(3)`, in the companion's
frozen lexicographic coefficient order.  A positional atlas cell has covector

```text
sigma=(sign h_1,...,sign h_15) in {-1,0,+1}^15.         (4)
```

The containment comparisons plus `(4)` determine the complete ten-distance
preorder, and every positive gap vector determines one such covector.

## 3. Completeness without a bounded grid

### 3.1 Global Cramer-vertex certificate

Fix a proposed covector `sigma`.  Its homogeneous open stratum is nonempty if
and only if the following closed rational polyhedron is nonempty:

```text
x_i >= 1                                      for i=1,...,4,
sigma_j h_j(x) >= 1                           if sigma_j != 0,
h_j(x)=0                                      if sigma_j = 0.   (5)
```

Indeed, a point in the open stratum can be multiplied by a positive scalar so
that all four gaps and every prescribed nonzero wall value have absolute value
at least one.  The converse is immediate.

The polyhedron `(5)` is pointed: a vector in both directions of its recession
cone has every coordinate zero because of the four lower bounds.  Therefore a
nonempty `(5)` has a vertex.  Four linearly independent active equations at a
vertex are selected from

```text
x_i=1                         (4 equations),
h_j=-1,  h_j=0,  h_j=1        (3*15 equations),          (6)
```

for a global bank of `49` affine equations.  Thus every feasible covector is
seen among the `binom(49,4)=211,876` four-equation systems.  Conversely every
positive solution of one of these systems supplies an actual feasible
covector, whether or not that solution is a vertex for every covector it
realizes.

Exact integer Cramer arithmetic gives:

| gate | count |
|---|---:|
| four-equation systems | `211,876` |
| nonsingular systems | `123,706` |
| positive rational solutions, with multiplicity | `17,143` |
| distinct feasible covectors | `477` |

This proves completeness of the atlas independently of any bounded search for
small integer gaps.

### 3.2 Independent Fourier--Motzkin route

The companion separately begins with the empty wall word and inserts the
fifteen signs one at a time.  For every partial word it encodes positive and
negative signs by margin-one inequalities and a zero sign by the two opposite
weak inequalities.  It then eliminates all four variables by exact
Fourier--Motzkin pairing, normalizing positive scalar multiples and retaining
the strongest duplicate inequality.

This route makes `7,053` exact feasibility decisions.  Its surviving prefix
counts are

```text
1,3,9,15,19,25,75,125,175,213,253,283,339,389,427,477. (7)
```

It returns the same `477` terminal covectors as the Cramer route.  With the
walls lexicographically sorted and the covectors serialized as compact sorted
JSON, their SHA-256 is

```text
38c1e9f857ebfe1ab5473bcab15695500fa5de7293cba436924199b4f0024eff. (8)
```

The two routes share only the derived wall matrix and exact arithmetic.  One
uses global active-equation vertices; the other decides recursively by
projection.  A later independent audit reconstructed a compressed Cramer
solver, every Fourier--Motzkin prefix set, the tie/face/reversal/S5 counts,
both sharp representative bounds, and the Fibonacci hostile without modifying
the package.

## 4. Tie census and projective face vector

For each covector, evaluate the ten interval sums at its exact Cramer witness
and group equal values.  The covector determines the resulting blocks by
Section 2, so this is a classification rather than sampling.  Listing only
the nonsingleton block sizes gives:

| tie-block type | count |
|---|---:|
| strict | `114` |
| `2` | `162` |
| `2+2` | `124` |
| `2+2+2` | `36` |
| `2+2+2+2` | `5` |
| `3+2` | `14` |
| `3+2+2` | `16` |
| `3+2+2+2` | `3` |
| `3+3+2` | `2` |
| `4+3+2` | `1` |
| **total** | **`477`** |

The unique `4+3+2` row is `a=b=c=d`: the four length-one intervals tie, the
three length-two intervals tie, the two length-three intervals tie, and `04`
is the unique diameter.

If `Z(sigma)` is the set of zero wall rows, the relative projective dimension
of a stratum is

```text
dim(sigma)=3-rank_Q Z(sigma).                            (9)
```

Exact row reduction gives the face vector, from projective dimension zero to
three,

```text
(f_0,f_1,f_2,f_3)=(20,125,218,114).                     (10)
```

The independent Euler control is

```text
20-125+218-114=-1,                                      (11)
```

the compactly supported Euler characteristic of the open three-simplex.

## 5. Reversal, `S5`, and arbitrary labellings

Line reversal acts by

```text
(a,b,c,d) -> (d,c,b,a),
ij -> (4-j)(4-i).                                       (12)
```

A reversal-fixed total preorder must put each reversed pair of labels in the
same ordered block: an order-preserving permutation of a finite chain of
blocks fixes every block.  In particular `01=34` and `12=23`, so

```text
a=d=p,  b=c=q.                                          (13)
```

Conversely `(13)` makes every distance reversal-invariant.  On this slice the
only additional critical ratios are `p=q` and `p=2q`, giving exactly five
strata:

```text
p<q,  p=q,  q<p<2q,  p=2q,  p>2q.                      (14)
```

Thus reversal fixes exactly five of the `477` signatures, and Burnside gives

```text
(477+5)/2=241                                           (15)
```

reversal classes.

No arbitrary vertex permutation identifies two further positional classes.
The edge `04` is the unique largest distance.  After its unordered endpoints
are identified, choosing either endpoint makes its four incident distances a
strict chain; that chain recovers the complete order of the remaining
vertices.  The two choices give precisely the original order and full line
reversal.

Consequently `S5` stabilizers have order one on the `236` non-fixed reversal
classes and order two on the five fixed classes.  Their arbitrary-label orbits
have sizes `120` and `60`, respectively, and are pairwise disjoint.  The exact
number of distance-preorder signatures after arbitrary labelling is therefore

```text
236*120+5*60=28,620.                                    (16)
```

This classifies labelled preorders, not unlabelled distance multisets and not
point-set reconstruction from such a multiset.

## 6. Sharp primitive representative bounds

Only after the `477` covectors were certified by Sections 3.1 and 3.2, the
companion enumerated positive integer compositions to optimize one primitive
representative per known cell.  The cumulative total-bound coverage is

```text
total <=21: 475,     total <=22: 475,     total <=23: 477. (17)
```

The final two cells are reversals.  One has representative `(6,5,4,8)`, and
its signature forces

```text
c<b<a<d,       d<b+c,       a+b<c+d.                    (18)
```

For positive integers, `(18)` gives

```text
d >= a+b-c+1,       d <= b+c-1,       hence a <= 2c-2.
```

Since `c<b<a`, also `b>=c+1` and `a>=c+2`; therefore `c>=4`, `b>=5`,
`a>=6`, and then `d>=8`.  Every integer representative of this cell has total
at least `23`, attained by `(6,5,4,8)`.  Hence the universal primitive total
bound `23` is sharp.

Similarly, the cumulative coordinate-height coverage is

```text
max gap <=9: 465,       max gap <=10: 477.               (19)
```

The chamber represented by `(2,3,4,10)` forces

```text
a<b<c<a+b,       a+b+c<d.                               (20)
```

If `a=1`, the integer inequalities `b<c<a+b` are impossible.  Thus `a>=2`,
then `b>=3`, `c>=4`, and `d>=a+b+c+1>=10`.  The displayed representative
attains height `10`, proving sharpness.  A first representative at a minimal
total or height is automatically primitive, since dividing a common factor
would preserve its covector and improve the bound.

These finite sweeps optimize representatives inside an already complete exact
atlas; they are not the completeness argument.

## 7. Fibonacci two-wall lane and its hostile

For four consecutive positive Fibonacci gaps

```text
(a,b,c,d)=(F_(k-3),F_(k-2),F_(k-1),F_k),   k>=5,
```

one has

```text
c=a+b,       d=b+c,       a<b.                           (21)
```

They occupy atlas row `C054`, with signature

```text
01<12<23=02<34=13<03<24<14<04.                          (22)
```

At the initial boundary `(1,1,2,3)`, the extra wall `a=b` gives row `C057`
and tie type `2+2+2`.  The converse to `(22)` fails sharply: `(1,3,4,7)` has
the same two recurrence walls and the same preorder, while its Cassini--Pell
value

```text
b^2-ab-a^2=5
```

is not a Fibonacci unit `+1` or `-1`.  As in THM-3454/3457, the preorder
retains additive walls but forgets recurrence origin and the Pell sidecar.

## 8. Exact companion and controls

The companion uses only standard-library integer and `Fraction` arithmetic.
It prints all `477` covectors with stable IDs, projective dimensions, tie
types, sharp representatives, reversal partners, and complete labelled
distance signatures.  Its main gates are:

- derivation of the ten interval indicators, the `25+20` comparison split,
  and the fifteen disjoint-block walls;
- all `211,876` exact Cramer systems and all `7,053` independent
  Fourier--Motzkin decisions;
- equality of both `477`-covector sets and the frozen digest `(8)`;
- the complete tie census, face vector, Euler control, and equal-gap maximum;
- all `120` vertex permutations on each reversal representative, with no
  collapse beyond reversal;
- post-certification sharp representative sweeps and both hostile chambers;
  and
- an AST self-audit allowing only named standard-library imports, with zero
  `assert`, dynamic-evaluation, write, or network calls.

Reproduction commands are

```text
python -B 04-computation/five_point_line_metric_interval_sum_preorder_atlas_thm3462.py
python -B -O 04-computation/five_point_line_metric_interval_sum_preorder_atlas_thm3462.py
```

Normal, optimized, and stored standard output agree byte for byte with LF
newlines.  Source, output, covector, signature, and semantic hashes are frozen
in the frontmatter.  The script also freezes the LF hash of THM-3457 at
`1e0acb07...c01d5`.

## 9. Scope and non-consequences

| field | exact boundary |
|---|---|
| source | five distinct labelled points on the real line in one fixed order |
| target | total preorder on their ten labelled distances |
| quotient | common translation and positive scale; reflection only after orbiting |
| preserved | every distance comparison, equality block, wall sign, and edge label |
| destroyed | magnitudes, scale, origin, ancestry, Pell norm, phase, owner, and current |
| sharp total hostile | `(6,5,4,8)` and its reversal |
| sharp height hostile | `(2,3,4,10)` |

The strict rows are transitive `T10` comparison tournaments, but the boundary
rows are honest preorders with missing or bidirected tie arcs.  Neither object
is “a tournament of size five.”  The theorem does not solve the classical
turnpike reconstruction problem, classify all `n`-point line metrics, or
identify a recurrence from its edge order.  It supplies no LRC owner, physical
clock, spectral current, Rule 30 trace, Keller map, Jacobian degree exclusion,
or harmonic-density theorem.  LRC(14) and `JC(2)` remain open.
