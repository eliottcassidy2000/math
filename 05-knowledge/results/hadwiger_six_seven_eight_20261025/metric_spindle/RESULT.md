# The unique seven-image unit cube, and why the Moser seven is different

2026-10-25. **PROVED elementary; independently audited.** The
finite receipts are **FINITE-EXACT**. The [independent audit](../audits/METRIC_AUDIT.md)
uses a separate two-circle/K2,3 obstruction and Gram-matrix checks.
No improvement of a plane chromatic
bound is claimed. “Unit-edge drawing” below requires every graph edge to
have Euclidean length one, but permits vertex collisions and extra unit
distances on nonedges unless injectivity or inducedness is explicitly stated.

## 0. Inheritance and the live board

The closest proved mechanism is the square-face/character theorem in
[ports_regular_actions_hn, geometry §§2--4](../../ports_regular_actions_hn_20261025/geometry/README.md):
every injective planar unit cube is a subset-sum drawing, and its centered
coordinates lie in the three singleton Walsh modes. The canonical hostile
is a dependent cubelike direction set, which cannot be drawn injectively
with unit edges. The corrected near miss is treating a graph quotient as
preserving chromatic number or Euclidean metric. The least-used sidecar is
the **collision fibre**, not just the number of image points.

We also inherit the
[Boolean empty/full quotient](../../reflection_bridge_20261025/seven.md),
the [Cremona geometric repair](../../six_seven_geometry_20261025/support_reflections/RESULT.md),
and Session61's [A2 metric and flower](../../session61_hadwiger_six_seven_20261025/hexagonal/RESULTS.md).
Session61 already proves that the flower has twelve unit edges and chromatic
number three. The new-to-route assertion here is the **classification of
all unit-edge cube drawings with exactly seven distinct images**, which
forces that flower and its collision/action data.

Anchor: classify the metric eight-to-seven transition. Niche: retain the
K3 conic-cube relation modes through the collapse. Wildcard: compare the
recovered Moser spindle with the flower using metric field/rank and complete
colour constraints. The board is cube square faces / collision fibres /
A2 Gram matrix / descended actions / two-diamond hinge / palette relation.
No “seven vertices therefore seven colours” inference is entertained.

![The metric cube-to-flower quotient](cube_to_flower.svg)

The source cube in this figure is schematic; the target uses actual regular
hexagon coordinates. `figure.py` regenerates the SVG and optional PNG. The
PNG was visually inspected, and unsupported subscript glyphs were replaced
with plain A2/Q3 labels before delivery.

## 1. A sharp seven-image theorem

Write `Q3` on `{0,1}³`, with Hamming-distance-one edges.

**Theorem M1.** Let `f:{0,1}³ -> R²` carry every cube edge to
a segment of length one. If f has exactly seven distinct values, then,
after a cube automorphism and a Euclidean isometry, it has the form

```
f(epsilon)=epsilon_1*u+epsilon_2*v+epsilon_3*w,
|u|=|v|=|w|=1,           u+v+w=0.                    (1)
```

Thus exactly one opposite pair is identified, and the seven-point image
is the center and six vertices of a regular unit-side hexagon. Conversely
(1) is such a drawing. This is a uniqueness theorem for the **whole plane**
drawing class, not an enumeration inside a chosen lattice or coordinate ansatz.

**Proof.** Exactly seven values on eight inputs means one fibre has size
two and all other fibres are singletons. Its two vertices cannot be adjacent,
because an edge cannot collapse. Their Hamming distance is two or three.

First recall the planar unit-quadrilateral lemma. If four distinct points
`p0,p1,p2,p3` occur cyclically on a unit four-cycle, p1 and p3 are the two
distinct intersections of unit circles around p0 and p2. They are symmetric
about the midpoint of p0,p2. Hence

```
p0+p2=p1+p3.                                         (2)
```

The statement does not require a crossing-free drawing. It does require
four distinct points and the Euclidean plane.

Suppose the one collision has Hamming distance two. Normalize it to
`f(100)=f(010)`. The two faces
`000,100,101,001` and `000,010,011,001` each have four distinct images,
since neither contains both members of the unique collision. Applying (2)
to both gives

```
f(101)=f(100)+f(001)-f(000)
      =f(010)+f(001)-f(000)=f(011),
```

a second collision. This contradicts the seven-image hypothesis.

The collision is therefore opposite; normalize it to `000,111` and translate
their common value to zero. Every cube face now has four distinct images,
so (2) holds on every face. Set `u=f(100),v=f(010),w=f(001)`.
The three faces through zero give `f(110)=u+v`, `f(101)=u+w`,
`f(011)=v+w`; any remaining face gives `f(111)=u+v+w`. Therefore their
collision says `u+v+w=0`. All three vectors have unit length. Consequently

```
u.v=u.w=v.w=-1/2,
Gram(u,v,w)=[ 1  -1/2 -1/2 ]
            [-1/2  1  -1/2 ] .                      (3)
            [-1/2 -1/2  1  ]
```

The directions are 120 degrees apart. Their six signed values are the
vertices of a regular hexagon, and all proper nonempty subset sums are
precisely those six values. This proves the claim and converse. QED.

Equation (3) has eigenvalues `3/2,3/2,0`, with null vector `(1,1,1)`.
It explains the six equiangular directions and seventh center through one
specific linear relation. Six arbitrary direction labels do not suffice.

## 2. The image graph and the exact information lost

Use the Eisenstein basis with norm `Q(a,b)=a²-ab+b²`, and choose

```
u=(1,0), v=(0,1), w=(-1,-1).
```

Then f identifies only `000` with `111`, and its other values are the six
Eisenstein units. Cube complement `epsilon -> 1-epsilon` becomes the actual
Euclidean central inversion `p -> -p`. The empty/full Boolean quotient is
therefore realized **metrically**, with its former 3+1+3 action unchanged.
As before, it is not an ordered-lattice quotient and does not restore the
four Fano lines destroyed by modified complement.

The graph quotient `Q3/(000=111)` is exactly the wheel consisting of a
six-cycle and its center. Every one of the twelve cube edges has a different
edge image: three spokes arise from 000, three from 111, and six rim edges
from the two middle Boolean levels. Among the 21 unordered pairs of distinct
image points, the squared-distance distribution is

```
1 : 12 pairs,      3 : 6 pairs,      4 : 3 pairs.      (4)
```

Thus the quotient graph is the **full induced unit-distance graph** of the
image. There are no additional image edges hidden by the construction.

However f does not reflect adjacency on its original eight labelled inputs.
For example `000,110` is a cube nonedge, but their images are unit-separated;
that image edge is already the image of `111,110`. On the 28 labelled input
pairs the squared-distance distribution is instead

```
0 : 1 pair,   1 : 18 pairs,   3 : 6 pairs,   4 : 3 pairs.  (5)
```

The cube has chromatic number two. The flower has chromatic number three:
each center/rim-edge triple is a triangle, and the even rim cycle uses two
colours while the center uses a third. The quotient has identified opposite
parities of the bipartition, so the old two-colouring does not descend.
This is a finite metric counterpart to Session61's three-colour triangular
cover with seven-colour quotient. Quotient colourings pull back; the converse
is not valid without fibre-constant colourings.

## 3. Exactly D12 descends, and the K3 relation space remains visible

Every graph automorphism of Q3 is a coordinate permutation followed by a
bit translation. One proof: the image of 000 has eight choices and its
three neighbours have six orders; all remaining vertices are determined by
their incidences to those neighbours. The resulting 48 maps are exactly
`F2³ semidirect S3`.

Such a map descends through the unique collision iff it preserves the pair
`{000,111}`. Its translation must be 000 or 111, while any coordinate
permutation is allowed. The descended group is therefore

```
C2 x S3 ≅ D12,               |D12|=12.               (6)
```

It acts by all Euclidean symmetries of the regular hexagonal flower.
The six signed directions are preserved, coordinate permutations permute
`u,v,w`, and the nontrivial translation becomes central inversion.
Only **two of the eight** regular cube translations descend. The other
six move the collapsed fibre and cannot act on this seven-point quotient.
Regularity of a full deck translation action is consequently lost.

For the inherited eight-conic intersection matrix
`M=-2I+2A(Q3)`, the three singleton Walsh modes span its nullspace.
Since `u+v+w=0`, equation (1) may also be written

```
f(epsilon)=-1/2 sum_j (-1)^epsilon_j u_j.             (7)
```

Each of the two coordinate functions is therefore in `ker M`, as in the
injective cube theorem. The extra coefficient relation picks the unique
seven-image degeneration. This transports a labelled coefficient-space
relation and cube adjacency. It does not identify the indefinite conic
intersection matrix with the positive-semidefinite Gram matrix (3), nor
embed a K3 surface in the Euclidean plane. Intersection multiplicity,
curve classes and complex structure remain outside this metric map.

A faithful eight-point control makes the boundary explicit. The subset sums
of the three rational unit vectors

```
(1,0), (3/5,4/5), (5/13,12/13)
```

are eight distinct points, whose induced unit graph is exactly Q3. All
28 pair distances are checked. Planarity of the ambient Euclidean space
does not force the eight-to-seven collision; relation (3) does.

## 4. The recovered Moser seven escapes through a different operation

The actual Moser mechanism is already proved and authenticated in
[spindle_tournament_interfaces, §2](../../spindle_tournament_interfaces_20261011g.md)
and [unit_field_escape_recursion, §§2--4](../../unit_field_escape_recursion_20261011f.md).
We recover it as a metric hostile to unqualified claims about seven points,
six directions, or dense edge counts, not as a new construction.

Let

```
alpha=(1+i*sqrt(3))/2,    beta=(5+i*sqrt(11))/6.
```

The seven points

```
0, 1, alpha, 1+alpha, beta, alpha*beta, beta*(1+alpha)  (8)
```

have exactly eleven unit edges: two copies of the five-edge diamond
`K4 minus one edge`, joined at 0, and the edge between their outer tips
`1+alpha` and `beta*(1+alpha)`. The hinge equation is exact:

```
|1+alpha|²=3,     |1-beta|²=1/3.
```

In a proper three-colouring, the two tips of a diamond have the same colour:
their adjacent common neighbours already use the other two colours.
Both outer tips in (8) would therefore match 0, contradicting their joining
edge. A four-colouring in the displayed point order is
`(0,1,2,0,1,2,3)`. Hence its exact chromatic number is four.

The lost coordinate is also arithmetic. Put `K=Q(i*sqrt(3))` and
`Lambda=Z[alpha]`. Each diamond lies in one copy of the rank-two triangular
lattice, but the rotated copy is `beta*Lambda`. Since `beta notin K`,

```
K intersect beta*K={0},
Lambda + beta*Lambda is a direct sum of Z-rank four.   (9)
```

Indeed a nonzero common element would give beta as a ratio of two elements
of K. The four elements `1,alpha,beta,alpha*beta` are independent over Q;
their coefficient determinant in the basis `1,s,t,st`, with `s²=-3,t²=-11`,
is `1/144`. Thus this point set cannot lie in any translated discrete
rank-two lattice. This is an actual escape from the single-lattice cap,
not merely an irrational angle: Session61's
[dense-flower control](../../session61_hadwiger_six_seven_20261025/dense-flower/DENSE_FLOWER_CONTROL.md)
has an infinite-order rotation still inside K and remains three-colourable.

The six equiangular steps of either diamond live in its own A2 lattice.
The rotation beta and the closing edge break a global three-direction cube
translation rule. Local unit rhombi still satisfy (2), but their overlap
pattern is no longer the six square faces of one cube. This is the precise
operation distinguishing the seven-image flower from the Moser seven.

Nor does field escape by itself force ever more colours: the inherited
all-length spindle family remains four-colourable even outside prescribed
finite collections of residue screens. The diamond is an equality gadget
at palette three; at palette four or six its tips can differ. That palette
boundary is retained, not extrapolated to the plane's six/seven question.

| Exact unit graph | Vertices | Edges | Chromatic number | Mechanism |
|---|---:|---:|---:|---|
| Faithful rational cube | 8 | 12 | 2 | Generic unit subset sums |
| Unique seven-image cube fold | 7 | 12 | 3 | Opposite-parity collision; A2 relation |
| Moser spindle | 7 | 11 | 4 | Two diamond equalities plus metric hinge |

For classical context, de Grey's
[*The chromatic number of the plane is at least 5*, v3](https://arxiv.org/html/1804.02385),
§§1--2 and §4, treats the regular seven-point hexagon and Moser spindle as
different ingredients. That primary source was checked 2026-10-25. The
present comparison does not use its larger five-chromatic construction or
make any assertion that its historical lower bound is the current one.

## 5. Connection contract, controls and stopping reason

| Source -> target | Map / preserved predicate | Loss / necessary sidecar |
|---|---|---|
| Eight cube vertices -> seven A2 points | (1); all twelve edges unit, quotient graph induced | One opposite fibre; six extra labelled unit nonedges in pullback |
| Boolean complement -> plane involution | Subset sum; complement becomes p↦-p | Inclusion and Fano incidence do not follow |
| Cube automorphisms -> flower isometries | Stabilizer of the collision fibre | Six regular translations excluded |
| K3 conic relation modes -> metric coordinates | (7), singleton Walsh span | Intersection form/multiplicity and complex geometry forgotten |
| Two A2 diamonds -> Moser graph | beta rotation plus one actual unit edge | One common lattice and global cube rule lost; palette-three equality needed |

The decisive improvement is M1's all-drawing classification and the exact
group/metric loss ledger. The next step toward stronger chromatic constraints
would need attachments with a proved joint boundary-colouring relation;
changing seven into eight, retaining a free group action, or adding another
rotation does not supply one. No new high-palette equality gadget is claimed.

Reproduce from this directory:

```
python3 check.py > normal.json
python3 -O check.py > optimized.json
```

The standard-library-only checker gives 297 explicit, optimization-safe
checks. Its complete finite universes are: all 28 possible cube collision
pairs; all square faces remaining injective under each candidate collision;
all 21 distances on the seven-image set and all 28 on the faithful cube;
all 48 cube automorphisms; all 21 distances of the Moser set; and every
labelled two-/three-/four-colouring used for the stated small graphs.
Rational row reduction independently verifies collision propagation; exact
arithmetic in `Q(sqrt(-3),sqrt(-11))` checks Moser distances and rank.
The checker does not assume a numerical tolerance or infer M1 from samples.
