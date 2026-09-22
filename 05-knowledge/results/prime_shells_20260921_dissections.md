# Exact simplex dissections, parity, and the two four-cycles

**Status:** PROVED elementary geometry and finite conjugacy; FINITE-EXACT
companion controls. The finite conjugacy is not a Collatz convergence theorem,
a Berggren-tree conjugacy, or an arithmetic group isomorphism. No novelty claim.

## 1. Inheritance and the five-concept board

The closest inherited construction is the alternating-vertex tetrahedron in
[`simplex_cuboid_packing.out`](simplex_cuboid_packing.out), §§1–4: its central
volume is `1/3`, its four corner volumes are `1/6`, and dimension four already
destroys the simplex identification. Those claims are re-proved below.
The current proved sidecar is **THM-3991, periodic-unimodular-toric-cusp-factorial-euler-obstruction**
([file](../../01-canon/theorems/THM-3991-periodic-unimodular-toric-cusp-factorial-euler-obstruction.md), §4):
the five cells have normalized volumes `(2,1,1,1,1)`, so counting cells loses
the lattice-volume weight. Its toric conclusions are not imported here.

The corrected near miss is **THM-1400, dsat-is-the-known-diameter-and-the-halfcube-correction**
([file](../../01-canon/theorems/THM-1400-dsat-is-the-known-diameter-and-the-halfcube-correction.md), §II):
distance-two moves preserve two parity sectors, but that label need not descend
through an additional quotient. The hostile in §6 below makes this loss concrete.
The arithmetic input is **OS-T4–T6** in
[`odd_square_20260921_triangles.md`](odd_square_20260921_triangles.md):
`D_u(v)=|u−2v|` is multiplication by two on unit residues modulo sign.

| Concept | Retained invariant / operation | Missing coordinate / cheap test |
|---|---|---|
| Central triangle and its ears | Area, rigid half-copy congruence | Rectangle aspect; test `h/w=√3/2` |
| Five tetrahedra | Incidence and volume ratios | Shape and lattice volume; compare squared edges |
| Demicube | Vertex parity, corner cuts | Simplex status; test dimension four |
| `u=17` shell | Permutation `D`, quadratic character | Full multiplier action and metric; test multiplier three |
| Three tetrahedral matchings | Quotient `S4/V4=S3` | No parent or height; compare with Berggren generation |

Anchor: the exact dissection. Niche: its all-dimensional corner-cut extension.
Wildcard: an explicit cube realization of the arithmetic two-four-cycle system.
The successful bridge is a labelled dynamical conjugacy with a precise operation
boundary, rather than a numerical coincidence between piece counts.

## 2. The equilateral picture is a rectangle; the square has an affine repair

**PS-D1.** For `w,h>0`, the rectangle `[0,w]×[0,h]` is the union of

\[
T=\operatorname{conv}\{(0,0),(w,0),(w/2,h)\},
\quad E_0=\operatorname{conv}\{(0,0),(0,h),(w/2,h)\},
\quad E_1=\operatorname{conv}\{(w,0),(w,h),(w/2,h)\}.
\]

Their interiors are disjoint. The altitude halves of `T` and both ears are
congruent right triangles with legs `w/2,h`; hence

\[
\operatorname{area}(T)=wh/2,\qquad
\operatorname{area}(E_i)=wh/4.
\]

Thus “one central triangle plus two halves” is literally correct.
The central triangle is equilateral exactly when
`(w/2)²+h²=w²`, or `h/w=√3/2`. It is therefore inside a rectangle of aspect
`√3/2`, not a square. In the square case `h=w`, the central sides are
`w,√5 w/2,√5 w/2`. An invertible diagonal affine map sends the equilateral
rectangle to a square and retains the entire dissection and its area ratios;
it changes the central triangle's angles and relative side lengths.
If desired the map can preserve area: with `s=√(wh)`, use
`(x,y)↦(sx/w,sy/h)`, whose determinant is one.

**PS-D2 (strong rigid obstruction).** No rigid rearrangement of an equilateral
triangle of side `a` together with two copies of its altitude-half triangle
can tile any square, with these pieces left uncut.

Every piece edge has length in `a·Q(√3)`. Any piece edge that contributes
positive length to a side of the square lies wholly on that supporting line
and therefore wholly on that side. The square side is consequently a sum
of such edge lengths, so `s/a=r+t√3` for rational `r,t`.
The total area requires `(s/a)²=√3/2`. But the rational part of this square is
`r²+3t²`, which is zero only when `r=t=0`; that would give zero area.
This proves the obstruction without an assumption on placement or orientation.
The two altitude halves alone have the same obstruction, now with total area
`√3 a²/4`; they do form a rectangle of aspect `√3` up to reciprocation.
Additional cuts can introduce new edge lengths, so the proof does not rule
out general scissors congruences. Its exact lost coordinate is the edge-length
field; area equality by itself is insufficient.

## 3. The central tetrahedron and four volume-halves really fill the cube

**PS-D3.** Put

\[
P=\operatorname{conv}\{000,110,101,011\}\subset[0,1]^3.
\]

The corner simplex at each unused, odd-parity cube vertex consists of that
vertex and its three cube neighbours. Each has volume `1/6`. The determinant
of the three nonzero vertices of `P` has absolute value two, giving
`vol(P)=1/3`; all six squared edge lengths of `P` equal two.

This is a dissection, not merely a volume identity. Within the cube,

\[
P=\{x\le y+z,\ y\le x+z,\ z\le x+y,\ x+y+z\le2\}.
\]

The reverse of each of these four inequalities describes the corresponding
corner simplex. Two strict reversals cannot hold simultaneously. Every point
outside `P` lies in one of the four corners, proving coverage and disjoint
interiors. Hence `1/3+4(1/6)=1` has its required geometric witness.

Here “half” means half the volume. The corner edge-square multiset is
`{1,1,1,2,2,2}`. A planar bisection of a regular tetrahedron into two
tetrahedra must pass through an edge and the midpoint of the opposite edge;
its half has edge-square multiset `{1/2,3/2,3/2,2,2,2}` at this scale.
The corner is not congruent to either such half. Indeed, for both cut pieces
to be tetrahedra the cutting plane must contain an original edge and meet
the opposite edge; equality of volumes then forces its midpoint.

Every nonsingular affine map transports this dissection to a parallelepiped
and preserves incidences and all five volume ratios. For a rectangular
cuboid with sides `a,b,c`, the central edge-squares are
`a²+b²,a²+c²,b²+c²`, each twice, so regularity survives exactly when
`a=b=c`. The four corner pieces remain congruent for a cuboid; a general
shear can destroy even their mutual congruence. Lattice-unimodularity is
another independent condition: the central normalized volume is two,
although each corner's is one.

**A second hostile.** The even- and odd-vertex central tetrahedra do not fill
the cube together. Their intersection is the octahedron

\[
|x-1/2|+|y-1/2|+|z-1/2|\le1/2,
\]

of volume `1/6`; their union has volume `1/2`. The cube-edge midpoint
`(1/2,0,0)` lies in neither. This repairs the literal “two tetrahedra fill
the cube” sentence in the exploratory
[`simplex_cuboid_eisenstein.out`](simplex_cuboid_eisenstein.out), before that
output's later recomputation. A union of vertex sets is not a union of convex
hulls. Both parities contain four vertices, while their convex hulls overlap.

## 4. The exact higher-dimensional continuation is a demicube

**PS-D4.** Let `n≥2`, let `E_n` be the even-parity vertices of `[0,1]^n`, and
put `P_n=conv(E_n)`. For each odd vertex `v`, remove the open corner
`||x−v||_1<1`. The remaining polytope is exactly `P_n`. The closures of
these `2^(n−1)` corners are congruent orthogonal simplices, with disjoint
interiors and volume `1/n!`. Therefore

\[
\operatorname{vol}(P_n)=1-\frac{2^{n-1}}{n!}.                 \tag{1}
\]

Here is an elementary proof of the polytope assertion. Impose all cube
inequalities and `||x−v||_1≥1` for odd vertices `v`; these are affine
inequalities on the cube. Every even vertex satisfies them, and every odd
vertex fails its own. A feasible nonintegral point cannot have just one
fractional coordinate: it would lie strictly within distance one of the
odd endpoint of its cube edge. Suppose it has `f≥2` fractional coordinates.
If two distinct corner inequalities are tight, their odd vertices have
Hamming distance at least two, while the sum of their distances from `x`
equals two. Thus they differ in exactly two positions, and every other
coordinate of `x` equals their common coordinate. This forces `f=2`, and
the two tight equalities restrict the same diagonal in that square face.
When `f≥3` at most one corner equality is tight. In either case a nonzero
motion within the minimal cube face preserves every active equality, so
the point is not a polytope vertex. All vertices are therefore precisely
the even cube vertices. Finally, two open odd corners cannot overlap by
the triangle inequality and Hamming distance at least two. This proves (1).

For `n≥3` the central polytope is full dimensional and all its
`2^(n−1)` selected cube vertices remain extreme. Consequently it is a
simplex exactly when `n=3`. The ratio of one corner volume to its central
volume is `1/(n!−2^(n−1))`, which equals `1/2` exactly at `n=3`.
At `n=4` there are eight central vertices, central volume `2/3`, and each
corner has volume `1/24`, one sixteenth of the central volume. Dimension
two gives the degenerate diagonal segment, not the equilateral picture.
This retains the inherited corner-cut mechanism and prevents an incorrect
all-dimensional extrapolation of “a simplex with half-copies around it.”

## 5. The `u=17` cycles have an explicit isometric cube model

**PS-D5.** On the eight admissible odd roots `v∈{1,3,…,15}`, `D(v)=|17−2v|`
has the two cycles `(1,15,13,9)` and `(3,11,5,7)`. The cube isometry

\[
G(x,y,z)=(1-y,x,1-z)
\]

has order four, determinant `−1`, and preserves vertex parity. The following
bijection `Φ` exactly satisfies `Φ(D(v))=G(Φ(v))`:

| v | 1 | 15 | 13 | 9 | 3 | 11 | 5 | 7 |
|---|---|---|---|---|---|---|---|---|
| Φ(v) | 000 | 101 | 110 | 011 | 001 | 100 | 111 | 010 |

Moreover `v^8 mod17` is `+1` on the first row-cycle and `−1` on the second,
so the arithmetic quadratic character is exactly `(-1)^(x+y+z)` in this
model. This is a character-preserving finite dynamical conjugacy. Choosing
the images of roots one and three fixes its phases; there are sixteen
character-preserving conjugacies of these two labelled permutation systems.

The cube reflection `C(x,y,z)=(x,y,1−z)` commutes with `G` and swaps parity.
Transporting it back gives the involution

\[
S=(1\ 3)(15\ 11)(13\ 5)(9\ 7),\qquad SD=DS.                 \tag{2}
\]

This is not multiplication by one scalar modulo sign. On
`(Z/17Z)^×/{±1}≅C8`, the only scalar involutions are classes one and four
(odd representatives `1,13`), which both preserve the two sectors.
Multiplier three has order eight and satisfies `M_3²=D³`; it swaps sectors
but is not an involution. Thus the split action `<G,C>≅C4×C2` differs from
the arithmetic multiplier action `C8`, even though both extend the same `D`.

The distinction is geometrically unavoidable: a cube isometry is a signed
coordinate permutation. A signed coordinate cycle of length `j` has order
`j` or `2j`; the partitions of three therefore give only orders
`1,2,3,4,6`. No relabelling realizes the arithmetic order-eight multiplier
as a Euclidean cube isometry. Under the displayed table, multiplier three
sends the cube-edge pair of roots `(1,7)` to `(3,13)`, opposite cube
vertices: squared distance changes from one to three.

## 6. Where the three-way quotient fits, and where it stops

**PS-D6.** The forty-eight cube isometries act on vertex parity by the
character “number of coordinate flips modulo two.” Its kernel has order
twenty-four and acts faithfully on the four even vertices; hence it is `S4`.
Its action on the three perfect matchings of these four vertices has kernel
`V4`, consisting of the identity and the three double transpositions, and
quotient `S3`. This gives exact two-sector and three-matching operations.
Only twelve of these twenty-four parity-preserving isometries preserve
ambient orientation. Vertex parity and orientation sign are different bits.

| Source → target | Map and preserved predicate | Lost data / decisive boundary |
|---|---|---|
| Rectangle dissection → square dissection | Diagonal affine map; incidences, area ratios | Equilateral metric; PS-D2 forbids the uncut rigid version |
| Cube → cuboid | Diagonal affine map; five cells and volume weights | Central regularity unless all side lengths agree |
| `u=17` shell → cube vertices | `Φ`; `D` and quadratic character | Multiplication by three is not a cube isometry |
| Cube symmetry → three matchings | `S4→S3`, kernel `V4` | Four-vertex labels and phase; no height or parent |
| Berggren tree → any finite picture | No conjugacy asserted | Three increasing children and unique decreasing parent cannot become one periodic successor |

The inherited Berggren maps on odd roots are `(u+2v,v)`, `(2u+v,u)`, and
`(2u−v,u)`. They change `u` and generate an infinite rooted tree.
`D_u` keeps `u` fixed and permutes a finite shell. The matching quotient
permutes three finite objects. These operations can coexist on the same
research board, but their cardinalities and recurrence predicates prohibit
an invertible identification of their complete directed graphs. Any proposed
elliptic-tree bridge must likewise provide its own map, height, guards, and
lost-data sidecar. A three-element quotient does not supply three children.

## 7. Replay and scope

Run `python 04-computation/experiments/prime_shells_20260921_dissections.py`
and the same command with `python -O`. The sibling JSON declares the finite
universe and source hash. Checks use exceptions, not removable assertions.
They enumerate polytope vertices exactly in dimensions two through four,
test 8,280 rational cube-grid points with repetition against independent
barycentric and inequality descriptions, check 64 rectangles, enumerate all
48 cube isometries and all 40,320 vertex permutations for the centralizer,
and audit the full eight-state conjugacy and its metric hostile.
The all-dimensional and quadratic-field statements are paper proofs above;
the dimension-two-through-ten table is not their proof. The comparison of
ordinary and optimized execution must match before freezing the output.
