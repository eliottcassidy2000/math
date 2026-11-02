# Hyperbolic binary branching: a bounded curvature control

**Status:** PROVED elementary construction; independent audit pending. The
finite rational checks are corroboration, not the proof for infinitely many
levels. This is an abstract binary-tree witness, not a construction of the full
Collatz graph, a one-step Collatz embedding, or a regular hyperbolic hexagonal
mesh. No novelty or priority claim is made.

## Statement and standard input

Let `T` be the rooted binary tree with vertices `(k,j)`, where `k >= 0` and
`0 <= j < 2^k`; the children of `(k,j)` are `(k+1,2j)` and `(k+1,2j+1)`.
Place its vertices in the upper half-plane by

\[
 z_{k,j}=((j+\tfrac12)2^{-k},\;2^{-k}).
\]

Equip the half-plane with the curvature `-1` metric
`ds^2=(dx^2+dy^2)/y^2`. The standard distance formula is

\[
 \cosh d((x,y),(x',y'))
 =1+\frac{(x-x')^2+(y-y')^2}{2yy'}.
\]

The formula and the upper-half-plane description of geodesics as vertical
lines or circle arcs orthogonal to `y=0` are standard inputs: see Anna
Felikson, *Geometry III/V*, Durham University, 2023–2024, Section 6.2 and
Theorem 6.22, printed pp. 95–96 ([author's university-hosted lecture
notes](https://www.maths.dur.ac.uk/users/anna.felikson/Geometry/Geometry23/notes.pdf),
checked 2026-11-02). All construction-specific claims below are derived here.

Join each parent to each child by the hyperbolic geodesic segment. This gives a
proper noncrossing embedding of `T` into the hyperbolic plane such that:

1. every edge has the same hyperbolic length
   `L=arcosh(21/16)=log((21+sqrt(185))/16)`;
2. any two distinct vertices have hyperbolic distance at least `log 2`;
3. every compact subset of the hyperbolic plane meets only finitely many
   vertices and edges;
4. the ball of radius `rL` around the root contains all `2^r` depth-`r`
   vertices, for every integer `r >= 0`.

Thus exponential branching, a uniform vertex separation, and uniformly short
noncrossing routes coexist in a space homeomorphic to the plane. The contrast
with a fixed Euclidean lattice is metric capacity, not planar topology.

## Proof

### Constant edge length

Write a parent as `P=(a,h)`, with `h=2^(-k)`. Its children are
`Q_s=(a+s h/4,h/2)`, where `s` is `-1` or `1`. Consequently

\[
 \cosh d(P,Q_s)
 =1+\frac{h^2/16+h^2/4}{h^2}=\frac{21}{16}.
\]

### Uniform separation

For levels `k` and `l`, dropping the nonnegative horizontal term gives

\[
 \cosh d(z_{k,j},z_{l,m})
 \ge \frac{2^{k-l}+2^{l-k}}2
 =\cosh(|k-l|\log 2).
\]

Because `cosh` is increasing on the nonnegative reals, different levels have
distance at least `|k-l| log 2`. On the same level, distinct vertices have
horizontal distance at least their common height, so `cosh d >= 3/2`.
Since `cosh(log 2)=5/4 < 3/2`, all distinct vertices are separated by at
least `log 2`.

The sharp separation is actually `L`. For adjacent levels, if the coarse
index is `j` and the fine index is `m`, then

\[
 x_{k+1,m}-x_{k,j}=h\frac{2m-4j-1}{4}.
\]

The numerator is odd, so the magnitude is at least `h/4`, with equality
exactly when `m=2j` or `m=2j+1`. Therefore adjacent-level distances have
`cosh d >= 21/16`, with equality precisely on tree edges. Same-level
distances have `cosh d >= 3/2 > 21/16`; levels differing by at least two have
`cosh d >= cosh(2 log 2)=17/8 > 21/16`. Hence the minimum distance is `L`,
and precisely the parent–child pairs attain it.

### Noncrossing geodesics: the height-band step is checked explicitly

The Euclidean circle defining the edge from `P=(a,h)` to
`Q_s=(a+s h/4,h/2)` has center and squared radius

\[
 c_s=a-s\frac{11h}{8},\qquad R^2=\frac{185h^2}{64}.
\]

Indeed, at `P` the squared radius is `(11h/8)^2+h^2`, and at `Q_s`
it is `(13h/8)^2+(h/2)^2`; both equal `185h^2/64`. Both endpoints
are on the same side of `c_s` and strictly above `y=0`. On that side the
upper semicircle `y=sqrt(R^2-(x-c_s)^2)` is strictly monotone. Its segment
between the two endpoints therefore stays in their horizontal interval and
in the height band `[h/2,h]`; the band boundary heights are attained only
at the respective endpoints. In particular, no unverified assumption that
an arbitrary hyperbolic segment stays below both endpoint heights is used.

At a fixed level, the two outgoing edge intervals for one parent meet only
at `x=a`, where both arcs meet only at the parent. The union of those
intervals is `[a-h/4,a+h/4]`. Distinct parents have horizontal centers
separated by at least `h`, so their outgoing intervals are disjoint.
Edges starting at different levels have disjoint open height bands. Only
consecutive bands share a boundary height, where the arcs have only their
endpoints. Distinct vertex coordinates are distinct, so such a common
point is a genuine incident vertex. Thus edges intersect only at their
shared endpoints, and no edge interior contains another vertex.

### Properness, exponential capacity, and the Euclidean boundary

A compact set `K` in the upper half-plane has a strictly positive minimum
height `epsilon`. Vertices in `K` must satisfy `2^(-k) >= epsilon`, so
there are only finitely many possible levels and vertices. An edge meeting
`K` has top height `2^(-k) >= epsilon`, giving only finitely many possible
edges as well. Together with the injective geodesic edge maps and their
intersection description, this proves a proper topological embedding of the
locally finite tree. The triangle inequality along a depth-`r` tree path
puts each of its `2^r` endpoints within hyperbolic distance `rL` of the root.

The apparent Euclidean crowding is not hyperbolic crowding. Within the open
half-plane there are no accumulation points. In the closed Euclidean
half-plane, the accumulation set of the vertices is exactly `[0,1] x {0}`:
all coordinates have `0 < x < 1`, unbounded levels have heights tending to
zero, and the depth-`k` horizontal coordinates are the midpoints of the
`2^k` dyadic subintervals of `[0,1]`. Every point of `[0,1]` is within
`2^(-k-1)` of one of those midpoints. The boundary `y=0` is not part of the
hyperbolic plane, so this does not violate properness there.

## Negative boundaries and connection contract

- **Not Euclidean uniform separation.** Same-level neighboring Euclidean
  distances are `2^(-k)`, tending to zero. Merely drawing the points in a
  Euclidean coordinate chart does not change the metric used in the theorem.
- **Not a quasi-isometric tree embedding.** For every `k >= 1`, the two
  depth-`k` vertices with indices `2^(k-1)-1` and `2^(k-1)` lie in different
  root-child subtrees, so their graph distance is `2k`; their hyperbolic
  distance is the constant `arcosh(3/2)`. There is no coarse inverse bound
  on graph distance from hyperbolic distance for this placement.
- **Not a Collatz completion.** The source is the abstract binary tree.
  A separate proven macro-address correspondence could transport this
  drawing to that selected branching structure only. No map from every
  Collatz state, no routing of every one-step edge, no new descent invariant,
  and no termination theorem are provided here.
- **Not a regular hexagonal mesh.** No tessellation, mesh adjacency, hexagon
  shape, or mesh-routing theorem is constructed. A hyperbolic hexagonal
  model requires its own construction and audit.

The map preserves vertex distinction, binary adjacency via the drawn edges,
planarity, uniform separation in the stated metric, and bounded route
length. It does not preserve the tree's large-scale intrinsic distance or
Euclidean lattice occupancy. The indispensable sidecar is the ambient
hyperbolic metric together with the address labels. The cheapest hostile
tests are Euclidean same-level crowding and the fixed-hyperbolic-distance
pairs across the root split; both are included in `verify.py`.

## Reproduction and scope of checks

Run `python verify.py` and `python -O verify.py` from this directory. The
checker uses rational arithmetic on all vertices through depth seven,
including every vertex pair and every edge pair, and raises explicit errors
rather than relying on optimization-sensitive `assert` statements. It checks
the edge constants, the sharp separation/equality classification, the circle
equations, monotone-side signs, the noncrossing band/interval certificates,
and the two negative controls. These finite checks do not replace the
all-level arguments above.

`draw.py` produces the small illustrative SVG and PNG from the stated exact
dyadic coordinates and circle centers. The plot is a Euclidean coordinate
view of the hyperbolic model, not a Euclidean equal-edge drawing and not a
hexagonal lattice. Only levels zero through four are displayed.
