# A boundary gap repairs the hyperbolic tree metric

**PROVED elementary extension, same-agent audit.** This is an explicit
all-depth construction, with finite rational controls. No external priority
claim is made for Cantor-boundary tree embeddings. Incoming Session141
independently supplies the broader two-parameter gap construction; the exact
nearest-distance phase classification below compares that construction with
our W=1 sharp-separation proof.

## Inherited model and the missing property

[Session140's curvature control](../session140_collatz_hexagonal_geometry_20261102/curvature_control/PROOF.md)
places a binary tree at dyadic interval centres and heights 2^-k. It is
proper, noncrossing and uniformly separated in the hyperbolic metric;
its edges have length arcosh(21/16). But two vertices across the root split
can have tree distance 2k and hyperbolic distance arcosh(3/2), independently
of k. The inherited theorem explicitly excludes a quasi-isometry.

We retain its interval-centre construction and introduce a positive gap.
The standard metric input is the upper-half-plane formula

    cosh d((x,y),(x',y')) = 1+[(x-x')^2+(y-y')^2]/(2yy'),       (1)

with geodesics vertical lines or circles orthogonal to y=0. See the source
already used by Session140, [Felikson's Geometry III/V notes, section6.2](https://www.maths.dur.ac.uk/users/anna.felikson/Geometry/Geometry23/notes.pdf).
The estimates specific to this construction are proved below.

## The all-depth theorem

Fix a REAL b>2. For a finite binary word w=e1...ek put

    a_w = sum_(r=1..k) (b-1)e_r b^-r,
    I_w = [a_w,a_w+b^-k],
    z_w = (a_w+b^-k/2, b^-k).

Include the empty word, at (1/2,1). The children occupy the leftmost and
rightmost subintervals of relative length 1/b. Their intervening gap has
relative length delta=(b-2)/b>0.

If D(u,v) is binary-tree distance, then for every pair of finite words,

    D(u,v) log b - 2 log[b/(b-2)]
          <= d(z_u,z_v) <= D(u,v) log b + log 4.             (2)

Thus, after the fixed scale log b, the embedding preserves tree distance
up to a uniform additive error. In particular it is a quasi-isometric
embedding, strengthening the inherited b=2 construction. The rough isometry
is onto its image; coarse surjectivity onto the hyperbolic plane is neither
asserted nor needed.

Moreover it is proper and noncrossing, and its exact minimum pair distance is

    L_b = arcosh[1+5(b-1)^2/(8b)],                          (3)

attained precisely on parent-child pairs. Thus its distance-L_b graph is
exactly the binary tree. Rescaling the hyperbolic metric makes L_b unit
length; this does not turn the geometry into the Euclidean plane.

## Proof of the distance estimate

Let k=|u|, l=|v|, and h be the length of their longest common prefix.
Then D=k+l-2h. Both points lie horizontally inside I_(common prefix),
whose length is b^-h; both heights are at most b^-h. Formula (1) gives

    cosh d <= 1+b^D.

Using arcosh X<=log(2X), and b^D>=1, gives the upper bound in (2).

If one word is an ancestor of the other, discarding the horizontal term
in (1) gives d>=|k-l|log b=Dlog b. Otherwise the two words occupy different
children of their common prefix. Their horizontal separation is at least
delta*b^-h. Write E for their Euclidean coordinate distance. Equation (1)
is equivalently

    d=2 asinh[E/(2 sqrt(y_u y_v))].

Since asinh t>=log(2t) for t>0,

    d >= log[E^2/(y_u y_v)] >= D log b+2 log delta.

This proves the lower bound, including the cases where it is negative.
The proof uses the full LCA and the gap, not just a level count.

## Edges, separation and noncrossing

A parent (x,h) has children (x plus/minus alpha*h,h/b), where
alpha=(b-1)/(2b). Substitution in (1) gives (3).

At a fixed level, distinct centres are separated horizontally by at least
(b-1) times their common height. To see this, if their last common prefix
has length h<k, the closest such centres have horizontal separation

    (b-2)b^(-h-1)+b^-k >= (b-1)b^-k.

Their cosh distance is therefore at least 1+(b-1)^2/2, greater than (3)'s
cosh value. Between adjacent levels, the child of another parent has
horizontal offset at least [(b-1)-alpha]h>alpha*h; only the two genuine
children attain (3). For levels differing by at least two, vertical
separation gives cosh d>=cosh(2 log b), and

    cosh(2 log b)-cosh L_b
      = (b-1)^2(4b^2+3b+4)/(8b^2)>0.

Thus (3) is sharp and only tree edges attain it.

For the child on side s=+1 or -1 the geodesic circle has centre

    c=x-s*gamma*h,       gamma=(3b+5)/(4b),

and squared radius h^2(gamma^2+1). The identity
gamma^2+1=(gamma+alpha)^2+b^-2 checks both endpoints. They are on the
same side of c, so the arc is monotone and stays in its horizontal interval
and height band [h/b,h]. At one level, different parents' outgoing intervals
are disjoint; across levels, open height bands are disjoint. Shared boundary
points are exactly shared tree vertices. This proves noncrossing.

A compact hyperbolic set has height bounded below by a positive number,
so only finitely many levels and edges meet it. This proves properness.

## The triangle inequality and the boundary

Negative curvature does NOT violate or weaken the triangle inequality.
It supplies more uniformly separated locations inside a metric ball.
The inequality still gives d(z_root,z_w)<=|w|L_b. What was missing at b=2
was a reverse estimate; equation (2) supplies one when a boundary gap exists.

For infinite binary addresses define rho(u,v)=b^(-h), where h is their
common-prefix length. Then

    rho(u,w)<=max(rho(u,v),rho(v,w)).                       (4)

Indeed the first and third addresses agree for at least the shorter of
the other two agreement lengths. This is the ultrametric triangle inequality.
The limit map to the real boundary is

    x(e)=sum_(r>=1) (b-1)e_r b^-r,
    delta*rho(u,v)<=|x(u)-x(v)|<=rho(u,v).                 (5)

The lower bound is precisely the child gap. At b=2 it disappears, and
the infinite addresses 0111... and 1000... have the same boundary point.
Their finite truncations are the inherited constant-distance hostile.
For every b>2 the boundary map is injective with the explicit estimates (5).

The level supply is 2^k while the total length of the 2^k boundary intervals
is (2/b)^k. This identifies a critical b=2 regime and a separated b>2 regime
without confusing binary source clocks with the ambient metric.

## Incoming construction and the exact nearest-distance phase boundary

The incoming [Session141 curvature theorem, section5](../session141_metric_aliquot_square_paths_20261102/curvature_triangle/PROOF.md)
uses q=1/b in (0,1/2) and horizontal width 0<W<=1: multiply every horizontal
coordinate above by W and leave the height q^k unchanged. Its proof gives
the same lower estimate with gap W(1-2q), and the stronger upper constant
log(W^2+2). Our displayed log4 bound is valid but not sharp. Its packing
argument proves the necessary edge-scale bound L>=log2 in curvature -1;
all strict values are realized by the two-parameter family, while attainment
at equality is open. Those incoming results are credited, not claimed anew.

Here is a further COMPLETE classification within this parameter family.
Write E for the parent-child distance and S for the sibling distance:

    cosh E = 1+(1+W^2/4)(1-q)^2/(2q),
    cosh S = 1+W^2(1-q)^2/(2q^2),
    Wcrit(q)^2 = 4q/(4-q).                              (6)

The exact minimum distance between distinct vertices is min(E,S), and the
graph of pairs at that minimum distance has three regimes:

- W>Wcrit: exactly the binary tree;
- W=Wcrit: the binary tree together with every sibling edge;
- W<Wcrit: exactly the disjoint sibling pairs, with the root isolated.

**Proof.** At one level the minimum horizontal separation is
W(1-q)q^(k-1), attained exactly by siblings, since q<1/2. Between adjacent
levels the true children have offset W(1-q)q^k/2; a child of another parent
has offset at least W(1-q)(1/q-1/2)q^k, strictly larger. Thus E is the
adjacent-level minimum, attained exactly by the tree edges. Levels differing
by at least two have distance >=2log(1/q), which exceeds E: at fixed q the
edge length increases with W, so the W=1 strict inequality proved above
applies to all W<=1. These three cases prove the minimum and its entire
equality graph. Comparing the two cosh expressions gives

    E<S iff W^2(1-q/4)>q,

and equality or reversal gives the other regimes. QED.

At the boundary, each parent and its two children is an equilateral
hyperbolic triangle at the minimum distance. This is an actual three-colour
graph, not a Euclidean Hadwiger-Nelson lower bound. Above the boundary the
minimum-distance graph needs two colours; below it there are only matching
edges and the isolated root. All three placements still satisfy the same
positive-gap coarse tree-metric estimates: metric recovery alone does not
determine the nearest-distance relation.

A pair closer than E is not automatically an edge of the distance-E graph.
The classification concerns the MINIMUM-distance relation. When S<E it does
not classify every pair at the parent-child length E, and supplies no lower
bound for arbitrary induced unit-distance realizations.

**Sharp edge-scale tradeoff in this family.** Requiring the minimum-distance
graph to be exactly the tree forces

    E > arcosh(9/7) > log2.                             (7)

The infimum arcosh(9/7) is sharp WITHIN this two-parameter construction,
not asserted as a universal lower bound for all hyperbolic tree drawings.
Indeed at Wcrit the edge cosh is

    1+2(1-q)^2/[q(4-q)],

and subtracting 9/7 gives

    2(2q-1)(4q-7)/[7q(4-q)] > 0.

For W>Wcrit the edge length is larger still. Let q approach1/2 from below
and W approach Wcrit(q) from above to attain the infimum as a limit;
Wcrit(q)<1 leaves room for this choice. The bound is not attained inside
the stated strict parameter range. This explains precisely why the incoming
near-log2 construction loses exact minimum-distance adjacency, and why our
W=1 construction retains it. The triangle inequality is unchanged in every
regime.

## Scope

The source here is the abstract binary tree. Transport to a selected Collatz
inverse subtree requires its actual source-address map and valuation clocks.
This is not a realization of every ordinary Collatz edge, a regular hyperbolic
honeycomb, or a convergence proof. Trees remain two-colourable despite their
exponential growth; this construction supplies no higher chromatic lower bound.

The exact checker covers all pairs through depth6 for b=2,9/4,5/2,3,4,5:
48,006 pairs, with rational cosh comparisons, circle equations and the
critical b=2 hostile. The all-real, all-depth claim is the proof above.
The incoming comparison adds23 rational (q,W) configurations through depth5,
including three EXACT equality cases, checking44,919 pairs and each full
minimum-distance graph. These finite checks support, not replace, (6)--(7).
