# A boundary gap repairs the hyperbolic tree metric

**PROVED elementary extension, same-agent audit.** This is an explicit
all-depth construction, with finite rational controls. No external priority
claim is made for Cantor-boundary tree embeddings.

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

## Scope

The source here is the abstract binary tree. Transport to a selected Collatz
inverse subtree requires its actual source-address map and valuation clocks.
This is not a realization of every ordinary Collatz edge, a regular hyperbolic
honeycomb, or a convergence proof. Trees remain two-colourable despite their
exponential growth; this construction supplies no higher chromatic lower bound.

The exact checker covers all pairs through depth6 for b=2,9/4,5/2,3,4,5:
48,006 pairs, with rational cosh comparisons, circle equations and the
critical b=2 hostile. The all-real, all-depth claim is the proof above.
