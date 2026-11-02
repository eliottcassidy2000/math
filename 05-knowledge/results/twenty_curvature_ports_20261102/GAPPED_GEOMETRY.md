# Gaps preserve binary tree distances in the hyperbolic plane

**Status: PROVED** elementary all-depth construction and bounds below;
**FINITE-EXACT** rational corroboration at 24,003 vertex pairs.
Hyperbolic tree embeddings are classical; no priority claim is made.
The contribution here is an explicit repair of the recovered dyadic
placement, with source-labelled Collatz distance bounds.

**Concurrent credit:** incoming Session141 at9e1ae5aaf7 independently
supplies the same interval-gap mechanism in a more general two-parameter
family, including the sharp edge-length infimum log2. See its
[curvature proof, Section5](../session141_metric_aliquot_square_paths_20261102/curvature_triangle/PROOF.md).
This packet's width-one case additionally retains the induced
minimum-distance property and proves the native metric and LCA-detour
bounds below. The shared gap construction is not claimed as exclusive
new work here.

## G1. The interval construction

Fix a real b>=2 and put r=1/b. Give the root the interval [0,1].
Inside every interval retain the leftmost and rightmost subintervals of
relative length r. A word w in {0,1}^k selects one such interval I_w.
Place its vertex at (midpoint(I_w), length(I_w)) in the upper half-plane
with metric ds^2=(dx^2+dy^2)/y^2.

Use the standard formula

    cosh d((x,y),(X,Y)) = 1 + ((x-X)^2+(y-Y)^2)/(2yY).

The inherited [curvature control](../session140_collatz_hexagonal_geometry_20261102/curvature_control/PROOF.md)
gives its primary-source reference and proves the b=2 case. General
hyperbolic tree embedding context is supplied by
[Sala, De Sa, Gu and Re, 2018](https://proceedings.mlr.press/v80/sala18a.html).

Write a parent as (a,h); its children are
(a +/- (1-r)h/2, rh). Every geodesic parent-child edge has length L_b where

    cosh L_b = 1 + 5(1-r)^2/(8r) = (5b^2-2b+5)/(8b).      (G1)

The geodesic circle for the right child has centre a-Ah on y=0, with
A=(3+5r)/4, B=(5+3r)/4 and

    R^2/h^2 = A^2+1 = B^2+r^2.

Both endpoints are on the right half of this circle, so the segment is
monotone in both x and y and stays in its endpoint horizontal interval
and height band [rh,h]. The left-child statement is its reflection.
Outgoing horizontal intervals of different same-level parents are disjoint;
different levels have disjoint open height bands. The resulting geodesic
embedding is noncrossing. Every compact set has a positive minimum height,
so it meets only finitely many levels and edges: the embedding is proper.

All distinct vertices are at least L_b apart, with equality exactly on
parent-child pairs. To see this, same-level centres are separated by at
least (b-1) times their height; therefore cosh d>=1+(b-1)^2/2>cosh L_b.
On adjacent levels the closest coarser point is the actual parent, with
horizontal offset (1-r)h/2; all other coarser points are strictly farther.
For a level difference at least two, cosh d>=cosh(2 log b)>cosh L_b,
because the positive numerator of the latter difference is

    4b^4-5b^3+2b^2-5b+4 = (b-1)^2(4b^2+3b+4).

This is an induced fixed-distance tree in the stated hyperbolic metric.

## G2. The missing gap is exactly the coarse-distance obstruction

Let D=d_T(u,v) be ordinary unweighted binary-tree distance.
For b>2 define g=1-2/b>0. Then for every pair,

    D log b + 2 log g <= d_H(u,v) <= L_b D.                 (G2)

The upper bound is the triangle inequality along the tree path. For the
lower bound, let q be the address LCA at depth l. If u,v lie in opposite
child regions of q, their horizontal distance is at least g b^-l.
Their heights are b^-k and b^-m, with D=k+m-2l. The distance formula gives

    2 sinh(d_H/2) >= g b^(D/2).

Since exp(d_H/2)>=2 sinh(d_H/2), (G2) follows. If one vertex is the LCA,
the vertical ratio already gives d_H>=D log b. Thus the construction is
a quasi-isometric embedding for every b>2.

The equality/failure boundary b=2 matters. For the words 01...1 and
10...0 at depth h, the horizontal gap divided by their common height is
g b^h+1. Exactly,

    cosh d_H = 1 + (g b^h+1)^2/2,   d_T=2h.               (G3)

At b=2 this is the constant 3/2 at every depth, so there is no coarse
inverse bound. At b=3 it is 1+(3^(h-1)+1)^2/2. In particular the simple
ternary construction has

    (D-2) log 3 <= d_H <= arcosh(11/6) D.

Negative curvature supplies room, but a specific drawing must still keep
different branches separated at the boundary. The dyadic counterexample
has sufficient volume and constant edge lengths while losing that property.

## G3. Triangle inequality and the LCA detour

There is a uniform upper bound on the excess caused by routing through
the LCA in the gapped construction:

    0 <= d_H(u,q)+d_H(q,v)-d_H(u,v) <= 2 log(5/(4g)).       (G4)

For an ancestor q and descendant u at depth difference a, the horizontal
offset is at most (1-b^-a)/2 times q's height. Put t=b^-a. Then

    cosh d_H(u,q) <= 5/(8t)+5t/8-1/4
                  <= cosh(a log b + log(5/4)).

The difference on the last line's right minus left is
1/4-9t/40>0 for 0<t<=1. Along with the vertical lower bound this gives

    a log b <= d_H(u,q) <= a log b + log(5/4).

Apply this to the two arms, then subtract G2. The left inequality of G4
is the ordinary triangle inequality. At b=3 the upper bound is
2 log(15/4); at b=2 no uniform bound is possible by G3.

The same ancestor estimate sharpens G2's upper bound to

    d_H(u,v) <= min{L_b D, D log b + 2 log(5/4)}.           (G4a)

This explicit additive upper bound is another direct use of triangle
inequality. The concurrent two-parameter proof gives a general additive
bound; G4a uses the narrower width-one ancestor calculation.

An exactly isometric copy of a positive-edge-length branching tree cannot
sit in H^2. For a degree-three vertex, tree distances make the triangle
inequality an equality for each pair of distinct outgoing neighbours
through the centre. Equality in H^2 forces the centre onto their common
geodesic segment. Three distinct rays cannot all be pairwise opposite.
Quasi-isometry accommodates a bounded detour; it does not violate or
remove triangle inequality.

## G4. Transport to actual selected Collatz sources and clocks

Use ordinary C(n)=n/2 for even n and C(n)=3n+1 for odd n. The inherited
SF5 binary family starts at 5. For an odd n coprime to 3 select the two
children

    u=(2^a n-1)/3,
    a in {2,4,6} if n=1 mod3, or {3,5,7} if n=2 mod3,
    discarding the unique a with 2^a n-1=0 mod9.

Each selected source is odd, coprime to 3, larger than its parent, and
returns to it in exactly a+1 ordinary steps. The source address, integer
label and exponent a are retained. This is the inherited guarded family,
not a claim that every positive integer has such an address.

Let d_C denote distance in the actual ordinary inverse tree rooted at 5,
restricted to these selected odd sources. For two different branches at
the odd LCA q, the actual meeting is 2^min(a,b) q, before reaching q.
Consequently the two adjacent macro edges cost |a-b|+2 in the native
tree, not a+b+2. The exponents are distinct and have the same parity,
so this cost is at least four. All other selected edges have clocks
between three and eight. Thus for all selected source pairs

    3D-2 <= d_C <= 8D.                                   (G5)

For ancestor pairs the stronger lower bound 3D holds. No shortcut through
another ordinary vertex can improve this distance: the inverse basin
stopped at 5 is a tree. For example 13 and 53 meet at 40, with clocks
1 and 3, so their native distance is 4 although their odd LCA is 5.

Composing G5 with the b=3 placement yields an explicit quasi-isometric
map of the selected source metric:

    (d_C/8-2) log3 <= d_H <= (arcosh(11/6)/3)(d_C+2).       (G6)

G4a also bounds the right side by ((d_C+2)/3)log3+2log(5/4).

The finite audit compares 8,001 source pairs by literal ordinary trajectories,
including earlier even meetings. The proof of G5 is the clock calculation
above, not extrapolation from that sample.

## Arithmetic constants and the geometric limit of the connection

The recovered pairs 11 and 13 occur in the b=2 circle exactly as

    11^2+8^2 = 13^2+4^2 = 185.

Here 11/8 and 13/8 are the horizontal radii at the parent and child,
whose heights have ratio 1:1/2. These are the concrete 12-1 and 12+1.
At b=3 the analogous identity is 7^2+6^2=9^2+2^2=85.

Vertical dilation by 2^k has distance k log2. Its Lorentz coordinates are

    cosh(k log2)=(4^k+1)/2^(k+1),
    sinh(k log2)=(4^k-1)/2^(k+1).

Thus 3,4,5 and 15,8,17 arise at k=1,2 by the same composition law.
The boost (1/4)[[5,3],[3,5]] squares to (1/8)[[17,15],[15,17]];
it preserves X^2-Y^2. This is a genuine map and invariant, rather than
an inference from the repeated appearance of 3,5 or 16+/-1.

The metric remains hyperbolic. A conformal Euclidean picture does not
preserve lengths and cannot transfer a unit-distance chromatic bound.
No regular hyperbolic honeycomb mesh or embedding of every ordinary
Collatz state has been constructed here. The all-source convergence
problem retains its missing coverage and descent obligations.
