# An exact twelve-edge optimum in a specified theta class

**Status: PROVED, using a CITED classical fixed-side polygon theorem.** This is a scoped optimality result for a known Friedman configuration, not a global product optimum and not a novelty claim about polygonal isoperimetry.

## Statement and boundary

Consider legal endpoint-only unit-edge product configurations whose graph is a theta with three internally vertex-disjoint paths of lengths 5,5,2 between distinct terminals u,v. Require that:

1. The two five-edge paths form the outer simple ten-edge polygon, and the two-edge path is the common separator of its two bounded fields.
2. The open straight segment uv lies in the interior of that outer polygon. This holds, in particular, when the outer polygon is convex. The separator itself is initially allowed to bend.

Then the maximum product of the two field areas is exactly

    (89+26 sqrt(13))/16.

It is attained by the construction in [PROOF.md, Section 4](../PROOF.md). Equality forces a straight two-edge separator and the reflected cyclic hexagons described there, up to Euclidean isometry and interchange of the two fields.

The second hypothesis is explicit. A straight chord between vertices of a nonconvex polygon can leave the polygon or intersect its boundary. That case is not silently covered by the following two-simple-polygons argument. Nor are the other planar theta embeddings, in which a five-edge path is the common separator, covered.

## Fixed-side cap

Let g(l) be the largest area of a simple hexagon with cyclic side-length list

    (1,1,1,1,1,l),  0<l<=2.

The classical secant-polygon theorem makes the maximum the unique convex cyclic polygon with those side lengths. The exact primary source used here is Kouřimská, Skuppin, and Springborn, *A Variational Principle for Cyclic Polygons with Prescribed Edge Lengths* (2016), [Theorems 1.1 and 1.5](https://link.springer.com/chapter/10.1007/978-3-662-50447-5_5). The first gives existence/uniqueness under the strict polygon inequalities, satisfied here; the second gives maximal area for prescribed side lengths. Both were directly inspected on 2026-11-03. This is classical input, not a new theorem of this session.

Give each unit side central half-angle alpha. The remaining side has half-angle beta=pi-5 alpha. Set

    l(alpha)=sin(5 alpha)/sin(alpha),
    g(l(alpha))=[5 sin(2 alpha)-sin(10 alpha)]/[8 sin(alpha)^2].

On pi/10<alpha<pi/5, beta lies in (0,pi/2), and l decreases strictly from 1+sqrt(5)>2 to zero. Indeed

    l'(alpha) = [5 cos(5 alpha) sin(alpha)
                 - sin(5 alpha) cos(alpha)] / sin(alpha)^2 < 0.

Consequently every 0<l<=2 has exactly one alpha in this interval, and these central angles construct the unique cyclic hexagon. Differentiating its area gives

    dg/dl = -cos(5 alpha)/(2 sin(alpha)) > 0.

For an entirely algebraic derivative check, put z=sin(alpha)^2 and

    L(z)=16z^2-20z+5,
    C(z)=16z^2-12z+1,
    B(z)=5-L(z)C(z)=-256z^4+512z^3-336z^2+80z.

Then g=sqrt((1-z)/z) B(z)/4, l=L(z), and the derivative identity reduces to the polynomial identity

    2z(1-z) B'(z) - B(z) + 4z(1-z) C(z) L'(z) = 0.

The checker verifies every coefficient exactly. Thus g(l)<=g(2), with strict inequality for l<2.

## From the cap to the two-field product

The terminal distance l=|u-v| is at most two by the triangle inequality along the shared two-unit-edge path. The interior-chord hypothesis makes each outer five-edge path together with uv a simple hexagon. Their areas sum to the area of the outer ten-edge polygon, independently of the position of the common separator. Hence, for the two field areas A_1,A_2,

    A_1+A_2 <= 2g(l) <= 2g(2),
    A_1 A_2 <= ((A_1+A_2)/2)^2 <= g(2)^2.

At l=2 the equalities along the two-edge path force its two edges to be collinear in the same direction. The equations l(alpha)=2 and z=sin(alpha)^2 give

    16z^2-20z+3=0,  z=(5-sqrt(13))/8.

Substitution, or the independently verified shoelace coordinates in PROOF.md, gives

    g(2)^2=(89+26 sqrt(13))/16.

The known reflected construction attains every inequality. Conversely equality requires l=2, both chord-closed hexagons to be the fixed-side cyclic maximizers, and equal field areas. The unique cyclic maximizers lie on opposite sides of uv and are reflections of each other. This proves the statement and its equality case.

## What changed, and what did not

The cactus theorem says shared cycle edges are necessary to exceed the restricted cycle-partition maximum. This result goes one step further: within the stated theta embedding class, the **metric terminal distance** controls the optimum and forces the shared path to straighten.

Abstract theta incidence or cycle rank alone is insufficient. Retain the path lengths, which path separates the two fields, and the straight terminal chord's position relative to the outer polygon. General twelve-edge product optimization and nonconvex terminal-chord placements remain outside this theorem.
