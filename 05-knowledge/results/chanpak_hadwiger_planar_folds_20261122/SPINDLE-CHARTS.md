# Unit-distance geometry and polynomial charts on physical tournament flips

**Proved.** There is a labelled strong 7-vertex tournament cube with twelve literal Hamilton-gap-two edges and an exact planar unit-distance realization. A non-affine polynomial automorphism over the rationals deforms its nonedge lengths while preserving all twelve seam lengths. A finite sufficiently replicated seam test forces bounded-degree polynomial maps to be affine isometries. Separately, the strong order-six **isomorphism quotient** of all one-arc flips contains an induced Moser spindle: all eleven edges are genuine one-arc adjacencies after relabelling, but only **four** have Hamilton count gap two. These operations do not improve the chromatic bounds for the plane or settle the planar Jacobian conjecture.

## 1. A literal unit cube with a lossy unmarked read

Use outgoing masks \((80,73,35,5,14,27,60)\), on vertices numbered 0 through 6. Independently reverse the three named arcs at unordered pairs \((0,1),(2,6),(3,4)\), indexed by the bits \(b=(b_0,b_1,b_2)\). Let \(Z_T(2)\) count ordered spanning covers by two possibly empty paths. Exact subset-endpoint DP gives
\[
(H(T_b),Z_{T_b}(2))=(141,1656)+(2,4)(b_0+b_1+b_2\bmod2).
\]
All eight tournaments are strong. Their labelled one-arc graph is the cube \(Q_3\), since the Hamming distance of the masks is the number of differing arcs. Every edge has absolute Hamilton count gap two. They are only **two unmarked isomorphism classes**, one per parity. Both the labels of the switched arcs and the distinction between geometric unit distance and Hamilton height are load-bearing.

There is no injective planar unit-distance realization of this cube in which \(H\) is an affine coordinate: the height layers would lie on two parallel lines, and the unit circle about a point on either line meets the other line in at most two points. Every vertex needs three opposite-layer neighbors. Using \((H,Z_T(2))\) as the source coordinates is worse: these coincide four-to-one, and a polynomial automorphism cannot separate a collision already present in its input.

## 2. Exact inverse and a genuinely non-affine deformation

Take the rational unit vectors
\[
u=(1,0),\quad v=(0,1),\quad w=(3/5,4/5),\quad w'=(5/13,12/13).
\]
The points \(P_b=b_0u+b_1v+b_2w\) and \(P'_b=b_0u+b_1v+b_2w'\) each realize exactly \(Q_3\). More generally with \(w=(a,b)\) in the open first quadrant of the unit circle and neither coordinate \(1/2\), a two-term signed sum has squared norm 2 or \(2\pm2a\) or \(2\pm2b\), none 1 or 0. A three-term sum has squared norm \(3\pm2a\pm2b\), never 1 since \(a+b>1\) and \(|a-b|<1\), and never 0 since \(w\ne(\pm1,\pm1)\). Differences of one chosen bit alone have length one.

The determinant-one linear coordinates
\[
L(x,y)=(4x+7y,y/4)
\]
give the same first coordinate \(m_b=4b_0+7b_1+8b_2\) in both frameworks, distinct for all eight masks. Let \(R\) be the Lagrange polynomial solving \(R(m_b)=b_2\). Its leading coefficient is \(17/300960\), so its degree is exactly seven. In the coordinates of \(L\), the vertical shear
\[
(m,z)\longmapsto(m,z+(2/65)R(m))
\]
sends \(L(P_b)\) to \(L(P'_b)\); its inverse subtracts the same term, and its Jacobian determinant is exactly one. Conjugate by \(L\) to get the polynomial automorphism of the physical unit framework. It is not a Euclidean isometry: the squared length between \(u\) and \(w\) changes from \(4/5\) to \(16/13\) while every seam length stays one.

There is also a det-one polynomial decoder from the natural **address/height** records \((m_b,\epsilon_b)\), \(\epsilon_b=b_0+b_1+b_2\bmod2\). Interpolate \(C(m_b)=(b_1+(4/5)b_2)/4-\epsilon_b\) and use the chart
\[
(m,e)\longmapsto L^{-1}(m,e+C(m)).
\]
The inverse applies \(L\) and subtracts \(C(m)\). It recovers \(H=141+2e\) and all three marked arc bits by the eight-entry lookup for \(m\). This chart transports the abstract edge relation to exact Euclidean edges, not initial Euclidean distances between the address records.

For an arbitrary cube \(Q_d\), one can also obtain an induced rational unit-distance embedding by successively adjoining rational unit vectors, excluding the finitely many choices that make any already nonzero signed sum plus the new vector have length zero or one. This applies directly to induced labelled physical families in which all Hamming-one edges are unit seams, but not to a selective subgraph with some Hamming-one edges deleted.

## 3. A finite physical rigidity test for bounded degree

Fix \(D\ge1\). Place translated copies of the cube framework on a spaced rectangular grid with \(2D-1\) distinct translation coordinates in each axis (for example spacing 10). The physical tournaments in each copy are obtained by attaching a number of common ordinal singleton sinks that is different for each copy. This preserves both the counts and the marked flips; distinct copies have different orders, so there are no new physical cross-copy adjacencies. Spacing keeps all geometric cross-copy distances greater than one.

Suppose a two-variable polynomial map \(F\) of total degree at most \(D\) preserves the squared unit length of the three edges out of the base point in each copy, in directions \(u,v,w\) from section 2. For each direction \(a\) the polynomial
\[
\big(F_1(x+a)-F_1(x)\big)^2+\big(F_2(x+a)-F_2(x)\big)^2-1
\]
has degree at most \(2(D-1)\) and vanishes throughout that rectangular grid. Applying the ordinary univariate root bound twice makes it identically zero. Even over the complex numbers, factoring the sum of squares as \((g_1+i g_2)(g_1-i g_2)=1\) in the polynomial ring forces both factors, hence both differences, to be constant. Two independent directions now force \(F\) affine: after subtracting its corresponding linear difference map, both coordinates of the residual are polynomials periodic in two independent directions, hence constant. The linear part has squared length one on three distinct real unit lines; a binary homogeneous quadratic vanishing on three directions is zero, so its matrix \(A\) satisfies \(A^{\mathsf T}A=I\). This is an isometry over the reals; over the complex numbers it preserves the algebraic squared-distance form. If the given Jacobian determinant is 1, it is orientation-preserving in the real case. The finite non-affine deformation above is compatible with this degree-bounded replicated rigidity.

## 4. An induced Moser spindle in the strong order-six quotient

Keep one labelled canonical mask for each strong six-vertex tournament isomorphism class. Bits are the lexicographic pairs \(i<j\), with bit 1 meaning \(i\to j\). A quotient edge means that some (equivalently, the selected) representative of the first class has a one-arc flip to some representative of the other.

Use the ordered mask list \((20,16,22,26,81,156,89)\), with Hamilton counts \((19,17,15,19,23,33,25)\). All are strong and pairwise nonisomorphic. Index them \(0,\dots,6\). The complete adjacency and labelled certificates follow. Each witness is a flip from the first vertex's canonical mask to the listed mask in the second vertex's orbit.

| Quotient edge | Flipped target mask | Absolute Hamilton gap |
|---|---:|---:|
| 01 | 16 | 2 |
| 02 | 22 | 4 |
| 04 | 84 | 4 |
| 05 | 276 | 14 |
| 12 | 8208 | 2 |
| 13 | 2064 | 2 |
| 23 | 4118 | 4 |
| 36 | 90 | 6 |
| 45 | 337 | 10 |
| 46 | 89 | 2 |
| 56 | 4252 | 8 |

No other quotient edges occur: take every single flip of each listed canonical mask and test membership in the other six explicit relabelling orbits. This also verifies the eleven witnesses. An exhaustive scout can canonically enumerate all 56 six-vertex classes then backtrack on the seven-vertex incidence pattern, matching edge and nonedge constraints at every added role. Among strong induced copies at order six this finds a maximum of four gap-two edges. Order five has only six strong classes, so order six is minimal for this strong seven-role fixture.

For a fully geometric replay use two unit diamonds sharing \(O=(0,0)\). The first has arms \(P=(\sqrt3/2,1/2)\), \(Q=(\sqrt3/2,-1/2)\) and tip \(T=(\sqrt3,0)\). Rotate the three noncentral points about \(O\) by the rotation with cosine \(5/6\) and sine \(\sqrt{11}/6\) to obtain \(P',Q',T'\). Assign roles \(O,P,Q,T,P',Q',T'\) in order. The squared distances in \(\mathbb Q(\sqrt3,\sqrt{11})\) give exactly the listed eleven unit edges. In each diamond, a proper 3-coloring forces its central point and tip to have the same color; the unit edge between tips defeats this, and a fourth color on one tip suffices. Thus the quotient graph's actual selected adjacencies form a 4-chromatic unit-distance graph. They are **not** eleven Hamilton-gap-two seams.

This distinction is necessary, not cosmetic. On every pure gap-two graph the phase \((H-1)/2\bmod2\) is a 2-coloring. Every odd cycle of quotient arc reversals therefore needs an edge whose gap is 0 modulo 4. The first triangle here records label holonomy explicitly: the labelled path \(20\to16\to8208\) needs the odd relabelling \((1,2,3,4,5,0)\) to mask 22 before the literal arc \(22\to20\) closes it; the labelled flip cube itself has no triangle.

Placing the quotient geometry into a polynomial chart creates no additional obstruction. For instance start at records \((\text{mask},H)\), whose masks are distinct, and shear the target linearly by \((x,y)\mapsto(x,y+x)\), whose seven resulting second coordinates are distinct (checked exactly). A vertical Lagrange shear sends each starting height to its desired second coordinate; a horizontal Lagrange shear interpolated at those distinct second coordinates sends the masks to the target first coordinates. Undo the linear shear. The inverse reverses the steps, and the chart has determinant one over \(\mathbb Q(\sqrt3,\sqrt{11})\). It decodes the named class, not a consistent global labelled seven-vertex lift; the odd holonomy proves the latter impossible.

## Replay

check_spindle_charts.py is standalone Python using fractions. It checks the actual strong labelled unit cube, both rational unit frameworks, the reciprocal shears, lost unmarked coordinates, and the complete strong quotient spindle in its exact algebraic metric. The [maximum-unit spindle census](check_spindle_max.py) independently enumerates the six-vertex quotient and all induced role assignments.
