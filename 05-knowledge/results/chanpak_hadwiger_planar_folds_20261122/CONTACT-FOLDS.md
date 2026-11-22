# Actual two-contact cuts: complete planar fold classification and a metric obstruction

**PROVED.** Start with a marked literal unit pair \(T_-,T_+\), with Hamilton-path counts differing by two. A common apex beats a contact set \(S\) in both hosts; write \(q_\pm(S)=(H(T_\pm+x_S)-H(T_\pm))/2\). These are the exact normalized quadratic cuts of the [previous packet](../chanpak_mirimanoff_unit_lifts_20261122/CUT-GAUGE.md). Pin every contact except two, denoted \(x,y\). In this face the cut pair is
\[
             F(x,y)=P+A x+B y+Cxy,\qquad A,B,C,P\in\mathbb Z^2.
\]
On Boolean inputs the values are literal cut counts. At fractional inputs they are **expected** cuts with independent random choices for the two contacts, not graph witnesses. The vectors \(C_-,C_+\) are both strictly negative: each is minus the positive normalized sum of the two rooted endpoint weights for that contact pair. Thus all actual two-free-contact faces have \(C\ne0\). The results below also describe arbitrary real constants in this form.

## Bounded inverse, fold and collapsed line

In output coordinates, take one linear combination eliminating \(C\); after subtracting a constant this is \(\alpha x+\beta y\). Take a second independent combination to retain a nonzero \(xy\) coefficient.

* If \(\alpha\beta\ne0\), prescribe the first coordinate and use it to eliminate \(y\). The second read is now quadratic in \(x\), with nonzero leading coefficient independent of the target. Therefore the complex fiber has at most two points, and real solutions in the square are recovered by checking both real roots and contact bounds. The target change removing a polynomial in the first coordinate, and affine completion of the square in the remaining input variable, reduce the map to \((u,v^2)\) up to nonzero rescalings and these coordinate changes.
* If exactly one of \(\alpha,\beta\) vanishes, the first coordinate fixes a single input coordinate. The other read is linear in the leftover coordinate with coefficient a nonconstant affine function of the fixed coordinate. This is polynomially equivalent to \((u,uv)\) by translations, scales, and subtraction of an affine function of \(u\); off the exceptional line it has one rational inverse, on the line it has either no solution or an entire input line.
* If the first coordinate is constant, the output lies in a line. This case cannot be injective on an open square.

The determinant is the affine function
\[
   J(x,y)=\det(A,B)+x\det(A,C)+y\det(C,B).
\]
If it is a nonzero constant then \(C=0\): otherwise \(A,B\) would both be proportional to \(C\) and \(\det(A,B)=0\). Thus **no directly interpolated two-contact tournament cut face is a Keller map of the plane**. This is a restricted statement about bilinear cut maps, not a conclusion about arbitrary polynomial Keller maps. The fold and collapsed-line alternatives retain the fiber information that the [earlier higher-dimensional fold revision](../chanpak_two_tape_banks_jacobian_anatomy_20261115/JACOBIAN.md) emphasized under a different construction.

## A sharp square-domain injection criterion

For every real quadratic map, hence for this one,
\[
          F(z)-F(w)=DF((z+w)/2)(z-w).
\]
It is the exact symmetric difference identity: the even Taylor terms cancel. Consequently an interior zero of \(J\) gives a genuine collision of two distinct points in the square: take a short nonzero segment in the kernel centered at that zero. Conversely, if \(J\) is nonzero in the open square, any collision must have its midpoint on the boundary. Two points with boundary midpoint must lie on the same boundary edge, and an edge restriction is affine. Therefore **the square is injective exactly when there is no interior zero of \(J\) and no collapsed edge**. Both tests are finite: \(J\) is affine so its corner readings specify its sign and zero set; compare the adjacent endpoint values to detect edge collapse. Zero determinant at a boundary corner need not itself break injection.

This gives geometric certificates for interpreting relaxations in the [shared-contact selector](../chanpak_middle_out_square_sum_20261122/README.md). For absolute placement, after elimination every surviving \(x,y\) still has to be Boolean. An injective fractional face does not assert a target lies in its four-point physical image.

## Why a paired cut face cannot itself be a unit-distance square

Suppose the four cut vectors are distinct and, in these output coordinates, all four one-contact changes have the same positive Euclidean length. The two opposite intermediate corners are then the distinct intersections of two equal-radius circles centered at the other two corners. The perpendicular-bisector equation gives
\[
         F(1,0)+F(0,1)=F(0,0)+F(1,1).
\]
That forces \(C=0\), contrary to the two positive contact weights. The same argument applies after **any invertible affine change of the two read coordinates**, because vanishing of the mixed difference is invariant. Thus a literal unit-distance depiction of the combinatorial contact cube is not the same as using its two Hamilton-count readings as planar geometric coordinates. Both metrics cannot be imposed even on one two-dimensional face. This is the smallest obstruction, unlike a global chromatic argument; it coexists with the actual [Hadwiger–Nelson port amplifiers](../session254_four_frontiers_20261111/hadwiger_nelson/PROOF.md), which keep boundary states rather than imposing those two coordinates on their vertices.

## Physical fixtures

Below the vertex labels are zero-based. The outside contact mask is the set of fixed wins, with both free bits initially off. The cuts are in order \(00,10,01,11\). In every row the free pair is the **same pair used for the old single-arc reversal**.

| Lower outgoing rows | Flipped pair | Old counts | Outside mask | Cut vectors | Corner Jacobians |
|---|---|---|---:|---|---|
| 16,1,3,7,14 | 0,3 | 9,11 | 0 | (0,0),(4,6),(12,9),(14,13) | -36,-32,-30,-26 |
| 8,17,3,6,13 | 0,1 | 13,15 | 16 | (12,9),(12,14),(15,9),(12,11) | -15,0,-6,9 |
| 2,24,3,5,13 | 0,3 | 11,13 | 16 | (13,12),(13,16),(17,12),(13,12) | -16,0,0,16 |
| 2,36,17,7,11,29 | 4,5 | 29,31 | 14 | (43,43),(44,52),(50,42),(43,43) | -64,0,0,64 |

The first face is injective throughout the square. The second has four different physical labels but folds the relaxation: \((16/25,2/5)\) and \((19/25,3/5)\) both map to \((1554/125,1429/125)\). The third has a **literal Boolean collision** between two opposite masks whose Jacobians are both nonzero. Their common cut pair is \((13,12)\). The fourth face goes further: its two opposite masks 14 and 62 both preserve the literal +2 gap and produce the same ordered pair **115,117**, despite determinants −64 and 64. The two pairs are not simultaneously isomorphic (the replay tests all seven-vertex relabellings), even though the lower hosts are isomorphic separately and so are the upper hosts. Retaining both unmarked classes still loses how the common flip is glued. [Replay the four rows](check_contact_folds.py) by direct Hamilton-path DP.

## A finite metric rigidity test for polynomial charts

A polynomial \(G:\mathbb R^2\to\mathbb R^2\) of degree at most \(d\ge1\) is forced to be affine if for both standard unit vectors \(e_j\) the identity \(\|G(z+e_j)-G(z)\|=1\) is verified on \(z\in\{0,\ldots,2d-2\}^2\). The squared norm minus one has degree at most \(2d-2\) and vanishes on the product grid; univariate interpolation in each coordinate makes it identically zero. Each coordinate of each difference is thus a bounded polynomial, hence constant. Subtract the linear map with those two constant columns. The remainder is polynomial and periodic in both standard directions, so constant. If \(G\) also has \(|\det DG|=1\), the two unit columns must be orthogonal by equality in the area bound, hence \(G\) is an isometry. Alternatively a single extra \(3/5,4/5\) unit direction tests orthogonality once affine rigidity is known. This is a finite **metric** test complementary to a determinant test; a shear shows that the determinant alone does not certify distances.
