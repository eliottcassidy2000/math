# Exact order-nine attachment dimension

This result builds on Session 145 `POOLED_PARITY_CARRIER.md` and the separately certified order-seven and order-eight results, including their lifting theorem. We consider rational linear mean responses to all fixed source-oblivious attachments, not merely an S/B-only port language. No statement is made about source-dependent attachments, conditioning, or shared-source products.

## Theorem

At nine named old vertices, the span of all fixed attachment mean responses has dimension exactly **65908**. For normalized laws, the constant coordinate is known, leaving exactly **65907** varying coordinates.

The pooled-parity carrier gives the upper bound 65908. For the lower bound, all 56296 coordinates with at least one isolated old vertex are individually observable by the saturated-vertex lifting of the exact order-eight certificate. The quotient on the remaining 9612 full-support signatures is shown full rank below. Exterior coefficient rows lie in the rational span of actual fixed attachment responses by exact Boolean multilinear interpolation.

## Four independent diagonal blocks

For any old forest with no isolated vertex, its degree-one set E has size 2r, for 1 ≤ r ≤ 4. The nine−2r remaining vertices are degree two. A fixed exterior context that touches every vertex of a set E gives a zero coefficient on any full-support forest whose degree-one set E' fails to contain E: a vertex of E that is not in E' already has old degree two and acquires degree at least three.

Arrange the full-support signature columns and these context rows by decreasing |E|, grouping equal E. A row for E is zero outside column blocks E' containing E, so the resulting square matrix is block triangular. At equal size, E ⊆ E' forces E=E'. The diagonal block depends, up to row/column permutations and nonzero sign gauges, only on |E|: relabel old vertices and carry the exterior context along. Sign gauges arise from the antisymmetry of centered edge variables and the canonical orientation of each old path. Thus it suffices to certify one diagonal block for E={0,...,2r−1} for each r.

The corresponding numbers of signatures, and determinant checks on the supplied square context matrices, are:

    r   dimension   determinant mod 101   determinant mod 103
    1       1               41                    43
    2       6               66                    74
    3      60               90                    64
    4     420                6                    65

The independent verifier reconstructs the representatives using a simple-path exact-cover enumeration, checks that each context touches every vertex of E, and recomputes both modular ranks and determinants. Each determinant is nonzero. Coefficient denominators are powers of two, so this proves rational invertibility. Deleting a row or duplicating a row lowers rank by one; for r=1 deletion gives the zero-rank empty matrix.

There are respectively C(9,2), C(9,4), C(9,6), C(9,8) endpoint sets. The full-support block dimensions sum to

    C(9,2)·1 + C(9,4)·6 + C(9,6)·60 + C(9,8)·420 = 9612.

Combining with the 56296 already observable coordinates proves the claimed dimension 65908. The quotient projection is only a device for establishing independence; it is not itself asserted to be a single attachment.

## Consequence at every higher order

By the same saturation and component-count interpolation lifting used for order seven, exactness at order nine makes every signature supported on at most nine old vertices individually observable at order n ≥ 9. Thus

    d(n) ≥ 1 + C(n,2) + 3 C(n,3) + 9 C(n,4)
        + 40 C(n,5) + 120 C(n,6) + 546 C(n,7)
        + 2233 C(n,8) + 9612 C(n,9),

while the pooled-parity count c(n) remains the upper bound. In particular 299701 ≤ d(10) ≤ 347626 and 1154726 ≤ d(11) ≤ 1902231. These counts include the constant; subtract one for normalized laws. Sharpness at order ten and beyond is not claimed here.
