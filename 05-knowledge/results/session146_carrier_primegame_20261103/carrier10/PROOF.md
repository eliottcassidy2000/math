# A certified order-ten lower bound

This result builds on the previous pooled-parity carrier theorem and the independently certified Session 146 exact results through order nine. We consider rational linear mean responses to all fixed source-oblivious attachments, not an asserted S/B-only query bank. No claim is made about source-dependent attachments, conditioning, or shared-source products.

## Theorem

At ten named old vertices, the response dimension d(10) satisfies

    345106 ≤ d(10) ≤ 347626.

For normalized source laws the corresponding bounds on varying coordinates are 345105 and 347625. The upper bound is the pooled-parity carrier count. The lower bound combines the exact order-nine lifting with independently verified endpoint-set block ranks. The exact order-ten minimum and all-order sharpness remain open.

## Proof

All 299701 signatures having at least one isolated old vertex are individually observable by lifting the exact order-nine response span. Among full-support signatures, group by the named degree-one set E, with |E|=2r. As in the order-nine proof, an exterior context touching every vertex of E is zero on a full-support forest with endpoint set E' unless E⊆E'. The resulting row/column matrix is triangular by endpoint-set inclusion. Different named endpoint sets of equal size give diagonal blocks equivalent up to relabelling and nonzero sign gauges.

For E={0,...,2r−1}, the independent verifier reconstructs the signed full-support representatives by an exact-cover census, checks endpoint saturation and evaluates the coefficient law. The number of columns, the number of certified independent rows over each of GF(101) and GF(103), and, where square, the determinants are:

    r   columns   independent rows   determinant mod 101 / mod 103
    1       1            1                    4 / 4
    2       6            6                   88 / 10
    3      60           60                   12 / 47
    4     735          679                 rectangular
    5     945          945                   21 / 18

The rectangular r=4 bank proves only a lower bound of 679 on that block's rational row rank. In particular its shortfall from 735 is not a proof of a relation or a dimension upper bound. Deletion and duplication controls reduce each nontrivial certified row rank by one in both prime fields. All modular denominators are powers of two, so the ranks imply rational independence.

Across all endpoint sets these independent full-support directions contribute at least

    C(10,2)·1 + C(10,4)·6 + C(10,6)·60
      + C(10,8)·679 + C(10,10)·945 = 45405.

Projecting onto the full-support quotient proves their independence modulo the 299701 individually observable smaller-support coordinates. Hence d(10)≥299701+45405=345106. Boolean multilinear interpolation puts each exterior-coefficient row in the rational span of actual fixed attachments. □

## Lifting the partial result to all orders

Let

    L9(n) = 1 + C(n,2) + 3 C(n,3) + 9 C(n,4)
          + 40 C(n,5) + 120 C(n,6) + 546 C(n,7)
          + 2233 C(n,8) + 9612 C(n,9).

For every n ≥ 10,

    L9(n) + 45405 C(n,10) ≤ d(n) ≤ c(n),

where c(n) is the pooled-parity carrier count. To prove the additional term, choose each ten-element old subset A and lift the 45405 independent full-support quotient rows. Give every old vertex outside A two exterior leaves, which kills any forest using that outside vertex; adding isolated new vertices and rational interpolation in the component count restores each order-ten row on forests supported in A. The projections of these lifted rows onto full-support signatures on A are independent, and the projections for distinct A have disjoint supports. All smaller-support coordinates are individually observable from exactness at order nine. This proves the displayed sum without asserting that the partially observed order-ten coordinates are individually observable.

Examples are 1654181 ≤ d(11) ≤ 1902231 and 6796879 ≤ d(12) ≤ 10777900. These bounds include the constant coordinate; subtract one for normalized laws.
