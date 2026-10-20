# Cokernel form of the odd-primary Walsh correspondence

**Status: PROVED.** This strengthens the presentation of the theorem in `PRIME_VALUATION_THEOREM.md`; it is equivalent to the complete family of prime-valuation multiset identities there, and stronger than equality of the products of automorphism orders.

Let R=Z[1/2], and suppose

    M D_b M^t = c D_a,

where M belongs to GL_N(R), c is a unit of R, and D_a,D_b are diagonal matrices of positive integer weights. Then

    M(D_b R^N) = c D_a M^(-t) R^N = D_a R^N.

Consequently the explicit R-linear map

    [x] -> [Mx]

is an isomorphism

    R^N / D_b R^N  ->  R^N / D_a R^N,

with inverse [y] -> [M^-1 y]. This uses the equality of the two image submodules, not merely an inclusion.

Write odd(d)=d/2^v_2(d). Since powers of 2 are invertible in R, the additive group R/dR is naturally isomorphic to Z/odd(d)Z. Thus for the tournament Walsh matrix:

    direct sum over even graph types X of Z/odd(|Aut X|)Z
      is isomorphic to
    direct sum over tournament types T of Z/|Aut T|Z.

Tournament automorphism groups have odd order, so the odd-part notation is unnecessary on that side. The isomorphism sends the graph-coordinate vector to the tournament-coordinate vector by the signed Walsh matrix. Its inverse contains dyadic fractions, interpreted modulo odd moduli by inverting 2. In particular it does not pair the summands one-by-one when their joint odd parts differ.

Equivalently, the diagonal matrices D_a and D_b have identical Smith invariant factors over R. Decomposing the cokernels into their p-primary cyclic factors recovers, for every odd prime p, the equality of the multisets of v_p(a_i) and v_p(b_j), including zero valuations after padding to the common dimension N.

The five-vertex contrast is especially transparent. The tournament group is

    (Z/3Z)^4 direct-sum Z/5Z.

The graph group is

    (Z/3Z)^3 direct-sum Z/15Z.

These groups are isomorphic, since Z/15Z is Z/3Z direct-sum Z/5Z. Their displayed cyclic summands cannot be paired by order. This is exactly the distinction between a canonical linear transport of odd-primary torsion and a set bijection simultaneously preserving all prime valuations.

The two-by-two matrix in the main note realizes the same mechanism: Z/15Z is transported by [[3,1],[-5,1]] to Z/3Z direct-sum Z/5Z. The intermediate direct summand Z/1Z is trivial.
