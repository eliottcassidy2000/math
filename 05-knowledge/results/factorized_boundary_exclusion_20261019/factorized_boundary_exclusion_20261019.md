# A factorized-boundary exclusion beyond the independent-generator cone

Packet 144, 2026-10-19. Mathematical certificate reviewed and approved; frozen repository-ready packet.

## Exact scoped result

Let Q have predecessor masks (0,28,0,0,12,4,13,93,125), and let M=(2,3,0,4,5,6,1,7,8). No nonnegative factorized boundary model with U_empty=V_empty=1 can give probability strictly greater than 2/3 to every M-forward incomparable pair.

This excludes all factorized models of the specified core and order, not merely the particular seven-generator cone witness. It does not claim that every core/order remaining in the (3,3,1) catalogue is excluded. It imposes no exterior-degree bound, disjunctive common-vertex condition, or realizability condition on the nonempty increments. Thus it is distinct from the simpler actual-exterior split certificate in packet145.

## Boundary model and variables

For each core extension w, a generator G_(u,v) requires w to begin with u and end with v, with both roots 0,1 between these two boundary words. E_(u,v) is its extension count. For a forward pair a<b in M, write

    C_ab(U,V) = sum_(u,v) U_u V_v [3 N_(u,v)(a<b) - 2 E_(u,v)].

Every retained nonempty word gets an independent nonnegative variable; U_empty=V_empty=1. There are 16 possible prefixes and 16 possible suffixes (including empty), with 175 nonempty generator supports among their 256 combinations. Unsupported combinations contribute zero, not an additional constraint. All 175 supports are retained. The core degrees are (5,5,2,3,3,6,3,3,2), so there are no saturated-vertex omissions.

The 540 core extensions and all generators are independently reconstructed by replay.py. The total mass is positive because the empty-empty term contributes 540. Consequently all forward probabilities exceeding 2/3 is equivalent to all 16 polynomials C_ab being positive.

For actual exterior attachment, ambient counts have the usual factorization L(pre(w))R(suf(w)); writing each L and R as a sum of its nonnegative word-tree increments gives exactly the model above. The proof below only needs these nonnegative increments and their empty normalization. It does not assume the increments can vary independently in an actual poset.

## Explicit certificate

Set a=U_(2) and b=U_(2,3,4,5). The following polynomial is coefficientwise nonpositive in all 30 nonempty U,V variables:

    F = C_30
      + (58 + 24a + 25b) C_04
      + (62 + 16a + 25b) C_61
      + 2 C_17
      + (7 + 3a + 4b) C_23
      + (67 + 22a + 26b) C_45
      + (62 + 17a + 27b) C_56
      + (8 + 2a + 4b) C_78.

Here C_30 denotes the oriented event 3 before 0, and likewise for every subscript. There are eight pair multipliers, represented by 20 monomial terms in small_certificate.json. Every multiplier is strictly positive throughout the nonnegative orthant.

Expanding F yields 486 strictly negative coefficients and 9 zero coefficients among the monomials encountered; absent monomials also have zero coefficient. The largest coefficient is zero. Therefore F<=0 for every nonnegative boundary array. If every forward incomparable probability exceeded 2/3, all eight displayed C polynomials would be positive and F>0. This contradiction proves the scoped result.

## Reproduction and trust boundary

Run python replay.py from this result folder, or invoke that script by its full path. Verification remains active under python -O. The replay uses Python standard library integer arithmetic only: it enumerates all 9! permutations, checks the predecessor masks directly, reconstructs generator margins, forms the displayed polynomial, and checks every coefficient. It imports neither the scan engine nor NumPy/SciPy nor an optimizer. It prints a JSON receipt and writes expanded_polynomial.json. verification.json and verification_optimized.json record matching normal and optimized runs. The verifier explicitly pins the two named multiplier words and all oriented pair labels.

Numerical search was used only for discovery and is excluded from this packet. The exact replay is the authority and reads only small_certificate.json, with no fallback certificate.

The independent-generator cone witness remains valid: coefficients of separate generators can be chosen to make all margins positive. This result identifies a genuine limitation of that relaxation. Retaining the product identities U_u V_v eliminates the specific obstruction, even without invoking exterior-degree slack.
