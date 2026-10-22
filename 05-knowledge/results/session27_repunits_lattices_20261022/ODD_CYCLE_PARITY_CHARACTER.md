# A Boolean reconstruction coefficient of odd-cycle parity

Status: PROVED, using the classical odd-cycle collection formula. This is a second proof of the fixed-endpoint response contrast, not a separate claim of priority for finite-difference reconstruction.

For a tournament T, let b(T) be the parity of the number of all directed odd simple cycles of length at least three, with cyclic rotations identified. The odd-cycle collection formula gives

    H(T) = 1 + 2 a1(T) + 4 a2(T) + ...,
    b(T) = ((H(T)-1)/2) mod 2.

Thus b is exactly the second least significant bit of H, not the constant parity bit. Here a1 counts ALL odd cycles, not only triangles.

Fix distinct labelled vertices u,v in a k-vertex core, k>=2. Adjoin k-1 ordered transitive base vertices. Arrange every vertex exceptu along the alternating path

    R=(v,b1,w1,b2,...,w_(k-2),b_(k-1)),

where the w's are the remaining core vertices. Fix b_(k-1)->u. Independently vary the orientations of the 2k-3 cross edges underlying R. Keep every other attachment arc arbitrary and fixed.

**Theorem.** The XOR of b(T) over these 2^(2k-3) extensions equals the indicator of the original core arc u->v.

**Direct proof.** Work over F2 and interchange the sum over edge assignments with the sum over cyclically ordered odd cycles. If a candidate cycle omits one varied edge of R, its contribution occurs in equal pairs under reversal of that unused edge and cancels. A cycle containing every underlying R edge must traverse R as a single block. Closing R directly from its last vertex back to v gives a cycle of length2k-2, which is even and is not counted. The only remaining vertex is u; using it closes a cycle of odd length2k-1. The fixed edge b_(k-1)->u excludes the reversed traversal, so exactly one assignment contributes precisely when u->v is present. This proves the XOR identity. ∎

For k=2 there is one varied edge, and the test reduces to whether a three-vertex extension is cyclic. The construction retains vertex labels, the endpoint direction, and every varied edge's address. It makes no assertion that a bare odd-cycle count classifies tournaments without contexts.

The query family is a high-order Boolean finite difference. Its retained coefficient records an actual arc, not a fitted numerical correlation. Over all unordered vertex pairs it yields a complete labelled tournament invariant from a single bit-valued observable under finitely many contexts. Relabelling and converse must act consistently on the contexts to obtain the corresponding isomorphism/converse quotient.

**Why triangles alone are insufficient.** At order five, the two explicit same-score examples in SCORE_ROOT_COLLISION.md both have four cyclic triples but have respectively two and one directed five-cycles. Their H values are13 and11, with different second bits. Replacing all odd cycles by just triples would erase exactly the signal used here.

**Source.** The OCF identity is pinned in repository THM-002-ocf.md and follows from Grinberg–Stanley, *The Redei–Berge symmetric function of a directed graph*, with the mod-four refinement in Theorem7.1 of [arXiv:2307.05569](https://arxiv.org/abs/2307.05569). The response packet provides independent full-path permutation and subset-DP tests; this proof supplies a different direct cancellation mechanism.
