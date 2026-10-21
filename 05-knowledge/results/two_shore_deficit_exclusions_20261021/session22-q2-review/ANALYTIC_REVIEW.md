# Independent two-shore insertion theorem review

Date: 2026-10-21 UTC. Verdict: PASS. This separate structural-family theorem does not edit or enlarge the frozen degree47 proof packet.

For strong self-converse T of order 2m, odd m>=5, with exact marked-cycle index q=2, the theorem Delta>=2m-2 is sound.

The full shore scores differ because their difference is 2k-m, nonzero for odd m. They sum to 2m-1. Thus automorphisms preserve shores and anti-automorphisms exchange them. After choosing the high-cross shore as A, k>m/2; marked-cycle cross arcs in both directions force k<=m-1, so the terminal v in B has both types of neighbors in A.

The punctured regular-shore lemma is correct: the regular connectivity averaging argument with one deleted vertex forces both nonempty surviving SCC blocks to have size at most one, incompatible with m>=5. Hence B minus v has a Hamilton cycle and a Hamilton word starting at any chosen outneighbor of the final A vertex.

On a directed Hamilton cycle of A, the binary classification of vertices according to their orientation against v has both values. A cyclic transition from a vertex beating v to one beaten by v supplies consecutive a,a' with a→v→a'. Every cut positioning a in position i, 1<=i<=m-2, yields an actual insertion word. The last A vertex has at least k-1>=2 outneighbors remaining in B minus v, permitting an individually chosen suffix. The separate P0 construction starts at v and then an outneighbor in A. No unsupported common-suffix assumption occurs.

The m-1 patterns are pairwise distinct. Their reverse-complement mates all begin with exactly m-1 copies of A, unlike the originals, which either start in B or have initial A-run at most m-2. Mate patterns are themselves distinct. As automorphisms preserve shores, these 2m-2 patterns certify different Hamilton-word orbits. None can be fixed under the converse/reversal involution, so all are black. The argument requires self-converseness only when turning reverse-complement patterns into actual mate orbits.

Both stated hostile boundaries are genuine: the omitted i=m-1 pattern is reverse-complement fixed, and a three-vertex regular shore loses strongness after one deletion. No claim of additional black orbits or m=3 exclusion is justified by this mechanism.

At n=10 or14, q<=8 and odd quotient n/q force q=2. Thus the theorem analytically excludes both orders for deficit six. This is a necessary family exclusion, not a classification of all remaining centers or metagraph adjacency.

Fresh copied structural producer runs under normal Python and -O each pass 21,504 always-active checks. Their nonvolatile receipts agree at source SHA256 e847e544503d7f16d41f519e174076c4bb2e1aee72c8c6d9d163f17628e55cda. The 128 m=5 and 2,048 m=7 frames are exactly the scoped marked construction universe, not a census of unmarked isomorphism classes. An independent literal consumer validates all twelve retained sample receipts, actual directed adjacency, exact vertex sets, intrinsic shore scores, original/mate patterns, and their disjointness. A negative-label corruption is independently rejected.

The order10 counting audit is recorded separately in ORDER10_REVIEW.md. Source pins and concise replay results accompany this review. Only session22-q2-review was written, with no Git or remote operations.
