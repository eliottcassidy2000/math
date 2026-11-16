# Sharp one-vertex growth for all span-3 ribbons

**PROVED.** If T is any nonempty span-3 ribbon and x is the next vertex of a span-3 ribbon, then
`H(T) ≤ H(T+x) ≤ 39 H(T)`.
The factor 39 is attained. No transitivity assumption is needed. This closes the bounded-growth part of [HYP-9802(ii)](https://github.com/eliottcassidy2000/math/blob/a2ec2ffce0d1a7ab3f57b976f1368570c35b7e0d/05-knowledge/hypotheses/HYP-9802-literal-span-3-certificates-need-a-bounded-growth-lemma-and-rank-adaptive-memory-short-cheap-cycles-are-cut-within-nine-letters.md), not its proposed lower-growth rank certificate.

## Delete the obstructing tail, not just the new vertex

Write x=n and use old labels 0,…,n−1. Let A=N⁺(x), contained in the last three old sites. By the odd-cycle packing identity of [Grinberg–Stanley](https://arxiv.org/abs/2307.05569), H is monotone under extension, and for every old subset S,
`H(T−S) H(T[S]) ≤ H(T)`
(with H(empty)=1). Indeed disjoint packings on the two separate supports give distinct packings in T with the same weight.

Paths ending at x contribute at most H(T), by deleting x. Otherwise write a Hamiltonian path as P₁,x,P₂, and let b∈A be the first vertex of P₂. Remove from P₁ its maximal terminal block of vertices all beaten by b; keep this block in its original order, and denote its vertex set by S. If a vertex of P₁ remains, its new last vertex beats b. Thus the remainder followed by P₂ is a Hamiltonian path in T−S, also when no vertex of P₁ remains. For fixed b and S, the Hamiltonian path on T−S and the ordered block recover the original by locating b and inserting the block and x immediately before it. The number for these fixed data is at most H(T−S)H(T[S])≤H(T). If S is nonempty its last vertex beats x, so S contains a vertex outside A.

## The local count and its sole overcount

For b=n−j in A, let d be its outdegree inside the last three old sites and c its outdegree inside A. There are at most j additional out-neighbors of b to the left of the last three sites: its backward arcs can span at most three. Consequently there are at most `1+2^(j+d)−2^c` possible sets S, including the empty set. If b has fewer than j such earlier out-neighbors, this bound is strictly larger than the true subset bound.

Add these quantities over b and add 1 for the paths ending at x. For each A mask, the maximum over the eight orientations of the last three sites is:

| A, mask on offsets 1,2,3 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Counting bound | 1 | 9 | 17 | 20 | 33 | 36 | 40 | 39 |

The sole 40 case has A={n−2,n−3} with the final old triple pointing forward. If any earlier out-neighbor counted at the maximum is missing the integer bound is already ≤39. Otherwise for b=n−3 the subset S={n−2,n−6} has been counted but is impossible. Its only directed order is n−6→n−2 (span four forces that orientation), which would put n−2∈A immediately before x. This contradicts the required arc into x. Remove this term, leaving 39. With fewer than three old vertices the whole new graph has at most 3! possible paths anyway.

The earlier sharpness fixture uses old transitive order (3,0,4,1,2,5) and new vertex 6 beating {3,4}. Its outgoing rows are (118,100,96,55,38,64,24), with H(T)=1 and H(T+x)=39 by exact subset DP. The companion checker replays the local table, its unique maximum, the excluded pair and this fixture.

## Consequences for the potential

Every finite suffix readout of length ℓ is now bounded by 39^ℓ H on physical prefix states. Thus a finite linear combination of continuation readouts is capped in absolute value by the corresponding sum of absolute coefficients weighted by 39^ℓ. The three-step positive sum in [the earlier envelope](POSITIVE_ENVELOPE.md) is capped by (1+312+312²+312³)H=30,468,985H, uniformly on the actual orbit. In particular its lack of a physical cap is no longer an open problem. The independently certified positive frontier row r there is pointwise ≤(11/208) times this sum row, so it too has a uniform physical cap.

The proof works more generally at fixed span w: the same deletion allows b in the last w places and at most w+j−1 out-neighbors (w−1 in the window and at most j earlier ones). Thus `H(T+x)≤(1+w·2^(2w−1))H(T)` for any fixed w≥1; if fewer than w old sites exist, use the actual window and its old out-neighbor counts, which give the same bound. The span-3 enumeration is what improves this crude bound to the sharp 39.

None of these upper comparisons proves a lower growth inequality on rank-preserving literal stretches, supplies the missing image states, or proves a power saving for missed counts. It isolates those as separate proof tasks.
