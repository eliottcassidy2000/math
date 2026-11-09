# Independent audit of rich unsaturated modules

Verdict: PASS for the explicit modules, their full path-cover coordinates, ordinal composition, star closure, rigidity, and all-length consequences stated below. This is not an exhaustive classification or an all-target coverage theorem.

## Independent computation

`check.py` uses vertex-subset dynamic programming recording BOTH endpoints of each directed spanning path of every induced subset. Ordered two-path covers are then counted by products over complementary subsets. The producer instead enumerates full permutations and cuts. No producer code is imported. Explicit exceptions rather than assertions make the checks effective under Python optimization.

Both `python check.py` and `python -O check.py` pass with byte-identical receipts. All three bases are checked to be tournaments. Literal permutation actions independently show their automorphism groups have order one. For each gap-two pair, all 49 contexts of up to three singleton vertices, over every binary cut word and every prefix/suffix split, are checked by fresh literal graph DP against composition. Finite context tests are controls; the all-length claim follows from the formula proof, not extrapolation.

The rich six-vertex base [26,44,57,16,34,9] has H=29,Q=320. Its entire A=81 fiber is masks 9,10,17,18, with G=129,117,93,81 respectively. The full coordinate lists and independently counted strong-flank Hamilton path totals 673,691,719,729 are in CHECK.json.

The seven-vertex base [64,121,35,37,13,81,28] with masks 12,20 gives (H,Q,A,G,S,E)=(95,1276,312,270,24,18), T=38,36, and star counts 681,683. The eight-vertex base [20,85,72,195,204,159,33,71] with masks 47,181 gives (389,5600,3369,3456,250,271), T=1858,1860, and star counts 5269,5267.

## Historical convention comparison

Compared against git HEAD:05-knowledge/results/session216_overnight_constructions_20261107/three-core/PROOF.md and the sibling numeration/PROOF.md. Those proofs use (h,q,A,G,T,S0,E); this packet uses (H,Q,A,G,S,E,T). These are a coordinate permutation, not different statistics. In particular Q includes both possible empty-path cases, A allows the first path empty, G allows the second path empty, and T requires both paths nonempty. S and E count full Hamilton paths beginning and ending in the selected set.

The old closure is exactly L→v→x, with x→selected vertices and all unselected vertices→x. Its count is 2A+H−T, exactly the closure tested here. The old formula N=(qL−G)P−h f(P)−S(P mod 2)+AL+h(2R+1−f(R))−T−E first_b(R) is compatible: a change solely in T changes N by its negative, for every prefix/suffix word. No claim should conflate the strong-flank construction with this star closure.

## Proof review

For nonempty cut-labeled K1⇒K2, every cover path is its K1 portion followed by its K2 portion. Choosing these portions establishes Q'=Q1Q2 and H'=H1H2. Distinguishing whether the second path has a K1 portion gives A'=A1Q2+H1A2. The corresponding first-path endpoint split gives G'=G2Q1+H2G1. Full-path endpoints give S'=S1H2 and E'=H1E2. The four cases of whether the first path uses K2 and the second path uses K1 give T'=A1G2+T1H2+H1T2+E1S2. All four classes are disjoint and exhaustive, including empty portions.

Thus replacing a block while preserving its first six coordinates changes only T. In any fixed left/right ordinal context its T difference is multiplied by H(left)H(right). Singleton contexts have H=1, so preserve an exact gap two at arbitrary lengths. Both H95 and H389 are odd and coprime to their Q; adding t singletons changes Q to 2^t Q and leaves H fixed, so coprimality persists at every length.

The seven-vertex core is strongly connected (this follows directly by reachability from each vertex), and both cuts are nonempty. Its star closure is therefore strong. A singleton suffix of arbitrary length, with any selected subset, preserves this strongness: the selected core is reachable from x, the strong core reaches the suffix, and every body vertex reaches v→x. Consequently there is a strong gap-two pair arising from unsaturated cores at every order n≥9. Unsaturation is established for the deletion bodies, not necessarily for the star-closure outputs. Their precise counts vary with the context. Prefixes also preserve strongness when their first vertex is selected. An arbitrary prefix whose first vertex is unselected is NOT automatically strong. The eight-vertex example similarly starts at n=10. These are adjacent attainable-count pairs at every order, not coverage of every adjacent odd pair.

For the four-choice module, the length-k ordinal word has H=29^k, Q=320^k, A=81 sum_{i=0}^{k-1}29^i320^{k-1-i}, independent of its choices. Its G equals sum_{i=0}^{k-1}29^{k-1-i}320^i G_i. Equality of two such values modulo 320 forces equality of the first digits, since 29 is invertible modulo 320 and the four G_i have distinct residues. Subtract those terms, divide by 320, and induct. Therefore all 4^k words have distinct G for EVERY k, not merely the tested small lengths. Coprimality also persists. This is an invariant-separation result, not automatically a scalar-count separation result for every closure.

Finally Aut(B) acts freely on Hamilton paths: fixing an ordered vertex sequence fixes every vertex. It likewise acts freely on ordered spanning two-path covers, even when one path is empty. Orbit-stabilizer implies |Aut(B)| divides both H and Q, so gcd(H,Q)=1 implies rigidity. This argument does not depend on strong connectivity. It explains why a genuinely unsaturated base cannot hide a nontrivial automorphism factor.

## Limits

No search logs are used to prove absence outside a precisely specified finite universe. No minimal-order claim, exhaustive rich-module classification, count-coverage theorem, or new universal transfer theorem follows from these witnesses. The local gap-two repair solves a concrete interface problem; the all-target inverse/numeration problem remains separate.

## Final producer proof and central-defect addendum

Reviewed the finalized producer PROOF.md, including its central-defect and bilateral-rigidity sections. PASS. With categories as in the inherited six-state convention, the eight weights are S,H−S,H−E,E,A−S−T,Q−2H−A+S−G+E+T,T,G−E−T in order (01,02,10,20,11,12,21,22). Thus a T-only difference δ changes (11,12,21,22) by (−δ,δ,δ,−δ). Direct substitution in the six-state update gives exactly F=−δ(e3−e5)(e0*+e2*+e4*), including the sign claimed by the producer. The independent checker now verifies this matrix identity and both eigenvector identities on every witness and both singleton cuts.

The universal identity Mw=Hw holds because nonempty first-path pieces erase the previous endpoint, while empty first-path pieces retain it and count H. The identity rho M=H rho counts precisely covers with second path empty. Hence the central defect commutes up to multiplication by H and squares to zero. With k identical switches, fixed first-six-coordinate interface and equal local δ, the total T shift is δ times the number of switched positions times H^(k−1). The scalar gap magnitude is therefore 2mH^(k−1), not 2^k distinct values. Strongness of both switch cores is also checked explicitly by directed reachability.

Bilateral rigidity follows because both A and G local coordinates lie in [0,Q−H], a subset of [0,Q). The modular argument fixes A digits from one end and G digits from the other when gcd(H,Q)=1. Thus cancellation of different local G words cannot disguise their full interface. No issue found with these statements.

The late exploratory `extend_triples.py` output is outside this independent audit. No negative assertion from that additional extension search is certified here.


## Release correction: unsaturation does not survive arbitrary star closure

The initial version of this review incorrectly called the output pairs themselves unsaturated. This was an overstatement: the proof preserves gcd(H,Q)=1 under singleton ordinal padding of the body, not under the final star closure. The wording above is corrected. Independent endpoint-subset DP now verifies the closure triples (H,Q,gcd(H,Q)): B95 masks12/20 give (681,15600,3) and (683,15664,1); B389 masks47/181 give (5269,111760,11) and (5267,111608,1). Thus one output in each displayed pair is saturated. The strongness and exact gap-two conclusions are unaffected. These hostile controls are included in both normal and optimized receipts. No claim that every output is unsaturated should be inferred from this audit.
