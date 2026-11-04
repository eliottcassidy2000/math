# Independent review: incoming deletion census and named-response exclusions

## Verdict and boundary

PASS for the precise finite conclusions in the producer's frozen PROOF.md (SHA-256 f40ce81caf1b8f8bbf39ac374c95ba3f7e7991168a8ea32062f2d3445785a2b8).

All three reported order-nine equal unweighted converse-merged deletion-deck pairs are excluded as full named two-port response collisions under every old-name bijection and each of the four independent converse choices. This is twelve orientation comparisons, not an enumeration or classification of all order-nine response collisions. The response keeps the port edge x→y fixed. No implication extracting the weighted deleted-card observer W from the full named response F is assumed or proved.

The incoming full census and Gram ranks have been reproduced from the exact original C source. This is not a new independent implementation of the full census or rank calculation. The response polynomials and all cards of the six reported examples do receive separate independent checks below.

## 1. Independent literal full polynomial

The new literal_polynomial.cpp does not use either of the producer's subset Hamilton-path recurrences. For each of the twelve physical graphs it enumerates all 9! permutations of the old names and all possible positions of the two ports.

A flattened old permutation can have at most two backward consecutive edges: each port can interrupt at most one such adjacency. If both ports occupy one gap, their only permitted order is x,y. If they occupy distinct gaps, the two intervening port-edge factors are

    a_predecessor (1-a_successor) c_predecessor (1-c_successor),

with nonexistent endpoint factors omitted. A common gap instead contributes a_predecessor (1-c_successor). These terms are expanded directly, with their signs, into physical named coefficients. Every augmented Hamilton path has a unique old permutation and placement, and every accepted placement gives exactly that path. Thus there is no factorial normalization or quotient-name aggregation.

At most two a variables and two c variables occur. The full basis has 46 masks in each family and 2,116 coefficients per graph, including all zeros. The independent literal expansion agrees with every coefficient and every basis evaluation in the producer bank:

- 12 graphs and 25,392 coefficients in total;
- 7,776 degree-three coefficients;
- 15,552 degree-four coefficients.

literal_bank.json retains these complete independent source/converse polynomials, evaluations, and named old-path start and end counts. This engine is different from predecessor.cpp, which was copied from a prior audit and must not be counted as an additional new independent algorithm here.

## 2. Independent deck identities and physical labels

canonical.cpp exhausts all permutations using row-major adjacency coding, rather than the producer's column-major coding or the incoming color-refinement canonicalizer. Minimizing both a code and its bit complement gives the converse-merged representative. The multiplicity of the minimum code equals the automorphism order; equality of the ordinary and converse minima detects self-converseness.

All twelve whole graphs and all 108 named deleted cards are independently canonicalized. The audit checks a bijection between the two different encodings' observed card-class identifiers, preserving each deleted physical name and its separately enumerated Hamilton endpoint weights. It verifies:

- the six original representative codes decode to exactly the retained physical rows;
- each paired orientation is the exact edge converse on those same names;
- each reported pair consists of distinct whole isomorphism-and-converse classes;
- the unweighted merged decks agree, with multiplicities;
- the primitive weighted decks differ and are unchanged by converse;
- all reported Hamilton counts, automorphism orders, and self-converse flags agree.

row_major_cards.json and exclusions.json retain the full records and the cross-encoding class correspondence.

## 3. All old-name bijections, all converse orientations

For a fixed port orientation, an old-name bijection merely permutes the nine coefficients [c_v]F. Therefore unequal multisets of these coefficients rule out equality of the full polynomials under every such bijection. This logical implication requires no enumeration of the 9! name bijections.

The audit uses multisets with multiplicities, not sets, and checks the Cartesian product of both converse choices on both sides for each of the three pairs. All twelve comparisons fail the necessary linear-multiset equality. exclusions.json retains every comparison and both complete multisets. This suffices even where the scalar Hamilton count is the same. Conversely, equality of these multisets would not be treated as a collision.

The displayed canonical H=645 example originally omitted the value 1070 from its second list. The author corrected that prose-only omission before this final audit; the retained code and all coefficient data were already correct.

## 4. Source and rank claims

The exact incoming commit is e58a09c18bc7536e2388c85df239afd9d10a9c1d. Targeted GitHub retrieval independently matched all three retained Git blobs:

- C source: 98953760287aa68650fb91728adea9686ea03e6d;
- original output: 550a9cf686bdf6f130cc5db71fd2bb8e74a285f0;
- HYP-9455: bd2b6ba5fa5fc446414c7726ab158451a3d1ee2b.

The untouched C source is compiled and rerun through order nine. Every non-timing output line matches. The exact source export regenerates the retained representatives and graph input byte-for-byte. The reported finite counts are order-eight 6,880 ordinary / 3,528 merged classes and order-nine 191,536 ordinary / 97,144 merged classes, with exactly three two-node unweighted-deck collisions at order nine.

The weighted observer is W(T)=sum_v (S_v+E_v)e_[T-v], and its total is 2H(T). Converse exchanges S and E and leaves a converse-merged card class fixed, so W is converse invariant. This differs from the named fixed-port response. Neither distinct weighted rows nor their rank supplies an unproved W-from-F extraction.

The stated full modular Gram rank does imply full rational column rank: the integer Gram determinant is nonzero modulo the prime, hence nonzero over the rationals. The audit independently checks that the modulus 2^31-1 is prime. No reverse inference from deficient modular Gram rank is made. The finite row uniqueness and Gram computation themselves remain source reproductions.

The later reversal-cube preregistration is not evidence for a proved separator. Its reserved theorem numbers and proposed outcomes are not dependencies of any accepted result here. The retained source metadata documents the producer's latest-source check; this review's independently retrieved source assertions are pinned to the exact commit above, not a claim that later source states were exhausted.

## 5. Reproducibility and release

The independent audit replays the producer verifier in normal and optimized Python modes in a clean temporary directory. Both outputs match its frozen receipt. The independent audit itself is run in both modes; its stdout and receipts agree byte-for-byte. All input producer files are pinned before and after execution. No producer source or data is modified.

PUBLICATION_FILES.json is the exact release allowlist, and SHA256SUMS binds every listed file. It includes the independent complete named bank, all card and exclusion records, source replay, and normal/optimized receipts. Binaries, temporary directories, and redundant raw remote fetch copies are excluded.
