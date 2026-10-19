# Independent review: pendant-pair parity reductions

2026-10-19. **PASS** for the two-cherry theorem, its named counterexample, the eight-vertex padding theorem, and the complete order-nine classification. The D8 argument below was independently proved before the final author exposition was checked in Section 5.

Reviewed artifacts: [TWO_CHERRY_PARITY.md](../TWO_CHERRY_PARITY.md) and [NINE_VERTEX_GOOD_TREES.md](../NINE_VERTEX_GOOD_TREES.md) (final SHA256 `dc21a674f2998b6c5dc797ef9cc2e4ceec2fcc68f33f2b1ce1ff7773d31a26fb`). This review supplements, and does not replace, the earlier core review in `INDEPENDENT_REVIEW.md`.

## 1. Two pendant pairs

Deletion–reversal on one source edge is an exact integer identity for injective embeddings into a tournament: each injection preserving the remaining arcs satisfies exactly one orientation of that edge. The identity remains valid when the source is not spanning the host.

For two pendant leaves at one center, equal directions give a transposition preserving every source arc. For opposite directions, reversing one edge gives that equal-direction transposition, so modulo two the original graph can be replaced by the graph with that edge deleted. The resulting isolated vertex must be retained.

With two disjoint pendant pairs at distinct centers, either a pair has equal directions already, or these reductions can be performed successively. The second pair is unaffected by the first deletion. The final source has two isolated vertices, whose transposition acts freely on injections. The universal-even conclusion is valid for every tournament host and arbitrary remaining source graph structure.

For the author's tree with edges 01,02,03,04,15,16,27,78, deleting vertices 4 and 6 leaves branches of lengths 1,2,3 at 0. This is an asymmetric subtree. The underlying automorphism group has order four, consisting of the two independent leaf swaps. Direct permutation enumeration independently confirms 64 asymmetric labelled orientations, hence 16 asymmetric orientation classes. Thus the tree is good outside both hereditary symmetry and the automorphism-forced-zero class. Minimal order nine follows from the separately verified complete nonconstancy certificates through order eight.

The named mask 40, reversing 04 and 16 relative to increasing orientation, has 260 embeddings into the transitive tournament and exactly one into its own directed host. These values were independently reproduced. The author's complete two-cherry census through order eleven was rerun in both normal and optimized Python; both outputs match. Its interpretation as filter counts, not a goodness classification, is correct.

## 2. Independent proof of the D8 plus isolate extension

Let D8 have edges 01,02,03,14,15,26,47, and retain a ninth isolated source vertex. Deleting a leaf edge gives two isolated source vertices, hence zero embedding parity. Deleting either edge 02 or 14 gives component sizes 6,2,1. Deleting 01 gives sizes 4,4,1.

All oriented trees on at most six vertices have constant tournament parity, by the inherited hereditary-symmetry theorem together with the exhaustively checked absence of asymmetric subtrees at these orders. Consequently, for any fixed orientations of the deleted-edge components, their forest embedding parity is a fixed product of the component parity constants times the number of ordered image-set partitions. The corresponding multinomial factors are 9!/(6!2!)=252 and 9!/(4!4!)=630, both even. This argument does not assume the two four-vertex component orientations coincide.

Deletion–reversal therefore shows that reversing any edge of D8 does not change the parity function of D8 plus an isolate. The orientation hypercube is connected, so it suffices to evaluate one anchor.

For that anchor, orient each rooted P4 half identically as a directed path, with the root at the same internal position; orient the connecting root edge arbitrarily. Fix the isolate image. Partition the remaining eight host vertices into an unordered pair of four-sets. For each pair of embeddings of the directed P4 halves, the two assignments of those identical rooted copies to the four-sets reverse the central-edge requirement. Exactly one assignment is allowed. Each four-set admits an odd number of directed Hamiltonian paths, so each unordered partition contributes odd parity. There are (1/2) binom(8,4)=35 such partitions. Finally there are nine choices of isolate image. Thus the anchor, and every orientation of D8 plus an isolate, has odd embedding count in every nine-tournament.

Add a new leaf at a D8 vertex already adjacent to a leaf. Equal directions on the resulting pair give zero. Opposite directions reduce to the corresponding D8 plus isolate and give one. Up to underlying tree isomorphism there are two choices: attach at 0 or at 2. Both resulting nine-vertex trees contain S(1,2,3), and both have automorphism group of order two. Each therefore has 128 odd labelled orientations and 64 odd orientation-isomorphism classes. Direct independent permutation and linear-extension calculations confirm these numbers.

This does not say D8 itself is good at order eight. The isolate extension changes the source and host order and is essential.

## 3. Exact structural reformulation

More generally, let H have a pendant leaf adjacent to v, and form U by adding another leaf adjacent to v. Then U is good if and only if H disjoint union an isolated vertex is good, both evaluated at order |H|+1. The equal-direction U orientations vanish; the opposite-direction orientations are in bijection with orientations of H and have precisely the corresponding H-plus-isolate parity function by deletion–reversal. This equivalence retains the isolate and does not assert heredity of goodness.

Any constant spanning n-source parity c also determines its injective count into every larger N-tournament: the parity is c binom(N,n), by partitioning embeddings according to their n-element image. In particular, the D8-plus-isolate odd statement extends with this binomial factor, not with an unqualified odd claim for all larger host orders.

## 4. Independent controls and scope

`independent_cherry_extensions.py` imports only the independent permutation-based tree utilities from this review directory, not the author's implementation. It performs 1,169 explicit checks using direct degree-class permutation automorphisms and an independently written subset-DP linear-extension counter for transitive tournaments. The normal and optimized outputs are byte-identical in `cherry_extension_receipt.json` and `cherry_extension_receipt.optimized.json`.

Checks include all 256 orientations of the two-cherry counterexample, the exact deletion/reversal/swap certificate for its 64 asymmetric orientations, all 128 orientations of D8 plus an isolate, and all 256 orientations of each of the two leaf extensions. Finite transitive-host counts are controls only; the universal tournament quantifiers follow from the proofs above.

No files inside the shared mathematical checkout were modified.

## 5. Final audit of the complete order-nine classification

The final [NINE_VERTEX_GOOD_TREES.md](../NINE_VERTEX_GOOD_TREES.md) has now been read directly. Its eight-vertex padding theorem, universal affine-orientation proposition, three-part structural criterion, critical-leaf selector, and two concrete odd-orientation examples all pass mathematical review.

The padding theorem correctly handles every edge-cut size. For symmetric eight-trees, a nontrivial vertex-fixing rooted automorphism group supplies an involution by the core branch recursion. Such an involution reverses no edge, so an invariant orientation exists and is universally even. If no such automorphism exists, the only symmetry exchanges two rooted four-vertex halves, both rigid as rooted trees. Those halves must be P4, rooted either at an endpoint or internally. Hence precisely P8 and D8 have odd padded parity; every other symmetric eight-tree has even padded parity. For an asymmetric eight-tree, adjoining one unique isolate retains asymmetry and the highest-degree seven-edge coefficient, so the common parity function is nonconstant.

The delicate necessity step for a 2+7 cut is valid. In the embedding polynomial of an asymmetric seven-tree A disjoint union P2, fix the monomial supported on a connected six-edge copy of A. Any six-edge source selection including the separate P2 edge is disconnected; injectivity prevents it from producing that connected monomial. Thus the selected edges must be exactly those of A. Its asymmetric underlying structure gives a unique map onto the selected host copy, and the two orders of the P2 endpoints contribute opposite literals whose sum is one. The coefficient is therefore one for every source orientation. This proves nonconstancy without misapplying the highest-degree obstruction to the entire disconnected forest.

For sufficiency, each one-edge deletion is rigid by the padding theorem or by the smaller good components and even shuffle factors. A symmetric odd-order tree has a vertex-fixing involution and hence an invariant zero anchor. Edge flips add constants, so all orientations are rigid. No converse stronger than the order-nine statement is assumed.

The universal affine-orientation proposition also passes. Every two-edge deletion gives three positive component orders summing to nine. A component of order seven forces two isolates. Otherwise all components have order at most six and are good, while their multinomial coefficient is even: three positive binary-disjoint summands would require at least three set bits in their sum, but nine has only two. Every mixed second orientation derivative therefore vanishes as a function of the tournament. This implies affineness for every nine-tree, including nongood ones, and does not wrongly turn the possibly nonconstant coefficient functions into constants.

For a good tree, non-leaf derivatives vanish by even shuffle factors. A leaf derivative equals one precisely when its eight-tree is P8 or D8. The critical-leaf parity formula relative to a zero invariant anchor follows and is independent of the choice of that anchor.

`independent_nine_criterion.py` independently regenerates all 47 unlabelled nine-trees using all-root canonicalization and computes automorphism orders by direct degree-class permutations. It confirms 36 good trees under the proved criterion, 24 hereditary trees, and 12 additional good trees. It independently checks orientation independence on all 23 eight-trees in the transitive host, identifies 22 symmetric eight-trees and exactly two odd padded types, and exhausts all 512 orientation/host pairs for the P4 input. Normal and optimized receipts agree.

The author's final verifier also reproduces its receipt in normal and optimized modes, including all 1,316 two-edge deletions, all 9,216 good-tree orientation-selector cases, and the 3+3+4+2 decomposition of the twelve nonhereditary good shapes. Explicit bad-host certificates are separately checked by `independent_nine_witnesses.py`, which uses all 9! permutations rather than the author's tree-embedding backtracker.

**Conclusion:** the complete 36-of-47 order-nine classification and its orientation selectors are approved within their stated order and inherited theorem dependencies. The earlier statement that all goodness classification above eight remains open is superseded only at order nine; the unrestricted classification remains open.
