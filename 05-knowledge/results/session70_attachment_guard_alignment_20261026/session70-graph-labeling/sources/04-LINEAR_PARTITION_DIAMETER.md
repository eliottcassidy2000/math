# A cycle traversal gives linear diameter for partition-matroid upper-tail fibers

2026-10-19. **PROVED, scoped; independently [reviewed](review-algebra/REVIEW.md).** This improves the quadratic bound on a structural subclass and supplies an explicit polynomial routing algorithm for that subclass. It makes no sharpness claim and no linear-diameter claim for arbitrary matroids.

## 1. Scope and the retained actual elements

A rank-r partition matroid in this note has r nonloop parts and capacity one in each part; loops may also be present. An actual element belongs to one part and to one color class. Different actual elements in the same part remain distinguishable, including different realizations of the same color.

Let A_1,...,A_r be disjoint color classes in such a matroid, with

    rank(A_i) >= r-i+1.

Assume their union has rank r. A state is an independent choice of one actual element per color, hence a basis meeting every part exactly once. Allowed moves change at most two actual color coordinates. Every such legal move preserves the full basis closure.

The preceding [polynomial rank-staircase theorem](../graceful_polynomial_diameter_20261019/POLYNOMIAL_STAIRCASE_DIAMETER.md) proves that all tight subsets are suffixes. Consecutive tight suffixes give blocks. In the quotient by the inner tight flat, a block of b colors has b parts, its first color meets every part, and every nonempty subset J of its remaining colors meets strictly more than |J| parts. These statements follow directly from the rank staircase and do not use the cited general-matroid path theorem. Contraction and restriction preserve the capacity-one partition structure.

The new step is a linear route inside each such block.

## 2. The full-buffer block lemma

**Lemma 1.** Suppose b color classes have an independent transversal meeting b capacity-one parts. One distinguished buffer color meets every part, and every nonempty subset J of the other colors meets more than |J| parts. Then any two states can be joined by at most

    f(1)=1;
    f(2)=1;
    f(b)=2b-3 for b>=3

moves changing at most two actual coordinates. The route is constructive, including exact parallel realizations.

**Proof.** Relabel the colors and parts as 0,...,b-1 so color 0 is the buffer and the target element of color i belongs to part i. The source assignment is a permutation sigma: its color i occupies part sigma(i). Regard the part occupied by color 0 as a movable hole; a move of another color into that hole is implemented by swapping that color's assignment with color 0, choosing actual elements in their new parts.

If sigma is the identity, only actual realizations may differ. Change them to target in pairs, using at most ceil(b/2) moves. For b=1,2 this proves the stated bound; for b>=3 it is at most 2b-3. The rest assumes sigma is not the identity, which entails b>=2.

Partition the indices into cycles of sigma. Form a directed graph on these cycles. There is an arc C to D when some part a in C is available to a nonbuffer color j in D. Record one such witness (a,j) for each chosen arc. The root is the cycle containing color 0.

Every cycle is reachable from the root. Otherwise let J be the union of the indices in the unreachable cycles. It excludes color 0. No color in J can use a part outside J, because every outside cycle is reachable and that availability would give an arc into an unreachable cycle. Hence the union of J's classes meets at most |J| parts, contrary to strict Hall. Choose a rooted spanning arborescence in this cycle graph, recording its witnesses.

Here is an explicit depth-first traversal. When the hole is at a part a in a parent cycle, enter a child cycle through its witness color j: move j into a. Its original part sigma(j) becomes the hole. Now follow the sequence

    sigma(j), sigma^2(j), ..., j.

Whenever the hole reaches a part k, first process every child cycle whose entry witness is that part. Each such recursive call returns the hole to k. If k is not j, move color k into its target part k; the hole becomes sigma(k). When k=j, move the temporarily parked color j from the parent part a into its target part j. The hole returns to a and all colors of this cycle and its descendants are at target.

This child-cycle routine uses |C|+1 swaps for a cycle C, apart from its descendants. A one-element cycle is allowed: its token moves out and back, providing access to descendants. The root routine starts with the hole at sigma(0), traverses its cycle in the same way, and finishes with the hole at part 0. It uses |C_root|-1 swaps, because color 0 itself need not be moved into a parent part.

Thus, if q is the number of nonroot cycles, the total number of swaps is

    b-1+q.

For a nonidentity permutation, the total number of cycles is at most b-1, so q<=b-2. The route uses at most 2b-3 moves. At b=2 this is one move.

Actual elements cause no extra steps. A temporarily parked color uses any available realization in its temporary part. Whenever a color enters its target part for the last time, choose its exact target element. The buffer uses an arbitrary available realization in each intermediate hole and its exact target realization on its final return to part 0. Every nonbuffer color is visited, including initially correct singleton cycles. Consequently the route reaches the full actual-element target, not merely its assignment quotient. Every intermediate choice occupies the b parts exactly once, so independence and closure are preserved. QED.

The displayed bound counts a particular buffer-only route. It is not claimed to be the true distance. A target transposition between two nonbuffer colors can often be done directly in one allowed move even when a buffer-only route needs a long relay.

## 3. Linear staircase theorem

**Theorem 2.** In the fully tight rank-staircase partition system of Section 1, every two states have a path of length at most the sum of f(b) over consecutive tight-suffix block sizes b. Uniformly, the diameter is at most

    0 if r=0;
    1 if r=1 or r=2;
    2r-3 if r>=3.

**Proof.** Process the blocks from inside out. The already completed inner flat is spanned by its fixed target coordinates. Lemma 1 applies in each quotient block. Every block route remains a basis of the same fixed quotient flat, so adjoining unchanged inner and outer coordinates preserves full independence and the global closure. Earlier completed blocks are never changed.

For r>=3, let q1 count size-one blocks, q2 count size-two blocks, and q3 count larger blocks. If there is a larger block, then

    sum f(b) = 2r-q1-3q2-3q3 <= 2r-3.

If all blocks have size one or two, their sum is at most r<=2r-3. At r=1 or2, replace all differing coordinates at once. QED.

This is genuinely broader than the complete factorial assignment family: class supports may be unequal, source-to-target cycles may require relay colors, and initially correct coordinates may move. The existing five-state target-preservation obstruction is a partition-matroid system and is included. No monotone-target hypothesis has been added.

## 4. Arithmetic specialization

Fix an acyclic upper tail F containing one edge of every difference color r+1,...,m on labels 0,...,m. Contract F. Form the multigraph of all remaining actual difference-d edges for d=1,...,r, retaining their color and original endpoints. Remove loops and simplify parallel edges.

**Corollary 3.** If this simplified contact graph is a tree, the graceful completion fiber has diameter at most 2r-3 for r>=3, with the small-r bounds above.

Indeed, the graphic matroid of a multigraph whose simplification is a tree is exactly a capacity-one partition matroid: its parts are the parallel classes corresponding to the simple tree edges. The earlier arithmetic residue-chain argument gives rank(A_d)>=r+1-d, and color one spans the quotient. Theorem 2 therefore applies. Actual parallel contacts are retained throughout the route.

The condition is on the entire available contact graph, not merely on one chosen completion, which is a tree in every graceful fiber. If the simplified contact graph has a cycle, this corollary supplies no linear bound. The universal quadratic theorem still applies there. No claim is made that every prescribed partition support system is realized by integer-difference classes.

## 5. Algorithmic complexity

For an explicitly given partition system, let N be the total number of actual input elements, including loops and repeated parallel realizations. Part identifiers and actual-element indices are the input representation. The route can be found in O(N+r^2) elementary operations, with O(N+r^2) storage when full intermediate states are materialized.

Compute color supports and one representative actual element for every available color/part pair in one input scan. Compute suffix supports and identify all tight suffixes; straightforward set unions use O(N+r^2) work. In each block, label parts by target colors, compute permutation cycles, and scan at most b^2 color/part pairs to construct the cycle arborescence. Its depth-first route has O(b) moves. An explicit event stack, rather than language-level recursion, implements the traversal. Representative lookup makes each move constant work aside from copying an r-coordinate output state. Across all blocks, the scans and copied output use O(N+r^2) work. Exact target actual-element indices are installed on final visits.

The supplied `route` function implements this cycle construction and never searches the state graph. It additionally audits strict Hall by subset enumeration only on blocks of size at most eight; this fixed-size debug check is not an exponential subroutine in r. The surrounding test driver enumerates some small state sets to test endpoint coverage, which is separate from the routing algorithm.

This complexity statement assumes the partition representation is given. It is not a polynomial-oracle algorithm for arbitrary matroids and does not replace the open implementation task for the general quadratic theorem.

## 6. Controls and status

The companion `verify_partition_routes.py` supplies the route and finite tests, including distinguishable parallel realizations and loops, complete endpoint pairs in bounded small systems, and seeded larger systems. `arithmetic_partition_controls.py` independently enumerates bounded upper-tail choices, recognizes the required contact condition, and tests actual-element routes. Run both normally and with Python optimization; matching receipts are retained.

The mathematical proof, not finite extrapolation, establishes the universal subclass bounds. Independent [review has passed](review-algebra/REVIEW.md). No assertion of optimal diameter or historical priority is made.

The exhaustive arithmetic range m<=10 and r<=5 contained 79,446 recognized tree-contact fibers and 100,890 tested endpoint routes. No nonstar contact tree or unequal color-support example appeared in that range. Thus the strict extension beyond the complete factorial family is certified for abstract partition systems; the arithmetic statement remains the proved conditional corollary, with no claim that a new nonstar integer-difference realization was found. All recognized fibers in that exhaustive range had r=1 or r=2; 120 constructed members/routes from the existing factorial family, through rank 12, provide higher-rank arithmetic controls.

The main routing receipt records 2,737,046 explicit checks, 71,738 endpoint routes across 960 systems through rank 16, and 355,829 validated moves (including small-r moves). A separate rank-1,100 event-stack stress reaches a cycle-arborescence depth of 1,098 and emits 2,197 legal buffer swaps. That stress has full width-two distance one and is explicitly not a diameter lower bound. Each of the three producer control suites has a matching normal/optimized receipt.
