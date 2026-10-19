# Independent review: linear routing in staircase partition systems

2026-10-19. **PASS: mathematical statement and constructive proof.**

Reviewed `../LINEAR_PARTITION_DIAMETER.md`, `../verify_partition_routes.py`, and the scope of `../arithmetic_partition_controls.py`. This is a review of the stated capacity-one partition subclass, not a general-matroid linear-diameter theorem or an optimal-distance claim.

## 1. Block reduction and strict Hall

The reduction can be verified directly without the general matroid path theorem. Use zero-indexed colors, supports S_i, and suffix unions U_k. If consecutive tight suffixes are at l<h, then |U_l|=r-l and |U_h|=r-h. Every full independent transversal saturates both suffixes, so colors l,...,h-1 occupy exactly B=U_l minus U_h. This applies to both endpoint states.

The effective support of color i in this block has size at least

    |S_i minus U_h| >= (r-i)-(r-h)=h-i.

For i=l this proves the full-buffer property. For a nonempty subset J of the other block colors, let j=min J. Its effective union has at least h-j parts, while |J|<=h-j. If that union had size at most |J|, both inequalities would be equalities, J would be the entire suffix j,...,h-1, and U_j would have size r-j. This would be an intervening tight suffix, a contradiction. Thus every such J satisfies strict Hall.

## 2. Availability-cycle reachability

After target relabeling, source assignments are a permutation sigma. If some cycles were unreachable from the buffer cycle, their union J would contain no buffer color. Any availability of a color in J at a part outside J would provide an incoming arc from a reachable cycle, making its cycle reachable. Hence the available parts of J would be contained in J, violating strict Hall. The spanning arborescence therefore exists with exactly the witnesses required by the algorithm.

The direction of this graph matters: a parent PART is usable by a child COLOR. The manuscript and implementation use the correct direction.

## 3. Traversal, final visits, and actual elements

For a child cycle entered at color j from parent part a, its entry swap parks j at a and exposes sigma(j). Each child subtree is processed before continuing the current cycle, and returns the buffer to exactly the incoming part. Moving color k into exposed part k then fixes its target assignment and exposes its untouched original part sigma(k). At the end, the entry color j moves from a to j and restores the buffer to a. Thus all colors in the completed cycle and descendants are fixed, and subsequent calls do not touch them.

This proves the claimed induction invariant. Singleton cycles also obey it: their color moves out and back, allowing child access and correcting an initially wrong parallel realization. The root performs the same traversal without an entry/exit pair. If the root cycle is a singleton, it still processes its child subtrees before finishing. This covers the potentially troublesome initially correct buffer.

Every nonbuffer color is visited, including initially correct singleton cycles. Choosing the exact target element on its final visit therefore corrects all parallel realizations. The buffer is set to its exact target element on its final return. If sigma is the identity, changing actual realizations in pairs changes no part assignments and is legal. Loops are never selected.

## 4. Budget and complexity

A nonroot cycle of length L costs L+1 swaps; the root costs L_root-1. With q nonroot cycles the total is b-1+q. A nonidentity permutation has at most b-1 cycles, so q<=b-2 and the route costs at most 2b-3. This bound includes the initially correct cycles deliberately disturbed by the route. The identity-assignment case costs at most ceil(b/2), which satisfies the displayed block bounds.

Summing the block bounds gives the stated 2r-3 bound for r>=3. At r=1,2 the whole target can be installed in one width-two move, even when the sum of singleton-block bounds is larger. No optimality statement follows from these route counts.

With explicit part identifiers, support/representative extraction costs O(N), suffix set work and all cycle-availability scans cost O(r^2), and there are O(r) moves. Materializing each full r-coordinate state therefore also costs O(r^2). The stated O(N+r^2) time/storage upper bound is valid in this representation. The bounded-size Hall debug enumeration does not change the asymptotic bound. The producer test driver's state enumeration is outside the routing algorithm.

Two implementation issues found during review were corrected: recursive DFS was replaced by an explicit event stack, and the rank<=2 early return now increments the aggregate MOVES counter. The stack order was independently inspected: each cycle position executes its hole check, then child subtrees, then its cycle move, and the total pending event count is O(b). The producer's rank-1,100 relay control completes 2,197 moves at arborescence depth 1,098, beyond the ordinary Python recursion ceiling. Its note correctly distinguishes this long constructed route from the actual width-two distance, which is one for that stress endpoint pair.

## 5. Graphic/arithmetic scope

A connected multigraph whose simplification is a tree has a capacity-one partition graphic matroid, with one part per simple tree link. Parallel contacts remain distinct actual elements, and contracted loops are unavailable. The inherited arithmetic rank staircase therefore gives the corollary exactly under the stated whole-contact-graph condition.

The bounded arithmetic receipt reports no nonstar or noncomplete-support examples; every recognized fiber in that exhaustive range has free rank one or two. Accordingly, strict extension beyond the complete factorial-assignment family is currently demonstrated for abstract partition systems, not by a new arithmetic example in that census. The final producer control additionally constructs 120 endpoint routes in the already proved factorial family through rank twelve. The corollary itself is unconditional under its explicit contact-tree hypothesis.

## 6. Independent finite verification

`independent_partition_audit.py` imports no producer code. It independently reconstructs suffix blocks, verifies strict Hall for every subset, constructs a differently ordered availability arborescence, checks every swap and child return, and verifies exact actual-element targets and route budgets. It also builds the complete width-two state graph and computes all-source BFS distances.

Exact universes:

- All 857 labelled support systems satisfying the staircase through rank four, with all 45,643 endpoint pairs.
- All 32 such support systems through rank three, with two distinguishable actual realizations at every available color/part pair and one loop per color, with all 14,692 actual endpoint pairs.
- In total, 60,335 independently constructed routes, 193,916 validated moves, and 1,880,058 explicit checks.

Maximum exact assignment-state diameters at ranks 1,2,3,4 are 0,1,2,4. Maximum exact doubled-realization diameters through rank three are 1,1,2. These are finite controls only, not sharpness claims at larger ranks.

Reproduce from the packet directory:

    python3 -B review-algebra/independent_partition_audit.py > review-algebra/independent_receipt.json
    python3 -B -O review-algebra/independent_partition_audit.py > review-algebra/independent_receipt.optimized.json
    cmp review-algebra/independent_receipt.json review-algebra/independent_receipt.optimized.json

Normal and optimized receipts agree byte-for-byte. The universal claim rests on the proof audited above.

## Frozen reviewed proof and implementation

- `LINEAR_PARTITION_DIAMETER.md`: SHA256 `9adbbf668476ef42d00d681f1543eab497750d32be06f2d74504f20383591f6b`.
- `verify_partition_routes.py`: SHA256 `91b1aa356816c764bb4e21accb8d970b700ef258cfb4dc06d737bca3a164d171`.

Both implementation fixes and the explicit arithmetic-scope clarification are included in these reviewed bytes.
