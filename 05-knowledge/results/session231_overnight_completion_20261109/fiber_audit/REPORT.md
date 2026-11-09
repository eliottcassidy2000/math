> Later audit: ADDENDUM.md completes order-8/9 validation and weighted-loop checks. SHORT12_THEOREM_AUDIT.md proves the final unit-width construction. Read those before the historical boundaries below.

# Independent audit of Session 231 endpoint erasure and moving-core claims

## Outcome

PASS for the current PROOF.md Sections 1–5, with the scope qualifications below. No claim of unrestricted Chan–Pak completeness follows. The order-9 census and non-ordinal extension census were not independently rerun here.

The audit imports none of the producer's Python modules. `audit.py` counts Hamilton paths and ordered disjoint covering pairs of paths literally by permutations and splits, independently explores all moving-core edges, and independently verifies the complete order-6/7 representative banks by canonicalization and orbit mass. `physical.py` builds actual adjacency lists, counts all small star-closure Hamilton paths by forward subset DP, and checks strong connectivity from every vertex. Neither checker uses `assert` for substantive checks.

## All-length finite-state proof

For a physical core, 0 <= A,G <= q-h: A (respectively G) counts pairs of covering paths whose second path starts (respectively first path ends) in the cut; the q pairs include exactly h with the relevant path empty. Thus the asserted range is valid without any hidden positivity or coprimality assumption.

A move C,y to x,C' has the two equations

    2A + h*y = q*x + A',
    G + q*y = h*x + 2G'.

As h is odd and q even, integrality gives x = G mod 2. The audit instead tries all four bit pairs and tests the equations, avoiding dependence on that reduction.

Forward iteration of G yields a convex binary-weight expression with coefficients q and -h. From a physical initial G, it stays strictly in (-h,q). Backward iteration of A from a physical final A gives the same strict range. Consequently every actual path has integer states in {-h+1,...,q-1} squared, a finite set of (q+h-1)^2 states. This argument is for paths with both physical endpoints, not for arbitrary unrestricted trajectories.

For the common-context cancellation step, removing a common suffix recovers each suffix bit from A modulo 2 because h is odd. Removing a common prefix is justified by the exact interval 0 <= A_remainder <= q*2^b-h, which is strictly shorter than the prefix place value q*2^b. Removing singleton contexts also recovers G by its affine composition equation. The differing positions therefore reduce to the displayed finite graph. Conversely, each accepted bit-labelled path constructs the two actual ordinal deletion graphs; no intermediate bookkeeping pair is required to be physical.

Independent full graph results:

- B29: 64 physical pairs; 236 source-counted reachable states; 0 nontrivial physical-return ingress edges.
- B389: 240 physical pairs; 10,460 source-counted reachable states; 0 nontrivial physical-return ingress edges.
- B95: 96 physical pairs; 100,488 source-counted reachable states; 68 nontrivial physical-return ingress edges from 20 sources; shortest length 10.

All outgoing edges of every reachable state were checked, including destinations already visited. Each bank's compressed `.graph.jsonl.gz` records every source, its complete reachable state set, and every retained bit-labelled edge. Hashes covering states AND edges are in AUDIT.json. These hashes intentionally differ from the producer's state-only digest.

Important hostile control: the producer suppresses every endpoint-to-same-endpoint return when reporting its totals. The independent checker verifies that each suppressed ingress is exactly the zero-bit self-loop at (0,0) or the one-bit self-loop at (q-h,q-h). There are no hidden nonconstant endpoint return walks. This closes a potential gap in the stated all-length inference.

The B29 single-topology at-most-one-T and B389 at-most-two-T conclusions follow. Literal core counting additionally confirms B389's ambiguous A/G fibers have identical S,E. These are not conclusions for every topology with the same h.

## Physical witness checks

All small examples below were verified using literal path/split enumeration on the actual padded deletion graphs, forward Hamilton-path DP on the actual star closures, and graph reachability from every vertex. Full star adjacency lists are saved in PHYSICAL.json.

- H29 cross-topology unit switch: padded states (29,1280,898,287,174,29,0) and (29,1280,898,287,176,29,0); direct 10-vertex star counts 1651 and 1649. Independent canonicalization proves both cores and original B29 pairwise nonisomorphic.
- H45 five-value fiber: padded T values 164,166,168,172,174 and common (h,q,A,G,S,E)=(45,1728,1122,303,45,0); direct 10-vertex star counts 2125,2123,2121,2117,2115.
- H5 unsaturated pair: order-7 class codes/cuts (2050,50),(1114112,70), using row-major upper-triangle encoding as in `erase_classes.cpp`; padded common (h,q,A,G,S,E)=(5,1024,716,209,5,0), T=140,142. Direct strong 11-vertex star counts 1297,1295.
- H9 saturated triple: class codes/cuts (32772,49),(1056768,42),(262400,54); padded common (h,q,A,G,S,E)=(9,1344,908,245,9,0), T=152,154,156. Direct strong 11-vertex star counts 1673,1671,1669.

H95's four 31-vertex stars were explicitly constructed and their strong connectivity checked. Their seven-coordinate states agree exactly with PROOF.md, and their star counts calculated via the proved composition/closure laws are 5364196511,5364196519,5364196551,5364196559. Their count calculation is NOT a generic independent 31-vertex subset enumeration. Independent literal controls test singleton composition on either side of B95 for the extreme masks and all four participating masks. The full adjacency witnesses allow further independent graph-level verification if desired.

The common selected leading singleton guarantees strongness: r reaches it, it reaches every later ordinal vertex, every deletion vertex reaches v, and v reaches r. No inference from odd Hamilton count to strongness is used.

## Complete small census and exact H45 hole

For order 6 and 7, independent canonicalization enumerates all label permutations compatible with sorted degree classes. This is a genuine canonical invariant: every isomorphism preserves degree, and all degree-compatible bijections are included. Distinct canonical codes establish pairwise nonisomorphism; the number of canonical-minimum labelings gives the automorphism order. The exact orbit sums are:

- Order 6: 56 classes; sum 6!/|Aut| = 32768 = 2^15.
- Order 7: 456 classes; sum 7!/|Aut| = 2097152 = 2^21.

These equal all labeled tournaments, so completeness does not rest on the producer's individualization/refinement implementation or a remembered class-count sequence. Literal counting over every cut in these complete banks gives maximal consecutive step-two runs 3 overall and 2 in the gcd(h,q)=1 subuniverse at both orders. All cuts on one representative cover all labeled cut interfaces because relabeling preserves the seven counts.

At order 6, h=45,q=432,A=G=129 has exactly T+S in {35,37,39,43,45}; hence the particular padded target T=170 really is absent from the complete six-vertex bank. This is a fixed-order/fixed-interface hole, not a global missing scalar value theorem.

Minor producer receipt artifact: the order-7 census has 104432 nonempty fibers, while ERASE_CLASSES7.stats prints 104433. The final `fibers[key]` lookup inserts an empty H45 target key absent at odd order. It does not change the extremal run calculation, which occurs before that insertion. The same code should avoid reporting this artifact at other orders.

## Boundaries

- No all-walk weighted scalar classification for H95 or the complete H29/H45 banks was established. A shortest-path representative for each ingress does not cover repeated cycles, even though it is enough for reachability.
- The tested zero-drift loops remain examples, not a theorem on every cycle or every compatible context.
- No independent order-8/9 census or non-ordinal two-vertex extension census is claimed here.
- Fixed-position decoding requires fixed h,q and fixed split, but does not require gcd(h,q)=1.
- Endpoint erasure loses S,E only after preserving the induced T corrections. It does not turn arbitrary A/G coincidences into pure T switches automatically.

## Reproduction

    python3 -B session231-fiber-audit/audit.py
    python3 -O -B session231-fiber-audit/audit.py
    python3 -B session231-fiber-audit/physical.py

Run from /workspace/scratch/07d5532f2fce. Producer paths are explicit in audit.py. Normal and optimized audit summaries were compared. No producer file was modified, and no commit was made by this audit worker.
