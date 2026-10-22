# Independent final audit: fixed-core attachment-word counter

Verdict: PASS. No blocking mathematical, indexing, multiplicity, or complexity defect found. This verdict applies to the exact source hashes recorded in `check.out`, not to subsequently changed files. The producer's existing SHA256SUMS was verified; the review also hashes the witness and both producer outputs, which the producer manifest did not include in full.

## Mathematical review

A Hamiltonian vertex word uniquely determines its core permutation and the slot containing each distinct base label. Within one slot base labels must increase, and an assignment reconstructs exactly one candidate word. This is an injective and surjective encoding of candidates, not an assumption that equal attachment letters identify vertices. Even when masks repeat, the positions are distinct vertices and are independently assigned; every transition consumes one particular position.

For a noninitial occupied slot, only its first vertex meets the preceding core vertex. Thus the first-port condition is imposed exactly once, when the slot becomes occupied. Only its last vertex meets the following core vertex. Appending a larger base label preserves all internal base edges and overwrites that exit port. Retaining temporarily bad exits is essential and correctly implemented. For a singleton interior slot both entry and exit conditions apply to the same vertex; the implementation does not accidentally treat them as independent vertices.

Empty interior slots contribute an actual consecutive core arc and must be tested against Q. Nonempty slots do not require that core arc. The leading empty slot needs no incoming edge, and the terminal slot has no exit condition. These are exactly the final predicate and transition rules in counter.py. For m=0, only the empty state survives, so acceptance sums the directed Hamiltonian core permutations. For k=0, validation still enforces masks equal to zero, and the public counting API returns one, including the empty graph. The internal cached transitions helper is used only for positive k; its lack of a k=0 transition table is not a public counting defect.

Each eligible slot choice has multiplicity one, but different choices can reach the same state. The implementation stores repeated target indices rather than deduplicating. A direct hostile checks that for k=1, mask zero, state four, the row is (4,4). This distinguishes weighted counting from Boolean reachability. A separate direct hostile checks repair of the bad exit: k=1, mask zero, state two moves to state one.

## State and computational claims

For each fixed core permutation there are at most 2*3^k states, hence at most 2*k!*3^k in the disjoint-block realization. This is an upper bound on states of a nonnegative-integer weighted automaton, not a minimal deterministic-automaton count and not a lower bound. The machine may contain unused/unreachable states without harming the upper bound.

The stated O(m*k!*(k+1)*3^k) propagation additions and separate O(k!*k*3^k) final checks are valid. The implementation's local-mask conversion adds O(k) per letter and permutation, within the stated propagation bound. The dense cached table has O(2^k*(k+1)*3^k) parameter-only preparation/storage in the stated indexing/arithmetic model. Sequential permutation processing allows O(3^k) live count entries, beyond the cached table and materialized word. Count bit length is bounded by log2(k!)+m*log2(k+1)+O(1), because all partial and final counts count assignments. The proof correctly excludes zero-cost m=0 processing and constant-bit-time implications. It does not solve presentation discovery, minimum transitive deletion, canonicalization, or unrestricted tournament/A000568 enumeration.

## Independent reusable controls

`check.py` imports only the producer counter, not its audit helpers. Its oracle directly enumerates admissible vertex words by recursive extension with no memoization or subset DP. It defines the edge relation directly from vertex types, core matrix and word positions. Exhaustive oracle ranges are k=0,m<=6; k=1,m<=6; k=2,m<=4; k=3,m<=3; k=4,m<=2, over every labelled core and every word in these ranges.

Results:

- 22,968 independent backtracking-oracle comparisons.
- 3,598 simultaneous core relabelling/word-bit permutation comparisons. These exhaust all core permutations for the tested k<=3,m<=2 subrange, including k=0.
- 22,968 reversal-duality comparisons. Reverse core arcs, reverse the attachment word, and complement its bits. This reverses all graph arcs after reversing the transitive base order and must preserve the Hamiltonian count.
- 685 explicit partial-slot-assignment mass comparisons, for k=1,2,3 and all words through length three. This checks the complete intermediate state distribution, including invalid exits retained for later repair, first-port filtering and colliding destinations, not merely the accepted total.
- 11 malformed-input rejections, plus acceptance of Boolean core entries and a generator word.
- Direct singleton/first-port, temporary bad-exit, repeated-target, empty-core-word and cyclic-core controls.
- H13 witness Q=[[0,1],[0,0]], word=[1,2,0] yields 13.
- Holding word=[1] fixed while reversing only the core arc yields 1 versus 3. Correct simultaneous relabelling also permutes attachment bits; core path-cover-profile equality is not enough for fixed-word equality.

The producer's source audit was also read. Its reported 25,146 generic subset-DP comparisons, 244 exact means and 8,191 one-apex recurrence controls have internally consistent range totals. Its path-cover computation enumerates core orders and allowable cuts, and exact inclusion-exclusion handles q=0,m=0 via 0**0=1, as required. The expectation theorem and inherited one-apex recurrence are dependencies already independently passed upstream; this review does not relabel the producer's 244 mean checks as newly independent theorem verification. Producer normal and optimized outputs are byte-identical. No producer script was rerun or source modified by this reviewer.

## Reproduction and source binding

Run from any directory:

    python3 /workspace/shared/session26-counter-review/check.py
    python3 -O /workspace/shared/session26-counter-review/check.py

The check reads the fixed source directory `/workspace/shared/session26-counter/`. If distributing elsewhere, change only S in check.py to the copied producer source directory, then refresh the review manifest. Any source or witness change is visible in the output hashes. All failures use explicit exceptions and remain active under Python optimization.

The review is confined to this final bounded instrument. No new research lane, repository write or Git operation was undertaken.
