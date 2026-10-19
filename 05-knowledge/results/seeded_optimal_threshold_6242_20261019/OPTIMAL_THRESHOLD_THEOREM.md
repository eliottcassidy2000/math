# Exact threshold6242 for square-sum cycles containing482–3362–6242

## Theorem

For an integer N, a Hamiltonian cycle on {1,…,N} with square-sum adjacent pairs can contain482,3362,6242 consecutively if and only if N≥6242.

The cycle may be oriented so that the three vertices occur in the displayed order. This is a constructive computer-assisted theorem. No claim of historical novelty or priority is made.

## Proof

Necessity is immediate: the required vertex6242 is absent when N<6242.

For sufficiency, the sibling bundle `../seeded_eventual_cycle_theorem_20261019/` proves the result for every N≥128453, using checked nice-pair bases,294 universal pair transitions and392 universal final split-lift states.

The present bundle independently certifies every remaining integer6242≤N≤128452. It starts with an explicit6242-cycle and supplies one exact extension record for each of the122210 subsequent lengths. The replay checker validates the initial cycle and proves the validity of every intermediate cycle by the invariant below. Thus there is no gap between the finite interval and the eventual theorem. This proves sufficiency and the exact threshold. □

## The local certificate invariant

Represent a cycle by a single ordered list q0,…,q_(n−1), with the last vertex followed by the first. Keep the two undirected protected edges {482,3362} and {3362,6242}.

An extension record first inserts the new vertex n+1 after a declared position k. This removes only the boundary edge {q_k,q_(k+1)} and adds the two edges through the new vertex. The checker refuses to remove a protected edge. Since the new vertex is the next integer and all old entries are retained, the new sequence is still a permutation of[1,n+1] and represents one cyclic order.

Each following repair specifies two distinct cut positions i<j. Reverse exactly the contiguous sublist with indices i+1,…,j. If the four boundary vertices before reversal are a=q_i, b=q_(i+1), c=q_j, d=q_(j+1), with cyclic indexing, this removes the two undirected edges {a,b},{c,d} and adds {a,c},{b,d}. Every internal edge of the reversed sublist remains the same undirected edge. Reversing a sublist preserves the permutation and, importantly, a single cyclic order; it cannot create two disconnected cycles. The checker again refuses to remove either protected edge.

The checker maintains the exact set of nonsquare boundary edges. It starts empty because the initial cycle is exhaustively checked. Insertion and reversal change only the edges just listed, so deleting their old nonsquare edges and adding their new nonsquare edges preserves exactness of this set by induction. Every edge test uses integer square root. The trace happens to maintain at most one nonsquare edge during a repair, and this additional property is checked too.

At the end of every extension record the bad-edge set must be empty. Therefore the current ordering is a square-sum Hamiltonian cycle on the next integer interval, and both seed edges remain present. Since3362 has exactly two cyclic neighbors, the preserved edges force the consecutive seed, up to reversal of the cycle's orientation.

This is a certificate lemma for the supplied moves. It does not assert that every seeded cycle admits such an extension or that the randomized search always succeeds. Existence here is proved by replaying every supplied finite record.

## Efficient exact replay

The standard-library checker represents the single sequence by an implicit binary tree. Subsequence split, merge and lazy reversal implement exactly the list operations above, so no orientation inference from an unordered adjacency graph is required. Tree priorities affect only the representation; they do not change any accepted trace operation or square test.

For extra checks beyond the induction, replay exhaustively traverses and checks every cycle through6500, every subsequent8192-vertex checkpoint, and the final128452-cycle. Together these are273 full intermediate audits, plus initial and final audits. The final sequence is also compared with the independently generated stored final cycle.

A second, independent C++ checker replays the same trace using ordinary explicit vector insertion and reversal, without the implicit tree. Its passing receipt agrees on all122211 lengths,504321 moves and the exact final sequence. This is an additional implementation cross-check; the standard-library Python replay remains sufficient.

The entire finite interval passes7147854 explicit checks in both ordinary and optimized Python execution. The receipts are identical. Checks use explicit failures, not removable assertions.

## Frozen evidence

- `initial_cycle.txt`: the complete6242-cycle. Its initial list begins482,3362,6242. Replay independently checks interval coverage, every cyclic square edge and both seed edges.
- `full_incremental_trace.bin.gz`: the complete extension trace,2373041 compressed bytes. It encodes122210 insertions and504321 reversals. There are19548 direct insertions requiring no reversal; the maximum number of reversals in one extension is58.
- `verify_full_incremental_trace.py`: the exact independent standard-library checker. It reads the compressed trace directly.
- `final_cycle.txt.gz`: the resulting128452-cycle, for independent final-sequence agreement.
- `full_incremental_controls.json` and its `.optimized.json` companion: passing replay receipts, including source, trace and final-cycle hashes and all checkpoint cycle hashes.
- `full_search_report.json` and `extend_seed_cycle.cpp`: optional discovery metadata and source. Search took about17 seconds in the recorded environment. Search speed and solver trust are not assumptions of the proof.
- `RELEASE_MANIFEST.json` and `SHA256SUMS`: release inventory, dependency and file-integrity records.

The decompressed initial sequence has SHA-256 `a892ae08c5304c35fae319dfcd9429faa5ea37b2f21e89d3a715e17912ceca6d`. The decompressed final sequence has SHA-256 `220f24a78945369fc9a485cd047952e8714f17aab68649b35d755e9509a05509`.

## Binary trace format

All integers are unsigned little-endian. The decompressed file begins with the eight bytes `SSINC1\0\0`, followed by32-bit first and last cycle lengths. For each successive N from first+1 to last, there is:

1. a32-bit insertion index into the current length N−1 list;
2. a16-bit number of reversals;
3. that many pairs of32-bit indices(i,j), interpreted in the current length N list immediately before each reversal.

Every operation index is checked. The checker requires all declared records and rejects trailing bytes. No intermediate cycle is omitted: the next length is implicit in the record's position.

## Reproduction

From this folder run:

    python verify_full_incremental_trace.py > replay.json
    python -O verify_full_incremental_trace.py > replay.optimized.json
    cmp replay.json replay.optimized.json
    cmp replay.json full_incremental_controls.json
    sha256sum -c SHA256SUMS

The finite-gap checker is self-contained. The all-N conclusion additionally uses the sibling eventual-theorem bundle, whose own checker and dependencies replay its infinite induction. Neither package modifies the other.
