# Induced H25 to H51 requires rebuilding

**PROVED, independently audited in the stated scope.** No tournament with 51 directed Hamiltonian paths contains an induced subtournament with 25 directed Hamiltonian paths. The orders, strong-component structure, number of births, and source family are unrestricted.

A specific six-vertex H25 parent reaches H51 by reversing exactly one old arc and adding one vertex. Thus one old-arc reversal is the exact minimum for this parent and target when births are allowed. The proof and independent audit are in [REVIEW.md](session203-induced51-audit/REVIEW.md).

This is an ancestry obstruction and a physical repair, not a missing-count result. It does not rule out eventual binary induction starting above 25, and it proves neither an eventual odd-count threshold nor a rebuilding operation for general parents.

## Explicit positive repair

Vertices are numbered from zero. A row is an integer bitmask: bit j in row i is one exactly when i points to j.

- Parent: (8,13,17,20,35,15), with H=25 and two Hamiltonian cycles.
- Reverse the old arc 0 to 3. Intermediate: (0,13,17,21,35,15), with H=11 and no Hamiltonian cycle.
- Add vertex 6 beating only old vertex 1. Child: (64,13,81,85,99,79,2), with H=51 and three Hamiltonian cycles.

All other old arcs stay unchanged. The intermediate is weak, which is allowed in this repair. Hamiltonian cycles are directed cyclic words modulo rotation.

## Why the all-order statement has a finite certificate

SCC multiplicativity and induced monotonicity reduce a possible H25 subset of an H51 tournament to one strong H51 component. The alternative factorization 51=3 times 17 cannot supply an induced count 25. The CITED Busch lower bound confines a strong H51 component to order at most 8. The complete bank through 8 has exactly 12 strong H51 classes, and none contains an H25 subset.

The bank is complete by replayed exhaustive one-vertex growth and full-adjacency canonicalization, not merely by a class-count comparison. Independent ending-at-vertex subset DP checks all 7,412 representatives. Literal simple-path DFS independently agrees on all 6,144 subset counts of the 28 H51 representatives. Literal permutations verify the displayed repair.

## Reproduction

Run from this directory:

    g++ -O2 -std=c++17 session194-weak-parent-certificate/verify_completeness.cpp -o /tmp/s203_induced_growth
    /tmp/s203_induced_growth session194-weak-parent-certificate/inherited_tournaments_through8.txt
    python3 session203-induced51-audit/check.py
    python3 -O session203-induced51-audit/check.py
    python3 session203-induced51-audit/check_repair.py
    python3 -O session203-induced51-audit/check_repair.py
    sha256sum -c SHA256SUMS

The original audit scripts and receipts are retained byte-for-byte, including their sibling directory layout. [DEPENDENCIES.md](DEPENDENCIES.md) separates the cited external order bound from the replayed finite certificate. RELOCATION.json records release-local replays. No exploratory unused counters or unaudited growth probes are part of the evidence.
