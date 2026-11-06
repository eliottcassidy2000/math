# Audited finite five-level rebuilding certificate

Status: **FINITE-EXACT, independently audited.** This is not an all-order proof, not a new finite scalar coverage range, and not a proof of an eventual threshold.

Every seed count23,25,...,45 has an actual two-apex tournament representative with a complete depth-five strategy. At every internal node, both children2H-1 and2H+1 exist in the displayed rebuilding menu and are recursively ready for the remaining levels. The choice may depend on the actual current word and remaining depth.

The certificate contains744 rewrite edges and755 distinct literal tournaments, all of order at most12. Its vertices cover every odd count23..1471. That interval is far below the project's previously certified scalar range; the result tests a physical rebuilding mechanism rather than extending numerical coverage.

A word w encodes two fixed apices0->1 and a transitive chain. Each letter0..3 says which apices that chain vertex beats. The rewrites globally change incidence-letter types, optionally erase a letter class, and then add at most one letter at each end or insert one letter at a specified chain position. These operations construct actual tournaments but can change old edges and delete vertices. They do not retain the original tournament as an induced subgraph.

## Portable verification

All required inputs are included. The two C++ programs require only a C++17 compiler and its standard library.

    g++ -O3 -std=c++17 verify_strategy.cpp -o verify_strategy
    ./verify_strategy STRATEGY5.tsv STRATEGY5_SEEDS.txt
    g++ -O3 -std=c++17 audit/independent.cpp -o independent
    ./independent STRATEGY5.tsv

The first program rebuilds every rewrite literally, counts every tournament by Hamilton-path subset DP, and verifies both-sign recursive readiness. The independently written second program uses a reverse-path subset recurrence to recount every listed graph. The release verifier explicitly caps inputs at the certified order12, ensuring signed64 count safety because12! is below its limit. This narrows the original author's generic guard in response to the independent audit's portability caution; it changes no certified graph or count.

`STRATEGY5.tsv` has eight whitespace-separated fields: parent word, parent H, child word, child H, operation kind, four map images, parameter1, parameter2. Map digit4 means erasure. For `ends`, parameters are prefix/suffix letters with-1 meaning empty. For `interior`, they are the insertion position in the mapped word and the inserted letter. `STRATEGY5_SEEDS.txt` identifies the twelve roots.

## Remaining proof obstruction

There is no proved infinite invariant family, closed guard system, or rule guaranteeing indefinite rebuilding. Finite depth does not imply arbitrary depth, and changing to another representative of the same count requires a justified physical construction.

`FINITE_MENU_BRIDGE.md` is a separately reviewed **unsatisfied sufficient criterion**, not a successful compiler. It explains how regular guarded rewrites, or a fixed finite menu checked by a tensor-product zero identity, could establish an all-word doubling law. No menu passes that test here. In particular, freely choosing an interior position in a finite certificate does not instantiate a fixed finite all-order rule.

`audit/REVIEW.md` records the independent mathematical and computational review. No eventual odd-count threshold follows from this packet.
