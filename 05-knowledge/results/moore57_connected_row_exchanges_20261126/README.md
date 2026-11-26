# Connected row exchanges for the fixed six-partition Moore inputs

2026-11-26. Follow-up to the [2,477-vertex incidence fitting packet](../moore57_singer_repunit_six_partitions_20261126/README.md). All changes concern ordinary six-column rows on those same boards. Distinguished rows, prism shell, and allowed pair cells remain fixed. No complete degree-57 Moore graph is claimed.

## New record on the same boards

All substitutions are recorded as literal rows in **moves.txt**. The first stage replayed the twelve supplied trades of size up to three, reaching 2,479. A complete search through size four then found twelve more gains across inputs 1–3. The initial four-row-stable equal-triangles input and the improved paired-squares input both admitted five-row gains; five of those gains are in the final replay. Each recorded trade increases the count by one. All four final witnesses admit **no** positive exchange with at most four removals.

| Input | Predecessor total | Ordinary now | Distinguished retained | Total now | Left per pair board |
|---|---:|---:|---:|---:|---:|
| Equal triangles, concurrent | 2,475 | 2,318 | 163 | 2,481 | 436 |
| Equal triangles, separated | 2,470 | 2,312 | 162 | 2,474 | 443 |
| Paired squares, concurrent | 2,477 | 2,324 | 163 | **2,487** | **430** |
| Paired squares, separated | 2,466 | 2,313 | 162 | 2,475 | 442 |

The largest verified fit is therefore ten above the predecessor. For the best final witness the exhaustive four-removal pass checked 949,725 connected four-owner supports in addition to the smaller ones. Neither the four-row proof nor the bounded five-row improvement run claims maximality under five removals, let alone a full incidence fit.

## Exact reduction from global packing to local exchanges

Every occupied pair cell has one **owner**, the index of the current ordinary row using it. For a candidate row, collect the distinct owners of its 15 pairs. It is admissible after removing an old set D exactly when all its pairs are allowed and its owners are in D. Two new rows can coexist exactly when they agree in at most one coordinate. Thus a positive exchange removing D exists exactly when the candidates blocked by subsets of D contain at least |D|+1 mutually compatible rows. It suffices to use exactly that many. Any old row removed and then reinserted could instead have been left in place, so we exclude old rows as candidates for a minimum-size gain.

For a smallest positive exchange, view the old rows removed as vertices and each new row's blocker set as a hyperedge. A candidate with no blockers would have been insertable immediately. If the blocker hypergraph is disconnected, one component gets more replacements than it removes, contradicting minimality of the original set of removals. A stronger Hall test follows: every nonempty proper subset of deleted owners must touch at least one more replacement row than its own size; otherwise discarding that owner subset and every replacement touching it leaves a smaller gain. In particular every deleted owner is needed by at least two new rows in a minimal gain of size at least two. This is a necessary filter, not a replacement for compatibility checking.

The replay constructs all allowed rows with at most a specified number of blockers (through five). Assign slots left to right, intersecting the currently allowed symbols. Once the blocker limit is full, the new symbol in each pair must either be free or have one of the retained owners. The code then explicitly checks the owners before advancing. Candidate keys index the distinct blocker sets; from each key it grows only connected unions through other keys sharing a blocker. Every support of a minimal gain is among these unions. Within a support it tests compatibility by enumerating a clique in the candidate pool. This is an exhaustive **local** assertion, not global maximality. A neutral trade can change its answer.

Retaining owners, rather than only pair-count totals, is what enables the connectedness reduction. This complements the [tournament port comparison](../moore_route_colors_tournament_ports_20261126/README.md): a path splice must keep the vertices used by a partial path, and this incidence splice must keep the actual owners of conflicting cells. The [114-neighborhood bound](../moore_route_colors_tournament_ports_20261126/TRADE.md) applies to transitions between finished Moore graphs, not to these ordinary-row edits in a partial incidence fit.

## Replay

1. Run **python3 replay.py** from this directory, or pass the predecessor directory as its sole argument. It starts from the predecessor witness for each input, applies the given moves in order, and verifies every allowed cell and absence of repeated pairs after each move. Intermediate outputs are local regenerated files.
2. Compile **g++ -O3 -std=c++17 check.cpp -o check**. For input i, run **./check ../moore57_singer_repunit_six_partitions_20261126/bi.txt afteri.b64 4**, replacing i by 0,1,2,3. Cap five is also supported; if an improving exchange is printed, the search stops at the first witness and does not certify the whole cap on that starting witness.
3. **extend.py** is the optional deterministic search driver for finding more moves. Proof certification does not require rerunning its heuristic choice of improving iterations.
