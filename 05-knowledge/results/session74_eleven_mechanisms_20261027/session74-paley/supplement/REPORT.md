# A fifteen-complement geodesic from Paley11 to a transitive tournament

**Status:** PROVED exact distance, using the inherited score-Lipschitz theorem and a fully verified explicit witness. FINITE-EXACT certificate and bounded search controls. This is a separate continuation of the frozen six-move forest packet; none of that packet's bytes were changed.

## 1. Result and exact scope

Let P11 be the labelled tournament on 0,...,10 with `x->y` iff `(y-x) mod 11` is in `{1,3,4,5,9}`. One native move reverses every arc OFF the currently named directed Hamilton path. Then

    d_native(P11, transitive class) = 15.

The lower bound is inherited from THM-4810 / S63: `L(P11)=15`, a transitive tournament has L=0, and a legal native complement changes L by at most one. The new contribution is a fifteen-move word attaining this bound.

The word needs no converse, relabelling, or other free operation. Thus it also proves distance exactly fifteen to its particular explicitly labelled transitive endpoint. If full converse, relabelling and legal path remarking are allowed for free, the distance to the transitive class is still fifteen because those operations preserve L. **We do not claim fifteen to every fixed labelled transitive tournament.**

The first six moves are byte-for-byte the previously verified forest geodesic. The remaining nine take that exact literal forest state to the transitive tournament. Hence its own distance to the transitive class is exactly nine. Every contiguous part of the fifteen-move word is a geodesic between its literal endpoints: its length equals the drop in L. The first forest frame occurs at time six, the earliest possible time proved by the parent packet.

No source weights, shared-port identifications, or vertex names are changed; no insertion, substitution, copy, or reconstruction is counted as free. The result concerns this native tournament grammar, not an arithmetic or LRC clock.

## 2. Nine additional named paths

The first six moves, every intermediate row and every rank witness are in the parent `certificate.json`. The self-contained `certificate15.json` repeats those six paths and appends:

| Total move | Target L | Hamilton path of the CURRENT tournament |
|---:|---:|---|
| 7 | 8 | 2,0,4,6,5,9,1,7,3,10,8 |
| 8 | 7 | 3,1,7,4,5,8,6,10,2,0,9 |
| 9 | 6 | 2,0,9,5,3,1,6,7,4,8,10 |
| 10 | 5 | 3,2,0,4,8,5,7,9,10,1,6 |
| 11 | 4 | 2,0,4,5,7,9,6,10,1,8,3 |
| 12 | 3 | 3,1,8,2,0,4,5,7,10,9,6 |
| 13 | 2 | 6,2,9,7,10,3,1,8,0,4,5 |
| 14 | 1 | 5,10,3,1,8,0,2,9,7,4,6 |
| 15 | 0 | 6,0,2,9,7,4,8,5,10,3,1 |

The final forward order is

    6,0,2,9,7,4,8,5,10,3,1.

Every one of its 55 pairs is oriented from earlier to later. Its labelled out-neighbour bitmasks are

    [1982,0,1978,2,1322,1034,1983,1338,1066,1466,10].

The checker validates the actual intermediate orientations, not just score sequences. Scores alone would not certify legal paths or a native word.

## 3. Feedback is not the monotone clock

Along this exact score-geodesic the independently computed feedback values after the fifteen moves are

    14,14,12,12,11,10,10,11,8,7,7,5,4,5,0.

L drops by exactly one at every move, but f rises at moves 8 and 14. In particular f=4 at L=2 before rising to f=5 at L=1; the last native complement then reaches f=0. This is a concrete reason not to replace the inherited score clock by greedy minimization of feedback edges. The parent cycle-rank invariant is likewise not 1-Lipschitz; its first step already drops intrinsic beta from 10 to 5.

## 4. A near miss explains the final guard

The initial purely greedy continuation reached a different L=1 tournament with rows

    [4,2037,512,1015,709,725,517,1093,757,129,893].

Its sorted scores are exactly `[1,1,2,3,4,5,6,7,8,9,9]`, the score profile of a one-move predecessor of a transitive eleven-tournament. Nevertheless it has f=4 rather than the transitive predecessor type's f=5. Independently enumerating ALL 243 legal Hamilton paths of this literal state finds zero L=0 successors. The independent auditor also obtained the full successor histogram `L=1:8, L=2:235` and a separate subset-DP path count of 243.

Thus the first failed implication was 'the right almost-transitive score profile gives a legal one-move landing'. The missing coordinate was orientation and path incidence. The repaired algorithm backtracked at L=3 and verified actual path complements. This hostile is not a new obstruction to the final endpoint: the corrected nine-move continuation succeeds.

There is also a useful safe pruning lemma.

**Lemma (PROVED).** A tournament with a source or sink has no native complement that strictly decreases L.

A source must be the first vertex of every Hamilton path; a sink must be its last. In the converse, their scores are respectively 0 and n-1. Restoring the path therefore moves the corresponding extreme inward by one. Moving a unique extreme score 0 to 1, or n-1 to n-2, increases the absolute score-assignment cost by exactly one. The other endpoint's unit change can reduce that cost by at most one. Thus the total cost cannot decrease, and neither can L. For n=2 the native move is the identity, also satisfying the conclusion. If both extremes occur and n>=3, the two inward changes act at opposite ends and increase L by one.

The continuation therefore discarded source/sink states before the final L=0 landing. The lemma is a necessary gate, not a guarantee of a decreasing successor: the L=1 trap above has no source or sink and still fails. The checker verifies the lemma on every legal move from every source-or-sink tournament of orders 2 through 5.

## 5. Search and independent validation

`extend_greedy.cpp` starts at the exact frozen forest rows, considers legal Hamilton paths in lexicographic order, and takes the first L-decreasing successor with no source or sink while L>0. It reaches L=1, then exhausts all 243 paths there without reaching L=0. Only its first six successful extension moves (target L=8,...,3) are used in the final certificate.

`backtrack.cpp` restarts at that L=3 literal state and uses a bounded depth-first search, with an explicit cap of 5,000,000 Hamilton-path candidate evaluations. It finds the final three-move landing after 17,766 such evaluations and 39 deduplicated states. A constant-size score-tie filter is used at L=1; its successful final move is independently checked, and no completeness claim depends on the filter. This was not a breadth-first search of the eleven-vertex tournament graph.

`check.py` does not import or execute either finder. It reconstructs Paley11 from its quadratic-residue definition using Boolean adjacency matrices, checks each path for legality at its own stage, computes native complements by independent unordered-pair logic, compares every literal row, checks the same-path involution and the whole L sequence `15,14,...,0`, and verifies all 55 final forward pairs. It additionally computes f by an independent exact subset DP and checks the complete trap universe and the small source/sink universe. It records the SHA256 of the parent six-move certificate.

Normal and optimized Python receipts are byte-identical. A separate audit agent independently reconstructed and verified all fifteen moves and the final total order without importing producer code. Its C2 report is supplied separately to the integrating session, with the C1 audit unchanged.

    g++ -std=c++17 -O2 extend_greedy.cpp -o /tmp/session74-extend-greedy
    /tmp/session74-extend-greedy > /tmp/session74-extend-greedy.out
    cmp extend_greedy.out /tmp/session74-extend-greedy.out
    g++ -std=c++17 -O2 backtrack.cpp -o /tmp/session74-backtrack
    /tmp/session74-backtrack > /tmp/session74-backtrack.out
    cmp backtrack.out /tmp/session74-backtrack.out
    python3 check.py > normal.json
    python3 -O check.py > optimized.json
    cmp normal.json optimized.json
    sha256sum -c SHA256SUMS

The lower bound is a theorem, not an exhaustive computation. The upper bound is a witness, not a claim that every greedy branch succeeds. P7 remains the inherited warning: it has L=6 but native distance seven to the transitive class, so attainment by P11 does not imply universal attainment by regular or Paley tournaments.
