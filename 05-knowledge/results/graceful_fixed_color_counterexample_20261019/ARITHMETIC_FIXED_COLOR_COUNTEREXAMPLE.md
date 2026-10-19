# An integer-difference fixed-color fiber disconnected under two-color moves

On labels 0,...,12, leave differences 6,9,11 free. Fix the following nine edges, indexed by their difference:

    1: 5--6       2: 1--3       3: 8--11
    4: 3--7      5: 0--5       7: 4--11
    8: 1--9     10: 0--10     12: 0--12.

These edges form a forest F with four components

    A={0,5,6,10,12}, B={1,3,7,9}, C={2}, D={4,8,11}.

The complete nonloop contact supports of the free differences are

    6: AD, CD
    9: AB, CD
    11: AB, AD.

Their union is the quotient path B--A--D--C, and each color omits a different edge of that path. Every spanning-tree completion must occupy all three path edges. Consequently only two type assignments are possible:

- difference 6 occupies CD, difference 9 occupies AB, difference 11 occupies AD;
- difference 6 occupies AD, difference 9 occupies CD, difference 11 occupies AB.

All three free colors change between these assignments, so no move changing at most two colors connects them while preserving F.

## Actual realizations and exact fiber

The first type assignment has exactly three realizations:

    6:2--8, 9:0--9,  11:0--11;
    6:2--8, 9:1--10, 11:0--11;
    6:2--8, 9:3--12, 11:0--11.

The second has exactly two:

    6:4--10, 9:2--11, 11:1--12;
    6:5--11, 9:2--11, 11:1--12.

Each row, together with F, is a tree with one edge of every difference 1,...,12. The entire completion fiber therefore has five states, with its width-two graph equal to a disjoint triangle and edge. Its connected-component sizes are exactly 3 and 2.

For transparency, the remaining free-color edges are loops after contracting F: difference 6 has 0--6,1--7,3--9,6--12; difference 9 has no loops; difference 11 has no loops. Such loops cannot occur in a completion containing F.

## Consequences and scope

This refutes the assertion that every nonempty arbitrary fixed-color graceful completion fiber is connected under at-most-two-color moves. Upper-tail structure in the earlier theorem cannot simply be dropped.

Choose source and target from different listed components. Their nine fixed edges are already correct target edges. Global width-two connectivity is guaranteed by the upper-tail theorem with empty tail, but every such path must temporarily remove at least one of those already-correct edges. Thus target-preserving reconfiguration fails for actual integer-difference graceful trees as well as for the abstract multistar.

This does not contradict the upper-tail theorem: the frozen set of colors is not an upper interval. It does not imply an exponential diameter lower bound or settle global one-color connectivity. Intermediate tree shapes may change.

## Discovery and direct verification

The witness was found by a targeted four-component-map CSP at m=12, rather than enumeration of all graceful trees. The map assigns labels to the vertices of a four-vertex tree, constrains every free-difference contact to one of its two prescribed cyclic support edges or a loop, requires both supported edges to occur, and searches for a same-component rainbow forest using all other colors. Search results for m=9,10,11 were negative within that exact model; no universal impossibility claim follows from them.

The independent `verify_m12_witness.py` reads only the explicit witness, directly enumerates all choices for the three free colors, checks connectivity, and constructs the exact width-two graph. The counterexample proof does not depend on trusting the CSP or its bounded negative results.

The direct enumeration has exactly 7*4*2=56 candidate free-edge tuples. No claim that m=12 is minimal is made. The frozen witness and direct proof are sufficient regardless of the discovery search.

This supplies a subsequent counterexample to the open refinement recorded in [packet97](../graceful_target_monotonicity_obstruction_20261019/), whose earlier bounded negative controls remain valid. The global two-color connectivity used above is [packet93](../graceful_upper_tail_two_color_connectivity_20261019/). The abstract contact classification is [packet101](../graceful_three_color_contact_classification_20261019/). None is a runtime dependency of this direct counterexample check.
