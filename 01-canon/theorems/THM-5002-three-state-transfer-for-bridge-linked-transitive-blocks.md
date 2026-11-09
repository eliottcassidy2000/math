---
id: THM-5002
title: "Three-state packing transfer for bridge-linked transitive blocks and the canonical run-window obstruction"
status: PROVED transfer and family exclusion; FINITE-EXACT complete bounded-value census. Chan-Pak remains OPEN.
session: Session233-composition-20261109
source: 05-knowledge/results/session233_composition_20261109/transducer/RESULT.md
audit: 05-knowledge/results/session233_composition_20261109/transducer_audit/AUDIT.md
---

# Bridge-linked blocks

For a composition (m_1,...,m_k), put a transitive tournament on each consecutive
block. All between-block arcs point backwards except the single base-path
bridge between consecutive blocks. The odd-cycle packing polynomial has an
exact three-state transfer: closed collections and an unfinished terminal
path of either parity. The full polynomial matrix and exhaustive decomposition
are in the [proof](../../05-knowledge/results/session233_composition_20261109/transducer/RESULT.md).

At the Hamilton-count specialization x=2, a block m>=2 has transfer

    [[1,2^(m-1),2^(m-1)],
     [2^(m-2),3^(m-2),3^(m-2)-1],
     [2^(m-2),3^(m-2)+1,3^(m-2)]].

A singleton instead has matrix [[1,2,0],[0,0,1],[1,1,0]]. Start with (1,0,0)
and read the first coordinate after all blocks. The optional extension at
both ends of an unfinished path is essential; omitting the head extension
fails on block sizes (1,3,1), where the true count is 13 rather than 11.

This family contains every complementary endpoint of THM-4951's canonical
run-window construction with all runs at least two. The other endpoint is
an ordinal product of staircase counts. Both families, even closed under
ordinal sums, omit 11. This is a family obstruction, not an obstruction to
general tight tilings or all tournaments.

Every multiblock member of order at least three is strong. Busch's inherited
strong minimum exceeds 10000 at order19. Exhaustive enumeration of every
composition through order18 gives exactly1374 values at most10000. This is a
complete bounded-value census for arbitrary block compositions, not a census
of just the narrower original run-window family. One-block tournaments at
arbitrary orders are transitive and contribute only1.

Independent coefficient checks cover all255 compositions through order8;
Hamilton DP covers all511 through order9; the unpruned census covers262143
compositions through order18. [Audit](../../05-knowledge/results/session233_composition_20261109/transducer_audit/AUDIT.md).
