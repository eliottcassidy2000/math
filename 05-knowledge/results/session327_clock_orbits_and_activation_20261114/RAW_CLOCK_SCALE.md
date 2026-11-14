# A second 30,000 scale, when the number is a source clock

The number near 30,000 has two different meanings in this work. As an output N in the slope-99 chart it belongs to the six-vertex strong-component grade. As a *raw Q2 budget*, a different exact recursive family already reaches just below 30,000, at order eleven. These are not the same recurrence.

## An adjacent-reversal chain weighted by slots

Label n vertices 0,…,n−1. Orient adjacent pairs i→i+1, and every nonadjacent pair from larger to smaller. This is the classical diameter-(n−1) tournament described in [Savchenko, Theorems II–III](https://arxiv.org/abs/2311.13377). It is strong for n≥3: the adjacent forward chain closes by n−1→0. Every directed cycle is one forward interval of length at least three, closed by the unique return edge. Indeed for any cut between consecutive labels, only the adjacent edge crosses forward. A cycle crossing all cuts between its extrema must use every one of those consecutive edges. Conversely any such interval closes.

In the [Grinberg–Stanley positive expansion](https://arxiv.org/abs/2307.05569), a collection of vertex-disjoint odd directed cycles contributes 2^k z^b, where k counts nontrivial cycles and b counts all blocks, including singleton vertices. For this chain the packing exactly tiles the line by singletons of weight z and odd intervals of length ≥3 and weight 2z. Therefore

    Σ_(n≥0) Q_n(z) t^n = (1−t²)/(1−z t−t²−z t³).

Thus the same graph recursion gives tribonacci at z=1 and a different clock recursion at z=2. Its profiles at orders 10,11,12 are

| order | H | Q2 |
|---:|---:|---:|
| 10 | 193 | 11200 |
| 11 | 355 | 29780 |
| 12 | 653 | 79184 |

All entries are recovered both from the tiling recurrence and independent mask subset DP. The order-eleven budget is only 250 below 30,030, but no equivalence with the primorial is inferred.

## Exact maximum grade throughout a full budget interval

By the new [optional-slot merge floor](OPTIONAL_SLOT_FLOOR.md), any strong tournament of order at least twelve has Q2≥32·6^4=41472. Ordinal joins multiply Q2, so the same lower bound applies to any tournament that contains such a strong component. Consequently, for every raw budget X in [29780,41472), the largest possible strong-component order among tournaments with Q2≤X is exactly eleven. In particular this holds at X=30000. This is a budget statement; it does not claim the chain has Q2 exactly 30000.

The recurring scale comes from how a finite grammar consumes a budget. The apparent same number near 30,000 can mean an order-eleven raw clock, a six-component-grade 99-scaled output, or the primorial's modular period. Distinguishing them produces more testable predictions than tying all three meanings to the same factorization.
