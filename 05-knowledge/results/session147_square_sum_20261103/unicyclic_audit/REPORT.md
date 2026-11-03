# Independent audit: unicyclic paths, growth, and two edge families

**Verdict:** Both [unicyclic theorem](../unicyclic/THEOREM.md) and [growth proof](../growth/REPORT.md) are correct in the statements audited. No mathematical errors found. The independent exhaustive computations and explicit failure controls described below pass identically with ordinary Python and `python -O`.

## Unicyclic theorem

In a connected simple unicyclic graph, each non-cycle edge is a bridge. Any connected spanning subgraph contains all bridges; an n-vertex Hamilton path has n−1 edges versus the graph's n, so the unique omitted edge is on the cycle. Conversely, deleting a cycle edge leaves a connected tree, and a finite connected tree with maximum degree at most 2 is a path. Its only degree changes are reductions by one at the deleted edge's endpoints. Consequently any off-cycle vertex of degree at least 3, or any cycle vertex of degree at least 4, obstructs all paths; otherwise the deleted edge must cover every degree-3 cycle vertex. A cycle has at least 3 vertices, so the number of candidate edges is the cycle length when S is empty, 2 when |S|=1, 1 when |S|=2 and the two vertices are adjacent, and 0 otherwise. This remains correct for triangles: every pair of distinct cycle vertices is adjacent, but three degree-3 vertices cannot be covered by one edge. Off-cycle branching and multiple branches at one cycle vertex are included in the obstruction. The edge-forcing and edge-deletion component-count claims in THEOREM.md also follow as stated; no vacuous path claim is made when there are no paths.

## Independent finite audit

`audit.py` generates all labelled trees by Prüfer codes, adds each missing edge, and deduplicates. This is exhaustive because deleting any cycle edge from a connected unicyclic graph yields a labelled tree. It identifies cycle edges independently by testing connectivity after one-edge deletion, then computes all Hamilton path edge sets by visited-subset/endpoint dynamic programming, not vertex permutations or the original leaf-peeling method. It compares those path sets with the predicted deletions and tests every subset of deleted graph edges against the component identity. It also checks six malformed/disconnected inputs against the original validator, rejects a non-unicyclic test graph, and exercises an intentionally false condition using explicit exceptions rather than `assert`. Results:

| n | graphs | distribution: number of paths → number of graphs |
|---|---:|---|
| 3 | 1 | 3 → 1 |
| 4 | 15 | 2 → 12, 4 → 3 |
| 5 | 222 | 0 → 30, 1 → 60, 2 → 120, 5 → 12 |
| 6 | 3660 | 0 → 1440, 1 → 1080, 2 → 1080, 6 → 60 |

The total is 3,898 graphs. Cycle lengths 3–6 all occur; 1,140 order-six graphs exhibit a degree obstruction, and other zero-path cases arise from too many or nonadjacent degree-3 cycle vertices. `audit.out` and `audit_O.out` match byte for byte. Both original checkers also ran in ordinary and optimized modes with matching outputs.

## Square-sum growth

The older neighbors of n correspond exactly to integer k with n < k² ≤ 2n−1. Hence d(n)=⌊√(2n−1)⌋−⌊√n⌋. For n=5–8, 9 lies in the interval. For n≥9, r=⌊√n⌋≥3 and (r+1)²≤2r²−1≤2n−1 since r²−2r−2≥0. Thus d(n)≥1 for all n≥5.

The independent direct pair enumeration reports G1 connected with rank 0, G2–G13 disconnected, and G14 connected with 13 edges, hence a tree. The older-neighbor sets for n=15–19 are {1,10}, {9}, {8}, {7}, {6,17}, respectively. The first four therefore have connected cycle rank 1 and G19 has rank 2. Since adding n≥15 preserves connectedness and changes rank by d(n)−1≥0, exactly G15–G18 among all connected prefix graphs are unicyclic. The full G1–G19 table appears in `audit.out`. The infinite conclusion relies on the inequality and induction, not the finite tests.

## Two infinite square-sum edge families

For positive k, k and k+1 are adjacent exactly when 2k+1 is an odd square (2t+1)². Positivity forces t≥1, and solving gives k=2t(t+1), t≥1. Examples are 4–5, 12–13, 24–25.

Likewise, k and k+2 are adjacent exactly when 2k+2 is an even square (2t)², equivalently k=2t²−1 for t≥1. Examples are 1–3, 7–9, 17–19. In particular, an even k makes 2k+2 congruent to 2 modulo 4, so every such gap-two pair is a nonedge, including 2–4.

Neither type of edge can belong to a triangle, in any prefix or the infinite graph: if a third vertex w were adjacent to both endpoints, the two positive squares k+w and k+w+1 (or k+w+2) would differ by 1 (or 2). Distinct positive squares differ by at least 3. `edge_families.py` checks both membership formulas for k=1,…,100000 and searches G1000 for common neighbors of the relevant pairs; its normal and optimized outputs agree. These finite checks support the exact proof but do not establish it on their own.
