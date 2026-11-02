# Independent count certificates for the two recovered prime witnesses

**Scope:** exact validation of the displayed graphs, not a new grammar, minimum-order claim, or proof of spectrum completeness.

The independent checker reconstructs every edge from the block recipe and checks that its Boolean adjacency agrees with both frozen representations. It verifies strong connectivity and primality separately. The initial four-vertex block is decoded as outgoing rows `(8,1,3,6)`; all later blocks are singletons ordered after the old parent. No arbitrary numerical state is promoted to a graph.

## A direct native-coordinate certificate

For a parent tournament A and apex x with cut S, let `(h,F,q,g)` denote its Hamiltonian count, the full A+x count, the number of ordered path covers `(P,Q)` with P possibly empty and nonempty Q starting in S, and the number of parent Hamilton paths ending in S.

Append a singleton z after every old parent vertex.

- If `z -> x`, then `(h,F,q,g)` becomes `(h,F+q+g,2q,0)`.
- If `x -> z`, it becomes `(h,F,2q+h,h)`.

For the first law, full paths ending at z are old full paths not ending at x, giving `F-h+g`. Every other full path uses `z -> x`: a nonempty suffix after x contributes q, while an empty suffix contributes h. The total is `F+q+g`. For q, z must end one of the two paths, giving two choices and `2q`. For the second law, z is a universal sink in the full graph, so F is unchanged. The same two choices contribute `2q`; the extra cover with second path `(z)` contributes h. Endpoint counts are immediate.

The checker obtains the initial states by literal permutations and subset covers, then uses only these exact graph laws:

| Target | Order | Literal initial `(h,F,q,g)` | Final `(h,F,q,g)` |
|---:|---:|---|---|
| 80557 | 19 | `(5,11,7,1)` | `(5,80557,120658,0)` |
| 100000007 | 29 | `(1,1,1,1)` | `(1,100000007,162678086,0)` |

All 2,196 singleton transitions from every labelled parent of orders 1–4 and every cut also agree with literal counting. Reversing one apex edge gives genuine hostile tournament counts 135233 and 17494899, respectively, rather than the requested targets.

## A second, direct count of the entire 29-vertex graph

Deleting its apex leaves the transitively ordered vertices `0,...,m-1`, with `m=28`. Write `b_i=1` exactly when the apex beats vertex i. Every Hamilton path consists of an increasing set P, the apex, and its increasing complement Q.

Empty P contributes `b_0`; empty Q contributes `1-b_(m-1)`. Otherwise let `i=max(P)` and `j=min(Q)`. The apex incidences require `b_i=0,b_j=1`. If `i<j`, necessarily `j=i+1`, with one partition. If `j<i`, the vertices strictly between j and i may independently lie in either part, giving `2^(i-j-1)` partitions. Therefore the exact full-graph count is

```
b_0 + 1-b_(m-1)
+ sum_(i=0)^(m-2) (1-b_i)b_(i+1)
+ sum_(0<=j<i<m) (1-b_i)b_j 2^(i-j-1).
```

The checker verifies the transitive-parent condition from the actual adjacency and evaluates this formula as **100000007**. This is not a recurrence-only certificate and is not represented as subset DP.

For the 19-vertex graph, the auditor separately rebuilt and ran the producer's checked-uint64 full endpoint subset-DP program; it returned **80557**. Its order limit excludes the 29-vertex graph deliberately. The independent Python checks pass identically in normal and optimized modes.

The inherited large-prefix computation was not rerun. Its source identity is checked separately; these two witness validations make no all-value coverage inference.
