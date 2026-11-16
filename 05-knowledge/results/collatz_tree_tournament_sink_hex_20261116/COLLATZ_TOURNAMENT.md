# The Collatz tree as a tournament built one node at a time

Session klein-2026-11-16-hex. Tools `scripts/treetour.c` (the tree tournament of the Collatz inverse tree, any build
order, Hamiltonian counts by subset DP to 20 nodes, inflation), `scripts/dfs_check.py` (random trees and orders, brute
force); receipts `treetour_*.out`, `dfs_check.out`.

## 1. The construction (the owner's reading made precise)

Collatz read globally is an ordering: distance from 1 in the inverse tree (height), ties broken by an even/odd
preference and then by value. Read locally it is a comparison of two numbers. The tournament on all positive integers
in which every pair is compared by the global order is transitive (one Hamiltonian path), and the Collatz tree is a
subset of its arcs (each parent is closer to 1 than its children). The content appears when the tree is **built one
node at a time** in the order `ℓ` and the new node is treated as the challenger: it loses to its ancestors, which are
already present, and beats every other node already present. Formally, for a rooted tree `P` (ancestor order `<_P`)
and a linear extension `ℓ` (parents before children), the **tree tournament** `T(P, ℓ)` has

    u → v   iff   u <_P v  (ancestor beats descendant),  or  u, v incomparable and ℓ(u) > ℓ(v).

Adding the next node `v` to `T_(N−1)` gives `T_N` with in-neighbourhood exactly the ancestor chain of `v`
(THM-4876's path-pairs formula applies with the cut (ancestors, rest)). The tournament depends only on `(P, ℓ)`: the
labels are irrelevant, as the owner's relabelling remark asks (replacing nodes by sets of equal cardinality is the
inflation of §4).

## 2. Theorem (tree tournaments; THM-5257)

Let `P` be a rooted tree, `ℓ` a linear extension, `T = T(P, ℓ)`.

1. The 3-cycles of `T` are exactly the triples `{x, y, z}` with `x <_P y`, `z` incomparable to `x` (hence to `y`), and
   `ℓ(x) < ℓ(z) < ℓ(y)`: an ancestor, a descendant, and an interloper built between them.
2. `T` is transitive iff `ℓ` is a preorder (depth-first) traversal: every subtree is built contiguously. Hence
   `H(T) = 1` iff the build is depth-first, and otherwise `H(T) ≥ 1 + 2·c_3(ℓ)` with `c_3` the number of triples in 1.
3. The global Collatz order (breadth-first by distance from 1) is as far from depth-first as a linear extension can be:
   each level is built before the next, so every pair (ancestor at level `i`, descendant at level `j > i + 1`) has as
   interlopers all nodes of levels `i + 1 … j − 1` outside the ancestor's subtree.

*Proof of 1.* Three vertices span 0, 1, 2 or 3 comparable pairs. With 0 the orientation is by `ℓ` alone (transitive).
With 3 they form a chain (transitive). With 2 the comparable pairs share a vertex that is an ancestor of the other two
(two ancestors of one vertex are comparable), and it is a source of the triple. With 1, say `x <_P y` and `z`
incomparable to both: `x → y`, and `y → z → x` iff `ℓ(y) > ℓ(z) > ℓ(x)`. ∎ *Proof of 2.* A tournament is transitive
iff it has no 3-cycle; a node `z` incomparable to `x` is incomparable to every descendant of `x` (an ancestor of a
descendant of `x` is comparable to `x`), so "no triple of 1" says: no node incomparable to `x` is built between `x`
and one of its descendants, i.e. the subtree of `x` occupies a contiguous block of `ℓ` starting at `x`, for every `x` —
the preorder property. `H = 1` iff transitive (a tournament with one Hamiltonian path is transitive), and the bound is
THM-002 with one triangle and singletons. ∎ Checked on 406 random (tree, order) pairs by brute force (`dfs_check.out`).

## 3. The Collatz sequence of counts (FINITE-EXACT)

Build orders: BFS by height with evens before odds (`even-first`), BFS with odds first (`odd-first`), numeric order of
the labels (not a linear extension; the ancestor rule still applies), and depth-first (even child first).

| `N` | even-first node | `H(T_N)` | `c_3` | odd-first node | `H(T_N)` | `c_3` |
|---|---|---|---|---|---|---|
| 1–7 | 1, 2, 4, 8, 16, 32, 5 | 1 | 0 | 1, 2, 4, 8, 16, 5, 32 | 1 | 0 |
| 8 | 10 | 1 | 0 | 10 | 3 | 1 |
| 9 | 64 | 5 | 2 | 64 | 5 | 2 |
| 10 | 20 | 11 | 4 | 3 | 15 | 5 |
| 11 | 128 | 41 | 8 | 21 | 41 | 8 |
| 12 | 3 | 131 | 12 | 20 | 157 | 13 |
| 13 | 21 | 517 | 18 | 128 | 517 | 18 |
| 14 | 6 | 2,107 | 25 | 6 | 2,793 | 28 |
| 15 | 40 | 10,977 | 35 | 40 | 12,359 | 37 |
| 16 | 42 | 58,679 | 47 | 42 | 70,429 | 50 |
| 17 | 256 | 334,877 | 62 | 256 | 332,969 | 62 |
| 18 | 80 | 2,351,809 | 79 | 13 | 2,166,681 | 79 |
| 19 | 84 | 13,778,505 | 99 | 85 | 11,604,367 | 95 |
| 20 | 512 | 87,545,175 | 121 | 12 | 105,295,407 | 121 |

- Even-first: `H = 1` for `N ≤ 8` and the first cycle `5 → 32 → 64 → 5` appears with the node 64 (pre-registered P-D
  **confirmed**, `H(T_9) = 5`, brute-force checked to `N = 9`); odd-first: the first cycle `5 → 10 → 32 → 5` appears one
  node earlier. Both sequences are odd (Rédei) and grow by a factor about 5–7 per node from `N = 12` on; the two
  preferences agree at `N = 9, 11, 13` (5, 41, 517) on different vertex sets and differ elsewhere. At `N = 20` the
  even-first count 87,545,175 and the odd-first count 105,295,407 both exceed 100 (P-D).
- Depth-first build: `H = 1` at every `N` (theorem 2). Numeric order (1, 2, 3, 4, 5, 8, 10, …): `H(T_14) = 79`.
  Orienting incomparable pairs by the labels instead of the build order gives `H(T_14) = 69` (larger beats smaller)
  and 3 (smaller beats larger) on the even-first node set.
- `c_3(N)`: even-first 0 (×8), 2, 4, 8, 12, 18, 25, 35, 47, 62, 79, 99, 121; odd-first 0 (×7), 1, 2, 5, 8, 13, 18, 28,
  37, 50, 62, 79, 95, 121. The bound `H ≥ 1 + 2 c_3` is far from tight (`N = 20`: 243 against 10^8): the count is
  carried by packings with many disjoint cycles.

## 4. Groups of nodes and relabelling by sets

Replacing every node by a transitive block of `c` vertices (the "sets of the same cardinality") multiplies the count
by a factor that is not a power of two: with `c = 2` the even-first counts become 1 (×8), 393, 6,069 at `N = 9, 10`
against 5, 11 (pre-registered P-E **confirmed**: the tree tournament is not a pure cone of THM-5109's kind). Adding a
whole level at once changes nothing in `T_N` (the tournament depends on the order, not on the grouping), but it changes
the *sequence* of intermediate tournaments: the counts above are the one-node-at-a-time sequence, and the level sums
are its subsequence at the level boundaries.

## 5. What this says about Collatz

Nothing about the conjecture: the tree is given, and the counts measure the build order, not the arithmetic. What it
gives is a clean tournament invariant of a linear extension of a tree — the Hamiltonian count of `T(P, ℓ)` — that is 1
exactly for depth-first builds and grows exponentially for the breadth-first build that the Collatz distance order
imposes. The even/odd preference is a second-order effect (same growth rate, different values), and the DFS/BFS
dichotomy is first-order.
