# Translucent Lines: The Internal Geometry of Tournament Classes

**Session**: opus-2026-03-23-S271

---

## The Two-Level Hierarchy

The hypercube Q_m of all labeled tournaments partitions its edges into two types:

**Opaque lines**: Arc flips that CHANGE the iso class → form the metagraph edges
**Translucent lines**: Arc flips that PRESERVE the iso class → form the internal fiber structure

```
Q_m = translucent subgraph ∪ opaque subgraph   (disjoint edge partition)
```

This creates a two-level hierarchy:
- **Level 1 (coarse)**: G_n/Z_2 with V_merged nodes, connected by opaque lines
- **Level 2 (fine)**: Each node contains n!/|Aut| labeled tournaments, connected by translucent lines

---

## Lines Per Metagraph Edge

Each metagraph edge is formed by multiple opaque lines (labeled tournament pairs). The distribution at n=6:

| Weight (lines) | # edges | % of edges |
|----------------|---------|-----------|
| 480 = 2n!/3 | 6 | 4.2% |
| 720 = n! | 10 | 7.0% |
| **1440 = 2n!** | **110** | **76.9%** |
| 2880 = 4n! | 15 | 10.5% |
| 4320 = 6n! | 2 | 1.4% |

**The dominant weight is 2n! = 1440.** This means 77% of metagraph edges carry exactly 2n! opaque lines. The weight is always a multiple of n!/|Aut|, reflecting the orbit structure.

**Lightest edges** connect classes with large ΔH (dissimilar tournaments — few paths between them). **Heaviest edges** connect classes with small ΔH (similar tournaments — many paths).

---

## Translucent Fraction and H

The translucent fraction (proportion of arc flips preserving the class) **anti-correlates** with H:

| n | Overall translucent % | Transitive (H=1) | Regular (H=max) |
|---|----------------------|------------------|-----------------|
| 4 | 37.5% | 33.3% | 33.3% |
| 5 | 17.2% | 25.0% | 0% |
| 6 | 10.4% | 20.0% | 0% |

The fraction converges to ~2^{2-n} (the universal twin probability from S263).

**Physical interpretation**: Low-H (hierarchical) tournaments are **flexible** — many arcs can be flipped without losing identity. High-H (cyclic) tournaments are **rigid** — every arc is essential to the structure.

---

## The Heap Connection

In Viennot's heap framework:
- Each arc = a **piece**
- Conflicting arcs (sharing a vertex) = pieces that can't commute
- **Translucent arcs = FREE PIECES** — they commute with everything and can be moved without changing the identity

A tournament's **heap** represents its conflict structure:
- **Wide heaps** (many free pieces) = low H, high translucent fraction
- **Tall heaps** (no free pieces) = high H, zero translucent fraction

The transitive tournament has the widest possible heap: n-1 "free levels" where each vertex can be independently scored. The regular tournament has the tallest possible heap: every arc is locked in a cycle with others.

---

## Arc Universality

**Every arc position generates ALL metagraph edges.** At n=5: all 10 arc positions produce exactly 21 edges each. At n=6: all 15 positions produce exactly 143 edges each.

This confirms kind-pasteur's edge-centric theorem (S20dz): E(G_n) can be computed from any single arc position. The metagraph is **isotropic** — no arc position is special.

In heap terms: every piece, when flipped, explores the entire metagraph. The heap structure is "ergodic" with respect to arc positions.

---

## The Translucent Subgraph Within Each Class

The translucent subgraph on the n!/|Aut| labeled tournaments of a class C has:

| Class type | Components | Diameter | Structure |
|-----------|------------|----------|-----------|
| Transitive (H=1) | 1 (connected!) | n(n-1)/2 | Cayley graph of S_n |
| Low H | Few, large | O(n) | Tree-like clusters |
| Mid H | Many, small | O(1) | Small isolated clusters |
| Regular (H=max) | n!/|Aut| (all isolated) | 0 | Zero translucent edges |

The transitive class is **uniquely connected** by translucent lines — you can reach any relabeling from any other through neutral arc flips alone. This makes it the **vacuum state** of tournament space.

---

## Binary Tree Interpretation

The translucent subgraph within a class can be viewed as a **forest of partial binary trees**:

1. **Root**: a canonical representative labeling
2. **Children**: labelings reachable by one neutral arc flip
3. **Branching factor**: number of neutral arcs (0 to ~n/2)
4. **Depth**: translucent diameter

For the transitive class:
- Branching factor = C(n,2) × neutral_fraction ≈ n(n-1)/2 × 2/(n-1) ≈ n
- Depth ≈ C(n,2) (the full staircase size)
- Total nodes = n! (all labeled tournaments)

This is like a **binary search tree** on permutations: each neutral flip is a "comparison swap" that doesn't change the sorted order (because the tournament is already transitive).

For regular classes:
- Branching factor = 0
- Every tiling is an isolated leaf
- The "tree" is a collection of n!/|Aut| disconnected points

This is like a **maximally balanced heap**: no element can be moved without violating the heap property.

---

## Open Questions

1. Is the translucent subgraph of the transitive class literally a Cayley graph of S_n?
2. What is the exact formula for the translucent diameter of each class?
3. Does the edge weight distribution (multiples of n!) have a closed form?
4. Can the 2n! dominant weight be derived from the Burnside factorization?
5. Does the translucent fraction converge to exactly 2^{2-n} or to a different limit?
6. What is the spectral gap of the translucent subgraph within each class?
