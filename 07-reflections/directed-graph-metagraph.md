# The Directed Graph Meta-Graph: A Generalization Hierarchy

*opus-2026-03-23-S250*

## The Four Levels

| Level | Object | State space | #States | Cell states | OEIS |
|-------|--------|-------------|---------|-------------|------|
| 0 | Tournament | {0,1}^{C(n,2)} | 2^m | 2 (binary) | A000568 |
| 1 | Oriented graph | {0,1,2}^{C(n,2)} | 3^m | 3 (ternary) | A001174 |
| 2 | Digraph | {0,1}^{n(n-1)} | 2^{2m} | 2×2 per pair | A000273 |
| 3 | Binary relation | {0,1}^{n^2} | 2^{n^2} | 2 per entry | A000595 |

Each level's meta-graph = Hamming graph / S_n:
- Level 0: Q_m / S_n
- Level 1: H(m, 3) / S_n
- Level 2: Q_{2m} / S_n
- Level 3: Q_{n^2} / S_n

## Why Oriented Graphs Are the Natural Home

**The 3 states per cell in an oriented graph are the 3 orbifold points of PSL(2,Z):**
- State 0 (no arc) = the **cusp** τ = i∞ (degenerate, no structure)
- State 1 (forward) = the **i-point** τ = i (order 2, complement fixed point)
- State 2 (backward) = the **ω-point** τ = ω (order 3, cycle fixed point)

**Tournaments are oriented graphs with the cusp excluded.** Every cell forced to be state 1 or 2. This is the "non-cuspidal" part of the modular curve.

This explains why tournament theory is governed by orders 2 and 3: removing the cusp from PSL(2,Z) leaves only the finite-order generators.

## What Survives and What Breaks

### Survives at all levels:
- Burnside formula for V(n)
- T_n = arc-orbit count via cycle index
- Case A of edge orbits = T_n/2 (the /2 from unordering)
- The S_n action structure

### Breaks at Level 1 (oriented):
- H is no longer always odd (empty graph has H=0)
- The staircase cells become ternary
- The recursive decomposition has 3 apex states
- The complement map only swaps 1↔2, fixing 0

### Breaks at Level 2 (digraph):
- The staircase becomes a SQUARE (not triangle)
- Pairs (i,j) and (j,i) are independent
- Two symmetry operations: complement AND transpose
- The merge group becomes Z_2 × Z_2 (Klein 4-group)

## At n=3

| Quantity | Tournament | Oriented | Digraph |
|----------|-----------|----------|---------|
| Iso classes V | 2 | 7 | 16 |
| Meta-graph edges E | 1 | 11 | 30 |
| Edge orbits | 3 | 15 | ? |
| Gap orbits | 2 | 4 | ? |

The oriented graph meta-graph is 3.5× larger in classes and 11× larger in edges than the tournament meta-graph.

## The Deep Question

Can the edge orbit formula edge_orbits = T_n/2 + (n-2)! generalize to:
- Oriented: edge_orbits = T_n^{orient}/2 + c_1 × (n-2)!
- Digraph: edge_orbits = T_n^{dig}/2 + c_2 × (n-2)!

where c_1, c_2 depend on the cell-state count?

If so: **what are c_1 and c_2?** And do they relate to the generator orders of the symmetry group?

Tournament c = 1 (from PSL(2,Z) with S of order 2).
Oriented c = ? (from the 3-state structure).
Digraph c = ? (from the Klein 4-group).

This is the research frontier.
