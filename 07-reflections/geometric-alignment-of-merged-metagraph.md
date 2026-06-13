# The Geometric Alignment of G_n/Z_2

**Session:** kind-pasteur-2026-03-23-S20cr
**Status:** FOUNDATIONAL — all future sessions must use this frame

---

## The Picture

G_n/Z_2 is a graph on V_merged = (A000568(n) + SC(n))/2 vertices.
It decomposes as A_C = A_B + A_K where B = blue, K = black, and C = combined.
This is not just a partition of edges — it is a **geometric decomposition with perpendicular axes**.

```
                    H = max (regular, Paley)
                         |
                    NS---SC---NS          <- black ribs, perpendicular
                         |
                    NS---SC---NS
                         |
                    NS---SC---NS---NS     <- more ribs higher up
                         |
                    NS---SC
                         |
                    H = 1 (transitive)
```

**The principal blue line** is the vertical axis: the SC blue spine from transitive (H=1, maximum hierarchy) to regular (H=max, zero hierarchy). All SC classes lie on or near this axis, connected by blue edges.

**The black subgraph** extends horizontally, perpendicular to the blue spine. Every black edge connects an SC node on the spine to an NS node in the halo. The black subgraph is BIPARTITE (SC vs NS) — proved in S20cr.

---

## Three Axes

### Axis 1: The Principal Blue Line (vertical — H direction)

The SC blue backbone, rooted at the transitive tournament:
- **Bottom:** Transitive (H=1, score=(0,1,...,n-1), maximum hierarchy, maximum score variance)
- **Top:** Regular/Paley (H=max, score=(k,k,...,k), zero hierarchy, zero score variance)
- **Structure:** A tree at n=3,4; develops genus at n≥5; fragments at n=8

Along this axis:
- H generally increases (but path to max is non-monotone at n≥6)
- Score variance decreases
- c3 (3-cycle count) increases
- |Aut| fluctuates

The first blue neighbor of transitive has Delta_H = 2^(n-2), matching the hypotenuse formula.

### Axis 2: The Black Perpendicular (horizontal — SC/NS direction)

Black edges extend perpendicular to the blue spine, connecting SC backbone to NS halo:
- Each black edge crosses the SC/NS boundary (bipartite)
- Black edges cluster around SC nodes — each SC node radiates black edges to its NS neighbors
- NS nodes often touch 2 or more SC nodes ("black=2 universality" for most NS classes)

**Bilateral direction:** Near the bottom (transitive), black ribs point UP (NS neighbors have higher H). Near the top (regular), ribs point DOWN (NS neighbors have lower H). The crossover occurs in the middle of the spine.

### Axis 3: The Blue Bulk (approaching completeness)

Beyond the spine, the full blue subgraph fills in NS-NS connections:
- NS nodes connect to each other via blue (same-type) edges
- This NS-NS blue subgraph becomes the dominant component as n grows
- Blue density within NS: avg_degree grows as ~2× per increment of n
- At n=8: 96% of all edges are blue; the blue graph is almost the full graph

---

## Left-Right Imbalance at Odd n

**The most structurally important asymmetry in G_n/Z_2.**

At odd n, the SC backbone has a c3-parity bipartition (c3 can be odd or even, and blue edges tend to change c3 parity). This creates a natural "left" and "right" to the black perpendicular:

| Property | Odd n | Even n |
|----------|-------|--------|
| First SC neighbor of transitive | H=3 | H=5 |
| Bilateral dominance | North > South | South >> North |
| Zero-halo SC classes | 2 (isolated) | 0 (none) |
| SC c3 parity | Mixed (odd+even) | All even |
| Black attachment direction | Asymmetric | More symmetric |

**At n=7:** North(91 nodes) > South(51 nodes) — the transitive side of the spine is larger.
**At n=8:** South(2191 nodes) >> North(282 nodes) — the regular side dominates.

The 2 zero-halo SC classes at odd n are SC nodes with NO black connections — they sit on the spine with no perpendicular ribs at all. These are the most "purely blue" nodes.

---

## Tiling Cell Uniformity

**Every cell in the staircase delta_{n-2} generates the same blue fraction of meta-graph edges.**

```
n=3: 1.000 (uniform across all 1 cells)
n=4: 0.200 (uniform across all 3 cells)
n=5: 0.547 (uniform across all 6 cells)
n=6: 0.691 (uniform across all 10 cells)
```

This means:
- The blue/black decomposition is ISOTROPIC on the staircase
- No cell is "more blue" or "more black"
- The (r,c) ↔ (c,r) reflection symmetry is PERFECT
- Row, column, and strip averages are all identical

BUT: the SC source blue fraction and NS source blue fraction are DIFFERENT and INVERT between odd and even n:
```
n=5: SC sources → 0.68 blue, NS sources → 0.20 blue
n=6: SC sources → 0.25 blue, NS sources → 0.80 blue
```

---

## The Approach to Tournament-ness

As n → ∞:
- E_blue/E_total → 1 (blue dominates)
- E_blue/C(V,2) → 0 (still sparse vs complete)
- Avg blue degree doubles each step: 1.0, 0.7, 2.6, 5.8, 11.6, 24.7
- The H-oriented blue subgraph is a TRANSITIVE partial tournament (0 directed 3-cycles)

The blue subgraph approaches a complete graph in the sense that it becomes the ONLY graph: black edges become negligible. The meta-graph IS the blue graph, and the blue graph is a sparse approximation to a complete graph (tournament underlying graph) on V_merged vertices.

---

## Summary for Future Agents

**Always think of G_n/Z_2 as having three perpendicular structures:**

1. **BLUE SPINE** (principal line): SC classes connected by blue, from transitive to regular. The axis of the whole graph. Carries H-gradient.

2. **BLACK RIBS** (perpendicular): Bipartite SC↔NS connections extending horizontally from the spine. Left-right imbalanced at odd n.

3. **BLUE BULK** (approaching complete): NS-NS blue connections filling the space around the spine. This is the "flesh" on the "skeleton."

When computing any invariant, always track it along all three structures. The spine tells you about the H-landscape. The ribs tell you about SC/NS mixing. The bulk tells you about the approach to completeness.
