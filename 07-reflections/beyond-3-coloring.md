> **CORRECTED HISTORICAL SYNTHESIS (2026-08-25).**  The former version said
> cubic graphs are always two-edge-colorable and identified three-edge-
> coloring with a `Z/3Z` flow.  Both claims are false; see MISTAKE-507.

# Cubic Edge Coloring, Flows, and Matchings

## Edge-coloring spectrum

Every nonempty cubic graph needs at least three edge colors because three
edges meet at each vertex.  For a simple cubic graph, Vizing's theorem gives
an upper bound of four, so its chromatic index is either three or four.
Snarks occupy the class-four side (with additional bridgelessness,
connectivity, or girth hypotheses depending on convention).  Thus two colors
never suffice, three colors are the exact obstruction threshold, and four
colors suffice for simple cubic graphs.

## Flow spectrum

- A cubic graph has a proper three-edge-coloring exactly when it has a
  nowhere-zero `F_2^2` flow, equivalently a nowhere-zero 4-flow.
- A nowhere-zero `Z/3Z` flow is a different condition.  In a cubic graph its
  edge orientations normalize to unit values, forcing every vertex to be a
  source or a sink and hence forcing bipartiteness.
- Tutte's nowhere-zero 5-flow conjecture for bridgeless graphs remains OPEN.
- Seymour's nowhere-zero 6-flow theorem gives a nowhere-zero 6-flow for every
  bridgeless graph.

The integers `3,4,5,6` label distinct algebraic constraints; they are not an
Eisenstein inert/split/ramified classification.

## Matching view

Petersen's matching theorem gives a perfect matching in every bridgeless
cubic graph.  Three disjoint perfect matchings partitioning the edges are
equivalent to a proper three-edge-coloring.  The Berge--Fulkerson conjecture,
asking for six perfect matchings covering each edge twice, remains OPEN.

These statements suggest a useful boundary question: for a small edge cut,
record not only whether a coloring or matching extends, but the complete
extension relation, its multiplicities, and the action of boundary
automorphisms.  THM-4116 implements the coloring-count version; resistance,
oddness, and critical deletion data remain needed to distinguish different
global snark obstructions.
