---
id: HYP-1796
status: EXPLORATORY
source: opus-2026-05-30-S4
related:
  - THM-356
  - HYP-1786
  - HYP-1787
  - HYP-1790
  - HYP-1791
  - HYP-1792
  - HYP-1793
---

# HYP-1796: Incidence Layer Necessity

## Statement

When a quotient support graph has the expected matching or adjacency shape but
an `F_2` rank statement fails, the theorem is not living in the scalar
invariant or the ordinary adjacency graph.  It is living in an incidence layer:

```text
transfer matrix,
collision hypergraph,
private-pivot structure,
or Smith torsion profile.
```

Endpoint transfer is the model case.  Support matching is necessary evidence,
but the actual proof object is a nonzero minor, private child system,
leaf-peelable collision hypergraph, or odd Smith profile.

## Evidence

- THM-356 proves private odd child witnesses give full row rank, while support
  matching alone does not.
- In unmerged endpoint transfer, every tested parent row has a private odd
  child through `6 -> 7`.
- In merged endpoint transfer, private witnesses stop covering all rows but
  full rank persists.
- In even-graph endpoint transfer, a full support matching can coexist with
  rank failure over `F_2`.
- HYP-1792 shows support-3 SC collision columns are incidence hyperedges, not
  necessarily triangles in the parent merged metagraph.

## Predictions

1. Future failures of naive projection-residue or metagraph-adjacency
   conjectures should often be repaired by recording the operation's incidence
   matrix.
2. Even-graph quotient failures should correlate with 2-primary Smith factors
   or automorphism 2-adics.
3. Merged tournament endpoint transfer should retain odd Smith invariant
   factors longer than the even-graph quotient.
4. Collision hypergraphs that are not clique shadows may still be
   leaf-peelable, giving triangular rank proofs after changing category.

## Test Plan

1. Compute endpoint-transfer Smith profiles one level higher if feasible.
2. Add incidence summaries to future TDA features: row rank, private columns,
   support matching size, zero rows, small Smith factors, and hypergraph
   peel-core size.
3. Re-examine projection-defect examples where scalar features succeed but
   rank or parity fails, and ask whether an incidence matrix was collapsed.

## Sources

- THM-356
- `07-reflections/endpoint-transfer-witnesses-s95.md`
- `07-reflections/endpoint-collisions-as-incidence-not-adjacency-s95.md`
- `07-reflections/residue-phase-incidence-synthesis.md`
