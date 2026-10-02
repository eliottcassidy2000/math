---
id: THM-4535
title: "The Cartesian product of the Petersen graph and a triangle is not the graph of any polytope (computer-assisted: an exact SAT search over all 7681 induced cycles of P x C3 as candidate 2-faces, with polyhedral vertex figures, edge-figure orientation, and lazily added exact cuts for facets and Z/2 homology, is UNSAT; dimension 5 is excluded by Pfeifle-Pilaud-Santos Thm 2.3)"
status: >
  PROVED (computer-assisted, exact). INDEPENDENT AUDIT: PENDING.
  Statement: P x C3 (30 vertices, 5-regular) is not polytopal.
  - dim <= 3: non-planar;
  - dim >= 6: the graph is 5-regular;
  - dim 5: a 5-polytope with a 5-regular graph is simple; PPS Thm 2.3 (a product is simply polytopal iff its factors
    are) and the non-polytopality of P exclude it;
  - dim 4: no 4-polytope boundary complex has this graph. The SAT model encodes only necessary conditions, over the
    complete list of 7681 induced cycles (lengths 3..12) as candidate 2-faces, and is UNSAT after 263 lazy iterations.
  PPS (2012, end of Section 2.4) cite P x C3 as the graph of a cellular S^1 x RP^2. It is an example where local
  combinatorics cannot decide; the decision here comes from facets and homology.
  Validation:
  - known 4-polytopes are accepted;
  - K33 x K2 and P x K2 are rejected;
  - without the homology condition, the pipeline recovers PPS's announced unique manifold with graph K33 x K3
    (6 triangular prisms + 9 cubes), and rejects it with that condition.
source: opus-2026-10-01-S15 (collatz-functional-uniqueness-20261001), the owner's prompt "try proving the product of two Petersen graphs is not 4-polytopal"
depends_on:
  - Pfeifle, Pilaud, Santos, Polytopality and Cartesian products of graphs, Israel J. Math. 192 (2012) (Thm 2.3)
  - Tutte (faces of 3-connected planar graphs = induced non-separating cycles); Steinitz
note: 05-knowledge/results/petersen_product_4polytope_attempt_20261001.md (section 6.2)
scripts: 04-computation/experiments/petersen_product_4polytope_20261001.py (subcommand `pc3 31`, also in `all-long`)
output: 04-computation/experiments/petersen_product_4polytope_20261001.out
---

# THM-4535 — Petersen × triangle is not polytopal

**Status: PROVED (computer-assisted, exact SAT over the complete candidate list); independent audit PENDING.**

## Statement

The Cartesian product `P×C₃` of the Petersen graph with a triangle is not the graph of a polytope of any dimension.

## Proof outline

1. **Dimensions ≠ 4.**
   - The graph is non-planar, which excludes dimension 3.
   - It is 5-regular, which excludes dimension 6 and up.
   - In dimension 5 the realization would be simple. By PPS Thm 2.3 the factors would then be polytopal, but `P` is not.
2. **Dimension 4.** Suppose `Π` were a 4-polytope with this graph. Its 2-faces are induced cycles of the graph, and all
   7681 of them are listed (lengths 3 to 12). The face lattice of `Π` satisfies:
   - **(I)** two 2-faces never share a non-adjacent pair of vertices;
   - **(V)** every vertex figure graph is polyhedral;
   - **(R)** the 2-faces around an edge are cyclically ordered the same way from both ends;
   - **(H)** the 2-faces span the `Z/2` cycle space;
   - **(F)** the facets, read off from the faces of the vertex figures and glued across edges, are 3-polytopes that
     meet in common faces.

   The SAT encoding of (I), (V) and (R), with exact lazy cuts for (H) and (F), is UNSAT. The soundness of each cut is
   argued in the note, §5.

## Reproduction

`python 04-computation/experiments/petersen_product_4polytope_20261001.py pc3 31` takes 20–40 minutes with CaDiCaL
1.9.5 via `python-sat`.
