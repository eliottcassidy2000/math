---
id: THM-4535
title: "The Cartesian product of the Petersen graph and a triangle is not the graph of any polytope (computer-assisted: an exact SAT search over all 7681 induced cycles of P x C3 as candidate 2-faces, with polyhedral vertex figures, edge-figure orientation, and lazily added exact cuts for facets and Z/2 homology, is UNSAT; dimension 5 is excluded by Pfeifle-Pilaud-Santos Thm 2.3)"
status: >
  PROVED (computer-assisted, exact) + INDEPENDENTLY REPRODUCED (audit 2026-10-01: own re-implementation, UNSAT after
  225 iterations; Glucose 4 and MapleChrono agree on its final formula; also UNSAT without the F6 escape-literal
  argument). No DRAT/LRAT proof certificate: the cut layer is custom code, and its validity is argued in the note.
  Statement: P x C3 (30 vertices, 5-regular) is not polytopal.
  - dim <= 3: non-planar.
  - dim >= 6: impossible, the graph is 5-regular.
  - dim 5: a 5-polytope with a 5-regular graph is simple; PPS Thm 2.3 (arXiv:1009.1499v1 numbering: a product is
    simply polytopal iff its factors are) and the non-polytopality of P exclude it.
  - dim 4: no 4-polytope boundary complex has this graph. The SAT model encodes only necessary conditions.
    - Candidate 2-faces: the complete list of 7681 induced cycles (lengths 3..12).
    - The recorded run (repository code, after the audit's F6 fix) is UNSAT after 184 lazy iterations (10 minutes).
      Its final formula (1,394,307 clauses) is re-solved from scratch by Glucose 4: UNSAT.
    - Earlier runs of the code before the F6 fix also ended UNSAT. They are reruns, not independent checks; the
      independent check is the audit's re-implementation.
  PPS (end of their Section 2.4) cite P x C3 as the graph of a cellular S^1 x RP^2, an example where local
  combinatorics cannot decide. Here both global layers are needed (audit-checked):
  - static constraints + homology are satisfiable without the facet checks;
  - static constraints + facets are satisfiable without homology.
  Validation:
  - known 4-polytopes are accepted;
  - K33 x K2 and P x K2 are rejected;
  - without the homology condition the pipeline recovers PPS's cellular S^1 x RP^2 with graph K33 x K3 (6 triangular
    prisms + 9 cubes) and rejects it with that condition. This is consistent with PPS's announced, unpublished claim
    that this complex is the unique strongly regular combinatorial manifold with that graph.
source: opus-2026-10-01-S15 (collatz-functional-uniqueness-20261001), the owner's prompt "try proving the product of two Petersen graphs is not 4-polytopal"
depends_on:
  - Pfeifle, Pilaud, Santos, Polytopality and Cartesian products of graphs, arXiv:1009.1499 (Israel J. Math. 192 (2012)), Thm 2.3
  - Tutte (faces of 3-connected planar graphs = induced non-separating cycles); Steinitz
note: 05-knowledge/results/petersen_product_4polytope_attempt_20261001.md (sections 5, 6.2, 9)
scripts: 04-computation/experiments/petersen_product_4polytope_20261001.py (subcommands `pc3 31 --crosscheck`, or `all-long --crosscheck`)
output: 04-computation/experiments/petersen_product_4polytope_20261001.out
---

# THM-4535 — Petersen × triangle is not polytopal

**Status: PROVED (computer-assisted, exact SAT over the complete candidate list); INDEPENDENTLY REPRODUCED by
the audit's re-implementation; no proof certificate.**

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

   The SAT encoding of (I), (V) and (R), with exact lazy cuts for (H) and (F), is UNSAT. The soundness of each cut,
   including the closure argument for whole-class cuts, is in the note, §5.

## Reproduction

`python 04-computation/experiments/petersen_product_4polytope_20261001.py pc3 31 --crosscheck` runs the lazy
CaDiCaL 1.9.5 loop and then the Glucose 4 re-solve of the final formula, via `python-sat`. The audit's independent
re-implementation is `run_audit.py pc3 0 path --xcheck` in the session's `audit22/` folder.
