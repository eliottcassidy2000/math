# Session140: Collatz geometry and faithful proof transport

**Checkpoint 2 · 2026-11-02.** Independently reviewed written proofs,
finite exact controls, and explicitly cited classical inputs; no new
Collatz or Hadwiger–Nelson solution and no novelty-priority claim.

## The result

The ordinary Collatz graph has a faithful directed **subdivision inside the
honeycomb lattice**, but its routing lengths must be **unbounded**. Every
root prime to three has at least `2^r` distinct predecessors at exactly `5r`
ordinary steps. This exponential growth excludes any bounded-step,
bounded-multiplicity coding into a polynomial-growth lattice. It does not
exclude a topological drawing or a Euclidean placement without uniform
separation. The known `1 -> 4 -> 2 -> 1` cycle becomes a literal hexagon
after two-edge subdivision; that local picture does not supply a global clock.

A complementary theorem realizes every countable locally finite simple
pseudoforest as an **induced Euclidean unit-distance graph**. Thus standard
Collatz admits exact unit-distance placements, despite admitting no bounded
local lattice coding. The realization can be proper: every compact set meets finitely many vertices.
Crossings and arbitrarily close nonadjacent vertices are allowed; the crowding
can escape to infinity rather than accumulating at a finite point.
Its abstract chromatic number is exactly three, independently of convergence;
this is not a lower-bound witness for the plane chromatic number.

The lattice quotient route also has a precise repair: modular edge labels
lift to exact positions precisely when the exact displacement around every
cycle vanishes. The seven-label triangular quotient does not itself give a
Euclidean unit-distance K7. A quotient coloring and a metric realization
preserve different predicates.

## Read next

- [Synthesis and boundaries](SYNTHESIS.md)
- [Concept board and connection contracts](BOARD.md)
- [First-checkpoint scope](CHECKPOINT1.md) and [second-checkpoint scope](CHECKPOINT2.md)
- [Exact first returns, site conjugacy and consumer-relative principles](principles/PRINCIPLES.md)
- [Proper induced unit-distance refinement](hadwiger_nelson_refinement/PROPER_UNIT_DISTANCE_REFINEMENT.md)
- [Hyperbolic binary-tree capacity control](curvature_control/PROOF.md) · [diagram](curvature_control/hyperbolic_binary_tree.svg)
- [Second independent audit](independent/checkpoint2/AUDIT.md)
- [Native shared-tail bridge and h=0 Gram-rank correction](independent/joint_ports/AUDIT.md)
- [Collatz growth and honeycomb proof](collatz_hex/PROOF.md)
- [Pseudoforest realization and cycle-lift proof](hadwiger_nelson/PROOF.md)
- [Seven recovered HN methods and corrections](hadwiger_nelson/HISTORICAL_RECOVERY.md)
- [Verification, source status and replay](VERIFICATION.md)

The honeycomb existence proof uses **CITED** Kühn Theorem 2, including its
disconnected-graph proviso; the input is not reproved or Lean-verified here.
Inherited HN arithmetic/descent material is recovered with individual status
labels. The external HN six-bound remains **CITED / not locally Lean-built**.
General Collatz convergence, the exact plane chromatic number, and LRC(14)
remain open in this work. The first-return/addressing continuation is independently accepted in
checkpoint2. The hyperbolic control concerns an abstract binary tree, not
the full Collatz graph or a regular hexagonal mesh, and is not quasi-isometric.
