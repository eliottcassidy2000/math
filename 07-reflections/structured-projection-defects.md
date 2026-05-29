# Structured Projection Defects

**Session:** kind-pasteur-2026-05-29-S2
**Computation:** `04-computation/projection_defect_structured_moves_s2.py`
**Results:** `05-knowledge/results/projection_defect_structured_moves_s2.out`

## The Place We Had Not Looked

The S1 projection-defect scan asked a global question: across all waggly lines at a fixed Hamming distance, do the tournament quotient `Q_m -> G_n/Z_2` and the even-graph quotient `Q_m -> E_n` change together?

The answer was mostly yes. The missed question was smaller:

> Which named geometric moves inside the same Hamming layer are responsible for the exceptions?

The S2 scan splits the same tiling hypercube motion into finite-difference probes:

- single tiles, grouped by range
- range flips
- upper and lower strips
- vertex-stars
- full complement-tiling moves

These are not interchangeable samples of a radius. They ask different derivative questions of the two quotient maps.

## Small Examples

At `n=5`, endpoint upper/lower strips and endpoint vertex-stars have tournament-minus-even defect `+0.3125`; their joint-change rate is only `59.38%`. The endpoint of the staircase is acting like a tournament-class amplifier.

At `n=6`, the same endpoint pattern persists: endpoint strips/stars have defect `+0.2109`, and the vertex-star family has defect `+0.1452`, much larger than the single-tile family defect `+0.0121`.

The opposite signal is local. At `n=6`, range-2 single tiles have defects `-0.0664` to `-0.0820`, but still have `91-93%` joint changes. So these are not independent moves; they are strongly coupled moves whose residual defect points toward the even-graph quotient.

This gives a polarity:

- boundary/endpoints/large structured supports -> tournament-only bias
- short local tiles and some middle strips -> even-only bias
- middle Hamming layers -> high synchronization

## Three Unrelated Threads That May Be One Thread

This resonates with three older-looking ideas that had not been forced into the same room.

First, the triangle doctrine says the three sides of the staircase carry different kinds of information: legs, hypotenuse, and interior are not decorative geometry. The endpoint-star defect makes that operational. Boundary motion is a different derivative of the tournament than interior motion.

Second, the SC blowup universal-score theorem says a construction can erase score variation while preserving and doubling other structure. Projection defects may be the local version of that same phenomenon: some moves are seen by score/tournament class and muted by even structure; other moves are seen by cycle structure and muted by tournament quotient symmetry.

Third, the engineering roadmap wants tournament-TDA features. A feature vector made from random Hamming-shell samples would miss the point. Structured probes are closer to tests in a lab: endpoint-star sensitivity, short-range sensitivity, range-flip sensitivity, and complement-tiling sensitivity measure different failure modes of a ranking system.

The deeper pattern may be:

> The project keeps finding structure by comparing two maps that almost agree, then studying the residual as a first-class object.

OCF compares Hamiltonian paths with independent odd-cycle collections. Blue/black mistakes came from conflating line color with class type. Waggly corrections came from conflating one layer with the whole hypercube. Projection defects continue the same lesson, but now make it quantitative.

## Meta-Hypothesis

The hypothesis collection itself has a shape: many good hypotheses are not statements about a single invariant, but about a commutator between two lenses.

Candidate schema:

1. Choose two quotient maps or summaries of the same tiling object.
2. Choose a structured perturbation family, not only a random one.
3. Measure where the two summaries agree, disagree, and silently absorb the move.
4. Interpret the disagreement geometrically: boundary, hypotenuse, interior, cut-space, cycle-space, score, or complement symmetry.

This suggests a search strategy. Instead of asking only "what invariant have we not computed?", ask:

> Which two already-computed invariants have never been differenced along a named geometric move?

## Engineering Reading

For ranking and preference data, the practical object is not just a tournament invariant. It is a perturbation response profile.

Endpoint-star tournament-only bias could detect rankings whose global class is sensitive to adding, removing, or reweighting extreme alternatives. Short-range even-only bias could detect local cyclic redundancy that survives quotienting by tournament isomorphism. High joint synchronization could detect robust structural changes that are visible to both feature channels.

So `tournament_tda.py` should eventually include structured projection-defect features:

- endpoint-star profile
- interior-star profile
- range-2 and range-3 single-tile profile
- range-flip profile
- complement-tiling profile
- random middle-layer baseline

The baseline says how noisy the system is. The structured probes say where the noise lives.

## Next Flashlight

The S2 scan did not yet condition on spine/ribs/sea location in `G_n/Z_2`. That is the next obvious negative space: the same move family may have different defect signs near the principal line, on ribs, or out in the NS sea.

If endpoint-star bias is strongest near the transitive end of the spine, it is probably score/hierarchy. If it persists uniformly in the sea, it is a more intrinsic quotient defect. That distinction matters both for proof search and for engineering use.
