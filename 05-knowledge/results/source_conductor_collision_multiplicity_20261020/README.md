# Source conductor follow-through

The proof packet [ACTUAL_SOURCE_CONDUCTOR_AND_COLLISION_MULTIPLICITY.md](ACTUAL_SOURCE_CONDUCTOR_AND_COLLISION_MULTIPLICITY.md) advances the prior planar-Jacobian source-slice and collision-locus target. It identifies the classical D-resultant with the actual normalization conductor and collision Fitting determinant, proves a global norm-square and exact slice collision length 2 delta, and classifies the minimum conductor-degree-two sector as a nodal-cylinder equalizer.

Status: scoped proofs and exact controls, [independently reviewed and approved](review/INDEPENDENT_REVIEW.md). The proof retains its original dated review-pending label; the separate completed audit records final approval and frozen source hashes. Classical complex JC remains open. The D-resultant machinery is explicitly credited to the classical literature. No generic-degree-four target is proposed.

Run `python3 verify_conductor.py` and `python3 -O verify_conductor.py`: 293 exact checks, 30 positive/hostile source families, identical outputs. All positive Keller examples are automorphisms; all collision-bearing examples are explicitly non-Keller. `normal.json` and `results.json` are the two identical replay results.
