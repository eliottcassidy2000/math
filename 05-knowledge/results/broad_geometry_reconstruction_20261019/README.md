# Broad geometry and reconstruction connections

2026-10-19. This packet contains six precise themes recovered from prior
repository work. Unrestricted algebraic results are separated from the
conditional geometry of the supplied S6 manuscript.

## Read first

- [Geometry overview](geometry/geometry_extension_overview.md): Seifert gauge
  torsion, exterior monodromy, and conductor deformation holonomy.
- [Three independently selected bridges](bridges/unprompted_bridges_session7.md):
  nonlinear conductor data, critical-group repair, and bounded-run moments.
- [Independent review](review/independent_review.md): proofs, separate exact
  controls, and all retained assumptions.

## Results and exact boundaries

1. The marked two-fibre Seifert gauge quotient is Z plus Z/gcd(a,b), with an
   exact primewise count of admissible classes. A cyclic clutch group can
   forget this finite marking information. This is not an analytic moduli
   classification.
2. Exterior degree two of a rank-two square-zero monodromy retains only the
   first Smith factor. Identical middle-degree data can coexist with an
   arbitrarily large second factor. Component/Euler consequences require
   the stated periodic toric quotient hypotheses.
3. The theta conductor has two degree-zero line-bundle holonomies. The
   manuscript's explicit t-preserving toric maps give both equal to one,
   conditional on the locally unramified quotient. This verifies one
   necessary smoothing compatibility, not a complex structure on S6.
4. The same linear tangent matroid can hide different higher conductor
   defects. Explicit rational line configurations give an unbounded loss;
   six branches are the first uniform-rank-three boundary, classified by
   whether the projective directions lie on a smooth conic.
5. The printed graph core needs exactly four edge duplications for its cut
   to become a critical-group character; preserving its involution needs
   six. These alter the graph's integral lattice and are not asserted to
   be lawful dynamical transitions.
6. The first 2k power moments reconstruct any integer set with at most k
   runs, by an exact endpoint ratio. The moment budget is sharp continuously
   and at discrete k=1,2; no all-k discrete sharpness claim is made.

Classical antecedents are cited in the individual notes. These are scoped
repository extensions; historical priority has not been certified.

## Reproduce

Python3 and SymPy are required. No network access is needed for these checks.

    cd geometry
    python3 verify_geometry_extensions.py
    python3 -O verify_geometry_extensions.py
    cd ../bridges
    python3 verify_unprompted_bridges.py
    python3 -O verify_unprompted_bridges.py
    cd ../review
    python3 verify_geometry_independent.py
    python3 -O verify_geometry_independent.py
    python3 verify_bridges_independent.py
    python3 -O verify_bridges_independent.py

Normal and optimized outputs agree. The publication copy of the geometry
checker replaces its original assert gates with explicit require calls,
so optimized mode exercises the same tests. The review retains this repair's
provenance. Original-source pins and final publication pins are distinct.
