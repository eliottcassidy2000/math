# Geometry research checkpoint 2

19 October 2026. Mostly-text checkpoint for the planar Jacobian research lane. No change to the main checkout or Git was made.

## Results

1. **Autonomous polynomiality beyond one fixed line.** A point-adically convergent scalar-time flow of a polynomial H in K+(x,y)^3 is polynomial at one nonzero time exactly when its derivation is locally nilpotent, equivalently H is a polynomial in a polynomial coordinate. The proof explicitly uses Shpilrain--Yu, Friedland--Milnor/Bedford--Smillie and Rentschler, with ground-field descent separated. An isolated critical point rules out every nonzero polynomial time. Common-invariant polynomial maps and finite words are included.
2. **Finite-dimensional normalized Hamiltonian Lie algebras are abelian.** A self-contained formal-power-series argument shows that finite-dimensional, or finite-step-nilpotent, Hamiltonian Lie structure in order at least three collapses to a common commuting family. Quadratic and infinite-dimensional metabelian examples show the boundaries.
3. **Fixed-direction jet completion and global boundary control.** Every symplectic N-jet tangent to identity has a completion by N+2 fixed-direction polynomial shears of degree at most N, with total degree at most N^(N+2). Existence is classical and explicitly attributed. A common finite-dimensional bound on the inverse coordinates, equivalently a suitable common inverse pole divisor at infinity, would make the formal limit globally invertible; no such bound is claimed in general.

The autonomous classification supplies a genuine global hypothesis beyond a selected source line. It does not prove arbitrary Keller maps admit a polynomial autonomous logarithm. The explicit Hénon-type automorphism has no such polynomial logarithm. None of these results resolves the planar Jacobian conjecture or improves a Hadwiger--Nelson chromatic bound.

## Exact controls

Run:

    python3 check_geometry_session2.py

Result: 100 exact checks pass. The JSON records rational shear profiles for orders N=2 through 7, and the degree-seven Hamiltonian logarithm of the Hénon control. The finite checks are supplemental; the general results are proved in the notes.

## Files

- planar_autonomous_polynomiality_beyond_fixed_line.md
- finite_dimensional_hamiltonian_lie_obstruction.md
- finite_symplectic_jets_and_boundary_control.md
- check_geometry_session2.py
- geometry_session2_exact_controls.json
- manifest.json

The main theorem was read by the coordinating agent before final additions; final common-invariant/Lie/jet text is frozen for its review. No historical-novelty claim is made.
