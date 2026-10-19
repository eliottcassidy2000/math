# Metabelian word and jet checkpoint

19 October 2026. A self-contained continuation of the normalized planar Hamiltonian Lie-algebra example.

For the two fixed Hamiltonians x²y and x³ over any characteristic-zero field:

- The exact finite scalar-flow group has only quadratic vertical shears as polynomial elements.
- Its point-adic closure contains all higher vertical shears, even though those flows cannot be obtained as finite words.
- Approximating the cubic vertical shear through degree N requires at least N−2 nonzero fourth-order pole sites.
- The exact minimum scalar-flow word length is 2N−3 for N≥4, with explicit rational coefficients attaining it. The low-order exceptions are recorded in the proof.
- The necessary denominator degree is at least 4(N−2), so increasingly accurate rational approximants cannot retain bounded pole complexity.

These are results about the stated generating pair, not a reduction or solution of the classical Jacobian conjecture. The finite poles are affine; the result does not apply the earlier inverse-at-infinity criterion for polynomial automorphisms. Earlier source principal-part and closure results are linked and distinguished. No historical-novelty claim is made.

## Files and verification

- `metabelian_word_closure_and_sharp_jet_length.md`: all-order proofs and scope.
- `check_metabelian_jets.py`: exact symbolic and rational finite controls.
- `metabelian_jets_exact_controls.json`: 166 passing checks and explicit shortest-word examples through N=14.

Run from this directory:

    python3 check_metabelian_jets.py

SymPy is required. The finite checks support the explicit identities; general minimality follows from the moment-functional and accumulated-time-walk proofs in the note.
