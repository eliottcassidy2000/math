# Blue forest constraints

Read `BLUE_FOREST_CONSTRAINTS.md` for the complete derivations and scope.

Main results:
- Exact blue and black score-energy and degree-two rooted-forest moments.
- A full blue forest polynomial expressed by compatible reflected-edge
  orbit constraints, not an unjustified mean-determinant substitution.
- An additional exact degree-three moment from reversal anticommutation.
- A minimal order-six blue equal-energy pair separated by its degree-four
  forest coefficient, with explicit adjacency matrices and permutation audit.
- Correct fundamental paired-black versus literal-black counting factors.

The endpoint/triangle parity and loop necessity laws are inherited from
THM-785, THM-790, and THM-648, with exact source slugs in the main note.
The new degree-two forest statistic is explicitly identified with the old
cyclic-triangle count plus a constant. No novelty or full-wiring claim.

Reproduce from the workspace root:

    python3 session31-blue-forests/check_blue_forests.py
    python3 session31-blue-forests/check_wiring.py

`probe.py` contains shared exact matrix routines. `verification.json` and
`wiring_verification.json` hold machine-readable finite-exact receipts.
`SHA256SUMS` binds the final mathematical note, code, and receipts.
