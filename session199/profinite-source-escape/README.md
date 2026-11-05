# Profinite source escape

The packet proves that one fixed triangle source already gives every
odd-compatible C2 residue modulo every modulus. The primitive nonphysical
profile `(1,0,6,1)` is a limit of actual primitive source profiles modulo all
moduli simultaneously, yet its formal target 67 is not in the C2 spectrum.
The physical approximants have an explicit ordinary-height escape bound.

This is an all-moduli local/global obstruction, not a coverage proof.
The finite-target inverse and a bounded-height compactness equivalence say
exactly what additional control is missing.

Run:

    python3 verify.py > replay-normal.json
    python3 -O verify.py > replay-optimized.json
    cmp replay-normal.json replay-optimized.json
    cmp replay-normal.json VALIDATION_NORMAL.json

All tests use explicit exception gates and Python's standard library.
The inherited proof snapshots are frozen under `dependencies/`; their hashes
are printed in the receipts. The historical Collatz and Mahler snapshots are
mechanism provenance, not open-conjecture dependencies.

New claims and scope are in `PROOF.md`. No files in the dirty math checkout
were modified. Publication is coordinated by the parent session.
