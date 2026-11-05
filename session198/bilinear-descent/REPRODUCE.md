# Reproduction

Python 3 standard library only. From this directory:

    python3 verify.py
    python3 -O verify.py
    python3 verify_growth_independently.py
    python3 -O verify_growth_independently.py

The primary normal and optimized receipts must match. The independent checker
reconstructs generator matrices from the old physical response laws and all
bilinear matrices from the inherited scalar formulas. It does not import the
new verifier's matrices.

`audit/AUDIT.md` records independent mathematical review and isolated replays.
The frozen C2 and C3 inputs and their graph helpers are in `dependencies/`.
`RELEASE_FILES.sha256` is the original audited source manifest;
`PUBLICATION_FILES.sha256` covers this complete publication bundle.

Scope: exact physical source membership, positive bilinear normal forms,
represented-value halving descent, and restricted C2/C3 NP membership.
Universal odd-target coverage remains open.
