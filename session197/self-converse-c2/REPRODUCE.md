# Reproduce the self-converse C2 classification

The main proof, code, receipts and three proof dependencies are unchanged from the independently audited snapshot. RELEASE_FILES.sha256 pins those files. PUBLICATION.sha256 pins the complete delivery. The independent audit is in audit/REPORT.md.

Run the main verifier in a temporary copy because it writes VALIDATION.json beside itself:

    cp -R . /tmp/session197-self-converse-c2-replay
    cd /tmp/session197-self-converse-c2-replay
    python3 verify.py
    python3 -O verify.py

Both runs must match VALIDATION_NORMAL.json and VALIDATION_OPTIMIZED.json. The source needs only the Python standard library. The exact inverse is converse.diagonal_inverse(target); a success includes actual graph rows.

Run the independent checkers:

    g++ -std=c++17 -O2 audit/check.cpp -o /tmp/session197-self-converse-c2-check
    /tmp/session197-self-converse-c2-check
    python3 audit/check_algebra.py
    python3 -O audit/check_algebra.py

The algebra checker has unchanged copies of its three comparison imports in audit/source/. It computes its own gadget constants and physical response bank before comparing with the source oracle/inverse. The C++ checker imports no source. Compare outputs with the audit receipts.

The proof uses ordered rooted graph classes for the converse criterion. Equality of response quadruples is sufficient for existential scalar realization, but not for self-converseness of a particular two-source graph; the H3321 hostile tests this distinction. The result gives an all-order class and inverse, not cofinite scalar coverage.
