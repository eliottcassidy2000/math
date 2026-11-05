# Reviewed fixed-port triangle gate

The [proof](C3_TRIANGLE_GATE.md) gives an all-order marked triangle replacement law, a physical two-parameter C3 family using the same two original ports, and an exact finite inverse for that named family. The [independent audit](audit/REPORT.md) passes without correction.

The family realizes67,533 and1861 explicitly. It does not realize821 and is not claimed to exhaust all C3/feedback-two tournaments. The separate H533 threshold proof is copied as referenced context; its complete negative-census reproduction and audit are at [the published H533 checkpoint](https://github.com/eliottcassidy2000/math/blob/231bda1d4b2f57506bfcd6a2b2b3decb5f5f089d/session197/lowcycle/feedback533/README.md).

From this directory run:

    python3 verify_c3_triangle_gate.py
    python3 -O verify_c3_triangle_gate.py
    g++ -O3 -std=c++17 audit/check.cpp -o /tmp/s197-triangle-gate-check
    /tmp/s197-triangle-gate-check > /tmp/s197-triangle-gate-report.txt 2> /tmp/s197-triangle-gate-image.txt

The two temporary outputs must match audit/independent-normal.txt and audit/image.txt. The C++ checker is standalone and independently counts literal path words, two-covers and increasing path slots; it does not import the Python implementation. The source verifier rewrites its adjacent certificate identically, so use a copy when preserving a pristine checkout matters. All dependencies use standard libraries. FILES.sha256 binds the frozen bytes; no navigation or index file is modified.
