# Reviewed square-root final-gadget inverse

The [proof](SQRT_FIXED_PREFIX_INVERSE.md) gives every final-gadget solution for a fixed response triple in O(sqrt(N/R) log(N/R)) arithmetic operations. It uses an exact binary-block identity and strictly increasing low-word branches. The [independent audit](audit/REPORT.md) passes without correction.

The prefix must already be physically available to turn a returned parameter into a graph. This does not prove cofinite scalar coverage or choose a reachable prefix for every target, and its complexity is not polynomial in the input bit length.

Run the source normally or with -O:

    python3 verify_sqrt_fixed_prefix.py
    python3 -O verify_sqrt_fixed_prefix.py

The verifier rewrites its certificate identically. Use a copy to preserve a pristine checkout. Its only dependencies are the two included, hash-pinned proof files.

The separate audit imports no source code. From audit/:

    g++ -std=c++17 -O3 sqrt-fixed-prefix/check.cpp -o /tmp/s197-sqrt-audit
    /tmp/s197-sqrt-audit
    g++ -std=c++17 -O1 -D_GLIBCXX_ASSERTIONS -fsanitize=undefined sqrt-fixed-prefix/check.cpp -o /tmp/s197-sqrt-audit-debug
    /tmp/s197-sqrt-audit-debug

It independently counts gadgets by before/after-root partitions, checks all three large exact inverses, and counts224 actual graphs by walk inclusion-exclusion. All executable dependencies use standard libraries. FILES.sha256 binds the released bytes.
