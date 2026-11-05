# Reviewed full-response binary-gadget inverse

The [proof](FULL_GADGET_INVERSE.md) gives every positive-integer solution for a fixed readout A g(L)+B L+C Q(L)+D Z(L)-E parity(L), under its stated coefficient inequalities. It uses dyadic shells and strictly increasing low-word branches despite the bit-reversal response Q. The [independent audit](audit/REPORT.md) passes without correction.

The exact bound is O(sqrt(K) log(2+K)) arithmetic operations, K=floor(T/B). It is not polynomial in input bit length. This algorithm does not choose physical sources, prove coefficient realizability or establish scalar coverage. Later graph applications are not silently included in this checkpoint.

Run the source normally or with -O:

    python3 verify_full_gadget_inverse.py
    python3 -O verify_full_gadget_inverse.py

It rewrites its certificate identically; use a copy for a pristine checkout. The two hash-pinned proof dependencies are included. The complete arithmetic argument is restated in the new proof.

The independent checker imports no source code. From audit/:

    g++ -std=c++17 -O2 full-gadget-inverse/check.cpp -o /tmp/s197-full-inverse-audit
    /tmp/s197-full-inverse-audit
    g++ -std=c++17 -O1 -D_GLIBCXX_ASSERTIONS -fsanitize=undefined full-gadget-inverse/check.cpp -o /tmp/s197-full-inverse-audit-debug
    /tmp/s197-full-inverse-audit-debug

The unchanged source snapshot used by that audit is preserved under audit/full-gadget-inverse/source/. All executable dependencies use standard libraries. FILES.sha256 binds every released file other than itself.
