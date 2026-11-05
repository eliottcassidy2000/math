# Reproduce the complete C2 two-source normal form

The main proof, source, support and frozen dependency files are unchanged from the independently reviewed snapshot. RELEASE_FILES.sha256 pins them. The independent audit is in audit/REPORT.md. PUBLICATION.sha256 pins all delivered files.

Use a temporary copy for the main replay, since it writes SUPPORT.json and VALIDATION.json beside the verifier:

    cp -R . /tmp/session197-bilinear-c2-replay
    cd /tmp/session197-bilinear-c2-replay
    python3 verify.py
    python3 -O verify.py

Both receipts must match VALIDATION_NORMAL.json and VALIDATION_OPTIMIZED.json. The verifier checks all four bundled proof dependency hashes. It uses only the Python standard library.

The independent audit checkers are standalone and import no main source:

    g++ -std=c++17 -O2 audit/check.cpp -o /tmp/session197-bilinear-c2-check
    /tmp/session197-bilinear-c2-check
    python3 audit/check_algebra.py
    python3 -O audit/check_algebra.py

Compare their JSON outputs with the corresponding audit receipts. The C++ output includes all bounded supported scalar values. The main support records contain actual source recipes and can be realized using bilinear.realize. bilinear.factor_graph performs the complete physical graph decomposition on promised strong C2 inputs.

The theorem covers all nonregular strong C2 graphs via a unique ordered pair of nontrivial rooted start-one sources and one of two physical forms; the regular-five graph is explicit. Its numerical-target inverse is complete for this class, but it is not a cofinite scalar-spectrum theorem or a returning start-two interface.
