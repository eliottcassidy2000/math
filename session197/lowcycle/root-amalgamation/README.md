# Reviewed three-cycle root-amalgamation templates

The [proof](ROOT_AMALGAMATION.md) gives exact ordered shared-root path/cover identities and four bilinear strong C=3 constructions on two physical start-one sources. It includes the precise singleton boundaries and a containment theorem for the canonical triangle-tail branch. It does not classify every C3 graph, prove uniqueness of these presentations, or establish scalar coverage.

The [independent audit](audit/REPORT.md) passes on the corrected snapshot, including the closing-only singleton suffix boundary. Its standalone implementation checks all 81,225 rooted labelled input pairs through order four, 72,092 allowed C3 templates and literal canonical tail reconstructions. The primary verifier independently enumerates directed words rather than subset-DP states.

Run in a copy because the primary verifier writes its certificate:

    python3 verify_root_amalgamation.py
    python3 -O verify_root_amalgamation.py

The independent audit is standalone C++20:

    cd audit
    g++ -std=c++20 -O2 check.cpp -o /tmp/s197-root-amalgamation
    /tmp/s197-root-amalgamation
    g++ -std=c++20 -O0 -D_GLIBCXX_ASSERTIONS check.cpp -o /tmp/s197-root-amalgamation-debug
    /tmp/s197-root-amalgamation-debug

Its unchanged source replay snapshot is included under audit/source/. All executable dependencies use standard libraries; no external graph or census file is needed. ROOTED.md and START_THREE.md are exact copies of the already-reviewed named mathematical context. The new identities, physical construction and branch-containment argument are restated in the new proof, and the copied context's broader enumeration claims are not being rerun by this packet. FILES.sha256 binds every released file except itself.
