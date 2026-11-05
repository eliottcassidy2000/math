# Reviewed backward start-one profile inverse

The [proof](BACKWARD_PROFILE_INVERSE.md) gives a complete all-order backward arithmetic membership test for rooted start-one profiles, with physical graph recipes. It includes the exact last-step endpoint-ratio test and a nonlinear recovery of the root-deletion count. The [independent audit](audit/REPORT.md) passes without correction.

The scalar C1 criterion still has to choose endpoint and cover responses. No cofinite scalar-coverage theorem, arbitrary-root inverse, or claim that all cone points are physical is made.

Run the frozen source normally or with -O:

    python3 verify_backward_profiles.py
    python3 -O verify_backward_profiles.py

For the independent physical/inverse audit, run from audit/:

    g++ -O3 -std=c++17 physical.cpp -o /tmp/s197-backward-physical
    /tmp/s197-backward-physical c1-inputs.csv physical-profiles.csv
    python3 check.py
    python3 -O check.py
    /tmp/s197-backward-physical c1-inputs.csv physical-profiles-replay.csv recipes.csv

The source verifier rewrites its certificate; the independent checker rewrites recipes.csv identically. Use a copy if preserving a pristine checkout matters. The complete small physical input bank is included, not fetched at runtime. The larger complete-bank generators and their completeness proofs are additionally included under audit/bank-source/: compile verify_c3_complete.cpp and run it with a fresh source CSV path; compile check_c3_complete.cpp and give it a fresh independent CSV path followed by the source CSV. Selecting its C1 rows of order at most9 reproduces c1-inputs.csv. The full bank need not be stored in the checkpoint.

The four large full-profile inversions are replayed through independently counted gadget constants, not by direct adjacency DP on their larger graphs. All executable dependencies use standard libraries. FILES.sha256 binds every bundled byte; no navigation or index is changed.
