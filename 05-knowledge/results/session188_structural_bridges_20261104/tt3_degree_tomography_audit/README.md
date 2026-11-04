# Reproduce the independent TT3 tomography review

Run `python3 audit.py > normal.json` and `python3 -O audit.py > optimized.json`, then compare. The script imports no source implementation. It reads pinned complete data from the sibling `tt3_degree_tomography` packet, uses the installed C++17 compiler temporarily, and installs nothing.

The review passes the exact coefficient identity, equal-log-concave positivity and the proposed source-relative root-family integration. It independently recovers all5,040 actual quartics, all4,941 exterior skew controls, all3,456 binomial records, both complete physical response banks and80 additional full root-family controls. It does not claim arbitrary-prime completeness.

Exact audit allowlist: `REVIEW.md`, `README.md`, `audit.py`, `count.cpp`, `run_covers.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`.
