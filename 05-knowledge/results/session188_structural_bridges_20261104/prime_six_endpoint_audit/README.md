# Reproduce the independent six-quotient audit

This auditor independently enumerates Hamilton words and all completions of their forced adjacent edges, rather than reusing the producer's per-tournament path recurrence. It checks every integer entry in the complete table, all module witnesses and prime relabeling orbits, then reconstructs every actual physical graph and full cubic/quartic response control.

Run `python3 audit.py > normal.json` and `python3 -O audit.py > optimized.json`, then compare. The frozen sibling `prime_six_reconnaissance` packet supplies pinned source data; no producer implementation is imported. Only Python's standard library and the installed C++17 compiler are required.

Exact allowlist: `REVIEW.md`, `README.md`, `audit.py`, `count.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`. This review validates the six-name theorem and its stated inherited consequences, not a general even-order claim.
