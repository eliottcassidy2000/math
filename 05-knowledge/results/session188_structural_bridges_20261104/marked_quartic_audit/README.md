# Marked-quartic independent audit

Run `python audit.py > normal.json` and `python -O audit.py > optimized.json`, then compare. A standard C++17 compiler is needed for the independent predecessor-DP literal controls; the build lives in a temporary directory. No software is installed.

Publication allowlist: REVIEW.md, README.md, audit.py, literal.cpp, normal.json, optimized.json, SHA256SUMS. The audit itself imports no producer code. Producer replay needs sibling marked_quartic and its covered_response/response helper dependencies. Exact source pins are in both receipts.
